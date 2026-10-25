"""Перекомбинация множителей по решётке следов — дорога ван Хоя.

Стена, названная §441 и сдвинутая §454: перебор подмножеств у
Зассенхауса экспоненциален по числу модульных множителей r, а дорога
LLL 1982 (`lllfactor.py`) полиномиальна, но на замеренном участке дороже
перебора на четыре порядка. Практическая дорога за порогом — ван Хой
[внешний источник: M. van Hoeij, «Factoring polynomials and the knapsack
problem», J. Number Theory 95 (2002); K. Belabas, «A relative van Hoeij
algorithm over number fields», 2004].

Устройство. Пусть f унитарен, свободен от квадратов, степени n, и
g₁…g_r — его множители по модулю p, поднятые Гензелем до p^a. Всякий
делитель f над ℤ есть произведение НЕКОТОРОГО подмножества S ⊆ {1…r}.
Ключ ван Хоя: у такого произведения степенные суммы корней —
ЦЕЛЫЕ и МАЛЫЕ, |P_k| ≤ n·B^k, где B — граница модуля корней f. По модулю
p^a они равны Σ_{i∈S} T_{i,k}, где T_{i,k} — степенные суммы корней g_i.

Отсюда решётка: вектор (x, y) ∈ ℤ^r × ℤ^s лежит в ней, если
Σ xᵢ·T_{i,k} + y_k·p^a — тот самый малый остаток. Базис:

    [ I_r | T ]
    [  0  | p^a·I_s ]

Приведение по LLL вытаскивает векторы с малым остатком; среди них
подмножества и живут. Перебора нет: решётка размерности r + s, а не 2^r.

Почему s берётся больше r. При s < r условий меньше, чем неизвестных, и
решение перестаёт быть единственным: замер 04.09 на степени 64 (r = 32,
s = 10) дал 0 годных векторов — не оттого, что их нет, а оттого, что
пространство решений было двадцатидвумерным и ни один вектор базиса не
оказался нулём-единицей. При s = r + 4 та же задача решается.

**Всякий найденный набор проверяется делением.** Ложная находка потому
безвредна: она просто не поделит f. Опасен лишь пропуск, и на этот
случай дорога возвращает `None`, а вызывающий берёт перебор — молчаливой
подмены нет.
"""
from __future__ import annotations

from .intfactor import _choose_prime, _hensel_lift, _mmul, _norm, _symmetric, _zdivmod
from .intpoly import IntPoly
from .lattice import lll

EXTRA_COLUMNS = 4        # s = r + EXTRA_COLUMNS: условий больше, чем неизвестных
SLACK_BITS = 20          # запас точности сверх границы малости остатков


def _iroot(x: int, k: int) -> int:
    """Целая часть корня k-й степени — двоичным поиском, без float."""
    if x <= 0:
        return 0
    lo, hi = 1, 1 << ((x.bit_length() // k) + 2)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if mid ** k <= x:
            lo = mid
        else:
            hi = mid - 1
    return lo


def root_bound(fc: list) -> int:
    """Граница модуля корней по Фудзиваре — целыми, всегда сверху.

    Коши (1 + max|aᵢ|) на этом семействе слишком груба: для
    Суиннертон-Дайера степени 8 она даёт 962 против 14 у Фудзивары, а
    граница входит в показатель s — то есть в разрядность p^a.
    """
    n = len(fc) - 1
    an = abs(fc[-1])
    best = 0
    for i in range(n):
        c = abs(fc[i])
        k = n - i
        best = max(best, _iroot((c * 2 if i == 0 else c) // an + 1, k) + 1)
    return 2 * best


def newton_sums(g: list, s: int, mod: int) -> list:
    """Первые s степенных сумм корней унитарного g по модулю mod.

    Рекуррента Ньютона: p_k + a_{d−1}p_{k−1} + … + k·a_{d−k} = 0.
    Коэффициенты с отрицательным номером — нули, и это не «на всякий
    случай»: при s > deg g они и встречаются.
    """
    d = len(g) - 1
    a = [c % mod for c in g]
    out: list = []
    for k in range(1, s + 1):
        acc = 0
        for i in range(1, k):
            j = d - i
            if j >= 0:
                acc += a[j] * out[k - i - 1]
        if k <= d:
            acc += k * a[d - k]
        out.append((-acc) % mod)
    return out


def knapsack_lattice(traces: list, mod: int) -> list:
    """Базис решётки следов: [I_r | T] сверху и [0 | p^a·I_s] снизу."""
    r, s = len(traces), len(traces[0])
    rows = [[1 if j == i else 0 for j in range(r)] + list(traces[i]) for i in range(r)]
    rows += [[0] * r + [mod if j == k else 0 for j in range(s)] for k in range(s)]
    return rows


def _subsets(reduced: list, r: int, limit: int) -> list:
    """Годные наборы из приведённого базиса: нуль-единичные с малым остатком."""
    out = []
    for v in reduced:
        head, tail = v[:r], v[r:]
        if max((abs(t) for t in tail), default=0) >= limit:
            continue
        if all(x == 0 for x in head):
            continue
        for sign in (1, -1):
            h = [sign * x for x in head]
            if all(x in (0, 1) for x in h):
                idx = tuple(i for i, x in enumerate(h) if x)
                if idx and idx not in out:
                    out.append(idx)
    return out


def factor_squarefree(p: IntPoly, extra: int = EXTRA_COLUMNS):
    """Множители f решёткой следов, либо None — тогда берётся перебор.

    None означает ровно одно: дорога не собрала разбиение целиком. Это
    не «неприводим» и не «ошибка»: вызывающий обязан взять другую дорогу.
    """
    fc = _norm(list(p.c))
    n = len(fc) - 1
    if n <= 1 or fc[-1] != 1:
        return None
    prime, modular = _choose_prime(fc)
    r = len(modular)
    if r == 1:
        return [IntPoly(fc)]
    B = root_bound(fc)
    s = r + extra
    limit = n * B ** s
    mod, lifted = _hensel_lift(fc, modular, prime, limit << SLACK_BITS)
    traces = [newton_sums(_symmetric(g, mod), s, mod) for g in lifted]
    reduced = lll(knapsack_lattice(traces, mod))
    found, used, rest = [], set(), list(fc)
    for idx in sorted(_subsets(reduced, r, limit), key=len):
        if set(idx) & used:
            continue
        g = [1]
        for i in idx:
            g = _mmul(g, lifted[i], mod)
        g = _symmetric(g, mod)
        if g[-1] != 1:
            continue
        q, rem = _zdivmod(rest, g)
        if rem == [0]:
            found.append(IntPoly(g))
            rest = q
            used |= set(idx)
    if len(used) == r:
        return found
    # Осталось неразобранное. Один неиспользованный множитель — остаток и
    # есть он сам, делить дальше нечего. Несколько — задача та же, но
    # меньше: зовём себя на остатке. Пустая находка означает, что дорога
    # не сдвинулась ни на шаг, и тогда честнее вернуть None, чем крутиться.
    remaining = [i for i in range(r) if i not in used]
    if len(rest) - 1 <= 0:
        return found if found else None
    if len(remaining) == 1:
        found.append(IntPoly(rest))
        return found
    if not found:
        return None
    sub = factor_squarefree(IntPoly(rest), extra)
    return None if sub is None else found + sub


# ---------- самоприменение ----------

# Окно приёмки для самопроверки: многочлен, число множителей над ℤ.
# `[вычислено здесь]` Все четыре — из приёмки proto/test_vanhoeij.py;
# первый неприводим над ℤ (дорога возвращает его самого), остальные
# дорога собирает целиком.
_WINDOW = (((1, 0, -10, 0, 1), 1), ((-1, 0, 0, 1), 2),
           ((-6, 11, -6, 1), 3), ((6, 5, 1), 2))
_ROOTS = (((-6, 11, -6, 1), (1, 2, 3)),      # (x−1)(x−2)(x−3)
          ((6, 5, 1), (-2, -3)),             # (x+2)(x+3)
          ((24, -50, 35, -10, 1), (1, 2, 3, 4)))   # (x−1)(x−2)(x−3)(x−4)


def self_check(deep: bool = False) -> list:
    """Самопроверка дороги ван Хоя ВСТРЕЧНЫМИ дорогами (ось самоприменения).

    `[вычислено здесь]` Дорога собрана из четырёх частей, и каждая
    сверяется не собой: (1) степенные суммы Ньютона — против ПРЯМОГО счёта
    Σ rᵢᵏ на многочленах с известными целыми корнями и против АДДИТИВНОСТИ
    (суммы произведения равны сумме сумм множителей — тот самый ключ ван
    Хоя, без которого решётки нет); (2) граница корней Фудзивары — что она
    СВЕРХУ настоящего наибольшего корня, и что `_iroot` есть целая часть
    корня: k-я степень не превосходит x, а у соседа превосходит; (3)
    решётка следов — после приведения всякий вектор лежит в ней:
    хвост сравним с Σ xᵢ·Tᵢ по модулю p^a (членство решается подстановкой,
    не приведением); (4) итог — против ПЕРЕБОРА `intfactor`: то же
    разбиение, произведение равно входу, каждый множитель унитарен, а на
    окне приёмки число множителей совпадает с записанным — отступление
    `None` там, где дорога прежде шла, есть находка.

    Граница объявлена: `None` на многочлене вне окна — честный переход на
    перебор, не находка; неунитарный вход и степень ≤ 1 обязаны давать
    `None` по договору. При `deep` окно корней расширяется степенью 6 и
    суммы сверяются до s = 12. Возвращает список находок; исключений не
    бросает.
    """
    from . import intfactor as _Z
    out: list = []
    M = 1 << 61
    s_max = 12 if deep else 6

    def sym(v):
        return [x - M if x > M // 2 else x for x in v]

    # (1) суммы Ньютона против прямого счёта и аддитивности
    roots_window = list(_ROOTS)
    if deep:
        roots_window.append(((720, -1764, 1624, -735, 175, -21, 1), (1, 2, 3, 4, 5, 6)))
    for fc, roots in roots_window:
        try:
            got = sym(newton_sums(list(fc), s_max, M))
        except Exception as e:                                  # noqa: BLE001
            out.append(f"суммы Ньютона отказали на {fc}: {type(e).__name__}: {e}")
            continue
        direct = [sum(r ** k for r in roots) for k in range(1, s_max + 1)]
        if got != direct:
            out.append(f"суммы Ньютона расходятся с прямым счётом по корням {roots}: "
                       f"{got} ≠ {direct}")
    g, h = [-1, 1], [-2, 1]                                   # (x−1), (x−2)
    gh = _mmul(g, h, M)
    try:
        sg, sh, sgh = (sym(newton_sums(q, s_max, M)) for q in (g, h, gh))
        if [a + b for a, b in zip(sg, sh)] != sgh:
            out.append(f"аддитивность сумм Ньютона нарушена: {sg} + {sh} ≠ {sgh}")
    except Exception as e:                                      # noqa: BLE001
        out.append(f"аддитивность сумм Ньютона не замерена: {type(e).__name__}: {e}")

    # (2) граница корней сверху и целый корень
    for fc, roots in roots_window:
        b = root_bound(list(fc))
        if b < max(abs(r) for r in roots):
            out.append(f"граница Фудзивары {b} ниже настоящего корня {max(abs(r) for r in roots)} у {fc}")
    for x, k in ((0, 2), (1, 3), (26, 3), (27, 3), (28, 3), (10 ** 12, 4), (2 ** 40, 5)):
        v = _iroot(x, k)
        if x > 0 and not (v ** k <= x < (v + 1) ** k):
            out.append(f"_iroot({x}, {k}) = {v} не есть целая часть корня")

    # (3) членство приведённых векторов в решётке следов
    traces = [[1, 2, 3], [4, 5, 6]]
    mod = 101
    try:
        red = lll(knapsack_lattice(traces, mod))
        r = len(traces)
        for v in red:
            head, tail = v[:r], v[r:]
            expect = [sum(head[i] * traces[i][k] for i in range(r)) % mod for k in range(len(tail))]
            if [t % mod for t in tail] != expect:
                out.append(f"вектор {v} после приведения не лежит в решётке следов")
        if len(red) != r + len(traces[0]):
            out.append(f"базис после приведения потерял размерность: {len(red)} ≠ {r + len(traces[0])}")
    except Exception as e:                                      # noqa: BLE001
        out.append(f"решётка следов не приведена: {type(e).__name__}: {e}")

    # (4) итог против перебора
    for fc, cnt in _WINDOW:
        p = IntPoly(list(fc))
        try:
            got = factor_squarefree(p)
        except Exception as e:                                  # noqa: BLE001
            out.append(f"дорога отказала исключением на {fc}: {type(e).__name__}: {e}")
            continue
        if got is None:
            out.append(f"дорога отступила (None) на {fc} из окна приёмки, где прежде шла")
            continue
        if len(got) != cnt:
            out.append(f"число множителей у {fc}: {len(got)} ≠ {cnt}")
        prod = IntPoly([1])
        for q in got:
            prod = prod * q
            if list(q.c)[-1] != 1:
                out.append(f"множитель {list(q.c)} у {fc} не унитарен")
        if list(prod.c) != list(fc):
            out.append(f"произведение множителей {list(prod.c)} ≠ входу {list(fc)}")
        by_enum = sorted(tuple(q.c) for q in _Z.factor_squarefree(p))
        if sorted(tuple(q.c) for q in got) != by_enum:
            out.append(f"разбиение {fc} расходится с перебором: {sorted(tuple(q.c) for q in got)} ≠ {by_enum}")

    # границы по договору
    if factor_squarefree(IntPoly([2, 0, 2])) is not None:
        out.append("неунитарный вход обязан давать None")
    if factor_squarefree(IntPoly([-3, 1])) is not None:
        out.append("степень 1 обязана давать None")
    return out
