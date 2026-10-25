"""Факторизация над ℤ через приведение решётки — дорога LLL 1982.

Стена, названная §441: перебор подмножеств у Зассенхауса экспоненциален
по числу модульных множителей, порог r ≤ 24 держится явным отказом, а
цена за ним замерена (r = 32 — около 2,4 часа). Дорога в обход известна с
1982 года и здесь написана.

Устройство [внешний источник: A. K. Lenstra, H. W. Lenstra Jr., L.
Lovász 1982; von zur Gathen & Gerhard, «Modern Computer Algebra», гл. 16].
Пусть f унитарен, свободен от квадратов, степени n; h — неприводимый по
модулю p множитель степени d, поднятый Гензелем до p^k. Тогда истинный
неприводимый множитель f, делящийся на h, есть КОРОТКИЙ вектор решётки

    L(m) = { g ∈ ℤ[x] : deg g ≤ m, h | g (mod p^k) },

порождённой строками p^k·x^i (i < d) и h·x^j (j ≤ m − d). Приведение по
LLL даёт короткий b₁, и если ‖b₁‖^n · ‖f‖^m < p^k, то НОД(b₁, f) —
нетривиальный множитель. Условие проверяется целыми числами, а не на
глаз: пока оно не выполнено, множитель не объявляется.

Перебора подмножеств здесь нет вовсе: множитель приходит прямо из решётки.
Цена, однако, своя, и она замерена (`proto/bench_lllfactor.py`): точная
арифметика приведения на больших размерностях дороже перебора, и место
перелома названо числом, а не обещанием.

Ненахождение здесь — вывод, а не незнание: точность p^k берётся такой,
что короткий вектор ОБЯЗАН найтись, если делитель есть. Потому «ни при
одном m» означает неприводимость, и это сказано в коде рядом с выкладкой,
а не подразумевается.
"""
from __future__ import annotations

from math import isqrt

from .errors import Unresolved
from .intfactor import (_choose_prime, _hensel_lift, _norm, _norm2_sq,
                        _zdivmod)
from .intpoly import IntPoly
from .lattice import lll, norm2


def _pseudo_rem(a: list, b: list) -> list:
    """Псевдоостаток: lc(b)^(deg a − deg b + 1)·a mod b, всё на целых."""
    a, b = list(a), list(b)
    if len(a) < len(b):
        return a
    lc = b[-1]
    r = [c * lc ** (len(a) - len(b) + 1) for c in a]
    for k in range(len(r) - len(b), -1, -1):
        c = r[len(b) - 1 + k]
        if c % lc == 0:
            q = c // lc
            for i, bc in enumerate(b):
                r[i + k] -= q * bc
    return _norm(r[: len(b) - 1] or [0])


def _content(a: list) -> int:
    from math import gcd
    g = 0
    for c in a:
        g = gcd(g, abs(c))
    return g or 1


def _primitive(a: list) -> list:
    g = _content(a)
    return _norm([c // g for c in a]) if g > 1 else _norm(a)


def poly_gcd(a: list, b: list) -> list:
    """НОД в ℤ[x] примитивной последовательностью псевдоостатков."""
    a, b = _primitive(a), _primitive(b)
    while b != [0]:
        r = _pseudo_rem(a, b)
        a, b = b, _primitive(r) if r != [0] else [0]
    return _primitive(a)


def _lattice(h: list, mod: int, m: int) -> list:
    """Базис решётки L(m): p^k·x^i (i < d) и h·x^j (j ≤ m − d)."""
    d = len(h) - 1
    rows = []
    for i in range(d):
        row = [0] * (m + 1)
        row[i] = mod
        rows.append(row)
    for j in range(m - d + 1):
        row = [0] * (m + 1)
        for i, c in enumerate(h):
            row[i + j] = c
        rows.append(row)
    return rows


def factor_one(fc: list, p: int | None = None, k_extra: int = 0) -> tuple:
    """Один неприводимый множитель f и остаток: (множитель, частное).

    f унитарен и свободен от квадратов. Точность p^k берётся по границе
    LLL: p^k > 2^(n²/2)·‖f‖^(2n); `k_extra` поднимает её ещё, если
    условие не выполнилось (это редкость, но молча не догадываемся).
    """
    fc = _norm(list(fc))
    n = len(fc) - 1
    if n <= 1:
        return list(fc), [1]
    p_used, modular = (_choose_prime(fc) if p is None else (p, None))
    if modular is None:
        from .intfactor import factor_mod_p
        modular = factor_mod_p(fc, p_used)
    if len(modular) == 1:
        return list(fc), [1]
    norm_f = isqrt(_norm2_sq(fc)) + 1
    # Точность обязана быть такой, чтобы НЕНАХОЖДЕНИЕ было доказательством
    # неприводимости, а не следствием малого p^k. Если u | f неприводим и
    # h | u (mod p^k), то u ∈ L(m) при m ≥ deg u, и LLL даёт
    # ‖b₁‖ ≤ 2^(m/2)·‖u‖ ≤ 2^(m/2)·2^n·‖f‖ (Миньотт). Условие
    # ‖b₁‖^n·‖f‖^m < p^k тогда выполнено, если
    # p^k > 2^(n·m/2 + n²)·‖f‖^(n+m); при m ≤ n − 1 хватает
    # p^k > 2^(3n²/2)·‖f‖^(2n). Запас на единицу в показателе — от
    # округлений вниз в isqrt.
    bound = (1 << (3 * n * n // 2 + n + 1)) * norm_f ** (2 * n)
    target = bound << k_extra
    mod, lifted = _hensel_lift(fc, modular, p_used, target)
    d = len(lifted[0]) - 1
    for m in range(d, n):
        basis = lll(_lattice(lifted[0], mod, m))
        b1 = _norm(list(basis[0]))
        if b1 == [0]:
            continue
        if isqrt(norm2(b1)) ** n * norm_f ** m < mod:
            g = poly_gcd(fc, b1)
            if 0 < len(g) - 1 < n:
                q, r = _zdivmod(fc, g)
                if r == [0]:
                    return g, q
            if len(g) - 1 == n:
                return list(fc), [1]
    # Ни при одном m ≤ n − 1 короткого вектора не нашлось. При взятой
    # точности это НЕ незнание, а вывод: делителя степени ≤ n − 1,
    # делящегося на h, нет — значит f неприводим.
    return list(fc), [1]


def self_check(deep: bool = False) -> list:
    """Самопроверка дороги LLL ВСТРЕЧНОЙ дорогой (ось самоприменения).

    `[вычислено здесь]` Сильнейшая встречная дорога здесь — ДРУГОЙ
    АЛГОРИТМ: тот же вход раскладывает `intfactor.factor_squarefree`
    (Зассенхаус с подъёмом Гензеля и перебором подмножеств), и два пути,
    не имеющие общего шага кроме выбора простого, обязаны дать ОДИН набор
    множителей — расхождение указывает на ошибку ровно одного из них.
    Вторая дорога — умножение обратно: произведение множителей обязано
    равняться входу (вход унитарен, потому с точностью до знака не нужно).
    Третья — степени: сумма степеней равна степени входа, и ни один
    множитель не имеет нулевой степени.

    Граница объявлена самой дверью: дорога написана для УНИТАРНОГО
    свободного от квадратов входа, и неунитарный сюда не входит — он
    отбит отказом `Unresolved` (проверяется здесь же ловлей отказа, а не
    молчанием). Пробы малой степени: цена приведения решётки растёт с
    числом модульных множителей, и полнота по всем многочленам не
    утверждается. Возвращает список находок; исключений не бросает.
    """
    from .intfactor import factor_squarefree as _zass
    out = []
    probes = [IntPoly([-1, 0, 0, 1]),                 # x³ − 1
              IntPoly([1, 0, 1]) * IntPoly([-1, 1]),  # (x²+1)(x−1)
              IntPoly([-2, 0, 1]),                    # x² − 2
              IntPoly([1, 1, 1, 1, 1])]               # Φ₅
    if deep:
        probes += [IntPoly([1, 0, 0, 0, 1]),          # x⁴ + 1
                   IntPoly([-1, 0, 0, 0, 0, 1])]      # x⁵ − 1
    for f in probes:
        try:
            mine = factor_squarefree(f)
        except Exception as e:                                  # noqa: BLE001
            out.append(f"дорога LLL отказала на {f.c}: {type(e).__name__}: {e}")
            continue
        prod = IntPoly([1])
        for g in mine:
            prod = prod * g
            if g.degree < 1:
                out.append(f"множитель нулевой степени на {f.c}: {g.c}")
        if prod != f:
            out.append(f"произведение множителей не равно входу {f.c}: "
                       f"{prod.c}")
        if sum(g.degree for g in mine) != f.degree:
            out.append(f"сумма степеней {sum(g.degree for g in mine)} ≠ "
                       f"степени входа {f.degree} на {f.c}")
        try:
            theirs = _zass(f)
        except Exception as e:                                  # noqa: BLE001
            out.append(f"встречный алгоритм отказал на {f.c}: "
                       f"{type(e).__name__}: {e}")
            continue
        key = lambda q: (len(q.c), q.c)                        # noqa: E731
        if sorted(mine, key=key) != sorted(theirs, key=key):
            out.append(f"две дороги факторизации разошлись на {f.c}: решёткой "
                       f"{[g.c for g in mine]}, Зассенхаусом "
                       f"{[g.c for g in theirs]}")
    try:
        factor_squarefree(IntPoly([1, 2, 3]))                  # старший 3 — не унитарен
        out.append("неунитарный вход принят вместо отказа")
    except Unresolved:
        pass
    except Exception as e:                                      # noqa: BLE001
        out.append(f"неунитарный вход отказал не Ф-отказом: {type(e).__name__}")
    return out


def factor_squarefree(p: IntPoly) -> list:
    """Все неприводимые множители унитарного свободного от квадратов f."""
    fc = _norm(list(p.c))
    if fc[-1] != 1:
        raise Unresolved("факторизация решёткой", repr(p),
                         hint="дорога LLL здесь написана для унитарного f")
    out = []
    rest = fc
    while len(rest) - 1 > 0:
        g, q = factor_one(rest)
        out.append(IntPoly(g))
        if len(q) - 1 <= 0:
            rest = [1]
            break
        rest = q
    return out
