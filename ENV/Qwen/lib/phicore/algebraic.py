"""Линейка класса 3 `POLY_ROOT`: точная арифметика алгебраических величин.

Величина = корень многочлена `P ∈ ℤ[x]`, изолированный рациональным
отрезком `[lo, hi]`. Явная запись: `Φ(C=3, P=…, отрезок=[lo; hi], s)`.
Знак — компонент записи: `−√3` строится с `sign=1` и НЕ равен `√3`;
`√3 − √3` даёт точный `Q(0)` (исправление замеренной потери
знака `(n^q)·radicand`).

Арифметика — результантами (через `intpoly`): сумма `Res_y(P(y), Q(x−y))`,
произведение `Res_y(P(y), y^m·Q(x/y))`; изоляция корня — цепочкой Штурма
по свободной от квадратов части; рациональные значения распознаются и
возвращаются как `Q`. Несведённое — `Unresolved`, не ложный класс.

Внутренний ход — целочисленный: цепочка Штурма — PRS на псевдоостатках
с контролем роста через content; знак p(a/b) — гомогенизацией
Σ c_i·a^i·b^(n−i) на целых; НОД многочленов — примитивная PRS с точным
целым делением. Q живёт только на границах: входы, выходы и
концы отрезков изоляции.
"""
from __future__ import annotations

from .rational import Q
from functools import lru_cache
from math import gcd
from .errors import Unresolved
from .intpoly import IntPoly, resultant_values, interpolate_int_poly

_MAX_REFINE = 200
# Предел перебора делителей при поиске рационального корня равен k —
# точности представления самого Ф-числа в этой разрядности (слово
# автора 02.09: это самый крайний случай). Разрядность здесь — битовая
# длина коэффициента; нижний этаж 64 покрывает мелкие разрядности.
_DIV_FLOOR = 64


# ---------- многочлены над ℤ: псевдоостаток, НОД, Штурм ----------

def _norm(p):
    """Снять хвостовые нули; ноль — всегда [0]. Правит список на месте."""
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    if not p:
        p.append(0)
    return p


# целый корень живёт в enclosure (сертификат корня из отрезка — там же)
from .enclosure import int_root_floor as _int_root_floor  # noqa: E402


def _content(p):
    g = 0
    for a in p:
        g = gcd(g, a)
    return g or 1


def _prim(p):
    """Деление на содержание: рост коэффициентов снят, знаки сохранены."""
    g = _content(p)
    return [a // g for a in p] if g != 1 else list(p)


def _prem(a, b):
    """Псевдоостаток над ℤ: r = lc(b)^s · (a mod b), дробей нет.
    Возвращает (r, s) — множитель нужен восстановлению знака остатка."""
    r = _norm(list(a))
    db = len(b) - 1
    lb = b[-1]
    steps = 0
    while len(r) - 1 >= db and any(r):
        lead = r[-1]
        shift = len(r) - 1 - db
        r = [c * lb for c in r]
        for i, bc in enumerate(b):
            r[shift + i] -= lead * bc
        r.pop()  # старший погашен точно
        _norm(r)
        steps += 1
    if not any(r):
        r = [0]
    return r, steps


def _int_poly_gcd(a, b):
    """НОД в ℤ[x] примитивной PRS: псевдоостаток + снятие content на
    каждом шаге; результат примитивный, старший коэффициент > 0.
    Замер 04.09 (набор Zech): 173 вызова на 59 различных пар — итог
    кэшируется по коэффициентам; субрезультантная PRS без content
    проверялась и оказалась медленнее (5,0 с против 2,2 с) — снята."""
    return list(_int_poly_gcd_cached(tuple(a), tuple(b)))


@lru_cache(maxsize=512)
def _int_poly_gcd_cached(a, b):
    a, b = _prim(_norm(list(a))), _prim(_norm(list(b)))
    while b != [0]:
        r, _ = _prem(a, b)
        a, b = b, _prim(_norm(r))
    if a[-1] < 0:
        a = [-c for c in a]
    return tuple(a)


def _divexact(a, g):
    """Точное деление в ℤ[x]: a = g·q с целым q (примитивный делитель
    примитивной части — гарантия Гаусса); каждый шаг делится нацело."""
    a, g = _norm(list(a)), _norm(list(g))
    q = [0] * (len(a) - len(g) + 1)
    r = list(a)
    lg = g[-1]
    for k in range(len(q) - 1, -1, -1):
        c = r[len(g) - 1 + k] // lg
        q[k] = c
        if c:
            for i, gc in enumerate(g):
                r[i + k] -= c * gc
    return q


def squarefree_part(p: IntPoly) -> IntPoly:
    """P / gcd(P, P′) — те же корни без кратности; примитивный из ℤ[x].
    Весь ход целочисленный: НОД — примитивная PRS, деление — точное."""
    g = _int_poly_gcd(p.c, p.derivative().c)
    if len(g) == 1:  # взаимно просты — уже свободен от квадратов
        return p.primitive()
    return IntPoly(_divexact(p.c, g)).primitive()


def _homog_sign(c, a, b):
    """Знак p(a/b) при b > 0 гомогенизацией: sign Σ c_i·a^i·b^(n−i).
    Горнер на целых, ни одной дроби."""
    v = c[-1]
    p = 1
    for i in range(len(c) - 2, -1, -1):
        p *= b
        v = v * a + c[i] * p
    return (v > 0) - (v < 0)


def _sturm_chain(p: IntPoly):
    """Целая цепочка Штурма: p, p′, дальше −остаток с точностью до
    положительного множителя (псевдоостаток с восстановлением знака,
    рост коэффициентов снимается content). Замер 04.09 (набор Zech):
    225 вызовов на 53 различных многочлена, до 25 повторов — цепочка
    кэшируется по коэффициентам (только читается: _sign_changes)."""
    return _sturm_chain_cached(tuple(p.c))


@lru_cache(maxsize=512)
def _sturm_chain_cached(coeffs):
    p = IntPoly(list(coeffs))
    chain = [list(p.c), list(p.derivative().c)]
    if chain[1] == [0]:
        raise ZeroDivisionError
    while True:
        a, b = chain[-2], chain[-1]
        r, steps = _prem(a, b)
        if r == [0]:
            break
        # истинный остаток = r / lc(b)^steps; элемент цепочки обязан быть
        # −остаток·(положительное) ⇒ при lc(b)^steps < 0 знак r уже нужный
        if b[-1] < 0 and steps % 2 == 1:
            nxt = _prim(r)
        else:
            nxt = [-c for c in _prim(r)]
        chain.append(nxt)
        if len(nxt) == 1:
            break
    return chain


def _sign_changes(chain, x: Q) -> int:
    a, b = x.numerator, x.denominator
    signs = []
    for p in chain:
        s = _homog_sign(p, a, b)
        if s != 0:
            signs.append(s)
    return sum(1 for u, v in zip(signs, signs[1:]) if u != v)


def count_roots(p_sf: IntPoly, lo: Q, hi: Q) -> int:
    """Число различных вещественных корней p_sf на (lo, hi]. Штурм."""
    chain = _sturm_chain(p_sf)
    return _sign_changes(chain, lo) - _sign_changes(chain, hi)


# ---------- поиск рационального корня на отрезке ----------

def _divisors(n: int):
    """Делители через факторизацию: предел перебора — разрядность k
    самого числа, не его величина (дефект найден тождеством
    (p/q)^(a/b·b/a): 5¹² — большое число с 13 делителями). Пробное
    деление до 10⁴. Остаток m > 1, чью простоту пробное деление не
    установило (порог 10⁴ достигнут при d² ≤ m), перечня НЕ даёт — ответ
    None, как и при переполнении перечня. До 18.09 такой остаток входил
    «атомом», и перечень 10002200057 = 100003·100019 выходил из 2
    делителей вместо 4 — ложный класс у корня x² − 200022x + 10002200057
    на [0; 100010]: as_fraction() = None при == 100003 → True
    `[вычислено здесь]`. При None вызывающие стороны идут дорогой
    факторизации над ℤ (_rational_root_by_factoring)."""
    n = abs(n)
    if n == 0:
        return None
    limit = max(_DIV_FLOOR, n.bit_length())
    fac = {}
    d, m = 2, n
    while d * d <= m and d <= 10 ** 4:
        while m % d == 0:
            fac[d] = fac.get(d, 0) + 1
            m //= d
        d += 1 if d == 2 else 2
    if m > 1:
        if d * d <= m:
            return None          # порог 10⁴ достигнут, простота m не установлена
        fac[m] = fac.get(m, 0) + 1
    divs = [1]
    for p, e in fac.items():
        divs = [dv * p ** k for dv in divs for k in range(e + 1)]
        if len(divs) > limit:
            return None
    return sorted(divs)


def iroot(n: int, q: int) -> int:
    """⌊n^(1/q)⌋ для n ≥ 0, q ≥ 1 — ТОЧНО, целым Ньютоном.

    Ни одного float: дорога через `n ** (1.0/q)` на 10^20 и q = 2 врёт на
    ±1 (мантисса IEEE короче числа), а здесь ответ проверяем тождеством
    x^q ≤ n < (x+1)^q. Начальное приближение — 2^⌈bits/q⌉ ≥ корня, сверху,
    потому итерация монотонно убывает и останов по y ≥ x точен."""
    if n < 0:
        raise ValueError("iroot: отрицательное подкоренное")
    if q < 1:
        raise ValueError("iroot: показатель корня < 1")
    if q == 1 or n < 2:
        return n
    x = 1 << (-(-n.bit_length() // q))
    while True:
        y = ((q - 1) * x + n // x ** (q - 1)) // q
        if y >= x:
            return x
        x = y


def exact_root(r: Q, q: int):
    """Точный q-й корень из ПОЛОЖИТЕЛЬНОЙ пары — парой, или None.

    Прямая формула вместо поиска: пара в каноне (gcd = 1), потому
    (n/d)^(1/q) рационально тогда и только тогда, когда n и d — полные
    q-е степени. До 12.09 тот же вопрос решался построением многочлена
    d·x^q − n и ЛИНЕЙНЫМ сканом целых до самой величины: Q(10**8)**Q(1,3)
    стоило 49 076 мс, Q(10**20)**Q(1,2) не досчитывалось вовсе."""
    num, den = r.numerator, r.denominator
    if num <= 0:
        raise ValueError("exact_root: подкоренное неположительно")
    rn, rd = iroot(num, q), iroot(den, q)
    if rn ** q == num and rd ** q == den:
        return Q(rn, rd)
    return None


# Бюджет линейного скана целых в _try_rational — ОБЪЯВЛЕННОЕ число шагов.
# Замер 12.09 на всём наборе соседей (12 приёмок): ДО правки 20 479
# вызовов, 331 021 шаг, самый длинный один скан 10 003 шага; ПОСЛЕ —
# 17 339 вызовов, 226 484 шага, самый длинный 2 057 (и 569 внутри десяти
# приёмок phicore). Бюджет взят 10^5 — в 48 раз выше замеренного
# максимума, цена ≤ 0,04 с при замеренных 2,4·10^6 шагов/с.
_SCAN_BUDGET = 10 ** 5


def _strip_zero_root(c):
    """p = x^m·g, где g(0) ≠ 0 — младшие нули сняты. Ненулевые корни p и
    g те же, а к g применима теорема о рациональных корнях."""
    i = 0
    while i < len(c) - 1 and c[i] == 0:
        i += 1
    return c[i:]


def _integer_root_by_theorem(p_sf: IntPoly, lo: Q, hi: Q):
    """Наименьший ЦЕЛЫЙ корень p_sf на [lo, hi] — теоремой о рациональных
    корнях (целый корень делит свободный член), а не сканом. Отказ
    поимённо, если перечень делителей не собирается."""
    ln, ld = lo.numerator, lo.denominator
    hn, hd = hi.numerator, hi.denominator
    g = _strip_zero_root(list(p_sf.c))
    if len(g) == 2:
        # степень 1 — прямая формула, перечень делителей не нужен:
        # единственный корень g₁x + g₀ есть −g₀/g₁, и он цел ровно когда
        # g₁ делит g₀ (эта ветка снимает дорогу from_root при q = 1)
        g0, g1 = g
        if g0 % g1 == 0:
            k = -(g0 // g1)
            if ln <= k * ld and k * hd <= hn and p_sf.eval(k) == 0:
                return True, Q(k)
        return True, None
    ds = _divisors(g[0])
    if ds is None:
        return False, None
    cands = sorted({s * dv for dv in ds for s in (1, -1)})
    for k in cands:
        if ln <= k * ld and k * hd <= hn and p_sf.eval(k) == 0:
            return True, Q(k)
    return True, None


def _rational_root_by_factoring(p_sf: IntPoly, lo: Q, hi: Q, steps=None):
    """Наименьший рациональный корень p_sf на [lo, hi] дорогой
    факторизации над ℤ (intfactor, Зассенхаус): рациональный корень p/q
    есть линейный множитель q·x − p, перечень делителей не нужен.

    Зовётся там, где перечень делителей не собрался (третья дорога,
    18.09). До неё в этом месте стояло два ответа: за бюджетом скана —
    отказ «сканов больше» при том, что устройство intfactor в доме есть
    и корень 10⁴ у x² − 10⁸ на [0; 10⁸ + 1] оно находит; ниже бюджета —
    `return None`, то есть ЛОЖНЫЙ КЛАСС: узел (√2 + 1/30030) − √2 отдавал
    величину класса 3 с as_fraction() = None при == 1/30030 → True
    (старший коэффициент аннулятора 27081081027000 = (2·3·5·7·11·13)³
    несёт 4096 делителей при пороге 64) `[вычислено здесь]`. Замер 18.09
    по 13 наборам, 249 вызовов _try_rational: на эту дорогу ушли бы 5,
    все — x² − 10⁸ из приёмки rational; остальные снимает монический
    короткий путь (187) либо перечень делителей.

    Отказ факторизации — отказ поимённо, не None: несведённое не
    выдаётся за иррациональное (обещание шапки модуля). Текст отказа
    несёт число сканов и бюджет, когда дорога пришла из-за бюджета
    (steps), и оба пути уточнения — refine и intfactor."""
    from .intfactor import factor_squarefree
    try:
        factors = factor_squarefree(p_sf)
    except Unresolved as e:
        budget = (f"сканов больше {_SCAN_BUDGET} (нужно {steps} целых на "
                  f"отрезке), " if steps is not None else "")
        raise Unresolved(
            "рациональный корень",
            f"Φ(C=3, P=[{','.join(map(str, p_sf.c))}], отрезок=[{lo}; {hi}])",
            hint=f"{budget}перечень делителей свободного члена ({p_sf.c[0]}) "
                 f"и старшего коэффициента ({p_sf.c[-1]}) не собрался, и "
                 f"факторизация над ℤ устройством intfactor отказала "
                 f"({e}) — несведённое не выдаётся за иррациональное. "
                 f"Путь уточнения: сузить отрезок изоляции (refine) либо "
                 f"спросить величину прямой формулой (exact_root для "
                 f"радикала)") from e
    roots = []
    for f in factors:
        if f.degree == 1:
            r = Q(-f.c[0], f.c[1])
            if lo <= r <= hi:
                roots.append(r)
    return min(roots) if roots else None


def _try_rational(p_sf: IntPoly, lo: Q, hi: Q):
    """Рациональный корень на [lo, hi] или None. Перебор целочисленный:
    сравнения — перекрёстным умножением, ноль — гомогенизацией;
    Q строится один раз, на найденном корне (граница).

    Целая дорога со СТОРОЖЕМ (12.09). Скан целых стоит ⌊hi⌋−⌈lo⌉+2 шагов,
    то есть растёт как САМА ВЕЛИЧИНА, а не как её разрядность: замерено
    2,4·10^6 шагов/с, отсюда Q(10)**Q(8) — 39 972 мс, а сумма над
    основанием 10 с размахом показателей 11 не досчитывалась вовсе.
    Решение по замеру — обе дороги, а не одна:

    * прямая формула там, где она есть, — вызывающие стороны
      (`from_root` через `exact_root`, `Q.__pow__` через целую степень);
      после них самый длинный скан на всём наборе соседей упал с 10 003
      шагов до 2 057, а общее число шагов с 331 021 до 226 484;
    * СТОРОЖ здесь, потому что общий многочлен из узлов +/× формулы не
      имеет: выше бюджета `_SCAN_BUDGET` = 10^5 скан заменяется теоремой
      о рациональных корнях, а если перечень делителей не собирается —
      отказ ПОИМЁННО с числом шагов, а не молчаливый счёт без конца.

    Ниже бюджета дорога побитово та же, что до 12.09 (тот же порядок
    перебора, тот же наименьший корень), потому прежние ответы целы.

    Третья дорога (18.09) — факторизация над ℤ там, где перечень
    делителей не собрался: до неё за бюджетом стоял отказ при живом
    устройстве, а ниже бюджета — `return None`, ложный класс (три случая
    замерены: неполный перечень «атомом», перечень за порогом 64 и
    свободный член 0 при нуле вне отрезка — x(30030x − 1) на
    [1/60060; 1] отдавал None при == 1/30030 → True). Перед перечнем
    снимается множитель x (теорема о рациональных корнях говорит о
    g(0) ≠ 0). Монический аннулятор (старший коэффициент ±1) перебора
    дробей не требует: его рациональные корни целы, а целые на отрезке
    уже исчерпаны сканом либо теоремой — ответ None без перечня."""
    ln, ld = lo.numerator, lo.denominator
    hn, hd = hi.numerator, hi.denominator
    if p_sf.c[0] == 0 and ln <= 0 <= hn:
        return Q(0)
    k = -((-ln) // ld)          # ⌈lo⌉
    kmax = hn // hd             # ⌊hi⌋
    steps = kmax + 2 - k
    if steps > _SCAN_BUDGET:
        enumerated, root = _integer_root_by_theorem(p_sf, lo, hi)
        if not enumerated:
            # перечень делителей свободного члена не собрался — дорога
            # факторизации (до 18.09 здесь стоял отказ «сканов больше»)
            return _rational_root_by_factoring(p_sf, lo, hi, steps)
        if root is not None:
            return root
    else:
        while k <= kmax + 1:
            if ln <= k * ld and k * hd <= hn and p_sf.eval(k) == 0:
                return Q(k)
            k += 1
    g = _strip_zero_root(list(p_sf.c))      # ноль вне отрезка: корни g те же
    a0, an = g[0], g[-1]
    if abs(an) == 1:
        return None             # монический: рациональный корень цел, целые исчерпаны
    ps, qs = _divisors(a0), _divisors(an)
    if ps is None or qs is None:
        # до 18.09: return None — ложный класс при неполном перечне
        return _rational_root_by_factoring(p_sf, lo, hi)
    for q in qs:
        for pnum in ps:
            for num in (pnum, -pnum):
                if (ln * q <= num * ld and num * hd <= hn * q
                        and _homog_sign(p_sf.c, num, q) == 0):
                    return Q(num, q)
    return None


def _power_by_theorem(base, exponent: "AlgebraicNumber"):
    """base^θ при показателе класса 3 (дверь Гильберта 7, 04.09).

    Показатель, рациональный на линейке, уходит в обычную степень (2^θ при
    θ = √4 = 2 есть 4). Основания 1 и 0 сведены без теоремы: 1^θ = 1,
    0^θ = 0 при θ > 0. Прочее — Гельфонд–Шнайдер (1934): α ∉ {0, 1}
    алгебраическое, θ алгебраическое иррациональное ⇒ α^θ трансцендентно.
    Носителя для такого значения на линейках 1–3 нет, и класс записи здесь
    назначает теорема, а не замер: отказ поимённо, не TypeError (до 04.09
    2^√2 падал как «Q не строится из AlgebraicNumber»)."""
    r = exponent.as_fraction()
    if r is not None:
        return base ** r
    b = Q(base) if isinstance(base, int) else base
    if isinstance(b, AlgebraicNumber):
        rb = b.as_fraction()
        if rb is not None:
            b = rb
    if isinstance(b, Q):
        if b == 1:
            return Q(1)
        if b == 0:
            if exponent._separate_sign() > 0:
                return Q(0)
            raise ZeroDivisionError("0 в неположительной степени")
        brec = f"Φ(C=1 RATIONAL, value={b})"
    else:
        brec = b.record()
    raise Unresolved("^", brec, exponent.record(),
                     hint="Гельфонд–Шнайдер (1934, 7-я проблема Гильберта): "
                          "алгебраическое основание вне {0, 1} в алгебраической "
                          "иррациональной степени трансцендентно — класс записи "
                          "назначен теоремой, носителя на линейках 1–3 нет")


def _defers_to(other, division: bool = False, product: bool = False) -> bool:
    """Операнд, чья дверь сводит действие с линейкой 3 сама — отражённой
    операцией Python: пара (re, im), класс 19, корень единицы, класс 4,
    и — в × и / — моном π, класс 10 (коэффициент монома несёт класс 3).
    Замер 25.09 (сдача Γ ↔ π, проверка G) и 01.10 `[вычислено здесь]`:
    √2 × 3π² отказывал, 3π² × √2 давал моном с коэффициентом класса 3 —
    асимметрия порядка в коммутативном действии; снята уступкой двери
    монома (`PiMonomial.__rmul__`, `__rtruediv__`). В + и − отражённой
    двери для линейки 3 у монома нет (`__radd__` → NotImplemented), и
    уступка дала бы голый TypeError — потому сложение остаётся отказом
    линейки 3 по имени, обеими сторонами (очередь, строка 103).
    До 18.09 линейка отказывала первой: AN + пара, AN − пара, AN × пара,
    AN / пара отказывали 4 из 4, тогда как пара + AN и т. д. сводились
    4 из 4 `[вычислено здесь]`. Корень единицы с линейкой 3 в +, − и ×
    отказывает своей дверью (unity_root.py: «сумма корней единицы вне
    линейки», «правила сведения с этой линейкой нет») — обе стороны
    говорят одним текстом; в делении у корня единицы отражённой двери
    нет (__rtruediv__ → NotImplemented), уступка дала бы голый TypeError
    (замер 18.09), потому AN / ζ остаётся отказом линейки 3."""
    from .complex_pair import ComplexPair
    if isinstance(other, ComplexPair):
        return True
    from .pi_value import PiMonomial
    if isinstance(other, PiMonomial):
        return product or division
    if division:
        return False
    from .unity_root import UnityRoot
    return isinstance(other, UnityRoot)


# ---------- сама линейка ----------

class AlgebraicNumber:
    __slots__ = ("p", "lo", "hi")

    def __init__(self, p: IntPoly, lo: Q, hi: Q):
        self.p = squarefree_part(p)
        self.lo, self.hi = Q(lo), Q(hi)
        if count_roots(self.p, self.lo, self.hi) != 1:
            raise Unresolved("изоляция", self.record(),
                             hint="отрезок не изолирует ровно один корень")

    # --- построение ---
    @classmethod
    def from_root(cls, q: int, radicand: Q, sign: int = 0) -> "AlgebraicNumber":
        """√[q]{radicand} со знаком: корень den·x^q − num.

        Рациональность корня решается ПРЯМОЙ ФОРМУЛОЙ (`exact_root`: обе
        части канонической пары — полные q-е степени), а не поиском по
        изолирующему отрезку: отрезок здесь шириной самой величины, и
        скан по нему стоил 49 076 мс на Q(10**8)**Q(1,3) и не
        досчитывался на Q(10**20)**Q(1,2). Отрезок изоляции [0; hi]
        оставлен прежним — он входит в запись record().

        Подкоренное любого знака (18.09; до того дверь отказывала на 4
        входах из 4, где величина в доме есть, и подсказка звала нечётный
        корень «чётным» `[вычислено здесь]`): 0 → Q(0); отрицательное при
        нечётном q — вещественный корень −|r|^(1/q) на этой же линейке
        (так уже шла дверь Q.__pow__ в rational.py:210); при чётном q —
        пара (re, im) главного значения (ComplexPair.from_real_root), как
        у `_root` этой линейки и у Q.__pow__ (rational.py:209). sign=1
        обращает результат при любом знаке подкоренного."""
        radicand = Q(radicand)
        if radicand == 0:
            return Q(0)
        if radicand < 0:
            if q % 2 == 0:
                from .complex_pair import ComplexPair
                val = ComplexPair.from_real_root(radicand, q)
            else:
                val = -cls.from_root(q, -radicand)
            return -val if sign == 1 else val
        num, den = radicand.numerator, radicand.denominator
        r = exact_root(radicand, q)          # прямая формула, не скан
        if r is not None:
            return -r if sign == 1 else r
        poly = IntPoly([-num] + [0] * (q - 1) + [den])
        hi = max(Q(1), radicand) + 1
        val = cls(poly, Q(0), hi)
        return -val if sign == 1 else val

    # --- запись ---
    def record(self) -> str:
        return (f"Φ(C=3 POLY_ROOT, P=[{','.join(map(str, self.p.c))}], "
                f"отрезок=[{self.lo}; {self.hi}])")

    __repr__ = record

    def annihilator(self, reduce: bool = False) -> IntPoly:
        """Аннулятор корня: свободный от квадратов многочлен из ℤ[x], обращающийся
        в ноль на величине. НЕ обязательно минимальный: узлы +/× строят его
        результантом, и он может быть приводим (замер 04.09 на семействе √2,
        √3, √5, ∛2: приводимых 2 из 12 — (√2+√3)² даёт степень 6 при
        минимальной 2, (1+√2)·√2 — 4 при 2). reduce=True снимает рациональные
        линейные множители, не обращающиеся в ноль на отрезке изоляции
        (дёшево: перебор делителей); полное решение — minimal_polynomial()
        (факторизация над ℤ, intfactor.py, с 04.09)."""
        if not reduce:
            return self.p
        c = list(self.p.c)
        changed = True
        while changed and len(c) > 2:
            changed = False
            b = 1 + max(abs(x) for x in c[:-1]) // abs(c[-1])       # граница Коши для корней
            cand = _try_rational(IntPoly(c), Q(-b), Q(b))
            roots = [cand] if cand is not None else []
            for r in roots:
                r = Q(r)
                if self.lo <= r <= self.hi:
                    continue                                    # этот корень — наш отрезок; множитель оставить
                q = _divexact(c, [-r.numerator, r.denominator])   # (den·x − num) делит нацело
                if _norm(list(q)) != [0]:
                    c = _prim(_norm(q)); changed = True; break
        return IntPoly(c)

    def minimal_annihilator(self) -> IntPoly:
        """Историческое имя; обещает больше, чем даёт: см. annihilator() —
        минимальность не гарантирована (замер 04.09); минимальный многочлен
        по имени — minimal_polynomial()."""
        return self.p

    def minimal_polynomial(self) -> IntPoly:
        """Минимальный многочлен величины над ℚ (примитивный, lc > 0):
        неприводимый множитель аннулятора, имеющий корень на отрезке
        изоляции. Аннулятор свободен от квадратов и отрезок изолирует
        ровно один корень, поэтому такой множитель ровно один. Дорога —
        Зассенхаус над ℤ (intfactor.py); при числе модульных множителей
        за порогом — отказ Unresolved, не молчание."""
        from .intfactor import factor_squarefree
        factors = factor_squarefree(self.p)
        hits = [f for f in factors if count_roots(f, self.lo, self.hi) == 1]
        if len(hits) != 1:
            raise Unresolved("минимальный многочлен", self.record(),
                             hint=f"множителей с корнем на отрезке: {len(hits)}, ожидался 1")
        return hits[0]

    def is_minimal(self) -> bool:
        """Совпадает ли аннулятор с минимальным многочленом. Степень ≤ 3
        решается перебором рациональных корней (неприводима ⇔ нет
        рационального корня; степень 1 — всегда); степень ≥ 4 — факторизацией
        над ℤ (с 04.09; раньше здесь стояло None)."""
        deg = len(self.p.c) - 1
        if deg <= 1:
            return True
        if deg <= 3:
            c = self.p.c
            b = 1 + max(abs(x) for x in c[:-1]) // abs(c[-1])       # граница Коши
            return _try_rational(self.p, Q(-b), Q(b)) is None
        return self.minimal_polynomial().degree == deg

    def self_check(self, deep: bool = False) -> list:
        """Линейка спрашивает саму себя («AIAΣA спрашивает AIAΣA», 04.09):
        (1) аннулятор свободен от квадратов; (2) отрезок изолирует ровно
        один корень (Штурм); (3) знак согласован: P(lo)·P(hi) ≤ 0 либо
        P(hi) = 0; deep=True — (4) минимальный многочлен делит аннулятор и
        сам имеет ровно один корень на отрезке (факторизация над ℤ).
        Возвращает список находок; пусто — согласовано."""
        out = []
        if len(_int_poly_gcd(self.p.c, self.p.derivative().c)) != 1:
            out.append("аннулятор не свободен от квадратов")
        n = count_roots(self.p, self.lo, self.hi)
        if n != 1:
            out.append(f"на отрезке корней {n}, ожидался 1")
        s_lo = _homog_sign(self.p.c, self.lo.numerator, self.lo.denominator)
        s_hi = _homog_sign(self.p.c, self.hi.numerator, self.hi.denominator)
        if s_hi != 0 and s_lo * s_hi > 0:
            out.append("знак аннулятора на концах отрезка не меняется")
        if deep:
            try:
                m = self.minimal_polynomial()
            except Unresolved as e:
                out.append(f"минимальный многочлен не вычислен: {e}")
            else:
                if _norm(list(_prem(self.p.c, m.c)[0])) != [0]:
                    out.append("минимальный многочлен не делит аннулятор")
                if count_roots(m, self.lo, self.hi) != 1:
                    out.append("минимальный многочлен без единственного корня на отрезке")
        return out

    def mirror_sum(self):
        """Зазеркальная сумма θ + 1/θ (слово автора 04.09: палиндромы и
        симметрии — инструмент). Если аннулятор — палиндром чётной степени
        2h, сумма есть корень Q = palindromic_reduce степени h: степень
        делится пополам без результанта (√2+√3 ↔ x⁴−10x²+1 → y² − 12:
        (√2+√3) + 1/(√2+√3) = 2√3). Не палиндром — обычный узел θ + 1/θ."""
        p = self.p
        if p.degree % 2 == 0 and p.is_palindromic():
            q = p.palindromic_reduce()
            inv = self.inverse()
            lo, hi = self.lo + inv.lo, self.hi + inv.hi
            a, b = AlgebraicNumber(p, self.lo, self.hi), inv
            def refiners():
                a.refine(); b.refine()
                return a.lo + b.lo, a.hi + b.hi
            return AlgebraicNumber._from_annihilator(q, lo, hi, refiners)
        return self + self.inverse()

    # --- уточнение отрезка делением пополам по знаку ---
    def enclosure(self, terms: int = 30):
        """Сертификат сетки — тот же интерфейс, что у прочих носителей:
        отрезок изоляции после terms уточнений (каждое — деление пополам
        по Штурму); величина внутри точно, ширина ≤ исходной/2^terms."""
        from .enclosure import Interval
        self.refine(terms)
        return Interval(self.lo, self.hi)

    def refine(self, times: int = 1) -> None:
        """Деление пополам по знаку. Сравнение — с правым концом: корень
        лежит в (lo, hi], потому P(hi) либо ≠ 0, либо hi и есть корень;
        левый конец сравнением не годится — lo может быть соседним корнем
        аннулятора (P(lo) = 0), и знак там пуст (дефект найден 04.09 на
        ζ₆ − ζ₆⁻¹: отрезок [0; 4] у x³−3x не сужался, equals отказывал)."""
        for _ in range(times):
            mid = (self.lo + self.hi) / 2
            s = _homog_sign(self.p.c, mid.numerator, mid.denominator)
            if s == 0:
                self.lo = self.hi = mid
                return
            s_hi = _homog_sign(self.p.c, self.hi.numerator, self.hi.denominator)
            if s_hi == 0 or s != s_hi:      # корень в (mid, hi]
                self.lo = mid
            else:                           # корень в (lo, mid]
                self.hi = mid

    def as_fraction(self):
        """Q, если величина рациональна на этой линейке; иначе None."""
        return _try_rational(self.p, self.lo, self.hi)

    # --- сборка результата бинарного действия ---
    @classmethod
    def _from_annihilator(cls, raw: IntPoly, lo: Q, hi: Q,
                          refiners) -> "AlgebraicNumber | Q":
        sf = squarefree_part(raw)
        for _ in range(_MAX_REFINE):
            if lo == hi:
                # вырожденный отрезок: уточнение попало ровно в корень
                return Q(lo)
            n = count_roots(sf, lo, hi)
            if n == 1:
                # единственный корень на (lo, hi]; рациональный кандидат
                # принимается только внутри этого полуинтервала — корень
                # на левой границе принадлежит другому значению
                r = _try_rational(sf, lo, hi)
                if r is not None and lo < r <= hi:
                    return r
                return cls(sf, lo, hi)
            lo, hi = refiners()
        raise Unresolved("изоляция", f"Φ(C=3, P=[{','.join(map(str, sf.c))}])",
                         hint=f"корень не изолирован за {_MAX_REFINE} делений")

    # --- арифметика ---
    def __neg__(self) -> "AlgebraicNumber":
        mirrored = IntPoly([a if i % 2 == 0 else -a for i, a in enumerate(self.p.c)])
        return AlgebraicNumber(mirrored, -self.hi, -self.lo)

    def _shift(self, r: Q) -> "AlgebraicNumber":
        """Величина + r: целые коэффициенты rd^d·P(x−r), дробей в цикле нет."""
        rn, rd = r.numerator, r.denominator
        d = len(self.p.c) - 1
        out = [0] * (d + 1)
        for j, a in enumerate(self.p.c):
            if a == 0:
                continue
            # rd^d·P(x−r) = Σ_j a_j·rd^(d−j)·(rd·x − rn)^j; вклад в x^(j−k):
            # a_j·C(j,k)·rd^(d−k)·(−rn)^k
            binom = 1
            for k in range(j + 1):
                out[j - k] += a * binom * rd ** (d - k) * (-rn) ** k
                binom = binom * (j - k) // (k + 1)
        return AlgebraicNumber(IntPoly(out), self.lo + r, self.hi + r)

    def _scale(self, r: Q):
        """Величина · r (r — рациональное): корни умножаются на r,
        коэффициенты P(x/r)·r^d очищены в целые степенями rn и rd."""
        if r == 0:
            return Q(0)
        d = self.p.degree
        rn, rd = r.numerator, r.denominator
        out = [self.p.c[i] * rn ** (d - i) * rd ** i for i in range(d + 1)]
        poly = IntPoly(out)
        a, b = self.lo * r, self.hi * r
        return AlgebraicNumber(poly, min(a, b), max(a, b))

    def __add__(self, other):
        from .zero_state import Zero
        if isinstance(other, Zero):
            return NotImplemented
        if isinstance(other, int):
            other = Q(other)
        if isinstance(other, Q):
            return self._shift(other)
        if _defers_to(other):
            return NotImplemented
        if not isinstance(other, AlgebraicNumber):
            rec = other.record() if hasattr(other, "record") else repr(other)
            raise Unresolved("+", self.record(), rec)
        a, b = self, other
        da, db = a.p.degree, b.p.degree
        pts = []
        fa = a.p.c  # целые коэффициенты — результант считается на ℤ
        for t in range(da * db + 1):
            qy = b.p.compose_x_minus_y_coeff_in_y(t)
            pts.append((t, resultant_values(fa, qy)))
        raw = interpolate_int_poly(pts)

        def refiners():
            a.refine()
            b.refine()
            return a.lo + b.lo, a.hi + b.hi

        return AlgebraicNumber._from_annihilator(raw, a.lo + b.lo, a.hi + b.hi,
                                                 refiners)

    __radd__ = __add__

    def __sub__(self, other):
        if isinstance(other, (int, Q)):
            return self._shift(Q(-other))
        return self + (-other)

    def __rsub__(self, other):
        return (-self) + other

    def __mul__(self, other):
        from .zero_state import Zero
        if isinstance(other, Zero):
            return NotImplemented
        if isinstance(other, int):
            other = Q(other)
        if isinstance(other, Q):
            return self._scale(other)
        if _defers_to(other, product=True):
            return NotImplemented
        if not isinstance(other, AlgebraicNumber):
            rec = other.record() if hasattr(other, "record") else repr(other)
            raise Unresolved("×", self.record(), rec)
        a, b = self, other
        da, db = a.p.degree, b.p.degree
        fa = a.p.c  # целые коэффициенты — результант считается на ℤ
        pts = []
        for t in range(da * db + 1):
            # q_t(y) = y^db · P_b(t/y) = Σ b_j · t^j · y^(db−j) — целые
            qy = [0] * (db + 1)
            for j, bj in enumerate(b.p.c):
                qy[db - j] += bj * t ** j
            pts.append((t, resultant_values(fa, qy)))
        raw = interpolate_int_poly(pts)

        def bounds():
            cands = [a.lo * b.lo, a.lo * b.hi, a.hi * b.lo, a.hi * b.hi]
            return min(cands), max(cands)

        def refiners():
            a.refine()
            b.refine()
            return bounds()

        lo, hi = bounds()
        return AlgebraicNumber._from_annihilator(raw, lo, hi, refiners)

    __rmul__ = __mul__

    def eval_poly(self, coeffs):
        """p(θ) для рационального многочлена p — ОДНИМ результантом:
        Res_x(P(x), L·y − p̃(x)); отрезок — интервальным Хорнером по
        отрезку θ. Степень не растёт цепочкой умножений (это и есть
        путь для сумм над радикальным полем ℚ(α^{1/N}))."""
        from math import lcm
        cs = [Q(c) for c in coeffs]
        while cs and cs[-1] == 0:
            cs.pop()
        if not cs:
            return Q(0)
        if len(cs) == 1:
            return cs[0]
        L = 1
        for c in cs:
            L = lcm(L, c.denominator)
        pt = [c.numerator * (L // c.denominator) for c in cs]   # L·p, целые
        fa = self.p.c
        dP = self.p.degree
        pts = []
        for t in range(dP + 1):
            # L·t − p̃(x): коэффициенты млад→ст
            qx = [L * t - pt[0]] + [-a for a in pt[1:]]
            pts.append((t, resultant_values(fa, qx)))
        raw = interpolate_int_poly(pts)
        theta = self

        def bounds():
            lo, hi = theta.lo, theta.hi
            a, b = cs[-1], cs[-1]
            for c in reversed(cs[:-1]):
                cands = (a * lo, a * hi, b * lo, b * hi)
                a, b = min(cands) + c, max(cands) + c
            return a, b

        def refiners():
            theta.refine()
            return bounds()

        lo, hi = bounds()
        return AlgebraicNumber._from_annihilator(raw, lo, hi, refiners)

    def inverse(self) -> "AlgebraicNumber":
        # отрезок, несущий 0 на краю или внутри, сперва отделяется от нуля
        # (дефект найден свипом ФIOS: from_root даёт [0; hi])
        if self.p.c[0] == 0 and self.lo <= 0 <= self.hi:
            raise Unresolved("1/x", self.record(), hint="ноль не обратим")
        for _ in range(_MAX_REFINE):
            if not (self.lo <= 0 <= self.hi):
                break
            self.refine()
        else:
            raise Unresolved("1/x", self.record(),
                             hint="отрезок не отделился от нуля")
        rev = IntPoly(list(reversed(self.p.c)))
        a, b = 1 / self.hi, 1 / self.lo
        return AlgebraicNumber(rev, min(a, b), max(a, b))

    def __truediv__(self, other):
        if isinstance(other, (int, Q)):
            return self._scale(1 / Q(other))
        if isinstance(other, AlgebraicNumber):
            return self * other.inverse()
        if _defers_to(other, division=True):
            return NotImplemented
        rec = other.record() if hasattr(other, "record") else repr(other)
        raise Unresolved("/", self.record(), rec)

    # --- степень: сущность Ф-числа — дробная экспонента (слово автора
    # 02.09: p/q возводим в степень, q/p извлекаем корень; отказ не
    # одобрен) ---
    def _separate_sign(self) -> int:
        """Знак величины; отрезок отделяется от нуля уточнением."""
        for _ in range(_MAX_REFINE):
            if not (self.lo <= 0 <= self.hi):
                break
            self.refine()
        else:
            raise Unresolved("знак", self.record(),
                             hint="отрезок не отделился от нуля")
        return 1 if self.lo >= 0 else -1

    def sign_certificate(self) -> dict:
        """Знак двумя независимыми дорогами (Гильберт 17 в дискретной форме,
        04.09): «отрезок» — деление пополам до отделения от нуля; «Штурм» —
        счёт корней аннулятора на (0; hi] по ИСХОДНОМУ отрезку, без единого
        деления. Дороги не знают друг о друге; расхождение — отказ, а не
        выбор одной из них. Сертификат суммой квадратов (Артин 1927) не
        строится: он требует SDP с точным округлением, спрос корпуса — 0."""
        lo, hi = self.lo, self.hi
        zero_inside = self.p.c[0] == 0 and lo < 0 <= hi
        if zero_inside:
            sturm = 0
        elif hi <= 0:
            sturm = -1
        elif lo >= 0:
            sturm = 1
        else:
            sturm = 1 if count_roots(self.p, Q(0), hi) == 1 else -1
        interval = 0 if zero_inside else self._separate_sign()
        if interval != sturm:
            raise Unresolved("знак", self.record(),
                             hint=f"дороги расходятся: отрезок {interval}, "
                                  f"Штурм {sturm}")
        return {"отрезок": interval, "Штурм": sturm, "знак": sturm}

    def sign(self) -> int:
        """Знак величины по сертификату двух дорог."""
        return self.sign_certificate()["знак"]

    def _root(self, q: int):
        """q-й корень величины: подстановка x→x^q в аннулятор — корень
        живёт на той же линейке точно; нечётный q проходит знак."""
        r = self.as_fraction()
        if r is not None:
            return r ** Q(1, q)
        if self._separate_sign() < 0:
            if q % 2 == 0:
                # вещественного значения нет — пара (re, im): |θ|^(1/q) под
                # углом π/q (главное значение), cos/sin — класс 3 (04.09)
                from .complex_pair import ComplexPair
                return ComplexPair.from_angle((-self)._root(q), Q(1, 2 * q))
            return -((-self)._root(q))
        w = AlgebraicNumber(self.p, self.lo, self.hi)
        for _ in range(_MAX_REFINE):
            if w.lo >= 0:
                break
            w.refine()
        c2 = [0] * (q * w.p.degree + 1)
        for i, a in enumerate(w.p.c):
            c2[q * i] = a
        raw = IntPoly(c2)
        prec = [8]

        def bounds():
            g = prec[0]
            n, d = w.lo.numerator, w.lo.denominator
            zl = Q(_int_root_floor((n * d ** (q - 1)) << (q * g), q),
                   d << g)
            n, d = w.hi.numerator, w.hi.denominator
            zh = Q(_int_root_floor((n * d ** (q - 1)) << (q * g), q) + 1,
                   d << g)
            return zl, zh

        def refiners():
            w.refine()
            prec[0] += 8
            return bounds()

        lo, hi = bounds()
        return AlgebraicNumber._from_annihilator(raw, lo, hi, refiners)

    def __pow__(self, p):
        """Целый показатель — квадрат-и-умножение на резльтантах;
        пара p/q — степень |p| и q-й корень, всё точно на линейке."""
        if isinstance(p, bool):
            return NotImplemented
        if isinstance(p, float):
            raise TypeError("float-показатель запрещён построением")
        if isinstance(p, AlgebraicNumber):
            return _power_by_theorem(self, p)
        if isinstance(p, int):
            if p == 0:
                return Q(1)
            base = self.inverse() if p < 0 else self
            p = abs(p)
            if p == 1:
                return base
            # θ^p = значение многочлена x^p в θ — один результант в том же
            # поле, а не цепочка произведений (степень d, не d^p)
            return base.eval_poly([Q(0)] * p + [Q(1)])
        if not isinstance(p, Q):
            return NotImplemented
        base = self.inverse() if p.n < 0 else self
        w = base ** abs(p.n)
        if p.d == 1:
            return w
        if isinstance(w, Q):
            return w ** Q(1, p.d)
        return w._root(p.d)

    def __rpow__(self, base):
        """основание ** θ: int/Q в степени класса 3 — по теореме (см.
        _power_by_theorem); float — запрещён построением."""
        if isinstance(base, bool) or isinstance(base, float):
            raise TypeError("основание — только int или Q")
        if isinstance(base, (int, Q)):
            return _power_by_theorem(base, self)
        return NotImplemented

    # --- сравнение величин ---
    def equals(self, other) -> bool:
        """Тождество величин: разность распознаётся как точный ноль."""
        diff = self + (-other) if isinstance(other, AlgebraicNumber) \
            else self - Q(other)
        return isinstance(diff, Q) and diff == 0

    def __eq__(self, other):
        if isinstance(other, AlgebraicNumber):
            if self.p == other.p and self.lo == other.lo and self.hi == other.hi:
                return True
            return self.equals(other)
        if isinstance(other, (int, Q)):
            # точно, не best-effort: величина — единственный корень
            # аннулятора на (lo; hi], поэтому «= r» ⟺ P(r)=0 и r внутри
            r = Q(other)
            if self.lo == self.hi:
                return Q(self.lo) == r
            if not (self.lo < r <= self.hi):
                return False
            return _homog_sign(self.p.c, r.numerator, r.denominator) == 0
        return NotImplemented

    def __hash__(self):
        return hash((tuple(self.p.c),))
