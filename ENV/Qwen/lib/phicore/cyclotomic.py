"""Круговые величины на линейке класса 3: 2cos(2πk/n), cos(π·r), sin(π·r)
при любом рациональном аргументе — алгебраические числа степени φ(m)/2.

До 04.09 sin_pi знал таблицу знаменателей {1, 2, 3, 4, 6}, вне её — отказ
(пересмотр отказов 23221be3: «следующий шаг — минимальный многочлен
2cos(2πk/n)»). Дорога — из строения делителей, без факторизации:
круговой многочлен Φ_m(y) = (y^m − 1) / ∏_{d|m, d<m} Φ_d(y) (точное деление
над ℤ), он палиндромен, и Φ_m(y) = y^{φ(m)/2}·Ψ_m(y + 1/y), где Ψ_m —
минимальный многочлен 2cos(2π/m) [внешний источник: Lehmer 1933;
Watkins, Zeitlin 1993]; переход y^k + y^{−k} = V_k(x) — многочлены Люка
V₀ = 2, V₁ = x, V_k = x·V_{k−1} − V_{k−2} (2cos kθ = V_k(2cos θ)).
Корни Ψ_m — 2cos(2πi/m), gcd(i, m) = 1, i ≤ m/2, упорядочены по i;
нужный — по рангу, изоляция делением пополам по Штурму: без π, без float.

Первая версия (04.09) брала свободную от квадратов часть V_n(x) − 2 и
факторизовала над ℤ: на n = 360 (степень 181) модульных множителей 35 —
за порогом перебора Зассенхауса; дорога через Φ_m порога не имеет.
Факторизация остаётся независимой сверкой в приёмке.
"""
from __future__ import annotations

from functools import lru_cache
from math import gcd

from .algebraic import AlgebraicNumber, count_roots, _divexact
from .intpoly import IntPoly
from .rational import Q

_NUM = (int, Q)


def lucas_poly(n: int) -> IntPoly:
    """V_n: 2cos(nθ) = V_n(2cosθ); V₀ = 2, V₁ = x, V_n = x·V_{n−1} − V_{n−2}."""
    if n < 0:
        raise ValueError("n ≥ 0")
    a, b = IntPoly([2]), IntPoly([0, 1])
    if n == 0:
        return a
    x = IntPoly([0, 1])
    for _ in range(n - 1):
        a, b = b, x * b - a
    return b


def _divisors(m: int) -> list:
    return [d for d in range(1, m + 1) if m % d == 0]


@lru_cache(maxsize=4096)
def cyclotomic_poly(m: int) -> IntPoly:
    """Φ_m(y) точным делением: (y^m − 1) / ∏_{d|m, d<m} Φ_d(y)."""
    if m < 1:
        raise ValueError("m ≥ 1")
    num = [-1] + [0] * (m - 1) + [1]
    for d in _divisors(m):
        if d < m:
            num = _divexact(num, cyclotomic_poly(d).c)
    return IntPoly(num)


@lru_cache(maxsize=4096)
def psi(m: int) -> IntPoly:
    """Ψ_m — минимальный многочлен 2cos(2π/m): Ψ₁ = x − 2, Ψ₂ = x + 2,
    при m ≥ 3 из палиндрома Φ_m: Φ_m(y)/y^h = a_h + Σ_{k=1}^{h} a_{h+k}(y^k + y^{−k}),
    h = φ(m)/2, и y^k + y^{−k} = V_k(x)."""
    if m == 1:
        return IntPoly([-2, 1])
    if m == 2:
        return IntPoly([2, 1])
    a = cyclotomic_poly(m).c
    h = (len(a) - 1) // 2
    out = IntPoly([a[h]])
    for k in range(1, h + 1):
        out = out + lucas_poly(k).scale(a[h + k])
    check_psi(m, out)
    return out


def check_psi(m: int, p: IntPoly) -> None:
    """Самопроверка Ψ_m («AIAΣA спрашивает AIAΣA», 04.09): степень равна
    φ(m)/2 (половина степени Φ_m), многочлен унитарен и палиндром Φ_m
    восстанавливается: Φ_m(y) = y^h·Ψ_m(y + 1/y) — проверяется значением
    в y = 2: Φ_m(2) = 2^h·Ψ_m(5/2). Расхождение — отказ, не тихий ответ."""
    from .errors import Unresolved
    phi = cyclotomic_poly(m)
    h = phi.degree // 2
    if p.degree != h or p.c[-1] != 1:
        raise Unresolved("Ψ_m", f"m={m}", hint=f"самопроверка: степень {p.degree} ≠ φ(m)/2 = {h} или не унитарен")
    if phi.eval(2) != 2 ** h * p.eval(Q(5, 2)):
        raise Unresolved("Ψ_m", f"m={m}", hint="самопроверка: Φ_m(2) ≠ 2^h·Ψ_m(5/2)")


def self_check(deep: bool = False) -> list:
    """Самопроверка круговой линейки ВСТРЕЧНОЙ дорогой (ось самоприменения).

    `[вычислено здесь]` Φ_m строится ДЕЛЕНИЕМ (y^m − 1 на произведение
    младших Φ_d), поэтому встречная дорога — УМНОЖЕНИЕ обратно: тождество
    ∏_{d | m} Φ_d(y) = y^m − 1 собирается здесь перемножением многочленов и
    сверяется с y^m − 1 коэффициент в коэффициент; ни одного деления на этой
    дороге нет. Вторая дорога — степень: deg Φ_m обязана равняться функции
    Эйлера φ(m), посчитанной здесь прямым счётом взаимно простых остатков, а
    не через делители. Третья — палиндром: `check_psi` восстанавливает Φ_m
    из Ψ_m значением в точке (Φ_m(2) = 2^h·Ψ_m(5/2)), и молчание там есть
    условие, а находка — отказ, пойманный поимённо.

    Граница объявлена дважды. Первая — палиндромная дорога идёт с m ≥ 3:
    Ψ₁ = x − 2 и Ψ₂ = x + 2 заданы явно, у них φ(m)/2 = 0, и `check_psi` к ним
    неприменим по построению (`psi` их и не зовёт) — это условие линейки, не
    находка (замерено этой же самопроверкой 12.09). Вторая — числом: перебор
    идёт по m ≤ 24 (при `deep` — до 60);
    выше цена растёт как φ(m)·число делителей, и молчание на непроверенных m
    находкой не считается — самопроверка говорит о проверенном окне, а не о
    всех m. Возвращает список находок; исключений не бросает.
    """
    out = []
    top = 60 if deep else 24
    for m in range(1, top + 1):
        try:
            phi_poly = cyclotomic_poly(m)
        except Exception as e:                                  # noqa: BLE001
            out.append(f"Φ_{m} не построен: {type(e).__name__}: {e}")
            continue
        prod = IntPoly([1])
        for d in _divisors(m):
            prod = prod * cyclotomic_poly(d)
        target = IntPoly([-1] + [0] * (m - 1) + [1])
        if prod != target:
            out.append(f"∏ Φ_d(y) по делителям {m} не равно y^{m} − 1: "
                       f"{prod!r}")
        totient = sum(1 for k in range(1, m + 1) if gcd(k, m) == 1)
        if phi_poly.degree != totient:
            out.append(f"deg Φ_{m} = {phi_poly.degree} ≠ φ({m}) = {totient}")
        if phi_poly.c[-1] != 1:
            out.append(f"Φ_{m} не унитарен: старший {phi_poly.c[-1]}")
        if m >= 3:            # Ψ₁ = x − 2 и Ψ₂ = x + 2 заданы явно: φ(m)/2 = 0
            try:
                check_psi(m, psi(m))
            except Exception as e:                              # noqa: BLE001
                out.append(f"Ψ_{m} не сошёлся со своим палиндромом: "
                           f"{type(e).__name__}: {e}")
    return out


def _isolate_rank(p: IntPoly, rank: int, lo: Q, hi: Q):
    """Отрезок (lo, hi] ровно с rank-м сверху корнем p (rank = 0 —
    наибольший); все корни p лежат в (lo, hi]. Штурм, деление пополам."""
    top = hi
    while count_roots(p, lo, hi) != 1:
        mid = (lo + hi) / 2
        if count_roots(p, mid, top) > rank:      # выше mid больше rank корней ⇒ x_rank > mid
            lo = mid
        else:
            hi = mid
    return lo, hi


def two_cos(k: int, n: int):
    """2cos(2πk/n) точно: Q при m = n/gcd(k, n) ∈ {1, 2, 3, 4, 6}, иначе
    AlgebraicNumber с минимальным многочленом Ψ_m степени φ(m)/2."""
    if isinstance(k, bool) or isinstance(n, bool) or not isinstance(k, int) or not isinstance(n, int):
        raise TypeError("k, n — только int")
    if n < 1:
        raise ValueError("n ≥ 1")
    k %= n
    j = min(k, n - k)                       # cos чётен: 2cos(2πk/n) = 2cos(2πj/n), 0 ≤ j ≤ n/2
    d = gcd(j, n)
    m, i = n // d, j // d                   # 2cos(2πi/m), gcd(i, m) = 1
    if m == 1:
        return Q(2)
    if m == 2:
        return Q(-2)
    rank = sum(1 for t in range(1, i) if gcd(t, m) == 1)   # корни Ψ_m убывают с ростом i
    p = psi(m)
    lo, hi = _isolate_rank(p, rank, Q(-3), Q(2))
    val = AlgebraicNumber(p, lo, hi)
    r = val.as_fraction()
    return r if r is not None else val


def cos_pi(r):
    """cos(π·r), r ∈ ℚ: = 2cos(2π·p/(2q))/2."""
    if isinstance(r, float) or not isinstance(r, _NUM):
        raise TypeError("r — только int/Q, float запрещён")
    r = Q(r)
    v = two_cos(int(r.numerator), 2 * int(r.denominator))
    return v / 2


def sin_pi(r):
    """sin(π·r) = cos(π(r − 1/2)), r ∈ ℚ."""
    if isinstance(r, float) or not isinstance(r, _NUM):
        raise TypeError("r — только int/Q, float запрещён")
    return cos_pi(Q(r) - Q(1, 2))


def unity_sum(n: int, ks, offset=None):
    """Σ_k c_k·ζ_n^k (+ offset ∈ ℚ) как пара (re, im) БЕЗ роста степени:
    ks — засечки (последовательность, вес 1) либо {k: c_k} с c_k ∈ ℚ —
    элемент группового кольца ℚ[ζ_n]. cos(2πk/n) = V_{4k}(θ)/2 и
    sin(2πk/n) = V_{4k−n}(θ)/2 при θ = 2cos(2π/N), N = 4n — все слагаемые
    лежат в одном поле ℚ(θ) степени φ(N)/2, сумма — многочлен от θ по
    модулю Ψ_N, значение — одним результантом (AlgebraicNumber.eval_poly).
    Последовательное сложение узлами класса 3 перемножало бы степени
    (замер 04.09: Σζ₇^k не досчитывался за 120 с; произведение двух
    ζ₇-сумм узлами — степени 67/145, 14 с)."""
    from .intfactor import _zdivmod
    weights = dict(ks) if isinstance(ks, dict) else {}
    if not isinstance(ks, dict):
        for k in ks:
            weights[k % n] = weights.get(k % n, Q(0)) + 1
    N = 4 * n
    theta = two_cos(1, N)
    p_n = psi(N)
    d = p_n.degree

    def poly_of(idx: int) -> list:
        j = idx % N
        j = min(j, N - j)
        v = lucas_poly(j)
        if v.degree >= d:
            v = IntPoly(_zdivmod(v.c, p_n.c)[1])
        return v.c

    re_c, im_c = [Q(0)] * max(d, 1), [Q(0)] * max(d, 1)
    for k, c in weights.items():
        c = Q(c)
        if c == 0:
            continue
        for i, a in enumerate(poly_of(4 * k)):
            re_c[i] += c * a
        for i, a in enumerate(poly_of(4 * k - n)):
            im_c[i] += c * a

    def value(coeffs):
        coeffs = [c / 2 for c in coeffs]
        while len(coeffs) > 1 and coeffs[-1] == 0:
            coeffs.pop()
        if isinstance(theta, Q):
            acc = Q(0)
            for c in reversed(coeffs):
                acc = acc * theta + c
            return acc
        return theta.eval_poly(coeffs)

    re, im = value(re_c), value(im_c)
    if offset is not None and offset != 0:
        re = re + Q(offset)
    return re, im
