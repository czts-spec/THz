"""Факторизация над ℤ — стена, названная пересмотром отказов 1602cda7
(аннулятор класса 3 приводим при степени ≥ 4, а приводимость там не
решалась). Дорога Зассенхауса [внешний источник: von zur Gathen, Gerhard,
«Modern Computer Algebra», гл. 14–15]: факторизация по модулю простого
(различные степени + Кантор–Цассенхаус), подъём Гензеля до границы
Миньотта, перебор подмножеств модульных множителей с симметричным
остатком и пробным делением над ℤ.

Отсев подмножеств (04.09): свободный член произведения обязан делить
свободный член f (Гаусс). Он считается по подмножеству за O(s)
умножений — против O(deg²) у пробного деления, — и отсев ТОЧЕН, а не
вероятностен. Замер на Суиннертон-Дайере степени 32: 39 202 подмножества
просмотрено, до перемножения дошло 256 (0,65 %), время 5,089 → 0,137 с.
Переключатель `PRUNE_BY_CONSTANT` снимает отсев, чтобы выигрыш
замерялся, а не заявлялся; `LAST_SEARCH` хранит счёт последнего перебора.

Всё на целых: коэффициенты int, остатки по модулю — int, ни одного
float. Случайность Кантора–Цассенхауса — детерминированная (семя
фиксировано): один и тот же вход даёт один и тот же ход. Перебор
подмножеств экспоненциален по числу модульных множителей r — это
честная граница дороги. За порогом (04.09) идёт решётка следов
(`vanhoeij.py`): перебора там нет вовсе, а всякий найденный набор
проверяется делением. Если и она разбиения не собрала — явный отказ
(Unresolved), а не молчание. Самопроверка (04.09):
произведение множителей сверяется со входом внутри factor_squarefree и
factor_int (verify_product) — расхождение есть отказ.
"""
from __future__ import annotations

import random
from math import isqrt

from .errors import Unresolved
from .intpoly import IntPoly

# Счётчик последнего перебора — чтобы цена дороги была замеряема, а не
# рассказана: сколько подмножеств просмотрено и сколько дошло до
# перемножения множителей.
LAST_SEARCH: dict = {}

# Отсев по свободному члену можно снять — тогда перебор идёт как прежде.
# Держится переключателем, чтобы выигрыш замерялся, а не заявлялся.
PRUNE_BY_CONSTANT = True

# ---------- многочлены по модулю p: списки int от младшей степени ----------


def _norm(a):
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    if not a:
        a.append(0)
    return a


def _madd(a, b, m):
    n = max(len(a), len(b))
    return _norm([((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % m for i in range(n)])


def _msub(a, b, m):
    n = max(len(a), len(b))
    return _norm([((a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)) % m for i in range(n)])


def _mmul(a, b, m):
    if a == [0] or b == [0]:
        return [0]
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return _norm([c % m for c in out])


def _mscale(a, k, m):
    return _norm([(c * k) % m for c in a])


def _mdivmod(a, b, m):
    """Деление с остатком по модулю m; старший коэффициент b обратим
    (m простое, либо b — унитарный)."""
    a = _norm([c % m for c in a])
    b = _norm([c % m for c in b])
    if b == [0]:
        raise ZeroDivisionError
    inv = pow(b[-1], -1, m)
    r = list(a)
    if len(r) < len(b):
        return [0], r
    q = [0] * (len(r) - len(b) + 1)
    for k in range(len(q) - 1, -1, -1):
        c = (r[len(b) - 1 + k] * inv) % m
        q[k] = c
        if c:
            for i, bc in enumerate(b):
                r[i + k] = (r[i + k] - c * bc) % m
    return _norm(q), _norm(r[: len(b) - 1] or [0])


def _mmod(a, b, m):
    return _mdivmod(a, b, m)[1]


def _monic(a, p):
    inv = pow(a[-1], -1, p)
    return _mscale(a, inv, p)


def _mgcd(a, b, p):
    a, b = _norm([c % p for c in a]), _norm([c % p for c in b])
    while b != [0]:
        a, b = b, _mmod(a, b, p)
    return _monic(a, p) if a != [0] else [0]


def _mxgcd(a, b, p):
    """Расширенный Евклид в 𝔽_p[x]: (g, s, t), s·a + t·b = g, g унитарный."""
    r0, r1 = _norm([c % p for c in a]), _norm([c % p for c in b])
    s0, s1 = [1], [0]
    t0, t1 = [0], [1]
    while r1 != [0]:
        q, r = _mdivmod(r0, r1, p)
        r0, r1 = r1, r
        s0, s1 = s1, _msub(s0, _mmul(q, s1, p), p)
        t0, t1 = t1, _msub(t0, _mmul(q, t1, p), p)
    inv = pow(r0[-1], -1, p)
    return _mscale(r0, inv, p), _mscale(s0, inv, p), _mscale(t0, inv, p)


def _mpowmod(base, e, f, p):
    """base^e mod f в 𝔽_p[x]."""
    result = [1]
    base = _mmod(base, f, p)
    while e:
        if e & 1:
            result = _mmod(_mmul(result, base, p), f, p)
        base = _mmod(_mmul(base, base, p), f, p)
        e >>= 1
    return result


def _mderiv(a, p):
    if len(a) <= 1:
        return [0]
    return _norm([(i * c) % p for i, c in enumerate(a)][1:])


# ---------- факторизация в 𝔽_p[x] ----------


def _distinct_degree(f, p):
    """Разложение по различным степеням: список (g_d, d) — g_d есть
    произведение всех неприводимых делителей f степени ровно d."""
    out = []
    h = list(f)
    xpow = [0, 1]                      # x^{p^i} mod h, i = 0
    d = 0
    while len(h) - 1 >= 2 * (d + 1):
        d += 1
        xpow = _mpowmod(xpow, p, h, p)   # x^{p^d} mod h
        g = _mgcd(_msub(xpow, [0, 1], p), h, p)
        if len(g) > 1:
            out.append((g, d))
            h = _mdivmod(h, g, p)[0]
            xpow = _mmod(xpow, h, p)
    if len(h) > 1:
        out.append((h, len(h) - 1))
    return out


def _equal_degree(f, d, p, rng):
    """Кантор–Цассенхаус: f — унитарное произведение неприводимых
    степени d в 𝔽_p[x], p нечётное; возвращает список множителей."""
    n = len(f) - 1
    if n == d:
        return [f]
    exp = (p ** d - 1) // 2
    while True:
        a = [rng.randrange(p) for _ in range(n)]
        a = _norm(a)
        if len(a) <= 1:
            continue
        g = _mgcd(a, f, p)
        if 1 < len(g) < len(f):
            break
        b = _msub(_mpowmod(a, exp, f, p), [1], p)
        g = _mgcd(b, f, p)
        if 1 < len(g) < len(f):
            break
    other = _mdivmod(f, g, p)[0]
    return _equal_degree(g, d, p, rng) + _equal_degree(other, d, p, rng)


def factor_mod_p(f, p, seed=0):
    """Унитарные неприводимые множители унитарного свободного от квадратов
    f в 𝔽_p[x] (p нечётное простое)."""
    rng = random.Random(seed)
    out = []
    for g, d in _distinct_degree(f, p):
        out.extend(_equal_degree(g, d, p, rng))
    return sorted(out, key=lambda g: (len(g), g))


# ---------- подъём Гензеля ----------


def _hensel_step(f, g, h, s, t, m):
    """Квадратичный шаг [von zur Gathen, Gerhard, алг. 15.10]: из
    f ≡ g·h (mod m), s·g + t·h ≡ 1 (mod m), g и h унитарные — те же
    соотношения по модулю m². f унитарный."""
    m2 = m * m
    e = _msub(f, _mmul(g, h, m2), m2)
    q, r = _mdivmod(_mmul(s, e, m2), h, m2)
    g1 = _madd(_madd(g, _mmul(t, e, m2), m2), _mmul(q, g, m2), m2)
    h1 = _madd(h, r, m2)
    b = _msub(_madd(_mmul(s, g1, m2), _mmul(t, h1, m2), m2), [1], m2)
    c, d = _mdivmod(_mmul(s, b, m2), h1, m2)
    s1 = _msub(s, d, m2)
    t1 = _msub(_msub(t, _mmul(t, b, m2), m2), _mmul(c, g1, m2), m2)
    return g1, h1, s1, t1


def _hensel_lift(f, factors, p, target):
    """Многомножительный подъём деревом [алг. 15.17]: f унитарный над ℤ,
    factors — его унитарные попарно взаимно простые множители по модулю p;
    возвращает (m, список множителей по модулю m), m = p^(2^j) ≥ target."""
    m = p
    while m < target:
        m *= m
    return m, _lift_tree(f, factors, p, m)


def _lift_tree(f, factors, p, mtarget):
    if len(factors) == 1:
        return [_norm([c % mtarget for c in f])]
    half = len(factors) // 2
    left, right = factors[:half], factors[half:]
    g = [1]
    for fac in left:
        g = _mmul(g, fac, p)
    h = [1]
    for fac in right:
        h = _mmul(h, fac, p)
    _, s, t = _mxgcd(g, h, p)
    m = p
    fmod = _norm([c % (mtarget) for c in f])
    while m < mtarget:
        g, h, s, t = _hensel_step(_norm([c % (m * m) for c in fmod]), g, h, s, t, m)
        m *= m
    return _lift_tree(g, left, p, mtarget) + _lift_tree(h, right, p, mtarget)


# ---------- над ℤ ----------


def _symmetric(a, m):
    half = m // 2
    return _norm([c - m if c > half else c for c in a])


def _zdivmod(a, b):
    """Деление в ℤ[x] с остатком по унитарному b: (q, r) точно."""
    r = list(a)
    if len(r) < len(b):
        return [0], _norm(r)
    q = [0] * (len(r) - len(b) + 1)
    for k in range(len(q) - 1, -1, -1):
        c = r[len(b) - 1 + k]
        q[k] = c
        if c:
            for i, bc in enumerate(b):
                r[i + k] -= c * bc
    return _norm(q), _norm(r[: len(b) - 1] or [0])


def _norm2_sq(a):
    return sum(c * c for c in a)


def _mignotte(fc):
    """Граница Миньотта для унитарного f степени n: любой унитарный
    делитель имеет |коэффициенты| ≤ 2^n · ‖f‖₂ [внешний источник:
    Mignotte 1974; von zur Gathen, Gerhard, следствие 6.33]."""
    n = len(fc) - 1
    return (1 << n) * (isqrt(_norm2_sq(fc)) + 1)


def _small_primes(limit=2000):
    sieve = bytearray([1]) * (limit + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, isqrt(limit) + 1):
        if sieve[i]:
            sieve[i * i:: i] = bytearray(len(sieve[i * i:: i]))
    return [i for i in range(3, limit + 1) if sieve[i]]


_PRIMES = _small_primes()

MAX_MODULAR = 24    # порог перебора подмножеств: 2^(r−1) ≤ 2^23


def _choose_prime(fc, tries=5):
    """Простое, при котором f унитарный свободен от квадратов по модулю p;
    из нескольких пригодных — то, где модульных множителей меньше всего."""
    best = None
    seen = 0
    for p in _PRIMES:
        fp = _norm([c % p for c in fc])
        if len(fp) != len(fc):
            continue
        if len(_mgcd(fp, _mderiv(fp, p), p)) != 1:
            continue
        facs = factor_mod_p(fp, p)
        seen += 1
        if best is None or len(facs) < len(best[1]):
            best = (p, facs)
        if len(facs) == 1 or seen >= tries:
            break
    if best is None:
        raise Unresolved("факторизация", f"IntPoly({fc})", hint="ни одно простое < 2000 не даёт свободный от квадратов образ")
    return best


def _factor_monic_squarefree(fc):
    """Неприводимые унитарные множители унитарного свободного от квадратов
    f ∈ ℤ[x] (список списков коэффициентов)."""
    if len(fc) <= 2:
        return [list(fc)]
    p, modular = _choose_prime(fc)
    if len(modular) == 1:
        return [list(fc)]
    if len(modular) > MAX_MODULAR:
        # За порогом перебора идёт решётка следов (ван Хой). Ввоз здесь, а не
        # наверху: vanhoeij опирается на этот же модуль, и кольцевой ввоз при
        # загрузке пакета был бы платой за красоту строки.
        from .vanhoeij import factor_squarefree as _by_lattice
        got = _by_lattice(IntPoly(fc))
        if got is not None:
            return [list(g.c) for g in got]
        raise Unresolved("факторизация", f"IntPoly({fc})",
                         hint=f"модульных множителей {len(modular)} > {MAX_MODULAR}: "
                              f"перебор подмножеств за порогом, решётка следов "
                              f"разбиения не собрала")
    bound = 2 * _mignotte(fc)
    m, lifted = _hensel_lift(fc, modular, p, bound)
    remaining = list(range(len(lifted)))
    f = list(fc)
    found = []
    s = 1
    trials = tried = 0
    while 2 * s <= len(remaining):
        progressed = False
        # свободный член произведения считается один раз по подмножеству: он
        # обязан делить свободный член f (Гаусс). Отсев точный, не
        # вероятностный, и стоит O(s) умножений против O(deg²) у пробного
        # деления — на нём и держится вся цена перебора.
        c0 = [g0[0] for g0 in lifted]
        f0 = f[0]
        for combo in _combinations(remaining, s):
            trials += 1
            if PRUNE_BY_CONSTANT and f0:
                t = 1
                for i in combo:
                    t = t * c0[i] % m
                t = t - m if t > m // 2 else t
                if t == 0 or f0 % t != 0:
                    continue
            tried += 1
            g = [1]
            for i in combo:
                g = _mmul(g, lifted[i], m)
            g = _symmetric(g, m)
            if g[-1] != 1:
                continue
            q, r = _zdivmod(f, g)
            if r == [0]:
                found.append(g)
                f = q
                remaining = [i for i in remaining if i not in combo]
                progressed = True
                break
        if not progressed:
            s += 1
    if len(f) > 1:
        found.append(f)
    LAST_SEARCH.update({"подмножеств": trials, "дошло до умножения": tried,
                        "модульных множителей": len(lifted)})
    return found


def _combinations(items, s):
    from itertools import combinations
    return combinations(items, s)


def _monicize(pc):
    """f(x) с lc = a → F(y) = a^(n−1)·f(y/a) унитарный над ℤ."""
    n = len(pc) - 1
    a = pc[-1]
    return [pc[i] * a ** (n - 1 - i) for i in range(n)] + [1]   # старший: a·a^(−1) = 1, без дроби


def _demonicize(Gc, a):
    """Множитель G(y) унитарного F → примитивная часть G(a·x)."""
    return IntPoly([c * a ** i for i, c in enumerate(Gc)]).primitive()


def verify_product(p: IntPoly, factors: list, mults=None) -> None:
    """Самопроверка факторизации («AIAΣA спрашивает AIAΣA», 04.09):
    произведение множителей (с кратностями) обязано равняться примитивной
    части входа с точностью до знака; иначе — отказ Unresolved, не тихий
    ответ. Цена — одно произведение многочленов."""
    prod = IntPoly([1])
    for i, f in enumerate(factors):
        for _ in range(mults[i] if mults else 1):
            prod = prod * f
    want = p.primitive()
    if prod != want and prod != -want:
        raise Unresolved("факторизация", f"IntPoly({p.c})",
                         hint=f"самопроверка: произведение множителей {prod.c} ≠ вход {want.c}")


def self_check(deep: bool = False) -> list:
    """Самопроверка факторизации ВСТРЕЧНОЙ дорогой (ось самоприменения).

    `[вычислено здесь]` Дороги: (1) ОТВЕТ ЗАДАН ЗАРАНЕЕ — здесь берутся
    заведомо неприводимые множители, ПЕРЕМНОЖАЮТСЯ, и факторизация обязана
    вернуть ровно тот же набор (сравнение множеством многочленов, а не
    числом множителей): разложение сверяется со сборкой, которую оно не
    видело; (2) произведение множителей обязано равняться примитивной части
    входа (`verify_product` — умножение обратно, дорога, обратная делению
    и подъёму Гензеля); (3) степени: сумма степеней множителей с кратностями
    равна степени входа; (4) факторизация по модулю: произведение
    `factor_mod_p` обязано сойтись с УНИТАРНЫМ приведением входа по тому же
    модулю — дорога, не похожая ни на подъём Гензеля, ни на сборку.

    Граница объявлена: проверяются пробные многочлены малой степени
    (при `deep` добавляется степень 6 и кратности); полнота неприводимости
    множителей здесь НЕ доказывается — она есть свойство алгоритма, и
    самопроверка говорит о сборке, а не о неприводимости. Многочлен с
    нулевым содержанием и нулевой многочлен в дверь не входят — их отбивает
    сама факторизация; модульная дорога пропускается, когда старший
    коэффициент гасится модулем — там многочлен не унитарен, и это условие
    самой двери. Возвращает список находок; исключений не бросает.
    """
    out = []
    irreducibles = [IntPoly([1, 1]), IntPoly([-1, 1]), IntPoly([1, 0, 1]),
                    IntPoly([1, 1, 1]), IntPoly([-2, 0, 1])]
    pairs = [(0, 2), (1, 3), (2, 4), (0, 1)]
    if deep:
        pairs += [(0, 4), (3, 4), (1, 2)]
    for i, j in pairs:
        f, g = irreducibles[i], irreducibles[j]
        prod = f * g
        try:
            got = factor_squarefree(prod)
        except Exception as e:                                  # noqa: BLE001
            out.append(f"факторизация отказала на {prod.c}: "
                       f"{type(e).__name__}: {e}")
            continue
        want = sorted((f, g), key=lambda q: (len(q.c), q.c))
        if sorted(got, key=lambda q: (len(q.c), q.c)) != want:
            out.append(f"заданный ответ не воспроизведён: вход {prod.c} "
                       f"разложен в {[x.c for x in got]}, собран из "
                       f"{[x.c for x in want]}")
        try:
            verify_product(prod, got)
        except Exception as e:                                  # noqa: BLE001
            out.append(f"произведение множителей не сошлось с входом "
                       f"{prod.c}: {e}")
        if sum(x.degree for x in got) != prod.degree:
            out.append(f"сумма степеней множителей {sum(x.degree for x in got)} "
                       f"≠ степени входа {prod.degree}")
        for pmod in (5, 7):
            # дверь модульной дороги принимает СПИСОК коэффициентов унитарного
            # многочлена (не IntPoly) — контракт factor_mod_p; при старшем,
            # гасимом модулем, многочлен там не унитарен и дорога пропускается
            if prod.c[-1] % pmod == 0:
                continue
            monic = _monic([c % pmod for c in prod.c], pmod)
            try:
                mod_factors = factor_mod_p(monic, pmod)
            except Exception as e:                              # noqa: BLE001
                out.append(f"факторизация по модулю {pmod} отказала: "
                           f"{type(e).__name__}: {e}")
                continue
            acc = [1]
            for mf in mod_factors:
                acc = _mmul(acc, mf, pmod)
            if _norm(acc) != _norm(monic):
                out.append(f"произведение модульных множителей по {pmod} не "
                           f"сошлось с унитарным приведением входа: {acc} "
                           f"против {monic}")
    if deep:
        cube = irreducibles[0] * irreducibles[0] * irreducibles[2]
        cont, fs = factor_int(cube)
        if cont != 1:
            out.append(f"содержание унитарного входа {cont} ≠ 1")
        try:
            verify_product(cube, [x for x, _ in fs], [m for _, m in fs])
        except Exception as e:                                  # noqa: BLE001
            out.append(f"кратности не сошлись: {e}")
        if sum(x.degree * m for x, m in fs) != cube.degree:
            out.append("сумма степеней с кратностями ≠ степени входа")
    return out


def factor_squarefree(p: IntPoly) -> list:
    """Неприводимые над ℚ примитивные множители (lc > 0) свободного от
    квадратов p ∈ ℤ[x]; произведение равно primitive(p) с точностью до
    знака. Содержание снимается, множитель x выделяется отдельно."""
    pc = list(p.primitive().c)
    if len(pc) <= 1:
        return []
    out = []
    while pc[0] == 0 and len(pc) > 1:      # множитель x
        out.append(IntPoly([0, 1]))
        pc = pc[1:]
    if len(pc) <= 1:
        return out
    if len(pc) == 2:
        return out + [IntPoly(pc).primitive()]
    a = pc[-1]
    F = _monicize(pc)
    for Gc in _factor_monic_squarefree(F):
        out.append(_demonicize(Gc, a))
    verify_product(p, out)
    return sorted(out, key=lambda q: (len(q.c), q.c))


def factor_int(p: IntPoly) -> list:
    """Полная факторизация над ℤ: (содержание, [(множитель, кратность)…]).
    Кратности — по Юну [внешний источник: Yun 1976] через целые НОД."""
    from .algebraic import _int_poly_gcd, _divexact, _norm as _zn
    cont = p.content()
    if p.is_zero():
        return 0, []
    pc = list(p.primitive().c)
    sign = 1
    if p.c[-1] < 0:
        cont = -cont
    out = []
    # Юн: a = p, b = p′, c = gcd(a, b), w = a / c, y = b / c, z = y − w′
    a = pc
    b = IntPoly(a).derivative().c
    c = _int_poly_gcd(a, b)
    w = _zn(_divexact(a, c))
    y = _zn(_divexact(b, c))
    k = 1
    while w != [1]:
        z = _zn(_divexact_sub(y, IntPoly(w).derivative().c))
        g = _int_poly_gcd(w, z)
        if len(g) > 1:
            for fac in factor_squarefree(IntPoly(g)):
                out.append((fac, k))
        w = _zn(_divexact(w, g))
        y = _zn(_divexact(z, g))
        k += 1
    verify_product(p, [f for f, _ in out], [m for _, m in out])
    return cont * sign, out


def _divexact_sub(y, dw):
    n = max(len(y), len(dw))
    return [((y[i] if i < len(y) else 0) - (dw[i] if i < len(dw) else 0)) for i in range(n)]


def is_irreducible(p: IntPoly) -> bool:
    """Неприводим ли примитивный многочлен над ℚ (степень ≥ 1)."""
    from .algebraic import squarefree_part
    if p.degree < 1:
        return False
    if squarefree_part(p).degree != p.degree:
        return False
    return len(factor_squarefree(p)) == 1
