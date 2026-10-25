"""Приведение решётки по Ленстре–Ленстре–Ловасу — точно, на целых.

Стена, названная §441: перебор подмножеств у Зассенхауса экспоненциален
по числу модульных множителей r, и цена замерена (r = 24 — около 35 с,
r = 32 — около 2,4 часа). Место, где перебор заменяется решением задачи о
коротком векторе, известно с 1982 года: LLL. Здесь он и написан.

Всё точно. Ортогонализация Грама–Шмидта считается парами целых (`Q` —
собственная пара пакета, где дверь из float закрыта построением; у
стандартной `Fraction` она открыта: `Fraction(0.1)` проходит молча),
базис остаётся целочисленным, ни одного float:
приближённая ортогонализация у LLL — известный источник тихих ошибок, а
тихих ошибок здесь быть не должно.

Что даёт LLL и чего не даёт. Даёт: базис, у которого |μᵢⱼ| ≤ ½ и
выполнено условие Ловаса, а первый вектор не длиннее 2^((n−1)/2) от
кратчайшего. Не даёт: сам кратчайший вектор — это задача NP-трудная, и
здесь она не решается и не заявляется.

Проверяемость устроена так, что «привело» нельзя сказать без доказа́ния:
`is_reduced` проверяет оба условия приведённости; `gram_det` считает
определитель Грама точно (метод Барейса, без дробей), и он ОБЯЗАН
совпасть до и после — решётка та же, сменился только базис.

Авторство: A. K. Lenstra, H. W. Lenstra Jr., L. Lovász, «Factoring
polynomials with rational coefficients», Math. Ann. 261 (1982)
[внешний источник]; определитель Барейса — Bareiss 1968.
"""
from __future__ import annotations

from .rational import Q, nearest_int

DELTA = Q(3, 4)          # δ Ловаса; классическое значение 3/4


def _dot(u, v) -> Q:
    return sum((a * b for a, b in zip(u, v)), Q(0))


def gram_schmidt(basis) -> tuple:
    """Ортогонализация Грама–Шмидта точными дробями: (b*, μ, ‖b*‖²)."""
    star, mu, norms = [], [], []
    for i, b in enumerate(basis):
        v = [Q(x) for x in b]
        row = []
        for j in range(i):
            m = _dot([Q(x) for x in b], star[j]) / norms[j] if norms[j] else Q(0)
            row.append(m)
            v = [vi - m * sj for vi, sj in zip(v, star[j])]
        star.append(v)
        mu.append(row)
        norms.append(_dot(v, v))
    return star, mu, norms


def is_reduced(basis, delta: Q = DELTA) -> bool:
    """Приведён ли базис: |μᵢⱼ| ≤ ½ и условие Ловаса на каждом шаге."""
    _, mu, norms = gram_schmidt(basis)
    for i in range(len(basis)):
        for j in range(i):
            if abs(mu[i][j]) > Q(1, 2):
                return False
    for k in range(1, len(basis)):
        if norms[k] < (delta - mu[k][k - 1] ** 2) * norms[k - 1]:
            return False
    return True


#: Ближайшее целое; половина — от нуля. Правило одно на весь пакет и
#: живёт при паре (`rational.nearest_int`), а не переписывается здесь.
_round_half = nearest_int


def lll_rational(basis, delta: Q = DELTA) -> list:
    """Приведение точными дробями — ЭТАЛОН устройства, не рабочий путь.

    Ортогонализация считается один раз, дальше μ и ‖b*‖² правятся по
    месту замкнутыми формулами [внешний источник: H. Cohen, алг. 2.6.3].
    Держится ради сверки: целочисленный путь обязан давать то же самое.
    """
    b = [list(map(int, row)) for row in basis]
    n = len(b)
    if n < 2:
        return b
    _, mu, B = gram_schmidt(b)
    mu = [row + [Q(0)] * (n - len(row)) for row in mu]

    def red(k, l):
        if abs(mu[k][l]) <= Q(1, 2):
            return
        q = _round_half(mu[k][l])
        b[k] = [x - q * y for x, y in zip(b[k], b[l])]
        mu[k][l] -= q
        for i in range(l):
            mu[k][i] -= q * mu[l][i]

    def swap(k):
        b[k], b[k - 1] = b[k - 1], b[k]
        m = mu[k][k - 1]
        bb = B[k] + m * m * B[k - 1]
        if bb == 0:
            B[k - 1], B[k] = B[k], B[k - 1]
        else:
            mu[k][k - 1] = m * B[k - 1] / bb
            B[k] = B[k - 1] * B[k] / bb
            B[k - 1] = bb
        for j in range(k - 1):
            mu[k][j], mu[k - 1][j] = mu[k - 1][j], mu[k][j]
        for i in range(k + 1, n):
            t = mu[i][k]
            mu[i][k] = mu[i][k - 1] - m * t
            mu[i][k - 1] = t + mu[k][k - 1] * mu[i][k]

    k = 1
    while k < n:
        red(k, k - 1)
        if B[k] < (delta - mu[k][k - 1] ** 2) * B[k - 1]:
            swap(k)
            k = max(k - 1, 1)
        else:
            for l in range(k - 2, -1, -1):
                red(k, l)
            k += 1
    return b


def _round_div(a: int, b: int) -> int:
    """Ближайшее целое к a/b при b > 0; половина — от нуля."""
    if a < 0:
        return -_round_div(-a, b)
    q, r = divmod(a, b)
    return q + 1 if 2 * r >= b else q


def lll(basis, delta: Q = DELTA) -> list:
    """Приведённый по LLL базис той же решётки — ЦЕЛОЧИСЛЕННЫЙ путь.

    Дробей нет вовсе: вместо μᵢⱼ хранятся целые λᵢⱼ и определители Грама
    dᵢ, а все деления точны [внешний источник: H. Cohen, «A Course in
    Computational Algebraic Number Theory», алг. 2.6.7]. Знаменатели у
    дробного пути растут вместе с определителями, и на решётках ван Хоя
    это и было узким местом: замер 04.09 — приведение для степени 32
    занимало 54,6 с дробями и 10,1 с целыми, для степени 16 — 2,21 с и
    0,31 с. Эталон (`lll_rational`) остаётся и сверяется приёмкой.

    δ здесь только 3/4: условие Ловаса записано целыми как
    4·d_k·d_{k−2} < 3·d_{k−1}² − 4·λ², и другое δ потребовало бы дробей —
    то есть того самого, от чего целочисленный путь и уходит.
    """
    if delta != DELTA:
        raise ValueError("целочисленный путь написан для δ = 3/4; "
                         "иное δ — через lll_rational")
    n = len(basis)
    if n < 2:
        return [list(map(int, r)) for r in basis]
    b = [None] + [list(map(int, r)) for r in basis]          # 1-индексация Коэна
    d = [0] * (n + 2)
    d[0] = 1
    lam = [[0] * (n + 2) for _ in range(n + 2)]

    def dot(u, v):
        return sum(x * y for x, y in zip(u, v))

    def gram(k):
        for j in range(1, k + 1):
            u = dot(b[k], b[j])
            for i in range(1, j):
                u = (d[i] * u - lam[k][i] * lam[j][i]) // d[i - 1]
            if j < k:
                lam[k][j] = u
            else:
                d[k] = u

    def redi(k, l):
        if 2 * abs(lam[k][l]) <= d[l]:
            return
        q = _round_div(lam[k][l], d[l])
        b[k] = [x - q * y for x, y in zip(b[k], b[l])]
        lam[k][l] -= q * d[l]
        for i in range(1, l):
            lam[k][i] -= q * lam[l][i]

    def swapi(k, kmax):
        b[k], b[k - 1] = b[k - 1], b[k]
        for j in range(1, k - 1):
            lam[k][j], lam[k - 1][j] = lam[k - 1][j], lam[k][j]
        L = lam[k][k - 1]
        B = (d[k - 2] * d[k] + L * L) // d[k - 1]
        for i in range(k + 1, kmax + 1):
            t = lam[i][k]
            lam[i][k] = (d[k] * lam[i][k - 1] - L * t) // d[k - 1]
            lam[i][k - 1] = (B * t + L * lam[i][k]) // d[k]
        d[k - 1] = B

    d[1] = dot(b[1], b[1])
    k, kmax = 2, 1
    while k <= n:
        if k > kmax:
            kmax = k
            gram(k)
        while True:
            redi(k, k - 1)
            if 4 * d[k] * d[k - 2] < 3 * d[k - 1] * d[k - 1] - 4 * lam[k][k - 1] ** 2:
                swapi(k, kmax)
                k = max(2, k - 1)
            else:
                break
        for l in range(k - 2, 0, -1):
            redi(k, l)
        k += 1
    return [b[i] for i in range(1, n + 1)]


def self_check(deep: bool = False) -> list:
    """Самопроверка приведения решётки ВСТРЕЧНОЙ дорогой (ось самоприменения).

    `[вычислено здесь]` Приведение — целая процедура с обменами и
    редукциями; сверять её пересчётом той же процедурой бессмысленно,
    поэтому дороги внешние ей: (1) ПРЕДИКАТ приведённости `is_reduced`
    (|μᵢⱼ| ≤ ½ и условие Ловаса) считается по Граму–Шмидту заново и обязан
    быть истинным на выходе; (2) ОПРЕДЕЛИТЕЛЬ ГРАМА (метод Барейса, без
    дробей) не меняется при смене базиса решётки — переход унимодулярен, и
    совпадение до и после есть доказательство, что решётка ТА ЖЕ; (3)
    членство: каждый вектор нового базиса обязан лежать в старой решётке
    (`in_lattice` решает нормальными уравнениями — третья дорога, ни на
    приведение, ни на определитель не похожая); (4) целая дорога `lll` и
    дробный эталон `lll_rational` обязаны дать один базис.

    Граница объявлена: проверяются пробные базисы малой размерности
    (2…4, при `deep` — и вырожденный по длине случай); LLL гарантирует
    приведённость при своём δ, и другое δ здесь не утверждается. Решётка
    с линейно зависимым базисом в дверь не входит — это условие `lll`, и
    молчание о ней находкой не считается. Возвращает список находок;
    исключений не бросает.
    """
    out = []
    bases = [[[1, 1, 1], [-1, 0, 2], [3, 5, 6]],
             [[201, 37], [1648, 297]],
             [[1, 0, 0, 1000], [0, 1, 0, 500], [0, 0, 1, 250], [0, 0, 0, 1]]]
    if deep:
        bases.append([[15, 23, 11], [46, 15, 17], [32, 1, 0]])
    for b in bases:
        try:
            r = lll(b)
        except Exception as e:                                  # noqa: BLE001
            out.append(f"приведение отказало на {b}: {type(e).__name__}: {e}")
            continue
        if not is_reduced(r):
            out.append(f"выход не приведён по предикату Ловаса: {r}")
        d0, d1 = gram_det(b), gram_det(r)
        if abs(d0) != abs(d1):
            out.append(f"определитель Грама изменился: {d0} → {d1} — решётка "
                       f"не та же")
        for v in r:
            if not in_lattice(v, b):
                out.append(f"вектор выхода {v} не лежит в исходной решётке")
        if norm2(r[0]) > norm2(b[0]) and is_reduced(b):
            out.append(f"первый вектор стал длиннее на уже приведённом базисе: "
                       f"{norm2(b[0])} → {norm2(r[0])}")
        try:
            if [list(map(int, v)) for v in lll_rational(b)] != [list(map(int, v)) for v in r]:
                out.append(f"целая дорога и дробный эталон дали разные базисы "
                           f"на {b}")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"эталон lll_rational отказал: {type(e).__name__}: {e}")
    return out


def gram_det(basis) -> int:
    """Определитель матрицы Грама — целое, метод Барейса (без дробей).

    Он не меняется при смене базиса решётки: у унимодулярного перехода
    определитель ±1, а Грам квадратичен. Потому совпадение до и после —
    доказательство, что решётка та же.
    """
    n = len(basis)
    m = [[sum(a * c for a, c in zip(basis[i], basis[j])) for j in range(n)]
         for i in range(n)]
    prev = 1
    sign = 1
    for k in range(n - 1):
        if m[k][k] == 0:
            for r in range(k + 1, n):
                if m[r][k] != 0:
                    m[k], m[r] = m[r], m[k]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                m[i][j] = (m[i][j] * m[k][k] - m[i][k] * m[k][j]) // prev
        prev = m[k][k]
    return sign * m[n - 1][n - 1]


def norm2(v) -> int:
    return sum(int(x) * int(x) for x in v)


def in_lattice(v, basis) -> bool:
    """Лежит ли вектор в решётке базиса (базис линейно независим).

    Решается точно: коэффициенты находятся из нормальных уравнений
    Грама рациональной прогонкой и обязаны быть целыми.
    """
    n = len(basis)
    g = [[Q(sum(a * c for a, c in zip(basis[i], basis[j]))) for j in range(n)]
         + [Q(sum(a * c for a, c in zip(basis[i], v)))] for i in range(n)]
    for k in range(n):
        p = next((r for r in range(k, n) if g[r][k] != 0), None)
        if p is None:
            return False
        g[k], g[p] = g[p], g[k]
        pivot = g[k][k]
        g[k] = [x / pivot for x in g[k]]
        for r in range(n):
            if r != k and g[r][k] != 0:
                f = g[r][k]
                g[r] = [x - f * y for x, y in zip(g[r], g[k])]
    coef = [g[i][n] for i in range(n)]
    if any(c.denominator != 1 for c in coef):
        return False
    got = [sum(int(coef[i]) * basis[i][j] for i in range(n)) for j in range(len(v))]
    return got == list(v)
