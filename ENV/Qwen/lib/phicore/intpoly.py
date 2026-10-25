"""Целочисленные многочлены и результант — фундамент линейки класса 3.

Всё точно и всё на целых: коэффициенты int, определитель — настоящий
Барейс (fraction-free, деление нацело без остатка), результант — через
целочисленную матрицу Сильвестра, интерполяция Лагранжа — с очисткой
знаменателей до целых. Ни одного float. Рациональные значения
принимаются на границе по паре .numerator/.denominator (int этой парой
тоже обладает) и немедленно очищаются в целые; накопленный множитель
учитывается в результате.
"""
from __future__ import annotations

from math import gcd, lcm


def _as_int(a, i: int) -> int:
    """Коэффициент — только целое; всё прочее отвергается ПОИМЁННО.

    `[вычислено здесь]` Здесь стоял `int(a)`, и дверь не держала ничего:
    замер 05.09 — `IntPoly([0.5, 1]).c` давал `[0, 1]` (половина
    усечена молча), `IntPoly([Q(1, 2)]).c` — `[0]`, то есть НУЛЕВОЙ
    многочлен вместо x/2, `IntPoly(['1'])` — разбор строки, `True` — 1.
    Минимальный многочлен мог молча стать нулём. Усечение — тот самый
    молчаливый перенос, которого в пакете нет: теперь оно отказ.

    Целое, у которого есть точная пара с единичным знаменателем (`Q(3)`,
    `Fraction(6, 2)`), принимается — это то же целое, записанное иначе, а
    не приближение. Дробное с знаменателем ≠ 1 — отказ.
    """
    if isinstance(a, bool):
        raise TypeError(f"коэффициент {i}: bool не число")
    if isinstance(a, int):
        return a
    if isinstance(a, float):
        raise TypeError(f"коэффициент {i}: float запрещён построением "
                        f"({a!r} усекалось бы молча)")
    num, den = getattr(a, "numerator", None), getattr(a, "denominator", None)
    if isinstance(num, int) and isinstance(den, int) and not isinstance(num, bool):
        if den == 1:
            return num
        raise TypeError(f"коэффициент {i}: {a!r} — не целое (знаменатель {den}); "
                        f"очистите знаменатели явно")
    raise TypeError(f"коэффициент {i}: {type(a).__name__} не целое")


class IntPoly:
    """Многочлен из ℤ[x]; coeffs — от младшей степени к старшей."""

    __slots__ = ("c",)

    def __init__(self, coeffs):
        c = [_as_int(a, i) for i, a in enumerate(coeffs)]
        while len(c) > 1 and c[-1] == 0:
            c.pop()
        self.c = c

    # --- устройство ---
    @property
    def degree(self) -> int:
        return len(self.c) - 1 if any(self.c) else 0

    def is_zero(self) -> bool:
        return all(a == 0 for a in self.c)

    def __eq__(self, other) -> bool:
        return isinstance(other, IntPoly) and self.c == other.c

    def __repr__(self) -> str:
        terms = [f"{a}·x^{i}" for i, a in enumerate(self.c) if a]
        return " + ".join(terms) or "0"

    # --- кольцо ---
    def __add__(self, o: "IntPoly") -> "IntPoly":
        n = max(len(self.c), len(o.c))
        return IntPoly([(self.c[i] if i < len(self.c) else 0)
                        + (o.c[i] if i < len(o.c) else 0) for i in range(n)])

    def __neg__(self) -> "IntPoly":
        return IntPoly([-a for a in self.c])

    def __sub__(self, o: "IntPoly") -> "IntPoly":
        return self + (-o)

    def __mul__(self, o: "IntPoly") -> "IntPoly":
        out = [0] * (len(self.c) + len(o.c) - 1)
        for i, a in enumerate(self.c):
            if a:
                for j, b in enumerate(o.c):
                    out[i + j] += a * b
        return IntPoly(out)

    def scale(self, k: int) -> "IntPoly":
        return IntPoly([a * k for a in self.c])

    # --- значения и производная ---
    def eval(self, x):
        """Значение в точке схемой Горнера; точно для int и для дроби
        (дробь на входе — дробь на выходе, целое — целое)."""
        acc = 0
        for a in reversed(self.c):
            acc = acc * x + a
        return acc

    def derivative(self) -> "IntPoly":
        if len(self.c) <= 1:
            return IntPoly([0])
        return IntPoly([i * a for i, a in enumerate(self.c)][1:])

    # --- нормализация ---
    def content(self) -> int:
        g = 0
        for a in self.c:
            g = gcd(g, abs(a))
        return g or 1

    def primitive(self) -> "IntPoly":
        g = self.content()
        p = IntPoly([a // g for a in self.c])
        if p.c[-1] < 0:
            p = -p
        return p

    # --- палиндромы: зазеркалье x ↔ 1/x (слово автора 04.09) ---
    def reciprocal(self) -> "IntPoly":
        """x^n·P(1/x): коэффициенты в обратном порядке — корни переходят в
        обратные; аннулятор 1/θ по аннулятору θ без единого результанта."""
        return IntPoly(list(reversed(self.c)))

    def is_palindromic(self) -> bool:
        """P = x^n·P(1/x): множество корней замкнуто относительно θ ↦ 1/θ."""
        return self.c == list(reversed(self.c))

    def is_antipalindromic(self) -> bool:
        """P = −x^n·P(1/x)."""
        return self.c == [-a for a in reversed(self.c)]

    def palindromic_reduce(self) -> "IntPoly":
        """Палиндром чётной степени 2h: P(x) = x^h·Q(x + 1/x), deg Q = h —
        степень делится пополам подстановкой зазеркалья y = x + 1/x
        (x^k + x^{−k} = V_k(y), многочлены Люка). Нечётная степень или не
        палиндром — ValueError (у нечётного палиндрома корень −1 снимается
        делением на x + 1)."""
        if not self.is_palindromic():
            raise ValueError("не палиндром")
        n = self.degree
        if n % 2:
            raise ValueError("нечётная степень: сперва снять множитель x + 1")
        h = n // 2
        # V_0 = 2, V_1 = y, V_k = y·V_{k−1} − V_{k−2}
        v_prev, v_cur = IntPoly([2]), IntPoly([0, 1])
        lucas = [v_prev, v_cur]
        for _ in range(2, h + 1):
            v_prev, v_cur = v_cur, IntPoly([0, 1]) * v_cur - v_prev
            lucas.append(v_cur)
        q = IntPoly([self.c[h]])
        for k in range(1, h + 1):
            q = q + lucas[k].scale(self.c[h + k])
        return q

    # --- подстановки для результантов ---
    def compose_x_minus_y_coeff_in_y(self, x_val) -> list:
        """P(x_val − y) как список коэффициентов по y; целая точка даёт
        целые коэффициенты, дробная — дробные (граница)."""
        # (x−y)^k раскрывается биномом; собираем по степеням y.
        n = len(self.c)
        out = [0] * n
        for k, a in enumerate(self.c):
            if a == 0:
                continue
            # (x−y)^k = Σ C(k,i)·x^(k−i)·(−y)^i
            binom = 1
            for i in range(k + 1):
                out[i] += a * binom * (x_val ** (k - i)) * ((-1) ** i)
                binom = binom * (k - i) // (i + 1)
        return out

    def self_check(self, deep: bool = False) -> list:
        """Самопроверка многочлена из ℤ[x] независимой дорогой (очередь п. 18).

        `[вычислено здесь]` Проверяется не тем счётом, каким строилось:
        коэффициенты целы (bool и float отвергнуты дверью `_as_int`), канон
        записи без хвостовых нулей, `degree` и `is_zero` отвечают самим
        коэффициентам, содержание делит каждый коэффициент нацело, а
        примитив восстанавливает многочлен обратным умножением. При `deep` —
        пять встречных дорог: значение схемой Горнера против явной степенной
        суммы Σcᵢ·tⁱ; производная против тождества (P²)′ = 2·P·P′;
        зазеркалье против тождества R(1/t)·tⁿ = P(t); результант через
        целую матрицу Сильвестра и Барейса против Горнера —
        Res(P, y−a) = (−1)^deg P·P(a) (знак замерен здесь: у нечётной степени
        он −1, на чётных степенях порча знака не видна); у палиндрома чётной
        степени — подстановка Люка против тождества Q(t + 1/t)·tʰ = P(t).

        Граница: неприводимость и разложение здесь НЕ проверяются — это
        работа `intfactor` (сверка произведением стоит там). Значение
        сверяется в трёх рациональных точках, а не тождественно: их меньше,
        чем deg+1, и совпадение в них тождества не заменяет — граница
        объявлена числом, а не словом «проверено».
        """
        out = []
        if not self.c:
            out.append("пустая запись: у нулевого многочлена стоит [0], а не []")
            return out
        for i, a in enumerate(self.c):
            if isinstance(a, bool) or not isinstance(a, int):
                out.append(f"коэффициент {i} не целый: {type(a).__name__} ({a!r})")
        if out:
            return out                      # дальше считать нечем: дверь уже пробита
        if len(self.c) > 1 and self.c[-1] == 0:
            out.append(f"канон записи нарушен: хвостовой нуль при x^{len(self.c) - 1}")
        want = len(self.c) - 1 if any(self.c) else 0
        if self.degree != want:
            out.append(f"degree {self.degree} не отвечает коэффициентам (ждётся {want})")
        if self.is_zero() != (not any(self.c)):
            out.append("is_zero не отвечает коэффициентам")
        g = self.content()
        if g < 1:
            out.append(f"содержание {g} не положительно")
        elif any(a % g for a in self.c):
            out.append(f"содержание {g} не делит коэффициенты нацело")
        else:
            pr = self.primitive()
            if pr.scale(g) != self and (-pr).scale(g) != self:
                out.append(f"примитив обратным умножением даёт {pr.scale(g)!r}, а не сам многочлен")
        if not deep:
            return out
        from .rational import Q                     # дверь из целых в пары — только здесь
        pts = (Q(3, 5), Q(-2), Q(7, 2))
        n = len(self.c) - 1
        for t in pts:
            horner = self.eval(t)
            direct = Q(0)
            for i, a in enumerate(self.c):
                direct = direct + Q(a) * t ** i
            if horner != direct:
                out.append(f"Горнер в точке {t} даёт {horner}, степенная сумма — {direct}")
        left, right = (self * self).derivative(), (self * self.derivative()).scale(2)
        if left != right:
            out.append(f"(P²)′ ≠ 2·P·P′: {left!r} против {right!r}")
        rec = self.reciprocal()
        for t in pts:
            if rec.eval(Q(1) / t) * t ** n != self.eval(t):
                out.append(f"зазеркалье: R(1/t)·t^{n} ≠ P(t) при t = {t}")
        cs = list(self.c)                             # та же обрезка, что внутри результанта
        while len(cs) > 1 and cs[-1] == 0:
            cs.pop()
        m = len(cs) - 1
        sign = -1 if m % 2 else 1                     # Res(P, y−a) = (−1)^deg P·P(a)
        for a in (2, -3):
            res = resultant_values(list(cs), [-a, 1])
            if res != sign * self.eval(a):
                out.append(f"результант Сильвестра–Барейса Res(P, y−{a}) = {res} "
                           f"≠ (−1)^{m}·P({a}) = {sign * self.eval(a)}")
        if n > 0 and n % 2 == 0 and self.is_palindromic():
            try:
                red = self.palindromic_reduce()
            except ValueError as e:
                out.append(f"палиндром не свёлся подстановкой Люка: {e}")
            else:
                h = n // 2
                for t in pts:
                    if red.eval(t + Q(1) / t) * t ** h != self.eval(t):
                        out.append(f"подстановка Люка: Q(t + 1/t)·t^{h} ≠ P(t) при t = {t}")
        return out


def _row_cleared(row):
    """Строка, домноженная на общий знаменатель: (целая строка, множитель,
    класс дроби или None). Целые проходят как есть (знаменатель 1)."""
    d = 1
    frac_cls = None
    for x in row:
        d = lcm(d, x.denominator)
        if frac_cls is None and not isinstance(x, int):
            frac_cls = type(x)  # дробь на входе — дробь и на выходе
    if d == 1:
        return [x.numerator for x in row], 1, frac_cls
    return [x.numerator * (d // x.denominator) for x in row], d, frac_cls


def _bareiss_det_int(a: "list[list[int]]") -> int:
    """Определитель целой матрицы алгоритмом Барейса: без дробей,
    каждое деление — нацело и без остатка. Матрица портится на месте."""
    n = len(a)
    if n == 0:
        return 1
    sign = 1
    prev = 1
    for k in range(n - 1):
        if a[k][k] == 0:
            for r in range(k + 1, n):
                if a[r][k] != 0:
                    a[k], a[r] = a[r], a[k]
                    sign = -sign
                    break
            else:
                return 0
        pivot = a[k][k]
        for i in range(k + 1, n):
            aik = a[i][k]
            row_i, row_k = a[i], a[k]
            for j in range(k + 1, n):
                # тождество Барейса: делитель prev делит нацело
                row_i[j] = (row_i[j] * pivot - aik * row_k[j]) // prev
            row_i[k] = 0
        prev = pivot
    return sign * a[n - 1][n - 1]


def _bareiss_det(m):
    """Определитель точно. Основной путь — целая матрица и целый Барейс;
    резерв для рациональных входов — домножение строк на общие
    знаменатели с учётом накопленного множителя в результате."""
    n = len(m)
    if n == 0:
        return 1
    rows = []
    mult = 1
    frac_cls = None
    for row in m:
        r, d, fc = _row_cleared(row)
        rows.append(r)
        mult *= d
        frac_cls = frac_cls or fc
    det = _bareiss_det_int(rows)
    if frac_cls is not None:
        return frac_cls(det, mult)
    return det if mult == 1 else det // mult


def resultant_values(p_y: list, q_y: list):
    """Результант двух многочленов от y (коэффициенты — рациональные на
    границе, млад→ст). Внутри: очистка знаменателей в целую матрицу
    Сильвестра, точный целочисленный Барейс, деление результата на
    накопленный множитель dp^n·dq^m. Дробь на входе — дробь на выходе,
    целые на входе — целое на выходе."""
    while len(p_y) > 1 and p_y[-1] == 0:
        p_y = p_y[:-1]
    while len(q_y) > 1 and q_y[-1] == 0:
        q_y = q_y[:-1]
    m, n = len(p_y) - 1, len(q_y) - 1
    if m == 0:
        return p_y[0] ** n
    if n == 0:
        return q_y[0] ** m
    pi, dp, fp = _row_cleared(p_y)
    qi, dq, fq = _row_cleared(q_y)
    frac_cls = fp or fq
    size = m + n
    s = [[0] * size for _ in range(size)]
    for r in range(n):                      # строки из P
        for i, a in enumerate(reversed(pi)):
            s[r][r + i] = a
    for r in range(m):                      # строки из Q
        for i, b in enumerate(reversed(qi)):
            s[n + r][r + i] = b
    det = _bareiss_det_int(s)
    mult = dp ** n * dq ** m
    if frac_cls is not None:
        return frac_cls(det, mult)
    return det if mult == 1 else det // mult


def interpolate_int_poly(points: list) -> IntPoly:
    """Интерполяция Лагранжа по точкам (t, v), t ∈ ℤ, с приведением к
    примитивному ℤ[x]. Значения v — рациональные на границе; внутри они
    очищаются общим знаменателем D, базисные веса — общим кратным L
    произведений M_i = Π(t_i−t_j), и собирается целый многочлен
    D·L·P(x); primitive() снимает оба множителя."""
    ts = [int(t) for t, _ in points]
    d_common = 1
    for _, v in points:
        d_common = lcm(d_common, v.denominator)
    ws = [v.numerator * (d_common // v.denominator) for _, v in points]
    prods = []
    for i, ti in enumerate(ts):
        m = 1
        for j, tj in enumerate(ts):
            if j != i:
                m *= ti - tj
        prods.append(m)
    l_common = 1
    for m in prods:
        l_common = lcm(l_common, abs(m))
    coeffs = [0] * len(points)
    for i, ti in enumerate(ts):
        w = ws[i] * (l_common // prods[i])
        if w == 0:
            continue
        # базисный многочлен Π_{j≠i}(x − t_j) — целый
        base = [1]
        for j, tj in enumerate(ts):
            if j == i:
                continue
            new = [0] * (len(base) + 1)
            for k, b in enumerate(base):
                new[k] += b * (-tj)
                new[k + 1] += b
            base = new
        for k, b in enumerate(base):
            coeffs[k] += b * w
    return IntPoly(coeffs).primitive()
