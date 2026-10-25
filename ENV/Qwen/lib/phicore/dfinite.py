"""Замыкание D-конечных рядов: сумма и свёртка Коши — одной линейкой класса 7.

Слово автора 03.09 («названные отказы подлежат явному пересмотру…»,
о строках ФIOS для 7/8/9): «постановка вопроса и отказ не предполагали
наличия других методов решения; мы не можем себе позволить закрыть 99 %
математики потому, что мы её не знали». Метод есть, и он точный:

Последовательность (aₙ) D-конечна (голономна), если удовлетворяет
линейной рекурренте с многочленными коэффициентами
    Σ_{i=0}^{r} pᵢ(n)·a_{n+i} = 0,  pᵢ ∈ ℤ[n], p_r ≠ 0,  n ≥ n0.
Производящий ряд f(x) = Σ aₙxⁿ тогда удовлетворяет линейному ОДУ с
многочленными коэффициентами, и обратно (Стэнли 1980, «Differentiably
finite power series»; алгоритмы — Salvy, Zimmermann 1994, gfun). Класс
замкнут относительно сложения, умножения (свёртки Коши коэффициентов)
и умножения на число, а рекуррента результата вычисляется точно —
исключением Гаусса над полем ℚ(x): производные f⁽ᵏ⁾ сводятся по ОДУ к
f, …, f⁽ʳ⁻¹⁾ (плюс многочлен из начальных членов), производные
произведения по Лейбницу лежат в пространстве размерности r·s + r + s
+ 1, зависимость среди r·s + r + s + 2 векторов существует по счёту
размерности и находится в конечное число шагов. Ни одной гипотезы:
теорема замыкания и арифметика ℚ(x).

Перевод рекуррента → ОДУ: θ = x·d/dx, θ(xᵐ) = m·xᵐ, поэтому
    Σᵢ x^{r−i}·pᵢ(θ − i) f(x) = P(x),
где многочлен P (степени < n0 + r) собирается из начальных членов;
θᵏ = Σⱼ S(k, j)·xʲ·Dʲ (числа Стирлинга второго рода). Перевод ОДУ →
рекуррента: коэффициент при xⁿ у xᵐDʲf есть (n−m+j)_j·a_{n−m+j}
(убывающий факториал); правая часть многочлен ⇒ рекуррента однородна
при n ≥ deg P + 1 сверх сдвига. Носитель хранит n0 и достаточно
начальных членов: где старший коэффициент p_r(n) обращается в ноль,
рекуррента a_{n+r} не определяет — член берётся из начальных.

Сумма (10.09, очередь п. 41 «интерполяция для больших порядков»): рекуррента
суммы — зависимость сдвигов c_{n+k} = a_{n+k} + b_{n+k} в базисе
(a_n…a_{n+r−1}, b_n…b_{n+s−1}), порядок ≤ r + s НА ЭТОЙ ДОРОГЕ (за порогом
_SHIFT_LIMIT сумма идёт через ОДУ, и там граница другая — rs + r + s: замер
12.09 дал 17 против 9 на одной паре); она ищется ТОЛЬКО ПО
ТОЧКАМ (_dependency_points): строки при целых n — числа, масштабы —
известные произведения старших коэффициентов, коэффициенты λ_i —
целые определители в D + 1 точках (ядро m × (m+1): Барейсс над ℤ, либо по
модулю 61-битных простых с восстановлением CRT и точной проверкой B·v = 0 —
_nullvector, порог по замеру перелома 11.09, очередь п. 41), интерполяция по
последовательным целым узлам конечными разностями целиком над ℤ (11.09; Ньютон
над Q стоил 900 с из 1435 на сумме трёх Лежандра K = 20); тождество
Σ λ_j·row_j = 0 доказано счётом степеней (равенство многочленов ограниченной
степени в D_ver + 1 точках, значения λ_j в добавочных узлах — продолжением
таблицы разностей) — ни одной операции над ℚ(n) и ни одной над Q в цикле.
Дыры старшего коэффициента результата (многочлен степени сотен с огромными
коэффициентами) не ищутся перебором делителей, а закрываются прямым счётом
в момент встречи (DFiniteSeq(…, direct)).

Каждое замыкание ПРОВЕРЯЕТСЯ прогоном: члены результата по новой
рекурренте сверяются с прямым счётом (сумма членов, свёртка Коши) на
_VERIFY номерах сверх начальных; расхождение — отказ, не тихий ответ.

Сертификат сетки (enclosure) — из строения: сумма и произведение
отрезков сомножителей (ряды абсолютно сходятся там, где у линеек есть
мажоранта; произведение рядов Коши сходится к произведению сумм —
Мертенс 1875). У голого D-конечного ряда без строения из
линеек мажоранта хвоста доказана для рекуррент ЛЮБОГО порядка r ≥ 1
(`geometric_tail`, `tail_enclosure`: порядок 1 — 11.09, порядок r ≥ 2 — 01.10,
индукцией по окну): знаки коэффициентов и неравенство
`ρ^r·|p_r(n)| − Σ_{i<r} ρ^i·|p_i(n)| ≥ 0` проверяются сдвигом Тейлора, точно в ℚ;
голый ряд получает сертификат подбором (n*, ρ) по сетке (`_auto_tail`).
Остался честный отказ, и это три разных утверждения, а не одно:
  ПРЕДМЕТ: сертификат голого ряда, у которого геометрическая мажоранта не
  доказана ни при одной паре (n*, ρ) сетки — сходимость не геометрическая
  (Σ 1/n²: отношение соседей → 1, ρ < 1 недоказуемо в принципе) либо знак
  коэффициента на [n*, ∞) не доказан сдвигом Тейлора (условие достаточное).
  ПРИЧИНА: геометрическая дорога мажорирует члены прогрессией C·ρⁿ, ρ < 1, и
  ряды с полиномиальным убыванием ей не по силам. Это отказ ДОМА (общая
  дорога не построена), а не теорема: мажоранта у них существует (Меццаробба
  2019 — границы через особенности), её просто не строили здесь.
  АДРЕС: отказ поимённо — `geometric_tail`, dfinite.py:1400 (порядок 0: хвост
  точно нулевой, мажоранта не нужна) и `enclosure` → `_auto_tail`,
  dfinite.py:1659 («геометрическая мажоранта голого ряда не доказана ни при
  одной паре (n*, ρ) из сетки»); сертификат — у строения из линеек и по явному
  запросу `tail_enclosure(ряд, N, n*, ρ)`. Числа адресов СВЕРЯЕТ приёмка.
  ЧЕМ МЕРИЛОСЬ БЫ, ценой: мажоранта Меццаробба берёт радиус доминирующей
  особенности — наименьший по модулю корень старшего коэффициента p_r, то
  есть изоляцию корней целочисленного многочлена степени deg p_r, и затем
  явную постоянную перехода. Цена в числах ЭТОГО дома: поиск целых корней уже
  отказывает на степени 90 с ~400-битными коэффициентами за 0,15 с (`_TRIAL`
  = 1 000 000, `_ROOT_SCAN` = 4000, стр. 116–118), а свёртка Коши здесь доходит
  до порядка 22 — та же изоляция корней, что и у поиска дыр старшего
  коэффициента, и оплачивается ею же. Охват замерен `[вычислено здесь]`:
  сертификат есть у порядков 1, 2, 3, где условие доказано (2⁻ⁿ, 1/n!,
  cosh 1 = Σ 1/(2k)!, Σ 1/(3k)!), отказ — у расходящихся |a_n| ≡ 1 любого
  порядка (приёмка двумя дорогами). И «не доказано» ≠ «ложно»: условие
  `_nonneg_from` достаточное, не необходимое — для 1/n! та же ρ = 1/2 НЕ
  доказана при n* = 0 и доказана при n* = 1, ρ = 1/3 — при n* = 2, не при 0.

Ни одного float; Q, IntPoly, RatFunc. Обращение ряда и деление ряда
на ряд остаются отказом: 1/f D-конечна лишь в частных случаях, общего
замыкания нет (это теорема, не незнание).
"""
from __future__ import annotations

from math import gcd, isqrt, lcm

from .enclosure import Interval
from .errors import Unresolved
from .intpoly import IntPoly
from .ratfunc import RatFunc, _add, _divmod, _gcd, _mul, _scale, _strip
from .rational import Q
from .zero_state import Zero

_NUM = (int, Q)
#: членов сверки замыкания с прямым счётом сверх начальных
_VERIFY = 24
#: Порог замыкания суммы сдвигами: исключение над ℚ(n) в чистом Python при r + s > 8
#: медленнее пути через ОДУ (замер 05.09: 7 + 7 сдвигами не завершилось за 12 мин,
#: через ОДУ — порядок 22 за 170 с); выше порога — ОДУ, названо в очереди п. 41.
_SHIFT_LIMIT = 8      # замер 10.09: 7 + 7 по точкам — 15 с, но порог 14 удлинил приёмку dfinite сверх 900 с (случаи с большими степенями) — порог 8 по замеру; выше — ОДУ
#: порог перебора целых корней старшего коэффициента по оценке Коши
_ROOT_SCAN = 4000
#: предел пробного деления при поиске делителей свободного члена
_TRIAL = 1_000_000


def _rec_of(x) -> str:
    r = getattr(x, "record", None)
    if callable(r):
        return r()
    if isinstance(x, _NUM):
        return f"Φ(C=1 RATIONAL, {Q(x)})"
    return repr(x)


# ---------------------------------------------------------------------------
# многочлены над ℚ (кортежи Q по возрастанию степени; () — ноль)
# ---------------------------------------------------------------------------

def _qpoly(p: IntPoly):
    return _strip(Q(c) for c in p.c)


def _poly_shift(p, s):
    """p(x + s) точно: схема Горнера по множителю (x + s)."""
    s = Q(s)
    acc = ()
    for c in reversed(p):
        acc = _add(_mul(acc, (s, Q(1))), (c,))
    return _strip(acc)


def _poly_deriv(p):
    return _strip(Q(i) * c for i, c in enumerate(p) if i > 0)


def _poly_lcm(a, b):
    if not a or not b:
        return ()
    g = _gcd(a, b)
    q, r = _divmod(_mul(a, b), g)
    assert not r
    return q


def _poly_eval(p, x):
    acc = Q(0)
    for c in reversed(p):
        acc = acc * x + c
    return acc


def _rf(p, d=None) -> RatFunc:
    return RatFunc(p, d if d is not None else (Q(1),))


def _rf_deriv(f: RatFunc) -> RatFunc:
    """(N/D)' = (N'D − ND')/D²."""
    n, d = f.num, f.den
    return RatFunc(_add(_mul(_poly_deriv(n), d), tuple(-c for c in _mul(n, _poly_deriv(d)))),
                   _mul(d, d))


def _stirling2(kmax: int):
    """S(k, j), 0 ≤ j ≤ k ≤ kmax."""
    S = [[0] * (kmax + 1) for _ in range(kmax + 1)]
    S[0][0] = 1
    for k in range(1, kmax + 1):
        for j in range(1, k + 1):
            S[k][j] = j * S[k - 1][j] + S[k - 1][j - 1]
    return S


def _falling(j: int, shift: int):
    """(n + shift)_j = Π_{t=0}^{j−1} (n + shift − t) как многочлен от n."""
    acc = (Q(1),)
    for t in range(j):
        acc = _mul(acc, (Q(shift - t), Q(1)))
    return acc


def _int_roots(p: IntPoly) -> list:
    """Неотрицательные целые корни p ∈ ℤ[n], точно; при слишком большом
    свободном члене за пределом перебора — отказ, не пропуск."""
    if p.is_zero():
        raise ValueError("нулевой многочлен: корни — все n")
    c = list(p.c)
    roots = []
    if c[0] == 0:
        roots.append(0)
        while c[0] == 0:
            c.pop(0)
    if len(c) == 1:
        return roots
    bound = 1 + max(abs(x) for x in c[:-1]) // abs(c[-1]) + 1
    q = IntPoly(c)
    if bound <= _ROOT_SCAN:
        roots += [n for n in range(1, bound + 1) if q.eval(n) == 0]
        return sorted(set(roots))
    c0 = abs(c[0])
    divs = {1}
    m, d = c0, 2
    while d * d <= m and d <= _TRIAL:
        while m % d == 0:
            divs = {x * d ** k for x in divs for k in range(0, 40) if x * d ** k <= c0}
            m //= d
        d += 1
    if m > 1 and d > _TRIAL:
        raise Unresolved("int_roots", repr(p),
                         hint=f"целые корни старшего коэффициента не отделены: свободный член "
                              f"{c0} не разложен пробным делением до {_TRIAL}")
    if m > 1:
        divs |= {x * m for x in divs}
    roots += [n for n in sorted(divs) if n <= bound and q.eval(n) == 0]
    return sorted(set(roots))


# ---------------------------------------------------------------------------
# рекуррента и последовательность
# ---------------------------------------------------------------------------

class Recurrence:
    """Σ_{i=0}^{r} pᵢ(n)·a_{n+i} = 0 при n ≥ n0; pᵢ ∈ ℤ[n], p_r ≠ 0."""

    __slots__ = ("p", "n0")

    def __init__(self, p, n0: int = 0):
        p = [q if isinstance(q, IntPoly) else IntPoly(q) for q in p]
        while len(p) > 1 and p[-1].is_zero():
            p.pop()
        if not p or p[-1].is_zero():
            raise ValueError("рекуррента без старшего коэффициента")
        if not isinstance(n0, int) or isinstance(n0, bool) or n0 < 0:
            raise ValueError("n0 — целое ≥ 0")
        g = 0
        for q in p:
            g = gcd(g, q.content()) if not q.is_zero() else g
        self.p = tuple(IntPoly([c // g for c in q.c]) for q in p) if g > 1 else tuple(p)
        self.n0 = n0

    @property
    def order(self) -> int:
        return len(self.p) - 1

    def record(self) -> str:
        return "[" + "; ".join("[" + ",".join(map(str, q.c)) + "]" for q in self.p) + f"] n0={self.n0}"

    __repr__ = record

    def __eq__(self, other):
        return isinstance(other, Recurrence) and self.p == other.p and self.n0 == other.n0

    def __hash__(self):
        return hash(("rec", tuple(tuple(q.c) for q in self.p), self.n0))


class DFiniteSeq:
    """D-конечная последовательность: рекуррента + начальные члены (Q)."""

    __slots__ = ("rec", "_a", "_holes", "_direct")

    def __init__(self, rec: Recurrence, initial, direct=None):
        """direct — прямой счёт члена (у замыканий: сумма/свёртка членов операндов);
        с ним дыры старшего коэффициента не ищутся заранее (целые корни многочлена
        степени сотни с огромными коэффициентами — перебор делителей, 117 с на 7 + 7,
        замер 10.09), а закрываются прямым счётом в момент встречи."""
        if not isinstance(rec, Recurrence):
            raise TypeError("rec — только Recurrence")
        a = []
        for v in initial:
            if isinstance(v, float) or isinstance(v, bool):
                raise TypeError("float запрещён построением")
            a.append(Q(v))
        r = rec.order
        holes = [] if direct is not None else [k for k in _int_roots(rec.p[-1]) if k >= rec.n0]
        need = rec.n0 + r
        if holes:
            need = max(need, holes[-1] + r + 1)
        if len(a) < need:
            raise ValueError(f"начальных членов {len(a)}, нужно ≥ {need} (порядок {r}, n0={rec.n0}, "
                             f"нули старшего коэффициента {holes})")
        self.rec, self._a, self._holes, self._direct = rec, a, set(holes), direct
        # начальные члены обязаны согласоваться с рекуррентой там, где она определена
        for k in range(rec.n0, len(a) - r):
            if k in self._holes:
                continue
            s = sum((rec.p[i].eval(k) * a[k + i] for i in range(r + 1)), Q(0))
            if s != 0:
                raise ValueError(f"начальные члены противоречат рекурренте при n={k}: невязка {s}")

    def _extend(self, n: int) -> None:
        r, p, a = self.rec.order, self.rec.p, self._a
        while len(a) <= n:
            k = len(a) - r
            lead = p[-1].eval(k)
            if lead == 0:
                if self._direct is None:
                    raise Unresolved("term", self.record(),
                                     hint=f"p_r({k}) = 0 — член a_{k + r} рекуррентой не определён")
                a.append(Q(self._direct(k + r)))     # дыра закрыта прямым счётом
                continue
            s = sum((p[i].eval(k) * a[k + i] for i in range(r)), Q(0))
            a.append(-s / lead)

    def term(self, n: int) -> Q:
        if isinstance(n, bool) or not isinstance(n, int):
            raise TypeError("номер — целое")
        if n < 0:
            raise ValueError("номер ≥ 0")
        self._extend(n)
        return self._a[n]

    @property
    def initial(self):
        return tuple(self._a[:self.rec.n0 + self.rec.order + (max(self._holes) + 1 if self._holes else 0)])

    def record(self) -> str:
        return f"rec={self.rec.record()}, a=[{','.join(str(v) for v in self._a[:self.rec.n0 + self.rec.order])}]"

    __repr__ = record


# ---------------------------------------------------------------------------
# рекуррента ↔ ОДУ (коэффициенты — многочлены над ℚ, правая часть — многочлен)
# ---------------------------------------------------------------------------

def _rec_to_ode(seq: DFiniteSeq):
    """(qs, P): Σⱼ qⱼ(x)·f⁽ʲ⁾ = P(x) для f = Σ aₙxⁿ."""
    rec, r = seq.rec, seq.rec.order
    kmax = max(q.degree for q in rec.p)
    S = _stirling2(kmax)
    qs = [() for _ in range(kmax + 1)]
    for i, pi in enumerate(rec.p):
        ct = _poly_shift(_qpoly(pi), -i)            # pᵢ(t − i) как многочлен от t = θ
        for k, ck in enumerate(ct):
            if ck == 0:
                continue
            for j in range(0, k + 1):
                s_kj = S[k][j]
                if s_kj == 0:
                    continue
                mono = (Q(0),) * (j + r - i) + (ck * s_kj,)
                qs[j] = _add(qs[j], mono)
    while len(qs) > 1 and not qs[-1]:
        qs.pop()
    # правая часть: x^r · Σᵢ x^{−i} Σ_{m<n0+i} pᵢ(m−i)·aₘ·xᵐ
    P = ()
    for i, pi in enumerate(rec.p):
        for m in range(0, rec.n0 + i):
            c = Q(pi.eval(m - i)) * seq.term(m)
            if c != 0:
                P = _add(P, (Q(0),) * (m - i + r) + (c,))
    return [tuple(q) for q in qs], P


def _derivs(qs, P, K: int):
    """Аффинные сведения f⁽ᵏ⁾ = Σ_{i<s} αₖ,ᵢ f⁽ⁱ⁾ + βₖ (RatFunc), k = 0..K."""
    s = len(qs) - 1
    lead = _rf(qs[s])
    red_alpha = [-(_rf(qs[j]) / lead) for j in range(s)]      # f⁽ˢ⁾ через младшие
    red_beta = _rf(P) / lead
    rows = []
    for k in range(K + 1):
        if k < s:
            alpha = [_rf((Q(1),)) if i == k else _rf(()) for i in range(s)]
            beta = _rf(())
        elif k == s:
            alpha, beta = list(red_alpha), red_beta
        else:
            pa, pb = rows[-1]
            alpha = [_rf_deriv(pa[i]) for i in range(s)]
            for i in range(1, s):
                alpha[i] = alpha[i] + pa[i - 1]
            top = pa[s - 1] if s else None
            if top is not None:
                for i in range(s):
                    alpha[i] = alpha[i] + top * red_alpha[i]
                beta = _rf_deriv(pb) + top * red_beta
            else:
                beta = _rf_deriv(pb)
        rows.append((alpha, beta))
    return rows


def _dependency(vectors):
    """Первая линейная зависимость среди векторов над ℚ(x) (инкрементный
    ступенчатый вид); возвращает список коэффициентов λ₀..λₖ, Σλᵢvᵢ = 0."""
    basis = []   # (pivot, row, combo)
    zero = _rf(())
    for k, v in enumerate(vectors):
        v = list(v)
        combo = [zero] * k + [_rf((Q(1),))]
        for piv, row, cmb in basis:
            if not v[piv].is_zero:
                f = v[piv] / row[piv]
                v = [v[i] - f * row[i] for i in range(len(v))]
                combo = [combo[i] - f * cmb[i] if i < len(cmb) else combo[i] for i in range(len(combo))]
        nz = [i for i, x in enumerate(v) if not x.is_zero]
        if not nz:
            return combo
        basis.append((nz[0], v, combo))
    raise RuntimeError("зависимость не найдена — счёт размерности нарушен")


def _combine_odes(ode_f, ode_g, kind: str):
    """ОДУ для f+g ('+') или f·g ('×') по ОДУ сомножителей: (qs, P)."""
    qf, Pf = ode_f
    qg, Pg = ode_g
    r, s = len(qf) - 1, len(qg) - 1
    if kind == "+":
        D = r + s
        K = D
        df, dg = _derivs(qf, Pf, K), _derivs(qg, Pg, K)
        vecs, consts = [], []
        for k in range(K + 1):
            af, bf = df[k]
            ag, bg = dg[k]
            vecs.append(list(af) + list(ag))
            consts.append(bf + bg)
    else:
        D = r * s + r + s
        K = D
        df, dg = _derivs(qf, Pf, K), _derivs(qg, Pg, K)
        binom = [[1] * (K + 1) for _ in range(K + 1)]
        for n in range(K + 1):
            for m in range(1, n):
                binom[n][m] = binom[n - 1][m - 1] + binom[n - 1][m]
        vecs, consts = [], []
        for k in range(K + 1):
            vec = [_rf(()) for _ in range(D)]
            const = _rf(())
            for i in range(k + 1):
                c = Q(binom[k][i])
                af, bf = df[i]
                ag, bg = dg[k - i]
                for a_i in range(r):
                    if af[a_i].is_zero:
                        continue
                    for b_j in range(s):
                        if not ag[b_j].is_zero:
                            vec[a_i * s + b_j] = vec[a_i * s + b_j] + af[a_i] * ag[b_j] * c
                    if not bg.is_zero:
                        vec[r * s + a_i] = vec[r * s + a_i] + af[a_i] * bg * c
                if not bf.is_zero:
                    for b_j in range(s):
                        if not ag[b_j].is_zero:
                            vec[r * s + r + b_j] = vec[r * s + r + b_j] + bf * ag[b_j] * c
                    if not bg.is_zero:
                        const = const + bf * bg * c
            vecs.append(vec)
            consts.append(const)
    lam = _dependency(vecs)
    rhs = _rf(())
    for l, c in zip(lam, consts):
        rhs = rhs + l * c
    # очистка знаменателей: общий знаменатель всех λ и правой части
    L = (Q(1),)
    for f in list(lam) + [rhs]:
        L = _poly_lcm(L, f.den)
    qs = [_divmod(_mul(f.num, L), f.den)[0] for f in lam]
    P = _divmod(_mul(rhs.num, L), rhs.den)[0]
    while len(qs) > 1 and not qs[-1]:
        qs.pop()
    return qs, P


def _ode_to_rec(qs, P) -> Recurrence:
    """Рекуррента коэффициентов по ОДУ Σ qⱼ f⁽ʲ⁾ = P."""
    shifts = {}
    for j, q in enumerate(qs):
        for m, c in enumerate(q):
            if c == 0:
                continue
            s = j - m
            shifts[s] = _add(shifts.get(s, ()), _scale(_falling(j, j - m), c))
    shifts = {s: p for s, p in shifts.items() if p}
    if not shifts:
        raise ValueError("нулевой оператор")
    s_min, s_max = min(shifts), max(shifts)
    R = s_max - s_min
    polys = [_poly_shift(shifts.get(s_min + i, ()), -s_min) if shifts.get(s_min + i) else () for i in range(R + 1)]
    den = 1
    for p in polys:
        for c in p:
            den = lcm(den, c.denominator)
    ipolys = [IntPoly([int(c * den) for c in p] or [0]) for p in polys]
    degP = len(P) - 1 if P else -1
    n0 = max(0, degP + s_min + 1, s_min)
    return Recurrence(ipolys, n0)


def _series_div(N, D, K: int):
    """Коэффициенты f₀..f_K ряда N(s)/D(s) при D₀ ≠ 0 — точно (Q)."""
    N = list(N) + [Q(0)] * (K + 1)
    D = list(D) + [Q(0)] * (K + 1)
    if D[0] == 0:
        raise ValueError("D(0) = 0: ряд частного не определён")
    f = []
    for j in range(K + 1):
        acc = N[j] - sum((D[i] * f[j - i] for i in range(1, j + 1)), Q(0))
        f.append(acc / D[0])
    return f


def _rescale_rec(rec: Recurrence, h) -> Recurrence:
    """u_k = c_k·h^k: из Σ pᵢ(k)·c_{k+i} = 0 следует Σ pᵢ(k)·h^{r−i}·u_{k+i} = 0
    (умножено на h^{k+r}); знаменатели степеней h очищены до ℤ[k]."""
    h = Q(h)
    r = rec.order
    polys = [_scale(_qpoly(pi), h ** (r - i)) for i, pi in enumerate(rec.p)]
    den = 1
    for q in polys:
        for c in q:
            den = lcm(den, c.denominator)
    return Recurrence([IntPoly([int(c * den) for c in q] or [0]) for q in polys], rec.n0)


#: предел деления отрезка при поиске кусков со сходящимся разложением
_SPLIT_DEPTH = 60


def _integral_pieces(num, den, a, b, depth: int = 0):
    """Куски [l; r] отрезка, на каждом из которых разложение первообразной
    вокруг середины m сходится на концах: D(s) = den(m + s) не имеет корней
    в круге |s| ≤ h (h — полуширина): достаточно |D₀| > Σ_{i≥1} |Dᵢ|·hⁱ —
    точная оценка на Q; иначе отрезок делится пополам."""
    m, h = (a + b) / 2, (b - a) / 2
    D = _poly_shift(_qpoly(den), m)
    if abs(D[0]) > sum((abs(c) * h ** i for i, c in enumerate(D) if i > 0), Q(0)):
        return [(a, b)]
    if depth >= _SPLIT_DEPTH:
        raise Unresolved("from_integral", f"∫[{a}; {b}]",
                         hint=f"полюс подынтегральной слишком близко к отрезку: за {_SPLIT_DEPTH} "
                              f"делений сходящееся разложение не найдено")
    return _integral_pieces(num, den, a, m, depth + 1) + _integral_pieces(num, den, m, b, depth + 1)


def _piece_series(num, den, l, r) -> "DFiniteSeries":
    """∫_l^r num/den как D-конечный ряд: F(m+s) = Σ c_k s^k, c₀ = 0,
    D(s)·F′(s) = N(s) даёт рекурренту c_k из ОДУ; значение куска
    F(m+h) − F(m−h) = Σ c_k(h^k − (−h)^k) — разность двух перемасштабирований."""
    m, h = (l + r) / 2, (r - l) / 2
    N = _poly_shift(_qpoly(num), m)
    D = _poly_shift(_qpoly(den), m)
    rec = _ode_to_rec([(), D], N)                       # D·F′ = N
    holes = [k for k in _int_roots(rec.p[-1]) if k >= rec.n0]
    need = rec.n0 + rec.order + (holes[-1] + 1 if holes else 0) + _VERIFY
    f = _series_div(N, D, need)
    c = [Q(0)] + [f[k - 1] / k for k in range(1, need + 1)]
    u = DFiniteSeries(DFiniteSeq(_rescale_rec(rec, h), [ck * h ** k for k, ck in enumerate(c)]))
    v = DFiniteSeries(DFiniteSeq(_rescale_rec(rec, -h), [ck * (-h) ** k for k, ck in enumerate(c)]))
    return u - v


def _shift_vectors(seq: DFiniteSeq, K: int):
    """a_{n+k}, k = 0…K, в базисе a_n … a_{n+r−1} над ℚ(n): единичные векторы до
    r − 1, дальше — рекуррентой при номере n + k − r (сдвиг многочленов)."""
    rec, r = seq.rec, seq.rec.order
    ps = [_qpoly(q) for q in rec.p]
    zero, one = _rf(()), _rf((Q(1),))
    vecs = [[one if i == k else zero for i in range(r)] for k in range(r)]
    for k in range(r, K + 1):
        sh = k - r
        lead = _rf(_poly_shift(ps[r], sh))
        v = [zero] * r
        for i in range(r):
            c = _rf(_poly_shift(ps[i], sh)) / lead
            if c.is_zero:
                continue
            u = vecs[k - r + i]
            v = [v[j] - c * u[j] for j in range(r)]
        vecs.append(v)
    return vecs


def _poly_gcd_all(polys):
    g = ()
    for q in polys:
        if q:
            g = _gcd(g, q) if g else q
    return g


def _poly_lcm_list(polys):
    L = (Q(1),)
    for q in polys:
        L = _poly_lcm(L, q)
    return L


def _exact_div(u, v):
    q, r = _divmod(u, v)
    if r:
        raise ArithmeticError("деление многочленов не точное — инвариант исключения нарушен")
    return q


#: Исключение БЕЗ ДРОБЕЙ над ℚ(x) (заготовка `_dependency_ff`) СНЯТО 12.09.
#: Что обещало: ни одной дроби ℚ(x) внутри — строки ступенчатого вида держатся
#: многочленами, общий множитель снимается НОД. Чем кончилось замером (05.09):
#: на приёмке dfinite (суммы порядков ≤ 8) набор НЕ ЗАВЕРШИЛСЯ за 600 с против
#: ≈ 60 с у исключения над ℚ(x) — рост степеней в строках без полного снятия
#: множителей. Снято 39 строк (38 непустых; 41 — это расстояние до следующего
#: `def` вместе с двумя пустыми-разделителями, и «кодом» оно не было: замер
#: 17.09 поправил названное здесь число), вхождений имени было 2, вызовов 0
#: (глава 11 §20 п. 4, замер 12.09):
#: мёртвая дорога снята, а её замер оставлен здесь числом. Работающие дороги
#: зависимости — `_dependency` (ℚ(x)), `_dependency_interp` (интерполяция) и
#: `_dependency_points_multi` (только по точкам).


def _int_det(M):
    """Определитель целой матрицы исключением Барейсса — без дробей."""
    A = [list(r) for r in M]
    n = len(A)
    if n == 0:
        return 1
    sign, prev = 1, 1
    for k in range(n - 1):
        if A[k][k] == 0:
            sw = next((i for i in range(k + 1, n) if A[i][k] != 0), None)
            if sw is None:
                return 0
            A[k], A[sw] = A[sw], A[k]
            sign = -sign
        akk = A[k][k]
        for i in range(k + 1, n):
            aik = A[i][k]
            row_i, row_k = A[i], A[k]
            for j in range(k + 1, n):
                row_i[j] = (row_i[j] * akk - aik * row_k[j]) // prev
            row_i[k] = 0
        prev = akk
    return sign * A[n - 1][n - 1]


def _solve_ff(A, b):
    """Барейсс на [A | b] (A — m × m целая) и обратная подстановка без дробей:
    (det(A)·x, det(A)) для A·x = −b; None при det(A) = 0. Целые по Крамеру,
    деление на диагональ точное — ни одного НОД над огромными числами."""
    m = len(A)
    M = [list(r) + [bi] for r, bi in zip(A, b)]
    sign, prev = 1, 1
    for k in range(m - 1):
        if M[k][k] == 0:
            sw = next((i for i in range(k + 1, m) if M[i][k] != 0), None)
            if sw is None:
                return None
            M[k], M[sw] = M[sw], M[k]
            sign = -sign
        akk = M[k][k]
        for i in range(k + 1, m):
            aik = M[i][k]
            row_i, row_k = M[i], M[k]
            if aik == 0:                       # строка лишь масштабируется — деление точное
                for j in range(k + 1, m + 1):
                    row_i[j] = row_i[j] * akk // prev
            else:
                for j in range(k + 1, m + 1):
                    row_i[j] = (row_i[j] * akk - aik * row_k[j]) // prev
            row_i[k] = 0
        prev = akk
    det = sign * M[m - 1][m - 1]
    if det == 0:
        return None
    y = [0] * m
    for i in range(m - 1, -1, -1):
        acc = -det * M[i][m]
        for j in range(i + 1, m):
            acc -= M[i][j] * y[j]
        q, r = divmod(acc, M[i][i])
        if r:
            raise ArithmeticError("обратная подстановка без дробей: деление не точное — инвариант Барейсса нарушен")
        y[i] = q
    return y, det


def _nullvector_minors(B):
    """Канонический целый вектор ядра матрицы m × (m+1): v_i = (−1)^i·det(B без
    столбца i) — одним исключением Барейсса вместо m + 1 определителей (очередь
    п. 41, 10.09). Исключается столбец p с det ≠ 0 (сначала последний); ядро
    одномерно, и решение с исключённым столбцом p равно (−1)^p·v. При всех
    минорах нулевых (ранг < m в этой точке — законный нуль многочленов миноров)
    — нулевой вектор."""
    m = len(B)
    for p in range(m, -1, -1):
        A = [[row[j] for j in range(m + 1) if j != p] for row in B]
        b = [row[p] for row in B]
        res = _solve_ff(A, b)
        if res is None:
            continue
        y, det = res
        w = y[:p] + [det] + y[p:]
        return [(-x if p % 2 else x) for x in w]
    return [0] * (m + 1)


#: основания Миллера–Рабина, при которых тест ДЕТЕРМИНИРОВАН для n < 3,3·10²⁴
#: (Sorenson, Webster 2015, «Strong pseudoprimes to twelve prime bases») [внешний источник];
#: 61-битные модули лежат далеко внутри
_MR_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)


def _is_prime_det(n: int) -> bool:
    """Детерминированный Миллер–Рабин для n < 3,3·10²⁴ (основания _MR_BASES)."""
    if n < 2:
        return False
    for p in _MR_BASES:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in _MR_BASES:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


#: машинные простые счёта по модулю: 61-битные, вниз от 2⁶¹ − 1 (простое Мерсенна);
#: список растёт по требованию (_prime61) и держится в памяти модуля
_PRIMES61 = [(1 << 61) - 1]


def _prime61(i: int) -> int:
    """i-е (с нуля) 61-битное простое вниз от 2⁶¹ − 1 — по требованию, с памятью."""
    while len(_PRIMES61) <= i:
        c = _PRIMES61[-1] - 2
        while not _is_prime_det(c):
            c -= 2
        _PRIMES61.append(c)
    return _PRIMES61[i]


def _nullvector_mod(B, p: int):
    """Вектор миноров (тот же, что у _nullvector_minors) ПО МОДУЛЮ простого p: Гаусс над
    ℤ/p с тем же выбором исключаемого столбца (первый с det ≢ 0 mod p от конца),
    v_i ≡ (−1)^i·det(B без столбца i); все миноры ≡ 0 — нулевой вектор (законный
    вычет). Значения в [0, p). Одно исключение стоит O(m³) операций интерпретатора."""
    m = len(B)
    for pc in range(m, -1, -1):
        M = [[row[j] % p for j in range(m + 1) if j != pc] + [row[pc] % p] for row in B]
        sign = 1
        singular = False
        for k in range(m):
            if M[k][k] == 0:
                sw = next((i for i in range(k + 1, m) if M[i][k]), None)
                if sw is None:
                    singular = True
                    break
                M[k], M[sw] = M[sw], M[k]
                sign = -sign
            row_k = M[k]
            inv = pow(row_k[k], -1, p)
            for i in range(k + 1, m):
                f = M[i][k] * inv % p
                if f:
                    row_i = M[i]
                    for j in range(k, m + 1):
                        row_i[j] = (row_i[j] - f * row_k[j]) % p
        if singular:
            continue
        det = sign % p
        for k in range(m):
            det = det * M[k][k] % p
        # A·x = −b, y = det·x — числители Крамера по модулю p (как в _solve_ff)
        y = [0] * m
        for i in range(m - 1, -1, -1):
            acc = -det * M[i][m]
            for j in range(i + 1, m):
                acc -= M[i][j] * y[j]
            y[i] = acc * pow(M[i][i], -1, p) % p
        w = y[:pc] + [det] + y[pc:]
        return [(-x) % p for x in w] if pc % 2 else w
    return [0] * (m + 1)


def _hadamard_bits(B) -> int:
    """Оценка Адамара: любой минор матрицы B по модулю < 2^bits (произведение длин строк)."""
    bits = 0
    for row in B:
        s = sum(x * x for x in row)
        bits += (s.bit_length() + 1) // 2
    return bits


class _CRTNotConverged(ArithmeticError):
    """Названный отказ: ядро по модулю простых не сошлось в отведённое число простых."""


def _nullvector_crt(B, max_primes=None, primitive: bool = False):
    """Вектор миноров m × (m+1) (тот же, что у _nullvector_minors) счётом ПО МОДУЛЮ
    61-битных простых с восстановлением по КТО (шаг Гарнера, симметричный
    представитель в (−M/2, M/2)); простые добавляются до сходимости: представитель
    не изменился при добавлении простого И B·v = 0 ТОЧНО над ℤ — либо модуль превысил
    удвоенную оценку Адамара (тогда представитель точен по оценке, проверка B·v = 0
    остаётся). primitive — снять целое содержание (НОД компонент; знак сохранён).
    Несходимость за max_primes (по умолчанию — до границы Адамара плюс запас) —
    _CRTNotConverged, названный отказ; запасной путь у вызывающего — Барейсс над ℤ.

    Замер 11.09 [вычислено здесь] (Лежандр 7 + 7, m = 14, элементы ≈ 700 бит, миноры
    ≈ 3000 бит): Гаусс mod p 0,17 мс × 52 простых = 9 мс на точку против Барейсса
    0,6 мс — в чистом Python интерпретируемый цикл по модулю проигрывает арифметике
    больших целых в C; перелом — см. _CRT_MIN_BITS."""
    m = len(B)
    if m == 0:
        return [1]
    hb = _hadamard_bits(B)
    limit = hb // 60 + 3 if max_primes is None else max_primes
    M = 1
    v = [0] * (m + 1)
    prev = None
    for i in range(limit):
        p = _prime61(i)
        vp = _nullvector_mod(B, p)
        inv = pow(M % p, -1, p)
        v = [vi + M * (((vpi - vi) * inv) % p) for vi, vpi in zip(v, vp)]
        M *= p
        half = M // 2
        sym = [x - M if x > half else x for x in v]
        by_bound = M.bit_length() > hb + 2
        if by_bound or sym == prev:
            if all(sum(r[j] * sym[j] for j in range(m + 1)) == 0 for r in B):
                if primitive:
                    g = gcd(*sym)
                    if g > 1:
                        sym = [x // g for x in sym]
                return sym
            if by_bound:
                raise ArithmeticError("ядро по модулю: представитель точен по Адамару, но B·v ≠ 0 — инвариант нарушен")
        prev = sym
    raise _CRTNotConverged(f"ядро по модулю простых не сошлось за {limit} простых "
                           f"(модуль {M.bit_length()} бит, оценка Адамара {hb} бит)")


#: порог (биты наибольшего элемента), от которого ядро считается по модулю простых с CRT,
#: ниже — Барейсс над ℤ. Замер перелома 11.09 [вычислено здесь], плотные случайные
#: матрицы, отношение времени CRT/Барейсс: m = 14 — 3,11 (700 бит), 1,71 (2500),
#: 1,23 (10 000), 1,03 (40 000); m = 20 — 0,86 (2500), 0,67 (10 000). На СТРУКТУРНЫХ
#: матрицах замыкания (единичные столбцы младших сдвигов, содержание миноров 2/3 битов)
#: Барейсс на два порядка быстрее плотного случая тех же битов (Лежандр 7 + 7: 0,6 мс
#: против 60 мс), а CRT структуры не использует (K = 20: 0,03–0,05 с против ≈ 0,5 с) —
#: порог консервативный, по перелому m = 14; реальные матрицы лежат ниже него.
_CRT_MIN_BITS = 40_000


def _nullvector(B):
    """Ядро m × (m+1): CRT при элементах от _CRT_MIN_BITS битов, иначе Барейсс над ℤ;
    несходимость CRT — Барейсс (запасной путь). Оба пути точны (B·v = 0 над ℤ)."""
    if max((abs(x).bit_length() for r in B for x in r), default=0) >= _CRT_MIN_BITS:
        try:
            return _nullvector_crt(B)
        except _CRTNotConverged:
            pass
    return _nullvector_minors(B)


def _rank_q(rows):
    """Ранг семейства числовых строк над ℚ и столбцы опор (порядок строк сохранён)."""
    A = [list(r) for r in rows]
    pivots, r = [], 0
    ncol = len(A[0]) if A else 0
    for c in range(ncol):
        piv = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c] / A[r][c]
                A[i] = [A[i][j] - f * A[r][j] for j in range(ncol)]
        pivots.append(c)
        r += 1
        if r == len(A):
            break
    return r, pivots


def _newton_interp(xs, ys):
    """Многочлен степени ≤ len(xs) − 1 через точки (xs, ys) над ℚ — разделённые разности."""
    n = len(xs)
    coef = [Q(y) for y in ys]
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            coef[i] = (coef[i] - coef[i - 1]) / Q(xs[i] - xs[i - j])
    poly = ()
    for i in range(n - 1, -1, -1):
        poly = _add(_mul(poly, (Q(-xs[i]), Q(1))), (coef[i],))
    return _strip(poly)


def _fwd_diffs(ys):
    """Таблица конечных разностей целых значений в последовательных целых узлах:
    (Δ^k y_0)_k — верхняя кромка, (Δ^k y_{n−1−k})_k — нижняя; O(n²) вычитаний."""
    a = list(ys)
    top, last = [a[0]], [a[-1]]
    while len(a) > 1:
        a = [y - x for x, y in zip(a, a[1:])]
        top.append(a[0])
        last.append(a[-1])
    return top, last


def _extend_diffs(last, count: int):
    """Значения многочлена степени ≤ len(last) − 1 в count СЛЕДУЮЩИХ целых узлах по
    нижней кромке таблицы разностей (Δ^{len} ≡ 0): O(deg) сложений на узел — вместо
    Горнера над Q (11.09: проверка тождества на K = 20 стоила 405 с в _poly_eval)."""
    d = list(last)
    out = []
    for _ in range(count):
        for k in range(len(d) - 2, -1, -1):
            d[k] += d[k + 1]
        out.append(d[0])
    return out


def _newton_interp_int(x0: int, ys, top=None):
    """Интерполяция по ПОСЛЕДОВАТЕЛЬНЫМ целым узлам x0, x0 + 1, … целых значений —
    целиком над ℤ: P(x) = Σ_k Δ^k y_0·C(x − x0, k), значит D!·P ∈ ℤ[x]; схема Горнера
    в базисе падающих факториалов, ни одного НОД внутри цикла (11.09: разделённые
    разности над Q на D = 1320 стоили 900 с). top — готовая верхняя кромка таблицы
    разностей (_fwd_diffs), если она уже посчитана для проверки тождества. Возвращает
    (целые коэффициенты N, знаменатель den > 0): P = N/den, содержание снято; нулевой
    P — ((), 1)."""
    if top is None:
        top, _ = _fwd_diffs(ys)
    D = len(top) - 1
    while D >= 0 and top[D] == 0:
        D -= 1
    if D < 0:
        return (), 1
    fact = [1] * (D + 1)                      # D!/k!
    for k in range(D - 1, -1, -1):
        fact[k] = fact[k + 1] * (k + 1)
    N = [top[D] * fact[D]]
    for k in range(D - 1, -1, -1):
        s = x0 + k
        new = [0] * (len(N) + 1)
        for i, a in enumerate(N):
            new[i + 1] += a
            new[i] -= a * s
        new[0] += top[k] * fact[k]
        N = new
    den = fact[0]
    g = gcd(den, *N)
    if g > 1:
        N = [c // g for c in N]
        den //= g
    while N and N[-1] == 0:
        N.pop()
    return N, den


def _ipoly_mul(a, b):
    """Произведение целых многочленов (списки коэффициентов по возрастанию степени)."""
    if not a or not b:
        return []
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def _ipoly_shift(p, s: int):
    """p(x + s) для целого многочлена — схема Горнера по множителю (x + s), над ℤ."""
    acc = []
    for c in reversed(p):
        new = [0] * (len(acc) + 1)
        for i, a in enumerate(acc):
            new[i + 1] += a
            new[i] += a * s
        new[0] += c
        acc = new
    while acc and acc[-1] == 0:
        acc.pop()
    return acc


def _dependency_interp(vectors, point: int = 100003):
    """Первая линейная зависимость векторов над ℚ(n) ИНТЕРПОЛЯЦИЕЙ (очередь п. 41,
    «интерполяция для больших порядков»): знаменатели каждого вектора снимаются —
    строки становятся целочисленными многочленами; порядок зависимости m и опорные
    компоненты находятся рангом при одной целой точке; коэффициенты
    λ_i = (−1)^i·det(B_{∖i}) — многочлены степени ≤ Σ deg строк — считаются целыми
    определителями (Барейсс) в D + 1 точках и восстанавливаются точной интерполяцией
    Ньютона; общий множитель снимается НОД; ответ СВЕРЯЕТСЯ символьно
    (Σ λ_j v_j = 0 над ℚ(n)); при любом сбое — исключение над ℚ(n) (_dependency).
    В исключении нет ни одного деления многочленов — замер 05.09 показал, что именно
    они (рост степеней в ℚ(n)) не давали 7 + 7 завершиться за 12 мин."""
    zero = ()
    rows, scales = [], []
    for v in vectors:
        L = _poly_lcm_list([f.den for f in v])
        row = [_exact_div(_mul(f.num, L), f.den) if not f.is_zero else zero for f in v]
        den = 1
        for q in row:
            for c in q:
                den = lcm(den, c.denominator)
        row = [tuple(c * den for c in q) for q in row]
        rows.append(row)
        scales.append(_mul(L, (Q(den),)))
    ncomp = len(rows[0])
    degs = [max((len(q) - 1 for q in row if q), default=0) for row in rows]
    # порядок зависимости и опорные компоненты — при одной точке
    num = [[_poly_eval(q, Q(point)) if q else Q(0) for q in row] for row in rows]
    m = None
    for k in range(len(rows)):
        rk, _ = _rank_q(num[:k + 1])
        if rk <= k:
            m = k
            break
    if m is None or m == 0:
        return _dependency(vectors)
    # m компонент, на которых первые m строк независимы: опорные столбцы исключения
    _, chosen = _rank_q([[num[j][c] for c in range(ncomp)] for j in range(m)])
    if len(chosen) != m:
        return _dependency(vectors)
    D = sum(degs[:m + 1])
    xs = list(range(D + 1))
    lam = []
    for i in range(m + 1):
        cols = [j for j in range(m + 1) if j != i]
        ys = []
        for x in xs:
            B = [[int(_poly_eval(rows[j][c], Q(x))) if rows[j][c] else 0 for j in cols] for c in chosen]
            ys.append(_int_det(B) * (1 if i % 2 == 0 else -1))
        lam.append(_newton_interp(xs, ys))
    if all(not q for q in lam):
        return _dependency(vectors)
    g = _poly_gcd_all(lam)
    if g and len(g) > 1:
        lam = [_exact_div(q, g) if q else q for q in lam]
    # символьная сверка на всех компонентах
    for c in range(ncomp):
        acc = ()
        for j in range(m + 1):
            if lam[j] and rows[j][c]:
                acc = _add(acc, _mul(lam[j], rows[j][c]))
        if _strip(acc):
            return _dependency(vectors)
    full = [_mul(lam[j], scales[j]) if lam[j] else () for j in range(m + 1)]
    g = _poly_gcd_all(full)                       # общий множитель, внесённый снятием знаменателей
    if g and len(g) > 1:
        full = [_exact_div(q, g) if q else q for q in full]
    out = [_rf(q) if q else _rf(()) for q in full]
    out += [_rf(())] * (len(vectors) - m - 1)
    return out


def _shift_rows_at(seq: DFiniteSeq, K: int, n: int):
    """Векторы сдвигов a_{n+k}, k = 0…K, в базисе a_n…a_{n+r−1} ПРИ ЦЕЛОМ n, УЖЕ
    умноженные на масштабы S_k(n) = Π_{j=0}^{k−r} p_r(n+j) — целые (V_k = S_k·u_k,
    V_k = −Σ_i p_i(n+k−r)·(S_{k−1}/S_{k−r+i})·V_{k−r+i}, деление точное: масштабы
    вложены); ни одной Q (11.09: rows_at на K = 20 стоил 51 с в Q). Возвращает
    (целые векторы, целые масштабы)."""
    rec, r = seq.rec, seq.rec.order
    ps = [list(q.c) for q in rec.p]

    def ev(cs, x):
        acc = 0
        for c in reversed(cs):
            acc = acc * x + c
        return acc

    vecs = [[1 if i == k else 0 for i in range(r)] for k in range(r)]
    scales = [1] * r
    for k in range(r, K + 1):
        sh = k - r
        lead = ev(ps[r], n + sh)
        if lead == 0:
            raise ZeroDivisionError("старший коэффициент рекурренты обращается в нуль в точке")
        s_prev = scales[k - 1]
        v = [0] * r
        for i in range(r):
            c = ev(ps[i], n + sh)
            if c == 0:
                continue
            idx = k - r + i
            f = c * (s_prev // scales[idx])
            u = vecs[idx]
            v = [v[j] - f * u[j] for j in range(r)]
        vecs.append(v)
        scales.append(s_prev * lead)
    return vecs, scales


def _dependency_points_multi(seqs, K: int):
    """Зависимость сдвигов СУММЫ n рядов ТОЛЬКО ПО ТОЧКАМ (очередь п. 41, 10.09):
    ни одной операции над ℚ(n). Строки [a_{n+k} | b_{n+k} | …] в базисе
    (a_n…, b_n…, …) считаются числами при целых n и умножаются на известные масштабы
    S_k(n) = Π_i Π_j p_{r_i}^{(i)}(n+j) — целые многочлены степени ≤ Σ_i (k−r_i+1)·deg_i.
    Порядок зависимости m и опорные компоненты — рангом в одной точке; λ_i =
    (−1)^i·det(B_{∖i}) — целые определители в D + 1 точках (ядро m × (m+1): Барейсс
    над ℤ или по модулю простых с CRT — _nullvector, порог по замеру), интерполяция
    по ПОСЛЕДОВАТЕЛЬНЫМ целым узлам конечными разностями целиком над ℤ (11.09; было:
    Ньютон над Q — 900 с из 1435 на K = 20). Доказательство тождества Σ λ_j·row_j = 0 —
    счётом степеней: обе части — многочлены степени ≤ D_ver, равенство проверено в
    D_ver + 1 точках на всех компонентах; значения λ_j в D_ver − D добавочных узлах —
    продолжением нижней кромки таблицы разностей (точно, O(D) сложений на узел; было:
    Горнер над Q — 405 с на K = 20). Сумма трёх и более рядов замыкается сразу (без
    промежуточного ряда, чьи коэффициенты раздувают степени: Лежандр 7 + 7 → 272,
    замер 10.09). Возвращает целые многочлены λ_0…λ_m (tuple Q со знаменателем 1,
    общее целое содержание снято) или None."""
    seqs = list(seqs)
    orders = [sq.rec.order for sq in seqs]
    degs_i = [max(len(q.c) - 1 for q in sq.rec.p) for sq in seqs]
    roots = [x for sq in seqs for x in _int_roots(sq.rec.p[-1])]
    n_start = max([sq.rec.n0 for sq in seqs] + [x + 1 for x in roots])
    degs = [sum(max(k - r + 1, 0) * d for r, d in zip(orders, degs_i)) for k in range(K + 1)]
    ncomp = sum(orders)

    def rows_at(n: int):
        parts = [_shift_rows_at(sq, K, n) for sq in seqs]
        out = []
        for k in range(K + 1):
            row = []
            for i, (vecs, _) in enumerate(parts):
                f = 1
                for i2, (_, sc2) in enumerate(parts):
                    if i2 != i:
                        f *= sc2[k]
                row.extend(x * f for x in vecs[k])
            out.append(row)
        return out

    probe = rows_at(n_start + 100003)
    m = None
    for k in range(K + 1):
        rk, _ = _rank_q([[Q(x) for x in probe[j]] for j in range(k + 1)])
        if rk <= k:
            m = k
            break
    if m is None or m == 0:
        return None
    _, chosen = _rank_q([[Q(x) for x in probe[j]] for j in range(m)])
    if len(chosen) != m:
        return None
    D = sum(degs[:m + 1])
    dmax = max(degs[:m + 1])
    # замер 10.09 (Лежандр 7 + 7): фактическая степень 272 при оценке 288 — удвоение D′ от 8
    # с доказательством на каждом шаге стоило 32 с против 15 с по оценке; берётся оценка
    D_ver = D + dmax
    xs = list(range(n_start, n_start + D_ver + 1))
    table = [rows_at(x)[:m + 1] for x in xs]
    # одним исключением на точку: вектор миноров (v_0…v_m) — ядро матрицы (компоненты × сдвиги)
    ys_all = []
    for t in range(D + 1):
        B = [[table[t][j][c] for j in range(m + 1)] for c in chosen]
        ys_all.append(_nullvector(B))
    if all(x == 0 for v in ys_all for x in v):
        return None
    # значения λ_j во всех D_ver + 1 узлах: в узлах интерполяции — сами миноры,
    # дальше — продолжение разностями; тождество проверяется над ℤ
    vals = []
    tops = []
    for i in range(m + 1):
        col = [ys_all[t][i] for t in range(D + 1)]
        top, last = _fwd_diffs(col)
        tops.append(top)
        vals.append(col + _extend_diffs(last, D_ver - D))
    for t in range(D_ver + 1):
        row_t = table[t]
        for c in range(ncomp):
            if sum(vals[j][t] * row_t[j][c] for j in range(m + 1)) != 0:
                return None
    # коэффициенты λ_j над ℤ с общим знаменателем, затем общее целое содержание
    lam = [_newton_interp_int(n_start, None, tops[i]) for i in range(m + 1)]
    L = 1
    for _, den in lam:
        L = lcm(L, den)
    ilam = [[c * (L // den) for c in N] for N, den in lam]
    # Σ μ_k·(S_k·row_k) = 0  ⇒  λ_k = μ_k·S_k, S_k = Π_i Π_j p_{r_i}(n+j) — многочлены сдвигов
    leads = [list(sq.rec.p[-1].c) for sq in seqs]
    full_i = []
    for k in range(m + 1):
        S = [1]
        for r, pl in zip(orders, leads):
            for j in range(k - r + 1):
                S = _ipoly_mul(S, _ipoly_shift(pl, j))
        full_i.append(_ipoly_mul(ilam[k], S) if ilam[k] else [])
    g = gcd(*[c for q in full_i for c in q])
    if g > 1:
        full_i = [[c // g for c in q] for q in full_i]
    full = [_strip(Q(c) for c in q) for q in full_i]
    # общий многочленный множитель — при умеренных степенях (Евклид над ℚ)
    g = _poly_gcd_all(full) if max(len(q) for q in full) <= 60 else ()
    if g and len(g) > 1:
        full = [_exact_div(q, g) if q else q for q in full]
    return full


def _dependency_points(seq_a: DFiniteSeq, seq_b: DFiniteSeq, K: int):
    """Двухоперандный случай _dependency_points_multi."""
    return _dependency_points_multi((seq_a, seq_b), K)


def _close_shift(seq_a: DFiniteSeq, seq_b: DFiniteSeq, direct) -> DFiniteSeq:
    """Замыкание СУММЫ на уровне рекуррент (05.09): c_{n+k} = a_{n+k} + b_{n+k}
    в базисе (a_n…a_{n+r−1}, b_n…b_{n+s−1}) над ℚ(n); первая линейная
    зависимость сдвигов даёт рекуррентту порядка ≤ r + s — без перехода к ОДУ,
    который раздувал порядок (замер 05.09: 7 + 7 → 22 через ОДУ). Начальные
    члены прямым счётом, сверка прогоном; n0 — за нулями старших коэффициентов."""
    r, s = seq_a.rec.order, seq_b.rec.order
    K = r + s
    polys = _dependency_points(seq_a, seq_b, K)      # только по точкам (10.09) — без ℚ(n)
    if polys is None:
        va, vb = _shift_vectors(seq_a, K), _shift_vectors(seq_b, K)
        lam = _dependency_interp([va[k] + vb[k] for k in range(K + 1)])   # интерполяция над строками ℚ(n)
        L = (Q(1),)
        for f in lam:
            L = _poly_lcm(L, f.den)
        polys = [_divmod(_mul(f.num, L), f.den)[0] for f in lam]
    while len(polys) > 1 and not polys[-1]:
        polys.pop()
    # нулевые младшие λ: зависимость начинается со сдвига t — порядок ниже, n0 выше на t
    t_shift = 0
    while len(polys) > 1 and not polys[0]:
        polys.pop(0)
        t_shift += 1
    if t_shift:
        polys = [_poly_shift(q, -t_shift) if q else () for q in polys]
    den = 1
    for q in polys:
        for c in q:
            den = lcm(den, c.denominator)
    ipolys = [IntPoly([int(c * den) for c in q] or [0]) for q in polys]
    n0 = max(seq_a.rec.n0, seq_b.rec.n0)
    for sq in (seq_a, seq_b):
        roots = _int_roots(sq.rec.p[-1])
        if roots:
            n0 = max(n0, roots[-1] + 1)          # рекуррента слагаемого применялась при номерах ≥ n
    rec = Recurrence(ipolys, n0 + t_shift)
    return _seq_from_rec(rec, direct, seq_a, seq_b)


def _close_sum_multi(seqs, direct):
    """Сумма n D-конечных последовательностей одним замыканием по точкам (10.09):
    порядок ≤ Σ r_i; при отказе точек — попарно (прежний путь)."""
    seqs = list(seqs)
    K = sum(sq.rec.order for sq in seqs)
    polys = _dependency_points_multi(seqs, K)
    if polys is None:
        return None
    while len(polys) > 1 and not polys[-1]:
        polys.pop()
    t_shift = 0
    while len(polys) > 1 and not polys[0]:
        polys.pop(0)
        t_shift += 1
    if t_shift:
        polys = [_poly_shift(q, -t_shift) if q else () for q in polys]
    den = 1
    for q in polys:
        for c in q:
            den = lcm(den, c.denominator)
    ipolys = [IntPoly([int(c * den) for c in q] or [0]) for q in polys]
    n0 = max(sq.rec.n0 for sq in seqs)
    for sq in seqs:
        roots = _int_roots(sq.rec.p[-1])
        if roots:
            n0 = max(n0, roots[-1] + 1)
    rec = Recurrence(ipolys, n0 + t_shift)
    return _seq_from_rec(rec, direct, seqs[0], seqs[-1])


def _seq_from_rec(rec: Recurrence, direct, seq_a, seq_b):
    need = rec.n0 + rec.order
    initial = [direct(n) for n in range(need)]
    if all(v == 0 for v in initial) and all(direct(n) == 0 for n in range(need, need + _VERIFY)):
        return None
    seq = DFiniteSeq(rec, initial, direct)       # дыры — прямым счётом при встрече
    for n in range(need, need + _VERIFY):
        if seq.term(n) != direct(n):
            raise Unresolved("замыкание", seq_a.record(), seq_b.record(),
                             hint=f"рекуррента результата разошлась с прямым счётом при n={n}")
    return seq


def _close(seq_a: DFiniteSeq, seq_b: DFiniteSeq, kind: str, direct) -> DFiniteSeq:
    """Замыкание: сумма — сдвигами либо через ОДУ, свёртка Коши — через ОДУ;
    начальные члены прямым счётом; сверка прогоном.

    ГРАНИЦА ПОРЯДКА У КАЖДОЙ ДОРОГИ СВОЯ, и это не мелочь: прежде здесь стояло
    «порядок ≤ r + s» без условия, а дорога ОДУ его не держит (глава 11 §20
    п. 6, замер 12.09 — 17 против 9 на одной и той же паре). По замеру:

    * сдвиги (`_close_shift`, выбираются при r + s ≤ _SHIFT_LIMIT) — порядок
      ≤ r + s, минимальный в базисе сдвигов;
    * ОДУ (`_combine_odes` + `_ode_to_rec`, за порогом и всегда для свёртки) —
      порядок ≤ rs + r + s: производные произведения по Лейбницу лежат в
      пространстве размерности rs + r + s + 1.

    Замер 12.09 [вычислено здесь]: (порядок 8) + (порядок 1) даёт **9**
    сдвигами и **17** через ОДУ; члены обоих путей совпали на 0…20."""
    if kind == "+" and seq_a.rec.order + seq_b.rec.order <= _SHIFT_LIMIT:
        return _close_shift(seq_a, seq_b, direct)   # эта дорога и только она держит порядок ≤ r + s
    ode = _combine_odes(_rec_to_ode(seq_a), _rec_to_ode(seq_b), kind)
    rec = _ode_to_rec(*ode)
    need = rec.n0 + rec.order
    initial = [direct(n) for n in range(need)]
    if all(v == 0 for v in initial) and all(direct(n) == 0 for n in range(need, need + _VERIFY)):
        return None          # определяющие и проверочные члены — нули: последовательность тождественно нулевая
    seq = DFiniteSeq(rec, initial, direct)       # дыры — прямым счётом при встрече
    for n in range(need, need + _VERIFY):
        if seq.term(n) != direct(n):
            raise Unresolved("замыкание", seq_a.record(), seq_b.record(),
                             hint=f"рекуррента результата разошлась с прямым счётом при n={n}")
    return seq


# ---------------------------------------------------------------------------
# ряд-величина: Σ tₙ (точка r вложена в члены)
# ---------------------------------------------------------------------------

def _nonneg_from(poly, n_star: int) -> bool:
    """Доказательство `p(n) ≥ 0` при всех `n ≥ n_star` сдвигом Тейлора.

    `[вычислено здесь]` Условие достаточное, не необходимое: если у `p(n_star + m)`
    все коэффициенты неотрицательны, то при `m ≥ 0` значение неотрицательно —
    сумма неотрицательных. Обратное неверно, потому «не доказано» здесь не
    значит «ложно»: ноль инструмента не есть ноль предмета.
    """
    return all(c >= 0 for c in _poly_shift(poly, n_star))


def _sign_from(poly, n_star: int) -> int:
    """Знак `p` на `[n_star, ∞)`: +1, −1 либо 0 (не доказан ни один)."""
    if _nonneg_from(poly, n_star):
        return 1
    if _nonneg_from(tuple(-c for c in poly), n_star):
        return -1
    return 0


def geometric_tail(obj, n_star: int, rho) -> dict:
    """Доказанная геометрическая мажоранта хвоста для рекуррент ЛЮБОГО порядка r ≥ 1.

    Слово корпуса: у голого D-конечного ряда мажоранты хвоста не было — стоял
    честный отказ (границы Меццаробба 2019 не строились). 11.09 закрыт порядок 1
    (`p₁(n)·a_{n+1} + p₀(n)·a_n = 0`, откуда `a_{n+1}/a_n = −p₀(n)/p₁(n)`);
    01.10 (вечер) — любой порядок r ≥ 1 ИНДУКЦИЕЙ ПО ОКНУ, без отношения соседей:

    `[вычислено здесь]` Пусть знаки s_i многочленов p_i на `[n_star, ∞)`
    доказаны сдвигом Тейлора (нулевой многочлен — знак 0, слагаемого нет) и
    при всех n ≥ n*
        ρ^r·|p_r(n)| − Σ_{i<r} ρ^i·|p_i(n)| ≥ 0                          (*)
    — многочлен Σ ±ρ^i·s_i·p_i(n), неотрицательность которого на `[n*, ∞)`
    проверяется тем же сдвигом. Тогда из |a_{n+i}| ≤ C·ρ^{n+i} (i < r) следует
    |a_{n+r}| = |Σ_{i<r} p_i(n)·a_{n+i}| / |p_r(n)| ≤ C·ρ^n·Σ|p_i(n)|ρ^i / |p_r(n)|
    ≤ C·ρ^{n+r} по (*), то есть |a_n| ≤ C·ρ^n держится для всех n от любого
    окна из r подряд идущих членов, и хвост мажорируется геометрической
    прогрессией. При r = 1 условие (*) есть прежнее ρ·|p₁(n)| − |p₀(n)| ≥ 0,
    и множитель хвоста прежний: ρ/(1−ρ) к |a_N|; при r ≥ 2 база — окно
    a_{N+1}, …, a_{N+r}: хвост Σ_{n>N}|a_n| ≤ max_{j<r}(|a_{N+1+j}|·ρ^{−j}) / (1−ρ).
    Ни одного float; всё — в ℚ.

    Возвращает словарь с доказательством и множителем хвоста; при
    недоказуемости — «доказано: False» либо отказ поимённо (знак не доказан,
    старший коэффициент — ноль в n*, порядок 0); ложного класса не возникает.
    """
    seq = obj.seq if isinstance(obj, DFiniteSeries) else obj
    if not isinstance(seq, DFiniteSeq):
        raise TypeError("geometric_tail — от DFiniteSeq или DFiniteSeries")
    rho = Q(rho)
    if not (0 <= rho < 1):
        raise ValueError("ρ — рациональное из [0, 1)")
    if not isinstance(n_star, int) or isinstance(n_star, bool) or n_star < seq.rec.n0:
        raise ValueError(f"n_star — целое ≥ n0 = {seq.rec.n0}")
    rec = seq.rec
    r = rec.order
    if r == 0:
        # Порядок 0 построением возможен (`Recurrence([p0])`): a_n = 0 при
        # n ≥ n0 — хвост точно нулевой, мажоранта ему не нужна. Текст обязан
        # быть верен о порядке 0 (глава 11 §20 п. 2, замер 12.09).
        raise Unresolved("мажоранта хвоста", _rec_of(obj),
                         hint="порядок 0: рекуррента порядка 0 задаёт a_n = 0 при "
                              "n ≥ n0 — хвост точно нулевой, и мажоранта ему не нужна")
    polys = [_qpoly(q) for q in rec.p]
    signs = []
    for i, q in enumerate(polys):
        if all(c == 0 for c in q):
            signs.append(0)
            continue
        si = _sign_from(q, n_star)
        if si == 0:
            raise Unresolved("мажоранта хвоста", _rec_of(obj),
                             hint=f"знак коэффициента p_{i} рекурренты на [n_star, ∞) не "
                                  "доказан сдвигом Тейлора — граница не строится")
        signs.append(si)
    if signs[r] == 0 or _poly_eval(polys[r], Q(n_star)) == 0:
        raise Unresolved("мажоранта хвоста", _rec_of(obj),
                         hint="старший коэффициент обращается в ноль в n_star — деления нет")
    cond = _mul((rho ** r,), tuple(signs[r] * c for c in polys[r]))
    for i in range(r):
        if signs[i]:
            cond = _add(cond, tuple(-signs[i] * rho ** i * c for c in polys[i]))
    proved = _nonneg_from(_strip(cond), n_star)
    out = {"порядок": r, "n*": n_star, "ρ": rho, "знаки": tuple(signs),
           "доказано": proved,
           "условие": ("ρ·|p₁(n)| − |p₀(n)| ≥ 0 при n ≥ n*" if r == 1 else
                       f"ρ^{r}·|p_{r}(n)| − Σ_{{i<{r}}} ρ^i·|p_i(n)| ≥ 0 при n ≥ n*"),
           "множитель хвоста": ((rho / (1 - rho)) if r == 1 else 1 / (1 - rho)) if proved else None,
           "база": "|a_N|" if r == 1 else f"max_{{j<{r}}} |a_{{N+1+j}}|·ρ^(−j)"}
    if r == 1:
        out["знак p1"], out["знак p0"] = signs[1], signs[0]
    return out


def tail_enclosure(series: "DFiniteSeries", N: int, n_star: int, rho) -> Interval:
    """Сертификат суммы через доказанную геометрическую мажоранту хвоста.

    `[вычислено здесь]` Порядок 1: при `|a_{n+1}/a_n| ≤ ρ` для `n ≥ n*` и
    `N ≥ n*` хвост не больше `|a_N|·ρ/(1−ρ)` — сумма геометрической прогрессии,
    мажорирующей члены поимённо. Порядок r ≥ 2 (01.10): по условию (*) из
    `geometric_tail` оценка |a_n| ≤ C·ρ^n распространяется с окна
    a_{N+1}, …, a_{N+r} на все n > N, и хвост не больше
    max_{j<r}(|a_{N+1+j}|·ρ^{−j}) / (1−ρ). Обрамление возвращается отрезком;
    ширина названа числом.
    """
    if not isinstance(series, DFiniteSeries):
        raise TypeError("tail_enclosure — от DFiniteSeries")
    if N < n_star:
        raise ValueError("N ≥ n_star: мажоранта доказана только от n_star")
    proof = geometric_tail(series, n_star, rho)
    if not proof["доказано"]:
        raise Unresolved("сертификат суммы", _rec_of(series),
                         hint=f"ρ = {proof['ρ']} не доказана при n ≥ {n_star}: возьмите большее ρ "
                              "или больший n*")
    r = proof["порядок"]
    if r == 1:
        bound = abs(series.term(N)) * proof["множитель хвоста"]
    else:
        rho = Q(rho)
        bound = max(abs(series.term(N + 1 + j)) / rho ** j for j in range(r)) * proof["множитель хвоста"]
    s = series.partial_sum(N)
    return Interval(s - bound, s + bound)


class DFiniteSeries:
    """Класс 7: Σ tₙ, (tₙ) D-конечна; строение помнит, из чего собрано."""

    __slots__ = ("seq", "parts", "_sums")

    def __init__(self, seq: DFiniteSeq, parts=("raw",)):
        if not isinstance(seq, DFiniteSeq):
            raise TypeError("seq — только DFiniteSeq")
        self.seq, self.parts = seq, parts
        self._sums = [seq.term(0)]

    # --- построение --------------------------------------------------------

    @classmethod
    def from_power(cls, ps) -> "DFiniteSeries":
        """Ряд линейки Σ aₙrⁿ с рекуррентой отношения → D-конечная
        последовательность членов tₙ = aₙrⁿ: den(n)·r_d·t_{n+1} − num(n)·r_n·tₙ = 0."""
        from .series import PowerSeriesValue
        if not isinstance(ps, PowerSeriesValue):
            raise TypeError("from_power — только PowerSeriesValue")
        rn, rd = ps.r.numerator, ps.r.denominator
        p0 = (-ps.ratio_num).scale(rn)
        p1 = ps.ratio_den.scale(rd)
        holes = _int_roots(p1)
        if holes:
            # тот же текст, что у линейки: один факт — один диагноз
            from .series import hole_hint
            raise Unresolved("from_power", ps.record(), hint=hole_hint(holes))
        return cls(DFiniteSeq(Recurrence([p0, p1], 0), [ps.a0]), ("leaf", ps))

    @classmethod
    def from_integral(cls, iv) -> "DFiniteSeries":
        """Класс 9 (∫_a^b num/den, f рациональна) как D-конечный ряд: отрезок
        делится на куски со сходящимся разложением первообразной, куски
        складываются замыканием; сертификат — от самого интеграла (обрамление
        Римана Interval-арифметикой), члены ряда — точные рациональные."""
        from .integral import IntegralValue
        if not isinstance(iv, IntegralValue):
            raise TypeError("from_integral — только IntegralValue")
        total = None
        for l, r in _integral_pieces(iv.num, iv.den, iv.a, iv.b):
            piece = _piece_series(iv.num, iv.den, l, r)
            total = piece if total is None else total + piece
        if isinstance(total, Q):
            return cls.const(total)
        return cls(total.seq, ("leaf", iv))

    @classmethod
    def from_product(cls, iv) -> "DFiniteSeries":
        """Класс 8 (Π fₙ, fₙ = P(n)/Q(n), n ≥ n0) как D-конечный ряд ТЕЛЕСКОПОМ —
        двойственность 7↔8 без логарифма: Π = p_{n0} + Σ_{N>n0} (p_N − p_{N−1}),
        d_N = p_{N−1}·(f_N − 1) — член с рациональным отношением
        d_{N+1}/d_N = P(N)·(P(N+1) − Q(N+1)) / (Q(N+1)·(P(N) − Q(N))),
        то есть рекуррента порядка 1 с многочленными коэффициентами (класс 7).
        Дыры (f_N = 1: старший коэффициент обращается в 0) закрываются прямым
        счётом начальных членов; рекуррента сверяется прогоном; сертификат —
        от самого произведения (обрамление хвоста Вейерштрасса)."""
        from .product import InfiniteProductValue
        if not isinstance(iv, InfiniteProductValue):
            raise TypeError("from_product — только InfiniteProductValue")
        P, Qd, n0 = iv.factor_num, iv.factor_den, iv.n0
        Dm = P - Qd                                          # числитель f − 1
        head = cls.const(iv.partial_product(n0))
        if Dm.is_zero():
            return cls(head.seq, ("leaf", iv))               # все множители 1 — ряд из одной головы

        def sh(poly: IntPoly, s: int) -> IntPoly:
            cs = _poly_shift(_qpoly(poly), s)
            return IntPoly([int(c.numerator) for c in cs] or [0])

        # t_k = d_{n0+1+k}, k ≥ 0: p1(k)·t_{k+1} + p0(k)·t_k = 0
        p1 = sh(Qd, n0 + 2) * sh(Dm, n0 + 1)
        p0 = -(sh(P, n0 + 1) * sh(Dm, n0 + 2))
        rec = Recurrence([p0, p1], 0)
        holes = [k for k in _int_roots(rec.p[-1]) if k >= 0]
        need = 1 + (holes[-1] + 1 if holes else 0)

        def d(k: int) -> Q:
            N = n0 + 1 + k
            return iv.partial_product(N - 1) * (iv.factor(N) - 1)

        initial = [d(k) for k in range(need)]
        if all(v == 0 for v in initial):
            return cls(head.seq, ("leaf", iv))
        tail_seq = DFiniteSeq(rec, initial)
        for k in range(need, need + _VERIFY):
            if tail_seq.term(k) != d(k):
                raise Unresolved("from_product", iv.record(),
                                 hint=f"телескоп: рекуррента разностей частных произведений разошлась с прямым счётом при k={k}")
        total = head + cls(tail_seq, ("raw",))
        if isinstance(total, Q):
            return cls.const(total)
        # самоприменение: частные суммы ряда суть частные произведения — S(K) = p_{n0+K+1}
        for K in (0, 1, 2, 5, 8):
            if total.partial_sum(K) != iv.partial_product(n0 + K + 1):
                raise Unresolved("from_product", iv.record(),
                                 hint=f"телескоп: S({K}) ≠ p_{n0 + K + 1} — двойственность не сошлась")
        return cls(total.seq, ("leaf", iv))

    @classmethod
    def const(cls, c) -> "DFiniteSeries":
        """Число как ряд c + 0·x + 0·x² + …: a_{n+1} = 0."""
        if isinstance(c, (float, bool)):
            raise TypeError("float запрещён построением")
        return cls(DFiniteSeq(Recurrence([IntPoly([0]), IntPoly([1])], 0), [Q(c)]), ("const", Q(c)))

    @staticmethod
    def _lift(x):
        from .series import PowerSeriesValue
        if isinstance(x, DFiniteSeries):
            return x
        if isinstance(x, PowerSeriesValue):
            return DFiniteSeries.from_power(x)
        from .integral import IntegralValue
        if isinstance(x, IntegralValue):
            return DFiniteSeries.from_integral(x)
        from .product import InfiniteProductValue
        if isinstance(x, InfiniteProductValue):
            return DFiniteSeries.from_product(x)
        if isinstance(x, bool):
            return None
        if isinstance(x, float):
            raise TypeError("float запрещён построением")
        if isinstance(x, _NUM):
            return DFiniteSeries.const(x)
        return None

    # --- члены и суммы -----------------------------------------------------

    def term(self, n: int) -> Q:
        return self.seq.term(n)

    def partial_sum(self, N: int) -> Q:
        if isinstance(N, bool) or not isinstance(N, int):
            raise TypeError("N — целое")
        if N < 0:
            raise ValueError("N ≥ 0")
        while len(self._sums) <= N:
            self._sums.append(self._sums[-1] + self.seq.term(len(self._sums)))
        return self._sums[N]

    @property
    def order(self) -> int:
        return self.seq.rec.order

    # --- сертификат сетки из строения ----------------------------------

    def enclosure(self, N: int = 40) -> Interval:
        kind = self.parts[0]
        if kind == "leaf":
            return self.parts[1].enclosure(N)
        if kind == "const":
            return Interval.point(self.parts[1])
        if kind == "scale":
            return self.parts[2].enclosure(N) * Interval.point(self.parts[1])
        if kind in ("+", "-", "×"):
            a, b = self.parts[1].enclosure(N), self.parts[2].enclosure(N)
            return a + b if kind == "+" else (a - b if kind == "-" else a * b)
        if kind == "sum":
            # Строение `sum_of` — то же сложение, только сразу многих слагаемых
            # (одно замыкание вместо цепочки попарных). Сертификат у него был
            # потерян: разбор знал "+" и не знал "sum", и ОДНА И ТА ЖЕ величина,
            # собранная двумя дорогами, в одном случае имела сертификат, в
            # другом нет (глава 11 §20 п. 3, замер 12.09). Сумма отрезков
            # конечного числа слагаемых — тот же приём, что у "+".
            acc = self.parts[1].enclosure(N)
            for part in self.parts[2:]:
                acc = acc + part.enclosure(N)
            return acc
        return self._auto_tail(N)

    _RHO_GRID = (Q(1, 2), Q(1, 4), Q(3, 4), Q(1, 8), Q(7, 8), Q(1, 16), Q(15, 16))
    _NSTAR_STEPS = (0, 1, 2, 4, 8, 16, 32)

    def _auto_tail(self, N: int) -> Interval:
        """Сертификат ГОЛОГО ряда (без строения из линеек) автоматическим
        подбором пары (n*, ρ) по сетке: n* ∈ n0 + {0, 1, 2, 4, 8, 16, 32},
        ρ ∈ {1/2, 1/4, 3/4, 1/8, 7/8, 1/16, 15/16} (01.10, вечер; до того —
        отказ «границы Меццаробба не строились»). Первая доказанная пара даёт
        `tail_enclosure(ряд, max(N, n*), n*, ρ)`; порядок 0 — точка частичной
        суммы (хвост нулевой). Ни одна пара не доказана — честный отказ по
        имени: сходимость не геометрическая (Σ1/n²) либо знак коэффициента не
        доказан; явный запрос `tail_enclosure` с своими n*, ρ остаётся."""
        rec = self.seq.rec
        if rec.order == 0:
            return Interval.point(self.partial_sum(max(N, rec.n0)))
        for step in self._NSTAR_STEPS:
            n_star = rec.n0 + step
            for rho in self._RHO_GRID:
                try:
                    proof = geometric_tail(self, n_star, rho)
                except Unresolved:
                    continue
                if proof["доказано"]:
                    return tail_enclosure(self, max(N, n_star), n_star, rho)
        raise Unresolved("enclosure", self.record(),
                         hint="геометрическая мажоранта голого ряда не доказана ни при одной "
                              "паре (n*, ρ) из сетки n* ∈ n0 + {0, 1, 2, 4, 8, 16, 32}, "
                              "ρ ∈ {1/2, 1/4, 3/4, 1/8, 7/8, 1/16, 15/16}: сходимость не "
                              "геометрическая либо знак коэффициента не доказан — сертификат "
                              "у строения из линеек, либо по явному запросу "
                              "tail_enclosure(ряд, N, n*, ρ)")

    # --- арифметика: замыкание ------------------------------------------

    def _binary(self, other, kind: str, swap: bool = False):
        o = self._lift(other)
        if o is None:
            return NotImplemented
        a, b = (o, self) if swap else (self, o)
        if kind == "+":
            seq = _close(a.seq, b.seq, "+", lambda n: a.term(n) + b.term(n))
        elif kind == "-":
            nb = -b
            seq = _close(a.seq, nb.seq, "+", lambda n: a.term(n) - b.term(n))
        else:
            seq = _close(a.seq, b.seq, "×",
                         lambda n: sum((a.term(k) * b.term(n - k) for k in range(n + 1)), Q(0)))
        if seq is None:
            return Q(0)      # взаимное уничтожение — точный ноль (как у линейки по a0)
        return DFiniteSeries(seq, (kind, a, b))

    @classmethod
    def sum_of(cls, series) -> "DFiniteSeries":
        """Сумма нескольких рядов ОДНИМ замыканием по точкам (порядок ≤ Σ rᵢ): без
        промежуточных рядов, чьи коэффициенты раздувают степени (очередь п. 41, 10.09);
        при отказе точек — попарное сложение."""
        items = [x if isinstance(x, cls) else cls.const(x) for x in series]
        if not items:
            raise ValueError("сумма пустого набора рядов")
        if len(items) == 1:
            return items[0]
        seq = _close_sum_multi([x.seq for x in items],
                               lambda n: sum((x.term(n) for x in items), Q(0)))
        if seq is None:
            acc = items[0]
            for x in items[1:]:
                acc = acc + x
            return acc
        return cls(seq, ("sum",) + tuple(items))

    def __add__(self, other):
        if isinstance(other, Zero):
            return NotImplemented
        return self._binary(other, "+")

    def __radd__(self, other):
        if isinstance(other, Zero):
            return NotImplemented
        return self._binary(other, "+", swap=True)

    def __sub__(self, other):
        if isinstance(other, Zero):
            return NotImplemented
        return self._binary(other, "-")

    def __rsub__(self, other):
        if isinstance(other, Zero):
            return NotImplemented
        return self._binary(other, "-", swap=True)

    def __neg__(self) -> "DFiniteSeries":
        return self.scale(-1)

    def scale(self, c) -> "DFiniteSeries":
        c = Q(c)
        if c == 0:
            return Q(0)
        seq = DFiniteSeq(self.seq.rec, [c * v for v in self.seq.initial])
        return DFiniteSeries(seq, ("scale", c, self))

    def __mul__(self, other):
        if isinstance(other, Zero):
            return NotImplemented
        if isinstance(other, bool):
            return NotImplemented
        if isinstance(other, float):
            raise TypeError("float запрещён построением")
        if isinstance(other, _NUM):
            return self.scale(other)
        return self._binary(other, "×")

    __rmul__ = __mul__

    def __truediv__(self, other):
        if isinstance(other, _NUM) and not isinstance(other, bool):
            q = Q(other)
            if q == 0:
                raise Unresolved("/", self.record(), _rec_of(other), hint="x/0 здесь не сведено")
            return self.scale(1 / q)
        if self._lift(other) is not None:
            raise Unresolved("/", self.record(), _rec_of(other),
                             hint="1/f D-конечна лишь в частных случаях — общего замыкания "
                                  "деления рядов нет (теорема, не незнание)")
        return NotImplemented

    def __rtruediv__(self, other):
        if isinstance(other, _NUM) or self._lift(other) is not None:
            raise Unresolved("/", _rec_of(other), self.record(),
                             hint="обращение ряда: общего замыкания нет")
        return NotImplemented

    def __pow__(self, n):
        if isinstance(n, bool) or not isinstance(n, int):
            return NotImplemented
        if n < 1:
            raise Unresolved("^", self.record(), _rec_of(n),
                             hint="степень ряда сведена для целых ≥ 1 (свёртками); "
                                  "0 не решается строением, отрицательная — обращение")
        out = self
        for _ in range(n - 1):
            out = out * self
        return out

    # --- запись и тождество -------------------------------------------

    # --- самопроверка носителя (ось самоприменения, 08.10) ----------------

    def self_check(self, deep: bool = False) -> list:
        """Самопроверка носителя встречными дорогами (ось «AIAΣA спрашивает
        AIAΣA», слово автора 04.09; очередь строка 18; 08.10). Список
        находок, пусто — здоров; исключений не бросает.

        `[вычислено здесь]` (1) слоты: seq — DFiniteSeq, рекуррента в каноне
        (старший коэффициент не нуль, содержание 1), строение из известных
        видов; (2) невязка рекурренты на 24 номерах от n0 вне дыр старшего
        коэффициента — нуль до последнего бита; (3) ВСТРЕЧНАЯ дорога у замыканий
        (+, −, ×, sum, scale): члены, выданные рекуррентой результата, против
        прямого счёта из членов операндов (`_direct`) на 24 номерах — две
        независимые дороги (исключение Гаусса против арифметики ℚ); у листа
        линейки — члены против самого листа (`term` носителя класса 7 линейки);
        (4) сертификат: enclosure(40) и enclosure(80) пересекаются и длинный не
        шире; отказ сертификата (мажоранта хвоста вне сетки доказательств) — объявленная граница,
        не находка. При deep — 48 номеров."""
        out = []
        seq = self.seq
        if not isinstance(seq, DFiniteSeq):
            return [f"слот seq чужого носителя: {type(seq).__name__}"]
        rec = seq.rec
        if not isinstance(rec, Recurrence) or rec.p[-1].is_zero():
            return ["рекуррента без старшего коэффициента"]
        g = 0
        for q in rec.p:
            g = gcd(g, q.content()) if not q.is_zero() else g
        if g > 1:
            out.append(f"рекуррента вне канона: общее содержание {g}")
        if not self.parts or self.parts[0] not in ("raw", "leaf", "const", "scale", "+", "-", "×", "sum"):
            out.append(f"строение неизвестного вида: {self.parts[:1]!r}")
        n_max = 48 if deep else 24
        try:
            r = rec.order
            for n in range(rec.n0, rec.n0 + n_max):
                if n in seq._holes or rec.p[-1].eval(n) == 0:
                    continue
                s = sum((rec.p[i].eval(n) * seq.term(n + i) for i in range(r + 1)), Q(0))
                if s != 0:
                    out.append(f"невязка рекурренты при n = {n}: {s}")
                    break
        except Unresolved as e:
            out.append(f"рекуррента не продлена: {e}")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"невязка не проверена: {type(e).__name__}: {e}")
        try:
            if seq._direct is not None:
                for n in range(n_max):
                    if seq.term(n) != Q(seq._direct(n)):
                        out.append(f"член замыкания t_{n} = {seq.term(n)} против прямого счёта {Q(seq._direct(n))}")
                        break
            elif self.parts[0] == "leaf" and hasattr(self.parts[1], "term"):
                leaf = self.parts[1]
                for n in range(n_max):
                    if seq.term(n) != leaf.term(n):
                        out.append(f"член t_{n} = {seq.term(n)} против листа линейки {leaf.term(n)}")
                        break
        except Unresolved as e:
            out.append(f"встречная дорога членов отказала: {e}")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"встречная дорога членов: {type(e).__name__}: {e}")
        try:
            e1, e2 = self.enclosure(40), self.enclosure(80)
            if e1.intersect(e2) is None:
                out.append(f"сертификаты на 40 и 80 членах не пересекаются: {e1.record()} против {e2.record()}")
            if e2.width() > e1.width():
                out.append("сертификат на 80 членах шире, чем на 40")
        except Unresolved:
            pass                                                # граница: мажоранта хвоста вне сетки доказательств — отказ законный
        except Exception as e:                                  # noqa: BLE001
            out.append(f"сертификат: {type(e).__name__}: {e}")
        return out

    def _shape(self) -> str:
        k = self.parts[0]
        if k == "leaf":
            return self.parts[1].record()
        if k == "const":
            return str(self.parts[1])
        if k == "scale":
            return f"{self.parts[1]}·({self.parts[2]._shape()})"
        if k in ("+", "-", "×"):
            return f"({self.parts[1]._shape()}) {k} ({self.parts[2]._shape()})"
        if k == "sum":
            # Вторая половина той же дыры, что и в `enclosure` (глава 11 §20
            # п. 3): разбор строения научили «sum», а ЗАПИСЬ — нет, и результат
            # `sum_of` печатал «строение=raw», то есть объявлял, что строения у
            # него нет, — при том, что сертификат берётся ИМЕННО из строения.
            # Замер 17.09 `[вычислено здесь]`: попарная сумма писала
            # «(…) + (…)», а sum_of тех же слагаемых — «raw».
            return " + ".join(f"({p._shape()})" for p in self.parts[1:])
        return "raw"

    def record(self) -> str:
        return (f"Φ(C=7 POWER_SERIES D-finite, r={self.order}, {self.seq.record()}, "
                f"строение={self._shape()})")

    __repr__ = record

    def __eq__(self, other):
        o = self._lift(other) if not isinstance(other, DFiniteSeries) else other
        if o is None:
            return NotImplemented
        if self.seq.rec == o.seq.rec and self.seq.initial == o.seq.initial:
            return True
        raise Unresolved("=", self.record(), o.record(),
                         hint="тождество сумм разных строений не сведено (канона записи нет)")

    def __hash__(self):
        return hash(("dfinite", self.seq.rec, self.seq.initial))


# ---------- самопроверка линейки (ось самоприменения, 08.10) ----------

def self_check(deep: bool = False) -> list:
    """Самопроверка линейки класса 7 (замыкание D-конечных рядов) встречными
    дорогами; список находок, пусто — здорова; исключений не бросает.

    `[вычислено здесь]` Листы: Σ 1/n! (экспонента) — сертификат обязан нести e
    дверью `exp_enclosure` (дорога ряда против двери), Σ (1/2)ⁿ — нести 2
    точно [внешний источник: геометрическая прогрессия]; замыкания: сумма
    экспоненты и прогрессии несёт e + 2, тождество cos²(1) + sin²(1) = 1 — сертификат
    суммы двух свёрток Коши несёт единицу (замыкание против внешнего тождества), произведение прогрессии
    на саму себя — члены (n + 1)/2ⁿ; `const` — частная сумма есть число;
    `Recurrence` сокращает общее содержание в канон. Носители окна проверяют
    себя (`DFiniteSeries.self_check`). Граница — окно."""
    out = []
    try:
        from .series import PowerSeriesValue
        from .enclosure import exp_enclosure
        exp_s = DFiniteSeries.from_power(PowerSeriesValue(1, IntPoly([1]), IntPoly([1, 1]), 1))
        geo = DFiniteSeries.from_power(PowerSeriesValue(1, IntPoly([1]), IntPoly([1]), Q(1, 2)))
        cos_s = DFiniteSeries.from_power(PowerSeriesValue(1, IntPoly([-1]), IntPoly([2, 6, 4]), 1))
        sin_s = DFiniteSeries.from_power(PowerSeriesValue(1, IntPoly([-1]), IntPoly([6, 10, 4]), 1))
        window = [("Σ 1/n!", exp_s), ("Σ (1/2)^n", geo), ("exp + geo", exp_s + geo), ("geo × geo", geo * geo)]
        if deep:
            window += [("cos² + sin²", cos_s * cos_s + sin_s * sin_s), ("exp − geo", exp_s - geo)]
        for name, s in window:
            out += [f"{name}: {f}" for f in s.self_check(deep)]
        e_iv = exp_enclosure(Q(1), 40)
        if exp_s.enclosure(40).intersect(e_iv) is None:
            out.append(f"Σ 1/n!: сертификат {exp_s.enclosure(40).record()} не несёт e")
        if not geo.enclosure(40).contains(Q(2)):
            out.append(f"Σ (1/2)^n: сертификат {geo.enclosure(40).record()} не несёт 2")
        if (exp_s + geo).enclosure(40).intersect(e_iv + Interval.point(Q(2))) is None:
            out.append("exp + geo: сертификат не несёт e + 2")
        one = cos_s * cos_s + sin_s * sin_s
        if not one.enclosure(20).contains(Q(1)):
            out.append(f"cos² + sin²: сертификат {one.enclosure(20).record()} не несёт 1")
        gg = geo * geo
        if any(gg.term(n) != Q(n + 1, 2 ** n) for n in range(12)):
            out.append("geo × geo: члены не (n + 1)/2^n")
        if DFiniteSeries.const(Q(5)).partial_sum(7) != Q(5):
            out.append("const(5): частная сумма ≠ 5")
        if Recurrence([IntPoly([2]), IntPoly([4])]) != Recurrence([IntPoly([1]), IntPoly([2])]):
            out.append("Recurrence: общее содержание не сокращено в канон")
    except Exception as e:                                      # noqa: BLE001
        out.append(f"окно: {type(e).__name__}: {e}")
    return out
