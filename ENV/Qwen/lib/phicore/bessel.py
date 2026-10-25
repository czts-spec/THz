"""Линейка класса 16 BESSEL_ROOT: k-й положительный корень J_ν.

Величина хранится структурно: порядок ν (целый либо рациональный) и
номер корня k ≥ 1.
Правил сведения корня с независимыми величинами нет — сумма с числом,
произведение, частное, степень вне {0, 1} получают отказ Unresolved
с явными Ф-записями обоих операндов, не ложный класс и не молчание.

ДВОЙСТВЕННОСТЬ с классом 10: полуцелые порядки приземляются в мономы
от π ещё фабрикой — J_{1/2}(x) = √(2/(π·x))·sin x, корни в k·π;
J_{−1/2}(x) = √(2/(π·x))·cos x, корни в (2k−1)·π/2 [внешний источник,
замкнутые формы полуцелых порядков]. Прочие нецелые ν направленного
правила не имеют — честный отказ.

Десятичная сторона — только сертификатом сетки (enclosure), по явному
запросу: восходящий ряд J_ν знакочередуется с убывающими членами после
старта мажоранты (частные суммы обрамляют — j_enclosure); счёт ряда —
целыми аккумуляторами на масштабе 2^−(g+запас) с направленными наружу
делениями, итог — на сетке 2^−g (приём этапа C, enclosure.py); корень
изолируется k-й сменой сертифицированного знака на сетке шага 1/8 и
уже́няется делением пополам до ширины 2^−bits; при нецелом ν < −1/2
(ν ∈ (−M−1; −M), M = ⌊−ν⌋ ≥ 0) скан начинается с узла под лейбницевой
нижней границей первого корня j_{ν,1} > 2√((M+1)(ν+M+1)) — первые M
членов ряда одного знака с лейбницевым хвостом, при M = 0 это 2√(ν+1)
[вычислено здесь]; соседние корни разнесены больше чем на π (ν² > 1/4,
сравнение Штурма в нормальной форме DLMF 10.13.1 [вычислено здесь]); у
самой границы −(M+1) шаг скана до первого корня — степень двойки под
границей (круг шестой 07.10), отказов по ν у линейки нет. Значение J_ν
при нецелом ν < −1 — спуском рекуррентой по порядку J_{μ−1} =
(2μ/x)·J_μ − J_{μ+1} от пары порядков в (−1; 1) [внешний источник:
Абрамовиц—Стиган 9.1.27] (круг пятый 07.10); (x/2)^ν при знаменателе ν
выше 64 — композицией exp(ν·ln(x/2)) (круг шестой), цена не зависит от
знаменателя. Рациональная точка не
есть корень J_ν целого порядка (положительные корни трансцендентны
[внешний источник, теорема Зигеля]) — сертификация знака завершима.
Ни одного float: float на входе — TypeError построением.
"""
from __future__ import annotations

import math

from .rational import Q

from .enclosure import (Interval, _cdiv, _grid_out, _GUARD, _scale_pair,
                        cos_enclosure, exp_enclosure, ln_enclosure, pi_enclosure,
                        root_enclosure, sin_enclosure, sqrt_enclosure)
from .errors import Unresolved
from .pi_value import PiMonomial
from .zero_state import Zero

#: Допустимые скаляры. float исключён построением.
_NUM = (int, Q)

#: `[вычислено здесь]` 10.10, очередь строка 105, шаг шестой: дверь оператора
#: носителя = дверь ФIOS, как у π, e^y/ln, arctrig, ζ и Γ. Для ЧУЖОГО класса с
#: записью операторы корня J_ν зовут `dispatch.combine`; класс 1 (int, Q) и
#: одноимённая линейка — прежний отказ носителя («направленного правила
#: арифметики нет»). Ввоз ленивый, внутри функций.
_IN_DOOR: set = set()
#: Прибитые ячейки — proto/test_phicore_bessel_dual.py:219–230 (проверка B:
#: восемь форм π ∘ J и J ∘ π, текст «не сведено» с записями C=16 и C=10);
#: со стороны монома те же ячейки прибиты как (16, все четыре).
_PINNED: frozenset = frozenset({(10, "+"), (10, "-"), (10, "*"), (10, "/")})


def _pinned(other, op: str) -> bool:
    from .Фios import classify, EXACT      # ленивый ввоз
    k = classify(other)
    return k is not None and k[0] == EXACT and (k[1], op) in _PINNED


def _door(a, b, op: str):
    from .dispatch import combine          # ленивый ввоз
    key = (id(a), id(b), op)
    if key in _IN_DOOR:
        return NotImplemented
    _IN_DOOR.add(key)
    try:
        return combine(a, b, op)
    finally:
        _IN_DOOR.discard(key)


def _foreign_record(x) -> bool:
    """Чужой класс с записью: не класс 1 и не одноимённая линейка."""
    return (not isinstance(x, (_NUM, BesselRootValue))) and callable(getattr(x, "record", None))
#: Одна запись отказа на bool во всех дверях модуля: bool числом в доме не
#: считается и называется bool'ом (канон `phicore/rational.py`; сосед по
#: главе отвечает тем же: `Notch.from_q(True, 8)` → «bool не число»).
#: Двух копий одного канона здесь нет намеренно: разошедшиеся копии — уже
#: пройденная поломка (`affine.py:75`).
_NOT_BOOL = "bool не число"

#: Потолок углубления сертификата знака: глубже — честный отказ.
_MAX_TERMS = 1 << 13
#: Знаменатель порядка, до которого (x/2)^ν берётся целым корнем на сетке;
#: выше — композицией exp(ν·ln(x/2)). Замер 07.10 (круг шестой): q-й корень
#: на сетке 200·q бит — 0,012 с при q = 64, 0,2 с при 256, 4 с при 1024,
#: 141 с при 4096; exp∘ln — < 1 мс при любом q [вычислено здесь].
_ROOT_DENOM_MAX = 64
#: Члены ряда ln с приведением к [1; 2): аргумент atanh ≤ 1/3, хвост ≤ 3^−2N.
_LN_TERMS = 64


def _rec(x: object) -> str:
    """Ф-запись операнда для сообщений отказа."""
    r = getattr(x, "record", None)
    if callable(r):
        return r()
    if isinstance(x, _NUM):
        return f"Φ(C=1 RATIONAL, value={Q(x)})"
    return f"Φ(?={x!r})"


def _fact(n: int) -> int:
    """n! точной рекуррентой k! = k·(k−1)!; таблиц значений нет."""
    out = 1
    for k in range(2, n + 1):
        out *= k
    return out


# ---------- сертификат ряда: вложение J_ν(x) ----------

def _half_order_enclosure(nu: Q, x: Q, terms: int) -> Interval:
    """J_{±1/2}(x) замкнутой формой: J_{1/2}(x) = √(2/(π·x))·sin x,
    J_{−1/2}(x) = √(2/(π·x))·cos x [внешний источник: Абрамовиц—Стиган
    10.1.11, 10.1.12]. Три сертификата дома перемножаются: π (Мэчин),
    √ (целый корень на сетке, концами отрезка 2/(π·x)) и sin/cos
    (знакочередующийся ряд при |x| ≤ 4, дальше — поворот периода с
    поправкой Липшица: `enclosure._rotate_period`, 01.10).
    x = 0: J_{1/2}(0) = 0 точно; J_{−1/2}(0) — полюс, отказ."""
    if x == 0:
        if nu > 0:
            return Interval.point(0)
        raise Unresolved("J_ν", _rec(nu), _rec(x),
                         hint="J_{−1/2}(0) — полюс √(2/(π·x)): значения нет")
    if x < 0:
        raise Unresolved("J_ν", _rec(nu), _rec(x),
                         hint="x < 0 у полуцелого порядка: √x вне ℝ — "
                              "пара (re, im) этой линейкой не несётся")
    bits = 8 * terms
    trig = sin_enclosure(x, terms) if nu > 0 else cos_enclosure(x, terms)
    inner = (pi_enclosure(max(30, bits)) * x).inverse() * Q(2)   # 2/(π·x) ⊂ ℚ, lo > 0
    root = Interval(sqrt_enclosure(inner.lo, bits).lo, sqrt_enclosure(inner.hi, bits).hi)
    return root * trig


_INV_GAMMA_MEMO: dict = {}


def _inv_gamma(nu: Q, m: int, terms: int) -> Interval:
    """1/Γ(m+ν+1) отрезком — КЭШ по (ν, m, terms) (02.10, круг третий).

    `[вычислено здесь]` Γ(m+ν+1) не зависит от x, а скан корня зовёт
    `j_enclosure` на сотнях узлов при одном ν: до кэша каждый узел заново
    гнал Стирлинг на всех m — ≈2,4 с на корень при bits = 32; с кэшем
    Γ считается один раз на (ν, m, terms). Ключ несёт terms: глубина
    Стирлинга — часть сертификата, и отрезок иной глубины — иной отрезок.
    Память: по одному отрезку на член ряда; `_MAX_TERMS` ограничивает рост."""
    key = (nu, m, terms)
    iv = _INV_GAMMA_MEMO.get(key)
    if iv is None:
        from .gamma_value import GammaValue
        g = GammaValue.of(m + nu + 1)
        g_iv = Interval.point(g) if isinstance(g, Q) else g.enclosure(terms)
        iv = g_iv.inverse()
        _INV_GAMMA_MEMO[key] = iv
    return iv


# ln с приведением к [1; 2) — одна копия в доме: `gamma_value._ln_reduced` (та же
# подпись; замер 07.10: собственная копия здесь разошлась бы с ней при первой же
# правке — пройденная поломка, `affine.py:75`).
from .gamma_value import _ln_reduced  # noqa: E402


def _power_enclosure(half: Q, nu: Q, terms: int, road: str = "") -> Interval:
    """(x/2)^ν отрезком, ν = p/q ∉ ℤ, x/2 > 0 — две дороги (07.10, круг шестой):
    при q ≤ `_ROOT_DENOM_MAX` — целый корень на сетке (`root_enclosure`,
    ширина 2^−190 при terms = 24); выше — композиция дверей дома
    exp(ν·ln(x/2)) (точное тождество), цена которой от q не зависит: ln с
    приведением к [1; 2) (`_ln_reduced`), exp с числом членов по |ν·ln(x/2)|
    (хвост |a|^N/N!, `exp_enclosure`), обе монотонны — концы переносятся
    концами. `road` ∈ {«root», «expln»} — принудительный выбор дороги для
    приёмки (две дороги обязаны пересекаться); пусто — по знаменателю.
    До 07.10 стоял только корень, и ν со знаменателем 256 стоил 0,3 с на
    точку, 4096 — 434 с; теперь — миллисекунды [вычислено здесь]."""
    pn, qn = nu.numerator, nu.denominator
    if road == "root" or (not road and qn <= _ROOT_DENOM_MAX):
        bits = 8 * terms + _GUARD
        if pn > 0:
            return root_enclosure(Interval.point(half ** pn), qn, bits)
        return root_enclosure(Interval.point(half ** (-pn)), qn, bits).inverse()
    a = _ln_reduced(half, max(_LN_TERMS, terms)) * nu
    mag = max(abs(a.lo), abs(a.hi))
    te = terms + 24 + 16 * (mag.numerator // mag.denominator)
    return Interval(exp_enclosure(a.lo, te).lo, exp_enclosure(a.hi, te).hi)


def _rational_order_enclosure(nu: Q, x: Q, terms: int) -> Interval:
    """J_ν(x) при рациональном ν ∉ ℤ, ν > −1, x > 0 (02.10): восходящий ряд
    с СЕРТИФИКАТАМИ Γ класса 14 вместо точных факториалов.

    J_ν(x) = (x/2)^ν · Σ_{m≥0} (−1)^m·u_m,  u_m = (x²/4)^m / (m!·Γ(m+ν+1))
    [внешний источник: Абрамовиц—Стиган 9.1.10]. Γ(m+ν+1) — дверь класса
    14 (`GammaValue.of`: рекуррента к базе ∈ (0, 1), Стирлинг с остатком),
    обратная — `Interval.inverse` (отрезок Γ при аргументе > 0 нуля не
    содержит); (x/2)^ν = ((x/2)^p)^{1/q} — целый корень на сетке
    (`root_enclosure`), при ν < 0 — обратный. Множители (x²/4)^m/m! точны
    в ℚ. Знакочередование и убывание |u_m|: u_{m+1}/u_m = (x²/4)/((m+1)(m+ν+1)),
    и при ν > −1 знаменатель растёт по m монотонно — как только
    (x²/4) < (m+1)(m+ν+1), ряд с этого места лейбницев, и предел лежит
    между S_m и S_{m+1}; их внешняя оболочка (отрезками) — сертификат.
    До 02.10 здесь стоял отказ «Γ(m+ν+1) при знаменателе ν вне {1, 2} дом
    точно не несёт» — точно не несёт и теперь, но сертификат есть.
    Границы по имени: x < 0 — (x/2)^ν вне ℝ, пара (re, im) не несётся;
    x = 0 — 0 при ν > 0, полюс при ν < 0. ν < −1 — Γ(m+ν+1) первых членов
    отрицательна, знакочередования нет; до 07.10 здесь стоял отказ «дорога
    рекуррентой по ν не построена» — с круга пятого она построена:
    `_below_minus_one_enclosure`, спуск от пары порядков в (−1; 1).
    """
    if x == 0:
        if nu > 0:
            return Interval.point(0)
        raise Unresolved("J_ν", _rec(nu), _rec(x),
                         hint="J_ν(0) при ν < 0 — полюс (x/2)^ν: значения нет")
    if x < 0:
        raise Unresolved("J_ν", _rec(nu), _rec(x),
                         hint="x < 0 у нецелого порядка: (x/2)^ν вне ℝ — "
                              "пара (re, im) этой линейкой не несётся")
    if nu < -1:
        return _below_minus_one_enclosure(nu, x, terms)
    r2 = x * x / 4
    coef = Q(1)                      # (x²/4)^m / m!
    s = Interval.point(0)
    m = 0
    while True:
        u = _inv_gamma(nu, m, terms) * coef
        s_next = s + u if m % 2 == 0 else s - u
        if m + 1 >= terms and r2 < (m + 1) * (m + nu + 1):
            hull = Interval(min(s.lo, s_next.lo), max(s.hi, s_next.hi))
            break
        s = s_next
        coef = coef * r2 / (m + 1)
        m += 1
    power = _power_enclosure(x / 2, nu, terms)
    out = power * hull
    return Interval.on_grid(out.lo, out.hi, 8 * terms)


def _below_minus_one_enclosure(nu: Q, x: Q, terms: int) -> Interval:
    """J_ν(x) при рациональном ν < −1, ν ∉ ℤ, x > 0 (07.10, круг пятый):
    спуск рекуррентой по порядку J_{μ−1}(x) = (2μ/x)·J_μ(x) − J_{μ+1}(x)
    [внешний источник: Абрамовиц—Стиган 9.1.27] от пары порядков
    ν₀ = ν + M ∈ (−1; 0) и ν₀ + 1 ∈ (0; 1), M = ⌊−ν⌋ ≥ 1 — обе несутся
    сертификатом ряда (`_rational_order_enclosure`; ±1/2 — замкнутыми
    формами). Рекуррента — точное тождество, потому отрезок верен при
    любой ширине; множитель 2μ/x точен в ℚ, спуск — M шагов направленной
    арифметикой отрезков. Ширина за шаг растёт не быстрее, чем в
    |2μ/x| + 1 раз — при x порядка единиц это считанные двоичные разряды
    на шаг; цена лечится числом членов terms, не новой дорогой. Это
    композиция дверей дома, не новая аналитика: знакочередования ряда при
    ν < −1 по-прежнему нет (Γ(m+ν+1) первых M членов отрицательна) — и оно
    не нужно. Граница по имени: x < 0 (та же, что у ряда).
    """
    big_m = (-nu.numerator) // nu.denominator      # ⌊−ν⌋ ≥ 1 точно в ℤ
    nu0 = nu + big_m                                # ∈ (−1; 0), нецелое
    a = j_enclosure(nu0, x, terms)                  # J_{ν₀}
    b = j_enclosure(nu0 + 1, x, terms)              # J_{ν₀+1}
    mu = nu0
    for _ in range(big_m):
        a, b = a * (2 * mu / x) - b, a              # J_{μ−1}, J_μ
        mu -= 1
    return Interval.on_grid(a.lo, a.hi, 8 * terms)


def j_enclosure(nu, x, terms: int = 24) -> Interval:
    """Сертификат сетки [lo, hi] ⊂ ℚ для J_ν(x): ν ∈ ℤ (любого знака),
    ν = ±1/2 (класс 10 — замкнутые формы, `_half_order_enclosure`) и, с
    02.10, любое рациональное ν > −1 (`_rational_order_enclosure`:
    восходящий ряд с сертификатами Γ класса 14), x ∈ ℚ.

    Отрицательные целые — зеркалом J_{−n}(x) = (−1)^n·J_n(x) [внешний
    источник: Абрамовиц—Стиган 9.1.5]; до 01.10 дверь принимала только
    ν ∈ ℤ, ν ≥ 0 (наблюдение сдачи 25.09, очередь строка 104); до 02.10
    нецелые ν вне ±1/2 отвечали отказом «Γ(m+ν+1) при знаменателе ν вне
    {1, 2} дом точно не несёт» — точно не несёт и теперь, сертификатом
    несёт; до 07.10 нецелые ν < −1 отвечали отказом «ряд не
    знакочередуется» — с круга пятого несутся спуском рекуррентой по
    порядку от пары порядков в (−1; 1) (`_below_minus_one_enclosure`).
    Остаток по имени: x < 0 при нецелом ν.

    Ряд при целом ν ≥ 0 — ниже, как и прежде:

    Восходящий ряд J_ν(x) = Σ (−1)^m (x/2)^{2m+ν}/(m!·(m+ν)!); члены
    считаются рекуррентой t_{m+1} = −t_m·(x/2)²/((m+1)(m+ν+1)),
    без таблиц. Ряд знакочередуется, и модули членов строго убывают,
    как только (x/2)² < (m+2)(m+ν+2) — знаменатель рекурренты растёт
    монотонно, потому убывание, раз начавшись, не прерывается.
    Суммирование ведётся до старта этой мажоранты И не короче terms
    членов; после этого две соседние частные суммы S_m и S_{m+1}
    обрамляют предел (теорема Лейбница) — они и есть сертификат.

    Счёт — целыми аккумуляторами на масштабе 2^−G, G = g + запас
    (приём этапа C, enclosure.py): член и сумма — пары целых границ,
    деления направлены к краям (низ — вниз, верх — вверх), потому
    пара [s_lo, s_hi] несёт точную частную сумму на каждом шаге.
    Итог — наружу на сетку 2^−g, g = 8·terms: сертификат несёт точный
    Q-отрезок [min(S_m, S_{m+1}); max(S_m, S_{m+1})], концы —
    на двоичной сетке (знаменатели — степени двойки).
    """
    if isinstance(nu, bool):
        raise TypeError(f"ν — {_NOT_BOOL}; дверь принимает int")
    if isinstance(nu, float):
        raise TypeError("ν — только int/Q, float запрещён")
    if not isinstance(nu, _NUM):
        raise TypeError("ν — только int/Q")
    if isinstance(x, float):
        raise TypeError("x — только int/Q, float запрещён")
    if not isinstance(x, _NUM):
        raise TypeError("x — только int/Q")
    if isinstance(terms, bool):
        raise TypeError(f"terms — {_NOT_BOOL}; дверь принимает int")
    if not isinstance(terms, int):
        raise TypeError("terms — только int")
    if terms < 1:
        raise ValueError("terms ≥ 1")
    x = Q(x)
    if isinstance(nu, Q) and nu.denominator != 1:
        if nu.denominator == 2 and abs(nu) == Q(1, 2):
            return _half_order_enclosure(nu, x, terms)
        return _rational_order_enclosure(nu, x, terms)
    nu = int(nu)
    if nu < 0:
        # зеркало J_{−n} = (−1)^n·J_n [внешний источник: Абрамовиц—Стиган 9.1.5]
        mirrored = j_enclosure(-nu, x, terms)
        return mirrored if nu % 2 == 0 else -mirrored
    g = 8 * terms
    G = g + _GUARD
    p, q = x.numerator, x.denominator
    # t_0 = (x/2)^ν / ν! = p^ν / ((2q)^ν·ν!) — границы на масштабе 2^−G
    num0, den0 = p ** nu, (2 * q) ** nu * _fact(nu)
    t_lo, t_hi = (num0 << G) // den0, _cdiv(num0 << G, den0)
    s_lo, s_hi = t_lo, t_hi             # частная сумма S_0
    r2_n, r2_d = p * p, 4 * q * q       # (x/2)² = p²/(4q²) — пара целых
    m = 0
    while True:
        # t_{m+1} = −t_m·(x/2)²/((m+1)(m+ν+1)) — границы к краям
        nxt_lo, nxt_hi = _scale_pair(t_lo, t_hi, -r2_n,
                                     r2_d * (m + 1) * (m + nu + 1))
        n_lo, n_hi = s_lo + nxt_lo, s_hi + nxt_hi
        if m + 1 >= terms and r2_n < r2_d * (m + 2) * (m + nu + 2):
            # старт мажоранты обеспечен: |t_{m+2}| < |t_{m+1}| и далее
            # монотонно — потому J_ν(x) лежит между S_m и S_{m+1};
            # берём внешние границы обеих сумм и выходим на сетку 2^−g
            return _grid_out(min(s_lo, n_lo), max(s_hi, n_hi), g)
        s_lo, s_hi, t_lo, t_hi, m = n_lo, n_hi, nxt_lo, nxt_hi, m + 1


def _exact_trig_enclosure(v, bits: int) -> Interval:
    """Отрезок точной величины cos πr / sin πr (класс 1 — точка; класс 3 —
    изолирующий отрезок носителя, суженный уточнением до 2^−bits)."""
    if isinstance(v, _NUM):
        return Interval.point(Q(v))
    refine = getattr(v, "refine", None)
    if refine is None or not hasattr(v, "lo") or not hasattr(v, "hi"):
        raise Unresolved("cos πν", _rec(v), hint="носитель без изолирующего отрезка")
    target = Q(1, 2 ** bits)
    while v.hi - v.lo > target:
        refine()
    return Interval(v.lo, v.hi)


def j_pair_enclosure(nu, x, terms: int = 24) -> "tuple[Interval, Interval]":
    """Пара отрезков (re, im) ⊂ ℚ² для J_ν(x), x ∈ ℚ ЛЮБОГО знака (07.10, круг
    седьмой) — двойственность класса 16 с парой (re, im) класса 4.

    x ≥ 0 — (J_ν(x), 0) вещественной дверью `j_enclosure`. x < 0 при нецелом ν:
    J_ν(x) = e^{iπν}·J_ν(|x|) — главная ветвь (x/2)^ν при arg x = π [внешний
    источник: DLMF 10.11.1 при m = 1, https://dlmf.nist.gov/10.11.E1], то есть
    re = cos(πν)·J_ν(|x|), im = sin(πν)·J_ν(|x|); cos πν и sin πν при ν ∈ ℚ —
    точные величины класса 1 либо 3 (`cyclotomic.cos_pi`, `sin_pi`: корни
    многочленов Ψ, без π и без ряда), их отрезки — уточнением носителя до
    2^−(8·terms); J_ν(|x|) — сертификат ряда, спуска или замкнутой формы. При
    целом ν зеркало J_n(−x) = (−1)^n·J_n(x) даёт im = 0 точно (sin πn = 0),
    при полуцелом ν — re = 0 точно (cos πν = 0): пара приземляется на ось.
    Вещественная дверь `j_enclosure` при x < 0 и нецелом ν по-прежнему
    отказывает по имени «вне ℝ» — она несёт ℝ, пара — эта дверь. x = 0 при
    ν < 0 — полюс, отказ той же двери."""
    if isinstance(nu, bool) or isinstance(x, bool):
        raise TypeError(f"ν, x — {_NOT_BOOL}; дверь принимает int/Q")
    if isinstance(nu, float) or isinstance(x, float):
        raise TypeError("ν, x — только int/Q, float запрещён")
    if not isinstance(nu, _NUM) or not isinstance(x, _NUM):
        raise TypeError("ν, x — только int/Q")
    x = Q(x)
    nu_q = Q(nu)
    if x >= 0 or nu_q.denominator == 1:
        return j_enclosure(nu, x, terms), Interval.point(0)
    from .cyclotomic import cos_pi, sin_pi
    mag = j_enclosure(nu_q, -x, terms)
    bits = 8 * terms
    c = _exact_trig_enclosure(cos_pi(nu_q), bits)
    sn = _exact_trig_enclosure(sin_pi(nu_q), bits)
    return mag * c, mag * sn


def _certified_sign(nu, x: Q, terms: int) -> "tuple[int, int]":
    """Знак J_ν(x) по сертификату: (знак ±1, достигнутая глубина terms).

    0 внутри сертификата — знак не разрешён: terms углубляется вдвое.
    Рациональный x не есть корень J_ν целого или рационального порядка
    (корни трансцендентны [внешний источник: теорема Зигеля 1929]),
    потому углубление завершается; потолок _MAX_TERMS даёт честный
    отказ, не догадку о знаке.
    """
    while True:
        iv = j_enclosure(nu, x, terms)
        if iv.lo > 0:
            return 1, terms
        if iv.hi < 0:
            return -1, terms
        if terms >= _MAX_TERMS:
            raise Unresolved("sign(J_ν)", _rec(x),
                             hint=f"знак не сертифицирован до глубины "
                                  f"{terms} — дальнейшее уже́ние отклонено")
        terms *= 2


def _bound2(nu: Q) -> Q:
    """Квадрат лейбницевой границы первого корня 4(M+1)(ν+M+1), M = ⌊−ν⌋ (M = 0
    при −1 < ν < −1/2): левее 2√((M+1)(ν+M+1)) корней J_ν нет (`_scan_start`)."""
    big_m = (-nu.numerator) // nu.denominator if nu < -1 else 0
    return 4 * (big_m + 1) * (nu + big_m + 1)


def _scan_step(nu) -> Q:
    """Шаг скана корней: 1/8, а при границе 2√((M+1)(ν+M+1)) ≤ 1/8 (ν → −(M+1)) —
    первая степень двойки ниже границы, чтобы узел под границей был (07.10,
    круг шестой; до него такие ν отвечали отказом «старт скана не отделён от
    ноля»). Цена — число узлов до первого корня, а после него шаг снова 1/8
    (разнос корней > π). Целые и ν ≥ −1/2 — 1/8."""
    nu = Q(nu)
    if nu >= Q(-1, 2) or nu.denominator == 1:
        return Q(1, 8)
    b2 = _bound2(nu)
    step = Q(1, 8)
    while step * step >= b2:
        step = step / 2
    return step


def _scan_start(nu) -> "Q | None":
    """Первый узел скана корней J_ν — левее него корней нет.

    ν ≥ −1/2: узел 1/8 — j_{ν,1} растёт по ν [внешний источник: DLMF 10.21(iv)]
    и при ν = −1/2 равен π/2 > 1/8. −1 < ν < −1/2 (02.10, круг четвёртый):
    лейбницева нижняя граница — J_ν(x) = (x/2)^ν·Σ(−1)^m·u_m, u_{m+1}/u_m =
    (x²/4)/((m+1)(m+ν+1)) ≤ (x²/4)/(ν+1) < 1 при x² < 4(ν+1), все u_m > 0
    (Γ(m+ν+1) > 0 при ν > −1), члены строго убывают — сумма > u_0 − u_1 > 0:
    корней на (0; 2√(ν+1)) нет [вычислено здесь]. Граница сравнима с корнем:
    j_{−3/4,1} = 1.0585… при 2√(1/4) = 1, j_{−254/255,1} = 0.12537… при
    2√(1/255) = 0.12524… [внешний источник: mpmath 1.3, besselj + findroot].

    Нецелое ν < −1, ν ∈ (−M−1; −M), M = ⌊−ν⌋ ≥ 1 (07.10, круг пятый) — та же
    граница с обобщением [вычислено здесь]: u_m = r^m/(m!·Γ(m+ν+1)), r = x²/4;
    при m ≥ M аргумент Γ положителен, при m < M он лежит в (m−M; m−M+1), и
    знак Γ есть (−1)^{M−m} — потому член (−1)^m·u_m имеет знак (−1)^M при
    m < M и (−1)^m при m ≥ M. Значит (−1)^M·Σ = H + (a_M − a_{M+1} + …),
    H = Σ_{m<M}|u_m| > 0, a_m = |u_m|; при r < (M+1)(M+ν+1) отношение
    a_{m+1}/a_m = r/((m+1)(m+ν+1)) < 1 для всех m ≥ M (знаменатель растёт по
    m), хвост лейбницев и лежит в (0; a_M), и (−1)^M·Σ > H > 0: корней на
    (0; 2√((M+1)(ν+M+1))) нет. При M = 0 это прежняя граница 2√(ν+1).

    Старт — наибольший узел n·h, h = `_scan_step(ν)`, с (n·h)² < 4(M+1)(ν+M+1);
    h выбран так, что h² < границы², потому n ≥ 1 всегда (07.10, круг шестой —
    до него при границе ≤ 1/8 стоял отказ «старт не отделён от ноля»: ν ≤
    −255/256 при M = 0, ν + M + 1 ≤ 1/(256·(M+1)) при M ≥ 1). Граница сравнима
    с корнем: j_{−3/2,1} = 2.7984… при 2; j_{−255/256,1} = 0.12512… при 1/8;
    j_{−1023/1024,1} = 0.062515… при 1/16 [внешний источник: mpmath 1.3].
    Целые ν < 0 сюда не приходят (зеркало у фабрики) — None.
    """
    nu = Q(nu)
    if nu >= Q(-1, 2):
        return Q(1, 8)
    if nu.denominator == 1:
        return None
    step = _scan_step(nu)
    t = _bound2(nu) / (step * step)           # n² < t, t > 1 по выбору шага
    top = -((-t.numerator) // t.denominator)  # ceil(t)
    n = math.isqrt(top - 1)                   # наибольшее n с n² < t
    return n * step


# ---------- сама линейка ----------

class BesselRootValue:
    """k-й положительный корень J_ν: ν ≥ 0 (int) либо любое рациональное
    ν ∉ ℤ (с 02.10 — сертификат Γ-ряда при ν > −1; с 07.10, круг пятый, —
    спуск рекуррентой при ν < −1; с круга шестого — и ν у самой границы
    −(M+1): старт скана под лейбницевой границей 2√((M+1)(ν+M+1)) с шагом,
    выбранным ниже неё, `_scan_step`), k ≥ 1 (int). Отказов по ν у носителя
    нет; цена — число узлов до первого корня и терms сертификата."""

    __slots__ = ("nu", "k")

    def __init__(self, nu, k: int):
        if isinstance(nu, bool):
            raise TypeError(f"ν — {_NOT_BOOL}; дверь принимает int/Q")
        if isinstance(nu, float) or not isinstance(nu, _NUM):
            raise TypeError("ν — только int/Q, float запрещён")
        if isinstance(k, bool):
            raise TypeError(f"k — {_NOT_BOOL}; дверь принимает int")
        if not isinstance(k, int):
            raise TypeError("k — только int")
        nu = Q(nu)
        if nu.denominator == 1:
            nu = int(nu)
            if nu < 0:
                raise ValueError("линейка несёт целые ν ≥ 0; отрицательный "
                                 "целый порядок — зеркалом у фабрики bessel_root")
        if k < 1:
            raise ValueError("нумерация положительных корней начинается "
                             "с k=1")
        self.nu, self.k = nu, k

    # --- арифметика: направленного правила нет, кроме ноля-состояния ---

    def _refuse(self, op: str, other: object, self_left: bool):
        """Чужой класс с записью — дверь ФIOS (порядок операндов как в
        выражении); класс 1, одноимённая линейка и прибитые ячейки —
        прежний отказ носителя."""
        if _foreign_record(other) and not _pinned(other, {"×": "*"}.get(op, op)):
            a, b = (self, other) if self_left else (other, self)
            return _door(a, b, {"×": "*"}.get(op, op))
        left, right = ((self.record(), _rec(other)) if self_left
                       else (_rec(other), self.record()))
        raise Unresolved(op, left, right,
                         hint="корень J_ν независим: направленного "
                              "правила арифметики нет")

    def __add__(self, other):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__radd__
        return self._refuse("+", other, True)

    def __radd__(self, other):
        if isinstance(other, Zero):
            return NotImplemented
        return self._refuse("+", other, False)

    def __sub__(self, other):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__rsub__
        return self._refuse("-", other, True)

    def __rsub__(self, other):
        if isinstance(other, Zero):
            return NotImplemented
        return self._refuse("-", other, False)

    def __mul__(self, other):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__rmul__
        return self._refuse("×", other, True)

    def __rmul__(self, other):
        if isinstance(other, Zero):
            return NotImplemented
        return self._refuse("×", other, False)

    def __truediv__(self, other):
        return self._refuse("/", other, True)

    def __rtruediv__(self, other):
        return self._refuse("/", other, False)

    def __pow__(self, n: int):
        if not isinstance(n, int) or isinstance(n, bool):
            return NotImplemented
        if n == 0:
            return Q(1)  # корень строго положителен: x⁰ = 1 точно
        if n == 1:
            return BesselRootValue(self.nu, self.k)
        raise Unresolved("^", self.record(), _rec(n),
                         hint="степень корня J_ν вне {0, 1} не сведена")

    def __neg__(self):
        raise Unresolved("neg", self.record(),
                         hint="линейка несёт k-й ПОЛОЖИТЕЛЬНЫЙ корень; "
                              "противоположной величины в ней нет")

    # --- запись, сравнение, сертификат ---------------------------------

    def record(self) -> str:
        return f"Φ(C=16 BESSEL_ROOT, ν={self.nu}, k={self.k})"

    __repr__ = record

    def __eq__(self, other):
        if isinstance(other, BesselRootValue):
            return self.nu == other.nu and self.k == other.k
        if isinstance(other, _NUM):
            # положительные корни J_ν целого и рационального порядка
            # трансцендентны [внешний источник: теорема Зигеля 1929 о
            # J_ν при рациональном ν] — рациональному не равны
            return False
        return NotImplemented

    def __hash__(self):
        return hash((BesselRootValue, self.nu, self.k))

    def self_check(self, deep: bool = False) -> list:
        """Самопроверка носителя встречной дорогой (ось «AIAΣA спрашивает AIAΣA»,
        слово автора 04.09; очередь строка 18): возвращает список находок,
        пусто — здоров; исключений не бросает.

        `[вычислено здесь]` Три дороги, ни одна не повторяет изоляцию корня:
        (1) слоты — ν рационален, k ≥ 1 целый, ν не есть ±1/2 (те приземлены в
        класс 10 фабрикой) и не отрицательное целое (зеркало); (2) сертификат
        enclosure(16) — на его концах знак J_ν сертифицирован и ПРОТИВОПОЛОЖЕН
        (корень внутри по промежуточному значению), ширина ≤ 2^−16, lo > 0;
        (3) старт скана `_scan_start` против теоремы: узел 1/8 при ν ≥ −1/2,
        наибольший узел под лейбницевой границей 2√((M+1)(ν+M+1)) при нецелом
        ν < −1/2 (арифметикой в ℚ), и строго левее сертификата: иначе k-я смена
        знака не есть k-й корень. При deep — четвёртая, по ряду: рекуррента
        J_{ν−1}(x) + J_{ν+1}(x) = (2ν/x)·J_ν(x) [внешний источник: А—С 9.1.27]
        в середине сертификата — отрезки обеих сторон обязаны пересекаться."""
        out = []
        nu, k = self.nu, self.k
        if isinstance(nu, bool) or not isinstance(nu, _NUM):
            return [f"слот ν не число: {nu!r}"]
        if isinstance(k, bool) or not isinstance(k, int) or k < 1:
            out.append(f"слот k не целое ≥ 1: {k!r}")
            return out
        nu_q = Q(nu)
        if nu_q.denominator == 2 and abs(nu_q) == Q(1, 2):
            out.append("ν = ±1/2 носителем: фабрика приземляет его в класс 10")
        if nu_q.denominator == 1 and nu_q < 0:
            out.append(f"отрицательный целый порядок {nu_q} носителем: зеркало у фабрики")
        try:
            iv = self.enclosure(16)
        except Exception as e:                                  # noqa: BLE001
            return out + [f"сертификат не построен: {type(e).__name__}: {e}"]
        if iv.lo <= 0:
            out.append(f"сертификат не положителен: lo = {iv.lo}")
        if iv.width() > Q(1, 2 ** 16):
            out.append(f"ширина сертификата {iv.width()} > 2^−16")
        try:
            s_lo, _ = _certified_sign(nu_q, iv.lo, 24)
            s_hi, _ = _certified_sign(nu_q, iv.hi, 24)
            if s_lo * s_hi != -1:
                out.append(f"знаки J_ν на концах сертификата не противоположны: {s_lo}, {s_hi}")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"знак на концах не сертифицирован: {type(e).__name__}: {e}")
        # старт скана — против теоремы, не против самого себя: при ν ≥ −1/2 старт
        # обязан быть узлом 1/8 (j_{ν,1} ≥ π/2 > 1/8, DLMF 10.21(iv)); при нецелом
        # ν < −1/2 — наибольшим узлом n·h под лейбницевой границей:
        # start² < 4(M+1)(ν+M+1) ≤ (start + h)², h = _scan_step(ν)
        start = _scan_start(nu_q)
        if nu_q >= Q(-1, 2):
            if start != Q(1, 8):
                out.append(f"старт скана {start} при ν ≥ −1/2 не есть узел 1/8")
        elif nu_q.denominator != 1:
            b2, h = _bound2(nu_q), _scan_step(nu_q)
            if start is None or not (start * start < b2 <= (start + h) ** 2):
                out.append(f"старт скана {start} вне лейбницевой границы: "
                           f"нужно start² < {b2} ≤ (start + {h})²")
        if start is not None and not start < iv.lo:
            out.append(f"старт скана {start} не левее сертификата [{iv.lo}; {iv.hi}]")
        if deep:
            mid = (iv.lo + iv.hi) / 2
            try:
                lhs = j_enclosure(nu_q - 1, mid, 24) + j_enclosure(nu_q + 1, mid, 24)
                rhs = j_enclosure(nu_q, mid, 24) * (2 * nu_q / mid)
                if lhs.intersect(rhs) is None:
                    out.append("рекуррента 9.1.27 в середине сертификата не сошлась")
            except Exception as e:                              # noqa: BLE001
                out.append(f"рекуррента в середине сертификата не построена: {type(e).__name__}: {e}")
        return out

    def enclosure(self, bits: int = 32) -> Interval:
        """Сертификат сетки [lo, hi] ⊂ ℚ ширины ≤ 2^−bits, по явному
        запросу.

        Изоляция: скан шага 1/8 по x от стартового узла (`_scan_start`:
        1/8 при ν ≥ −1/2; при −1 < ν < −1/2 — наибольший узел под
        лейбницевой границей 2√(ν+1), левее корней нет) до 4·(k+ν+2)
        собирает узлы с сертифицированным знаком J_ν (сертификат значения
        не несёт ноля на концах); k-я смена знака — стартовая вилка: между
        узлами со строгими противоположными знаками лежит нечётное
        число корней, а соседние корни J_ν разнесены не меньше чем на
        j_{ν,2} − j_{ν,1} ≥ π при |ν| ≥ 1/2 и > 3 при |ν| < 1/2 (ν ≥ −1/2,
        DLMF 10.21(ii)) ≫ 1/8 — значит ровно один; с 02.10 то же для
        рационального ν ≥ −1/2 — знак берётся сертификатом Γ-ряда; при
        −1 < ν < −1/2 (круг четвёртый) ν² > 1/4, и в нормальной форме
        u = √x·J_ν, u'' + (1 − (ν² − 1/4)/x²)·u = 0 [внешний источник:
        DLMF 10.13.1] коэффициент меньше 1 — по сравнению Штурма с
        v'' + v = 0 между соседними корнями u лежит корень всякого
        sin(x − a), то есть разнос > π [вычислено здесь]; при нецелом
        ν < −1 (круг пятый 07.10) ν² > 1 > 1/4 — тот же довод, старт — узел
        под обобщённой лейбницевой границей 2√((M+1)(ν+M+1)), знак —
        сертификатом спуска рекуррентой. У самой границы −(M+1) (круг
        шестой) граница ниже 1/8, и до первого корня шаг — степень двойки
        под ней (`_scan_step`), после него снова 1/8: разнос корней > π.
        Верхний предел скана 4·(k+|ν|+2)
        — не теорема, а честная граница: если смен знака не хватило,
        отказ «вилка не изолирована», не догадка. Дальше деление пополам
        по сертификату смены знака; конец с неразрешимым знаком (0 внутри
        сертификата) углубляет terms.
        """
        if isinstance(bits, bool):
            raise TypeError(f"bits — {_NOT_BOOL}; дверь принимает int")
        if not isinstance(bits, int):
            raise TypeError("bits — только int")
        if bits < 1:
            raise ValueError("bits ≥ 1")
        nu, k = self.nu, self.k
        step = _scan_step(nu)
        limit = 4 * (k + abs(nu) + 2)
        terms = 24
        changes = 0
        prev_x = prev_sign = None
        lo = hi = s_lo = None
        x = _scan_start(nu)
        while x <= limit:
            sign, terms = _certified_sign(nu, x, terms)
            if prev_sign is not None and sign != prev_sign:
                changes += 1
                if changes == k:
                    lo, s_lo, hi = prev_x, prev_sign, x
                    break
                step = Q(1, 8)                 # первый корень найден — дальше разнос > π
            prev_x, prev_sign = x, sign
            x += step
        if lo is None:
            raise Unresolved("enclosure", self.record(),
                             hint=f"скан шага 1/8 до {limit} дал только "
                                  f"{changes} смен знака — вилка k-го "
                                  "корня не изолирована")
        target = Q(1, 2 ** bits)
        while hi - lo > target:
            mid = (lo + hi) / 2
            s_mid, terms = _certified_sign(nu, mid, terms)
            if s_mid == s_lo:
                lo = mid
            else:
                hi = mid
        return Interval(lo, hi)


# ---------- самопроверка линейки (ось самоприменения) ----------

def self_check(deep: bool = False) -> list:
    """Самопроверка линейки класса 16 встречными дорогами (слово автора 04.09:
    «AIAΣA спрашивает AIAΣA»; очередь строка 18). Возвращает список находок,
    пусто — здорова; исключений не бросает.

    `[вычислено здесь]` Окно: ν ∈ {1/3, −3/4, 3/2, −4/3}, x ∈ {1, 2} (при deep —
    и ν ∈ {5/2, −7/3, −255/256}, x ∈ {1/2, 3}). Дороги: (1) рекуррента
    J_{ν−1}(x) + J_{ν+1}(x) = (2ν/x)·J_ν(x) [внешний источник: А—С 9.1.27] —
    отрезки сторон пересекаются; (2) две дороги степени (x/2)^ν — целый корень
    и exp∘ln — пересекаются (`_power_enclosure`); (3) пара при x < 0:
    re² + im² ∋ J_ν(|x|)² (DLMF 10.11.1); (4) полуцелые приземления:
    bessel_root(1/2, k) ∋ k·π, bessel_root(−1/2, k) ∋ (2k−1)·π/2 сертификатом
    класса 10 против скана носителем BesselRootValue(±1/2, k); (5) носитель
    первого корня каждого ν окна — `BesselRootValue.self_check`. Граница —
    окно: молчание на непроверенных ν находкой не считается."""
    out = []
    nus = [Q(1, 3), Q(-3, 4), Q(3, 2), Q(-4, 3)]
    xs = [Q(1), Q(2)]
    if deep:
        nus += [Q(5, 2), Q(-7, 3), Q(-255, 256)]
        xs += [Q(1, 2), Q(3)]
    for nu in nus:
        for x in xs:
            try:
                lhs = j_enclosure(nu - 1, x, 24) + j_enclosure(nu + 1, x, 24)
                rhs = j_enclosure(nu, x, 24) * (2 * nu / x)
                if lhs.intersect(rhs) is None:
                    out.append(f"рекуррента 9.1.27 при ν = {nu}, x = {x} не сошлась")
            except Exception as e:                              # noqa: BLE001
                out.append(f"рекуррента при ν = {nu}, x = {x}: {type(e).__name__}: {e}")
            try:
                a = _power_enclosure(x / 2, nu, 24, "root")
                b = _power_enclosure(x / 2, nu, 24, "expln")
                if a.intersect(b) is None:
                    out.append(f"дороги степени (x/2)^ν при ν = {nu}, x = {x} не пересекаются")
            except Exception as e:                              # noqa: BLE001
                out.append(f"дороги степени при ν = {nu}, x = {x}: {type(e).__name__}: {e}")
            try:
                re, im = j_pair_enclosure(nu, -x, 24)
                mag = j_enclosure(nu, x, 24)
                if (re * re + im * im).intersect(mag * mag) is None:
                    out.append(f"пара при ν = {nu}, x = {-x}: re² + im² не несёт J_ν(|x|)²")
            except Exception as e:                              # noqa: BLE001
                out.append(f"пара при ν = {nu}, x = {-x}: {type(e).__name__}: {e}")
        try:
            out += [f"носитель j_{{{nu},1}}: {f}" for f in BesselRootValue(nu, 1).self_check(deep)]
        except Exception as e:                                  # noqa: BLE001
            out.append(f"носитель j_{{{nu},1}}: {type(e).__name__}: {e}")
    for k in (1, 2):
        try:
            pi_k = pi_enclosure(40) * k
            if BesselRootValue(Q(1, 2), k).enclosure(16).intersect(pi_k) is None:
                out.append(f"скан носителем при ν = 1/2, k = {k} не несёт k·π класса 10")
            half = pi_enclosure(40) * Q(2 * k - 1, 2)
            if BesselRootValue(Q(-1, 2), k).enclosure(16).intersect(half) is None:
                out.append(f"скан носителем при ν = −1/2, k = {k} не несёт (2k−1)·π/2 класса 10")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"полуцелые приземления при k = {k}: {type(e).__name__}: {e}")
    return out


# ---------- фабрика линейки ----------

def bessel_root(nu, k: int):
    """k-й положительный корень J_ν по направленным правилам линейки.

    ν = 1/2  → PiMonomial(k, 1): J_{1/2}(x) = √(2/(π·x))·sin x, корни
               в k·π [внешний источник] — приземление в класс 10;
    ν = −1/2 → PiMonomial((2k−1)/2, 1): J_{−1/2}(x) = √(2/(π·x))·cos x,
               корни в (2k−1)·π/2 [внешний источник];
    ν ∈ ℤ, ν ≥ 0 → BesselRootValue(ν, k) — структурный носитель;
    ν ∈ ℤ, ν < 0 → тот же носитель с |ν|: J_{−n} = (−1)^n·J_n, корни
               общие [внешний источник: Абрамовиц—Стиган 9.1.5] (02.10);
    ν ∈ ℚ∖ℤ, ν > −255/256 → BesselRootValue(ν, k): изоляция знаком
               сертификата J_ν рационального порядка (02.10); при
               −1 < ν < −1/2 старт скана — лейбницева граница 2√(ν+1)
               (круг четвёртый 02.10, `_scan_start`);
    ν ∈ ℚ∖ℤ, ν < −1 → BesselRootValue(ν, k) (07.10, круг пятый): значение
               спуском рекуррентой по порядку, старт скана — обобщённая
               лейбницева граница 2√((M+1)(ν+M+1)), M = ⌊−ν⌋; у самой границы
               −(M+1) (круг шестой) шаг скана до первого корня — степень
               двойки под границей, отказа по ν больше нет; (x/2)^ν при
               знаменателе ν выше 64 — композицией exp∘ln, цена от
               знаменателя не зависит.
    """
    if isinstance(nu, float) or isinstance(k, float):
        raise TypeError("float запрещён построением")
    if isinstance(k, bool):
        raise TypeError(f"k — {_NOT_BOOL}; дверь принимает int")
    if not isinstance(k, int):
        raise TypeError("k — только int")
    if k < 1:
        raise ValueError("нумерация положительных корней начинается с k=1")
    if isinstance(nu, bool):
        raise TypeError(f"ν — {_NOT_BOOL}; дверь принимает int/Q")
    if not isinstance(nu, _NUM):
        raise TypeError("ν — только int/Q")
    nu = Q(nu)
    if nu == Q(1, 2):
        # ДВОЙСТВЕННОСТЬ J_{1/2} ↔ класс 10: корни синуса
        return PiMonomial(Q(k), Q(1))
    if nu == Q(-1, 2):
        # та же двойственность: корни косинуса
        return PiMonomial(Q(2 * k - 1, 2), Q(1))
    if nu.denominator == 1:
        # отрицательный целый порядок — зеркало: корни J_{−n} суть корни J_n
        return BesselRootValue(abs(int(nu)), k)
    return BesselRootValue(nu, k)
