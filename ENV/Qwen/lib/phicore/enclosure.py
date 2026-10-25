"""Сертификаты сетки: точные рациональные интервалы-вложения.

ОП-дискретность в устройстве: трансцендентная величина не «вычисляется
приближённо», а СЕРТИФИЦИРУЕТСЯ отрезком [lo, hi] — целым счётом на
выбранной ДВОИЧНОЙ сетке; уточнение — по явному запросу (уже́ние отрезка),
не скрытое округление. Ни одного float.

Концы сертификатов элементарных величин живут на сетке 2^−g: каждый конец —
пара целых вида m·2^−g (форма оси экспоненты контейнера v6 при основании 2;
Q здесь лишь обёртка совместимости, знаменатель — степень двойки).
Округление к сетке — всегда НАРУЖУ (низ — floor, верх — ceil), потому
вложение эталона сохраняется, а цена — не более одного шага сетки на конец.

Внутренняя арифметика рядов — целые аккумуляторы на общем масштабе
2^−(g+запас): масштабированное целое s ≈ x·2^G, произведения со сдвигами,
деления с округлением к обоим краям (вниз для нижней границы, вверх для
верхней). Дробей в горячих циклах нет. Оценки хвостов — строгие:
знакочередование (частные суммы обрамляют) либо мажоранта геометрической
прогрессией; обе — теоремы анализа, применённые в целой арифметике.
"""
from __future__ import annotations

from .rational import Q
from math import isqrt

from .errors import Unresolved

_NUM = (int, Q)


class Interval:
    """Отрезок [lo, hi] ⊂ ℚ, несущий величину; lo ≤ hi, оба Q."""

    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi):
        lo, hi = Q(lo), Q(hi)
        if lo > hi:
            raise ValueError(f"пустой отрезок: [{lo}; {hi}]")
        self.lo, self.hi = lo, hi

    # --- устройство ---
    def width(self) -> Q:
        return self.hi - self.lo

    def contains(self, x) -> bool:
        x = Q(x)
        return self.lo <= x <= self.hi

    def encloses(self, other: "Interval") -> bool:
        return self.lo <= other.lo and other.hi <= self.hi

    def intersect(self, other: "Interval"):
        lo, hi = max(self.lo, other.lo), min(self.hi, other.hi)
        return None if lo > hi else Interval(lo, hi)

    def record(self) -> str:
        return f"⟦{self.lo}; {self.hi}⟧"

    __repr__ = record

    def __eq__(self, other):
        return isinstance(other, Interval) \
            and (self.lo, self.hi) == (other.lo, other.hi)

    def __hash__(self):
        return hash((Interval, self.lo, self.hi))

    # --- самопроверка (ось «AIAΣA спрашивает AIAΣA», 04.09) ---
    def self_check(self, deep: bool = False) -> list:
        """Самопроверка отрезка ВСТРЕЧНОЙ дорогой (не пересчётом тем же кодом).

        `[вычислено здесь]` Дороги: (1) ширина — разностью концов здесь
        заново против `width()`; (2) принадлежность — сравнением ДРОБЕЙ
        (перекрёстным умножением внутри Q), а не сравнением концов между
        собой: концы и середина обязаны лежать внутри, а точки за концами —
        снаружи; (3) сетка: округление НАРУЖУ `on_grid(lo, hi, g)` обязано
        нести этот отрезок при любом g, а при g = `grid` (оба конца
        двоично-рациональны) округление точно и обязано дать тот же
        отрезок — это встречная дорога целыми к записи парами; (4)
        `encloses` и `intersect` сами с собой — тождества.

        Граница объявлена: арифметика отрезков монотонна по включению, но
        НЕ точна — зависимость операндов даёт избыток, поэтому `deep`
        сверяет ВКЛЮЧЕНИЕ ((I + J) − J ⊇ I, I² ⊆ I·I), а не равенство:
        равенства там нет и требовать его было бы ошибкой мерки. У отрезка
        с недвоичным концом `grid` равен None, и дорога (3) идёт одним
        включением. Возвращает список находок; исключений не бросает.
        """
        out = []
        if not isinstance(self.lo, Q) or not isinstance(self.hi, Q):
            out.append(f"концы не пары: {self.lo!r}, {self.hi!r}")
            return out
        if self.lo > self.hi:
            out.append(f"концы переставлены: {self.lo} > {self.hi}")
            return out
        if self.width() != self.hi - self.lo:
            out.append(f"ширина {self.width()} не равна разности концов "
                       f"{self.hi - self.lo}")
        mid = (self.lo + self.hi) / Q(2)
        for x in (self.lo, self.hi, mid):
            if not self.contains(x):
                out.append(f"величина внутри отрезка не признана своей: {x}")
        step = self.width() + Q(1)
        for x in (self.lo - step, self.hi + step):
            if self.contains(x):
                out.append(f"чужая величина признана своей: {x}")
        if not self.encloses(self):
            out.append("отрезок не несёт сам себя")
        if self.intersect(self) != self:
            out.append("пересечение с собой дало иное")
        for g in (0, 4, 10):
            try:
                wide = Interval.on_grid(self.lo, self.hi, g)
            except Exception as e:                              # noqa: BLE001
                out.append(f"округление на сетку {g} отказало: "
                           f"{type(e).__name__}: {e}")
                continue
            if not wide.encloses(self):
                out.append(f"округление наружу на сетку 2^−{g} потеряло "
                           f"отрезок: {wide.record()} ⊉ {self.record()}")
        g = self.grid
        if g is not None:
            exact = Interval.on_grid(self.lo, self.hi, g)
            if exact != self:
                out.append(f"на своей сетке 2^−{g} округление дало иное: "
                           f"{exact.record()}")
        if not deep:
            return out
        other = Interval(Q(1, 3), Q(2))
        if not ((self + other) - other).encloses(self):
            out.append("(I + J) − J не несёт I: арифметика не монотонна "
                       "по включению")
        if not (self * self).encloses(self ** 2):
            out.append("I² не вложен в I·I")
        if not (self * Interval.point(Q(1))).encloses(self):
            out.append("умножение на точку 1 потеряло отрезок")
        return out

    # --- двоичная сетка ---
    @classmethod
    def on_grid(cls, lo, hi, g: int) -> "Interval":
        """Отрезок [lo, hi], округлённый НАРУЖУ к сетке 2^−g.

        Низ — floor, верх — ceil на целых: результат несёт исходный
        отрезок и шире его не более чем на 2·2^−g.
        """
        if not isinstance(g, int) or isinstance(g, bool):
            raise TypeError(f"показатель сетки g — целое, не {type(g).__name__}")
        if g < 0:
            raise ValueError(f"показатель сетки g ≥ 0, получено {g}")
        lo, hi = Q(lo), Q(hi)
        scale = 1 << g
        n_lo = (lo.numerator * scale) // lo.denominator
        n_hi = -((-hi.numerator * scale) // hi.denominator)
        return cls(Q(n_lo, scale), Q(n_hi, scale))

    @property
    def grid(self):
        """Показатель g тончайшей сетки 2^−g, несущей оба конца.

        None, если хотя бы один конец не двоично-рационален
        (знаменатель — не степень двойки)."""
        g = 0
        for q in (self.lo.denominator, self.hi.denominator):
            if q & (q - 1):
                return None
            g = max(g, q.bit_length() - 1)
        return g

    # --- арифметика вложений (монотонна по включению) ---
    def __add__(self, other):
        if isinstance(other, Interval):
            return Interval(self.lo + other.lo, self.hi + other.hi)
        if isinstance(other, _NUM):
            q = Q(other)
            return Interval(self.lo + q, self.hi + q)
        return NotImplemented

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        if isinstance(other, (Interval,) + _NUM):
            return self + (-other if isinstance(other, Interval)
                           else -Q(other))
        return NotImplemented

    def __rsub__(self, other):
        if isinstance(other, _NUM):
            return (-self) + Q(other)
        return NotImplemented

    def __mul__(self, other):
        if isinstance(other, Interval):
            c = (self.lo * other.lo, self.lo * other.hi,
                 self.hi * other.lo, self.hi * other.hi)
            return Interval(min(c), max(c))
        if isinstance(other, _NUM):
            q = Q(other)
            a, b = self.lo * q, self.hi * q
            return Interval(min(a, b), max(a, b))
        return NotImplemented

    __rmul__ = __mul__

    def inverse(self) -> "Interval":
        if self.lo <= 0 <= self.hi:
            raise Unresolved("1/x", self.record(),
                             hint="отрезок несёт ноль — обращение не сведено: "
                                  "1/x на нём неограничено; сузить надо сам "
                                  "отрезок, а у Interval уточнения нет — оно "
                                  "у того, кто отрезок построил (довод terms "
                                  "сертификата)")
        return Interval(1 / self.hi, 1 / self.lo)

    def __truediv__(self, other):
        if isinstance(other, Interval):
            return self * other.inverse()
        if isinstance(other, _NUM):
            q = Q(other)
            if q == 0:
                raise Unresolved("/", self.record(), "Φ(C=1, value=0)",
                                 hint="x/0 здесь не сведено")
            return self * (1 / q)
        return NotImplemented

    def __rtruediv__(self, other):
        """q / отрезок = q · (1/отрезок): парная дверь к `__truediv__` — до
        25.09 её не было, и `2 / Interval` падал TypeError интерпретатора
        (замер приёмки `test_phicore_bessel_dual.py`); отрезок с нулём —
        отказ по имени из `inverse`."""
        if isinstance(other, _NUM):
            return self.inverse() * Q(other)
        return NotImplemented

    def __pow__(self, n: int):
        if not isinstance(n, int) or isinstance(n, bool):
            return NotImplemented
        if n == 0:
            return Interval(1, 1)
        if n < 0:
            return (self ** (-n)).inverse()
        out = Interval(1, 1)
        for _ in range(n):
            out = out * self
        return out

    @classmethod
    def point(cls, x) -> "Interval":
        x = Q(x)
        return cls(x, x)


# ---------- целый счёт на сетке: направленные округления ----------

_GUARD = 64            # запасные биты внутреннего масштаба G = g + _GUARD:
                       # накопленная за цикл потеря направленных делений
                       # (≤ 1 младший бит на шаг) тонет в запасе и не
                       # выходит за один шаг итоговой сетки 2^−g
#: `[вычислено здесь]` Здесь стоял порог 10^−40 — десятичная величина
#: внутри двоичного устройства, и мотив её не назывался. Порог теперь в
#: единицах САМОЙ сетки: счёт идёт, пока член не станет меньше одного узла
#: рабочего масштаба 2^−G, то есть пока он различим. Десятичное число из
#: двоичного устройства ушло, а порог стал свойством сетки, а не вкуса.


def _cdiv(a: int, b: int) -> int:
    """Деление с округлением ВВЕРХ, b > 0 (низ даёт обычный floor //)."""
    return -((-a) // b)


def _scale_pair(lo: int, hi: int, num: int, den: int) -> tuple:
    """Границы [lo, hi]·(num/den) на том же масштабе, den > 0:
    низ — вниз, верх — вверх; при num < 0 края меняются местами."""
    if num >= 0:
        return (lo * num) // den, _cdiv(hi * num, den)
    return (hi * num) // den, _cdiv(lo * num, den)


def _grid_out(lo_int: int, hi_int: int, g: int) -> Interval:
    """Итог с внутреннего масштаба 2^−(g+_GUARD) наружу к сетке 2^−g."""
    return Interval(Q(lo_int >> _GUARD, 1 << g),
                    Q(_cdiv(hi_int, 1 << _GUARD), 1 << g))


# ---------- строгие вложения элементарных величин ----------

def _fact(n: int) -> int:
    out = 1
    for k in range(2, n + 1):
        out *= k
    return out


def exp_enclosure(x, terms: int = 24) -> Interval:
    """e^x, x ∈ ℚ: частная сумма + мажоранта хвоста геометрией.

    Хвост: |R_N| ≤ |x|^N/N! · 1/(1 − |x|/(N+1)) при |x| < N+1.
    Счёт — целыми на масштабе 2^−G; итог — на сетке 2^−g, g = 8·N.
    """
    x = Q(x)
    n = max(terms, 2)
    while abs(x) >= n:  # гарантия сходимости мажоранты
        n += 8
    g = 8 * n
    G = g + _GUARD
    p, q = x.numerator, x.denominator
    t_lo = t_hi = 1 << G          # член k = 0 равен 1 — на сетке точен
    s_lo = s_hi = 0
    for k in range(n):
        s_lo += t_lo
        s_hi += t_hi
        # член k+1 = член k · x/(k+1); границы — с округлением к краям
        t_lo, t_hi = _scale_pair(t_lo, t_hi, p, q * (k + 1))
    # |x^n/n!| ≤ t_max·2^−G; множитель геометрии 1/(1−|x|/n) = qn/(qn−|p|)
    t_max = max(abs(t_lo), abs(t_hi))
    tail = _cdiv(t_max * q * n, q * n - abs(p))
    # `[вычислено здесь]` Знак остатка есть знак ПЕРВОГО ОТБРОШЕННОГО члена
    # x^n/n!: члены убывают по модулю (|x| < n обеспечено подбором n выше),
    # значит ряд с этого места знакочередующийся либо знакопостоянный, и
    # остаток не выходит за первый отброшенный член по знаку.
    #
    # Прежде мажоранта вычиталась снизу И прибавлялась сверху — засечка
    # была рыхла вдвое при том, что одна её сторона ТОЧНА. Замер 05.09:
    # exp(1, 3) давало ⟦9/4 − 2^−24; 11/4 + 2^−24⟧ при частной сумме 5/2,
    # хотя при x > 0 все отброшенные члены положительны и низ равен 5/2.
    if t_lo >= 0:                     # остаток ≥ 0 — точен низ
        return _grid_out(s_lo, s_hi + tail, g)
    if t_hi <= 0:                     # остаток ≤ 0 — точен верх
        return _grid_out(s_lo - tail, s_hi, g)
    return _grid_out(s_lo - tail, s_hi + tail, g)


def _alternating(x: Q, first_pow: int, step: int, denom_rule,
                 terms: int) -> Interval:
    """Σ (−1)^k · x^{first+step·k}/denom(k) при убывающих членах:
    частные суммы обрамляют предел (теорема Лейбница).

    denom_rule возвращает ЦЕЛОЕ. Степени x и частные суммы — пары целых
    границ на масштабе 2^−G; итог — наружу к сетке 2^−g,
    g = max(8·terms, 160): 160 — пол разрешения итога, чтобы при малом
    `terms` сетка не оказалась грубее того, что ряд уже насчитал.
    """
    # Пол разрешения итога: при малом `terms` сетка 8·terms была бы грубее
    # того, что ряд уже насчитал, и точность терялась бы на выходе, а не
    # в счёте. 160 — этот пол; порог остановки к нему больше не привязан.
    g = max(8 * terms, 160)
    G = g + _GUARD
    p, q = x.numerator, x.denominator
    pw_n, pw_d = p ** first_pow, q ** first_pow
    xs_n, xs_d = p ** step, q ** step
    p_lo = (pw_n << G) // pw_d
    p_hi = _cdiv(pw_n << G, pw_d)
    s_lo = s_hi = 0
    k = 0
    while True:
        d = denom_rule(k)
        t_lo, t_hi = p_lo // d, _cdiv(p_hi, d)
        if k % 2 == 0:
            n_lo, n_hi = s_lo + t_lo, s_hi + t_hi
        else:
            n_lo, n_hi = s_lo - t_hi, s_hi - t_lo
        t_max = max(abs(t_lo), abs(t_hi))
        # Порог — один узел ИТОГОВОЙ сетки 2^−g: член мельче него на итоге
        # уже неразличим. В единицах рабочего масштаба 2^−G это 2^_GUARD.
        # Ноль здесь взять нельзя: границы округляются наружу, и верхняя
        # никогда не станет нулём — замер 05.09: счёт упирался в предел
        # k > terms + 4000 и шёл 15 с вместо миллисекунд.
        if (k >= terms and t_max < 1 << _GUARD) or k > terms + 4000:
            # предел — между соседними частными суммами; берём обе границы
            return _grid_out(min(s_lo, n_lo), max(s_hi, n_hi), g)
        s_lo, s_hi = n_lo, n_hi
        p_lo, p_hi = _scale_pair(p_lo, p_hi, xs_n, xs_d)
        k += 1


def _rotate_period(x: Q, terms: int) -> "tuple[int, Interval]":
    """Поворот периода: x = r + 2π·k, k = ⌊x/(2π) + 1/2⌋ по середине
    сертификата π; остаток r — ОТРЕЗОК [x − 2k·π_hi; x − 2k·π_lo] ⊂ ℚ
    (π иррационально, точкой остаток не записывается), |r| ≤ π + ширина.

    `[вычислено здесь]` 01.10 (очередь строка 104): до того sin/cos за
    |x| ≤ 4 отказывали по имени, и сертификаты корней J_{±1/2} при k ≥ 2
    (2π, 3π/2 — приёмка Бесселя D) были закрыты 1 из 12. Ширина остатка —
    2k·ширина(π): при terms = 20 ширина π ≈ 2⁻¹⁵⁵, при 30 — 2⁻²³⁶.
    """
    piv = pi_enclosure(max(30, terms))
    mid = (piv.lo + piv.hi) / 2
    t = x / (2 * mid)
    k = (2 * t.numerator + t.denominator) // (2 * t.denominator)   # ⌊t + 1/2⌋
    if k == 0:
        return 0, Interval.point(x)
    if k > 0:
        return k, Interval(x - 2 * k * piv.hi, x - 2 * k * piv.lo)
    return k, Interval(x - 2 * k * piv.lo, x - 2 * k * piv.hi)


def _rotated(x: Q, terms: int, direct) -> Interval:
    """sin/cos за сеткой ряда: значение в левом конце остатка рядом, плюс
    ширина остатка по Липшицу (|sin′|, |cos′| ≤ 1) — сертификат строгий,
    не приближение."""
    _, r = _rotate_period(x, terms)
    base = direct(r.lo, terms)
    w = r.width()
    return Interval(base.lo - w, base.hi + w)


def sin_enclosure(x, terms: int = 20) -> Interval:
    """sin x, x ∈ ℚ: при |x| ≤ 4 — знакочередующийся ряд, суммы обрамляют;
    дальше — поворот периода (`_rotate_period`, 01.10) и тот же ряд на
    остатке с поправкой Липшица.

    `[вычислено здесь]` Откуда 4. Обрамление по Лейбницу требует убывания
    членов по модулю. У ряда синуса отношение соседних членов есть
    x²/((2k+2)(2k+3)); при k ≥ 1 знаменатель не меньше 20, значит члены
    убывают начиная со второго при x² < 20, то есть при |x| < 4,47…
    Граница 4 — целое внутри этого условия (16 < 20), и она названа, а не
    выбрана на глаз. Дальше её сводит поворот периода: до 01.10 его не
    было, и дверь отказывала по имени (стена снята, пересмотр отказов).
    """
    x = Q(x)
    if abs(x) > 4:
        return _rotated(x, terms, sin_enclosure)
    return _alternating(x, 1, 2, lambda k: _fact(2 * k + 1), terms)


def cos_enclosure(x, terms: int = 20) -> Interval:
    """cos x, x ∈ ℚ: при |x| ≤ 4 — знакочередующийся ряд, суммы обрамляют;
    дальше — поворот периода (`_rotate_period`, 01.10) с поправкой Липшица.

    `[вычислено здесь]` Граница та же, что у синуса, и по той же причине:
    отношение соседних членов x²/((2k+1)(2k+2)) при k ≥ 1 не превосходит
    x²/12, а условие убывания x² < 20 покрывает 16 с запасом. Мотив
    назван, чтобы число в подписи не выглядело вкусом.
    """
    x = Q(x)
    if abs(x) > 4:
        return _rotated(x, terms, cos_enclosure)
    return _alternating(x, 0, 2, lambda k: _fact(2 * k), terms)


def _atan_enclosure(x: Q, terms: int) -> Interval:
    """arctan x при |x| < 1: знакочередование, суммы обрамляют."""
    return _alternating(x, 1, 2, lambda k: 2 * k + 1, terms)


def ln_enclosure(y, terms: int = 40) -> Interval:
    """ln y, y ∈ ℚ, y > 0: через 2·atanh t, t = (y−1)/(y+1) ∈ (−1, 1);
    хвост — геометрическая мажоранта. Счёт целый, итог на сетке 2^−g."""
    y = Q(y)
    if y <= 0:
        raise Unresolved("ln", f"Φ(C=1, value={y})",
                         hint="ln на неположительном — не веществен")
    if y == 1:
        return Interval.point(0)
    g = max(8 * terms, 160)
    G = g + _GUARD
    ta, tb = y.numerator - y.denominator, y.numerator + y.denominator
    t2_n, t2_d = ta * ta, tb * tb          # t² = ta²/tb² < 1 при y > 0
    p_lo = (ta << G) // tb
    p_hi = _cdiv(ta << G, tb)
    s_lo = s_hi = 0
    for k in range(terms):
        d = 2 * k + 1
        s_lo += p_lo // d
        s_hi += _cdiv(p_hi, d)
        p_lo, p_hi = _scale_pair(p_lo, p_hi, t2_n, t2_d)
    # хвост: |t^{2N+1}|/((2N+1)·(1−t²)); 1−t² = (tb²−ta²)/tb²
    p_max = max(abs(p_lo), abs(p_hi))
    tail = _cdiv(p_max * t2_d, (2 * terms + 1) * (t2_d - t2_n))
    return _grid_out(2 * (s_lo - tail), 2 * (s_hi + tail), g)


def pi_enclosure(terms: int = 30) -> Interval:
    """π по Мэчину: 16·atan(1/5) − 4·atan(1/239); оба ряда обрамляют."""
    a = _atan_enclosure(Q(1, 5), terms)
    b = _atan_enclosure(Q(1, 239), terms)
    return a * 16 - b * 4


def int_root_floor(t: int, q: int) -> int:
    """⌊t^(1/q)⌋ на целых (Ньютон сверху), t ≥ 0 — ни одного float."""
    if t < 0:
        raise ValueError("подкоренное отрицательно")
    if t == 0 or q == 1:
        return t
    x = 1 << ((t.bit_length() + q - 1) // q + 1)
    while True:
        y = ((q - 1) * x + t // x ** (q - 1)) // q
        if y >= x:
            break
        x = y
    while x ** q > t:
        x -= 1
    return x


def root_enclosure(iv: Interval, n: int, bits: int = 64) -> Interval:
    """Корень n-й степени из ОТРЕЗКА на сетке 2^−bits — монотонность корня:
    [lo, hi] ↦ [⌊lo^(1/n)⌋_g, ⌈hi^(1/n)⌉_g]. Ни одного float: целый корень
    из ⌊lo·2^(n·bits)⌋ снизу и из ⌈hi·2^(n·bits)⌉ сверху. Нечётный n
    переносится на отрицательные симметрией; чётный n при lo<0 — отказ.
    Это сертификат дробной степени носителя: α^(p/q) = (α^(1/q))^p.
    """
    if not isinstance(n, int) or isinstance(n, bool) or n < 1:
        raise TypeError("степень корня — целое ≥ 1")
    if n == 1:
        return iv
    if iv.lo < 0:
        if n % 2 == 0:
            raise Unresolved("root", iv.record(),
                             hint=f"чётный корень ({n}) из отрезка с "
                                  f"отрицательным концом не веществен")
        if iv.hi <= 0:
            return -root_enclosure(-iv, n, bits)
        neg = root_enclosure(Interval(-iv.lo, -iv.lo), n, bits)
        pos = root_enclosure(Interval(iv.hi, iv.hi), n, bits)
        return Interval(-neg.hi, pos.hi)
    scale = 1 << bits
    big = 1 << (n * bits)
    lo, hi = iv.lo, iv.hi
    t_lo = (lo.numerator * big) // lo.denominator
    r_lo = int_root_floor(t_lo, n)
    t_hi = -((-hi.numerator * big) // hi.denominator)
    r_hi = int_root_floor(t_hi, n)
    if not (r_hi ** n == t_hi and t_hi * hi.denominator == hi.numerator * big):
        r_hi += 1
    return Interval(Q(r_lo, scale), Q(r_hi, scale))


def sqrt_enclosure(y, bits: int = 64) -> Interval:
    """√y, y ∈ ℚ, y ≥ 0: целый корень на сетке 2^−bits.

    lo = ⌊√y·2^g⌋, hi = lo либо lo+1 (точный квадрат распознаётся);
    ширина ≤ 2^−bits — как у прежнего деления отрезка пополам, но за
    один isqrt вместо цикла.
    """
    y = Q(y)
    if y < 0:
        raise Unresolved("sqrt", f"Φ(C=1, value={y})",
                         hint="чётный корень из отрицательного не веществен")
    if y == 0:
        return Interval.point(0)
    g = bits
    scale = 1 << g
    a, b = y.numerator, y.denominator
    lo_n = isqrt((a << (2 * g)) // b)      # ⌊√⌊y·4^g⌋⌋ = ⌊√(y·4^g)⌋
    hi_n = lo_n if lo_n * lo_n * b == a << (2 * g) else lo_n + 1
    return Interval(Q(lo_n, scale), Q(hi_n, scale))
