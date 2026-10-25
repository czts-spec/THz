"""Линейка класса 14 GAMMA: канонический вид mult·Γ(base), base ∈ (0, 1).

Величина хранится структурно: точный множитель mult (Q) и базис
base — аргумент Γ, сведённый рекуррентой Γ(z) = (z−1)·Γ(z−1) в интервал
(0, 1); накопление mult — в обе стороны сдвига. Целый аргумент z ≥ 1
вырождается в точный Q((z−1)!), факториал считается рекуррентой;
целый z ≤ 0 — полюс, честный отказ.

ДВОЙСТВЕННОСТЬ с классом 10: отражение Эйлера Γ(b)·Γ(1−b) = π/sin(π·b)
приземляет произведение в моном PiMonomial; sin(π·b) рационального
аргумента алгебраичен и берётся точной таблицей знаменателей (sin_pi):
класс 1 для {1, 2, 6}, класс 3 (POLY_ROOT) для {3, 4}. Прочие сведения
направленных правил не имеют: сумма Γ с чем угодно, кроме
ноля-состояния, произведение и частное разных базисов, 1/Γ — отказ
Unresolved с явными Ф-записями операндов, не ложный класс.

Десятичная сторона — только сертификатом сетки (enclosure) по явному
запросу: ряд Стирлинга для ln Γ(z) на вещественном z > 0 с остатком,
ограниченным первым отброшенным членом и того же знака [внешний
источник: DLMF 5.11.1, 5.11(ii); Уиттекер–Ватсон 12.33] — аргумент
сдвигается вверх рекуррентой (z = base + n), ln z и ln 2π берутся
сертификатами с приведением аргумента к [1, 2), числа Бернулли — точные
пары из zeta.py; назад к base — делением на точное произведение и
exp-сертификатом. Ни одного float; ширина отрезка убывает с terms.
"""
from __future__ import annotations

from .rational import Q

from .algebraic import AlgebraicNumber, _MAX_REFINE
from .enclosure import Interval, exp_enclosure, ln_enclosure, pi_enclosure
from .errors import Unresolved
from .pi_value import PiMonomial
from .zero_state import Zero

#: Допустимые скаляры. float исключён построением.
_NUM = (int, Q)

#: `[вычислено здесь]` 10.10, очередь строка 105, шаг пятый: дверь оператора
#: носителя = дверь ФIOS, как у π, e^y/ln, arctrig и ζ. Для ЧУЖОГО класса с
#: записью операторы Γ зовут `dispatch.combine`; класс 1 (int, Q) и
#: одноимённая линейка — прежние правила (отражение Эйлера, сокращение
#: Γ(base)) и собственные отказы носителя. Ввоз ленивый, внутри функций.
_IN_DOOR: set = set()
#: Места `[вычислено здесь]` 10.10 — proto/test_phicore_gamma.py: (10, ×) — «π·Γ(1/2) —
#: отказ (правила сведения линеек нет)» (ячейка (14, ×) прибита и со стороны монома);
#: (10, +) — «Γ(1/2)+π — отказ». Обе ячейки дверью дали бы запись класса 19; снятие —
#: с пересмотром приёмок.
_PINNED: frozenset = frozenset({(10, "*"), (10, "+")})


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
    return (not isinstance(x, (_NUM, GammaValue))) and callable(getattr(x, "record", None))


def _rec(x: object) -> str:
    """Ф-запись операнда для сообщений отказа."""
    r = getattr(x, "record", None)
    if callable(r):
        return r()
    if isinstance(x, _NUM):
        return f"Φ(C=1 RATIONAL, value={Q(x)})"
    return f"Φ(?={x!r})"


def _factorial(n: int) -> int:
    """n! точной рекуррентой k! = k·(k−1)!; таблиц значений нет."""
    out = 1
    for k in range(2, n + 1):
        out *= k
    return out


def sin_pi(r) -> Q | AlgebraicNumber:
    """Точное значение sin(π·r), r ∈ ℚ.

    Приведение: период 2 по r (r берётся по модулю 2) и знак —
    sin(π(1+t)) = −sin(π·t), остаток t ∈ [0, 1). Дальше значение
    зависит только от знаменателя несократимого t:
      1 → 0;  2 → 1;  6 → 1/2 — рациональные значения
      [внешний источник, точная таблица синусов: sin(π/6) = 1/2];
      3 → √3/2;  4 → √2/2 — корни, линейка класса 3.
    Иные знаменатели — с 04.09 тоже класс 3: sin(πr) = cos(π(r − ½)) =
    2cos(2πk/N)/2, минимальный многочлен степени φ(N)/2 через многочлены
    Люка и факторизацию над ℤ (cyclotomic.py); отказа больше нет.
    """
    if isinstance(r, float):
        raise TypeError("r — только int/Q, float запрещён")
    if not isinstance(r, _NUM):
        raise TypeError("r — только int/Q")
    t = Q(r) % 2          # период sin(π·r) по r равен 2
    sign = 1
    if t >= 1:                   # sin(π(1+t)) = −sin(π·t)
        sign, t = -1, t - 1
    den = t.denominator          # t ∈ [0, 1), Q несократима
    if den == 1:                 # t = 0
        return Q(0)
    if den == 2:                 # t = 1/2
        return Q(sign)
    if den == 6:                 # t ∈ {1/6, 5/6}: sin = 1/2, класс 1
        return Q(sign, 2)
    if den == 3:                 # t ∈ {1/3, 2/3}: sin = √3/2
        root = AlgebraicNumber.from_root(2, 3) / 2
        return -root if sign < 0 else root
    if den == 4:                 # t ∈ {1/4, 3/4}: sin = √2/2
        root = AlgebraicNumber.from_root(2, 2) / 2
        return -root if sign < 0 else root
    from .cyclotomic import sin_pi as _sin_pi_cyclotomic   # вне таблицы: класс 3 (04.09)
    val = _sin_pi_cyclotomic(t)
    return -val if sign < 0 else val


def _inv_sin(s):
    """Точное 1/sin(π·b) для b ∈ (0, 1): величина строго положительна.

    Изолирующий отрезок корня уточняется до отделения от ноля —
    уже́ние отрезка есть точный шаг сетки, не округление; обращение
    после него остаётся точным.

    Предел уточнений — тот же, что у `AlgebraicNumber.inverse`
    (`_MAX_REFINE`): замер 25.09 (приёмка `test_phicore_gamma_dual.py`,
    подмены дверей) — без предела цикл на неположительной величине не
    завершался, и ворота на такой двери не краснели, а стояли. Величина
    ≤ 0 (b вне (0, 1) либо подменённая дверь синуса) — отказ по имени.
    """
    if isinstance(s, AlgebraicNumber):
        for _ in range(_MAX_REFINE):
            if s.lo > 0:
                return s.inverse()
            if s.hi <= 0:
                break
            s.refine()
        raise Unresolved("1/sin", s.record(),
                         hint="sin(π·b) не положителен либо отрезок не отделился от нуля: "
                              "b вне (0, 1) или подменённая дверь синуса")
    if s <= 0:
        raise Unresolved("1/sin", str(s), hint="sin(π·b) не положителен: b вне (0, 1)")
    return 1 / s


class GammaValue:
    """mult·Γ(base): mult ≠ 0 (Q), base ∈ (0, 1) (Q)."""

    __slots__ = ("mult", "base")

    def __init__(self, mult, base):
        if isinstance(mult, float) or isinstance(base, float):
            raise TypeError("float запрещён построением")
        if not isinstance(mult, _NUM) or not isinstance(base, _NUM):
            raise TypeError("mult и base — только int/Q")
        mult, base = Q(mult), Q(base)
        if mult == 0:
            raise ValueError("mult=0 не есть величина линейки — ноль "
                             "живёт в zero_state")
        if not 0 < base < 1:
            raise ValueError("канонический вид требует base ∈ (0, 1); "
                             "сдвиг аргумента — фабрикой of()")
        self.mult, self.base = mult, base

    # --- построение -----------------------------------------------------

    @classmethod
    def of(cls, z) -> Q | GammaValue:
        """Γ(z) по рациональному z.

        Целые z ≥ 1 вырождаются в точный Q((z−1)!); целые z ≤ 0 —
        полюс, отказ. Прочие сводятся рекуррентой Γ(z) = (z−1)·Γ(z−1)
        к base ∈ (0, 1) с накоплением точного mult в обе стороны.
        """
        if isinstance(z, float):
            raise TypeError("z — только int/Q, float запрещён")
        if not isinstance(z, _NUM):
            raise TypeError("z — только int/Q")
        z = Q(z)
        if z.denominator == 1:
            if z <= 0:
                raise Unresolved("Γ", _rec(z),
                                 hint="полюс Γ: неположительное целое, "
                                      "значения нет")
            return Q(_factorial(z.numerator - 1))
        mult = Q(1)
        while z > 1:             # спуск: Γ(z) = (z−1)·Γ(z−1)
            z -= 1
            mult *= z
        while z < 0:             # подъём: Γ(z) = Γ(z+1)/z
            mult /= z
            z += 1
        return cls(mult, z)

    # --- умножение / деление -------------------------------------------

    def __mul__(self, other):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__rmul__
        if isinstance(other, GammaValue):
            if self.base + other.base == 1:
                # ОТРАЖЕНИЕ ЭЙЛЕРА: Γ(b)·Γ(1−b) = π/sin(π·b) —
                # приземление в класс 10 (двойственность линеек);
                # случай base = 1/2 даёт Γ(1/2)² = π этим же правилом
                s = sin_pi(self.base)          # любой знаменатель: таблица либо класс 3 (04.09)
                coeff = _inv_sin(s) * (self.mult * other.mult)
                return PiMonomial(coeff, Q(1))
            if self.base == other.base:
                raise Unresolved("×", self.record(), other.record(),
                                 hint="Γ(b)·Γ(b) вне b=1/2 не сведено")
            raise Unresolved("×", self.record(), other.record(),
                             hint="произведение Γ разных базисов без "
                                  "отражения не сведено")
        if isinstance(other, _NUM):
            q = Q(other)
            if q == 0:
                return Q(0)  # масштаб нолём — точный ноль класса 1
            return GammaValue(self.mult * q, self.base)
        if _foreign_record(other) and not _pinned(other, "*"):
            return _door(self, other, "*")      # чужой класс — дверь ФIOS
        r = getattr(other, "record", None)
        if callable(r):
            raise Unresolved("×", self.record(), r(),
                             hint="правила сведения с этой линейкой нет")
        return NotImplemented

    def __rmul__(self, other):
        """Скаляр и одноимённая линейка — коммутативно, тем же правилом; чужой
        класс — дверь с порядком операндов как в выражении."""
        if _foreign_record(other) and not _pinned(other, "*"):
            return _door(other, self, "*")
        return self.__mul__(other)

    def __truediv__(self, other):
        if isinstance(other, GammaValue):
            if other.base == self.base:
                return self.mult / other.mult  # Γ(base) сокращается точно
            raise Unresolved("/", self.record(), other.record(),
                             hint="частное Γ разных базисов не сведено")
        if isinstance(other, _NUM):
            q = Q(other)
            if q == 0:
                raise Unresolved("/", self.record(), _rec(other),
                                 hint="x/0 здесь не сведено")
            return GammaValue(self.mult / q, self.base)
        if isinstance(other, Zero) or (_foreign_record(other) and not _pinned(other, "/")):
            return _door(self, other, "/")      # ноль-состояние — таблица «x/0»
        r = getattr(other, "record", None)
        if callable(r):
            raise Unresolved("/", self.record(), r(),
                             hint="правила сведения с этой линейкой нет")
        return NotImplemented

    def __rtruediv__(self, other):
        if isinstance(other, _NUM):
            raise Unresolved("/", _rec(other), self.record(),
                             hint="обратная величина 1/Γ этой линейкой "
                                  "не несётся")
        if _foreign_record(other) and not _pinned(other, "/"):
            return _door(other, self, "/")
        return NotImplemented

    def __pow__(self, n: int):
        if not isinstance(n, int) or isinstance(n, bool):
            return NotImplemented
        if n == 0:
            return Q(1)
        if n == 1:
            return GammaValue(self.mult, self.base)
        if n == 2:
            return self * self  # сведётся лишь отражением (base = 1/2)
        raise Unresolved("^", self.record(), _rec(n),
                         hint="степень Γ вне {0, 1, 2} не сведена")

    # --- сложение / вычитание: направленного правила нет ----------------

    def _refuse(self, op: str, other: object, self_left: bool):
        left, right = ((self.record(), _rec(other)) if self_left
                       else (_rec(other), self.record()))
        raise Unresolved(op, left, right,
                         hint="Γ независима: сложение сведено только "
                              "с нолём-состоянием")

    def __add__(self, other):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__radd__
        if _foreign_record(other) and not _pinned(other, "+"):
            return _door(self, other, "+")      # чужой класс — дверь ФIOS
        self._refuse("+", other, True)

    def __radd__(self, other):
        if _foreign_record(other) and not _pinned(other, "+"):
            return _door(other, self, "+")
        self._refuse("+", other, False)

    def __sub__(self, other):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__rsub__
        if _foreign_record(other) and not _pinned(other, "-"):
            return _door(self, other, "-")
        self._refuse("-", other, True)

    def __rsub__(self, other):
        if _foreign_record(other) and not _pinned(other, "-"):
            return _door(other, self, "-")
        self._refuse("-", other, False)

    def __neg__(self) -> GammaValue:
        return GammaValue(-self.mult, self.base)

    # --- запись, сравнение, сертификат ---------------------------------

    def record(self) -> str:
        return f"Φ(C=14 GAMMA, mult={self.mult}, base={self.base})"

    __repr__ = record

    def __eq__(self, other):
        if isinstance(other, GammaValue):
            return self.mult == other.mult and self.base == other.base
        # равенство Γ(base) с рациональным или алгебраическим этой
        # линейкой не решается (иррациональность Γ(b) в общем случае
        # не доказана) — ответ не выдумывается
        return NotImplemented

    def __hash__(self):
        return hash((GammaValue, self.mult, self.base))

    def enclosure(self, terms: int = 30) -> Interval:
        """Сертификат сетки [lo, hi] ⊂ ℚ: mult · вложение Γ(base) рядом
        Стирлинга с доказанным остатком (см. gamma_enclosure)."""
        return gamma_enclosure(self.base, terms) * self.mult

    def self_check(self, deep: bool = False) -> list:
        """Самопроверка носителя встречными дорогами (ось «AIAΣA спрашивает
        AIAΣA», слово автора 04.09; очередь строка 18; 07.10). Список находок,
        пусто — здоров; исключений не бросает.

        `[вычислено здесь]` (1) слоты: mult — ненулевое ℚ, base ∈ (0; 1) ∩ ℚ;
        (2) отражение Эйлера Γ(b)·Γ(1−b) = π/sin(πb) [внешний источник: DLMF
        5.5.3, https://dlmf.nist.gov/5.5.E3] — сертификат Γ(b) этого носителя,
        умноженный на НЕЗАВИСИМЫЙ сертификат Γ(1−b) (свой Стирлинг в другой
        точке), обязан пересекать π/sin(πb): π — Мэчин, sin πb — точная
        величина класса 1 либо 3 (`sin_pi`), ни одна сторона не повторяет
        другую; (3) при deep — рекуррента Γ(b+1) = b·Γ(b) [DLMF 5.5.1]
        сертификатами двух точек. Граница: сертификаты на 30 членах — ширина,
        а не теорема; пересечение пустых не бывает у здоровых точек окна."""
        out = []
        m, b = self.mult, self.base
        if isinstance(m, bool) or not isinstance(m, _NUM) or Q(m) == 0:
            out.append(f"слот mult не ненулевое число: {m!r}")
        if isinstance(b, bool) or not isinstance(b, _NUM) or not (0 < Q(b) < 1):
            out.append(f"слот base вне (0; 1): {b!r}")
            return out
        b = Q(b)
        g = None
        try:
            g = gamma_enclosure(b, 30)
            if g.lo <= 0:
                out.append(f"сертификат Γ({b}) не положителен: lo = {g.lo}")
            g1 = gamma_enclosure(1 - b, 30)
            sv = sin_pi(b)
            s_iv = Interval.point(sv) if isinstance(sv, Q) else sv.enclosure(60)
            if s_iv.lo <= 0:
                out.append(f"sin(π·{b}) не отделён от нуля сверху нуля: {s_iv.record()}")
            else:
                rhs = pi_enclosure(40) * s_iv.inverse()
                if (g * g1).intersect(rhs) is None:
                    out.append(f"отражение Эйлера при b = {b}: Γ(b)·Γ(1−b) не несёт π/sin(πb)")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"отражение при b = {b} не построено: {type(e).__name__}: {e}")
        if deep and g is not None:
            try:
                if gamma_enclosure(b + 1, 30).intersect(g * b) is None:
                    out.append(f"рекуррента Γ(b+1) = b·Γ(b) при b = {b} не сошлась")
            except Exception as e:                              # noqa: BLE001
                out.append(f"рекуррента при b = {b} не построена: {type(e).__name__}: {e}")
        return out


# ---------- самопроверка линейки (ось самоприменения, 07.10) ----------

def self_check(deep: bool = False) -> list:
    """Самопроверка линейки класса 14 встречными дорогами; список находок,
    пусто — здорова; исключений не бросает.

    `[вычислено здесь]` Окно базисов {1/3, 1/4, 1/6, 2/5} (deep — ещё 3/7, 1/12,
    7/8): носители проверяют себя отражением Эйлера и рекуррентой
    (`GammaValue.self_check`); Γ(1/2)² несёт π (самодуальная точка отражения,
    DLMF 5.4.6); целые аргументы — точный факториал фабрикой `of()`, и
    сертификат Γ(n+1) несёт n! (n ∈ {1, 2, 3, 5}); канон `of(b+2)` — слоты
    mult = (b+1)·b и base = b (рекуррента спуском). Граница — окно."""
    out = []
    bases = [Q(1, 3), Q(1, 4), Q(1, 6), Q(2, 5)]
    if deep:
        bases += [Q(3, 7), Q(1, 12), Q(7, 8)]
    for b in bases:
        try:
            out += [f"Γ({b}): {f}" for f in GammaValue(1, b).self_check(deep)]
        except Exception as e:                                  # noqa: BLE001
            out.append(f"Γ({b}): {type(e).__name__}: {e}")
    try:
        g = gamma_enclosure(Q(1, 2), 30)
        if (g * g).intersect(pi_enclosure(40)) is None:
            out.append("Γ(1/2)² не несёт π")
    except Exception as e:                                      # noqa: BLE001
        out.append(f"Γ(1/2)²: {type(e).__name__}: {e}")
    for n in (1, 2, 3, 5):
        try:
            exact = GammaValue.of(n + 1)
            if not isinstance(exact, Q) or exact != _factorial(n):
                out.append(f"Γ({n + 1}) фабрикой не есть {n}!: {exact!r}")
            elif not gamma_enclosure(Q(n + 1), 30).contains(exact):
                out.append(f"сертификат Γ({n + 1}) не несёт {n}!")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"Γ({n + 1}): {type(e).__name__}: {e}")
    for b in bases[:2]:
        try:
            v = GammaValue.of(b + 2)
            if not (isinstance(v, GammaValue) and v.base == b and v.mult == (b + 1) * b):
                out.append(f"канон of({b + 2}): ожидались слоты mult = {(b + 1) * b}, base = {b}, получено {v!r}")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"канон of({b + 2}): {type(e).__name__}: {e}")
    return out


# ---------- сертификат Γ: Стирлинг с остатком ----------

def _ln_reduced(y: Q, terms: int) -> Interval:
    """ln y для y > 0 с приведением: y = 2^k·r, r ∈ [1, 2) — ряд atanh
    сходится как 9^−N (t ≤ 1/3) при любом порядке y."""
    y = Q(y)
    k = 0
    while y >= 2:
        y, k = y / 2, k + 1
    while y < 1:
        y, k = y * 2, k - 1
    iv = ln_enclosure(y, terms) if y != 1 else Interval.point(0)
    if k:
        iv = iv + ln_enclosure(2, terms) * k
    return iv


def _ln_interval(iv: Interval, terms: int) -> Interval:
    """ln на отрезке положительных концов — монотонность."""
    return Interval(_ln_reduced(iv.lo, terms).lo, _ln_reduced(iv.hi, terms).hi)


def ln_gamma_enclosure(b, terms: int = 30) -> Interval:
    """ln Γ(b), b ∈ ℚ, b > 0 — отрезок ⊂ ℚ.

    z = b + n, n = terms: Γ(b) = Γ(z)/∏_{i<n}(b+i). Для z > 0:
    ln Γ(z) = (z−½)·ln z − z + ½·ln 2π + Σ_{k=1}^{K} B_{2k}/(2k(2k−1)z^{2k−1}) + R_K,
    |R_K| ≤ |B_{2K+2}|/((2K+2)(2K+1)z^{2K+1}) [внешний источник: DLMF
    5.11(ii)]. Все члены — точные пары; ln z, ln 2, ln π — сертификаты."""
    from .zeta import bernoulli
    b = Q(b)
    if b <= 0:
        raise Unresolved("ln Γ", f"Φ(C=1, value={b})",
                         hint="сертификат Γ — только при положительном аргументе")
    n = max(0, terms - int(b.numerator // b.denominator))
    z = b + n
    prod = Q(1)
    for i in range(n):
        prod *= b + i
    K = max(2, terms // 2)
    series = Q(0)
    zp = z                                     # z^(2k−1)
    for k in range(1, K + 1):
        series += bernoulli(2 * k) / (2 * k * (2 * k - 1) * zp)
        zp *= z * z
    bound = abs(bernoulli(2 * K + 2)) / ((2 * K + 2) * (2 * K + 1) * zp)
    ln_z = _ln_reduced(z, terms)
    ln_2pi = ln_enclosure(2, terms) + _ln_interval(pi_enclosure(terms), terms)
    out = (ln_z * (z - Q(1, 2)) - z + ln_2pi * Q(1, 2)
           + Interval(series - bound, series + bound))
    if n:
        out = out - _ln_reduced(prod, terms)
    return out


def gamma_enclosure(b, terms: int = 30) -> Interval:
    """Γ(b), b ∈ ℚ, b > 0: exp сертификата ln Γ(b) — монотонность exp."""
    lg = ln_gamma_enclosure(b, terms)
    mag = max(abs(lg.lo), abs(lg.hi))
    t = max(terms, 3 * (int(mag.numerator // mag.denominator) + 1) + 10)
    return Interval(exp_enclosure(lg.lo, t).lo, exp_enclosure(lg.hi, t).hi)


def gamma_half_integer(z):
    """Γ(z) при ПОЛУЦЕЛОМ z ТОЧНО, как моном класса 10: q·π^{1/2}.

    Это самодуальная точка отражения Эйлера. Отражение Γ(b)·Γ(1−b) =
    π/sin(πb) сводит ПАРУ; в точке b = 1/2, где b = 1−b, пара вырождается
    в квадрат Γ(1/2)² = π, и ОДНО значение прежде оставалось несведённым
    (`GammaValue.of(1/2)` → Φ(C=14 GAMMA, mult=1, base=1/2)) — при том что
    величина π^{1/2} в доме есть и живёт мономом класса 10.

    Равенство Γ(1/2) = π^{1/2} здесь НЕ ввозится извне, а выводится двумя
    дорогами самого дома:
      * Γ(1/2)·Γ(1/2) = π — отражением Эйлера (`GammaValue.__mul__`);
      * Γ(1/2) > 0 — сертификатом `gamma_enclosure(1/2)` (нижний конец
        положителен), потому из двух корней π берётся положительный.
    Обе дороги проверяются приёмкой, а не объявляются.

    Дальше по лестнице работает та же рекуррента Γ(z) = (z−1)·Γ(z−1),
    которой живёт `of()`: Γ(m+1/2) = mult·π^{1/2} с точным рациональным
    mult в обе стороны (Γ(3/2) = π^{1/2}/2, Γ(−1/2) = −2·π^{1/2}).

    Граница названа: знаменатель z вне {2} сюда не идёт — у Γ(p/q) при
    q ∉ {1, 2} замкнутой формы дом не знает, и её не знает математика
    (Чоула—Сельберг даёт лишь произведения при мнимом квадратичном поле).
    Целое z сводит `of()` факториалом — здесь оно тоже отказ по адресу.
    """
    if isinstance(z, float):
        raise TypeError("z — только int/Q, float запрещён")
    if not isinstance(z, _NUM):
        raise TypeError("z — только int/Q")
    z = Q(z)
    if z.denominator != 2:
        raise Unresolved("Γ", _rec(z),
                         hint="точная форма q·π^{1/2} есть только у "
                              "полуцелого z; целое сводит of() факториалом, "
                              "а Γ(p/q) при q ∉ {1,2} замкнутой формы не "
                              "имеет — её не знает и математика")
    g = cls_of_half(z)
    return PiMonomial._make(g.mult, Q(1, 2))


def cls_of_half(z) -> GammaValue:
    """Каноничный носитель Γ(z) при полуцелом z: mult·Γ(1/2), база 1/2.

    Отдельным именем — чтобы дорога «через mult рекурренты» была видна
    и проверяема отдельно от подстановки Γ(1/2) = π^{1/2}.
    """
    g = GammaValue.of(Q(z))
    if not isinstance(g, GammaValue) or g.base != Q(1, 2):
        raise Unresolved("Γ", _rec(z),
                         hint="полуцелое z обязано канонизоваться в базу 1/2")
    return g
