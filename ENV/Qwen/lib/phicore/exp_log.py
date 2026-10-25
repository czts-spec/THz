"""Носители классов 5 EXPONENTIAL и 6 LOGARITHMIC: coeff·e^y и coeff·ln(arg).

Обе линейки независимы от ℚ: e^y при y ∈ ℚ, y ≠ 0 и ln(arg) при
arg ∈ ℚ₊, arg ≠ 1 трансцендентны (теорема Линдемана—Вейерштрасса
[внешний источник]). Сведено только то, что сводится тождествами самой
линейки:

  * e^a · e^b = e^{a+b} — произведение экспонент складывает показатели;
  * сложение экспонент — лишь при одноимённом y (сумма coeff);
  * c·ln a + c·ln b = c·ln(ab), c·ln a − c·ln b = c·ln(a/b) — при равных
    coeff; при равных arg складываются coeff;
  * канон логарифма: arg > 1, точная целая степень свёрнута в coeff
    (ln 8 = 3·ln 2) — распознаётся целым корнем Ньютона, не таблицей.

Несведённое действие получает отказ Unresolved с явными Ф-записями
обоих операндов, не ложный класс и не молчание. Десятичная сторона —
только сертификатом сетки (enclosure), по явному запросу. Ни одного
float: float на входе — TypeError построением.
"""
from __future__ import annotations

from .rational import Q
from math import gcd

from .enclosure import Interval, exp_enclosure, ln_enclosure
from .errors import Unresolved
from .zero_state import Zero, ZeroKind

#: Допустимые скаляры. float исключён построением.
_NUM = (int, Q)

#: `[вычислено здесь]` 10.10, очередь строка 105, шаг второй: дверь оператора
#: носителя = дверь ФIOS, как у монома π (pi_value.py, шаг первый). Для ЧУЖОГО
#: класса операторы e^y и ln x не отказывают сами и не отвечают NotImplemented
#: — зовут `dispatch.combine` и возвращают её ответ (запись класса 19 там, где
#: точного правила нет; отказ пустой ячейки — насквозь). Свой класс и класс 1
#: (int, Q) решаются точными правилами и собственными отказами носителя, как
#: прежде: ячейки (E5, E1, ±), (E6, E1, ±) таблицы — refuse с той же сутью, а
#: (E5, E5, ±), (E6, E6, ±) — post_check, который звал бы оператор обратно.
#: Ввоз ленивый, внутри функций: Фios ввозит exp_log сам.
_IN_DOOR: set = set()

#: ячейки (класс другого операнда, действие), в которых дверь НЕ зовётся —
#: приёмки прибили в них отказ оператора поимённо (правило 09.10: краснеет
#: существующая приёмка — снять правку ровно в ней, приёмку не трогать).
#: Места `[вычислено здесь]` 10.10 — proto/test_phicore_exp_log.py:
#:   (6, +) — «e¹ + ln 2 — отказ (независимые линейки)»;
#:   (6, ×) — «e¹ · ln 2 — отказ, несёт обе записи»;
#:   (5, ×) — «ln 2 · e¹ — отказ и с другой стороны».
#: Дверь в тех же ячейках даёт запись класса 19 (e − ln 2 уже идёт дверью);
#: расхождение названо, снятие — вместе с пересмотром приёмки.
_PINNED: frozenset = frozenset({(6, "+"), (6, "*"), (5, "*")})

#: метка «чужой класс» — отлична от NotImplemented (законный ответ для нуля)
_FOREIGN = object()


def _pinned(other, op: str) -> bool:
    """Ячейка (класс другого операнда, действие) прибита приёмкой — _PINNED."""
    from .Фios import classify, EXACT      # ленивый ввоз: Фios ввозит exp_log
    k = classify(other)
    return k is not None and k[0] == EXACT and (k[1], op) in _PINNED


def _door(a, b, op: str):
    """Дверь ФIOS для чужого класса; второй вход с той же парой — NotImplemented
    (кольцо delegate → оператор → дверь разорвано)."""
    from .dispatch import combine          # ленивый ввоз
    key = (id(a), id(b), op)
    if key in _IN_DOOR:
        return NotImplemented
    _IN_DOOR.add(key)
    try:
        return combine(a, b, op)
    finally:
        _IN_DOOR.discard(key)


def _foreign(a, b, op: str, pin_other, hint: str | None = None):
    """Чужой класс: прибитая ячейка — прежний отказ носителя (hint) либо
    NotImplemented, когда своего текста у ячейки нет; иначе дверь."""
    if _pinned(pin_other, op):
        if hint is None:
            return NotImplemented
        raise Unresolved(op, _rec(a), _rec(b), hint=hint)
    return _door(a, b, op)


def _is_zero_num(x) -> bool:
    """Нуль класса 1 — int 0 или Q(0), не bool: x ± 0 = x в любом классе.
    Это поглощение нейтрального элемента, не сведение разных классов;
    до 01.10 дверь отказывала «независим: правила нет» и на нуле (наблюдение
    сдачи 25.09, очередь строка 104; замер 01.10 `[вычислено здесь]`:
    π + 0, 0 + π, π − 0, 0 − π — четыре отказа из четырёх)."""
    return isinstance(x, _NUM) and not isinstance(x, bool) and Q(x) == 0
#: Одна запись отказа на bool во всех дверях модуля: bool числом в доме не
#: считается и называется bool'ом (канон `phicore/rational.py`; сосед по
#: главе отвечает тем же: `Notch.from_q(True, 8)` → «bool не число»).
#: Двух копий одного канона здесь нет намеренно: разошедшиеся копии — уже
#: пройденная поломка (`affine.py:75`).
_NOT_BOOL = "bool не число"


def _rec(x: object) -> str:
    """Ф-запись операнда для сообщений отказа."""
    r = getattr(x, "record", None)
    if callable(r):
        return r()
    if isinstance(x, _NUM):
        return f"Φ(C=1 RATIONAL, value={Q(x)})"
    return f"Φ(?={x!r})"


def _as_fraction(x, name: str) -> Q:
    """int/Q → Q; float и прочее — TypeError построением."""
    if isinstance(x, bool):
        # bool числом в доме не считается, и назван он обязан быть bool'ом:
        # канон `phicore/rational.py` («bool не число»), тем же отвечают
        # соседи по главе (`Notch.from_q(True, 8)`). Замер 12.09: одна
        # запись на два разных входа давала float'у и bool'у ОДИН диагноз,
        # а ветвление у вызывающего идёт по тексту отказа.
        raise TypeError(f"{name} — {_NOT_BOOL}; дверь принимает int/Q")
    if not isinstance(x, _NUM):
        raise TypeError(f"{name} — только int/Q, float запрещён")
    return Q(x)


# ---------- целые корни: свёртка совершенных степеней без таблиц ----------

def _iroot(n: int, k: int) -> int:
    """⌊n^{1/k}⌋ целым счётом (Ньютон на int), n ≥ 0, k ≥ 1."""
    if n < 2 or k == 1:
        return n
    # старт заведомо не ниже корня: n < 2^bits ⇒ n^{1/k} < 2^⌈bits/k⌉
    x = 1 << ((n.bit_length() + k - 1) // k)
    while True:
        y = ((k - 1) * x + n // x ** (k - 1)) // k
        if y >= x:
            return x
        x = y


def _power_index(n: int) -> int:
    """Наибольшее e ≥ 1 с n = a^e (n ≥ 2): перебор показателей до битовой
    длины, каждый кандидат проверяется точным возведением r^k == n."""
    best = 1
    for k in range(2, n.bit_length() + 1):
        r = _iroot(n, k)
        if r ** k == n:
            best = k
    return best


# ---------- класс 5: EXPONENTIAL ----------

class ExpValue:
    """coeff·e^y: coeff ≠ 0 (Q), y ≠ 0 (Q)."""

    __slots__ = ("coeff", "y")

    def __init__(self, coeff, y):
        coeff = _as_fraction(coeff, "coeff")
        y = _as_fraction(y, "y")
        if coeff == 0:
            raise ValueError("coeff=0 не есть экспонента — ноль живёт в "
                             "zero_state; стройте через make")
        if y == 0:
            raise ValueError("y=0 — обычное число, не экспонента; "
                             "стройте через make")
        self.coeff, self.y = coeff, y

    # --- сборка результата: вырождение в Q ---
    @classmethod
    def make(cls, coeff, y):
        """Фабрика с вырождением: y=0 → сам коэффициент (e⁰=1),
        coeff=0 → точный ноль."""
        coeff = _as_fraction(coeff, "coeff")
        y = _as_fraction(y, "y")
        if coeff == 0:
            return Q(0)
        if y == 0:
            return Q(coeff)
        return cls(coeff, y)

    # --- умножение / деление: точные правила впереди, чужой класс — дверь ---
    def _mul_exact(self, other):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__rmul__
        if isinstance(other, ExpValue):
            # e^a·e^b = e^{a+b}; одноимённая линейка сводится показателями
            return self.make(self.coeff * other.coeff, self.y + other.y)
        if isinstance(other, _NUM):
            q = Q(other)
            if q == 0:
                return Q(0)
            return self.make(self.coeff * q, self.y)
        return _FOREIGN

    _MUL_PIN_HINT = "произведение экспоненты и логарифма не сведено"

    def __mul__(self, other):
        r = self._mul_exact(other)
        if r is not _FOREIGN:
            return r
        return _foreign(self, other, "*", other, self._MUL_PIN_HINT)

    def __rmul__(self, other):
        r = self._mul_exact(other)
        if r is not _FOREIGN:
            return r
        return _foreign(other, self, "*", other, self._MUL_PIN_HINT)  # порядок — как в выражении

    def inverse(self) -> "ExpValue | Q":
        """1/(c·e^y) = (1/c)·e^{−y}."""
        return self.make(1 / self.coeff, -self.y)

    def __truediv__(self, other):
        if isinstance(other, ExpValue):
            return self * other.inverse()
        if isinstance(other, _NUM):
            q = Q(other)
            if q == 0:
                raise Unresolved("/", self.record(), _rec(other),
                                 hint="x/0 здесь не сведено")
            return self.make(self.coeff / q, self.y)
        if isinstance(other, Zero):
            return _door(self, other, "/")    # таблица: «x/0 здесь не сведено»
        return _foreign(self, other, "/", other)

    def __rtruediv__(self, other):
        """число / (c·e^y) = (число/c)·e^{−y}; 0/e^y = точный ноль."""
        if isinstance(other, _NUM):
            q = Q(other)
            if q == 0:
                return Q(0)
            return self.make(q / self.coeff, -self.y)
        return _foreign(other, self, "/", other)

    def __pow__(self, n: int):
        if not isinstance(n, int) or isinstance(n, bool):
            return NotImplemented
        if n == 0:
            return Q(1)
        return self.make(self.coeff ** n, self.y * n)

    def __neg__(self) -> "ExpValue":
        return ExpValue(-self.coeff, self.y)

    # --- сложение: только одноимённый показатель ---
    def _addsub(self, other, op: str, sign: int):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__radd__/__rsub__
        if _is_zero_num(other):
            return self             # x ± 0 = x: нуль класса 1 поглощается
        if isinstance(other, ExpValue) and other.y == self.y:
            s = self.coeff + sign * other.coeff
            if s == 0:
                return Q(0)
            return self.make(s, self.y)
        # класс 1 и одноимённый класс с другим показателем — отказ самого
        # носителя (ячейки таблицы refuse / post_check, дверь звала бы обратно)
        if isinstance(other, (_NUM, ExpValue)) or _pinned(other, op):
            raise Unresolved(op, self.record(), _rec(other),
                             hint="e^y независим: сложение сведено только "
                                  "для одноимённого показателя y")
        if callable(getattr(other, "record", None)):
            return _door(self, other, op)   # чужой класс — дверь ФIOS
        return NotImplemented

    def __add__(self, other):
        return self._addsub(other, "+", +1)

    def __sub__(self, other):
        return self._addsub(other, "-", -1)

    def __radd__(self, other):
        if _is_zero_num(other):
            return self             # 0 + x = x
        if isinstance(other, _NUM) or _pinned(other, "+"):
            raise Unresolved("+", _rec(other), self.record(),
                             hint="e^y независим: правила нет")
        if callable(getattr(other, "record", None)):
            return _door(other, self, "+")
        return NotImplemented

    def __rsub__(self, other):
        if _is_zero_num(other):
            return self.make(-self.coeff, self.y)   # 0 − x = −x
        if isinstance(other, _NUM) or _pinned(other, "-"):
            raise Unresolved("-", _rec(other), self.record(),
                             hint="e^y независим: правила нет")
        if callable(getattr(other, "record", None)):
            return _door(other, self, "-")
        return NotImplemented

    # --- запись, сравнение, сертификат ---
    def self_check(self, deep: bool = False) -> list:
        """Самопроверка экспоненты ВСТРЕЧНОЙ дорогой (ось самоприменения).

        `[вычислено здесь]` Дороги: (1) тождество САМОЙ линейки, не ряд:
        e^y · e^(−y) обязано свестись классом 5 в точную единицу (сложение
        показателей), то есть произведение проверяет запись, а не сертификат
        сам себя; (2) `[внешний источник]` неравенства Эйлера 1 + y ≤ e^y
        (всюду) и e^y ≤ 1/(1 − y) (при y < 1) — ряд Тейлора в
        `exp_enclosure` о них не знает, и сертификат обязан им отвечать;
        (3) запись называет класс 5 и оба поля, контракт eq/hash держится.

        Граница объявлена: верхнее неравенство годится ТОЛЬКО при y < 1 —
        при y ≥ 1 остаётся нижнее, и молчание верхнего находкой не
        считается; сертификат берётся при 30 членах, у большого |y| он
        широк, и ширина сама по себе тоже не находка (уточнение — явным
        запросом). Возвращает список находок; исключений не бросает.
        """
        out = []
        if not isinstance(self.coeff, Q) or not isinstance(self.y, Q):
            out.append(f"поля не пары: coeff={self.coeff!r}, y={self.y!r}")
            return out
        if self.coeff == 0:
            out.append("coeff = 0: ноль записан экспонентой — ноль живёт "
                       "состоянием (zero_state)")
        if self.y == 0:
            out.append("y = 0: обычное число записано экспонентой")
        if "C=5 EXPONENTIAL" not in self.record():
            out.append(f"запись не называет класс 5: {self.record()}")
        if hash(self) != hash(ExpValue(self.coeff, self.y)):
            out.append("контракт eq/hash нарушен: равные записи — разные ключи")
        if out:
            return out
        back = ExpValue(Q(1), self.y) * ExpValue(Q(1), -self.y)
        if not (isinstance(back, Q) and back == 1):
            out.append(f"e^y · e^(−y) не свелось в единицу: {back!r}")
        try:
            iv = self.enclosure(30)
        except Exception as e:                                  # noqa: BLE001
            out.append(f"сертификат отказал: {type(e).__name__}: {e}")
            return out
        lo_e, hi_e = Q(1) + self.y, None                       # [внешний источник]
        if self.y < 1:
            hi_e = Q(1) / (Q(1) - self.y)
        if self.coeff > 0:
            b_lo = self.coeff * lo_e
            b_hi = None if hi_e is None else self.coeff * hi_e
            if iv.hi < b_lo:
                out.append(f"сертификат ниже неравенства Эйлера 1 + y ≤ e^y: "
                           f"{iv.record()} против {b_lo}")
            if b_hi is not None and iv.lo > b_hi:
                out.append(f"сертификат выше неравенства e^y ≤ 1/(1 − y): "
                           f"{iv.record()} против {b_hi}")
        else:
            b_hi = self.coeff * lo_e
            b_lo = None if hi_e is None else self.coeff * hi_e
            if iv.lo > b_hi:
                out.append(f"сертификат выше неравенства (коэффициент "
                           f"отрицателен): {iv.record()} против {b_hi}")
            if b_lo is not None and iv.hi < b_lo:
                out.append(f"сертификат ниже неравенства (коэффициент "
                           f"отрицателен): {iv.record()} против {b_lo}")
        if not deep:
            return out
        fine = self.enclosure(44)
        if not iv.encloses(fine):
            out.append(f"уточнение не вложено в грубый сертификат: "
                       f"{fine.record()} ⊄ {iv.record()}")
        sq = self * self
        if isinstance(sq, ExpValue) and sq.y != 2 * self.y:
            out.append(f"квадрат не удвоил показатель: {sq.record()}")
        return out

    def record(self) -> str:
        return f"Φ(C=5 EXPONENTIAL, coeff={self.coeff}, y={self.y})"

    __repr__ = record

    def __eq__(self, other):
        if isinstance(other, ExpValue):
            return self.coeff == other.coeff and self.y == other.y
        if isinstance(other, _NUM):
            # e^y при y ≠ 0 трансцендентно [внешний источник: Линдеман]
            return False
        return NotImplemented

    def __hash__(self):
        return hash((ExpValue, self.coeff, self.y))

    def enclosure(self, terms: int = 30) -> Interval:
        """Сертификат сетки [lo, hi] ⊂ ℚ: coeff · вложение e^y."""
        return exp_enclosure(self.y, terms) * self.coeff


# ---------- класс 6: LOGARITHMIC ----------

class LogValue:
    """coeff·ln(arg): coeff ≠ 0 (Q), arg > 0, arg ≠ 1 (Q).

    Канон (arg > 1, совершенная степень свёрнута) держит фабрика make;
    в каноне совпадение величин распознаётся сравнением слотов."""

    __slots__ = ("coeff", "arg")

    def __init__(self, coeff, arg):
        coeff = _as_fraction(coeff, "coeff")
        arg = _as_fraction(arg, "arg")
        if arg <= 0:
            raise Unresolved("ln", f"Φ(C=1 RATIONAL, value={arg})",
                             hint="ln на неположительном не веществен")
        if coeff == 0:
            raise ValueError("coeff=0 не есть логарифм — ноль живёт в "
                             "zero_state; стройте через make")
        if arg == 1:
            raise ValueError("arg=1 — ln 1 есть обычный ноль; "
                             "стройте через make")
        self.coeff, self.arg = coeff, arg

    # --- сборка результата: вырождение и нормализация к канону ---
    @classmethod
    def make(cls, coeff, arg):
        """Фабрика с вырождением и каноном arg > 1.

        arg=1 или coeff=0 → точный ноль; arg < 1 → знак уходит в coeff
        (c·ln a = −c·ln(1/a)); точная целая степень arg = arg0^k
        (k ≥ 2, по числителю и знаменателю) свёрнута:
        c·ln(arg0^k) = (c·k)·ln(arg0)."""
        coeff = _as_fraction(coeff, "coeff")
        arg = _as_fraction(arg, "arg")
        if arg <= 0:
            raise Unresolved("ln", f"Φ(C=1 RATIONAL, value={arg})",
                             hint="ln на неположительном не веществен")
        if arg == 1 or coeff == 0:
            return Q(0)
        if arg < 1:
            return cls.make(-coeff, 1 / arg)
        n, d = arg.numerator, arg.denominator
        k = _power_index(n) if d == 1 else gcd(_power_index(n),
                                               _power_index(d))
        if k >= 2:
            return cls.make(coeff * k, Q(_iroot(n, k), _iroot(d, k)))
        return cls(coeff, arg)

    # --- сложение / вычитание: равные arg либо равные coeff ---
    def _addsub(self, other, op: str, sign: int):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__radd__/__rsub__
        if _is_zero_num(other):
            return self             # x ± 0 = x: нуль класса 1 поглощается
        if isinstance(other, LogValue):
            if other.arg == self.arg:
                s = self.coeff + sign * other.coeff
                if s == 0:
                    return Q(0)
                return self.make(s, self.arg)
            if other.coeff == self.coeff:
                # c·ln a ± c·ln b = c·ln(a·b^{±1}); канон соберёт make
                new_arg = (self.arg * other.arg if sign > 0
                           else self.arg / other.arg)
                return self.make(self.coeff, new_arg)
            if (Q(self.coeff).denominator == 1
                    and Q(other.coeff).denominator == 1):
                # ОБЩИЙ закон логарифма при ЦЕЛЫХ coeff (12.09, §3 строка 3):
                #     c₁·ln a ± c₂·ln b = ln(a^{c₁}·b^{±c₂}),
                # и аргумент остаётся РАЦИОНАЛЬНЫМ — значит величина не
                # выходит из класса 6 и правило законно. Прежде здесь стоял
                # отказ, и двойственность 7↔8 (§3, строка 3) ломалась на
                # накоплении: Σ ln fₙ произведения Валлиса доходила до N = 4
                # и упиралась на пятом члене, потому что канон свернул
                # совершенный квадрат в coeff = 2, а у нового члена coeff = 1.
                # Замер 12.09: 4 шага из 8 закрывались, 4 отказывали.
                # При ДРОБНОМ coeff правило НЕ действует, и отказ ниже верен
                # по ПРЕДМЕТУ, а не по бедности: ln 3 + ½·ln 2 = ln(3·√2)
                # требует класса 3, которого линейка 6 не несёт.
                a = self.arg ** Q(self.coeff).numerator
                b = other.arg ** (sign * Q(other.coeff).numerator)
                return self.make(Q(1), a * b)
        if isinstance(other, (_NUM, LogValue)) or _pinned(other, op):
            raise Unresolved(op, self.record(), _rec(other),
                             hint="ln сведён при равных arg (сумма coeff) "
                                  "или равных coeff (произведение "
                                  "аргументов); иного правила нет")
        if callable(getattr(other, "record", None)):
            return _door(self, other, op)   # чужой класс — дверь ФIOS
        return NotImplemented

    def __add__(self, other):
        return self._addsub(other, "+", +1)

    def __sub__(self, other):
        return self._addsub(other, "-", -1)

    def __radd__(self, other):
        if _is_zero_num(other):
            return self             # 0 + x = x
        if isinstance(other, _NUM) or _pinned(other, "+"):
            raise Unresolved("+", _rec(other), self.record(),
                             hint="ln независим: правила нет")
        if callable(getattr(other, "record", None)):
            return _door(other, self, "+")
        return NotImplemented

    def __rsub__(self, other):
        if _is_zero_num(other):
            return self.make(-self.coeff, self.arg)  # 0 − x = −x
        if isinstance(other, _NUM) or _pinned(other, "-"):
            raise Unresolved("-", _rec(other), self.record(),
                             hint="ln независим: правила нет")
        if callable(getattr(other, "record", None)):
            return _door(other, self, "-")
        return NotImplemented

    # --- умножение / деление: точные правила впереди, чужой класс — дверь ---
    def _mul_exact(self, other):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__rmul__
        if isinstance(other, _NUM):
            q = Q(other)
            if q == 0:
                return Q(0)
            return self.make(self.coeff * q, self.arg)
        return _FOREIGN

    _MUL_PIN_HINT = ("произведение логарифмов и произведение с экспонентой "
                     "не сведены")

    def __mul__(self, other):
        r = self._mul_exact(other)
        if r is not _FOREIGN:
            return r
        # произведение логарифмов — тоже чужая для линейки величина: дверь
        return _foreign(self, other, "*", other, self._MUL_PIN_HINT)

    def __rmul__(self, other):
        r = self._mul_exact(other)
        if r is not _FOREIGN:
            return r
        return _foreign(other, self, "*", other, self._MUL_PIN_HINT)

    def __truediv__(self, other):
        if isinstance(other, _NUM):
            q = Q(other)
            if q == 0:
                raise Unresolved("/", self.record(), _rec(other),
                                 hint="x/0 здесь не сведено")
            return self.make(self.coeff / q, self.arg)
        if isinstance(other, LogValue):
            if other.arg == self.arg:
                # c₁·ln g / (c₂·ln g) = c₁/c₂ — канон делает одно
                # основание распознаваемым (ln 8 / ln 2 = 3)
                return self.coeff / other.coeff
            # разные основания — величина вне линейки: дверь (запись класса 19)
            return _foreign(self, other, "/", other)
        if isinstance(other, Zero):
            return _door(self, other, "/")    # таблица: «x/0 здесь не сведено»
        return _foreign(self, other, "/", other)

    def __rtruediv__(self, other):
        if isinstance(other, _NUM):
            q = Q(other)
            if q == 0:
                return Q(0)
            # число / ln a — обращение вне линейки: дверь (запись класса 19)
            return _door(other, self, "/")
        return _foreign(other, self, "/", other)

    def __neg__(self) -> "LogValue":
        return LogValue(-self.coeff, self.arg)

    def __pow__(self, n: int):
        if not isinstance(n, int) or isinstance(n, bool):
            return NotImplemented
        if n == 0:
            return Q(1)
        if n == 1:
            return self
        raise Unresolved("^", self.record(), _rec(n),
                         hint="степень логарифма не сведена")

    # --- запись, сравнение, сертификат ---
    def self_check(self, deep: bool = False) -> list:
        """Самопроверка логарифма ВСТРЕЧНОЙ дорогой (ось самоприменения).

        `[вычислено здесь]` Дороги: (1) `[внешний источник]` неравенства
        Наперовой границы 1 − 1/a ≤ ln a ≤ a − 1 (всюду при a > 0, равенство
        только при a = 1) — ряд внутри `enclosure` о них не знает, и
        сертификат обязан лежать между ними, с переменой концов при
        отрицательном коэффициенте; (2) знак величины против положения
        аргумента: ln a > 0 ⟺ a > 1 — сверяется по сертификату, а не по
        слоту; (3) запись называет класс 6 и оба поля, контракт eq/hash.

        Граница объявлена: канон фабрики (`make`: arg > 1, совершенная
        степень свёрнута) здесь НЕ требуется — прямой конструктор законно
        строит запись с arg < 1, и это не находка, а другая запись той же
        величины; сведение записей — работа `make`, не самопроверки.
        Возвращает список находок; исключений не бросает.
        """
        out = []
        if not isinstance(self.coeff, Q) or not isinstance(self.arg, Q):
            out.append(f"поля не пары: coeff={self.coeff!r}, arg={self.arg!r}")
            return out
        if self.arg <= 0:
            out.append(f"arg = {self.arg} ≤ 0: ln на неположительном не веществен")
        if self.arg == 1:
            out.append("arg = 1: ln 1 есть обычный ноль, не логарифм")
        if self.coeff == 0:
            out.append("coeff = 0: ноль записан логарифмом — ноль живёт "
                       "состоянием (zero_state)")
        if "C=6 LOGARITHMIC" not in self.record():
            out.append(f"запись не называет класс 6: {self.record()}")
        if hash(self) != hash(LogValue(self.coeff, self.arg)):
            out.append("контракт eq/hash нарушен: равные записи — разные ключи")
        if out:
            return out
        try:
            iv = self.enclosure(60)
        except Exception as e:                                  # noqa: BLE001
            out.append(f"сертификат отказал: {type(e).__name__}: {e}")
            return out
        lo_l = Q(1) - Q(1) / self.arg                          # [внешний источник]
        hi_l = self.arg - Q(1)
        b_lo, b_hi = (self.coeff * lo_l, self.coeff * hi_l) if self.coeff > 0 \
            else (self.coeff * hi_l, self.coeff * lo_l)
        if iv.hi < b_lo or iv.lo > b_hi:
            out.append(f"сертификат вышел за границы 1 − 1/a ≤ ln a ≤ a − 1: "
                       f"{iv.record()} против ⟦{b_lo}; {b_hi}⟧")
        positive = (self.arg > 1) == (self.coeff > 0)
        if positive and iv.hi < 0:
            out.append(f"знак величины против положения аргумента: "
                       f"{self.record()} даёт {iv.record()}")
        if not positive and iv.lo > 0:
            out.append(f"знак величины против положения аргумента: "
                       f"{self.record()} даёт {iv.record()}")
        if not deep:
            return out
        fine = self.enclosure(90)
        if not iv.encloses(fine):
            out.append(f"уточнение не вложено в грубый сертификат: "
                       f"{fine.record()} ⊄ {iv.record()}")
        return out

    def record(self) -> str:
        return f"Φ(C=6 LOGARITHMIC, coeff={self.coeff}, arg={self.arg})"

    __repr__ = record

    def __eq__(self, other):
        if isinstance(other, LogValue):
            return self.coeff == other.coeff and self.arg == other.arg
        if isinstance(other, _NUM):
            # ln(arg) при arg ∈ ℚ₊, arg ≠ 1 трансцендентен
            # [внешний источник: Линдеман—Вейерштрасс]
            return False
        return NotImplemented

    def __hash__(self):
        return hash((LogValue, self.coeff, self.arg))

    def enclosure(self, terms: int = 60) -> Interval:
        """Сертификат сетки [lo, hi] ⊂ ℚ: coeff · вложение ln(arg)."""
        return ln_enclosure(self.arg, terms) * self.coeff


# ---------- композиция классов 5 ↔ 6: exp и ln как функции ----------

def exp(x):
    """e^x по классу аргумента (композиционные правила 5↔6, известное
    аппаратом): пара y → e^y (класс 5; y=0 → 1); c·ln a → a^c ТОЧНО
    (пара либо линейка 3 — запись Φα над рациональной базой a); ноль-
    состояние абсолютный → 1. Прочее (e^{e^y}, e^π, e^θ) носителя не
    имеет — отказ поимённо, не приближение."""
    if isinstance(x, Zero):
        if x.kind is ZeroKind.ABSOLUTE:
            return Q(1)
        raise Unresolved("exp", x.record(),
                         hint="e^0 для ноля-состояния с провенансом не сведено")
    if isinstance(x, _NUM):
        return ExpValue.make(1, Q(x))
    if isinstance(x, LogValue):
        from .alpha import PhiAlpha
        return PhiAlpha(1, x.arg, x.coeff).as_carrier()     # e^{c·ln a} = a^c
    raise Unresolved("exp", _rec(x),
                     hint="e^x для этого носителя линейкой не несётся "
                          "(e^{e^y}, e^π, e^θ — носителя нет)")


def ln(x):
    """ln x по классу аргумента: пара a>0 → ln a (класс 6; a=1 → 0);
    c·e^y при c=1 → y ТОЧНО (пара); запись Φα над рациональной базой
    m·a^E при m>0 → ln(m·a^E) (целый E — одним логарифмом; дробный —
    суммой ln m + E·ln a по правилам линейки 6); запись над e с
    рациональным показателем и мантиссой 1 → E. Прочее — отказ поимённо:
    y + ln c при c≠1, ln ln a, ln π, ln θ носителя не имеют."""
    from .alpha import PhiAlpha
    if isinstance(x, _NUM):
        return LogValue.make(1, Q(x))
    if isinstance(x, ExpValue):
        if x.coeff == 1:
            return x.y
        raise Unresolved("ln", x.record(),
                         hint=f"ln(c·e^y) = y + ln c при c={x.coeff} ≠ 1: "
                              "сумма пары и логарифма линейкой не несётся")
    if isinstance(x, PhiAlpha):
        a, e, m = x.alpha, x.e, x.m
        if isinstance(a, Q) and isinstance(e, Q) and m > 0:
            if e.d == 1:
                return LogValue.make(1, m * a ** e)
            if m == 1:
                return LogValue.make(e, a)                  # E·ln a
            return LogValue.make(1, m) + LogValue.make(e, a)  # правила линейки 6
        if isinstance(a, ExpValue) and a.coeff == 1 and a.y == 1 \
                and isinstance(e, Q) and m == 1:
            return e
    raise Unresolved("ln", _rec(x),
                     hint="ln для этого носителя линейкой не несётся "
                          "(ln ln a, ln π, ln θ, y + ln c — носителя нет)")
