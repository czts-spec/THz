"""Линейка класса 12 HYPERBOLIC: sinh/cosh/tanh рационального аргумента.

Носитель — СТРУКТУРНАЯ ПАРА над классом 5: величина определена парой
экспонент sinh x = (e^x − e^{−x})/2, cosh x = (e^x + e^{−x})/2 — и
хранится структурно слотами (fn, x, sign): имя функции fn ∈ {sinh, cosh,
tanh}, аргумент x ∈ ℚ в каноне x > 0, знак sign ∈ {+1, −1} — множитель
величины. Правила линейки выражаются через ExpValue, не через таблицы.

Направленные правила фабрики hyperbolic():
  * x = 0 — вырождение в точный класс 1: sinh 0 = 0, cosh 0 = 1,
    tanh 0 = 0;
  * чётность — каноном x > 0: sinh(−x) = −sinh(x) и tanh(−x) = −tanh(x)
    несут знак слотом sign (нечётность), cosh(−x) = cosh(x) знака не
    порождает (чётность).

ДВОЙСТВЕННОСТЬ с классом 5: пифагорово тождество cosh²x − sinh²x = 1
сводится функцией pythagorean_check(x) РАСКРЫТИЕМ через ExpValue —
квадраты раскрываются в экспоненты, одноимённые показатели гасятся
классом 5 сам (e^{2x} − e^{2x} = 0), а перекрёстный член e^x·e^{−x} = 1
сводит тоже класс 5 (сложение показателей, вырождение e⁰). Иных правил
арифметики нет: суммы и произведения гиперболик, сложение с рациональным,
масштаб, частное, степень вне {0, 1} — отказ Unresolved с явными
Ф-записями обоих операндов, не ложный класс.

Величина нетривиального слота не равна ни одному рациональному и слоту
другой функции линейки: sinh/cosh/tanh ненулевого рационального
трансцендентны, а совпадение разных функций свелось бы к исчезающей
целочисленной комбинации экспонент с различными рациональными
показателями [внешний источник: теорема Линдемана—Вейерштрасса].

Десятичная сторона — только сертификатом сетки (enclosure) по явному
запросу: sinh/cosh — интервальная полусумма/полуразность вложений
exp_enclosure(x) и exp_enclosure(−x); tanh — интервальное частное
sinh/cosh (cosh ≥ 1 — ноль в знаменатель не попадает). Ни одного float:
float на входе — TypeError построением.
"""
from __future__ import annotations

from .rational import Q

from .enclosure import Interval, exp_enclosure
from .errors import Unresolved
from .exp_log import ExpValue
from .zero_state import Zero

#: Допустимые скаляры. float исключён построением.
_NUM = (int, Q)

#: `[вычислено здесь]` 10.10, очередь строка 105, шаг седьмой: дверь оператора
#: носителя = дверь ФIOS (как у π, e^y/ln, arctrig, ζ, Γ, корня J_ν). Для ЧУЖОГО
#: класса с записью операторы зовут `dispatch.combine`; класс 1 (int, Q) и
#: одноимённая линейка — прежние правила и собственные отказы носителя.
#: Ввоз ленивый, внутри функций.
_IN_DOOR: set = set()
_PINNED: frozenset = frozenset()


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
    return (not isinstance(x, (_NUM, HyperbolicValue))) and callable(getattr(x, "record", None))
#: Одна запись отказа на bool во всех дверях модуля: bool числом в доме не
#: считается и называется bool'ом (канон `phicore/rational.py`; сосед по
#: главе отвечает тем же: `Notch.from_q(True, 8)` → «bool не число»).
#: Двух копий одного канона здесь нет намеренно: разошедшиеся копии — уже
#: пройденная поломка (`affine.py:75`).
_NOT_BOOL = "bool не число"

#: Имена функций линейки; иное имя — не слот.
_FNS = ("sinh", "cosh", "tanh")

#: Нечётные функции: знак аргумента уходит в слот sign.
_ODD = ("sinh", "tanh")


def _rec(x: object) -> str:
    """Ф-запись операнда для сообщений отказа."""
    r = getattr(x, "record", None)
    if callable(r):
        return r()
    if isinstance(x, _NUM):
        return f"Φ(C=1 RATIONAL, value={Q(x)})"
    return f"Φ(?={x!r})"


def _check_fn_x(fn, x) -> tuple:
    """Общая проверка входов фабрики и конструктора: имя из линейки,
    аргумент — только int/Q; float и bool — TypeError построением."""
    if not isinstance(fn, str):
        raise TypeError("fn — только str")
    if fn not in _FNS:
        raise ValueError("fn ∈ {sinh, cosh, tanh}; иного линейка не несёт")
    if isinstance(x, bool):
        raise TypeError(f"x — {_NOT_BOOL}; дверь принимает int/Q")
    if not isinstance(x, _NUM):
        raise TypeError("x — только int/Q, float запрещён")
    return fn, Q(x)


# ---------- двойственность с классом 5: пифагорово тождество ----------

def pythagorean_check(x) -> Q:
    """cosh(x)² − sinh(x)² = точная 1 раскрытием через ExpValue.

    Раскрытие структурной пары:
      cosh²x = e^{2x}/4 + e^{−2x}/4 + (1/2)·e^x·e^{−x},
      sinh²x = e^{2x}/4 + e^{−2x}/4 − (1/2)·e^x·e^{−x};
    в разности одноимённые показатели гасит сам класс 5
    (e^{2x} − e^{2x} = 0), а перекрёстный член e^x·e^{−x} класс 5
    сводит сложением показателей в e⁰ = 1 — двойственность линеек.
    Каждый шаг обязан свестись в Q; иначе — честный отказ,
    не подстановка ответа.
    """
    if isinstance(x, bool):
        raise TypeError(f"x — {_NOT_BOOL}; дверь принимает int/Q")
    if not isinstance(x, _NUM):
        raise TypeError("x — только int/Q, float запрещён")
    x = Q(x)
    if x == 0:
        # вырождение: cosh 0 = 1, sinh 0 = 0 — тождество рационально
        return Q(1)
    quarter = Q(1, 4)
    e_hi = ExpValue.make(quarter, 2 * x)     # e^{2x}/4
    e_lo = ExpValue.make(quarter, -2 * x)    # e^{−2x}/4
    # одноимённые показатели: разность сводит класс 5 в точный ноль
    d_hi = e_hi - e_hi
    d_lo = e_lo - e_lo
    # перекрёстный член: e^x·e^{−x} — сложение показателей, e⁰ = 1
    cross = ExpValue.make(1, x) * ExpValue.make(1, -x)
    if not (isinstance(d_hi, Q) and d_hi == 0
            and isinstance(d_lo, Q) and d_lo == 0
            and isinstance(cross, Q)):
        raise Unresolved("cosh²−sinh²", _rec(x),
                         hint="класс 5 не свёл раскрытие структурной "
                              "пары — тождество не приземлено")
    return d_hi + d_lo + cross


# ---------- сама линейка ----------

class HyperbolicValue:
    """sign·fn(x): fn ∈ {sinh, cosh, tanh}, x ∈ ℚ, x > 0, sign ∈ {±1}.

    Канон держится построением: x = 0 в слот не помещается (вырождение
    делает фабрика hyperbolic()), x < 0 не строится вовсе — чётность
    приводит аргумент к x > 0 со знаком у нечётных функций. Иного
    множителя линейка не несёт: слоты только (fn, x, sign).
    """

    __slots__ = ("fn", "x", "sign")

    def __init__(self, fn, x, sign=1):
        fn, x = _check_fn_x(fn, x)
        if x == 0:
            raise ValueError("x=0 — вырождение в класс 1; стройте через "
                             "фабрику hyperbolic()")
        if x < 0:
            raise ValueError("канон требует x > 0; чётность — фабрикой "
                             "hyperbolic()")
        if isinstance(sign, bool):
            raise TypeError(f"sign — {_NOT_BOOL}; дверь принимает int")
        if not isinstance(sign, int):
            raise TypeError("sign — только int ∈ {1, −1}")
        if sign not in (1, -1):
            raise ValueError("sign ∈ {1, −1}; иного множителя слот не несёт")
        self.fn, self.x, self.sign = fn, x, sign

    # --- запись ---------------------------------------------------------

    def self_check(self, deep: bool = False) -> list:
        """Самопроверка гиперболической величины ВСТРЕЧНОЙ дорогой.

        `[вычислено здесь]` Носитель хранится СТРУКТУРНО — тройкой
        (fn, x, sign), величины в нём нет; поэтому встречные дороги обе
        внешние этой тройке: (1) тождество Пифагора cosh²x − sinh²x = 1,
        приземлённое через класс 5 (`pythagorean_check`: раскрытие пары
        экспонент, где одноимённые показатели гасит сама линейка, а
        перекрёстный член сводится в e⁰) — тройка о нём не знает; (2)
        `[внешний источник]` неравенства анализа на сертификате: при x > 0
        sinh x ≥ x, cosh x ≥ 1 + x²/2, |tanh x| < 1 и cosh x > sinh x —
        ряд внутри `enclosure` их не использует; (3) запись называет класс
        12, знак несёт множитель, контракт eq/hash держится.

        Граница объявлена: канон хранит x > 0 и знак отдельным слотом,
        поэтому неравенства проверяются на модуле величины, а знак — на
        слоте; при x = 0 вырождение (sinh 0 = 0, tanh 0 = 0) в носитель не
        попадает вовсе — его отдаёт фабрика `hyperbolic` парой, и это
        объявлено там. Возвращает список находок; исключений не бросает.
        """
        out = []
        if self.fn not in ("sinh", "cosh", "tanh"):
            out.append(f"имя функции вне линейки: {self.fn!r}")
        if not isinstance(self.x, Q) or self.x <= 0:
            out.append(f"аргумент вне канона x > 0: {self.x!r}")
        if self.sign not in (1, -1):
            out.append(f"знак не ±1: {self.sign!r}")
        if out:
            return out
        if "C=12" not in self.record():
            out.append(f"запись не называет класс 12: {self.record()}")
        if hash(self) != hash(HyperbolicValue(self.fn, self.x, self.sign)):
            out.append("контракт eq/hash нарушен: равные записи — разные ключи")
        try:
            one = pythagorean_check(self.x)
            if one != 1:
                out.append(f"тождество cosh²−sinh² приземлилось не в единицу: "
                           f"{one}")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"тождество Пифагора не приземлилось: "
                       f"{type(e).__name__}: {e}")
        try:
            iv = self.enclosure(30)
        except Exception as e:                                  # noqa: BLE001
            out.append(f"сертификат отказал: {type(e).__name__}: {e}")
            return out
        lo, hi = (iv.lo, iv.hi) if self.sign > 0 else (-iv.hi, -iv.lo)
        if self.fn == "sinh" and hi < self.x:
            out.append(f"sinh x < x при x > 0: {iv.record()} против {self.x}")
        if self.fn == "cosh" and hi < Q(1) + self.x * self.x / Q(2):
            out.append(f"cosh x < 1 + x²/2: {iv.record()}")
        if self.fn == "tanh" and (lo <= -1 or hi >= 1):
            out.append(f"|tanh x| не меньше единицы: {iv.record()}")
        if not deep:
            return out
        if -(-self) != self:
            out.append("двойная смена знака не вернула запись")
        try:
            ch = HyperbolicValue("cosh", self.x).enclosure(30)
            sh = HyperbolicValue("sinh", self.x).enclosure(30)
            if ch.hi <= sh.lo:
                out.append(f"cosh x не больше sinh x: {ch.record()} против "
                           f"{sh.record()}")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"сверка cosh против sinh отказала: {type(e).__name__}")
        return out

    def record(self) -> str:
        return f"Φ(C=12 HYPERBOLIC, fn={self.fn}, x={self.x}, " \
               f"sign={self.sign})"

    __repr__ = record

    # --- равенство: канон делает слоты инвариантом ----------------------

    def __eq__(self, other):
        if isinstance(other, HyperbolicValue):
            if self.fn != other.fn:
                # разные функции линейки не совпадают на ℚ: равенство
                # свелось бы к исчезающей комбинации экспонент с
                # различными рациональными показателями [внешний
                # источник: Линдеман—Вейерштрасс]
                return False
            # одна функция: на каноне x > 0 каждая из sinh/cosh/tanh
            # строго монотонна — слоты и есть инвариант
            return self.x == other.x and self.sign == other.sign
        if isinstance(other, _NUM):
            # sinh/cosh/tanh ненулевого рационального трансцендентны
            # [внешний источник: Линдеман—Вейерштрасс]
            return False
        return NotImplemented

    def __hash__(self):
        return hash((HyperbolicValue, self.fn, self.x, self.sign))

    # --- знак: чётность в устройстве ------------------------------------

    def __neg__(self) -> "HyperbolicValue":
        # знак — компонент записи и не теряется: −(sign·fn(x)) несёт
        # тот же слот x с перевёрнутым sign; для нечётных fn это ровно
        # fn(−x), для cosh — честный носитель −cosh(x)
        return HyperbolicValue(self.fn, self.x, -self.sign)

    # --- сложение / вычитание: направленного правила нет ----------------

    def _refuse(self, op: str, other: object, self_left: bool):
        """Чужой класс с записью — дверь ФIOS (порядок операндов как в
        выражении); класс 1, одноимённая линейка и прибитые ячейки —
        прежний отказ носителя."""
        if _foreign_record(other) and not _pinned(other, op):
            a, b = (self, other) if self_left else (other, self)
            return _door(a, b, op)
        left, right = ((self.record(), _rec(other)) if self_left
                       else (_rec(other), self.record()))
        raise Unresolved(op, left, right,
                         hint="линейка независима: сложение сведено "
                              "только с нолём-состоянием; суммы "
                              "гиперболик не сведены")

    def __add__(self, other):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__radd__
        return self._refuse("+", other, True)

    def __radd__(self, other):
        return self._refuse("+", other, False)

    def __sub__(self, other):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__rsub__
        return self._refuse("-", other, True)

    def __rsub__(self, other):
        return self._refuse("-", other, False)

    # --- умножение / деление: направленного правила нет ------------------

    def __mul__(self, other):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__rmul__
        if _foreign_record(other) and not _pinned(other, "*"):
            return _door(self, other, "*")      # чужой класс — дверь ФIOS
        raise Unresolved("×", self.record(), _rec(other),
                         hint="слоты линейки — только (fn, x, sign): "
                              "множителя не несёт; произведения "
                              "гиперболик не сведены")

    # Сторона операнда не теряется: `__rmul__ = __mul__` печатал 2·x как
    # «x × 2» — операнды местами (замер 12.09, гл. 07 §19 п. 6). У соседей
    # по главе порядок верен (`bessel.py`, `de_solution.py`), и сложение
    # этого же модуля сторону уже различало. Отказ и подсказка те же —
    # названы лишь операнды как даны. Ноль-состояние слева сведёт
    # Zero.__mul__, потому NotImplemented остаётся.
    def __rmul__(self, other):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__mul__
        if _foreign_record(other) and not _pinned(other, "*"):
            return _door(other, self, "*")
        raise Unresolved("×", _rec(other), self.record(),
                         hint="слоты линейки — только (fn, x, sign): "
                              "множителя не несёт; произведения "
                              "гиперболик не сведены")

    def __truediv__(self, other):
        if isinstance(other, Zero) or (_foreign_record(other) and not _pinned(other, "/")):
            return _door(self, other, "/")      # ноль-состояние — таблица «x/0»
        raise Unresolved("/", self.record(), _rec(other),
                         hint="частное этой линейкой не сведено")

    def __rtruediv__(self, other):
        if _foreign_record(other) and not _pinned(other, "/"):
            return _door(other, self, "/")
        raise Unresolved("/", _rec(other), self.record(),
                         hint="обратная величина этой линейкой не несётся")

    def __pow__(self, n: int):
        if not isinstance(n, int) or isinstance(n, bool):
            return NotImplemented
        if n == 0:
            return Q(1)
        if n == 1:
            return HyperbolicValue(self.fn, self.x, self.sign)
        raise Unresolved("^", self.record(), _rec(n),
                         hint="степень вне {0, 1} не сведена; тождество "
                              "cosh²−sinh² — функцией pythagorean_check")

    # --- десятичная сторона: сертификат сетки ----------------------------

    def enclosure(self, terms: int = 30) -> Interval:
        """Сертификат сетки [lo, hi] ⊂ ℚ; уже́ние — явным запросом.

        sinh/cosh — интервальная полуразность/полусумма вложений
        exp_enclosure(x) и exp_enclosure(−x) (структурная пара над
        классом 5 дословно); tanh — интервальное частное sinh/cosh:
        cosh ≥ 1, потому ноль в знаменатель не попадает. Знак слота —
        точный множитель сертификата.
        """
        if isinstance(terms, bool):
            raise TypeError(f"terms — {_NOT_BOOL}; дверь принимает int")
        if not isinstance(terms, int):
            raise TypeError("terms — только int")
        if terms < 1:
            raise ValueError("сетка требует terms ≥ 1")
        half = Q(1, 2)
        e_pos = exp_enclosure(self.x, terms)
        e_neg = exp_enclosure(-self.x, terms)
        if self.fn == "sinh":
            return (e_pos - e_neg) * half * self.sign
        if self.fn == "cosh":
            return (e_pos + e_neg) * half * self.sign
        sinh_iv = (e_pos - e_neg) * half
        cosh_iv = (e_pos + e_neg) * half
        return sinh_iv / cosh_iv * self.sign


def hyperbolic(fn: str, x) -> Q | HyperbolicValue:
    """fn(x) по рациональному x: сведение направленными правилами.

    x = 0 — точное вырождение в класс 1 (sinh 0 = 0, cosh 0 = 1,
    tanh 0 = 0); чётность — каноном x > 0: у нечётных sinh/tanh знак
    аргумента уходит в слот sign, у чётного cosh — гасится; прочее —
    носитель HyperbolicValue.
    """
    fn, x = _check_fn_x(fn, x)
    if x == 0:
        return Q(1) if fn == "cosh" else Q(0)
    if x < 0:
        # чётность в фабрике: sinh(−x) = −sinh(x), tanh(−x) = −tanh(x),
        # cosh(−x) = cosh(x)
        return HyperbolicValue(fn, -x, -1 if fn in _ODD else 1)
    return HyperbolicValue(fn, x, 1)
