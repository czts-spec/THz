"""Линейка класса 11 ARCTRIG: arctan/arcsin/arccos рационального аргумента.

Величина хранится структурно: имя функции fn ∈ {arctan, arcsin, arccos}
и несводимый рациональный аргумент x. Направленные правила фабрики
arctrig():
  * табличные точки приземляются в класс 10 (моном q·π) либо в точный
    ноль класса 1; таблица исчерпывающая — рациональные значения
    тригонометрических функций в рациональных долях π суть
    {0, ±1/2, ±1} [внешний источник: теорема Нивена];
  * arcsin/arccos вне [−1, 1] — вещественного значения нет: честный
    отказ Unresolved («вне области»), не ложный класс;
  * прочие точки несводимы — носитель ArcTrigValue со слотами (fn, x).

ДВОЙСТВЕННОСТЬ с классом 10: одно направленное правило арифметики —
дополнение арктангенсов arctan x + arctan(1/x) = π/2 при x > 0
приземляет сумму в моном класса 10 ТОЧНО; его десятичная тень —
сертификат Мэчина 16·arctan(1/5) − 4·arctan(1/239), пересекающий
вложение π (приёмка proto/test_phicore_arctrig.py). Иных правил
сведения нет: сумма с рациональным, произведение, частное, степень
вне {0, 1} — отказ Unresolved с явными Ф-записями обоих операндов.
Вычитание сведено ровно на РАВНЫХ величинах и закрывается ИНВАРИАНТОМ
(четверть, sin²), не совпадением слотов: arcsin(3/5) − arctan(3/4) есть
точный ноль класса 1, потому что это один угол в двух записях (правка
12.09, оговорка к §19 п. 8 главы 07).

Равенство величин решается ТОЧНО инвариантом (четверть, sin²): квадрат
синуса рационален для всех трёх fn, а в пределах одной четверти синус
инъективен — так распознаётся и межфункциональное тождество вида
arcsin(3/5) = arctan(3/4) = arccos(4/5) (пифагорова тройка). Величина
нетабличного слота не равна ни одному рациональному: тангенс ненулевого
рационального иррационален [внешний источник: Ламберт], синус и косинус
ненулевого алгебраического трансцендентны [внешний источник: Линдеман].

Десятичная сторона — только сертификатом сетки (enclosure) по явному
запросу: arctan при |x| < 1 — знакочередующийся ряд (частные суммы
обрамляют предел, теорема Лейбница; счёт строгий — целыми
аккумуляторами на масштабе 2^−(g+запас) с направленными наружу
делениями, итог на сетке 2^−g — приём этапа C, enclosure.py), при
|x| ≥ 1 — дополнение ±π/2 − arctan(1/x) с вложением π по Мэчину;
arcsin — через arctan(x/√(1−x²)) корнем-сертификатом и интервальным
частным; arccos = π/2 − arcsin. Ни одного float: float на входе —
TypeError построением.
"""
from __future__ import annotations

from .rational import Q

from .enclosure import Interval, _atan_enclosure, pi_enclosure, sqrt_enclosure
from .errors import Unresolved
from .pi_value import PiMonomial
from .zero_state import Zero

#: Допустимые скаляры. float исключён построением.
_NUM = (int, Q)

#: `[вычислено здесь]` 10.10, очередь строка 105, шаг третий: дверь оператора
#: носителя = дверь ФIOS, как у монома π (шаг первый) и у e^y, ln (шаг второй).
#: Для ЧУЖОГО класса операторы линейки arctrig зовут `dispatch.combine` и
#: возвращают её ответ; класс 1 (int, Q) и одноимённая линейка — прежние
#: точные правила и собственные отказы носителя (ячейки таблицы refuse с той
#: же сутью). Ввоз ленивый, внутри функций: Фios ввозит arctrig сам.
_IN_DOOR: set = set()
#: Места `[вычислено здесь]` 10.10 — proto/test_phicore_arctrig.py:
#:   (10, ×) — «arctan(1/3)·π — отказ (правила сведения линеек нет)».
#: Дверь в той же ячейке даёт запись класса 19 (π·arctan(1/3) слева уже идёт
#: дверью монома); расхождение названо, снятие — с пересмотром приёмки.
_PINNED: frozenset = frozenset({(10, "*")})
_FOREIGN = object()


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


def _own(x) -> bool:
    """Класс 1 или одноимённая линейка — решает носитель сам."""
    return isinstance(x, _NUM) or isinstance(x, ArcTrigValue)


def _foreign_record(x) -> bool:
    return (not _own(x)) and callable(getattr(x, "record", None))
#: Одна запись отказа на bool во всех дверях модуля: bool числом в доме не
#: считается и называется bool'ом (канон `phicore/rational.py`; сосед по
#: главе отвечает тем же: `Notch.from_q(True, 8)` → «bool не число»).
#: Двух копий одного канона здесь нет намеренно: разошедшиеся копии — уже
#: пройденная поломка (`affine.py:75`).
_NOT_BOOL = "bool не число"

#: Имена функций линейки; иное имя — не слот.
_FNS = ("arctan", "arcsin", "arccos")

#: Таблица точных приземлений (fn, x) → q: значение равно q·π; q = 0 —
#: точный ноль класса 1. Таблица полна по теореме Нивена
#: [внешний источник: рациональные sin/cos/tan рациональных долей π].
_EXACT = {
    ("arctan", Q(0)): Q(0),
    ("arctan", Q(1)): Q(1, 4),
    ("arctan", Q(-1)): Q(-1, 4),
    ("arcsin", Q(0)): Q(0),
    ("arcsin", Q(1)): Q(1, 2),
    ("arcsin", Q(-1)): Q(-1, 2),
    ("arcsin", Q(1, 2)): Q(1, 6),
    ("arcsin", Q(-1, 2)): Q(-1, 6),
    ("arccos", Q(1)): Q(0),
    ("arccos", Q(0)): Q(1, 2),
    ("arccos", Q(-1)): Q(1),
    ("arccos", Q(1, 2)): Q(1, 3),
    ("arccos", Q(-1, 2)): Q(2, 3),
}


def _rec(x: object) -> str:
    """Ф-запись операнда для сообщений отказа."""
    r = getattr(x, "record", None)
    if callable(r):
        return r()
    if isinstance(x, _NUM):
        return f"Φ(C=1 RATIONAL, value={Q(x)})"
    return f"Φ(?={x!r})"


def _check_args(fn, x) -> tuple:
    """Общая проверка входов фабрики и конструктора: имя из линейки,
    аргумент — только int/Q; float и bool — TypeError построением;
    arcsin/arccos вне [−1, 1] — отказ («вне области»)."""
    if not isinstance(fn, str):
        raise TypeError("fn — только str")
    if fn not in _FNS:
        raise ValueError("fn ∈ {arctan, arcsin, arccos}; иного линейка "
                         "не несёт")
    if isinstance(x, bool):
        raise TypeError(f"x — {_NOT_BOOL}; дверь принимает int/Q")
    if not isinstance(x, _NUM):
        raise TypeError("x — только int/Q, float запрещён")
    x = Q(x)
    if fn != "arctan" and not -1 <= x <= 1:
        raise Unresolved(fn, _rec(x),
                         hint="вне области: |x| > 1 — вещественного "
                              "значения нет")
    return fn, x


# ---------- десятичная сторона: строгие вложения ----------

def _atan_series(x: Q, terms: int) -> Interval:
    """arctan x при |x| < 1: ряд x − x³/3 + x⁵/5 − … Члены строго
    убывают по модулю, потому соседние частные суммы обрамляют предел
    (теорема Лейбница). Счёт — целыми аккумуляторами на масштабе
    2^−(g+запас) с направленными наружу делениями (общий ход этапа C —
    _alternating из enclosure.py, тот же, что несёт π по Мэчину);
    итог — на сетке 2^−g: сертификат несёт точный Q-отрезок
    частных сумм, концы — со знаменателями-степенями двойки."""
    if x == 0:
        return Interval.point(0)
    return _atan_enclosure(x, terms)


def _atan_point(x: Q, terms: int) -> Interval:
    """arctan x для любого x ∈ ℚ: |x| < 1 — ряд; |x| = 1 — четверть
    вложения π; |x| > 1 — дополнение arctan x = ±π/2 − arctan(1/x)
    с вложением π по Мэчину (двойственность с классом 10)."""
    if -1 < x < 1:
        return _atan_series(x, terms)
    if x == 1:
        return pi_enclosure(terms) * Q(1, 4)
    if x == -1:
        return pi_enclosure(terms) * Q(-1, 4)
    half = pi_enclosure(terms) * Q(1, 2)
    inner = _atan_series(1 / x, terms)
    return half - inner if x > 1 else -half - inner


def _asin_point(x: Q, terms: int) -> Interval:
    """arcsin x при |x| < 1: arctan(x/√(1−x²)) — корень сертификатом
    сетки, интервальное частное, затем монотонность arctan по концам:
    нижний конец — низом вложения в t.lo, верхний — верхом в t.hi."""
    if x == 0:
        return Interval.point(0)
    root = sqrt_enclosure(1 - x * x, max(8 * terms, 160))
    t = Interval.point(x) / root
    return Interval(_atan_point(t.lo, terms).lo,
                    _atan_point(t.hi, terms).hi)


# ---------- встречная дорога арктангенса (ось самоприменения, 08.10) ----------

def _atan_road(x: Q, terms: int) -> Interval:
    """arctan x ВСТРЕЧНОЙ дорогой: формула сложения с единицей
    arctan y = π/4 + arctan((y − 1)/(y + 1)) при y = |x| > 0, аргумент ряда
    лежит в [0; 1) при всяком y > 0 — ряд и π входят ВСЕГДА, в отличие от
    сертификата (ряд при |x| < 1, дополнение ±π/2 при |x| > 1); знак — по
    нечётности [внешний источник: DLMF 4.24.2]. x = 0 — точка."""
    if x == 0:
        return Interval.point(0)
    y = abs(x)
    inner = _atan_series((y - 1) / (y + 1), terms)
    iv = pi_enclosure(terms) * Q(1, 4) + inner
    return iv if x > 0 else Interval(-iv.hi, -iv.lo)


def _atan_road_iv(t: Interval, terms: int) -> Interval:
    """Та же дорога на отрезке — монотонность arctan по концам."""
    return Interval(_atan_road(t.lo, terms).lo, _atan_road(t.hi, terms).hi)


# ---------- сама линейка ----------

class ArcTrigValue:
    """Несводимое fn(x): fn ∈ {arctan, arcsin, arccos}, x ∈ ℚ.

    Канон держится построением: табличные точки в слот не помещаются —
    их приземление делает фабрика arctrig(); arcsin/arccos вне [−1, 1]
    не строятся вовсе. Множителя линейка не несёт: слоты только (fn, x).
    """

    __slots__ = ("fn", "x")

    def __init__(self, fn, x):
        fn, x = _check_args(fn, x)
        if (fn, x) in _EXACT:
            raise ValueError("слот ArcTrigValue несёт только несводимую "
                             "точку; приземление табличных — фабрикой "
                             "arctrig()")
        self.fn, self.x = fn, x

    # --- запись ---------------------------------------------------------

    def record(self) -> str:
        return f"Φ(C=11 ARCTRIG, fn={self.fn}, x={self.x})"

    __repr__ = record

    # --- инвариант равенства --------------------------------------------

    def _quadrant(self) -> int:
        """Четверть угла: −1 — (−π/2, 0), 1 — (0, π/2), 2 — (π/2, π).

        Нетабличный слот на границы четвертей не попадает: точки
        x = 0 у всех трёх fn лежат в таблице."""
        if self.fn == "arccos":
            return 1 if self.x > 0 else 2
        return 1 if self.x > 0 else -1

    def _sin2(self) -> Q:
        """Квадрат синуса угла — рационален для всех трёх fn.

        Вместе с четвертью — точный инвариант равенства: в пределах
        четверти синус инъективен, знак его фиксирован."""
        if self.fn == "arctan":
            return self.x ** 2 / (1 + self.x ** 2)
        if self.fn == "arcsin":
            return self.x ** 2
        return 1 - self.x ** 2  # arccos: sin² = 1 − cos²

    def _same_value(self, other: "ArcTrigValue") -> bool:
        """Один ли это УГОЛ — по инварианту величины, не по слотам.

        Инвариант — пара (четверть, sin²). Он точен в обе стороны:
        равный sin² даёт равный |sin|, а в пределах одной четверти знак
        синуса постоянен и синус инъективен (на (−π/2; π/2) возрастает,
        на (π/2; π) убывает), потому равенство инварианта РАВНОСИЛЬНО
        равенству углов. Слоты (fn, x) — лишь одна из записей угла:
        одному углу их отвечает несколько (пифагорова тройка
        arcsin(3/5) = arctan(3/4) = arccos(4/5)), и совпадение слотов
        есть ДОСТАТОЧНОЕ, но не необходимое условие.

        Одна копия канона на модуль: ею отвечают и `__eq__`, и `__hash__`,
        и вычитание (правка 12.09, вердикт главы 07 §19 п. 8). Прежде
        `==` уже шло инвариантом, а вычитание отказывало при доказанном
        равенстве — два ответа на один вопрос.
        """
        return (self._quadrant() == other._quadrant()
                and self._sin2() == other._sin2())

    def __eq__(self, other):
        if isinstance(other, ArcTrigValue):
            # тождество углов, включая межфункциональное
            # (arcsin(3/5) = arctan(3/4)): одна четверть и один sin²
            return self._same_value(other)
        if isinstance(other, _NUM):
            # нетабличная точка не равна рациональному: тангенс
            # ненулевого рационального иррационален [внешний источник:
            # Ламберт], sin/cos ненулевого алгебраического
            # трансцендентны [внешний источник: Линдеман]
            return False
        return NotImplemented

    def __hash__(self):
        # хеш берётся с ИНВАРИАНТА, не со слотов: иначе равные величины
        # разошлись бы по корзинам (arcsin(3/5) и arctan(3/4))
        return hash((ArcTrigValue, self._quadrant(), self._sin2()))

    # --- сложение / вычитание -------------------------------------------

    def __add__(self, other):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__radd__
        if (isinstance(other, ArcTrigValue) and self.fn == "arctan"
                and other.fn == "arctan" and self.x > 0
                and other.x == 1 / self.x):
            # ДОПОЛНЕНИЕ АРКТАНГЕНСОВ: arctan x + arctan(1/x) = π/2
            # при x > 0 — приземление в класс 10 (двойственность линеек)
            return PiMonomial(Q(1, 2), Q(1))
        if _foreign_record(other) and not _pinned(other, "+"):
            return _door(self, other, "+")      # чужой класс — дверь ФIOS
        raise Unresolved("+", self.record(), _rec(other),
                         hint="линейка независима: сложение сведено только "
                              "с нолём-состоянием и дополнением "
                              "arctan x + arctan(1/x) при x > 0")

    def __radd__(self, other):
        if _foreign_record(other) and not _pinned(other, "+"):
            return _door(other, self, "+")      # порядок — как в выражении
        raise Unresolved("+", _rec(other), self.record(),
                         hint="линейка независима: сложение сведено только "
                              "с нолём-состоянием и дополнением "
                              "arctan x + arctan(1/x) при x > 0")

    _SUB_HINT = ("вычитание сведено только на РАВНЫХ величинах "
                 "(инвариант: четверть и sin²); иных правил линейка не несёт")

    def __sub__(self, other):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__rsub__
        if isinstance(other, ArcTrigValue) and self._same_value(other):
            # x − x закрывается ИНВАРИАНТОМ величины, а не совпадением
            # слотов (правка 12.09, оговорка к §19 п. 8 главы 07): слоты
            # у одного угла бывают разные. Замер — свод 8×8 = 64 пар на
            # носителях {atan ±3/4, asin ±3/5, acos ±4/5, atan 1/3,
            # asin 1/3}: совпадением слотов закрылось бы 8 пар, инвариант
            # закрывает 16 — вдвое больше; восемь добавленных суть
            # межфункциональные тождества пифагоровой тройки, и среди них
            # arcsin(3/5) − arctan(3/4), где
            # равенство величин доказано (`__eq__` отвечал True, а
            # вычитание отказывало). Итог — арифметический ноль класса 1
            # (метка 0{1} переписи `zeros.py`): сократилась ВЕЛИЧИНА,
            # а не форма записи.
            return Q(0)
        if _foreign_record(other) and not _pinned(other, "-"):
            return _door(self, other, "-")
        raise Unresolved("-", self.record(), _rec(other),
                         hint=self._SUB_HINT)

    def __rsub__(self, other):
        if isinstance(other, ArcTrigValue) and self._same_value(other):
            return Q(0)
        if _foreign_record(other) and not _pinned(other, "-"):
            return _door(other, self, "-")
        raise Unresolved("-", _rec(other), self.record(),
                         hint=self._SUB_HINT)

    def __neg__(self):
        if self.fn == "arccos":
            raise Unresolved("−", self.record(),
                             hint="arccos не нечётен: знак в слотах (fn, x) "
                                  "не несётся")
        # нечётность: −arctan x = arctan(−x), −arcsin x = arcsin(−x);
        # аргумент остаётся нетабличным — таблица симметрична по знаку
        return ArcTrigValue(self.fn, -self.x)

    # --- умножение / деление: направленного правила нет ------------------

    def __mul__(self, other):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__rmul__
        if _foreign_record(other) and not _pinned(other, "*"):
            return _door(self, other, "*")
        raise Unresolved("×", self.record(), _rec(other),
                         hint="слоты линейки — только (fn, x): множителя "
                              "не несёт, правила сведения нет")

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
                         hint="слоты линейки — только (fn, x): множителя "
                              "не несёт, правила сведения нет")

    def __truediv__(self, other):
        if isinstance(other, Zero) or (_foreign_record(other) and not _pinned(other, "/")):
            return _door(self, other, "/")   # ноль-состояние — таблица «x/0»
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
            return ArcTrigValue(self.fn, self.x)
        raise Unresolved("^", self.record(), _rec(n),
                         hint="степень вне {0, 1} не сведена")

    # --- десятичная сторона: сертификат сетки ----------------------------

    # --- самопроверка носителя (ось самоприменения, 08.10) ----------------

    def self_check(self, deep: bool = False) -> list:
        """Самопроверка носителя встречными дорогами (ось «AIAΣA спрашивает
        AIAΣA», слово автора 04.09; очередь строка 18; 08.10). Список
        находок, пусто — здоров; исключений не бросает.

        `[вычислено здесь]` (1) слоты: fn из трёх, x — ℚ, точка не табличная
        (иначе жила бы в классах 1/10), arcsin/arccos — |x| < 1; (2)
        инвариант: четверть согласна со знаком x, sin² ∈ (0; 1); (3)
        сертификат enclosure(40): концы упорядочены, ширина ≤ 2^−30, область —
        arctan и arcsin в (−π/2; π/2), arccos в (0; π); (4) ВСТРЕЧНАЯ дорога:
        формула сложения arctan y = π/4 + arctan((y−1)/(y+1)) (`_atan_road`),
        у которой ряд и π входят всегда, против сертификата (ряд либо
        дополнение ±π/2): для arctan — на самом x; для arcsin — на тангенсе
        x/√(1−x²) отрезком; для arccos — на √(1−x²)/x с добавкой π при x < 0,
        то есть мимо дороги π/2 − arcsin, которой идёт сертификат; отрезки
        обязаны пересекаться. При deep — ещё terms = 80 и ширина не растёт."""
        out = []
        fn, x = self.fn, self.x
        if fn not in _FNS:
            return [f"слот fn вне {{arctan, arcsin, arccos}}: {fn!r}"]
        if not isinstance(x, Q):
            return [f"слот x чужого носителя: {type(x).__name__}"]
        if (fn, x) in _EXACT:
            return [f"слот ({fn}, {x}) табличен — носителю не принадлежит (приземление фабрикой arctrig)"]
        if fn != "arctan" and not -1 < x < 1:
            return [f"слот x = {x} вне (−1; 1) при {fn}"]
        try:
            q, s2 = self._quadrant(), self._sin2()
            if fn == "arccos" and q != (1 if x > 0 else 2):
                out.append(f"четверть {q} не согласна со знаком x = {x} при arccos")
            if fn != "arccos" and q != (1 if x > 0 else -1):
                out.append(f"четверть {q} не согласна со знаком x = {x} при {fn}")
            if not 0 < s2 < 1:
                out.append(f"инвариант sin² = {s2} вне (0; 1)")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"инвариант не построен: {type(e).__name__}: {e}")
        try:
            iv = self.enclosure(40)
        except Exception as e:                                  # noqa: BLE001
            return out + [f"сертификат не построен: {type(e).__name__}: {e}"]
        if iv.lo > iv.hi:
            out.append(f"концы сертификата переставлены: {iv.record()}")
        if iv.width() > Q(1, 2 ** 30):
            out.append(f"ширина сертификата {iv.width()} > 2^−30")
        try:
            half = pi_enclosure(40) * Q(1, 2)
            if fn == "arccos":
                if iv.lo < 0 or iv.hi > half.hi * 2:
                    out.append(f"arccos вне (0; π): {iv.record()}")
            elif iv.lo < -half.hi or iv.hi > half.hi:
                out.append(f"{fn} вне (−π/2; π/2): {iv.record()}")
            if fn == "arctan":
                road = _atan_road(x, 40)
            else:
                root = sqrt_enclosure(1 - x * x, 320)
                if fn == "arcsin":
                    road = _atan_road_iv(Interval.point(x) / root, 40)
                else:
                    road = _atan_road_iv(root / Interval.point(x), 40)
                    if x < 0:
                        road = road + pi_enclosure(40)
            if road.intersect(iv) is None:
                out.append(f"встречная дорога сложения arctan ∌ сертификат при {fn}({x}): {road.record()} против {iv.record()}")
            if deep:
                iv2 = self.enclosure(80)
                if iv2.intersect(iv) is None or iv2.width() > iv.width():
                    out.append(f"сертификат на 80 членах не вложен или шире: {iv2.record()} против {iv.record()}")
        except Unresolved as e:
            out.append(f"встречная дорога отказала при {fn}({x}): {e}")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"встречная дорога не построена при {fn}({x}): {type(e).__name__}: {e}")
        return out

    def enclosure(self, terms: int = 40) -> Interval:
        """Сертификат сетки [lo, hi] ⊂ ℚ; уже́ние — явным запросом.

        arctan: |x| < 1 — знакочередующийся ряд (обрамление Лейбница);
        |x| > 1 — дополнение ±π/2 − arctan(1/x) с вложением π по
        Мэчину. arcsin: arctan(x/√(1−x²)) — корень сертификатом сетки
        и интервальное частное; ноль в интервале корня даёт честный
        отказ, не приближение. arccos: π/2 − arcsin.
        """
        if isinstance(terms, bool):
            raise TypeError(f"terms — {_NOT_BOOL}; дверь принимает int")
        if not isinstance(terms, int):
            raise TypeError("terms — только int")
        if terms < 1:
            raise ValueError("сетка требует terms ≥ 1")
        if self.fn == "arctan":
            return _atan_point(self.x, terms)
        if self.fn == "arcsin":
            return _asin_point(self.x, terms)
        return (pi_enclosure(terms) * Q(1, 2)
                - _asin_point(self.x, terms))


def arctrig(fn: str, x) -> Q | PiMonomial | ArcTrigValue:
    """fn(x) по рациональному x: сведение направленными правилами.

    Табличные точки — точное приземление (моном q·π класса 10 либо ноль
    класса 1); arcsin/arccos вне [−1, 1] — отказ («вне области»);
    прочее — носитель ArcTrigValue.
    """
    fn, x = _check_args(fn, x)
    q = _EXACT.get((fn, x))
    if q is not None:
        return Q(0) if q == 0 else PiMonomial(q, Q(1))
    return ArcTrigValue(fn, x)


# ---------- самопроверка линейки (ось самоприменения, 08.10) ----------

def self_check(deep: bool = False) -> list:
    """Самопроверка линейки класса 11 встречными дорогами; список находок,
    пусто — здорова; исключений не бросает.

    `[вычислено здесь]` Окно носителей: arctan 1/2, arctan 2, arctan −3,
    arcsin 3/5, arccos 4/5, arcsin −1/3 (deep — ещё arctan 1/5, arccos −2/3)
    проверяют себя (`ArcTrigValue.self_check`); таблица Нивена приземляет
    arctan 1, arcsin 1/2, arccos 1/2, arccos 0 в мономы π и arcsin 0 в нуль;
    arcsin 2 — отказ вне области; пифагорова тройка arcsin 3/5 = arctan 3/4 =
    arccos 4/5 — одно значение инвариантом; дополнение arctan 2 + arctan 1/2
    приземляется в π/2 точно; Мэчин: 16·arctan(1/5) − 4·arctan(1/239)
    сертификатами носителей несёт π [внешний источник: формула Мэчина] —
    дорога класса 11 против вложения π класса 10. Граница — окно."""
    out = []
    pts = [("arctan", Q(1, 2)), ("arctan", Q(2)), ("arctan", Q(-3)), ("arcsin", Q(3, 5)), ("arccos", Q(4, 5)), ("arcsin", Q(-1, 3))]
    if deep:
        pts += [("arctan", Q(1, 5)), ("arccos", Q(-2, 3))]
    for fn, x in pts:
        try:
            v = arctrig(fn, x)
            if isinstance(v, ArcTrigValue):
                out += [f"{fn}({x}): {f}" for f in v.self_check(deep)]
            else:
                out.append(f"{fn}({x}) приземлилось в {v!r} без носителя")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"{fn}({x}): {type(e).__name__}: {e}")
    try:
        table = [(("arctan", Q(1)), PiMonomial(Q(1, 4), Q(1))), (("arcsin", Q(1, 2)), PiMonomial(Q(1, 6), Q(1))),
                 (("arccos", Q(1, 2)), PiMonomial(Q(1, 3), Q(1))), (("arccos", Q(0)), PiMonomial(Q(1, 2), Q(1)))]
        for (fn, x), want in table:
            if arctrig(fn, x) != want:
                out.append(f"таблица: {fn}({x}) = {arctrig(fn, x)!r} ≠ {want!r}")
        if arctrig("arcsin", Q(0)) != Q(0):
            out.append("таблица: arcsin 0 не точный нуль")
        try:
            arctrig("arcsin", Q(2))
            out.append("arcsin 2 не отказал вне области")
        except Unresolved:
            pass
        a, b, c = arctrig("arcsin", Q(3, 5)), arctrig("arctan", Q(3, 4)), arctrig("arccos", Q(4, 5))
        if not (a == b == c):
            out.append("пифагорова тройка arcsin 3/5, arctan 3/4, arccos 4/5 не распознана одним углом")
        if arctrig("arctan", Q(2)) + arctrig("arctan", Q(1, 2)) != PiMonomial(Q(1, 2), Q(1)):
            out.append("дополнение arctan 2 + arctan 1/2 не приземлилось в π/2")
        m = arctrig("arctan", Q(1, 5)).enclosure(60) * 16 - arctrig("arctan", Q(1, 239)).enclosure(60) * 4
        if m.intersect(pi_enclosure(40)) is None:
            out.append(f"Мэчин: {m.record()} не несёт π")
    except Exception as e:                                      # noqa: BLE001
        out.append(f"таблица и тождества: {type(e).__name__}: {e}")
    return out
