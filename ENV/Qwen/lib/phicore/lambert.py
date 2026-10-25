"""Линейка класса 17 LAMBERT_W: обращение x·e^x = y, ветви 0 и −1.

Величина хранится структурно: аргумент y — ТОЧНАЯ величина дома, номер
ветви ∈ {0, −1}. Вход y объявлен тремя носителями, и только ими:

  * пара ℚ (int/Q);
  * класс 5 EXPONENTIAL — `ExpValue` = coeff·e^y;
  * класс 3 POLY_ROOT — `AlgebraicNumber` (рациональный корень
    приземляется в ℚ ещё в двери, `as_fraction`);
  * ноль-состояние абсолютный (класс 0) — то же, что y = 0.

Прочее — класс 6 LOGARITHMIC, класс 10 (π-моном), класс 14 GAMMA,
вложенный класс 17 — получает отказ Unresolved ПОИМЁННО, с путём
`phicore/lambert.py` и с явной Ф-записью носителя: точной дороги нет, и
приближения на её место не ставится. float/bool — TypeError построением.

ТОЧНАЯ ДОРОГА (`exact_w`, без единого сертификата):

  * W(0) = 0 на ветви 0;
  * W(k·e^k) = k при k ∈ ℚ — потому что k·e^k есть в точности величина
    класса 5 `ExpValue(k, k)`, и тождество двойственности замыкается
    арифметикой линейки 5: `w_exp_w(k)` == y как объекты, БЕЗ разрядов.
    Ветвь 0 несёт k ≥ −1, ветвь −1 несёт k ≤ −1 — это область значений
    ветвей, не соглашение.

Отсюда достижима точка ветвления: y = −1/e есть `ExpValue(−1, −1)`, и
W(−1/e) = −1 ТОЧНО на обеих ветвях (k = −1 — единственная точка, где
ветви сходятся). Рационального вырождения при y ∈ ℚ, y ≠ 0 линейка не
несёт: будь w = W(y) рационально и ненулево, произведение w·e^w было бы
трансцендентно [внешний источник, теорема Линдемана] и не равнялось бы
рациональному y. По той же причине точного W нет и при иррациональном
алгебраическом y (класс 3).

Правил сведения W с независимыми величинами нет — сумма с числом,
произведение, частное, степень вне {0, 1} получают отказ Unresolved
с явными Ф-записями обоих операндов, не ложный класс и не молчание.

ДВОЙСТВЕННОСТЬ с классом 5 EXPONENTIAL: определяющее тождество
W(y)·e^{W(y)} = y — обратный ход экспоненциальной линейки. Где точная
дорога есть, тождество проверяется ТОЧНО (`w_exp_w` против самого y);
где её нет — сертификатом: произведение вложений W и e^W обязано
нести y.

Граница вещественной области — точка ветвления −1/e. Где y не равен ей
точной дорогой, принадлежность разбирается сертификатом сетки
inv_e = exp_enclosure(−1, grid): заведомо ниже границы — отказ
«ниже -1/e», зазор сертификата — честный отказ «уточнить сетку»
(глубина сетки названа параметром grid и поднимается по явному
запросу), не догадка о стороне.

Десятичная сторона — только сертификатом сетки (enclosure), по явному
запросу: бисекция f(w) = w·e^w − y, где e^w на рациональном w несёт
exp_enclosure, y несёт свой сертификат (точка для ℚ), а знак f
сертифицирован, когда отрезок w·e^w не пересекает отрезок y.
Ни одного float: float на входе — TypeError построением.
"""
from __future__ import annotations

from .rational import Q

from .algebraic import AlgebraicNumber
from .enclosure import Interval, exp_enclosure, ln_enclosure
from .errors import Unresolved
from .exp_log import ExpValue
from .zero_state import Zero, ZeroKind

#: Допустимые скаляры. float исключён построением.
_NUM = (int, Q)

#: `[вычислено здесь]` 10.10, очередь строка 105, шаг девятый: дверь оператора
#: носителя = дверь ФIOS (как у восьми прежних носителей). Для ЧУЖОГО класса с
#: записью единый `_refuse` ведёт в `dispatch.combine` с порядком операндов
#: как в выражении; класс 1 (int, Q) и одноимённая линейка — прежний отказ
#: носителя. Ввоз ленивый, внутри функций.
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
    return (not isinstance(x, (_NUM, WValue))) and callable(getattr(x, "record", None))

#: Носители y, у которых есть точная либо сертификатная дорога. Список
#: закрыт: всё вне его отказывает поимённо с путём модуля.
_Y_CARRIERS = (int, Q, ExpValue, AlgebraicNumber)

#: Объявленная глубина сертификата границы −1/e. Зазор на этой глубине —
#: отказ «уточнить сетку»; поднять глубину можно только явно (grid=…).
_BOUND_GRID = 40

#: Потолок углубления сертификата знака: глубже — честный отказ.
_MAX_TERMS = 1 << 13

#: Потолок скана стартовой вилки ветви −1 (шаг 1/4 от −1 вниз).
_MAX_SCAN = 1 << 12

#: Путь модуля в сообщениях отказа: отказ называет своё место.
_HERE = "phicore/lambert.py"


def _rec(x: object) -> str:
    """Ф-запись операнда для сообщений отказа."""
    r = getattr(x, "record", None)
    if callable(r):
        return r()
    if isinstance(x, _NUM):
        return f"Φ(C=1 RATIONAL, value={Q(x)})"
    return f"Φ(?={x!r})"


# ---------- дверь аргумента: величины дома, прочее — поимённо ----------

def _admit_y(y, branch):
    """Вход y: точная величина дома → носитель ℚ / класс 5 / класс 3.

    Приземления в двери: int → Q; рациональный корень класса 3 → Q
    (`as_fraction`); абсолютный ноль-состояние → Q(0). Прочий носитель —
    отказ Unresolved поимённо, с путём модуля; float и bool — TypeError
    построением (bool числом в доме не считается)."""
    if isinstance(y, bool):
        raise TypeError("y — bool числом в доме не считается; дверь несёт "
                        "int/Q, класс 5 ExpValue, класс 3 AlgebraicNumber")
    if isinstance(y, float):
        raise TypeError("y — float запрещён построением; дверь несёт int/Q, "
                        "класс 5 ExpValue, класс 3 AlgebraicNumber")
    if isinstance(y, _NUM):
        return Q(y)
    if isinstance(y, AlgebraicNumber):
        r = y.as_fraction()
        return y if r is None else Q(r)
    if isinstance(y, ExpValue):
        return y
    if isinstance(y, Zero):
        if y.kind is ZeroKind.ABSOLUTE:
            return Q(0)
        raise Unresolved("lambert_w", _rec(y), _rec(branch),
                         hint=f"{_HERE}: ноль-состояние с провенансом за y не "
                              "идёт — точной дороги к W нет (абсолютный ноль "
                              "сведён как y = 0)")
    raise Unresolved("lambert_w", _rec(y), _rec(branch),
                     hint=f"{_HERE}: точной дороги к W для этого носителя "
                          "нет; y несут пара ℚ, класс 5 EXPONENTIAL "
                          "(ExpValue) и класс 3 POLY_ROOT (AlgebraicNumber)")


def _y_record(y) -> str:
    """Запись y внутри Ф-записи класса 17: пара ℚ — самим числом, иная
    точная величина дома — своей Ф-записью."""
    return str(y) if isinstance(y, Q) else y.record()


def _y_sign(y) -> int:
    """Знак y ТОЧНО, по его собственной линейке.

    ℚ — сравнением; класс 5 — знаком coeff (e^y > 0 при любом y ∈ ℚ);
    класс 3 — сертификатом знака самой линейки 3 (`sign`)."""
    if isinstance(y, Q):
        return (1 if y > 0 else -1) if y != 0 else 0
    if isinstance(y, ExpValue):
        return 1 if y.coeff > 0 else -1
    return y.sign()


def _y_enclosure(y, grid: int) -> Interval:
    """Сертификат сетки для самого y: точка для ℚ, вложение линейки — иначе.

    Класс 3 уточняет свой изолирующий отрезок на месте (`enclosure` →
    `refine`): величина от этого не меняется, отрезок только у́же — канон
    дома (`phicore/algebraic.py`)."""
    if isinstance(y, Q):
        return Interval.point(y)
    return y.enclosure(grid)


# ---------- точная дорога: тождество, а не сертификат ----------

def w_exp_w(w):
    """w·e^w ТОЧНОЙ дорогой дома — левая часть тождества двойственности.

    w = 0 → точный ноль; w ≠ 0 → величина класса 5 `ExpValue(w, w)`:
    коэффициент и показатель суть одна и та же пара ℚ. Ни одного
    сертификата — сравнение с y идёт объектами."""
    if isinstance(w, bool) or not isinstance(w, _NUM):
        raise TypeError("w — только int/Q, float запрещён")
    w = Q(w)
    return Q(0) if w == 0 else ExpValue.make(w, w)


def exact_w(y, branch: int = 0):
    """Точное W(y) ∈ ℚ, где точная дорога есть; иначе None.

    Дорог ровно две, и обе суть тождества:

      * y = 0, ветвь 0 → W = 0 (0·e⁰ = 0);
      * y = k·e^k, то есть класс 5 `ExpValue(k, k)` → W = k, при k ≥ −1 на
        ветви 0 и k ≤ −1 на ветви −1 (область значений ветвей). Точка
        ветвления k = −1 лежит на обеих ветвях: W(−1/e) = −1.

    Иных рациональных значений линейка не несёт: k·e^k при k ∈ ℚ, k ≠ 0
    трансцендентно [внешний источник, теорема Линдемана], потому ни
    рациональный, ни иррациональный алгебраический y точного W не даёт."""
    if isinstance(y, bool) or isinstance(y, float):
        return None
    if isinstance(y, _NUM):
        return Q(0) if (Q(y) == 0 and branch == 0) else None
    if isinstance(y, AlgebraicNumber):
        r = y.as_fraction()
        return Q(0) if (r is not None and r == 0 and branch == 0) else None
    if isinstance(y, ExpValue) and y.coeff == y.y:
        k = y.coeff
        if branch == 0 and k >= -1:
            return k
        if branch == -1 and k <= -1:
            return k
    return None


# ---------- сертификатный разбор границы −1/e ----------

def _certify_branch_point(y, branch: int, grid: int) -> None:
    """Разбор y относительно точки ветвления −1/e сертификатом сетки.

    Граница иррациональна, потому сравнение ведётся не с точкой, а с её
    вложением: inv_e = exp_enclosure(−1, grid) несёт 1/e отрезком
    [lo; hi], значит −1/e ∈ [−hi; −lo]; сам y несёт свой отрезок
    Y = _y_enclosure(y, grid) (точка при y ∈ ℚ). Разбор ровно таков:

      Y.hi < −hi        — заведомо ниже −1/e: отказ;
      Y.lo < −lo        — зона неразличимости сертификата: честный отказ
                          с подсказкой уже́ния (grid=…), не догадка о
                          стороне;
      иначе             — заведомо y > −1/e: допустимо.

    Точное равенство y = −1/e сюда не доходит: оно снято точной дорогой
    (`exact_w` на `ExpValue(−1, −1)`) ещё в фабрике."""
    inv_e = exp_enclosure(-1, grid)
    Y = _y_enclosure(y, grid)
    if Y.hi < -inv_e.hi:
        raise Unresolved("lambert_w", _rec(y), _rec(branch),
                         hint="ниже -1/e — вещественного значения на "
                              "ветвях 0 и −1 нет")
    if Y.lo < -inv_e.lo:
        raise Unresolved("lambert_w", _rec(y), _rec(branch),
                         hint=f"y в зоне неразличимости сертификата "
                              f"границы -1/e на сетке grid={grid} — "
                              f"уточнить сетку (grid > {grid})")


# ---------- сертификат знака f(w) = w·e^w − y ----------

def _certified_sign(w: Q, y, terms: int) -> "tuple[int, int]":
    """Знак f(w) = w·e^w − y по сертификату: (знак ±1, глубина terms).

    Вложение w·e^w — произведение точки w на exp_enclosure(w, terms);
    y несёт свой отрезок (точку при y ∈ ℚ). Знак сертифицирован, когда
    отрезки не пересекаются. При w = 0 значение f(0) = −y, и знак берётся
    точным знаком самого y (носитель гарантирует y ≠ 0).

    Углубление завершимо на всех трёх носителях y, и это заявление о
    предмете, не о вере [внешний источник, теорема Линдемана; сборка —
    вычислено здесь]:
      * y ∈ ℚ: w·e^w при w ≠ 0 трансцендентно, рациональному y не равно;
      * y = c·e^d (класс 5): w·e^w = c·e^d дало бы e^{w−d} = c/w ∈ ℚ,
        значит w = d и c = w = d — а такой y снят точной дорогой;
      * y алгебраический иррациональный: w·e^w трансцендентно.
    Потолок _MAX_TERMS даёт честный отказ, не догадку."""
    if w == 0:
        s = _y_sign(y)
        return (1 if s < 0 else -1), terms
    while True:
        p = exp_enclosure(w, terms) * w
        Y = _y_enclosure(y, terms)
        if p.lo > Y.hi:
            return 1, terms
        if p.hi < Y.lo:
            return -1, terms
        if terms >= _MAX_TERMS:
            raise Unresolved("sign(w·e^w − y)", _rec(w), _rec(y),
                             hint=f"знак не сертифицирован до глубины "
                                  f"{terms} — дальнейшее уже́ние отклонено")
        terms *= 2


# ---------- сама линейка ----------

class WValue:
    """W(y) на ветви branch: y — точная величина дома (ℚ, класс 5 или
    класс 3), ветвь ∈ {0, −1}, grid — объявленная глубина сертификата
    границы −1/e (в тождество величины не входит).

    Носитель держит ЛИШЬ те y, у которых точного W нет: точная дорога
    приземляется в ℚ фабрикой `lambert_w`, и прямое построение такого y
    здесь — ValueError с названным точным значением."""

    __slots__ = ("y", "branch", "grid")

    def __init__(self, y, branch: int = 0, grid: int = _BOUND_GRID):
        if isinstance(branch, bool) or not isinstance(branch, int):
            raise TypeError("branch — только int, float запрещён")
        if isinstance(grid, bool) or not isinstance(grid, int):
            raise TypeError("grid — только int")
        if grid < 1:
            raise ValueError("grid ≥ 1")
        y = _admit_y(y, branch)
        if branch not in (0, -1):
            raise Unresolved("lambert_w", _rec(y), _rec(branch),
                             hint="вещественны только ветви 0 и −1; "
                                  "иная ветвь не сведена")
        if branch == 0 and isinstance(y, Q) and y == 0:
            raise ValueError("y=0 — W(0) есть обычный ноль; стройте "
                             "через фабрику lambert_w")
        k = exact_w(y, branch)
        if k is not None:
            raise ValueError(f"W({_y_record(y)}) = {k} ТОЧНО — величина живёт "
                             "в ℚ, а не в носителе класса 17; стройте через "
                             "фабрику lambert_w")
        s = _y_sign(y)
        if branch == 0:
            if s < 0:
                _certify_branch_point(y, branch, grid)
        else:
            if s >= 0:
                raise Unresolved("lambert_w", _rec(y), _rec(branch),
                                 hint="ветвь −1 несёт лишь −1/e ≤ y < 0")
            _certify_branch_point(y, branch, grid)
        self.y, self.branch, self.grid = y, branch, grid

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
                         hint="W(y) независим: направленного правила "
                              "арифметики нет")

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
            # W ≠ 0 на носителе: ветвь 0 несёт y ≠ 0, ветвь −1 лежит
            # в w ≤ −1 — потому x⁰ = 1 точно
            return Q(1)
        if n == 1:
            return WValue(self.y, self.branch, self.grid)
        raise Unresolved("^", self.record(), _rec(n),
                         hint="степень W вне {0, 1} не сведена")

    def __neg__(self):
        raise Unresolved("neg", self.record(),
                         hint="противоположной величины линейка не несёт")

    # --- запись, сравнение, сертификат ---------------------------------

    def record(self) -> str:
        return (f"Φ(C=17 LAMBERT_W, ветвь={self.branch}, "
                f"y={_y_record(self.y)})")

    __repr__ = record

    def __eq__(self, other):
        if isinstance(other, WValue):
            # grid есть глубина сертификата, не часть величины
            return self.y == other.y and self.branch == other.branch
        if isinstance(other, _NUM):
            # W(y) на носителе иррационально (см. шапку модуля): точные
            # рациональные значения сняты фабрикой и здесь не живут
            return False
        return NotImplemented

    def __hash__(self):
        return hash((WValue, self.y, self.branch))

    def self_check(self, deep: bool = False) -> list:
        """Самопроверка носителя встречной дорогой (ось «AIAΣA спрашивает
        AIAΣA», слово автора 04.09; очередь строка 18; 08.10). Список находок,
        пусто — здоров; исключений не бросает.

        `[вычислено здесь]` (1) слоты: ветвь ∈ {0, −1}, y — ℚ, класс 5 либо
        класс 3, точной дороги у y нет (иначе величина жила бы в ℚ);
        (2) сертификат enclosure(16): ширина ≤ 2^−16, область ветви — W ≥ −1
        на ветви 0, W ≤ −1 на ветви −1 [внешний источник: DLMF 4.13(i)];
        (3) ВСТРЕЧНАЯ дорога логарифмом, не экспонентой: из W·e^W = y следует
        ln|W| + W = ln|y| [DLMF 4.13.1], и отрезок ln|W| + W по концам
        сертификата (ряд ln, `ln_enclosure`) обязан пересекать ln|y| — у ℚ
        и класса 3 ряд ln по модулю, у класса 5 c·e^d — ln|c| + d точно по
        показателю; сертификат же строился ТОЛЬКО экспонентой, потому дороги
        независимы. При deep — тождество фабрики: lambert_w(k·e^k, ветвь) = k
        для k своей ветви (точная дорога, без сертификата)."""
        out = []
        br, y = self.branch, self.y
        if isinstance(br, bool) or br not in (0, -1):
            return [f"слот branch вне {{0, −1}}: {br!r}"]
        if not isinstance(y, (Q, ExpValue, AlgebraicNumber)):
            return [f"слот y чужого носителя: {type(y).__name__}"]
        if exact_w(y, br) is not None:
            out.append(f"y имеет точное W = {exact_w(y, br)} — носителю не принадлежит")
        try:
            iv = self.enclosure(16)
        except Exception as e:                                  # noqa: BLE001
            return out + [f"сертификат не построен: {type(e).__name__}: {e}"]
        if iv.width() > Q(1, 2 ** 16):
            out.append(f"ширина сертификата {iv.width()} > 2^−16")
        if br == 0 and iv.lo < -1:
            out.append(f"ветвь 0 ниже −1: lo = {iv.lo}")
        if br == -1 and iv.hi > -1:
            out.append(f"ветвь −1 выше −1: hi = {iv.hi}")
        try:
            if iv.lo <= 0 <= iv.hi:
                out.append(f"сертификат несёт ноль: {iv.record()} (W = 0 лишь при y = 0)")
            else:
                a, b = abs(iv.lo), abs(iv.hi)
                lo_w, hi_w = min(a, b), max(a, b)
                lhs = Interval(ln_enclosure(lo_w, 40).lo, ln_enclosure(hi_w, 40).hi) + iv
                if isinstance(y, ExpValue):
                    c = abs(Q(y.coeff))
                    rhs = ln_enclosure(c, 40) + Q(y.y)
                else:
                    Y = _y_enclosure(y, 40)
                    if Y.lo <= 0 <= Y.hi:
                        out.append("y не отделён от нуля — ln|y| не строится")
                        rhs = None
                    else:
                        ya, yb = abs(Y.lo), abs(Y.hi)
                        rhs = Interval(ln_enclosure(min(ya, yb), 40).lo, ln_enclosure(max(ya, yb), 40).hi)
                if rhs is not None and lhs.intersect(rhs) is None:
                    out.append(f"встречная дорога ln|W| + W ∌ ln|y| на ветви {br}: {lhs.record()} против {rhs.record()}")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"встречная дорога ln не построена: {type(e).__name__}: {e}")
        if deep:
            for k in ((Q(1, 2), Q(2)) if br == 0 else (Q(-2), Q(-3, 2))):
                try:
                    if lambert_w(w_exp_w(k), br) != k:
                        out.append(f"тождество фабрики lambert_w(k·e^k, {br}) ≠ k при k = {k}")
                except Exception as e:                          # noqa: BLE001
                    out.append(f"тождество фабрики при k = {k}: {type(e).__name__}: {e}")
        return out

    def enclosure(self, bits: int = 32) -> Interval:
        """Сертификат сетки [lo, hi] ⊂ ℚ ширины ≤ 2^−bits, по явному
        запросу.

        Бисекция f(w) = w·e^w − y: e^w на рациональном w несёт
        exp_enclosure, y несёт свой отрезок, знак f сертифицирован, когда
        отрезки не пересекаются (_certified_sign). Стартовая вилка:

          ветвь 0:  [−1, max(1, Y.hi) + 1] — слева f(−1) = −1/e − y < 0
                    (y правее границы по сертификату носителя), справа
                    w·e^w > w > y при w = max(1, Y.hi) + 1 ≥ 2;
          ветвь −1: верх −1 (там f < 0), низ — скан шага 1/4 от −1
                    вниз до сертифицированной смены знака: w·e^w → 0
                    при w → −∞, потому f > 0 достигается на глубине
                    порядка −ln(−y); потолок скана — честный отказ.
        """
        if isinstance(bits, bool) or not isinstance(bits, int):
            raise TypeError("bits — только int")
        if bits < 1:
            raise ValueError("bits ≥ 1")
        y = self.y
        terms = 24
        if self.branch == 0:
            lo = Q(-1)
            hi = max(Q(1), _y_enclosure(y, self.grid).hi) + 1
            s_lo, terms = _certified_sign(lo, y, terms)
        else:
            hi = Q(-1)
            lo, s_lo = None, 1
            w, step = hi, Q(1, 4)
            for _ in range(_MAX_SCAN):
                w -= step
                sign, terms = _certified_sign(w, y, terms)
                if sign > 0:
                    lo = w
                    break
                hi = w
            if lo is None:
                raise Unresolved("enclosure", self.record(),
                                 hint=f"скан шага 1/4 до {w} не дал смены "
                                      "знака — вилка ветви −1 не "
                                      "изолирована")
        target = Q(1, 2 ** bits)
        while hi - lo > target:
            mid = (lo + hi) / 2
            s_mid, terms = _certified_sign(mid, y, terms)
            if s_mid == s_lo:
                lo = mid
            else:
                hi = mid
        return Interval(lo, hi)


# ---------- самопроверка линейки (ось самоприменения, 08.10) ----------

def self_check(deep: bool = False) -> list:
    """Самопроверка линейки класса 17 встречными дорогами; список находок,
    пусто — здорова; исключений не бросает.

    `[вычислено здесь]` Окно: y ∈ {1, 2, 1/2, 3/7} на ветви 0, y ∈ {−1/4, −1/10}
    на обеих ветвях (deep — ещё y = 5, −3/10 и класс 5 y = 2·e^{1/2}):
    носители проверяют себя логарифмической дорогой (`WValue.self_check`);
    точные точки — W(0) = 0, W(−1/e) = −1 на обеих ветвях, W(k·e^k) = k —
    фабрикой без сертификата [DLMF 4.13.1]. Граница — окно."""
    out = []
    pts = [(Q(1), 0), (Q(2), 0), (Q(1, 2), 0), (Q(3, 7), 0), (Q(-1, 4), 0), (Q(-1, 4), -1), (Q(-1, 10), 0), (Q(-1, 10), -1)]
    if deep:
        pts += [(Q(5), 0), (Q(-3, 10), 0), (Q(-3, 10), -1), (ExpValue.make(Q(2), Q(1, 2)), 0)]
    for y, br in pts:
        try:
            v = lambert_w(y, br)
            if isinstance(v, WValue):
                out += [f"W({_y_record(y)}, {br}): {f}" for f in v.self_check(deep)]
            else:
                out.append(f"W({_y_record(y)}, {br}) приземлилось в {v!r} без носителя")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"W({_y_record(y)}, {br}): {type(e).__name__}: {e}")
    try:
        if lambert_w(0) != Q(0):
            out.append("W(0) ≠ 0")
        bp = ExpValue.make(Q(-1), Q(-1))
        if lambert_w(bp, 0) != Q(-1) or lambert_w(bp, -1) != Q(-1):
            out.append("W(−1/e) ≠ −1 на одной из ветвей")
        for k, br in ((Q(3), 0), (Q(-1, 3), 0), (Q(-5, 2), -1)):
            if lambert_w(w_exp_w(k), br) != k:
                out.append(f"W(k·e^k) ≠ k при k = {k}, ветвь {br}")
    except Exception as e:                                      # noqa: BLE001
        out.append(f"точные точки: {type(e).__name__}: {e}")
    return out


# ---------- фабрика линейки ----------

def lambert_w(y, branch: int = 0, grid: int = _BOUND_GRID):
    """W(y) по направленным правилам линейки: обращение x·e^x = y.

    Вход y — точная величина дома: пара ℚ, класс 5 EXPONENTIAL, класс 3
    POLY_ROOT (рациональный корень приземляется в ℚ), абсолютный
    ноль-состояние. Иной носитель — отказ поимённо с путём модуля.

    точная дорога   → ТОЧНОЕ значение в ℚ (`exact_w`): W(0) = 0 и
                      W(k·e^k) = k, в том числе точка ветвления
                      W(−1/e) = −1 на обеих ветвях; тождество
                      W·e^W = y замыкается объектами (`w_exp_w`);
    ветвь 0         → допустимо y ≥ −1/e; где точного равенства с
                      границей нет, она разбирается сертификатом
                      (_certify_branch_point): заведомо ниже — отказ
                      «ниже -1/e», зазор сетки grid — отказ «уточнить
                      сетку»; допустимый y — носитель WValue(y, 0);
    ветвь −1        → нужен −1/e ≤ y < 0: тот же сертификатный разбор;
                      y ≥ 0 — отказ;
    иные ветви      → отказ Unresolved: вещественных значений нет.
    """
    if isinstance(branch, bool) or not isinstance(branch, int):
        raise TypeError("branch — только int, float запрещён")
    if isinstance(grid, bool) or not isinstance(grid, int):
        raise TypeError("grid — только int")
    if grid < 1:
        raise ValueError("grid ≥ 1")
    y = _admit_y(y, branch)
    k = exact_w(y, branch)
    if k is not None:
        return k
    return WValue(y, branch, grid)
