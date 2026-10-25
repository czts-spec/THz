"""Линейка класса 13 ZETA: ζ(s) рационального аргумента.

Величина хранится структурно: несводимый слот s (Q). Направленные
правила фабрики zeta():
  * s = 1 — полюс: честный отказ Unresolved, значения нет;
  * s = 0 — точное Q(−1/2); целые s < 0 — тривиальная сторона:
    ζ(−n) = −B_{n+1}/(n+1) из чисел Бернулли, чётные n дают точный ноль
    класса 1 [внешний источник: тривиальные нули ζ];
  * чётные s = 2n ≥ 2 — приземление в класс 10 (двойственность линеек):
    ζ(2n) = (−1)^{n+1}·B_{2n}·(2π)^{2n}/(2·(2n)!) — моном q·π^{2n}
    [внешний источник: формула Эйлера для ζ(2n)];
  * нечётные s ≥ 3 и нецелые s — несводимы: ZetaValue со слотом s.

Числа Бернулли считаются точной рекуррентой Σ_{k≤n} C(n+1, k)·B_k = 0
с кэшем; таблиц значений нет, биномы — мультипликативной рекуррентой.

ДВОЙСТВЕННОСТЬ с классами 14 и 10: функциональное уравнение
ζ(s) = 2^s·π^{s−1}·sin(πs/2)·Γ(1−s)·ζ(1−s) связывает тривиальную
сторону с чётной через sin_pi и GammaValue; тождество при s = −1
проверяется точно приёмкой (proto/test_phicore_zeta.py).

Правил сведения ζ с независимыми величинами нет: ζ(3)+1, 2·ζ(3) и т.п.
получают отказ Unresolved с явными Ф-записями операндов, не ложный
класс; сведена только разность одноимённого s — точный ноль класса 1.
Десятичная сторона — только сертификатом сетки (enclosure) по явному
запросу: частная сумма плюс строгое обрамление хвоста интегралом.

КРИТИЧЕСКАЯ ПРЯМАЯ s = 1/2 + it (12.09, строка 77 очереди). Точно:
множитель функционального уравнения χ(s) = ζ(s)/ζ(1−s) сводится в
величину дома при целом s ≤ 0 и при ПОЛУЦЕЛОМ s (functional_factor),
и χ(1/2) = 1 РОВНО — уравнение на критической прямой вырождается в
тождество и о величине не сообщает ничего. Сертификатом: ζ(1/2) и
вообще нецелое s > 0 — дорогой DLMF 25.2.8 (zeta_enclosure_integral),
сверка — знакочередующимся рядом (zeta_enclosure_eta) и при s > 1
прямым рядом (zeta_enclosure_dirichlet); прежний отказ «нецелый слот
сертификата не имеет» снят замером. Отказ об ИНСТРУМЕНТЕ: точка
1/2 + it при t ≠ 0 — слот линейки есть ℚ, цена дороги названа в
_pair_refusal. Отказа о ПРЕДМЕТЕ дом здесь не заявляет: ζ(1/2) есть
величина класса 13 построением, а несводимость её к классам 1, 3, 10
есть открытый вопрос математики, а не ноль инструмента.
"""
from __future__ import annotations

from .rational import Q

from .algebraic import AlgebraicNumber
from .enclosure import Interval, int_root_floor, pi_enclosure
from .errors import Unresolved
from .gamma_value import GammaValue, gamma_half_integer, sin_pi
from .pi_value import PiMonomial
from .zero_state import Zero

#: Допустимые скаляры. float исключён построением.
_NUM = (int, Q)

#: `[вычислено здесь]` 10.10, очередь строка 105, шаг четвёртый: дверь оператора
#: носителя = дверь ФIOS, как у π, e^y/ln и arctrig. Для ЧУЖОГО класса с
#: записью операторы ζ зовут `dispatch.combine`; класс 1 (int, Q) и
#: одноимённая линейка — прежние правила и собственные отказы носителя.
#: Ввоз ленивый, внутри функций: Фios ввозит zeta сам.
_IN_DOOR: set = set()
#: Места `[вычислено здесь]` 10.10 — proto/test_phicore_zeta.py:
#:   (10, ×) — «ζ(3)·π — отказ (правила сведения линеек нет)».
#: Дверь в той же ячейке даёт запись класса 19 (π·ζ(3) слева уже идёт дверью
#: монома); расхождение названо, снятие — с пересмотром приёмки.
_PINNED: frozenset = frozenset({(10, "*")})


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
    return (not isinstance(x, (_NUM, ZetaValue))) and callable(getattr(x, "record", None))
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


def _exact_zero(x: object) -> bool:
    """Точный ноль-ЗНАЧЕНИЕ класса 1 на входе умножения (не bool).

    `[вычислено здесь]` разнобой 12.09: сосед по главе уже держал это
    правило — `gamma_value.GammaValue.__mul__` на q == 0 отдаёт Q(0)
    («масштаб нолём — точный ноль класса 1»), а ζ отвечала отказом даже
    нолю. Довод тот же и для ζ: канон слота пускает только s ≠ 1, значит
    величина КОНЕЧНА, и произведение с нолём есть ноль, а не «не сведено».
    Без этого правила функциональное уравнение не доносит тривиальные
    нули: χ(−2) = 0 точно, но 0·ζ(3) упиралось в отказ.
    """
    return isinstance(x, _NUM) and not isinstance(x, bool) and Q(x) == 0


def _factorial(n: int) -> int:
    """n! точной рекуррентой k! = k·(k−1)!; таблиц значений нет."""
    out = 1
    for k in range(2, n + 1):
        out *= k
    return out


def _binomial(n: int, k: int) -> int:
    """C(n, k) мультипликативной рекуррентой C(n, j) = C(n, j−1)·(n−j+1)/j.

    Деление на каждом шаге нацело — точность держится построением,
    промежуточных дробей нет.
    """
    if k < 0 or k > n:
        return 0
    k = min(k, n - k)
    out = 1
    for j in range(1, k + 1):
        out = out * (n - j + 1) // j
    return out


#: Кэш рекурренты: B_0 = 1 — база; дальше список растёт по запросу.
_BERNOULLI: list[Q] = [Q(1)]


def bernoulli(n: int) -> Q:
    """B_n точной рекуррентой Σ_{k=0}^{n} C(n+1, k)·B_k = 0, с кэшем.

    База B_0 = 1; далее B_n = −(Σ_{k<n} C(n+1, k)·B_k)/(n+1). Конвенция
    B_1 = −1/2 не постулируется, а выводится этой же рекуррентой;
    таблицы значений нет.
    """
    if isinstance(n, bool):
        raise TypeError(f"n — {_NOT_BOOL}; дверь принимает int")
    if not isinstance(n, int):
        raise TypeError("n — только int")
    if n < 0:
        raise ValueError("индекс числа Бернулли неотрицателен")
    while len(_BERNOULLI) <= n:
        m = len(_BERNOULLI)
        acc = Q(0)
        for k in range(m):
            acc += _binomial(m + 1, k) * _BERNOULLI[k]
        _BERNOULLI.append(-acc / (m + 1))
    return _BERNOULLI[n]


class ZetaValue:
    """Несводимое ζ(s): слот s — нечётное целое ≥ 3 либо нецелое ≠ 1.

    Канон держится построением: сводимые точки (полюс s = 1, целые
    s ≤ 0, чётные s ≥ 2) в слот не помещаются — их сведение делает
    фабрика zeta(). Множителя линейка не несёт: слот только s.
    """

    __slots__ = ("s",)

    def __init__(self, s):
        if isinstance(s, float):
            raise TypeError("s — только int/Q, float запрещён")
        if hasattr(s, "re") and hasattr(s, "im"):
            _pair_refusal(s)       # 1/2 + it при t ≠ 0 — отказ поимённо
        if not isinstance(s, _NUM):
            raise TypeError("s — только int/Q")
        s = Q(s)
        if s.denominator == 1:
            n = s.numerator
            if n <= 0 or n == 1 or n % 2 == 0:
                raise ValueError("слот ZetaValue несёт только несводимое s "
                                 "(нечётное целое ≥ 3 или нецелое); "
                                 "сведение — фабрикой zeta()")
        self.s = s

    # --- запись ---------------------------------------------------------

    def record(self) -> str:
        return f"Φ(C=13 ZETA, s={self.s})"

    __repr__ = record

    def __eq__(self, other):
        if isinstance(other, ZetaValue):
            return self.s == other.s
        # равенство ζ(s) с рациональным или мономом π этой линейкой не
        # решается (арифметическая природа ζ нечётных и нецелых s в
        # общем случае не сведена) — ответ не выдумывается
        return NotImplemented

    def __hash__(self):
        return hash((ZetaValue, self.s))

    # --- сложение / вычитание -------------------------------------------

    def __add__(self, other):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__radd__
        if _foreign_record(other) and not _pinned(other, "+"):
            return _door(self, other, "+")      # чужой класс — дверь ФIOS
        raise Unresolved("+", self.record(), _rec(other),
                         hint="ζ независима: сложение сведено только "
                              "с нолём-состоянием")

    def __radd__(self, other):
        if _foreign_record(other) and not _pinned(other, "+"):
            return _door(other, self, "+")      # порядок — как в выражении
        raise Unresolved("+", _rec(other), self.record(),
                         hint="ζ независима: сложение сведено только "
                              "с нолём-состоянием")

    def __sub__(self, other):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__rsub__
        if isinstance(other, ZetaValue) and other.s == self.s:
            return Q(0)  # ζ(s) − ζ(s) — точный ноль класса 1
        if _foreign_record(other) and not _pinned(other, "-"):
            return _door(self, other, "-")
        raise Unresolved("-", self.record(), _rec(other),
                         hint="ζ независима: сведена только разность "
                              "одноимённого s")

    def __rsub__(self, other):
        if _foreign_record(other) and not _pinned(other, "-"):
            return _door(other, self, "-")
        raise Unresolved("-", _rec(other), self.record(),
                         hint="ζ независима: сведена только разность "
                              "одноимённого s")

    def __neg__(self):
        raise Unresolved("−", self.record(),
                         hint="слот линейки — только s: знака без "
                              "множителя эта запись не несёт")

    # --- умножение / деление: направленного правила нет ------------------

    def __mul__(self, other):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__rmul__
        if _exact_zero(other):
            return Q(0)            # масштаб нолём — точный ноль класса 1
        if _foreign_record(other) and not _pinned(other, "*"):
            return _door(self, other, "*")
        raise Unresolved("×", self.record(), _rec(other),
                         hint="слот линейки — только s: множителя ζ не "
                              "несёт, правила сведения нет")

    # Сторона операнда не теряется: `__rmul__ = __mul__` печатал 2·x как
    # «x × 2» — операнды местами (замер 12.09, гл. 07 §19 п. 6). У соседей
    # по главе порядок верен (`bessel.py`, `de_solution.py`), и сложение
    # этого же модуля сторону уже различало. Отказ и подсказка те же —
    # названы лишь операнды как даны. Ноль-состояние слева сведёт
    # Zero.__mul__, потому NotImplemented остаётся.
    def __rmul__(self, other):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__mul__
        if _exact_zero(other):
            return Q(0)            # масштаб нолём — точный ноль класса 1
        if _foreign_record(other) and not _pinned(other, "*"):
            return _door(other, self, "*")
        raise Unresolved("×", _rec(other), self.record(),
                         hint="слот линейки — только s: множителя ζ не "
                              "несёт, правила сведения нет")

    def __truediv__(self, other):
        if isinstance(other, Zero) or (_foreign_record(other) and not _pinned(other, "/")):
            return _door(self, other, "/")   # ноль-состояние — таблица «x/0»
        raise Unresolved("/", self.record(), _rec(other),
                         hint="частное с ζ не сведено")

    def __rtruediv__(self, other):
        if _foreign_record(other) and not _pinned(other, "/"):
            return _door(other, self, "/")
        raise Unresolved("/", _rec(other), self.record(),
                         hint="обратная величина 1/ζ этой линейкой "
                              "не несётся")

    def __pow__(self, n: int):
        if not isinstance(n, int) or isinstance(n, bool):
            return NotImplemented
        if n == 0:
            return Q(1)
        if n == 1:
            return ZetaValue(self.s)
        raise Unresolved("^", self.record(), _rec(n),
                         hint="степень ζ вне {0, 1} не сведена")

    # --- десятичная сторона: сертификат сетки ----------------------------

    def enclosure(self, terms: int = 60, bits: int = 128) -> Interval:
        """Сертификат сетки [lo, hi] ⊂ ℚ; уже́ние — явным запросом (рост N).

        Дороги по слоту, каждая со своей границей:
          * ЦЕЛОЕ s ≥ 3 — частная сумма Σ_{n≤N} n^{−s} плюс строгое
            обрамление остатка интегралом: R = Σ_{n>N} n^{−s} лежит в
            [1/((s−1)·(N+1)^{s−1}); 1/((s−1)·N^{s−1})] [внешний источник:
            интегральный признак]. Счёт весь в ℚ — корней не нужно;
          * НЕЦЕЛОЕ s > 0 (в том числе s = 1/2 — единственная точка
            критической прямой, которую несёт слот линейки) — с 02.10
            `zeta_enclosure_euler_maclaurin`: DLMF 25.2.9 с остатком,
            обрамлённым первым отброшенным членом (n = min(16, 2N/3)
            членов Бернулли). До 12.09 здесь стоял отказ «нецелый слот
            сертификата не имеет» — снят замером дорогой DLMF 25.2.8
            (`zeta_enclosure_integral`, ≈1,5 разряда на декаду N); с
            02.10 та дорога и `zeta_enclosure_eta` (при s > 1 ещё
            `zeta_enclosure_dirichlet`) — независимые сверки, а ответ
            даёт быстрая: при N = 60 — сетка 2^−128 против 2^−14;
          * НЕЦЕЛОЕ s < 0 — с 02.10 той же быстрой дорогой напрямую (DLMF
            25.2.9 держит ℜs > −2n, число членов Бернулли ≥ ⌊−s/2⌋ + 2);
            до того — точный множитель `functional_factor(s)` (при
            знаменателе 2 это ровно q·π^m) на сертификат зеркальной точки
            ζ(1−s), а знаменатель вне {1, 2} отвечал отказом «точного
            множителя нет» — теперь это сверка, не единственная дорога.
        `[вычислено здесь]` цена на s = 1/2: N = 10⁴ → 6 разрядов за
        0,03 с; N = 10⁶ → 9 разрядов за 3,1 с (≈1,5 разряда на декаду N).
        """
        if isinstance(terms, bool) or isinstance(bits, bool):
            raise TypeError(f"N — {_NOT_BOOL}; дверь принимает int")
        if not isinstance(terms, int) or not isinstance(bits, int):
            raise TypeError("N — только int")
        if terms < 1:
            raise ValueError("сетка требует N ≥ 1")
        if self.s.denominator != 1:
            # 02.10: и s < 0 — напрямую (DLMF 25.2.9 держит ℜs > −2n);
            # зеркало через точный χ осталось сверкой в приёмках
            return zeta_enclosure_euler_maclaurin(self.s, terms, bits)
        s = self.s.numerator
        n_grid = terms
        partial = Q(0)
        for n in range(1, n_grid + 1):
            partial += Q(1, n ** s)
        lo_tail = Q(1, (s - 1) * (n_grid + 1) ** (s - 1))
        hi_tail = Q(1, (s - 1) * n_grid ** (s - 1))
        return Interval(partial + lo_tail, partial + hi_tail)


    # --- самопроверка носителя (ось самоприменения, 08.10) ----------------

    def self_check(self, deep: bool = False) -> list:
        """Самопроверка носителя встречными дорогами (ось «AIAΣA спрашивает
        AIAΣA», слово автора 04.09; очередь строка 18; 08.10). Список
        находок, пусто — здоров; исключений не бросает.

        `[вычислено здесь]` (1) слот: s — ℚ, несводимое (нечётное целое ≥ 3
        либо нецелое ≠ 1), иначе величина жила бы в ℚ или классе 10;
        (2) сертификат enclosure(60): концы упорядочены, ширина ≤ 2^−14;
        (3) ВСТРЕЧНАЯ дорога по слоту: при s > 0 — знакочередующийся ряд
        η [внешний источник: DLMF 25.2.3] (`zeta_enclosure_eta`), независимый
        и от прямого ряда, которым идёт целое s, и от Эйлера—Маклорена,
        которым идёт нецелое; при s < 0 с знаменателем 2 — зеркало точным
        множителем χ(s) на η-сертификат точки 1 − s [DLMF 25.4.2]; при
        s < 0 с иным знаменателем точного χ нет — граница объявлена, дорога
        не строится. Отрезки обязаны пересекаться. При deep — сверка ещё
        дорогой DLMF 25.2.8 (`zeta_enclosure_integral`) при s > 0."""
        out = []
        s = self.s
        if not isinstance(s, Q):
            return [f"слот s чужого носителя: {type(s).__name__}"]
        if s == 1 or (s.denominator == 1 and (s.numerator <= 0 or s.numerator % 2 == 0)):
            return [f"слот s = {s} сводим — носителю не принадлежит (сведение фабрикой zeta)"]
        try:
            iv = self.enclosure(60)
        except Exception as e:                                  # noqa: BLE001
            return out + [f"сертификат не построен: {type(e).__name__}: {e}"]
        if iv.lo > iv.hi:
            out.append(f"концы сертификата переставлены: {iv.record()}")
        if iv.width() > Q(1, 2 ** 14):
            out.append(f"ширина сертификата {iv.width()} > 2^−14")
        try:
            road = None
            if s > 0:
                road = ("η", zeta_enclosure_eta(s, 400, 64))
            elif s.denominator == 2:
                chi = functional_factor(s)
                road = ("χ·η(1−s)", _factor_interval(chi, 64) * zeta_enclosure_eta(1 - s, 400, 64))
            if road is not None and road[1].intersect(iv) is None:
                out.append(f"встречная дорога {road[0]} ∌ сертификат при s = {s}: {road[1].record()} против {iv.record()}")
            if deep and s > 0:
                j = zeta_enclosure_integral(s, 2000, 64)
                if j.intersect(iv) is None:
                    out.append(f"дорога DLMF 25.2.8 ∌ сертификат при s = {s}: {j.record()} против {iv.record()}")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"встречная дорога не построена при s = {s}: {type(e).__name__}: {e}")
        return out


def zeta(s) -> Q | PiMonomial | ZetaValue:
    """ζ(s) по рациональному s: сведение направленными правилами.

    Полюс s = 1 — отказ; s = 0 и целые s < 0 — точные Q из
    Бернулли; чётные s ≥ 2 — моном класса 10; прочее — ZetaValue.
    """
    if isinstance(s, float):
        raise TypeError("s — только int/Q, float запрещён")
    if hasattr(s, "re") and hasattr(s, "im"):
        _pair_refusal(s)           # 1/2 + it при t ≠ 0 — отказ поимённо
    if not isinstance(s, _NUM):
        raise TypeError("s — только int/Q")
    s = Q(s)
    if s == 1:
        raise Unresolved("ζ(s)", _rec(s),
                         hint="полюс ζ: в точке s=1 значения нет")
    if s == 0:
        # ζ(0) = −1/2 [внешний источник]; сведено отдельной точкой:
        # формула ζ(−n) = −B_{n+1}/(n+1) при конвенции B_1 = −1/2
        # (её даёт рекуррента) в точке n = 0 меняет знак
        return Q(-1, 2)
    if s.denominator == 1:
        n = s.numerator
        if n < 0:
            m = -n
            # тривиальная сторона: ζ(−m) = −B_{m+1}/(m+1);
            # чётные m дают точный ноль класса 1 (B нечётных ≥ 3 нулевые)
            return -bernoulli(m + 1) / (m + 1)
        if n % 2 == 0:
            # приземление в класс 10: ζ(2k) = (−1)^{k+1}·B_{2k}·
            # (2π)^{2k}/(2·(2k)!) [внешний источник: формула Эйлера]
            k = n // 2
            sign = 1 if k % 2 == 1 else -1  # (−1)^{k+1}
            coeff = (sign * bernoulli(2 * k) * Q(2 ** (2 * k))
                     / (2 * _factorial(2 * k)))
            return PiMonomial(coeff, Q(n))
        # нечётное целое ≥ 3 — несводимое
    return ZetaValue(s)


# ══════════════════════════════════════════════════════════════════════
#  КРИТИЧЕСКАЯ ПРЯМАЯ s = 1/2 + it (строка 77 очереди, 12.09)
#
#  Разделение проведено замером, а не предположением:
#   ТОЧНО    — множитель функционального уравнения χ(s) при целом s ≤ 0
#              и при ПОЛУЦЕЛОМ s; в частности χ(1/2) = 1 — уравнение на
#              критической прямой есть тождество ζ(1/2) = ζ(1/2) и о
#              величине не говорит НИЧЕГО (отказ о ПРЕДМЕТЕ для этой
#              дороги, доказанный счётом, а не назначенный);
#   СЕРТИФИКАТОМ — сама величина ζ(1/2) и вообще ζ(s) при нецелом s > 0:
#              отрезок с объявленной разрядностью, две независимые дороги;
#   ОТКАЗ ОБ ИНСТРУМЕНТЕ — точка 1/2 + it при t ≠ 0: слот линейки есть ℚ,
#              комплексного s он не несёт (цена названа в _pair_refusal);
#   ОТКАЗ О ПРЕДМЕТЕ ЗАЯВЛЕН НЕ БУДЕТ — «ζ(1/2) не есть величина никакого
#              класса» дом НЕ доказывает и заявлять не станет: ζ(1/2) есть
#              величина класса 13 построением (слот 1/2 законен), а
#              несводимость её к классам 1, 3, 10 есть ОТКРЫТЫЙ вопрос
#              математики, а не ноль инструмента.
# ══════════════════════════════════════════════════════════════════════

def _pair_refusal(s) -> None:
    """Названный отказ на пару (re, im) в слоте вместо голой двери.

    Цена дороги названа числом, чтобы отказ не читался как «нельзя»:
    каждый член n^{−s} при s = σ + it есть n^{−σ}·(cos(t·ln n) −
    i·sin(t·ln n)), то есть НА ЧЛЕН — один сертификат ln n (класс 8,
    `enclosure.ln_enclosure`) и два тригонометрических
    (`sin_enclosure`/`cos_enclosure`, и те лишь при |x| ≤ 4, приведения
    периодом в доме нет), плюс своё обрамление хвоста для комплексного
    показателя, которого ни одна дорога ниже не даёт.
    """
    raise Unresolved("ζ(s)", _rec(s),
                     hint="слот класса 13 есть ℚ: точка критической "
                          "прямой 1/2 + it при t ≠ 0 записью этой линейки "
                          "не является. Вещественная её точка ровно одна — "
                          "s = 1/2 (неподвижная точка отражения s ↔ 1−s в "
                          "ℚ), и она в слот помещается. Цена комплексной "
                          "дороги: на каждый член ряда — сертификат ln n "
                          "и два тригонометрических, приведения периодом "
                          "у sin/cos дома нет, обрамления хвоста при "
                          "комплексном показателе нет вовсе")


def functional_factor(s):
    """χ(s) = ζ(s)/ζ(1−s) ТОЧНО, как величина дома (Q либо класс 10).

    χ(s) = 2^s·π^{s−1}·sin(πs/2)·Γ(1−s) [внешний источник: DLMF 25.4.2,
    http://dlmf.nist.gov/25.4.E2]. Точно сводится в двух полосах:

      * целое s ≤ 0: 2^s и Γ(1−s) = (−s)! рациональны, sin(πs/2) ∈
        {0, ±1}; чётное s ≤ −2 даёт ТОЧНЫЙ НОЛЬ — вторая, независимая
        от чисел Бернулли дорога к тривиальным нулям ζ;
      * s = m + 1/2 (знаменатель 2): гасятся ДВЕ пары сразу —
        2^s·sin(πs/2) = ±2^m (класс 3 схлопывается в ℚ: √2·(√2/2) = 1)
        и π^{s−1}·Γ(1−s) = q·π^m (Γ(1/2) = π^{1/2}); итог — моном q·π^m
        с ЦЕЛЫМ m = (2s−1)/2. При m = 0, то есть НА КРИТИЧЕСКОЙ ПРЯМОЙ,
        χ(1/2) = 1 ровно, и уравнение вырождается в ζ(1/2) = ζ(1/2).

    Граница названа, а не назначена:
      * целое s ≥ 1 — полюс Γ(1−s) (при чётном s ещё и ноль синуса:
        неопределённость 0·∞); уравнение идёт в обратную сторону;
      * s = 0 — χ(0) = 0 РОВНО как значение множителя (ноль синуса при
        Γ(1) = 1), но зеркальная точка ζ(1) — полюс: уравнение
        ζ(0) = χ(0)·ζ(1) здесь есть 0·∞ и ζ(0) им не сводится —
        ζ(0) = −1/2 берётся в `zeta` по Бернулли (наблюдение сдачи
        25.09, очередь строка 103: прежде точка s = 0 проходила под
        именем «тривиальный ноль» молча);
      * знаменатель s вне {1, 2} — 2^s есть класс 3 степени q, а Γ(1−s)
        остаётся несведённой величиной класса 14: пары для гашения нет.
        Это отказ О ПРЕДМЕТЕ настолько, насколько замкнутой формы Γ(p/q)
        при q ∉ {1, 2} не знает и математика (Chowla–Selberg даёт лишь
        произведения при мнимом квадратичном поле) — дом же не знает её
        ТОЧНО и потому молчит, а не выдумывает.
    """
    if isinstance(s, float):
        raise TypeError("s — только int/Q, float запрещён")
    if hasattr(s, "re") and hasattr(s, "im"):
        _pair_refusal(s)
    if not isinstance(s, _NUM):
        raise TypeError("s — только int/Q")
    s = Q(s)
    q = s.denominator
    if q == 1:
        n = s.numerator
        if n >= 1:
            raise Unresolved("χ(s)", _rec(s),
                             hint="целое s ≥ 1: Γ(1−s) — полюс"
                                  + (", а sin(πs/2) — ноль: 0·∞ не сведено"
                                     if n % 2 == 0 else "")
                                  + "; уравнение сводит обратную сторону")
        sn = sin_pi(Q(n, 2))
        if sn == 0:
            # чётное s ≤ −2: тривиальный ноль ζ — χ(s) = 0 вместе с ζ(s).
            # s = 0: ноль синуса ПРИ ПОЛЮСЕ ζ(1) — χ(0) = 0 как значение
            # множителя 2^s·π^{s−1}·sin(πs/2)·Γ(1−s) (Γ(1) = 1 конечна),
            # а уравнение ζ(0) = χ(0)·ζ(1) в этой точке — 0·∞ и ζ(0) им
            # не сводится (см. докстринг; ζ(0) = −1/2 — Бернулли в `zeta`).
            return Q(0)
        coeff = Q(1, 2 ** (-n)) * sn * Q(_factorial(-n))
        return PiMonomial._make(coeff, s - 1)
    if q == 2:
        m = (s.numerator - 1) // 2      # s = m + 1/2, m ∈ ℤ
        two = AlgebraicNumber.from_root(2, 2) * (Q(2 ** m) if m >= 0
                                                 else Q(1, 2 ** (-m)))
        two_part = two * sin_pi(s / 2)  # ±2^m — класс 3 схлопнулся в ℚ
        # π^{s−1}·Γ(1−s) = q·π^m. Подстановка Γ(1/2) = π^{1/2} живёт в ОДНОМ
        # месте — `gamma_value.gamma_half_integer`, где она и выведена двумя
        # дорогами (квадрат отражением Эйлера + положительность сертификатом).
        # Прежде та же подстановка стояла копией здесь; разошедшиеся копии —
        # пройденная поломка (`affine.py:75`).
        pi_part = PiMonomial(Q(1), s - 1) * gamma_half_integer(Q(1) - s)
        return two_part * pi_part
    raise Unresolved("χ(s)", _rec(s),
                     hint=f"знаменатель слота {q} вне {{1, 2}}: 2^s есть "
                          f"класс 3 степени {q}, Γ(1−s) остаётся величиной "
                          "класса 14 — пары для точного гашения нет; "
                          "замкнутой формы Γ(p/q) при q ∉ {1,2} дом не знает")


def functional_factor_square(s):
    """χ(s)² — ВТОРАЯ дорога к тому же множителю, независимая от первой.

    Первая дорога (`functional_factor`) пользуется равенством
    Γ(1/2) = π^{1/2}. Эта — не пользуется им вовсе: Γ(1−s)² берётся
    ОТРАЖЕНИЕМ ЭЙЛЕРА самого дома (`gamma_value.GammaValue.__mul__`:
    Γ(b)·Γ(1−b) = π/sin(πb) при b = 1/2 даёт ровно π), а sin²(πs/2) —
    возведением величины класса 3 в квадрат. Совпадение
    functional_factor(s)² == functional_factor_square(s) проверяется
    приёмкой: подстановка Γ(1/2) = π^{1/2} есть СЛЕДСТВИЕ устройства
    дома, а не внесённое допущение.
    """
    s = Q(s)
    if s.denominator == 1:
        chi = functional_factor(s)      # целая полоса: дорога одна
        # прежде здесь стояло условие с одинаковыми ветвями (мёртвое —
        # наблюдение сдачи 25.09, очередь строка 103); квадрат один для Q
        # и для монома π, оператор выбирает дорогу сам
        return chi * chi
    if s.denominator != 2:
        return functional_factor(s)     # отказ тем же текстом
    g = GammaValue.of(Q(1) - s)
    gg = g * g                          # отражение Эйлера: q·π^1
    sn = sin_pi(s / 2)
    two2 = Q(2) ** s.numerator if s.numerator >= 0 \
        else Q(1, 2 ** (-s.numerator))  # 2^{2s} = 2^p
    return gg * (sn * sn) * two2 * PiMonomial(Q(1), 2 * s - 2)


# ---------- сертификаты: целый счёт на сетке 2^−bits ----------

def _grid_pow_neg(n: int, p: int, q: int, bits: int) -> tuple:
    """Пара целых границ n^{−p/q} на сетке 2^{−bits}; p, q ≥ 1, n ≥ 1.

    Ни одного float: n^{−p/q}·2^bits = (2^{bits·q}/n^p)^{1/q}, целый
    корень — `enclosure.int_root_floor`. Низ есть floor, верх — floor+1,
    и верх ДОКАЗАН: будь (lo+1)^q < 2^{bits·q}/n^p, целое (lo+1)^q не
    превосходило бы ⌊2^{bits·q}/n^p⌋, откуда lo+1 ≤ ⌊t^{1/q}⌋ = lo —
    противоречие. Цена округления — не более узла сетки на конец.
    """
    lo = int_root_floor((1 << (bits * q)) // n ** p, q)
    return lo, lo + 1


def _grid_pow(n: int, p: int, q: int, bits: int) -> tuple:
    """Пара целых границ n^{p/q} на сетке 2^{−bits}; тот же довод верха."""
    lo = int_root_floor((n ** p) << (bits * q), q)
    return lo, lo + 1


def _grid_pow_signed(n: int, p: int, q: int, bits: int) -> tuple:
    """Пара целых границ n^{−p/q} при ЛЮБОМ знаке p: p > 0 — обратная
    степень (`_grid_pow_neg`), p ≤ 0 — прямая (`_grid_pow` с |p|); нужна
    дороге Эйлера—Маклорена под нулём (02.10)."""
    return _grid_pow_neg(n, p, q, bits) if p > 0 else _grid_pow(n, -p, q, bits)


def _partial_sum(n_from: int, n_to: int, p: int, q: int, bits: int) -> tuple:
    """Пара целых границ Σ_{n=n_from}^{n_to} n^{−p/q} на сетке 2^{−bits}
    (p любого знака)."""
    lo = hi = 0
    for n in range(n_from, n_to + 1):
        a, b = _grid_pow_signed(n, p, q, bits)
        lo += a
        hi += b
    return lo, hi


def _slot_doors(s, terms, bits) -> Q:
    """Общие двери сертификатных дорог: bool, float, тип, сетка."""
    if isinstance(terms, bool) or isinstance(bits, bool):
        raise TypeError(f"N и bits — {_NOT_BOOL}; дверь принимает int")
    if not isinstance(terms, int) or not isinstance(bits, int):
        raise TypeError("N и bits — только int")
    if terms < 1:
        raise ValueError("сетка требует N ≥ 1")
    if bits < 8:
        raise ValueError("разрядность сетки bits ≥ 8")
    if isinstance(s, float):
        raise TypeError("s — только int/Q, float запрещён")
    if hasattr(s, "re") and hasattr(s, "im"):
        _pair_refusal(s)
    if not isinstance(s, _NUM):
        raise TypeError("s — только int/Q")
    return Q(s)


def zeta_enclosure_integral(s, terms: int = 1000,
                            bits: int = 128) -> Interval:
    """ζ(s), s ∈ ℚ, s > 0, s ≠ 1 — сертификат, работающий ПОД единицей.

    Дорога первая. Тождество [внешний источник: DLMF 25.2.8,
    http://dlmf.nist.gov/25.2.E8]:
        ζ(s) = Σ_{k≤N} k^{−s} + N^{1−s}/(s−1) − s·J,
        J = ∫_N^∞ (x − ⌊x⌋)·x^{−s−1} dx,  ℜs > 0, s ≠ 1.
    Ряд Дирихле при s < 1 расходится — потому именно эта дорога, а не
    частная сумма, даёт величину НА КРИТИЧЕСКОЙ ПРЯМОЙ (s = 1/2).

    Обрамление J взято ПЕРВЫМИ ПРИНЦИПАМИ, без внешней оценки остатка:
    на [n, n+1) дробная доля есть u = x − n, множитель (n+u)^{−s−1}
    по u убывает, ∫_0^1 u du = 1/2, значит
        c_n = ∫_0^1 u(n+u)^{−s−1} du ∈ [½(n+1)^{−s−1}, ½ n^{−s−1}],
    и, суммируя, J ∈ [½·T_{N+1}, ½·T_N] при T_M = Σ_{n≥M} n^{−s−1};
    сам T_M обрамлён интегралом (показатель s+1 > 1):
        M^{−s}/s ≤ T_M ≤ M^{−s−1} + M^{−s}/s.
    Ширина итога — порядка s·N^{−s−1}/2, то есть на целый порядок N
    у́же прямого ряда с интегральным хвостом (тот даёт N^{−s}).
    `[вычислено здесь]` замер на s = 1/2: N = 1000 → ширина < 2·10⁻⁵,
    N = 10⁴ → < 6·10⁻⁷.
    """
    s = _slot_doors(s, terms, bits)
    if s <= 0 or s == 1:
        raise Unresolved("enclosure", _rec(s),
                         hint="дорога DLMF 25.2.8 держит ℜs > 0, s ≠ 1; "
                              "отрицательную сторону сводит точный "
                              "множитель functional_factor(s)")
    p, q, N = s.numerator, s.denominator, terms
    scale = 1 << bits
    s_lo, s_hi = _partial_sum(1, N, p, q, bits)
    part = Interval(Q(s_lo, scale), Q(s_hi, scale))
    # N^{1−s}/(s−1): знак делителя разный по сторонам единицы — его
    # разбирает арифметика отрезка, а не отдельная ветка кода
    if q > p:                                   # s < 1: показатель > 0
        a, b = _grid_pow(N, q - p, q, bits)
    else:                                       # s > 1: показатель < 0
        a, b = _grid_pow_neg(N, p - q, q, bits)
    part = part + Interval(Q(a, scale), Q(b, scale)) * (1 / (s - 1))
    # J ∈ [½(N+1)^{−s}/s, ½(N^{−s−1} + N^{−s}/s)]
    j_lo = Q(1, 2) * Q(_grid_pow_neg(N + 1, p, q, bits)[0], scale) / s
    hi_ss = Q(_grid_pow_neg(N, p + q, q, bits)[1], scale)   # N^{−s−1}
    hi_s = Q(_grid_pow_neg(N, p, q, bits)[1], scale)        # N^{−s}
    j_hi = Q(1, 2) * (hi_ss + hi_s / s)
    return _on_grid(Interval(part.lo - s * j_hi, part.hi - s * j_lo), bits)


def zeta_enclosure_dirichlet(s, terms: int = 1000,
                             bits: int = 128) -> Interval:
    """ζ(s), s ∈ ℚ, s > 1 — сертификат прямым рядом с хвостом-интегралом.

    Дорога вторая, НЕЗАВИСИМАЯ от DLMF 25.2.8: Σ_{n≤N} n^{−s} плюс
    строгое обрамление остатка интегральным признаком,
        R ∈ [1/((s−1)(N+1)^{s−1}); 1/((s−1)N^{s−1})].
    Та же теорема, что держит целый слот в `ZetaValue.enclosure`, — но
    при нецелом s члены n^{−p/q} берутся целым корнем на сетке.
    Служит сверкой первой дороги там, где обе применимы (s > 1).
    """
    s = _slot_doors(s, terms, bits)
    if s <= 1:
        raise Unresolved("enclosure", _rec(s),
                         hint="прямой ряд Дирихле сходится лишь при s > 1; "
                              "под единицей дорога — zeta_enclosure_integral")
    p, q, N = s.numerator, s.denominator, terms
    scale = 1 << bits
    s_lo, s_hi = _partial_sum(1, N, p, q, bits)
    hi_np1 = Q(_grid_pow(N + 1, p - q, q, bits)[1], scale)   # (N+1)^{s−1} ↑
    lo_n = Q(max(_grid_pow(N, p - q, q, bits)[0], 1), scale)  # N^{s−1} ↓
    tail_lo = 1 / ((s - 1) * hi_np1)
    tail_hi = 1 / ((s - 1) * lo_n)
    return _on_grid(Interval(Q(s_lo, scale) + tail_lo,
                             Q(s_hi, scale) + tail_hi), bits)


def eta_enclosure(s, terms: int = 1000, bits: int = 128) -> Interval:
    """η(s) = Σ_{n≥1} (−1)^{n−1}·n^{−s}, s ∈ ℚ, s > 0 — сертификат.

    Обрамление — теорема Лейбница, первые принципы: члены n^{−s} строго
    убывают к нулю, потому частные суммы чётного и нечётного номера
    зажимают предел с двух сторон. Ширина — порядка N^{−s}: дорога
    ЧЕСТНАЯ, но дорогая (цена названа в `zeta_enclosure_eta`).
    """
    s = _slot_doors(s, terms, bits)
    if s <= 0:
        raise Unresolved("η", _rec(s),
                         hint="знакочередующийся ряд сходится при s > 0")
    p, q = s.numerator, s.denominator
    scale = 1 << bits
    m = terms if terms % 2 == 0 else terms + 1   # чётный номер — низ
    lo = hi = 0
    for n in range(1, m + 1):
        a, b = _grid_pow_neg(n, p, q, bits)
        if n % 2:
            lo += a
            hi += b
        else:
            lo -= b
            hi -= a
    nxt_hi = _grid_pow_neg(m + 1, p, q, bits)[1]   # член m+1 — со знаком +
    return Interval(Q(lo, scale), Q(hi + nxt_hi, scale))


def zeta_enclosure_eta(s, terms: int = 1000, bits: int = 128) -> Interval:
    """ζ(s) через η: ζ(s) = η(s)/(1 − 2^{1−s}) [внешний источник:
    DLMF 25.2.3, http://dlmf.nist.gov/25.2.E3], s ∈ ℚ, s > 0, s ≠ 1.

    Дорога третья — единственная, что держит 0 < s < 1 НЕЗАВИСИМО от
    DLMF 25.2.8, и потому служит сверкой критической точки. Цена её
    названа числом: ширина падает как N^{−s}, то есть на s = 1/2 каждый
    десятичный разряд стоит СТА крат по N. `[вычислено здесь]` замер:
    N = 10⁴ → ширина ≈ 2,4·10⁻², N = 10⁶ → ≈ 2,4·10⁻³; шесть разрядов
    потребовали бы N ≈ 10¹², то есть дороги нет — есть цена.
    """
    s = _slot_doors(s, terms, bits)
    if s == 1:
        raise Unresolved("ζ(s)", _rec(s), hint="полюс ζ: значения нет")
    eta = eta_enclosure(s, terms, bits)
    p, q = s.numerator, s.denominator
    scale = 1 << bits
    if q > p:                                    # 1−s > 0
        a, b = _grid_pow(2, q - p, q, bits)
    else:                                        # 1−s < 0
        a, b = _grid_pow_neg(2, p - q, q, bits)
    denom = 1 - Interval(Q(a, scale), Q(b, scale))
    return _on_grid(eta / denom, bits)


def _pochhammer(s: Q, j: int) -> Q:
    """(s)_j = s·(s+1)·…·(s+j−1), точно в ℚ; (s)_0 = 1."""
    out = Q(1)
    for i in range(j):
        out *= s + i
    return out


def zeta_enclosure_euler_maclaurin(s, terms: int = 20, bits: int = 128,
                                   bernoulli_terms: int | None = None) -> Interval:
    """ζ(s), s ∈ ℚ, s ≠ 1 — быстрая дорога Эйлера—Маклорена, и ПОД нулём тоже.

    Дорога третья, независимая от двух прежних (DLMF 25.2.8 и η-ряд);
    названа в записке 12.09 «критическая прямая ζ» п. 4 как не построенная
    («удешевила бы десятичную сторону на пять порядков N») и построена
    02.10. Тождество [внешний источник: DLMF 25.2.9,
    http://dlmf.nist.gov/25.2.E9], N ≥ 1, n ≥ 1, ℜs > −2n:
        ζ(s) = Σ_{k=1}^{N−1} k^{−s} + N^{1−s}/(s−1) + N^{−s}/2
               + Σ_{k=1}^{n} B_{2k}/(2k)!·(s)_{2k−1}·N^{−s−2k+1} + R_n,
        R_n = −(s)_{2n+1}/(2n+1)! · ∫_N^∞ B̃_{2n+1}(x)·x^{−s−2n−1} dx.
    Остаток обрамлён ПЕРВЫМ ОТБРОШЕННЫМ ЧЛЕНОМ [внешний источник: DLMF
    §2.10(i), http://dlmf.nist.gov/2.10.i, и Hardy, Divergent Series,
    1949, §13.2]: когда f^{(2n+2)} и f^{(2n+4)} знакопостоянны и одного
    знака на [N, ∞), R_n лежит между 0 и B_{2n+2}/(2n+2)!·(s)_{2n+1}·N^{−s−2n−1};
    здесь берётся симметричное обрамление ±|этого члена| — шире вдвое,
    зато без довода о знаке. У f = x^{−s} это так при s > 0 для всех n
    (производные знакопостоянны и монотонны), а при s < 0 — как только
    2n + 2 > −s: тогда f^{(k)} = (s)_k·x^{−s−k} при k ≥ 2n+2 имеет
    показатель −s−k < 0 и множители (−s−j), j < k, все уже
    отрицательные за j > −s, то есть (s)_{2n+2} и (s)_{2n+4} одного знака;
    потому при s < 0 число членов n не меньше ⌊−s/2⌋ + 2 (02.10 — так
    отрицательная сторона ζ, включая знаменатели вне {1, 2}, получила
    сертификат НАПРЯМУЮ, без точного χ; зеркало через χ — сверка).

    `[вычислено здесь]` Весь счёт в ℚ, кроме одного корня: N^{−s} = N^{−p/q}
    на сетке 2^{−bits} (`_grid_pow_neg`), и каждый член суммы Σ_{k<N} k^{−s}
    — тем же корнем; множители B_{2k}/(2k)!·(s)_{2k−1}·N^{−2k+1} и
    N/(s−1), 1/2 — точные рациональные, общий множитель N^{−s} один на
    всех. Цена против дороги DLMF 25.2.8 (≈1,5 разряда на декаду N): при
    s = 3/2 N = 10, n = 8 дают ширину ≈ 2^−55, N = 20, n = 12 — ≈ 2^−90,
    N = 60, n = 16 — сетка 2^−128 (замер приёмкой). По умолчанию n =
    min(16, max(2, N·2/3)) — больше членов при малом N расходятся
    (ряд асимптотический: |B_{2k}|/(2k)! ~ 2/(2π)^{2k}, а (s)_{2k−1}
    растёт факториально), потому n привязан к N явно.
    """
    s = _slot_doors(s, terms, bits)
    if s == 1:
        raise Unresolved("enclosure", _rec(s), hint="полюс ζ: значения нет")
    n_min = 1 if s > 0 else int((-s).numerator // (-s).denominator) // 2 + 2
    if bernoulli_terms is None:
        n = max(min(16, max(2, terms * 2 // 3)), n_min)
    else:
        if isinstance(bernoulli_terms, bool) or not isinstance(bernoulli_terms, int):
            raise TypeError("bernoulli_terms — только int")
        if bernoulli_terms < n_min:
            raise ValueError(f"bernoulli_terms ≥ {n_min}: при s = {s} остаток обрамлён первым "
                             f"отброшенным членом лишь при 2n + 2 > −s")
        n = bernoulli_terms
    N = terms
    p, q = s.numerator, s.denominator
    # Внутренняя сетка тоньше объявленной на log₂N + 8 бит: N округлений
    # частной суммы иначе сдвигали бы концы на N узлов, и ответ при
    # большем N переставал бы вкладываться в ответ при меньшем (замер
    # 02.10: EM(400) ⊉ EM(800) при равной сетке; с запасом — равны).
    fine = bits + N.bit_length() + 8
    scale = 1 << fine
    # Σ_{k=1}^{N−1} k^{−s} — пара границ на сетке
    if N >= 2:
        s_lo, s_hi = _partial_sum(1, N - 1, p, q, fine)
    else:
        s_lo = s_hi = 0
    part = Interval(Q(s_lo, scale), Q(s_hi, scale))
    r_lo, r_hi = _grid_pow_signed(N, p, q, fine)
    r = Interval(Q(r_lo, scale), Q(r_hi, scale))          # N^{−s}, p любого знака
    coef = Q(N) / (s - 1) + Q(1, 2)                       # N^{1−s}/(s−1) + N^{−s}/2 = N^{−s}·(N/(s−1) + 1/2)
    for k in range(1, n + 1):
        coef += bernoulli(2 * k) / _factorial(2 * k) * _pochhammer(s, 2 * k - 1) / Q(N) ** (2 * k - 1)
    rem = abs(bernoulli(2 * n + 2) / _factorial(2 * n + 2) * _pochhammer(s, 2 * n + 1) / Q(N) ** (2 * n + 1))
    main = r * Interval(coef - rem, coef + rem)
    return _on_grid(part + main, bits)


def _factor_interval(chi, bits: int) -> Interval:
    """Отрезок точного множителя χ: Q — точкой, моном q·π^e — через
    сертификат π (показатель e у χ всегда ЦЕЛЫЙ: e = m при полуцелом s
    и e = s−1 при целом).

    `[вычислено здесь]` порядок ряда Мэчина берётся ОТ РАЗРЯДНОСТИ, а не
    от N: член atan(1/5) весит ≈ log₂25 ≈ 4,6 бита, потому bits//4 + 12
    членов покрывают сетку с запасом. Прежде сюда шло само N, и при
    N = 2000 знаменатели концов доходили до 4300+ десятичных знаков —
    счёт верный, но запись негодная.
    """
    if isinstance(chi, Q):
        return Interval.point(chi)
    if chi.e.denominator != 1:
        raise Unresolved("χ", chi.record(),
                         hint="дробный показатель π у множителя не ожидался")
    pi_iv = pi_enclosure(max(12, bits // 4 + 12))
    return Interval.point(chi.coeff) * (pi_iv ** chi.e.numerator)


def _on_grid(iv: Interval, bits: int) -> Interval:
    """Концы — наружу к объявленной сетке 2^−bits (канон enclosure.py:
    «концы сертификатов живут на сетке»). Вложение сохраняется, цена —
    не более узла на конец; знаменатели перестают расти произведениями.
    """
    return Interval.on_grid(iv.lo, iv.hi, bits)


# ---------- самопроверка линейки (ось самоприменения, 08.10) ----------

def self_check(deep: bool = False) -> list:
    """Самопроверка линейки класса 13 встречными дорогами; список находок,
    пусто — здорова; исключений не бросает.

    `[вычислено здесь]` Окно носителей: s ∈ {3, 5, 1/2, 3/2, −1/2, 7/3}
    (deep — ещё 7, 5/4, −3/2) проверяют себя (`ZetaValue.self_check`);
    точные точки фабрики против внешних фактов: ζ(0) = −1/2, ζ(−1) = −1/12,
    тривиальные нули ζ(−2) = ζ(−4) = 0, ζ(2) = π²/6, ζ(4) = π⁴/90 [внешний
    источник: формула Эйлера]; числа Бернулли B₂ = 1/6, B₄ = −1/30, B₃ = 0;
    функциональное уравнение на целой полосе: χ(−1)·ζ(2) = ζ(−1) и
    χ(−3)·ζ(4) = ζ(−3) точно — вторая дорога к тривиальной стороне, минуя
    рекурренту Бернулли в ζ(−n). Граница — окно."""
    out = []
    pts = [Q(3), Q(5), Q(1, 2), Q(3, 2), Q(-1, 2), Q(7, 3)]
    if deep:
        pts += [Q(7), Q(5, 4), Q(-3, 2)]
    for s in pts:
        try:
            v = zeta(s)
            if isinstance(v, ZetaValue):
                out += [f"ζ({s}): {f}" for f in v.self_check(deep)]
            else:
                out.append(f"ζ({s}) приземлилось в {v!r} без носителя")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"ζ({s}): {type(e).__name__}: {e}")
    try:
        if zeta(0) != Q(-1, 2):
            out.append(f"ζ(0) = {zeta(0)!r} ≠ −1/2")
        if zeta(-1) != Q(-1, 12):
            out.append(f"ζ(−1) = {zeta(-1)!r} ≠ −1/12")
        if not (_exact_zero(zeta(-2)) and _exact_zero(zeta(-4))):
            out.append("тривиальные нули ζ(−2), ζ(−4) не точный ноль")
        if zeta(2) != PiMonomial(Q(1, 6), 2) or zeta(4) != PiMonomial(Q(1, 90), 4):
            out.append(f"приземление Эйлера: ζ(2) = {zeta(2)!r}, ζ(4) = {zeta(4)!r}")
        if bernoulli(2) != Q(1, 6) or bernoulli(4) != Q(-1, 30) or bernoulli(3) != 0:
            out.append(f"числа Бернулли: B₂ = {bernoulli(2)}, B₃ = {bernoulli(3)}, B₄ = {bernoulli(4)}")
        for m in (1, 3):
            lhs = functional_factor(-m) * zeta(m + 1)
            if lhs != zeta(-m):
                out.append(f"функциональное уравнение при s = −{m}: χ·ζ({m + 1}) = {lhs!r} ≠ ζ(−{m}) = {zeta(-m)!r}")
    except Exception as e:                                      # noqa: BLE001
        out.append(f"точные точки: {type(e).__name__}: {e}")
    return out
