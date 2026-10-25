"""Линейка класса 9 INTEGRAL: x = ∫_a^b f(t) dt, f = num/den рациональна.

Величина хранится структурно: слоты num, den ∈ ℤ[t] (IntPoly) и пределы
a, b ∈ ℚ. Направленные правила фабрики integral() — двойственность
«первообразная ↔ 5/6/10»: сводимая первообразная приземляется в чужую
линейку, не считается заново:

  * den — константа: точный Q почленной первообразной
    Σ cᵢ·(b^{i+1} − a^{i+1})/(i+1) (класс 1);
  * f = c/t при a > 0, b > 0: ∫ c/t = c·ln(b/a) — класс 6 (LogValue);
  * f = c/(1+t²): c·(arctan b − arctan a) по табличным точным точкам
    arctan {0 → 0, 1 → π/4, −1 → −π/4} — класс 10 (PiMonomial);
    иные пределы этой таблицей не сведены;
  * иначе — несводимое IntegralValue со слотами (num, den, a, b);
    корень den на [a; b] — честный отказ: особая точка, значения нет.

Устранимая особенность отделена от настоящей (правка 12.09, вердикт
главы 07 §19 п. 11): num и den сокращаются на НОД в ℤ[t] ДО проверки
Штурмом, потому общий корень (0/0 в точке меры нуль) значения не
отменяет — `∫₀¹ (2t−1)/((2t−1)(t+1))` сводится к ∫ dt/(t+1). Отказывает
только корень den, которого у num нет.

Канон держит фабрика: a = b — точный ноль, a > b — переворот знаком
числителя (∫_a^b = ∫_b^a(−f)), так что в слотах всегда a < b.

Структурные правила: линейность над общим отрезком (общий den — сумма
числителей, разные — приведение к den·den′) и аддитивность по смежным
отрезкам одной подынтегральной (∫_a^b + ∫_b^c = ∫_a^c); разность равных —
точный ноль. Прочее — с 05.09 замыканием D-конечных рядов (dfinite.py):
интеграл рациональной f есть D-конечный ряд (разложение первообразной
вокруг середины куска, сходимость на концах — точной оценкой |D₀| >
Σ|Dᵢ|hⁱ, отрезок делится до сходимости), поэтому ∫ + ∫ разных слотов,
∫ ± число, ∫ · ∫, ∫ · ряд, ∫ⁿ — ряд класса 7 с точной рекуррентой и
сертификатом из строения; линейки без D-конечного ряда (ln 2 класса 6)
получают отказ поимённо; обращение и деление — отказ (замыкания нет).
Десятичная сторона — только сертификатом сетки (enclosure) по явному
запросу: обрамление f на каждом из N подотрезков Interval-арифметикой,
суммы нижних/верхних концов дают строгое вложение. Ни одного float:
float на входе — TypeError построением.
"""
from __future__ import annotations

from .rational import Q

from .algebraic import _divexact, _int_poly_gcd, count_roots, squarefree_part
from .enclosure import Interval
from .errors import Unresolved
from .exp_log import LogValue
from .intpoly import IntPoly
from .pi_value import PiMonomial
from .zero_state import Zero

#: Допустимые скаляры. float исключён построением.
_NUM = (int, Q)
#: Одна запись отказа на bool во всех дверях модуля: bool числом в доме не
#: считается и называется bool'ом (канон `phicore/rational.py`; сосед по
#: главе отвечает тем же: `Notch.from_q(True, 8)` → «bool не число»).
#: Двух копий одного канона здесь нет намеренно: разошедшиеся копии — уже
#: пройденная поломка (`affine.py:75`).
_NOT_BOOL = "bool не число"

#: Табличные точные точки arctan: только {0, 1, −1} → доли π
#: (arctan 1 = π/4 — определение четверти круга, не разложение).
_ATAN_TABLE = {
    Q(0): Q(0),
    Q(1): Q(1, 4),
    Q(-1): Q(-1, 4),
}


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
        raise TypeError(f"{name} — {_NOT_BOOL}; дверь принимает int/Q")
    if not isinstance(x, _NUM):
        raise TypeError(f"{name} — только int/Q, float запрещён")
    return Q(x)


def _slots_record(num: IntPoly, den: IntPoly, a: Q, b: Q) -> str:
    """Явная Ф-запись слотов линейки — одна форма для record() и отказов."""
    return (f"Φ(C=9 INTEGRAL, f=[{','.join(map(str, num.c))}]"
            f"/[{','.join(map(str, den.c))}], пределы=[{a}; {b}])")


def _reduce_common(num: IntPoly, den: IntPoly) -> "tuple[IntPoly, IntPoly]":
    """Сокращение num/den на НОД в ℤ[t] — снятие УСТРАНИМОЙ особенности.

    Общий корень num и den есть 0/0: подынтегральная в этой ОДНОЙ точке
    не определена, но величина интеграла от точки меры нуль не зависит и
    равна интегралу сокращённой дроби, непрерывной на всём отрезке.
    Настоящая особенность — корень den, которого у num нет, — НОДом не
    снимается: там значения нет, и отказ остаётся отказом.

    Ход целиком целочисленный: НОД — примитивная PRS
    (`algebraic._int_poly_gcd`), деление — точное (`_divexact`); по лемме
    Гаусса примитивный делитель делит num и den в ℤ[t] без остатка.
    НОД примитивен, потому константный НОД всегда есть [1] — общий
    целый множитель коэффициентов (content) сокращением не считается и
    слоты не трогает. Точность деления проверяется умножением обратно:
    структурная невозможность должна падать громко, не молча.

    Замер 12.09 (вердикт главы 07 §19 п. 11, строка 75 очереди):
    `∫₀¹ (2t−1)/((2t−1)(t+1))` до сокращения отказывал «особая точка» при
    существующем значении ln 2; после — den сокращается до (t+1), корней
    на [0; 1] нет. Свод 394 пар (num·g)/(den·g) со случайными
    коэффициентами: сокращений 378, нарушений тождества величины
    qn·D = qd·N — НИ ОДНОГО.

    ГРАНИЦА, названная 12.09, ЗАКРЫТА 17.09 ТОЧНОЙ ДОРОГОЙ. Прежде пара
    (0, den) возвращалась как есть, и `∫₀¹ 0/t` отказывал «особая точка»
    при существующем значении. Разбор показал, что это был ноль
    ИНСТРУМЕНТА, а не предмета: НОД(0, den) = den в ℤ[t] по определению
    (всякий делитель den делит и нуль, и den есть наибольший общий из
    них), а «не определён» был лишь примитивный PRS `_int_poly_gcd`,
    которому нужен ненулевой вход. Тот же корень дом уже закрывал в
    контейнере (`gcd(0, d) = d`, сторож в `_reduce_ratio`) `[в корпусе]`.
    Потому num ≡ 0 сокращается на весь den: 0/den → 0/1. Довод — тот же
    довод устранимости, что и выше: при num ≡ 0 ВСЯКИЙ корень den есть
    ОБЩИЙ корень пары, то есть 0/0 в точке меры нуль, а не полюс. На
    [a; b] минус конечное множество корней den подынтегральная тождественно
    нуль, и величина интеграла есть точный нуль — обеими дорогами:
    по Лебегу (равна нулю почти всюду) и несобственной по Риману
    (lim_{ε→0} ∫_{a+ε}^{b} 0 dt = 0). Цена дороги названа: den ≡ 0
    сокращением НЕ трогается — там 0/0 не в точке, а на всём отрезке,
    значения нет ни по одной дороге, и отказ фабрики остаётся отказом.
    """
    if den.is_zero():
        return num, den
    if num.is_zero():
        # НОД(0, den) = den: сокращение на весь den, 0/den → 0/1.
        # Ниже по ходу это приземляет ∫ в класс 1 точным нулём вместо
        # отказа «особая точка» (замер 17.09: `∫₀¹ 0/t` → Q(0)).
        return num, IntPoly([1])
    g = _int_poly_gcd(num.c, den.c)
    if len(g) == 1:  # НОД — константа: общего корня нет, сокращать нечего
        return num, den
    gp = IntPoly(g)
    qn, qd = IntPoly(_divexact(num.c, g)), IntPoly(_divexact(den.c, g))
    if qn * gp != num or qd * gp != den:
        raise ValueError("сокращение на НОД не точно — нарушена лемма "
                         f"Гаусса на num={num.c}, den={den.c}, НОД={g}")
    return qn, qd


def _canon_content(num: IntPoly, den: IntPoly) -> "tuple[IntPoly, IntPoly]":
    """Канон содержания слотов: общий целый множитель коэффициентов num и
    den сокращается, знак старшего коэффициента den — положителен.

    До 01.10 (наблюдение сдачи 25.09, очередь строка 104; приёмка D
    двойственности: слоты [9]/[9,0,0,9] против [1]/[1,0,0,1]) равные
    подынтегральные с разными слотами не были равны `==` — запись несла
    не величину, а историю приведения к общему знаменателю. Величина
    интеграла от общего множителя num/den не зависит, и знак den уходит
    в num без потери: сокращение — целыми, точность проверяется в
    `_reduce_common` тем же приёмом (умножением обратно)."""
    from math import gcd
    g = gcd(num.content(), den.content())
    if g > 1:
        num = IntPoly([c // g for c in num.c])
        den = IntPoly([c // g for c in den.c])
    if den.c[-1] < 0:
        num, den = -num, -den
    return num, den


def _den_root_on(num: IntPoly, den: IntPoly, a: Q, b: Q) -> bool:
    """Есть ли НАСТОЯЩАЯ особая точка на замкнутом [a; b].

    Сокращение на НОД(num, den) идёт ДО проверки Штурмом (правка 12.09):
    общий корень num и den есть устранимая особенность, значение
    существует — отказывать там нечем. Далее — Штурм по свободной от
    квадратов части СОКРАЩЁННОГО den: count_roots считает на (a; b],
    левая граница проверяется значением отдельно. Сокращённый den —
    константа (≠ 0 построением): корней нет вовсе, Штурм не зовётся
    (цепочка константы не строится — производная тождественно нуль).
    """
    _, den = _reduce_common(num, den)
    if len(den.c) == 1:
        return False
    if den.eval(a) == 0:
        return True
    return count_roots(squarefree_part(den), a, b) > 0


def _poly_on(p: IntPoly, x: Interval) -> Interval:
    """Вложение множества значений p на отрезке x: схема Горнера в
    Interval-арифметике (монотонна по включению)."""
    acc = Interval.point(p.c[-1])
    for c in reversed(p.c[:-1]):
        acc = acc * x + c
    return acc


class IntegralValue:
    """Несводимое ∫_a^b num/den: num ≢ 0, den — не константа, a < b.

    Канон (сводимые первообразные приземлены в классы 1/6/10, пределы
    упорядочены) держит фабрика integral(); конструктор отбрасывает лишь
    вырожденные слоты и указывает на фабрику.
    """

    __slots__ = ("num", "den", "a", "b")

    def __init__(self, num: IntPoly, den: IntPoly, a, b):
        if not isinstance(num, IntPoly) or not isinstance(den, IntPoly):
            raise TypeError("num/den — только IntPoly")
        a = _as_fraction(a, "a")
        b = _as_fraction(b, "b")
        if den.is_zero():
            raise ValueError("den ≡ 0 не есть подынтегральная величина")
        if len(den.c) == 1:
            raise ValueError("den-константа сводится точной первообразной — "
                             "стройте через integral()")
        if num.is_zero():
            raise ValueError("num ≡ 0 — точный ноль, не интеграл; "
                             "стройте через integral()")
        if a == b:
            raise ValueError("a=b — вырожденный отрезок, точный ноль; "
                             "стройте через integral()")
        if a > b:
            raise ValueError("ориентация a < b держится фабрикой integral()")
        self.num, self.den, self.a, self.b = num, den, a, b

    # --- запись ---------------------------------------------------------

    def record(self) -> str:
        return _slots_record(self.num, self.den, self.a, self.b)

    __repr__ = record

    def __eq__(self, other):
        if isinstance(other, IntegralValue):
            return (self.num == other.num and self.den == other.den
                    and self.a == other.a and self.b == other.b)
        # равенство несводимого интеграла с другими линейками этой
        # записью не решается — ответ не выдумывается
        return NotImplemented

    def __hash__(self):
        return hash((IntegralValue, tuple(self.num.c), tuple(self.den.c),
                     self.a, self.b))

    # --- сложение / вычитание: только одноимённые слоты -----------------

    def _addsub(self, other, op: str, sign: int):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__radd__/__rsub__
        if isinstance(other, IntegralValue):
            onum = other.num if sign > 0 else -other.num
            if other.a == self.a and other.b == self.b:
                # линейность над общим отрезком: общий den — сумма
                # числителей; разные den — приведение к общему den·den′
                # (корней на отрезке нет ни у одного, значит и у
                # произведения); фабрика доведёт вырождения (взаимное
                # уничтожение — точный ноль) и приземления
                if other.den == self.den:
                    return integral(self.num + onum, self.den, self.a, self.b)
                return integral(self.num * other.den + onum * self.den,
                                self.den * other.den, self.a, self.b)
            if sign > 0 and other.num == self.num and other.den == self.den:
                # аддитивность по отрезку: ∫_a^b f + ∫_b^c f = ∫_a^c f
                if self.b == other.a:
                    return integral(self.num, self.den, self.a, other.b)
                if other.b == self.a:
                    return integral(self.num, self.den, other.a, self.b)
            # иная пара — замыкание D-конечных рядов (05.09): оба интеграла —
            # D-конечные ряды с точными членами, сумма — одной линейкой
            return self._dfinite(other, "+" if sign > 0 else "-")
        if isinstance(other, _NUM) and not isinstance(other, bool):
            return self._dfinite(other, "+" if sign > 0 else "-")
        from .dfinite import DFiniteSeries
        if DFiniteSeries._lift(other) is not None:
            return self._dfinite(other, "+" if sign > 0 else "-")
        if callable(getattr(other, "record", None)):
            raise Unresolved(op, self.record(), _rec(other),
                             hint="интеграл с этой линейкой не сведён: у неё нет D-конечного ряда")
        return NotImplemented

    def _dfinite(self, other, op: str):
        """Действие через D-конечные ряды: ∫ → ряд (dfinite.from_integral),
        другой операнд поднимается тем же мостом; результат — ряд класса 7
        с сертификатом из строения (сертификат интеграла — обрамление Римана)."""
        from .dfinite import DFiniteSeries
        a = DFiniteSeries.from_integral(self)
        b = DFiniteSeries._lift(other)
        if op == "+":
            return a + b
        if op == "-":
            return a - b
        return a * b

    def __add__(self, other):
        return self._addsub(other, "+", +1)

    def __sub__(self, other):
        return self._addsub(other, "-", -1)

    def __radd__(self, other):
        if isinstance(other, _NUM) and not isinstance(other, bool):
            from .dfinite import DFiniteSeries
            return DFiniteSeries.const(other) + DFiniteSeries.from_integral(self)
        return NotImplemented

    def __rsub__(self, other):
        if isinstance(other, _NUM) and not isinstance(other, bool):
            from .dfinite import DFiniteSeries
            return DFiniteSeries.const(other) - DFiniteSeries.from_integral(self)
        return NotImplemented

    def __neg__(self):
        # −∫ f = ∫ (−f): знак уходит в числитель, слоты остаются каноном
        return integral(-self.num, self.den, self.a, self.b)

    # --- умножение / деление: сведён только рациональный масштаб --------

    def __mul__(self, other):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__rmul__
        if isinstance(other, _NUM):
            q = Q(other)
            if q == 0:
                return Q(0)
            # q·(num/den) = (q_num·num)/(q_den·den) — точно, без дробных
            # коэффициентов; корни den не меняются
            return integral(self.num.scale(q.numerator),
                            self.den.scale(q.denominator), self.a, self.b)
        from .dfinite import DFiniteSeries
        if DFiniteSeries._lift(other) is not None:
            # 05.09: произведение интегралов и произведение с рядом — свёртка
            # Коши D-конечных рядов (слово автора 03.09: «отказ не предполагал
            # других методов»)
            return self._dfinite(other, "*")
        if callable(getattr(other, "record", None)):
            raise Unresolved("×", self.record(), _rec(other),
                             hint="произведение с этой линейкой не сведено: у неё нет D-конечного ряда")
        return NotImplemented

    __rmul__ = __mul__

    def __truediv__(self, other):
        if isinstance(other, _NUM):
            q = Q(other)
            if q == 0:
                raise Unresolved("/", self.record(), _rec(other),
                                 hint="x/0 здесь не сведено")
            return self * (1 / q)
        if callable(getattr(other, "record", None)):
            raise Unresolved("/", self.record(), _rec(other),
                             hint="частное с другими линейками не сведено")
        return NotImplemented

    def __rtruediv__(self, other):
        if isinstance(other, _NUM):
            q = Q(other)
            if q == 0:
                return Q(0)
            raise Unresolved("/", _rec(other), self.record(),
                             hint="обращение интеграла не сведено")
        return NotImplemented

    def __pow__(self, n: int):
        """Целая степень: 0 — точная единица, 1 — сам носитель, ≥ 2 —
        ряд класса 7 (свёртка Коши), отрицательная — отказ.

        КАНОН НУЛЕВОЙ СТЕПЕНИ — один на дом (замерен 12.09, строка 75
        очереди; прежде здесь стояла 1, а класс 18 отказывал тем же
        доводом — один вопрос с двумя ответами). Нулевая степень есть
        ПУСТОЕ ПРОИЗВЕДЕНИЕ: факт о строении операции, не о величине
        основания, потому она равна 1 у всякого носителя, который несёт
        ВЕЛИЧИНУ, — включая носитель, чья величина есть точный ноль
        (`Q(0) ** 0` → 1, оговорено прямо в `rational.py:196—197`;
        `Hyper.const(0) ** 0` → 1; `ArcTrigValue ** 0` → 1). Здесь это
        важно поимённо: `∫₋₁¹ t/(1+t²)` есть ТОЧНЫЙ ноль (куски [−1; 0] и
        [0; 1] зеркальны, сертификат симметричен), и ответ всё равно 1.

        Отказ канон держит не за нулевой величиной, а за носителем,
        который величины не несёт вовсе: класс 0 (пустота) и класс 20
        (не вычислено) возвращают сам носитель — `container.py:1028`
        «пустота не становится единицей даже при n=0». Класс 9 величину
        несёт, потому единица.
        """
        if not isinstance(n, int) or isinstance(n, bool):
            return NotImplemented
        if n == 0:
            return Q(1)
        if n == 1:
            return self
        if n >= 2:
            from .dfinite import DFiniteSeries
            return DFiniteSeries.from_integral(self) ** n
        raise Unresolved("^", self.record(), _rec(n),
                         hint="отрицательная степень интеграла — обращение ряда, замыкания нет")

    # --- десятичная сторона: сертификат сетки ---------------------------

    def enclosure(self, n: int = 64) -> Interval:
        """Сертификат сетки [lo, hi] ⊂ ℚ обрамлением f на N подотрезках.

        Отрезок [a; b] делится на N равных частей; на каждой множество
        значений f обрамляется интервальной оценкой f([lo, hi]) =
        num([lo, hi])/den([lo, hi]) в Interval-арифметике; суммы нижних
        и верхних концов, умноженные на шаг, дают строгое вложение
        интеграла. Уже́ние — явным запросом (рост N). Корень den на
        [a; b] (Штурм по свободной от квадратов части) — честный отказ:
        особая точка.
        """
        if isinstance(n, bool):
            raise TypeError(f"N — {_NOT_BOOL}; дверь принимает int")
        if not isinstance(n, int):
            raise TypeError("N — только int")
        if n < 1:
            raise ValueError("сетка требует N ≥ 1")
        num, den = _reduce_common(self.num, self.den)
        if _den_root_on(num, den, self.a, self.b):
            raise Unresolved("enclosure", self.record(),
                             hint="особая точка: den имеет корень на "
                                  "отрезке [a; b]")
        h = (self.b - self.a) / n
        lo_sum = hi_sum = Q(0)
        for k in range(n):
            piece = Interval(self.a + k * h, self.a + (k + 1) * h)
            num_iv = _poly_on(num, piece)
            den_iv = _poly_on(den, piece)
            try:
                f_iv = num_iv * den_iv.inverse()
            except Unresolved:
                # корней den на [a; b] нет (Штурм выше) — ноль в
                # den([lo, hi]) есть переоценка вложения на грубой сетке
                raise Unresolved(
                    "enclosure", self.record(),
                    hint=f"сетка N={n} груба: den([lo, hi]) несёт ноль "
                         "переоценкой вложения; увеличить N явным "
                         "запросом") from None
            lo_sum += f_iv.lo
            hi_sum += f_iv.hi
        return Interval(lo_sum * h, hi_sum * h)


    # --- самопроверка носителя (ось самоприменения, 08.10) ----------------

    def self_check(self, deep: bool = False) -> list:
        """Самопроверка носителя встречными дорогами (ось «AIAΣA спрашивает
        AIAΣA», слово автора 04.09; очередь строка 18; 08.10). Список
        находок, пусто — здоров; исключений не бросает.

        `[вычислено здесь]` (1) слоты: num, den — IntPoly, num ≢ 0, den не
        константа, a < b в ℚ; канон — num и den взаимно просты в ℤ[t] и
        содержание с знаком канонические (иначе фабрика сократила бы);
        (2) сертификат enclosure(64): концы упорядочены, и сетка вдвое
        мельче вложена в него (вложенность интервальной арифметики — теорема,
        не договор); (3) ВСТРЕЧНАЯ дорога — ряд класса 7 (`dfinite.from_integral`:
        разложение первообразной вокруг середины куска, члены точны): частная
        сумма S_40 обязана лежать в сертификате с запасом |S_40 − S_20| —
        запас есть объявленная мера хвоста, не его сертификат; дорога
        независима от обрамления Римана по построению. При deep — S_80 и
        запас |S_80 − S_40|, и сетка 256."""
        out = []
        num, den, a, b = self.num, self.den, self.a, self.b
        if not isinstance(num, IntPoly) or not isinstance(den, IntPoly):
            return [f"слоты num/den чужого носителя: {type(num).__name__}, {type(den).__name__}"]
        if not isinstance(a, Q) or not isinstance(b, Q):
            return [f"пределы чужого носителя: {type(a).__name__}, {type(b).__name__}"]
        if num.is_zero() or den.is_zero() or len(den.c) == 1:
            return ["вырожденные слоты: num ≡ 0, den ≡ 0 либо den — константа (сведение фабрикой integral)"]
        if a >= b:
            return [f"пределы не упорядочены: a = {a} ≥ b = {b} (канон фабрики — a < b)"]
        try:
            n2, d2 = _canon_content(*_reduce_common(num, den))
            if n2 != num or d2 != den:
                out.append(f"слоты вне канона: сокращение даёт [{','.join(map(str, n2.c))}]/[{','.join(map(str, d2.c))}]")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"канон не проверен: {type(e).__name__}: {e}")
        try:
            iv = self.enclosure(64)
        except Exception as e:                                  # noqa: BLE001
            return out + [f"сертификат не построен: {type(e).__name__}: {e}"]
        if iv.lo > iv.hi:
            out.append(f"концы сертификата переставлены: {iv.record()}")
        try:
            fine = self.enclosure(256 if deep else 128)
            if fine.lo < iv.lo or fine.hi > iv.hi:
                out.append(f"мелкая сетка не вложена в грубую: {fine.record()} ⊄ {iv.record()}")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"мелкая сетка не построена: {type(e).__name__}: {e}")
        try:
            from .dfinite import DFiniteSeries
            s = DFiniteSeries.from_integral(self)
            n = 80 if deep else 40
            s_n, s_h = s.partial_sum(n), s.partial_sum(n // 2)
            gap = abs(s_n - s_h)
            if s_n < iv.lo - gap or s_n > iv.hi + gap:
                out.append(f"ряд класса 7: S_{n} = {s_n} вне сертификата {iv.record()} с запасом {gap}")
        except Unresolved as e:
            out.append(f"ряд класса 7 отказал: {e}")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"ряд класса 7 не построен: {type(e).__name__}: {e}")
        return out


def integral(num: IntPoly, den: IntPoly, a, b) \
        -> Q | LogValue | PiMonomial | IntegralValue:
    """∫_a^b (num/den) dt: сведение направленными правилами.

    den-константа — точный Q; c/t при a > 0, b > 0 — класс 6;
    c/(1+t²) по табличным пределам {0, 1, −1} — класс 10; прочее —
    IntegralValue; корень den на [a; b] — отказ (особая точка).
    """
    if not isinstance(num, IntPoly) or not isinstance(den, IntPoly):
        raise TypeError("num/den — только IntPoly")
    a = _as_fraction(a, "a")
    b = _as_fraction(b, "b")
    if den.is_zero():
        raise Unresolved("∫", _slots_record(num, den, a, b),
                         hint="den ≡ 0: подынтегральная величина не "
                              "определена ни в одной точке")
    # УСТРАНИМАЯ особенность снимается ДО всякого разбора и до Штурма:
    # общий множитель num и den даёт 0/0 в точках меры нуль, величина от
    # него не зависит. Слоты канона — сокращённые; потому сокращение
    # способно приземлить интеграл в класс 1 (den стал константой) — то
    # же правило 1, которому оно предшествует.
    num, den = _reduce_common(num, den)
    num, den = _canon_content(num, den)      # канон содержания и знака (01.10)
    if a == b:
        if den.eval(a) == 0:
            raise Unresolved("∫", _slots_record(num, den, a, b),
                             hint="особая точка: den имеет корень на "
                                  "отрезке [a; b]")
        return Q(0)  # вырожденный отрезок — точный ноль класса 1
    if a > b:
        # ∫_a^b f = ∫_b^a (−f): знак уходит в числитель, пределы — в канон
        return integral(-num, den, b, a)

    # --- правило 1: den — константа, точная почленная первообразная ---
    if len(den.c) == 1:
        total = Q(0)
        for i, c in enumerate(num.c):
            if c:
                total += Q(c, i + 1) * (b ** (i + 1) - a ** (i + 1))
        return total / den.c[0]

    # den не константа: сперва честная проверка особых точек (Штурм по
    # сокращённому den — устранимая особенность отказом уже не считается)
    if _den_root_on(num, den, a, b):
        raise Unresolved("∫", _slots_record(num, den, a, b),
                         hint="особая точка: den имеет корень на "
                              "отрезке [a; b]")
    if num.is_zero():
        return Q(0)  # f ≡ 0 на всём [a; b] — точный ноль

    # --- правило 2: f = c/t, 0 ∉ [a; b] — приземление в класс 6 ---
    if len(num.c) == 1 and len(den.c) == 2 and den.c[0] == 0:
        # Штурм выше уже исключил 0 из [a; b]: либо 0 < a < b, либо a < b < 0.
        # В обоих случаях ∫_a^b c/t dt = c·ln(b/a), b/a > 0 (сторона a < 0
        # приземляется с 01.10 — наблюдение сдачи 25.09, очередь строка 104;
        # прежде отрицательная сторона оставалась носителем класса 9).
        return LogValue.make(Q(num.c[0], den.c[1]), b / a)

    # --- правило 3: f = c/(k+k·t²) — приземление в класс 10 ---
    if len(num.c) == 1 and len(den.c) == 3 and den.c[1] == 0 \
            and den.c[0] == den.c[2]:
        if a in _ATAN_TABLE and b in _ATAN_TABLE:
            coeff = (Q(num.c[0], den.c[2])
                     * (_ATAN_TABLE[b] - _ATAN_TABLE[a]))
            if coeff == 0:  # страховка вырождения — ноль живёт в классе 1
                return Q(0)
            return PiMonomial(coeff, 1)
        # пределы вне таблицы arctan — несводимое, вниз к правилу 4

    # --- правило 4: несводимый слот линейки ---
    return IntegralValue(num, den, a, b)


# ---------- самопроверка линейки (ось самоприменения, 08.10) ----------

def self_check(deep: bool = False) -> list:
    """Самопроверка линейки класса 9 встречными дорогами; список находок,
    пусто — здорова; исключений не бросает.

    `[вычислено здесь]` Окно носителей: ∫₀¹ dt/(1+t³), ∫₀¹ t dt/(1+t⁴),
    ∫₁² dt/(t²+1) (deep — ещё ∫_{−1}^{1} dt/(2+t³) и ∫₀³ dt/(1+t²)) проверяют
    себя (`IntegralValue.self_check`); дорога через класс 6: ∫₀¹ 2t/(1+t²) —
    носитель, сертификат обязан нести ln 2 рядом логарифма (`ln_enclosure`)
    [внешний источник: первообразная ln(1+t²)]; приземления фабрики: den-
    константа — точный Q, ∫₁² dt/t = ln 2 классом 6, ∫₀¹ dt/(1+t²) = π/4
    классом 10, a = b — ноль, a > b — переворот знаком, устранимая
    особенность (2t−1)/((2t−1)(t+1)) сокращается до dt/(t+1), полюс на
    отрезке — отказ; аддитивность ∫₀¹ + ∫₁² = ∫₀² точно. Граница — окно."""
    out = []
    window = [(IntPoly([1]), IntPoly([1, 0, 0, 1]), 0, 1), (IntPoly([0, 1]), IntPoly([1, 0, 0, 0, 1]), 0, 1),
              (IntPoly([1]), IntPoly([1, 0, 1]), 1, 2)]
    if deep:
        window += [(IntPoly([1]), IntPoly([2, 0, 0, 1]), -1, 1), (IntPoly([1]), IntPoly([1, 0, 1]), 0, 3)]
    for num, den, a, b in window:
        try:
            v = integral(num, den, a, b)
            if isinstance(v, IntegralValue):
                out += [f"{v.record()}: {f}" for f in v.self_check(deep)]
            else:
                out.append(f"∫ {_slots_record(num, den, Q(a), Q(b))} приземлился в {v!r} без носителя")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"∫ {_slots_record(num, den, Q(a), Q(b))}: {type(e).__name__}: {e}")
    try:
        from .enclosure import ln_enclosure
        v = integral(IntPoly([0, 2]), IntPoly([1, 0, 1]), 0, 1)
        if not isinstance(v, IntegralValue):
            out.append(f"∫₀¹ 2t/(1+t²) приземлился в {v!r} без носителя")
        elif v.enclosure(128).intersect(ln_enclosure(Q(2), 40)) is None:
            out.append(f"дорога класса 6: ∫₀¹ 2t/(1+t²) = {v.enclosure(128).record()} не несёт ln 2")
    except Exception as e:                                      # noqa: BLE001
        out.append(f"дорога класса 6: {type(e).__name__}: {e}")
    try:
        if integral(IntPoly([0, 0, 1]), IntPoly([1]), 0, 1) != Q(1, 3):
            out.append("den-константа: ∫₀¹ t² dt ≠ 1/3")
        if integral(IntPoly([1]), IntPoly([0, 1]), 1, 2) != LogValue.make(Q(1), Q(2)):
            out.append(f"класс 6: ∫₁² dt/t = {integral(IntPoly([1]), IntPoly([0, 1]), 1, 2)!r} ≠ ln 2")
        if integral(IntPoly([1]), IntPoly([1, 0, 1]), 0, 1) != PiMonomial(Q(1, 4), 1):
            out.append(f"класс 10: ∫₀¹ dt/(1+t²) = {integral(IntPoly([1]), IntPoly([1, 0, 1]), 0, 1)!r} ≠ π/4")
        if integral(IntPoly([1]), IntPoly([1, 0, 0, 1]), Q(1, 2), Q(1, 2)) != Q(0):
            out.append("a = b: не точный ноль")
        flip = integral(IntPoly([1]), IntPoly([1, 0, 0, 1]), 1, 0)
        if not isinstance(flip, IntegralValue) or flip.num != IntPoly([-1]) or (flip.a, flip.b) != (Q(0), Q(1)):
            out.append(f"a > b: переворот не канонический: {flip!r}")
        if integral(IntPoly([-1, 2]), IntPoly([-1, 1, 2]), 0, 1) != integral(IntPoly([1]), IntPoly([1, 1]), 0, 1):
            out.append("устранимая особенность (2t−1)/((2t−1)(t+1)) не сокращена до dt/(t+1)")
        try:
            integral(IntPoly([1]), IntPoly([-1, 0, 1]), 0, 2)
            out.append("полюс t = 1 на [0; 2] не отказал")
        except Unresolved:
            pass
        f_num, f_den = IntPoly([1]), IntPoly([1, 0, 0, 1])
        if integral(f_num, f_den, 0, 1) + integral(f_num, f_den, 1, 2) != integral(f_num, f_den, 0, 2):
            out.append("аддитивность ∫₀¹ + ∫₁² ≠ ∫₀² на одной подынтегральной")
    except Exception as e:                                      # noqa: BLE001
        out.append(f"приземления фабрики: {type(e).__name__}: {e}")
    return out
