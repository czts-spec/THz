"""Линейка класса 19 ITERATION_LIMIT: предел итерации lim fⁿ(x0).

Собственного носителя линейка НЕ хранит: её направленное правило —
перевод предела в ту линейку, где предел живёт. Дорог две, и обе
объявлены границей входа:

  * f = num/den ∈ ℚ(x) (`iteration_limit`) — предел идёт в класс 3
    POLY_ROOT или в класс 1 RATIONAL по §3: кандидаты суть корни
    G(x) = num(x) − x·den(x), ибо всякая неподвижная точка f(r) = r с
    den(r) ≠ 0 есть корень G, и обратно;
  * f = coeff·e^{rate·x} над величинами дома (`ExpMap`,
    `iteration_limit_map`) — предел идёт в класс 17 LAMBERT_W или в
    класс 1: перепись x = c·e^{a·x} ⟺ (−a·x)·e^{−a·x} = −a·c даёт
    x = −W(−a·c)/a, то есть обращение той же экспоненты, которым живёт
    класс 17. Точная дорога есть при a = −1 (x = W(c) без деления) и
    всюду, где W(−a·c) рационально (`phicore/lambert.py: exact_w`);
    иное — отказ поимённо, не приближение.

Потому модуль несёт фабрики, сертификаты и отказы, а не класс-носитель;
record() результата — запись той линейки, где предел живёт.

ДВОЙСТВЕННОСТЬ линейки — неподвижная точка: r = lim fⁿ(x0) обязан
удовлетворять f(r) = r, то есть num(r) = r·den(r) на рациональной
дороге и c·e^{a·r} = r на экспоненциальной. Это обращение хода
итерации — тот же ход, каким класс 17 LAMBERT_W обращает x·e^x = y:
у f = e^{−x} неподвижная точка ЕСТЬ W(1) (Ω-константа), и равенство
проверяется объектами, не разрядами. В приёмке тождество проверяется
точной подстановкой в арифметике класса 3 (рациональная дорога) и
класса 5 (экспоненциальная), не сравнением разрядов.

Сертификат сходимости — по последним 8 итерациям точного прогона:

  * итерация x_{i+1} = num(x_i)/den(x_i) считается в Q без
    единого округления; полюс den(x_i) = 0 — отказ Unresolved;
  * x_{i+1} = x_i — точная неподвижная точка (корень G): хвост
    постоянен, предел равен x_i и возвращается Q;
  * свободная от квадратов часть G, счёт корней Штурмом
    (phicore.algebraic) на отрезке [min хвоста − 1; max + 1]; требуется
    ровно один корень r там, а при нескольких — ровно один в
    сужающейся [x_k; x_{k+1}]-оболочке (последняя пара хвоста,
    раздвинутая на два её зазора);
  * |x_i − r| строго убывает по всему хвосту: для рационального r —
    точным сравнением Q, для алгебраического — через уточняемый
    изолирующий отрезок (AlgebraicNumber.refine): пара сертифицирована,
    когда наибольшее расстояние до отрезка у следующей итерации меньше
    наименьшего у предыдущей; неразличимое за потолок уже́ний — отказ,
    не догадка о стороне.

Нет сертификата за steps шагов — отказ «сходимость не сертифицирована
за N шагов» с явными Ф-записями f и x0, не правдоподобный предел и не
молчание. Десятичная сторона — только сертификатом сетки
(limit_enclosure → Interval), по явному запросу. Ни одного float:
float на входе — TypeError построением.
"""
from __future__ import annotations

from .rational import Q

from .algebraic import AlgebraicNumber, count_roots, squarefree_part
from .enclosure import Interval, exp_enclosure
from .errors import Unresolved
from .exp_log import ExpValue
from .intpoly import IntPoly
from .lambert import WValue, exact_w, lambert_w

#: Допустимые скаляры. float исключён построением.
_NUM = (int, Q)

#: Путь модуля в сообщениях отказа: отказ называет своё место.
_HERE = "phicore/iteration.py"

#: Сертификат читает последние 8 итераций хвоста.
_TAIL = 8

#: Потолок уже́ния изолирующего отрезка на одно сравнение расстояний.
_MAX_HALVING = 512

#: Решётка разрядности итерации: рост дроби сверх — честный отказ,
#: не зависание счёта (числитель + знаменатель, в битах).
_MAX_BITS = 1 << 16


def _rec(x: object) -> str:
    """Ф-запись операнда для сообщений отказа."""
    r = getattr(x, "record", None)
    if callable(r):
        return r()
    if isinstance(x, _NUM):
        return f"Φ(C=1 RATIONAL, value={Q(x)})"
    return f"Φ(?={x!r})"


def _phi_f(num: IntPoly, den: IntPoly) -> str:
    """Явная Ф-запись отображения f = num/den для сообщений отказа.

    Носителя класс 19 не держит — запись собирается на месте отказа."""
    ns = ",".join(str(c) for c in num.c)
    ds = ",".join(str(c) for c in den.c)
    return f"Φ(C=19 ITERATION_LIMIT, num=[{ns}], den=[{ds}])"


def _check_poly(p: object, name: str) -> None:
    """num/den — только IntPoly (ℤ[x]); прочее — TypeError построением."""
    if not isinstance(p, IntPoly):
        raise TypeError(f"{name} — только IntPoly (ℤ[x]); float и голые "
                        "списки запрещены построением")


def fixed_point_poly(num: IntPoly, den: IntPoly) -> IntPoly:
    """G(x) = num(x) − x·den(x) — направленное правило к классу 3.

    Корни G — в точности неподвижные точки f = num/den вне полюсов:
    num(r)/den(r) = r ⟺ num(r) − r·den(r) = 0 при den(r) ≠ 0.
    Построение точное: сдвиг den на одну степень и вычитание в ℤ[x]."""
    _check_poly(num, "num")
    _check_poly(den, "den")
    return num - IntPoly([0] + list(den.c))


# ---------- сертификат хвоста ----------

def _dist_lo(x: Q, lo: Q, hi: Q) -> Q:
    """Наименьшее расстояние от точки x до отрезка [lo; hi]."""
    if x < lo:
        return lo - x
    if x > hi:
        return x - hi
    return Q(0)


def _dist_hi(x: Q, lo: Q, hi: Q) -> Q:
    """Наибольшее расстояние от точки x до величины, несомой [lo; hi]."""
    return max(hi - x, x - lo)


def _certify(sf: IntPoly, tail: "list[Q]"):
    """Сертификат сходимости хвоста к одному изолированному корню sf.

    Возвращает пару (точный предел | None, причина отказа). Предел —
    Q, когда корень рационален (распознан самой линейкой 3),
    иначе AlgebraicNumber(sf, lo, hi). sf — свободная от квадратов
    часть G, степень ≥ 1."""
    lo_seg = min(tail) - 1
    hi_seg = max(tail) + 1
    n = count_roots(sf, lo_seg, hi_seg)
    if n == 0:
        return None, "корней G на [min хвоста − 1; max + 1] нет"
    lo, hi = lo_seg, hi_seg
    if n > 1:
        # сужающаяся оболочка последней пары: [x_k; x_{k+1}] ± два зазора
        a, b = tail[-2], tail[-1]
        d = abs(b - a)
        lo = max(lo_seg, min(a, b) - 2 * d)
        hi = min(hi_seg, max(a, b) + 2 * d)
        if count_roots(sf, lo, hi) != 1:
            return None, "в [x_k; x_{k+1}]-оболочке не ровно один корень G"
    root = AlgebraicNumber(sf, lo, hi)
    r = root.as_fraction()
    if r is not None:
        # рациональный корень: расстояния точны, сравнение — Q
        prev = None
        for x in tail:
            cur = abs(x - r)
            if prev is not None and not cur < prev:
                return None, "|x_i − r| не убывает строго (r рационален)"
            prev = cur
        return r, ""
    # алгебраический корень: сравнение через уточняемый изолирующий отрезок
    for u, v in zip(tail, tail[1:]):
        m = (u + v) / 2
        if root.lo < m <= root.hi and sf.eval(m) == 0:
            # r — точная середина пары: расстояния равны (цикл 2 вокруг r)
            return None, "r — середина пары итераций: расстояния равны"
        for _ in range(_MAX_HALVING):
            if _dist_lo(u, root.lo, root.hi) > _dist_hi(v, root.lo, root.hi):
                break  # пара сертифицирована: |x_{i+1} − r| < |x_i − r|
            if _dist_hi(u, root.lo, root.hi) < _dist_lo(v, root.lo, root.hi):
                return None, "|x_i − r| строго растёт на паре хвоста"
            root.refine()
        else:
            return None, (f"сравнение |x_i − r| не сертифицировано за "
                          f"{_MAX_HALVING} уже́ний отрезка")
    return root, ""


# ---------- фабрика линейки ----------

def iteration_limit(num: IntPoly, den: IntPoly, x0, steps: int = 60):
    """lim fⁿ(x0) для f(x) = num(x)/den(x) — по направленным правилам.

    Точная итерация Q до steps шагов; полюс den(x_i) = 0 — отказ.
    Точное попадание в неподвижную точку (x_{i+1} = x_i) возвращает
    Q немедленно: хвост постоянен, G(x_i) = 0. Иначе после
    каждого шага с полным хвостом читается сертификат (_certify);
    при успехе — ТОЧНЫЙ предел: Q для рационального корня G,
    AlgebraicNumber(sf, lo, hi) для алгебраического. Нет сертификата
    за steps шагов — отказ Unresolved с явными Ф-записями f и x0
    и причиной последней неудачи; результат никогда не приближается."""
    _check_poly(num, "num")
    _check_poly(den, "den")
    if isinstance(x0, bool) or not isinstance(x0, _NUM):
        raise TypeError("x0 — только int/Q, float запрещён")
    if isinstance(steps, bool) or not isinstance(steps, int):
        raise TypeError("steps — только int")
    if steps < 1:
        raise ValueError("steps ≥ 1")
    if den.is_zero():
        raise ValueError("den=0 не задаёт рациональной функции")
    x0 = Q(x0)
    g = fixed_point_poly(num, den)
    if g.degree == 0 and not g.is_zero():
        raise Unresolved("lim fⁿ", _phi_f(num, den), _rec(x0),
                         hint="G = num − x·den — ненулевая константа: "
                              "неподвижных точек нет, кандидатов предела нет")
    # G ≡ 0 ⟺ f = id вне полюсов: первый же шаг даёт x₁ = x₀
    sf = None if g.is_zero() else squarefree_part(g)
    x = x0
    xs = [x0]
    reason = f"хвост короче {_TAIL} итераций — сертификат не читается"
    for i in range(steps):
        d = den.eval(x)
        if d == 0:
            raise Unresolved("lim fⁿ", _phi_f(num, den), _rec(x0),
                             hint=f"полюс: den(x_{i}) = 0 при x_{i} = {x} — "
                                  "итерация не определена")
        x_next = num.eval(x) / d
        if x_next == x:
            # точная неподвижная точка: G(x) = 0, предел достигнут
            return x
        bits = x_next.numerator.bit_length() + x_next.denominator.bit_length()
        if bits > _MAX_BITS:
            raise Unresolved("lim fⁿ", _phi_f(num, den), _rec(x0),
                             hint=f"сходимость не сертифицирована за {i + 1} "
                                  f"шагов: разрядность итерации превысила "
                                  f"решётку 2^16 битов — рост без сжатия")
        xs.append(x_next)
        x = x_next
        if sf is not None and len(xs) >= _TAIL:
            val, reason = _certify(sf, xs[-_TAIL:])
            if val is not None:
                return val
    raise Unresolved("lim fⁿ", _phi_f(num, den), _rec(x0),
                     hint=f"сходимость не сертифицирована за {steps} шагов; "
                          f"{reason}")


# ---------- десятичная сторона: сертификат сетки по явному запросу ----------

def limit_enclosure(value, bits: int = 32) -> Interval:
    """Сертификат сетки [lo, hi] ⊂ ℚ для предела, по явному запросу.

    Предел живёт в классах 1/3: рациональный — точечный отрезок;
    алгебраический — его изолирующий отрезок, суженный уточнением
    (AlgebraicNumber.refine — уже́ние отрезка носителя) до ширины
    ≤ 2^−bits. Иной вход сертификата не несёт — отказ Unresolved."""
    if isinstance(bits, bool) or not isinstance(bits, int):
        raise TypeError("bits — только int")
    if bits < 1:
        raise ValueError("bits ≥ 1")
    if isinstance(value, bool):
        raise TypeError("value — Q или AlgebraicNumber")
    if isinstance(value, float):
        raise TypeError("float запрещён построением")
    if isinstance(value, _NUM):
        return Interval.point(Q(value))
    if isinstance(value, AlgebraicNumber):
        target = Q(1, 2 ** bits)
        while value.hi - value.lo > target:
            value.refine()
        return Interval(value.lo, value.hi)
    raise Unresolved("enclosure", _rec(value),
                     hint="сертификат сетки сведён только для классов 1 и 3")


# ---------- f из величин дома: отображение класса 19 над классом 5 ----------

def _admit_coeff(coeff):
    """Вход coeff: пара ℚ\\{0} либо класс 5 EXPONENTIAL; прочее — поимённо."""
    if isinstance(coeff, bool):
        raise TypeError("coeff — bool числом в доме не считается; дверь несёт "
                        "int/Q и класс 5 ExpValue")
    if isinstance(coeff, float):
        raise TypeError("coeff — float запрещён построением; дверь несёт "
                        "int/Q и класс 5 ExpValue")
    if isinstance(coeff, _NUM):
        q = Q(coeff)
        if q == 0:
            raise ValueError("coeff=0 не есть отображение класса 19: "
                             "f ≡ 0 неподвижна в самом нуле, и предел "
                             "живёт в ℚ без итерации")
        return q
    if isinstance(coeff, ExpValue):
        return coeff
    raise Unresolved("ExpMap", _rec(coeff),
                     hint=f"{_HERE}: coeff несут пара ℚ и класс 5 "
                          "EXPONENTIAL; для иного носителя переписи "
                          "x = c·e^{a·x} ⟺ (−a·x)·e^{−a·x} = −a·c нет, "
                          "потому и точной неподвижной точки нет")


def _admit_rate(rate) -> Q:
    """Вход rate: пара ℚ\\{0}. Показатель вне ℚ·x линейка не несёт."""
    if isinstance(rate, bool) or not isinstance(rate, _NUM):
        raise TypeError("rate — только int/Q, float запрещён построением")
    q = Q(rate)
    if q == 0:
        raise ValueError("rate=0 даёт постоянное f = coeff: предел есть сам "
                         "coeff и живёт в классах 1/5 без итерации")
    return q


class ExpMap:
    """f(x) = coeff·e^{rate·x} — отображение класса 19 над величинами дома.

    Объявленная граница входа:

      * coeff — пара ℚ\\{0} либо класс 5 EXPONENTIAL (`ExpValue`);
      * rate — пара ℚ\\{0};
      * float/bool — TypeError построением; иной носитель coeff (класс 6
        LOGARITHMIC, класс 10, класс 17 …) — отказ Unresolved поимённо с
        путём модуля.

    Чем граница держится: неподвижная точка x = c·e^{a·x} переписывается
    ТОЧНО как (−a·x)·e^{−a·x} = −a·c, то есть x = −W(−a·c)/a — обращение
    той же экспоненты, которым живёт класс 17. Сумма c·e^{a·x} + b,
    показатель вне ℚ·x, основание вне e такой переписи не имеют, и
    точной неподвижной точки в доме за ними нет.

    f от рациональной точки — ТОЧНАЯ величина класса 5 (`eval_exact`),
    f от отрезка — сертификат сетки (`eval_interval`). Точной итерации по
    точкам здесь быть не может: уже второй шаг требует показателя
    a·c·e^{a·x} — носителя такому показателю дом не несёт, потому
    сходимость сертифицируется на отрезках (`iteration_limit_map`)."""

    __slots__ = ("coeff", "rate")

    def __init__(self, coeff, rate):
        self.coeff = _admit_coeff(coeff)
        self.rate = _admit_rate(rate)

    # --- запись ---
    def record(self) -> str:
        c = (str(self.coeff) if isinstance(self.coeff, Q)
             else self.coeff.record())
        return (f"Φ(C=19 ITERATION_LIMIT, f={c}·e^({self.rate}·x))")

    __repr__ = record

    def __eq__(self, other):
        if isinstance(other, ExpMap):
            return self.coeff == other.coeff and self.rate == other.rate
        return NotImplemented

    def __hash__(self):
        return hash((ExpMap, self.coeff, self.rate))

    # --- значение: точное в классе 5, сертификат на отрезке ---
    def eval_exact(self, x):
        """f(x) при x ∈ ℚ — ТОЧНО, величиной класса 5 (или ℚ при
        вырождении показателя): coeff·e^{rate·x} собирается арифметикой
        линейки 5, где e^a·e^b = e^{a+b}. Ни одного разряда."""
        if isinstance(x, bool) or not isinstance(x, _NUM):
            raise TypeError("x — только int/Q, float запрещён построением")
        return self.coeff * ExpValue.make(1, self.rate * Q(x))

    def eval_interval(self, iv: Interval, terms: int = 60) -> Interval:
        """f на отрезке — сертификат сетки: e^{rate·[lo; hi]} монотонен,
        потому концы берутся наружу (низ от низа, верх от верха), а
        coeff умножается своим сертификатом (точкой при coeff ∈ ℚ)."""
        if not isinstance(iv, Interval):
            raise TypeError("iv — только Interval")
        if isinstance(terms, bool) or not isinstance(terms, int):
            raise TypeError("terms — только int")
        if terms < 1:
            raise ValueError("terms ≥ 1")
        e = iv * self.rate
        ex = Interval(exp_enclosure(e.lo, terms).lo,
                      exp_enclosure(e.hi, terms).hi)
        if isinstance(self.coeff, Q):
            return ex * self.coeff
        return ex * self.coeff.enclosure(terms)

    # --- направленное правило: неподвижная точка в класс 17 либо в ℚ ---
    def fixed_point(self, branch: int = 0):
        """Точная неподвижная точка f: класс 17 LAMBERT_W либо ℚ.

        x = c·e^{a·x} ⟺ (−a·x)·e^{−a·x} = −a·c ⟺ x = −W(−a·c)/a, где
        y = −a·c есть точная величина дома (ℚ либо класс 5):

          * W(y) рационально точной дорогой линейки 17 (`exact_w`:
            y = k·e^k → W = k, в том числе точка ветвления y = −1/e) →
            x = −k/a ТОЧНО в ℚ, и подстановка это подтверждает в
            арифметике класса 5 (`check_fixed`);
          * a = −1 → x = W(c) — НОСИТЕЛЬ класса 17 без деления: у
            f = e^{−x} это ровно W(1), Ω-константа;
          * иначе — отказ поимённо: величины −W(y) и W(y)/a линейка 17 не
            несёт по построению (`phicore/lambert.py`: `WValue.__neg__`,
            `WValue.__truediv__` отказывают), потому точного носителя у
            неподвижной точки в доме нет.

        branch выбирает ветвь W: ветвь 0 несёт притягивающую точку, ветвь
        −1 — вторую, отталкивающую (она неподвижна точно так же, но
        пределом итерации не бывает)."""
        y = -self.rate * self.coeff          # точная величина дома
        k = exact_w(y, branch)
        if k is not None:
            return -k / self.rate
        if self.rate == -1:
            return lambert_w(self.coeff, branch)
        raise Unresolved("неподвижная точка", self.record(), _rec(y),
                         hint=f"{_HERE}: точная неподвижная точка есть "
                              f"−W(y)/a при a={self.rate}, y = −a·c; "
                              "величин −W и W/a линейка 17 не несёт "
                              "(phicore/lambert.py: WValue.__neg__, "
                              "WValue.__truediv__). Точная дорога есть при "
                              "a = −1 либо при рациональном W(y) (exact_w)")

    def check_fixed(self, x) -> bool:
        """ДВОЙСТВЕННОСТЬ подстановкой: f(x) = x ТОЧНО, для x ∈ ℚ.

        Сравнение идёт объектами класса 5: f(x) есть coeff·e^{rate·x}, и
        равенство с рациональным x возможно лишь при вырождении
        показателя в ноль — иначе величина трансцендентна [внешний
        источник, теорема Линдемана] и равной x не бывает.

        Для носителя класса 17 точной подстановки нет: e^{−W} линейка 5
        не несёт (`phicore/exp_log.py: exp` отказывает на классе 17) —
        отказ поимённо; сертификатная дорога — `iteration_limit_map`."""
        if isinstance(x, bool):
            raise TypeError("x — bool числом в доме не считается")
        if isinstance(x, _NUM):
            return self.eval_exact(x) == Q(x)
        if isinstance(x, WValue):
            raise Unresolved("f(x) = x", self.record(), x.record(),
                             hint=f"{_HERE}: точной подстановки носителя "
                                  "класса 17 нет — e^{−W} линейка 5 не "
                                  "несёт (phicore/exp_log.py: exp на классе "
                                  "17 отказывает); сертификатная дорога — "
                                  "iteration_limit_map")
        raise Unresolved("f(x) = x", self.record(), _rec(x),
                         hint=f"{_HERE}: подстановка сведена для пары ℚ")


    # --- самопроверка отображения (ось самоприменения, 08.10) --------------

    def self_check(self, deep: bool = False) -> list:
        """Самопроверка отображения встречными дорогами (ось «AIAΣA
        спрашивает AIAΣA», слово автора 04.09; очередь строка 18; 08.10).
        Список находок, пусто — здорово; исключений не бросает.

        `[вычислено здесь]` (1) слоты: coeff — ℚ без нуля либо класс 5, rate —
        ℚ без нуля; (2) неподвижная точка (`fixed_point`): рациональная —
        подстановка f(x) = x ТОЧНО в арифметике класса 5 (`check_fixed`),
        носитель класса 17 — его собственная самопроверка (линейка 17) и
        сертификатная дорога f(X) ∩ X ≠ ∅ на отрезке X = enclosure(16)
        носителя (`eval_interval` дверью exp против сертификата W, строенного
        бисекцией той же двери, но по другой формуле); отказ точной
        неподвижной точки — объявленная граница, не находка. При deep —
        `iteration_limit_map` от нуля обязан сойтись к той же точке."""
        out = []
        c, a = self.coeff, self.rate
        if not isinstance(c, (Q, ExpValue)) or (isinstance(c, Q) and c == 0):
            return [f"слот coeff вне границы входа: {c!r}"]
        if not isinstance(a, Q) or a == 0:
            return [f"слот rate вне границы входа: {a!r}"]
        fp = None
        try:
            fp = self.fixed_point(0)
            if isinstance(fp, Q):
                if not self.check_fixed(fp):
                    out.append(f"рациональная неподвижная точка {fp} не выдерживает подстановки f(x) = x")
            elif isinstance(fp, WValue):
                out += [f"носитель W: {f}" for f in fp.self_check()]
                X = fp.enclosure(16)
                fX = self.eval_interval(X, 60)
                if fX.intersect(X) is None:
                    out.append(f"сертификатная дорога: f(X) = {fX.record()} не пересекает X = {X.record()}")
            else:
                out.append(f"неподвижная точка чужого носителя: {type(fp).__name__}")
        except Unresolved:
            pass                                                # граница: точной неподвижной точки в доме нет
        except Exception as e:                                  # noqa: BLE001
            out.append(f"неподвижная точка не построена: {type(e).__name__}: {e}")
        if deep and fp is not None:
            try:
                lim = iteration_limit_map(self, 0)
                same = (lim == fp) if isinstance(fp, Q) else (isinstance(lim, WValue) and lim == fp)
                if not same:
                    out.append(f"предел итерации {lim!r} не совпал с неподвижной точкой {fp!r}")
            except Unresolved:
                pass                                            # сходимость не сертифицирована — граница инструмента
            except Exception as e:                              # noqa: BLE001
                out.append(f"предел итерации: {type(e).__name__}: {e}")
        return out


def iteration_limit_map(f: ExpMap, x0, steps: int = 60, terms: int = 60,
                        shrink: int = 10 ** 6, bits: int = 60):
    """lim fⁿ(x0) для f = coeff·e^{rate·x} — сертификат плюс ТОЧНЫЙ итог.

    Точной итерации по точкам нет и быть не может: f(x) при x ∈ ℚ есть
    величина класса 5, и второй шаг требует показателя a·c·e^{a·x} —
    носителя такому показателю дом не несёт. Потому сходимость
    сертифицируется на ОТРЕЗКАХ, а итог возвращается точный:

      * вилка I — последняя пара хвоста из _TAIL интервальных шагов от
        ⟦x0; x0⟧, РАЗДВИНУТАЯ на два её зазора в каждую сторону (тот же
        канон, каким рациональная дорога изолирует корень G: хвост даёт
        место, зазор даёт запас на недобранный остаток);
      * f(I) ⊆ I — неподвижная точка на I существует (Брауэр на
        отрезке); иначе отказ «f(I) ⊄ I» — так отказывает и монотонно
        расходящийся ход, и ход с |f′| > 1, ушедший в цикл 2;
      * steps пересечений I ← I ∩ f(I) обязаны сузить отрезок не менее
        чем в shrink раз (по умолчанию 10⁶) — иначе отказ «сжатия нет».
        Это заявление о ПРЕДМЕТЕ: при |f′| = 1 (нейтральная точка, ровно
        случай c = e, rate = −1) и при |f′| > 1 сжатия нет и предела
        итерации нет, хотя точная неподвижная точка есть;
      * точный итог (`ExpMap.fixed_point`) обязан лежать внутри
        сертификата: Q — точкой, носитель класса 17 — пересечением
        своего сертификата ширины ≤ 2^−bits. Расхождение — отказ, не
        подгонка.

    Ни одного float: float на входе — TypeError построением."""
    if not isinstance(f, ExpMap):
        raise TypeError("f — только ExpMap (класс 19 над величинами дома); "
                        "рациональная дорога живёт в iteration_limit")
    if isinstance(x0, bool) or not isinstance(x0, _NUM):
        raise TypeError("x0 — только int/Q, float запрещён")
    for name, v in (("steps", steps), ("terms", terms), ("shrink", shrink),
                    ("bits", bits)):
        if isinstance(v, bool) or not isinstance(v, int):
            raise TypeError(f"{name} — только int")
        if v < 1:
            raise ValueError(f"{name} ≥ 1")
    x0 = Q(x0)
    cur = Interval.point(x0)
    tail = []
    for _ in range(_TAIL):
        cur = f.eval_interval(cur, terms)
        tail.append(cur)
    a, b = tail[-2], tail[-1]
    lo, hi = min(a.lo, b.lo), max(a.hi, b.hi)
    gap = hi - lo
    cur = Interval(lo - 2 * gap, hi + 2 * gap)
    if not cur.encloses(f.eval_interval(cur, terms)):
        raise Unresolved("lim fⁿ", f.record(), _rec(x0),
                         hint=f"f(I) ⊄ I на вилке {cur.record()} — "
                              "неподвижная точка на отрезке не изолирована, "
                              "Брауэр не применим; сходимость не "
                              "сертифицирована")
    start = cur.width()
    for i in range(steps):
        nxt = cur.intersect(f.eval_interval(cur, terms))
        if nxt is None:
            raise Unresolved("lim fⁿ", f.record(), _rec(x0),
                             hint=f"пересечение I ∩ f(I) пусто на шаге {i} — "
                                  "сертификат противоречив, сходимость не "
                                  "сертифицирована")
        cur = nxt
    if not (start == 0 or cur.width() * shrink < start):
        raise Unresolved("lim fⁿ", f.record(), _rec(x0),
                         hint=f"сжатия нет: за {steps} пересечений отрезок "
                              f"{cur.record()} не стал у́же исходного в "
                              f"{shrink} раз — при |f′| ≥ 1 предела итерации "
                              "нет, хотя точная неподвижная точка есть "
                              "(ExpMap.fixed_point)")
    val = f.fixed_point()
    inside = (cur.contains(val) if isinstance(val, Q)
              else cur.intersect(val.enclosure(bits)) is not None)
    if not inside:
        raise Unresolved("lim fⁿ", f.record(), _rec(x0),
                         hint=f"точный итог {_rec(val)} вне сертификата "
                              f"{cur.record()} — расхождение дорог, не "
                              "подгонка")
    return val


# ---------- самопроверка линейки (ось самоприменения, 08.10) ----------

def self_check(deep: bool = False) -> list:
    """Самопроверка линейки класса 19 встречными дорогами; список находок,
    пусто — здорова; исключений не бросает.

    `[вычислено здесь]` Рациональная дорога: Ньютон для √2 — f = (x² + 2)/(2x)
    от x0 = 2 — предел в классе 3, его сертификат в квадрате обязан нести 2
    (дорога класса 3 против внешнего числа), а предел как корень G = num −
    x·den обязан обращать G в нуль по изолирующему отрезку (Штурм на нём
    даёт ровно один корень); f = (x + 6)/2 сходится к 6 точно; f = x
    неподвижна в x0; G — ненулевая константа — отказ. Экспоненциальная
    дорога: f = e^{−x} — неподвижная точка есть носитель W(1) (Ω),
    отображение проверяет себя; f = e^{−1}·e^{x} — точная точка 1 с
    подстановкой; deep — ещё f = 2·e^{−x} и Ньютон для √3. Граница — окно."""
    out = []
    try:
        r = iteration_limit(IntPoly([2, 0, 1]), IntPoly([0, 2]), Q(2))
        if not isinstance(r, AlgebraicNumber):
            out.append(f"Ньютон √2: предел {r!r} не в классе 3")
        else:
            iv = limit_enclosure(r, 40)
            sq = iv * iv
            if not sq.contains(Q(2)):
                out.append(f"Ньютон √2: сертификат в квадрате не несёт 2: {sq.record()}")
            g = fixed_point_poly(IntPoly([2, 0, 1]), IntPoly([0, 2]))
            if count_roots(squarefree_part(g), r.lo, r.hi) != 1:
                out.append("Ньютон √2: на изолирующем отрезке предела не ровно один корень G")
        if iteration_limit(IntPoly([6, 1]), IntPoly([2]), Q(0)) != Q(6):
            out.append("f = (x + 6)/2: предел ≠ 6")
        if iteration_limit(IntPoly([0, 1]), IntPoly([1]), Q(5)) != Q(5):
            out.append("f = x: предел ≠ x0")
        try:
            iteration_limit(IntPoly([1, 1]), IntPoly([1]), Q(0))
            out.append("f = x + 1: не отказал при G — ненулевой константе")
        except Unresolved:
            pass
        if deep:
            r3 = iteration_limit(IntPoly([3, 0, 1]), IntPoly([0, 2]), Q(2))
            iv3 = limit_enclosure(r3, 40)
            if not (iv3 * iv3).contains(Q(3)):
                out.append(f"Ньютон √3: сертификат в квадрате не несёт 3")
    except Exception as e:                                      # noqa: BLE001
        out.append(f"рациональная дорога: {type(e).__name__}: {e}")
    maps = [ExpMap(Q(1), Q(-1)), ExpMap(ExpValue.make(Q(1), Q(-1)), Q(1))]
    if deep:
        maps.append(ExpMap(Q(2), Q(-1)))
    for m in maps:
        try:
            out += [f"{m.record()}: {f}" for f in m.self_check(deep)]
        except Exception as e:                                  # noqa: BLE001
            out.append(f"{m.record()}: {type(e).__name__}: {e}")
    try:
        omega = ExpMap(Q(1), Q(-1)).fixed_point(0)
        if not isinstance(omega, WValue) or omega != lambert_w(Q(1), 0):
            out.append(f"f = e^{{−x}}: неподвижная точка {omega!r} не есть W(1)")
        one = ExpMap(ExpValue.make(Q(1), Q(-1)), Q(1)).fixed_point(0)
        if one != Q(1):
            out.append(f"f = e^{{−1}}·e^x: неподвижная точка {one!r} ≠ 1")
    except Exception as e:                                      # noqa: BLE001
        out.append(f"экспоненциальная дорога: {type(e).__name__}: {e}")
    return out
