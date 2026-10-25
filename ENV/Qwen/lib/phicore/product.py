"""Носитель класса 8 INFINITE_PRODUCT: x = Π fₙ с рациональными множителями.

Линейка держит произведение СТРОЕНИЕМ, а не готовым пределом: множитель
fₙ = factor_num(n)/factor_den(n) — отношение целочисленных многочленов
номера (intpoly.IntPoly), счёт от старта n0 ∈ ℤ. Каждый множитель и
каждое частное произведение точны: Q, накопление с кэшем, ни
одного float — float на входе даёт TypeError построением.

Горячий цикл накопления — на парах целых с отложенной нормализацией
(этап F выселения Q): множитель — пара (числитель, знаменатель
> 0), цепочка частных произведений перемножает числители и знаменатели
как есть, gcd берётся один раз при выдаче Q наружу — форма
контейнера v6 (сокращение по запросу, не на каждом действии), результат
бит-в-бит тот же. Страховка от разбухания: при росте пары за _PAIR_BITS
бит — разовое сокращение gcd, и счёт продолжается (нормализация
отложенная, не отменённая).

Нули проверяются ЛЕНИВО, в точке запроса: factor_den(n) = 0 — полюс
множителя, factor_num(n) = 0 — нулевой множитель, гасящий произведение;
оба получают отказ Unresolved при обращении к factor(n), не сплошным
перебором номеров при построении.

Двойственность 7↔8, приземлённая в класс 6: лог произведения есть
сумма логов, Σ ln fₙ = ln Π fₙ. Для конечного N свёртка суммы
LogValue.make(1, fₙ) по правилам класса 6 обязана дать
LogValue.make(1, partial_product(N)) либо точный ноль — метод
log_partial(N). Свёртка ведётся от старшего номера к младшему: канон
класса 6 сворачивает совершенную степень аргумента в coeff, и частное
произведение, оказавшееся точной степенью, при свёртке от младших
номеров уводит coeff от единицы — дальнейшее сложение получает честный
отказ (замер на Валлисе: P₄ = (128/105)², свёртка слева гаснет на
пятом члене). Порядок — часть правила сведения; отказ класса 6, если
он всё же случился, не подменяется.

Десятичная сторона — только сертификатом сетки (enclosure), по явному
запросу: частное произведение до m ≥ N плюс строгое двустороннее
обрамление хвоста мажорантой |fₙ − 1| ≤ 2A/(B·n²); при
deg(factor_num − factor_den) > deg(factor_den) − 2 мажоранта хвоста
не сходится — честный отказ Unresolved, не приближение.

Арифметика произведений (05.09, двойственность 7↔8 телескопом): слияние
строений при одном старте n0 — как прежде (Π fₙ · Π gₙ = Π fₙgₙ); всё
остальное — сложение и вычитание с числом, рядом, интегралом и другим
произведением, умножение на число ≠ 1, на ряд, на интеграл и на произведение
с иным стартом, деление на число, q/Π, знак — идёт через ряд класса 7:
Π = p_{n0} + Σ (p_N − p_{N−1}), разности частных произведений — D-конечны
порядка 1 (dfinite.DFiniteSeries.from_product), логарифм не нужен; итог —
ряд класса 7 с сертификатом от произведения. Деление произведений разных
стартов и Π⁰ — отказ по теореме (обращение ряда не замкнуто). Отказ —
Unresolved с явными Ф-записями обоих операндов, не ложный класс.
"""
from __future__ import annotations

from .rational import Q
from math import gcd

from .enclosure import Interval, exp_enclosure
from .errors import Unresolved
from .exp_log import LogValue
from .intpoly import IntPoly
from .zero_state import Zero

#: Допустимые скаляры. float исключён построением.
_NUM = (int, Q)

#: Порог страховки: пара целых длиннее этого числа бит сокращается
#: разовым gcd прямо в цикле — отложенная нормализация, не отменённая.
_PAIR_BITS = 1 << 20


def _rec(x: object) -> str:
    """Ф-запись операнда для сообщений отказа."""
    r = getattr(x, "record", None)
    if callable(r):
        return r()
    if isinstance(x, _NUM):
        return f"Φ(C=1 RATIONAL, value={Q(x)})"
    return f"Φ(?={x!r})"


def _as_int(n, name: str) -> int:
    """Номер множителя: только int; bool и float исключены построением."""
    if isinstance(n, bool) or not isinstance(n, int):
        raise TypeError(f"{name} — только int")
    return n


class InfiniteProductValue:
    """Π fₙ: fₙ = factor_num(n)/factor_den(n) ∈ ℚ, n = n0, n0+1, …

    factor(n) и partial_product(N) точны (Q) и кэшируются;
    канона записи у линейки нет — разные строения могут нести один
    предел, поэтому тождество решается только совпадением строения,
    иное сравнение — отказ."""

    __slots__ = ("factor_num", "factor_den", "n0", "_factors", "_partials")

    def __init__(self, factor_num, factor_den, n0: int):
        if not isinstance(factor_num, IntPoly) \
                or not isinstance(factor_den, IntPoly):
            raise TypeError("factor_num и factor_den — только intpoly.IntPoly")
        if factor_den.is_zero():
            raise ValueError("factor_den — нулевой многочлен: множитель "
                             "fₙ не определён ни в одной точке")
        if factor_num.is_zero():
            raise ValueError("factor_num — нулевой многочлен: каждый "
                             "множитель гасит произведение — ноль живёт "
                             "в zero_state")
        self.factor_num, self.factor_den = factor_num, factor_den
        self.n0 = _as_int(n0, "n0")
        #: кэш множителей: {n: (числитель, знаменатель > 0)} — пары
        #: целых; заполняется лениво, по точке запроса.
        self._factors: dict[int, tuple[int, int]] = {}
        #: кэш частных произведений: пары целых (числитель, знаменатель
        #: > 0); _partials[i] несёт Π_{n=n0}^{n0+i} fₙ, gcd отложен.
        self._partials: list[tuple[int, int]] = []

    # --- множители и частные произведения: точный счёт с кэшем ----------

    def _factor_pair(self, n: int) -> tuple:
        """Множитель fₙ парой целых (числитель, знаменатель > 0), кэш.

        Проверка нулей — та же ленивая, что у factor; gcd пары не
        берётся — нормализация отложена до выдачи Q наружу."""
        f = self._factors.get(n)
        if f is None:
            den = self.factor_den.eval(n)
            if den == 0:
                raise Unresolved("factor", self.record(),
                                 hint=f"factor_den({n})=0 — полюс множителя: "
                                      f"произведение в этой точке не "
                                      f"определено")
            num = self.factor_num.eval(n)
            if num == 0:
                raise Unresolved("factor", self.record(),
                                 hint=f"factor_num({n})=0 — нулевой "
                                      f"множитель гасит произведение; ноль "
                                      f"живёт в zero_state")
            if den < 0:
                num, den = -num, -den
            f = (num, den)
            self._factors[n] = f
        return f

    def factor(self, n: int) -> Q:
        """Точный множитель fₙ = factor_num(n)/factor_den(n) (Q).

        Проверка нулей — ленивая, в самой точке n: полюс (нуль
        знаменателя) и нулевой множитель (нуль числителя) получают
        отказ Unresolved здесь, а не перебором номеров при построении.
        Внутри — пара целых; gcd один раз здесь, при выдаче наружу."""
        n = _as_int(n, "n")
        if n < self.n0:
            raise ValueError(f"n должен быть ≥ n0={self.n0}, получено {n}")
        return Q(*self._factor_pair(n))

    def partial_product(self, N: int) -> Q:
        """Точное частное произведение Π_{n=n0}^{N} fₙ (Q).

        Цепочка — на парах целых: числители и знаменатели перемножаются
        как есть, gcd один раз здесь, при выдаче Q наружу
        (бит-в-бит то же значение); страховка: пара длиннее _PAIR_BITS
        бит сокращается разовым gcd, и счёт продолжается."""
        N = _as_int(N, "N")
        if N < self.n0:
            raise ValueError(f"N должен быть ≥ n0={self.n0}: пустой "
                             f"отрезок произведения")
        while len(self._partials) <= N - self.n0:
            k = self.n0 + len(self._partials)
            fn, fd = self._factor_pair(k)
            pn, pd = self._partials[-1] if self._partials else (1, 1)
            pn *= fn
            pd *= fd
            if max(pn.bit_length(), pd.bit_length()) > _PAIR_BITS:
                g = gcd(pn, pd)
                if g > 1:
                    pn, pd = pn // g, pd // g
            self._partials.append((pn, pd))
        return Q(*self._partials[N - self.n0])

    # --- двойственность 7↔8: сумма логов сворачивается в лог произведения

    def log_partial(self, N: int):
        """Σ_{n=n0}^{N} LogValue.make(1, fₙ) по правилам класса 6.

        Двойственность структурная и точная: свёртка обязана дать
        LogValue.make(1, partial_product(N)) либо точный ноль. Порядок
        свёртки — от старшего номера к младшему — часть правила: канон
        класса 6 сворачивает совершенную степень аргумента в coeff, и
        частное произведение, оказавшееся точной степенью, при свёртке
        слева уводит coeff от единицы, после чего сложение честно
        гаснет. Точный ноль из вырождения make (fₙ = 1, ln 1 = 0)
        нейтрален и пропускается; отказ класса 6, если он случился на
        самом сложении, не подменяется."""
        N = _as_int(N, "N")
        if N < self.n0:
            raise ValueError(f"N должен быть ≥ n0={self.n0}: пустой "
                             f"отрезок произведения")
        acc = LogValue.make(1, self.factor(N))
        for n in range(N - 1, self.n0 - 1, -1):
            term = LogValue.make(1, self.factor(n))
            if isinstance(term, Q):
                continue  # fₙ = 1: ln 1 — точный ноль, нейтрален
            if isinstance(acc, Q):
                acc = term  # накопитель — точный ноль из вырождения make
                continue
            acc = term + acc
        return acc

    # --- сертификат сетки -----------------------------------------------

    def enclosure(self, N: int = 40) -> Interval:
        """Сертификат сетки [lo, hi] ⊂ ℚ: частное произведение + строгое
        двустороннее обрамление хвоста.

        Мажоранта отклонения множителя от единицы: d = factor_num −
        factor_den; при n ≥ 1 верно |d(n)| ≤ A·n^{dd} и
        |factor_den(n)| ≥ n^{dq−1}·(B·n − C), где A — сумма модулей
        коэффициентов d, B — модуль старшего коэффициента знаменателя,
        C — сумма модулей остальных; значит при dd ≤ dq − 2 и B·n > 2C

            |fₙ − 1| ≤ A/(n·(B·n − C)) ≤ 2A/(B·n²).

        Берётся m ≥ N с B·m > 2C и S̄ = 2A/(B·m) < 1; телескопом
        1/n² ≤ 1/(n−1) − 1/n сумма Σ_{n>m} |fₙ−1| ≤ S̄, и хвост
        обрамлён двусторонне: Π_{n>m} fₙ ≥ 1 − S̄ (неравенство
        Вейерштрасса) и Π_{n>m} fₙ ≤ exp(S̄) — верхняя сторона берётся
        сертификатом exp_enclosure. За m знаменатель отделён от нуля
        той же мажорантой — хвост определён без перебора номеров.

        При dd > dq − 2 мажоранта Σ|fₙ−1| не сходится — честный отказ
        Unresolved (расходимость мажоранты), не приближение."""
        N = _as_int(N, "N")
        if N < self.n0:
            raise ValueError(f"N должен быть ≥ n0={self.n0}: пустой "
                             f"отрезок произведения")
        d = self.factor_num - self.factor_den
        q = self.factor_den
        dq = q.degree
        if not d.is_zero() and d.degree > dq - 2:
            raise Unresolved("enclosure", self.record(),
                             hint="deg(factor_num − factor_den) > "
                                  "deg(factor_den) − 2: мажоранта Σ|fₙ−1| "
                                  "не сходится — расходимость мажоранты, "
                                  "хвост произведения не сведён")
        A = sum(abs(c) for c in d.c)
        B = abs(q.c[-1])
        C = sum(abs(c) for c in q.c[:-1])
        # m: с этой точки мажоранта строга (B·m > 2C) и мала (S̄ < 1)
        m = max(N, self.n0, 2, (2 * C) // B + 1, (2 * A) // B + 1)
        s_bar = Q(2 * A, B * m)
        head = self.partial_product(m)
        tail = Interval(1 - s_bar, exp_enclosure(s_bar).hi)
        return tail * head

    # --- арифметика: слияние строений либо дорога через ряд класса 7 -----

    def _series(self):
        """Π как D-конечный ряд (телескоп разностей частных произведений)."""
        from .dfinite import DFiniteSeries
        return DFiniteSeries.from_product(self)

    @staticmethod
    def _lift(other):
        from .dfinite import DFiniteSeries
        if isinstance(other, (bool, float)):
            return None
        return DFiniteSeries._lift(other)

    def _refuse(self, op: str, other, hint: str):
        if isinstance(other, _NUM) or callable(getattr(other, "record", None)):
            raise Unresolved(op, self.record(), _rec(other), hint=hint)
        return NotImplemented

    def __add__(self, other):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__radd__
        o = self._lift(other)
        if o is None:
            return self._refuse("+", other, "слагаемое не поднимается к D-конечному ряду")
        return self._series() + o

    def __sub__(self, other):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__rsub__
        o = self._lift(other)
        if o is None:
            return self._refuse("-", other, "вычитаемое не поднимается к D-конечному ряду")
        return self._series() - o

    def __radd__(self, other):
        o = self._lift(other)
        if o is None:
            return NotImplemented
        return o + self._series()

    def __rsub__(self, other):
        o = self._lift(other)
        if o is None:
            return NotImplemented
        return o - self._series()

    def _merge(self, other, invert: bool):
        """Π fₙ · Π gₙ = Π (fₙ·gₙ) при одном старте n0 — слияние строений
        (произведение отношений многочленов номера); invert — Π (fₙ/gₙ).
        Нули и полюсы множителей по-прежнему проверяются лениво."""
        gn, gd = ((other.factor_den, other.factor_num) if invert
                  else (other.factor_num, other.factor_den))
        return InfiniteProductValue(self.factor_num * gn,
                                    self.factor_den * gd, self.n0)

    def __mul__(self, other):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__rmul__
        if isinstance(other, InfiniteProductValue):
            if other.n0 == self.n0:
                return self._merge(other, invert=False)
            return self._series() * other._series()      # разные старты — свёртка Коши рядов
        if isinstance(other, _NUM) and not isinstance(other, bool):
            if Q(other) == 1:
                return self
            return self._series().scale(other)           # масштаб — в ряду (у строения слота нет)
        o = self._lift(other)
        if o is None:
            return self._refuse("×", other, "множитель не поднимается к D-конечному ряду")
        return self._series() * o

    __rmul__ = __mul__

    def inverse(self) -> "InfiniteProductValue":
        """1/Π fₙ = Π (1/fₙ): обращение строением — обмен числителя и
        знаменателя множителя."""
        return InfiniteProductValue(self.factor_den, self.factor_num, self.n0)

    def __truediv__(self, other):
        if isinstance(other, InfiniteProductValue):
            if other.n0 == self.n0:
                return self._merge(other, invert=True)
            return self._refuse("/", other, "слияние строений — только при "
                                            "одном старте n0; частное рядов "
                                            "не замкнуто (теорема)")
        if isinstance(other, _NUM) and not isinstance(other, bool):
            q = Q(other)
            if q == 1:
                return self
            if q == 0:
                raise Unresolved("/", self.record(), _rec(other), hint="Π/0 не сведено")
            return self._series().scale(1 / q)
        return self._refuse("/", other, "деление произведения на носитель: "
                                        "обращение ряда не замкнуто (теорема)")

    def __rtruediv__(self, other):
        if isinstance(other, _NUM) and not isinstance(other, bool):
            q = Q(other)
            if q == 1:
                return self.inverse()
            if q == 0:
                return Q(0)
            return self.inverse()._series().scale(q)     # q/Π = q·Π(1/fₙ)
        return NotImplemented

    def __pow__(self, n):
        if isinstance(n, bool) or not isinstance(n, int):
            return NotImplemented
        if n == 1:
            return self
        if n >= 2:
            # (Π fₙ)^n = Π fₙ^n — степень многочленов множителя
            num, den = self.factor_num, self.factor_den
            for _ in range(n - 1):
                num, den = num * self.factor_num, den * self.factor_den
            return InfiniteProductValue(num, den, self.n0)
        if n <= -1:
            return self.inverse() ** (-n)
        # n=0 не сводится: предел произведения может быть нулевым
        # (Π n/(n+1) гаснет), x⁰ не решается строением записи.
        raise Unresolved("^", self.record(), _rec(n),
                         hint="степень 0 произведения не сведена")

    def __neg__(self):
        return self._series().scale(-1)                  # знак — в ряду: строение слота не имеет

    # --- запись и тождество ---------------------------------------------

    # --- самопроверка носителя (ось самоприменения, 08.10) ----------------

    def self_check(self, deep: bool = False) -> list:
        """Самопроверка носителя встречными дорогами (ось «AIAΣA спрашивает
        AIAΣA», слово автора 04.09; очередь строка 18; 08.10). Список
        находок, пусто — здоров; исключений не бросает.

        `[вычислено здесь]` (1) слоты: строение — два IntPoly, ни один не
        нулевой, n0 — целое; (2) кэш против свежей свёртки: каждая кэшированная
        пара множителя обязана равняться factor_num(n)/factor_den(n) по
        прямому счёту, каждое кэшированное частное произведение — свёртке
        точных Q от n0 ВНЕ кэша (дорога без отложенной нормализации);
        (3) двойственность 7↔8 на первых трёх номерах: log_partial(N) —
        LogValue.make(1, partial_product(N)) либо точный ноль; (4) при
        сходящейся мажоранте — сертификаты enclosure(N) и enclosure(2N) оба
        несут предел, потому пересекаются, и ширина не растёт; расходимость
        мажоранты — законный отказ, не находка. При deep — ещё N = 40 и
        сертификат накрывает дальнейшее частное произведение P_{4N}."""
        out = []
        fn, fd, n0 = self.factor_num, self.factor_den, self.n0
        if not isinstance(fn, IntPoly) or not isinstance(fd, IntPoly):
            return [f"слоты строения чужого носителя: {type(fn).__name__}, {type(fd).__name__}"]
        if fn.is_zero() or fd.is_zero():
            return ["нулевой многочлен в строении (нуль живёт в zero_state)"]
        if isinstance(n0, bool) or not isinstance(n0, int):
            return [f"слот n0 не целое: {n0!r}"]
        try:
            for n, (a, b) in sorted(self._factors.items()):
                if b <= 0:
                    out.append(f"кэш множителя f_{n}: знаменатель {b} ≤ 0")
                elif Q(a, b) != Q(fn.eval(n), fd.eval(n)):
                    out.append(f"кэш множителя f_{n} = {Q(a, b)} против прямого счёта {Q(fn.eval(n), fd.eval(n))}")
            acc = Q(1)
            for i, (a, b) in enumerate(self._partials):
                acc = acc * Q(fn.eval(n0 + i), fd.eval(n0 + i))
                if b <= 0 or Q(a, b) != acc:
                    out.append(f"кэш частного произведения P_{n0 + i} = {Q(a, b)} против свежей свёртки {acc}")
                    break
        except Exception as e:                                  # noqa: BLE001
            out.append(f"кэш не сверен: {type(e).__name__}: {e}")
        for N in range(n0, n0 + 3):
            try:
                lp = self.log_partial(N)
                want = LogValue.make(1, self.partial_product(N))
                if lp != want:
                    out.append(f"двойственность 7↔8 при N = {N}: Σ ln f_n = {lp!r} ≠ ln P_N = {want!r}")
            except Unresolved as e:
                out.append(f"двойственность 7↔8 при N = {N}: отказ {e}")
            except Exception as e:                              # noqa: BLE001
                out.append(f"двойственность 7↔8 при N = {N}: {type(e).__name__}: {e}")
        d = fn - fd
        if d.is_zero() or d.degree <= fd.degree - 2:
            for N in ((20, 40) if deep else (20,)):
                try:
                    e1, e2 = self.enclosure(N), self.enclosure(2 * N)
                    if e1.intersect(e2) is None:
                        out.append(f"сертификаты enclosure({N}) и enclosure({2 * N}) не пересекаются: {e1.record()} против {e2.record()}")
                    if e2.width() > e1.width():
                        out.append(f"ширина сертификата растёт с N: {e1.width()} → {e2.width()}")
                    if deep and not e1.contains(self.partial_product(4 * N)):
                        out.append(f"enclosure({N}) не накрывает P_{4 * N}")
                except Exception as e:                          # noqa: BLE001
                    out.append(f"сертификат при N = {N}: {type(e).__name__}: {e}")
        return out

    def record(self) -> str:
        num = ",".join(map(str, self.factor_num.c))
        den = ",".join(map(str, self.factor_den.c))
        return (f"Φ(C=8 INFINITE_PRODUCT, f=[{num}]/[{den}], n0={self.n0})")

    __repr__ = record

    def __eq__(self, other):
        if isinstance(other, InfiniteProductValue):
            if (self.factor_num == other.factor_num
                    and self.factor_den == other.factor_den
                    and self.n0 == other.n0):
                return True
            # Канона записи нет: различие строений не решает различие
            # пределов — тождество величин не сведено, отказ вместо
            # ложного «не равно».
            raise Unresolved("=", self.record(), other.record(),
                             hint="тождество пределов разных строений "
                                  "не сведено")
        if isinstance(other, _NUM):
            raise Unresolved("=", self.record(), _rec(other),
                             hint="предел произведения может быть "
                                  "рациональным — тождество с числом не "
                                  "решается строением")
        return NotImplemented

    def __hash__(self):
        return hash((InfiniteProductValue, tuple(self.factor_num.c),
                     tuple(self.factor_den.c), self.n0))


# ---------- самопроверка линейки (ось самоприменения, 08.10) ----------

def self_check(deep: bool = False) -> list:
    """Самопроверка линейки класса 8 встречными дорогами; список находок,
    пусто — здорова; исключений не бросает.

    `[вычислено здесь]` Окно: произведение Валлиса Π 4n²/(4n² − 1) — сертификат
    enclosure(60) обязан нести π/2 дорогой класса 10 (`pi_enclosure`), чужой
    линейки [внешний источник: формула Валлиса]; Π (1 − 1/n²) от n = 2 —
    частное произведение равно (N + 1)/(2N) телескопом, сертификат несёт 1/2;
    телескоп Π (n + 1)/n — конечная сторона точна, обрамление хвоста законно
    отказывает (мажоранта расходится). Носители окна проверяют себя
    (`InfiniteProductValue.self_check`). Граница — окно."""
    out = []
    try:
        from .enclosure import pi_enclosure
        wallis = InfiniteProductValue(IntPoly([0, 0, 4]), IntPoly([-1, 0, 4]), 1)
        out += [f"Валлис: {f}" for f in wallis.self_check(deep)]
        enc = wallis.enclosure(60)
        if enc.intersect(pi_enclosure(30) / 2) is None:
            out.append(f"Валлис: enclosure(60) не несёт π/2 — {enc.record()} против pi_enclosure(30)/2")
    except Exception as e:                                      # noqa: BLE001
        out.append(f"Валлис: {type(e).__name__}: {e}")
    try:
        half = InfiniteProductValue(IntPoly([-1, 0, 1]), IntPoly([0, 0, 1]), 2)
        out += [f"Π(1 − 1/n²): {f}" for f in half.self_check(deep)]
        for N in (5, 17):
            if half.partial_product(N) != Q(N + 1, 2 * N):
                out.append(f"Π(1 − 1/n²): P_{N} = {half.partial_product(N)} ≠ (N + 1)/(2N) телескопом")
        if not half.enclosure(30).contains(Q(1, 2)):
            out.append(f"Π(1 − 1/n²): enclosure(30) не несёт 1/2: {half.enclosure(30).record()}")
    except Exception as e:                                      # noqa: BLE001
        out.append(f"Π(1 − 1/n²): {type(e).__name__}: {e}")
    try:
        tele = InfiniteProductValue(IntPoly([1, 1]), IntPoly([0, 1]), 1)
        if tele.partial_product(9) != Q(10):
            out.append(f"телескоп: P_9 = {tele.partial_product(9)} ≠ 10")
        try:
            tele.enclosure(10)
            out.append("телескоп: enclosure не отказал при расходящейся мажоранте")
        except Unresolved:
            pass
    except Exception as e:                                      # noqa: BLE001
        out.append(f"телескоп: {type(e).__name__}: {e}")
    return out
