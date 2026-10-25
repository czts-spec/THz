"""Носитель класса 4 UNITY_ROOT: корни единицы x^n = 1.

Образец корпуса [в корпусе]: группа ℤ/Lℤ и поле порядка — живой одометр
ОП-геометрии, цикл и счёт одним числом. Величина хранится засечкой
колеса: n — число засечек (n ≥ 1), k — номер засечки с каноном k mod n;
значение — e^{2πi·k/n}. Сокращение gcd(n, k) построением НЕ выполняется:
носитель — засечка именно этого колеса; переход к младшему колесу —
только явным запросом reduced().

Вещественные вырождения снимает фабрика of(): k ≡ 0 (mod n) — точная
Q(1), 2k = n — точная Q(-1); всякий живой экземпляр
линейки не веществен. Умножение на одном колесе — одометр: of(n, k1+k2),
переполнение есть цикл, точно; разные колёса сводятся на общее колесо
lcm(n1, n2). Вещественная и мнимая части — по явному запросу re()/im()
через sin_pi: таблица знаменателей {1, 2, 3, 4, 6}, вне её — класс 3
(cyclotomic.py, 04.09: 2cos(2πk/n) с минимальным многочленом).

Сложение корней единицы сведено целиком (04.09): вещественный итог —
ζ^a + ζ^b при b ≡ −a (mod n) есть 2cos(2πa/n) (Q или класс 3), ζ^a − ζ^a
— точный ноль; комплексный итог (в том числе ζ + рациональное) —
структурная запись «пара (re, im)» над линейками 1/3
(complex_pair.ComplexPair, спрос замерен по корпусу: PhiComplex).
Десятичная сторона: величина не вещественна, отрезок
[lo, hi] ⊂ ℚ её не несёт — enclosure() отвечает честным отказом.
Ни одного float: float на входе — TypeError построением.
"""
from __future__ import annotations

from .rational import Q
from math import gcd

from .algebraic import AlgebraicNumber
from .enclosure import Interval
from .errors import Unresolved
from .gamma_value import sin_pi
from .zero_state import Zero

#: Допустимые скаляры бинарных действий. float исключён построением.
_NUM = (int, Q)


def _rec(x: object) -> str:
    """Ф-запись операнда для сообщений отказа."""
    r = getattr(x, "record", None)
    if callable(r):
        return r()
    if isinstance(x, _NUM):
        return f"Φ(C=1 RATIONAL, value={Q(x)})"
    return f"Φ(?={x!r})"


def _as_int(x, name: str) -> int:
    """int → int; bool, float и прочее — TypeError построением."""
    if isinstance(x, bool) or not isinstance(x, int):
        raise TypeError(f"{name} — только int (целая засечка колеса), "
                        "float запрещён")
    return x


class UnityRoot:
    """Засечка k колеса n: значение e^{2πi·k/n}, канон k mod n.

    Живой экземпляр не веществен: k ≡ 0 (единица) и 2k = n (минус
    единица) сняты построением — вещественные вырождения возвращает
    фабрика of() точными Q."""

    __slots__ = ("n", "k")

    def __init__(self, n, k):
        n = _as_int(n, "n")
        k = _as_int(k, "k")
        if n < 1:
            raise ValueError("колесо требует n ≥ 1 засечек")
        k %= n                   # канон: засечка берётся по модулю колеса
        if k == 0:
            raise ValueError("k ≡ 0 (mod n) — вещественная единица, не "
                             "засечка линейки; стройте через of()")
        if 2 * k == n:
            raise ValueError("2k = n — вещественная минус единица, не "
                             "засечка линейки; стройте через of()")
        self.n, self.k = n, k

    # --- построение -----------------------------------------------------

    @classmethod
    def of(cls, n, k) -> Q | UnityRoot:
        """Корень e^{2πi·k/n} с каноном k mod n и снятием вещественных
        вырождений: k ≡ 0 → Q(1), 2k = n → Q(-1)."""
        n = _as_int(n, "n")
        k = _as_int(k, "k")
        if n < 1:
            raise ValueError("колесо требует n ≥ 1 засечек")
        k %= n
        if k == 0:
            return Q(1)
        if 2 * k == n:
            return Q(-1)
        return cls(n, k)

    def reduced(self) -> UnityRoot:
        """Та же точка окружности на младшем колесе n/g, g = gcd(n, k).

        Только явным запросом: построением сокращение не выполняется —
        носитель есть засечка своего колеса. Для живой засечки итог
        всегда живой (вырождения делителем не создаются)."""
        g = gcd(self.n, self.k)
        return self.of(self.n // g, self.k // g)

    # --- умножение / деление: одометр -----------------------------------

    def __mul__(self, other):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__rmul__
        if isinstance(other, UnityRoot):
            if other.n == self.n:
                # одометр: сложение засечек; переполнение = цикл, точно
                return self.of(self.n, self.k + other.k)
            # разные колёса: общее колесо lcm(n1, n2)
            common = self.n * other.n // gcd(self.n, other.n)
            return self.of(common, self.k * (common // self.n)
                           + other.k * (common // other.n))
        if isinstance(other, _NUM):
            q = Q(other)
            if q == 0:
                return Q(0)  # масштаб нолём — точный ноль класса 1
            if q == 1:
                return UnityRoot(self.n, self.k)
            if q == -1:
                return -self        # −1 = e^{iπ}: поворот на полоборота
            raise Unresolved("×", self.record(), _rec(other),
                             hint="модуль корня единицы равен 1: масштаб "
                                  "q ∉ {0, ±1} линейкой не несётся")
        r = getattr(other, "record", None)
        if callable(r):
            raise Unresolved("×", self.record(), r(),
                             hint="правила сведения с этой линейкой нет")
        return NotImplemented

    # умножение коммутативно; чужие линейки получают тот же отказ
    __rmul__ = __mul__

    def inverse(self) -> UnityRoot:
        """Обратный элемент группы: e^{−2πi·k/n} = of(n, −k)."""
        return self.of(self.n, -self.k)

    def __truediv__(self, other):
        if isinstance(other, UnityRoot):
            return self * other.inverse()
        if isinstance(other, _NUM):
            q = Q(other)
            if q == 0:
                raise Unresolved("/", self.record(), _rec(other),
                                 hint="x/0 здесь не сведено")
            if q == 1:
                return UnityRoot(self.n, self.k)
            if q == -1:
                return -self
            raise Unresolved("/", self.record(), _rec(other),
                             hint="модуль корня единицы равен 1: масштаб "
                                  "q ∉ {±1} линейкой не несётся")
        r = getattr(other, "record", None)
        if callable(r):
            raise Unresolved("/", self.record(), r(),
                             hint="правила сведения с этой линейкой нет")
        return NotImplemented

    def __rtruediv__(self, other):
        """q/ζ = q·ζ^{−1}: сведено для q ∈ {0, ±1}."""
        if isinstance(other, _NUM):
            q = Q(other)
            if q == 0:
                return Q(0)
            if q == 1:
                return self.inverse()
            if q == -1:
                return -self.inverse()
            raise Unresolved("/", _rec(other), self.record(),
                             hint="модуль корня единицы равен 1: масштаб "
                                  "q ∉ {±1} линейкой не несётся")
        return NotImplemented

    def __pow__(self, m: int):
        if not isinstance(m, int) or isinstance(m, bool):
            return NotImplemented
        # одометр: m оборотов по k засечек; цикл замыкается точно,
        # m = 0 вырождается фабрикой в точную единицу
        return self.of(self.n, self.k * m)

    def __neg__(self) -> UnityRoot:
        """−ζ = e^{iπ}·ζ: полоборота на общем с колесом 2 колесе."""
        common = self.n if self.n % 2 == 0 else 2 * self.n
        return self.of(common, self.k * (common // self.n) + common // 2)

    def order(self) -> int:
        """Порядок в группе колеса: наименьшее m ≥ 1 с ζ^m = 1."""
        return self.n // gcd(self.n, self.k)

    # --- сложение / вычитание: направленного правила нет ----------------

    def _refuse(self, op: str, other: object, self_left: bool):
        left, right = ((self.record(), _rec(other)) if self_left
                       else (_rec(other), self.record()))
        raise Unresolved(op, left, right,
                         hint="сумма корней единицы вне линейки "
                              "(циклотомическая сумма не несётся); "
                              "сложение сведено только с "
                              "нолём-состоянием")

    def _real_sum(self, other: "UnityRoot", sign: int):
        """ζ^a ± ζ^b на общем колесе N: веществен ⇔ мнимые части гасятся —
        для «+» b ≡ −a (mod N): итог 2cos(2πa/N); для «−» b ≡ a: итог 0.
        Иначе None (комплексная сумма — не на линейке)."""
        n = self.n * other.n // gcd(self.n, other.n)
        a, b = self.k * (n // self.n), other.k * (n // other.n)
        if sign > 0 and (a + b) % n == 0:
            from .cyclotomic import two_cos
            return two_cos(a, n)
        if sign < 0 and (a - b) % n == 0:
            return Q(0)
        return None

    def _pair(self):
        from .complex_pair import ComplexPair
        return ComplexPair.from_wheel(self.n, (self.k,))

    def __add__(self, other):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__radd__
        if isinstance(other, UnityRoot):
            r = self._real_sum(other, +1)
            if r is not None:
                return r
            return self._pair() + other._pair()                 # комплексная сумма — пара (re, im)
        if isinstance(other, _NUM):
            return self._pair() + Q(other)
        self._refuse("+", other, True)

    def __radd__(self, other):
        if isinstance(other, _NUM):
            return self._pair() + Q(other)
        self._refuse("+", other, False)

    def __sub__(self, other):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__rsub__
        if isinstance(other, UnityRoot):
            r = self._real_sum(other, -1)
            if r is not None:
                return r
            return self._pair() - other._pair()
        if isinstance(other, _NUM):
            return self._pair() - Q(other)
        self._refuse("-", other, True)

    def __rsub__(self, other):
        if isinstance(other, _NUM):
            return Q(other) - self._pair()
        self._refuse("-", other, False)

    # --- вещественная и мнимая части: по явному запросу -----------------

    def re(self) -> Q | AlgebraicNumber:
        """Вещественная часть: cos(2πk/n) = sin(π(1/2 − 2k/n)) — таблица
        sin_pi, вне её класс 3 (cyclotomic, 04.09); отказа нет."""
        return sin_pi(Q(1, 2) - Q(2 * self.k, self.n))

    def im(self) -> Q | AlgebraicNumber:
        """Мнимая часть: sin(2πk/n) — та же дорога."""
        return sin_pi(Q(2 * self.k, self.n))

    # --- запись, сравнение, сертификат ----------------------------------

    def record(self) -> str:
        return f"Φ(C=4 UNITY_ROOT, n={self.n}, k={self.k})"

    __repr__ = record

    def __eq__(self, other):
        if isinstance(other, UnityRoot):
            # тождество величин — точка окружности k/n; колесо
            # (носитель) — компонент записи, не величины
            return Q(self.k, self.n) == Q(other.k, other.n)
        if isinstance(other, _NUM) or isinstance(other, AlgebraicNumber):
            return False  # живая засечка линейки не вещественна
        return NotImplemented

    def __hash__(self):
        return hash((UnityRoot, Q(self.k, self.n)))

    def self_check(self, deep: bool = False) -> list:
        """Самопроверка засечки колеса независимой дорогой (очередь п. 18).

        `[вычислено здесь]` Проверяется не тем счётом, каким строилось:
        канон засечки (0 < k < n, 2k ≠ n — живая засечка не вещественна),
        порядок равен n/gcd(n, k), сокращение на младшее колесо возвращает
        ту же точку, хеш согласован с равенством по величине (та же точка на
        удвоенном колесе — та же величина). При `deep` — три встречные
        дороги: одометр `z^порядок = 1`, обратная `z·z⁻¹ = 1` и точный
        Пифагор `re² + im² = 1` (класс 3, таблица sin_pi).

        Граница: десятичной стороны у величины нет — `enclosure` отказывает
        поимённо, и здесь ничего не приближается.
        """
        out = []
        if not (0 < self.k < self.n):
            out.append(f"засечка вне канона: k={self.k}, n={self.n}")
        if 2 * self.k == self.n:
            out.append("2k = n: вещественная минус единица под видом засечки")
        g = gcd(self.n, self.k)
        if self.order() != self.n // g:
            out.append(f"порядок {self.order()} ≠ n/gcd(n,k) = {self.n // g}")
        try:
            red = self.reduced()
            if red != self:
                out.append(f"сокращение на младшее колесо дало другую точку: {_rec(red)}")
        except Exception as e:                              # noqa: BLE001
            out.append(f"сокращение отказало: {type(e).__name__}: {e}")
        try:
            twin = UnityRoot(2 * self.n, 2 * self.k)
            if twin != self or hash(twin) != hash(self):
                out.append("та же точка на удвоенном колесе не равна себе либо хеш не согласован")
        except Exception as e:                              # noqa: BLE001
            out.append(f"удвоение колеса отказало: {type(e).__name__}: {e}")
        if "UNITY_ROOT" not in self.record():
            out.append("запись не читается как засечка колеса")
        if not deep:
            return out
        m = self.order()
        try:
            if self ** m != Q(1):
                out.append(f"одометр: z^{m} ≠ 1")
        except Exception as e:                              # noqa: BLE001
            out.append(f"одометр отказал: {type(e).__name__}")
        try:
            if self * self.inverse() != Q(1):
                out.append("обратная засечка не даёт единицы")
        except Exception as e:                              # noqa: BLE001
            out.append(f"обратная отказала: {type(e).__name__}")
        try:
            re, im = self.re(), self.im()
            if re * re + im * im != Q(1):
                out.append(f"Пифагор точным счётом не даёт единицы: {_rec(re * re + im * im)}")
        except Exception as e:                              # noqa: BLE001
            out.append(f"вещественная и мнимая части не сведены: {type(e).__name__}: {e}")
        return out

    def enclosure(self, terms: int = 30) -> Interval:
        """Десятичная сторона — только сертификатом сетки; величина
        не вещественна, отрезок [lo, hi] ⊂ ℚ её не несёт — отказ."""
        raise Unresolved("enclosure", self.record(),
                         hint="величина не вещественна: отрезок "
                              "[lo, hi] ⊂ ℚ её не несёт; вещественная "
                              "и мнимая части — re()/im() по явному "
                              "запросу")


#: Мнимая единица — засечка четверти колеса: i² = −1 одометром.
I = UnityRoot(4, 1)
