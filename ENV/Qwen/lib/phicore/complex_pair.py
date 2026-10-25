"""Пара (re, im) — комплексная величина как структурная запись над двумя
вещественными линейками (Декартова форма re + im·i).

Спрос замерен 04.09: в корпусе класс `PhiComplex` [в корпусе:
libphi v3.9–v3.10, 19 версий; 171 упоминание в 10 наборах тестов; 23
файла docs/notes] — «мост к результатам, которые PhiNumber отвергает»:
чётный корень отрицательного, опора для UNITY_ROOT («существовал в enum
без арифметики»). Решения корпуса, взятые как есть: Декартова форма;
re²+im² точно; деление на ноль — отказ («комплексная ∞ — единая точка
сферы Римана, не пара (ω, ω)»). Что здесь сверх корпуса: угол π·p/q при
любом q точен — cos/sin класса 3 (cyclotomic.py), а не float APPROX.

Компоненты — любые вещественные носители phicore (Q, класс 3, мономы π…),
лишь бы несли +, −, ×; пара с im = 0 сворачивается в re (is_real()).
Запись: Φ(C=19 STRUCT, пара re=(…) im=(…)).

С 18.09 — целая степень `**` и модуль `magnitude`/`abs` (спрос в корпусе:
PhiComplex.__pow__ и .magnitude; до того у пары двери `**` не было —
TypeError на 3 входах из 3, а abs(3 + 4i) падал TypeError при |z| = 5), и
линейка 3 отдаёт паре отражённые действия: AN + пара, AN − пара, AN × пара,
AN / пара сводятся так же, как пара + AN (до 18.09 — отказ 4 из 4)
`[вычислено здесь]`.
"""
from __future__ import annotations

from .errors import Unresolved
from .rational import Q

_NUM = (int, Q)


def _rec(x) -> str:
    r = getattr(x, "record", None)
    if callable(r):
        return r()
    if isinstance(x, _NUM):
        return f"Φ(C=1 RATIONAL, value={Q(x)})"
    return f"Φ(?={x!r})"


def _real(x):
    if isinstance(x, bool) or isinstance(x, float):
        raise TypeError("float/bool запрещены построением")
    return Q(x) if isinstance(x, int) else x


def _is_zero(x) -> bool:
    if isinstance(x, Q):
        return x == 0
    f = getattr(x, "as_fraction", None)
    if callable(f):
        r = f()
        return r is not None and r == 0
    return False


class ComplexPair:
    """re + im·i над вещественными носителями. Пара, собранная из корней
    единицы, помнит колесо: элемент группового кольца ℚ[ζ_n] — (n, {k: c_k});
    +, −, × таких пар идут на колесе (lcm), re/im считаются лениво, один раз,
    в поле ℚ(2cos(2π/4n)) без роста степени (cyclotomic.unity_sum)."""

    __slots__ = ("_re", "_im", "_wheel")

    def __init__(self, re=None, im=None, wheel=None):
        self._wheel = wheel                     # (n, {k: Q}) либо None
        if wheel is None:
            self._re, self._im = _real(re), _real(im)
        else:
            self._re = None if re is None else _real(re)
            self._im = None if im is None else _real(im)

    def _materialize(self):
        if self._re is None or self._im is None:
            from .cyclotomic import unity_sum
            n, w = self._wheel
            self._re, self._im = unity_sum(n, w)

    @property
    def re(self):
        self._materialize()
        return self._re

    @property
    def im(self):
        self._materialize()
        return self._im

    @classmethod
    def from_wheel(cls, n: int, ks, offset=Q(0)) -> "ComplexPair":
        """Σ c_k ζ_n^k (+ offset): ks — засечки (вес 1) или {k: c_k}; re/im — лениво."""
        w = {}
        items = ks.items() if isinstance(ks, dict) else ((k, 1) for k in ks)
        for k, c in items:
            k = int(k) % n
            w[k] = w.get(k, Q(0)) + Q(c)
        if offset:
            w[0] = w.get(0, Q(0)) + Q(offset)
        w = {k: c for k, c in w.items() if c != 0}
        return cls(wheel=(n, w))

    @staticmethod
    def _lift(wheel, n):
        """Колесо (m, w) на колесо n (m | n)."""
        m, w = wheel
        f = n // m
        return {k * f: c for k, c in w.items()}

    @staticmethod
    def _merge_wheels(a, b, op: str):
        """a ± b, a × b в ℚ[ζ_n] на общем колесе lcm."""
        from math import gcd
        na, nb = a[0], b[0]
        n = na * nb // gcd(na, nb)
        wa, wb = ComplexPair._lift(a, n), ComplexPair._lift(b, n)
        out = {}
        if op in "+-":
            for k, c in wa.items():
                out[k] = out.get(k, Q(0)) + c
            for k, c in wb.items():
                out[k] = out.get(k, Q(0)) + (c if op == "+" else -c)
        else:
            for ka, ca in wa.items():
                for kb, cb in wb.items():
                    k = (ka + kb) % n
                    out[k] = out.get(k, Q(0)) + ca * cb
        return n, {k: c for k, c in out.items() if c != 0}

    # --- построение ---------------------------------------------------------
    @classmethod
    def from_real(cls, x) -> "ComplexPair":
        if isinstance(x, _NUM):
            return cls.from_wheel(1, {}, Q(x))
        return cls(x, Q(0))

    @classmethod
    def from_polar_exact(cls, magnitude, quarter_turns: int) -> "ComplexPair":
        """|z| под углом quarter_turns·π/2: cos/sin ∈ {0, ±1} точно [в корпусе]."""
        q = quarter_turns % 4
        m = _real(magnitude)
        zero = Q(0)
        return (cls(m, zero), cls(zero, m), cls(-m, zero), cls(zero, -m))[q]

    @classmethod
    def from_angle(cls, magnitude, turns) -> "ComplexPair":
        """|z|·e^{2πi·turns}, turns ∈ ℚ: cos/sin — класс 3 через cyclotomic
        (точно при любом знаменателе; корпус здесь давал float APPROX)."""
        from .cyclotomic import cos_pi, sin_pi
        t = Q(turns) * 2                     # угол π·t
        m = _real(magnitude)
        return cls(m * cos_pi(t), m * sin_pi(t))

    @classmethod
    def from_real_root(cls, radicand, q: int, p: int = 1) -> "ComplexPair":
        """radicand^(p/q) при radicand < 0 [в корпусе: PhiComplex.from_real_root]:
        |r|^(p/q) под углом π·p/q — главное значение; q = 2 даёт i·√|r| точно,
        любое q — точно через класс 3 (сверх корпуса)."""
        from .algebraic import AlgebraicNumber
        r = Q(radicand)
        if r >= 0:
            raise ValueError("radicand < 0 — иначе корень веществен, пара не нужна")
        mag = AlgebraicNumber.from_root(q, -r)
        if p != 1:
            mag = mag ** p
        return cls.from_angle(mag, Q(p, 2 * q))

    @classmethod
    def from_unity(cls, z) -> "ComplexPair":
        """Корень единицы как пара (re, im) класса 1/3 с памятью о колесе."""
        return cls.from_wheel(z.n, (z.k,))

    # --- запись, равенство ---------------------------------------------------
    def record(self) -> str:
        return f"Φ(C=19 STRUCT, пара re=({_rec(self.re)}) im=({_rec(self.im)}))"

    __repr__ = record

    def is_real(self) -> bool:
        return _is_zero(self.im)

    def real_or_self(self):
        """re, если im = 0 точно; иначе сама пара."""
        return self.re if self.is_real() else self

    def __eq__(self, other):
        if isinstance(other, ComplexPair):
            return self.re == other.re and self.im == other.im
        if isinstance(other, _NUM) or hasattr(other, "record"):
            return self.is_real() and self.re == _real(other)
        return NotImplemented

    def __hash__(self):
        return hash((ComplexPair, repr(self.re), repr(self.im)))

    # --- арифметика ----------------------------------------------------------
    @staticmethod
    def _wheel_of(other):
        """Колесо операнда: пара с колесом, корень единицы (n, {k: 1}), рациональное — (1, {0: q})."""
        if isinstance(other, ComplexPair):
            return other._wheel
        if hasattr(other, "n") and hasattr(other, "k") and callable(getattr(other, "re", None)):
            return (other.n, {other.k % other.n: Q(1)})          # UnityRoot
        if isinstance(other, _NUM):
            return (1, {0: Q(other)}) if Q(other) != 0 else (1, {})
        return None

    @staticmethod
    def _pair(other) -> "ComplexPair":
        if isinstance(other, ComplexPair):
            return other
        if hasattr(other, "re") and hasattr(other, "im") and callable(other.re):
            return ComplexPair.from_wheel(other.n, (other.k,))     # UnityRoot
        return ComplexPair.from_real(other)

    def _combine(self, other, op: str):
        wa, wb = self._wheel, self._wheel_of(other)
        if wa is not None and wb is not None:
            n, w = self._merge_wheels(wa, wb, op)
            return ComplexPair(wheel=(n, w))
        o = self._pair(other)
        if op == "+":
            return ComplexPair(self.re + o.re, self.im + o.im)
        if op == "-":
            return ComplexPair(self.re - o.re, self.im - o.im)
        a, b, c, d = self.re, self.im, o.re, o.im
        return ComplexPair(a * c - b * d, a * d + b * c)

    def __add__(self, other):
        return self._combine(other, "+")

    __radd__ = __add__

    def __neg__(self):
        if self._wheel is not None:
            n, w = self._wheel
            return ComplexPair(wheel=(n, {k: -c for k, c in w.items()}))
        return ComplexPair(-self.re, -self.im)

    def __sub__(self, other):
        return self._combine(other, "-")

    def __rsub__(self, other):
        return (-self) + other

    def __mul__(self, other):
        return self._combine(other, "*")

    __rmul__ = __mul__

    def conjugate(self) -> "ComplexPair":
        if self._wheel is not None:
            n, w = self._wheel
            return ComplexPair(wheel=(n, {(-k) % n: c for k, c in w.items()}))
        return ComplexPair(self.re, -self.im)

    def magnitude_squared(self):
        """re² + im² — точно на носителях (без корня) [в корпусе]."""
        return self.re * self.re + self.im * self.im

    def __truediv__(self, other):
        o = self._pair(other)
        den = o.magnitude_squared()
        if _is_zero(den):
            raise Unresolved("/", self.record(), o.record(),
                             hint="деление на ноль: комплексная ∞ — единая точка сферы Римана, "
                                  "не пара (ω, ω) [в корпусе: PhiComplex, БАГ-C]")
        num = self * o.conjugate()
        inv = (1 / den) if isinstance(den, Q) else den.inverse()
        return ComplexPair(num.re * inv, num.im * inv)

    def __rtruediv__(self, other):
        return self._pair(other) / self

    def __pow__(self, n):
        """Целая степень — квадрат-и-умножение на своём же `*` (на колесе,
        когда колесо есть) [в корпусе: PhiComplex.__pow__ — «целая степень —
        точное повторное умножение»]. n = 0 — вещественная единица парой;
        n < 0 — степень обратной 1/z (ноль отказывает сферой Римана).
        Показатель Q со знаменателем 1 — та же целая дорога; дробный —
        отказ поимённо: корень из пары есть выбор ветви, главное значение
        живёт у from_angle / from_real_root; float и bool — запрещены
        построением."""
        if isinstance(n, bool) or isinstance(n, float):
            raise TypeError("показатель пары — только int или Q; float/bool запрещены построением")
        if isinstance(n, Q):
            if n.denominator != 1:
                raise Unresolved("^", self.record(), f"Φ(C=1 RATIONAL, value={n})",
                                 hint="дробная степень пары — выбор ветви корня; "
                                      "главное значение строится from_angle / "
                                      "from_real_root, здесь ветвь не назначается")
            n = n.numerator
        if not isinstance(n, int):
            return NotImplemented
        if n == 0:
            return ComplexPair.from_real(Q(1))
        base = self if n > 0 else (1 / self)
        e = -n if n < 0 else n
        result = None
        while e > 0:
            if e & 1:
                result = base if result is None else result * base
            e >>= 1
            if e:
                base = base * base
        return result

    def magnitude(self):
        """|z| = (re² + im²)^(1/2) точно на линейках 1 и 3 [в корпусе:
        PhiComplex.magnitude — EXACT лишь при полном квадрате, иначе APPROX;
        здесь — класс 3 при любом re² + im² из ℚ или класса 3]. Иные
        носители компонент — отказ поимённо; квадрат модуля точен всегда
        (magnitude_squared)."""
        from .algebraic import AlgebraicNumber
        m2 = self.magnitude_squared()
        if isinstance(m2, (Q, AlgebraicNumber)):
            return m2 ** Q(1, 2)
        raise Unresolved("|z|", self.record(),
                         hint="корень из re² + im² сведён на линейках 1 и 3; на этом "
                              "носителе точен только квадрат модуля (magnitude_squared)")

    __abs__ = magnitude

    def self_check(self, deep: bool = False) -> list:
        """Самопроверка пары (re, im) независимой дорогой (очередь п. 18).

        `[вычислено здесь]` Проверяется не тем счётом, каким строилось:
        компоненты вещественны (float/bool дверью `_real` не проходят), канон
        колеса (n ≥ 1, засечки приведены по модулю n, нулевых весов нет),
        `is_real` отвечает ТОЧНОМУ нулю мнимой части, `real_or_self`
        согласован с ним, хеш и равенство согласованы с пересборкой из своих
        же компонент, запись читается как C=19. При `deep` — четыре встречные
        дороги: у пары с колесом материализованные re/im сверяются с
        ПЕРЕСЧЁТОМ `unity_sum` из самого колеса (две разные дороги к одной
        величине — групповое кольцо ℚ[ζ_n] и пара на линейках, и подложенная
        re эту сверку не проходит, хотя запись остаётся с виду законной);
        Пифагор ТОЧНЫМ счётом — z·conj(z) обязано быть вещественным и равным
        re² + im², а у одиночной засечки колеса равным единице; сложение и
        умножение НА КОЛЕСЕ сверяются с теми же действиями по паре (re, im);
        деление обратным умножением — (z/w)·w = z.

        Граница: модуль |z| берётся дорогой линейки 3 (magnitude, с 18.09) и
        сверяется своим квадратом; на носителях вне линеек 1 и 3 модуль —
        отказ, и там проверяется только КВАДРАТ модуля (так и в корпусе).
        Тождество величин на РАЗНЫХ колёсах решается сравнением
        пар, а не сокращением колеса: сокращение живёт у засечки
        (`unity_root.reduced`). Глубокая сверка МАТЕРИАЛИЗУЕТ ленивую пару —
        цена названа, она не бесплатна.
        """
        out = []
        if self._wheel is not None:
            try:
                n, w = self._wheel
                items = list(w.items())
            except Exception:                           # noqa: BLE001
                out.append(f"колесо не пара (n, {{k: c}}): {self._wheel!r}")
                return out
            if isinstance(n, bool) or not isinstance(n, int) or n < 1:
                out.append(f"колесо {n!r} — не целое ≥ 1")
            else:
                for k, c in items:
                    if isinstance(k, bool) or not isinstance(k, int) or not (0 <= k < n):
                        out.append(f"засечка {k!r} не приведена по модулю {n}")
                    if c == 0:
                        out.append(f"нулевой вес при ζ_{n}^{k}: член не член")
        for name, x in (("re", self._re), ("im", self._im)):
            if x is None:
                if self._wheel is None:
                    out.append(f"{name} не задана и колеса нет: пары нет")
                continue
            if isinstance(x, bool) or isinstance(x, float):
                out.append(f"{name} — {type(x).__name__}: float/bool запрещены построением")
        if out:
            return out                      # дальше материализовать нечего
        try:
            re, im = self.re, self.im
        except Exception as e:                          # noqa: BLE001
            out.append(f"материализация пары отказала: {type(e).__name__}: {e}")
            return out
        if self.is_real() != _is_zero(im):
            out.append(f"is_real не отвечает точному нулю мнимой части: im = {_rec(im)}")
        ros = self.real_or_self()
        if self.is_real():
            if ros != re:
                out.append("real_or_self не вернул re у вещественной пары")
        elif ros is not self:
            out.append("real_or_self отдал не саму пару у невещественной")
        if "C=19" not in self.record():
            out.append("запись не читается как пара класса 19")
        try:
            twin = ComplexPair(re, im)
            if twin != self or hash(twin) != hash(self):
                out.append("пересборка из своих же компонент не равна себе либо хеш не согласован")
        except Exception as e:                          # noqa: BLE001
            out.append(f"пересборка из своих же компонент отказала: {type(e).__name__}: {e}")
        if not deep:
            return out
        if self._wheel is not None:
            n, w = self._wheel
            try:
                from .cyclotomic import unity_sum
                re2, im2 = unity_sum(n, w)
            except Exception as e:                      # noqa: BLE001
                out.append(f"пересчёт колеса unity_sum отказал: {type(e).__name__}: {e}")
            else:
                if re2 != re:
                    out.append(f"колесо и пара разошлись по re: unity_sum даёт {_rec(re2)}, "
                               f"записано {_rec(re)}")
                if im2 != im:
                    out.append(f"колесо и пара разошлись по im: unity_sum даёт {_rec(im2)}, "
                               f"записано {_rec(im)}")
        try:
            conj = self.conjugate()
            prod = self * conj
            m2 = self.magnitude_squared()
        except Exception as e:                          # noqa: BLE001
            out.append(f"сопряжение либо квадрат модуля отказали: {type(e).__name__}: {e}")
            return out
        if not _is_zero(prod.im):
            out.append(f"z·conj(z) не вещественно: im = {_rec(prod.im)}")
        if prod.re != m2:
            out.append(f"Пифагор точным счётом: z·conj(z) = {_rec(prod.re)} "
                       f"≠ re² + im² = {_rec(m2)}")
        if conj.conjugate() != self:
            out.append("двойное сопряжение не тождественно")
        try:
            mag = self.magnitude()
        except Unresolved:
            pass                                    # носитель вне линеек 1/3 — граница названа в magnitude
        else:
            if mag * mag != m2:
                out.append(f"модуль дорогой линейки 3: |z|² = {_rec(mag * mag)} "
                           f"≠ re² + im² = {_rec(m2)}")
        if self._wheel is not None and len(self._wheel[1]) == 1 \
                and next(iter(self._wheel[1].values())) == Q(1) and m2 != Q(1):
            out.append(f"одиночная засечка колеса: |ζ|² = {_rec(m2)} ≠ 1")
        if self._wheel is not None:
            if self + self != ComplexPair(re + re, im + im):
                out.append("сложение на колесе не сходится со сложением по паре (re, im)")
            if self * self != ComplexPair(re * re - im * im, re * im + im * re):
                out.append("умножение на колесе не сходится с умножением по паре (re, im)")
        w = ComplexPair(Q(1), Q(-2))                    # |w|² = 5 ≠ 0: полюса нет
        try:
            div = self / w
        except Exception as e:                          # noqa: BLE001
            out.append(f"деление на 1 − 2i отказало: {type(e).__name__}: {e}")
        else:
            if div * w != self:
                out.append(f"деление обратным умножением: (z/w)·w = {(div * w).record()} ≠ z")
        return out
