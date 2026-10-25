"""Поле ℚ(x) — нормальная форма записи над ТРАНСЦЕНДЕНТНОЙ базой.

Запись Φα = m·α^E с показателем класса 6 (zech.py) есть значение
m·α^q·Π(1 + kᵢ·α^{dᵢ})^{cᵢ}; при x = β = α^(1/N) это рациональная функция
от x с коэффициентами-парами. Над трансцендентной базой (e^y при
рациональном y≠0 — Эрмит 1873/Линдеман 1882; q·π^e — Линдеман 1882) β
трансцендентна, и ℚ(β) ≅ ℚ(x): два выражения равны как ЧИСЛА тогда и
только тогда, когда равны как рациональные функции. Нормальная форма —
несократимая дробь P/Q с унитарным знаменателем; сравнение записей с
разной историей сложений — сравнение форм, не структуры (план 03.09,
шаг 1). Над рациональной базой поле другое — ℚ[x]/(x^N−α) (alpha.py).

Ни одного float; знак значения — интервальным Горнером по сертификату
сетки β (корень из отрезка — enclosure.root_enclosure).
"""
from __future__ import annotations

from math import gcd

from .enclosure import Interval
from .errors import Unresolved
from .rational import Q

Poly = tuple  # коэффициенты по возрастанию степени, старший ≠ 0; () — ноль


def _strip(p) -> Poly:
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return tuple(p)


def _add(a: Poly, b: Poly) -> Poly:
    n = max(len(a), len(b))
    return _strip((a[i] if i < len(a) else Q(0)) + (b[i] if i < len(b) else Q(0))
                  for i in range(n))


def _neg(a: Poly) -> Poly:
    return tuple(-x for x in a)


def _mul(a: Poly, b: Poly) -> Poly:
    if not a or not b:
        return ()
    out = [Q(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j, y in enumerate(b):
            if y != 0:
                out[i + j] += x * y
    return _strip(out)


def _scale(a: Poly, k: Q) -> Poly:
    return () if k == 0 else tuple(x * k for x in a)


class ZeroInQx(Unresolved, ZeroDivisionError):
    """Обращение нуля или нулевой знаменатель в ℚ(x): отказ с Ф-записью
    (README, правило 1) и одновременно ZeroDivisionError — прежние обработчики
    остаются в силе (очередь п. 47, 10.09)."""


def _show_poly(p: Poly) -> str:
    if not p:
        return "0"
    parts = []
    for i, c in enumerate(p):
        if c == 0:
            continue
        parts.append(f"{c}" if i == 0 else (f"{c}·x" if i == 1 else f"{c}·x^{i}"))
    return " + ".join(parts) or "0"


def _divmod(u: Poly, v: Poly) -> tuple[Poly, Poly]:
    if not v:
        raise ZeroInQx("/", _show_poly(tuple(u)), "0", hint="деление многочлена на нулевой многочлен в ℚ(x)")
    u = list(u)
    q = [Q(0)] * max(1, len(u) - len(v) + 1)
    while len(u) >= len(v):
        c = u[-1] / v[-1]
        d = len(u) - len(v)
        q[d] = c
        for i, y in enumerate(v):
            u[i + d] -= c * y
        u = list(_strip(u))
        if not u:
            break
    return _strip(q), tuple(u)


def _gcd(a: Poly, b: Poly) -> Poly:
    while b:
        a, b = b, _divmod(a, b)[1]
    return a


def _monic(a: Poly) -> tuple[Poly, Q]:
    """(a/lead, lead)."""
    if not a:
        return (), Q(1)
    lead = a[-1]
    return tuple(x / lead for x in a), lead


def _horner(a: Poly, x):
    """Значение многочлена в точке x (пара или Interval)."""
    if not a:
        return Interval.point(0) if isinstance(x, Interval) else Q(0)
    acc = Interval.point(a[-1]) if isinstance(x, Interval) else a[-1]
    for c in reversed(a[:-1]):
        acc = acc * x + c
    return acc


class RatFunc:
    """Несократимая дробь P/Q ∈ ℚ(x): Q унитарен, gcd(P, Q) = 1; ноль — P=()."""

    __slots__ = ("num", "den")

    def __init__(self, num, den=(Q(1),)):
        num, den = _strip(Q(c) for c in num), _strip(Q(c) for c in den)
        if not den:
            raise ZeroInQx("P/Q", _show_poly(num), "0", hint="знаменатель — нулевой многочлен в ℚ(x)")
        if not num:
            den = (Q(1),)
        else:
            g = _gcd(num, den)
            if len(g) > 1:
                num, den = _divmod(num, g)[0], _divmod(den, g)[0]
            den, lead = _monic(den)
            num = tuple(x / lead for x in num)
        object.__setattr__(self, "num", num)
        object.__setattr__(self, "den", den)

    def __setattr__(self, *_):
        raise AttributeError("RatFunc неизменяема")

    # --- построение ------------------------------------------------------
    @classmethod
    def const(cls, c) -> "RatFunc":
        return cls((Q(c),))

    @classmethod
    def mono(cls, n: int, c=1) -> "RatFunc":
        """c·x^n; отрицательный n — c/x^(−n)."""
        if n >= 0:
            return cls((Q(0),) * n + (Q(c),))
        return cls((Q(c),), (Q(0),) * (-n) + (Q(1),))

    # --- устройство ------------------------------------------------------
    @property
    def is_zero(self) -> bool:
        return not self.num

    @property
    def is_monomial(self):
        """(c, n), если P/Q = c·x^n (n ∈ ℤ); иначе None."""
        nz_num = [i for i, c in enumerate(self.num) if c != 0]
        nz_den = [i for i, c in enumerate(self.den) if c != 0]
        if len(nz_num) != 1 or len(nz_den) != 1:
            return None
        return self.num[nz_num[0]] / self.den[nz_den[0]], nz_num[0] - nz_den[0]

    def reduce_index(self, N: int) -> tuple["RatFunc", int]:
        """x = β = α^(1/N): если встречаются лишь степени, кратные g | N, —
        та же функция от α^(1/(N/g)). Канон для хеша и сравнения."""
        g = N
        for p in (self.num, self.den):
            for i, c in enumerate(p):
                if c != 0:
                    g = gcd(g, i)
        if g <= 1:
            return self, N
        return RatFunc(self.num[::g], self.den[::g]), N // g

    def __eq__(self, other):
        return (isinstance(other, RatFunc) and self.num == other.num
                and self.den == other.den)

    def __hash__(self):
        return hash(("ℚ(x)", self.num, self.den))

    def record(self) -> str:
        def show(p):
            if not p:
                return "0"
            parts = []
            for i, c in enumerate(p):
                if c == 0:
                    continue
                parts.append(f"{c}" if i == 0 else
                             (f"{c}·x^{i}" if c != 1 else f"x^{i}"))
            return " + ".join(parts)
        return show(self.num) if self.den == (Q(1),) \
            else f"({show(self.num)}) / ({show(self.den)})"

    __repr__ = record

    # --- арифметика поля -------------------------------------------------
    def __add__(self, other):
        if not isinstance(other, RatFunc):
            return NotImplemented
        return RatFunc(_add(_mul(self.num, other.den), _mul(other.num, self.den)),
                       _mul(self.den, other.den))

    def __neg__(self):
        return RatFunc(_neg(self.num), self.den)

    def __sub__(self, other):
        if not isinstance(other, RatFunc):
            return NotImplemented
        return self + (-other)

    def __mul__(self, other):
        if isinstance(other, RatFunc):
            return RatFunc(_mul(self.num, other.num), _mul(self.den, other.den))
        if isinstance(other, (int, Q)) and not isinstance(other, bool):
            return RatFunc(_scale(self.num, Q(other)), self.den)
        return NotImplemented

    __rmul__ = __mul__

    def inverse(self) -> "RatFunc":
        if self.is_zero:
            raise ZeroInQx("1/x", self.record(), hint="обращение нуля в ℚ(x)")
        return RatFunc(self.den, self.num)

    def __truediv__(self, other):
        if isinstance(other, RatFunc):
            return self * other.inverse()
        return NotImplemented

    def __pow__(self, n: int) -> "RatFunc":
        if not isinstance(n, int) or isinstance(n, bool):
            return NotImplemented
        base = self.inverse() if n < 0 else self
        n = abs(n)
        out, acc = RatFunc.const(1), base
        while n:
            if n & 1:
                out = out * acc
            n >>= 1
            if n:
                acc = acc * acc
        return out

    # --- знак по сертификату сетки аргумента ---------------------------
    def sign_at(self, x: Interval):
        """Знак P(x)/Q(x) при x ∈ отрезке: ±1, если оба множителя
        отделились от нуля; None — отрезок слишком широк (уточнять)."""
        if not isinstance(x, Interval):
            raise TypeError(f"sign_at: нужен отрезок Interval, получено {type(x).__name__}")
        if self.is_zero:
            return 0
        s = 1
        for p in (self.num, self.den):
            iv = _horner(p, x)
            if iv.lo > 0:
                continue
            if iv.hi < 0:
                s = -s
                continue
            return None
        return s

    def value_at(self, x: Interval) -> Interval:
        return _horner(self.num, x) / _horner(self.den, x)

    def self_check(self, deep: bool = False) -> list:
        """Самопроверка дроби из ℚ(x) независимой дорогой (очередь п. 18).

        `[вычислено здесь]` Проверяется не тем счётом, каким строилось:
        коэффициенты — пары ℚ (float дверью `Q` не проходит), знаменатель
        не нулевой и унитарен, у нуля знаменатель ровно (1), старшие
        коэффициенты не нулевые, пересчитанный ЗАНОВО `gcd(P, Q)` есть
        единица (несократимость), `is_monomial` сходится со сборкой `mono`,
        хеш согласован с равенством. При `deep` — пять встречных дорог:
        ОБРАТНОЕ УМНОЖЕНИЕ (домножить P и Q на посторонний множитель g и
        сократить снова — обязана выйти та же форма, а построение шло
        делением на gcd, не умножением); обращение P/Q · Q/P = 1; деление
        с остатком P = Q·q + r; степень двоичной дорогой против линейной;
        знак по сертификату сетки (`sign_at`) против знака отрезка значения
        (`value_at`).

        Граница: тождество ЗНАЧЕНИЙ над РАЦИОНАЛЬНОЙ базой здесь не
        решается — там поле другое, ℚ[x]/(x^N−α) (`alpha.py`); равенство
        форм есть равенство чисел только над трансцендентной базой.
        `reduce_index` сверяется как канон ЗАПИСИ (идемпотентность), а не
        как перенос величины на другое колесо. Знак сверяется в трёх
        рациональных точках, а не на всей прямой: точка, где знаменатель
        обнулился, пропускается названно — полюс не знак.
        """
        out = []
        for name, p in (("числитель", self.num), ("знаменатель", self.den)):
            for i, c in enumerate(p):
                if not isinstance(c, Q):
                    out.append(f"{name}: коэффициент {i} не пара ℚ, а {type(c).__name__}")
            if p and p[-1] == 0:
                out.append(f"{name}: старший коэффициент нулевой — канон записи нарушен")
        if out:
            return out                      # дальше считать нечем: запись уже не запись
        if not self.den:
            out.append("знаменатель — нулевой многочлен: дроби нет")
            return out
        if self.den[-1] != 1:
            out.append(f"знаменатель не унитарен: старший коэффициент {self.den[-1]}")
        if self.is_zero and self.den != (Q(1),):
            out.append(f"у нуля знаменатель {_show_poly(self.den)}, а канон — 1")
        if self.num:
            g = _gcd(self.num, self.den)
            if len(g) > 1:
                out.append(f"дробь сократима: gcd(P, Q) = {_show_poly(g)}")
        mono = self.is_monomial
        if mono is not None:
            c, k = mono
            if RatFunc.mono(k, c) != self:
                out.append(f"is_monomial назвал {c}·x^{k}, а сборка mono даёт иную форму")
        try:
            if hash(self) != hash(RatFunc(self.num, self.den)) \
                    or self != RatFunc(self.num, self.den):
                out.append("хеш либо равенство не согласованы с пересборкой из своих же частей")
        except Exception as e:                          # noqa: BLE001
            out.append(f"пересборка из своих же частей отказала: {type(e).__name__}: {e}")
        if not deep:
            return out
        for g in ((Q(2), Q(-3), Q(7)), (Q(1), Q(1)), (Q(5),)):
            try:
                back = RatFunc(_mul(self.num, g), _mul(self.den, g))
            except Exception as e:                      # noqa: BLE001
                out.append(f"обратное умножение на {_show_poly(g)} отказало: {type(e).__name__}")
                continue
            if back != self:
                out.append(f"обратное умножение на {_show_poly(g)} дало иную форму: {back.record()}")
        if not self.is_zero:
            if self * self.inverse() != RatFunc.const(1):
                out.append(f"P/Q · Q/P не даёт единицы: {(self * self.inverse()).record()}")
            q, r = _divmod(self.num, self.den)
            if _add(_mul(self.den, q), r) != self.num:
                out.append("деление с остатком: Q·q + r ≠ P")
            lin = RatFunc.const(1)
            for _ in range(3):
                lin = lin * self
            if self ** 3 != lin:
                out.append(f"степень двоичной дорогой {(self ** 3).record()} "
                           f"не сходится с линейной {lin.record()}")
        red, n_red = self.reduce_index(12)
        again, n2 = red.reduce_index(n_red)
        if again != red or n2 != n_red:
            out.append(f"reduce_index не идемпотентен: {red.record()} (N={n_red}) → "
                       f"{again.record()} (N={n2})")
        for t in (Q(1, 2), Q(3), Q(-2)):
            iv = Interval.point(t)
            try:
                val = self.value_at(iv)
            except Exception:                           # noqa: BLE001 — полюс: точка не годна
                continue
            by_val = 0 if self.is_zero else (1 if val.lo > 0 else (-1 if val.hi < 0 else None))
            by_cert = self.sign_at(iv)
            if by_cert is not None and by_val is not None and by_cert != by_val:
                out.append(f"знак по сертификату {by_cert} не сходится со знаком "
                           f"значения {by_val} в точке {t}")
        return out
