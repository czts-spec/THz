"""Q — собственный носитель класса 1: пара целых (n, d), d > 0.

Слово автора (02.09): total_nonFraction — тотально, без исключений.
Чужой универсальный тип снят; вместо него — своя пара целых с
ЗАКРЫТЫМИ подвохами, которые чужой тип нёс скрыто:

* дробная степень у чужого типа тихо давала float — здесь она ТОЧНАЯ:
  полный корень возвращается парой (4^(1/2) → 2), прочее — переходом
  класса на линейку 3 (2^(1/2) → √2 с аннулятором); отказ
  `ExactRepresentationNotAvailable` остаётся только там, где
  вещественного значения нет (чётный корень из отрицательной пары);
* перевода в float нет вовсе (__float__ не определён): десятичная
  сторона — только сертификатом сетки по явному запросу;
* float/bool на входе — TypeError; единственная дверь из мира IEEE —
  явный `Q.from_float` (битово точный перевод, для замеров стороны Б).

Канон: d > 0, gcd(n, d) = 1. Хэш совместим с int (Q(5) и 5 — один
хэш), чтобы целые вырождения жили в множествах рядом с int.
"""
from __future__ import annotations

import sys
from math import gcd

from .errors import ExactRepresentationNotAvailable

_HASH_MOD = sys.hash_info.modulus
_HASH_INF = sys.hash_info.inf


class Q:
    """Точное рациональное: пара целых. Носитель класса 1."""

    __slots__ = ("n", "d")

    def __init__(self, n=0, d=None):
        if isinstance(n, bool) or isinstance(d, bool):
            raise TypeError("bool не число")
        if isinstance(n, float) or isinstance(d, float):
            raise TypeError("float запрещён построением; явная дверь — "
                            "Q.from_float")
        if d is None:
            if isinstance(n, Q):
                object.__setattr__(self, "n", n.n)
                object.__setattr__(self, "d", n.d)
                return
            if isinstance(n, int):
                object.__setattr__(self, "n", n)
                object.__setattr__(self, "d", 1)
                return
            if isinstance(n, str):
                num, den = _parse(n)
            else:
                raise TypeError(f"Q не строится из {type(n).__name__}")
        else:
            if not isinstance(n, int) or not isinstance(d, int):
                if isinstance(n, Q) or isinstance(d, Q):
                    nq, dq = Q(n), Q(d)
                    num = nq.n * dq.d
                    den = nq.d * dq.n
                else:
                    raise TypeError("пара Q(n, d) — только целые или Q")
            else:
                num, den = n, d
        if den == 0:
            raise ZeroDivisionError("знаменатель 0")
        if den < 0:
            num, den = -num, -den
        g = gcd(num, den)
        if g > 1:
            num //= g
            den //= g
        object.__setattr__(self, "n", num)
        object.__setattr__(self, "d", den)

    def __setattr__(self, *_):
        raise AttributeError("Q неизменяем")

    # --- совместимость по именам компонентов ---------------------------
    @property
    def numerator(self) -> int:
        return self.n

    @property
    def denominator(self) -> int:
        return self.d

    @classmethod
    def _canonical(cls, num: int, den: int) -> "Q":
        """Пара, канон которой ДОКАЗАН вызывающей стороной (den > 0,
        gcd(num, den) = 1) — без повторного gcd.

        Замер 12.09 на целой степени: gcd по уже канонической паре стоит
        7,6× самого возведения ((3/7)^20000: 3693 мкс против 483 мкс).
        Дверь закрытая: имя с подчерком, единственный вызов — целая
        степень, где канон сохраняется тождеством
        gcd(n, d) = 1 ⇒ gcd(n^p, d^p) = 1."""
        r = cls.__new__(cls)
        object.__setattr__(r, "n", num)
        object.__setattr__(r, "d", den)
        return r

    @classmethod
    def from_float(cls, x: float) -> "Q":
        """Явный битово точный перевод IEEE 754 → пара целых.

        Единственная дверь из float; используется замерами стороны Б —
        скрытых переводов в библиотеке нет."""
        if not isinstance(x, float):
            raise TypeError("from_float принимает только float")
        num, den = x.as_integer_ratio()
        return cls(num, den)

    # --- арифметика (int и Q; прочее — NotImplemented) ------------------
    def _coerce(self, other):
        if isinstance(other, Q):
            return other
        if isinstance(other, int) and not isinstance(other, bool):
            return Q(other)
        return None

    def __add__(self, other):
        o = self._coerce(other)
        if o is None:
            return NotImplemented
        return Q(self.n * o.d + o.n * self.d, self.d * o.d)

    __radd__ = __add__

    def __sub__(self, other):
        o = self._coerce(other)
        if o is None:
            return NotImplemented
        return Q(self.n * o.d - o.n * self.d, self.d * o.d)

    def __rsub__(self, other):
        o = self._coerce(other)
        if o is None:
            return NotImplemented
        return Q(o.n * self.d - self.n * o.d, self.d * o.d)

    def __mul__(self, other):
        o = self._coerce(other)
        if o is None:
            return NotImplemented
        return Q(self.n * o.n, self.d * o.d)

    __rmul__ = __mul__

    def __truediv__(self, other):
        o = self._coerce(other)
        if o is None:
            return NotImplemented
        if o.n == 0:
            raise ZeroDivisionError("деление на нулевую пару")
        return Q(self.n * o.d, self.d * o.n)

    def __rtruediv__(self, other):
        o = self._coerce(other)
        if o is None:
            return NotImplemented
        if self.n == 0:
            raise ZeroDivisionError("деление на нулевую пару")
        return Q(o.n * self.d, o.d * self.n)

    def __pow__(self, p):
        """Степень точная всегда (слово автора 02.09: «дробная — почти
        вся математика», отказ не одобрен). Целый показатель — пара;
        дробный — ТОЧНЫЙ переход на линейку 3: полный корень
        возвращается парой (4^(1/2) → 2), прочее — алгебраической
        величиной с аннулятором. Чужой тип здесь тихо давал float —
        мы даём точное значение на несущей линейке."""
        if isinstance(p, bool):
            return NotImplemented
        if isinstance(p, int):
            if p >= 0:
                return self._int_pow(p)
            if self.n == 0:
                raise ZeroDivisionError("0 в отрицательной степени")
            return Q(self.d, self.n)._int_pow(-p)
        if isinstance(p, float):
            raise TypeError("float-показатель запрещён построением")
        if not isinstance(p, Q):
            return NotImplemented
        num, den = p.n, p.d
        if self.n == 0:
            if num > 0:
                return Q(0)
            raise ZeroDivisionError("0 в отрицательной степени")
        if den == 1:
            # показатель ЦЕЛЫЙ, только в одежде пары: целая дорога, без
            # корня и без сведения на линейке 3. До 12.09 такой показатель
            # всё равно уходил алгебраической дорогой, и Q(10)**Q(8) стоило
            # 39 972 мс на линейном скане (algebraic.py:_try_rational);
            # тот же случай целой дорогой — 0,013 мс, отношение 3·10⁶.
            # Ветка стоит НИЖЕ разбора нулевого основания: 0^Q(0) остаётся
            # отказом (ZeroDivisionError), тогда как 0**0 дало бы 1.
            return self ** num
        base = self if num >= 0 else Q(self.d, self.n)
        radicand = base ** abs(num)          # целая степень — пара
        # ленивый импорт: линейка 3 стоит выше носителя в сборке пакета
        from .algebraic import AlgebraicNumber
        if radicand < 0:
            if den % 2 == 0:
                # вещественного значения нет — переход на пару (re, im)
                # (очередь 02.09: «представление через линейку 4 — в очереди»;
                # 04.09 — пара: |r|^(1/q) под углом π/q, главное значение, точно)
                from .complex_pair import ComplexPair
                return ComplexPair.from_real_root(radicand, den)
            return AlgebraicNumber.from_root(den, -radicand, sign=1)
        return AlgebraicNumber.from_root(den, radicand)

    def _int_pow(self, p: int) -> "Q":
        """Целая неотрицательная степень пары — БЫСТРЫМ возведением в
        квадрат, не линейным умножением.

        Решение принято замером 12.09 (четыре дороги, значения совпали
        побитово на всех шести случаях):

        | дорога                        | (3/7)^20000 |
        |-------------------------------|-------------|
        | линейное умножение (p шагов)  | 23 240 141 мкс |
        | квадрат-и-умножение на Python |        823 мкс |
        | int ** int, канон через gcd   |      3 693 мкс |
        | int ** int, канон без gcd     |        483 мкс |

        Линейная дорога проигрывает 48 000×. Квадрат-и-умножение,
        выписанное уровнем Python, проигрывает 1,7× тому же алгоритму
        внутри CPython: `int.__pow__` ЕСТЬ двоичное возведение (Handbook
        of Applied Cryptography, алг. 14.76), исполненное на C. Потому
        квадрат-и-умножение здесь не переписывается, а берётся у целых,
        и снимается единственная лишняя цена — повторный gcd (7,6×).
        Точность: ни одного float, ни одного деления; канон сохраняется
        тождеством gcd(n, d) = 1 ⇒ gcd(n^p, d^p) = 1."""
        if p < 0:
            raise ValueError("_int_pow: показатель отрицателен")
        return Q._canonical(self.n ** p, self.d ** p)

    def __rpow__(self, base):
        """основание ** пара: то же точное сведение (2 ** 1/2 = √2)."""
        if isinstance(base, bool) or isinstance(base, float):
            raise TypeError("основание — только int или Q")
        if isinstance(base, int):
            return Q(base) ** self
        return NotImplemented

    def __floordiv__(self, other):
        o = self._coerce(other)
        if o is None:
            return NotImplemented
        if o.n == 0:
            raise ZeroDivisionError("деление на нулевую пару")
        return (self.n * o.d) // (o.n * self.d)

    def __rfloordiv__(self, other):
        o = self._coerce(other)
        if o is None:
            return NotImplemented
        if self.n == 0:
            raise ZeroDivisionError("деление на нулевую пару")
        return (o.n * self.d) // (self.n * o.d)

    def __mod__(self, other):
        o = self._coerce(other)
        if o is None:
            return NotImplemented
        return self - o * (self // o)

    def __rmod__(self, other):
        o = self._coerce(other)
        if o is None:
            return NotImplemented
        return o - self * (o // self)

    def __divmod__(self, other):
        o = self._coerce(other)
        if o is None:
            return NotImplemented
        q = self // o
        return q, self - o * q

    def __rdivmod__(self, other):
        o = self._coerce(other)
        if o is None:
            return NotImplemented
        q = o // self
        return q, o - self * q

    def __neg__(self):
        return Q(-self.n, self.d)

    def __pos__(self):
        return self

    def __abs__(self):
        return Q(-self.n, self.d) if self.n < 0 else self

    # --- сравнения -------------------------------------------------------
    def _cmp(self, other):
        o = self._coerce(other)
        if o is None:
            return None
        return self.n * o.d - o.n * self.d

    def __eq__(self, other):
        c = self._cmp(other)
        return NotImplemented if c is None else c == 0

    def __lt__(self, other):
        c = self._cmp(other)
        return NotImplemented if c is None else c < 0

    def __le__(self, other):
        c = self._cmp(other)
        return NotImplemented if c is None else c <= 0

    def __gt__(self, other):
        c = self._cmp(other)
        return NotImplemented if c is None else c > 0

    def __ge__(self, other):
        c = self._cmp(other)
        return NotImplemented if c is None else c >= 0

    def __bool__(self):
        return self.n != 0

    def __hash__(self):
        # числовой хэш, совместимый с int: Q(5) и 5 — один хэш
        try:
            dinv = pow(self.d, -1, _HASH_MOD)
        except ValueError:
            h = _HASH_INF
        else:
            h = (abs(self.n) * dinv) % _HASH_MOD
        h = -h if self.n < 0 else h
        return -2 if h == -1 else h

    # --- вырождения и запись --------------------------------------------
    def __int__(self):
        """Усечение к нулю; точно, когда d == 1."""
        return -((-self.n) // self.d) if self.n < 0 else self.n // self.d

    __trunc__ = __int__

    def self_check(self, deep: bool = False) -> list:
        """Самопроверка пары ВСТРЕЧНОЙ дорогой (ось «AIAΣA спрашивает AIAΣA»).

        `[вычислено здесь]` Канон пары (d > 0, gcd(|n|, d) = 1) проверяется не
        тем счётом, каким строился: НОД считается здесь заново, а сама пара
        восстанавливается ЦЕПНОЙ ДРОБЬЮ — независимым алгоритмом (частные
        деления с остатком, свёртка подходящих дробей): у несократимой пары
        свёртка обязана дать её же. Порядок сверяется двумя дорогами —
        знаком разности и перекрёстным умножением. Контракт eq/hash с целым
        (Q(5) и 5 — один хэш) проверяется на самой паре.

        Граница: самопроверка говорит о ЗАПИСИ уже построенной пары, не о
        входе конструктора — знаменатель 0, float и bool отбиты построением
        (ZeroDivisionError и TypeError), сюда такая пара не доходит. У пары с
        огромными членами цепная дробь стоит O(log) делений — это и есть
        цена `deep=False`; `deep` добавляет тождества арифметики
        (x + (−x) = 0, x² = x·x, 1/x = x^−1), каждое — другой дорогой.
        Возвращает список находок; пустой список — молчание.
        """
        out = []
        if not isinstance(self.d, int) or self.d <= 0:
            out.append(f"знаменатель не положителен: {self.d!r}")
            return out
        if not isinstance(self.n, int):
            out.append(f"числитель не целое: {self.n!r}")
            return out
        g = gcd(abs(self.n), self.d)
        if g != 1:
            out.append(f"пара не сокращена: gcd(|{self.n}|, {self.d}) = {g}")
        # встречная дорога: цепная дробь и свёртка обратно
        a, b = self.n, self.d
        parts = []
        while b:
            q = a // b
            parts.append(q)
            a, b = b, a - q * b
        num, den = 1, 0
        for q in reversed(parts):
            num, den = q * num + den, num
        if den < 0:
            num, den = -num, -den
        if (num, den) != (self.n, self.d):
            out.append(f"свёртка цепной дроби дала {num}/{den}, а записано "
                       f"{self.n}/{self.d}")
        if self.d == 1 and hash(self) != hash(self.n):
            out.append("контракт eq/hash с целым нарушен: Q(n) и n — разные хэши")
        for other in (Q(self.n + 1, self.d), Q(self.n - 1, self.d)):
            cross = self.n * other.d - other.n * self.d
            byorder = (self > other) - (self < other)
            if (cross > 0) - (cross < 0) != byorder:
                out.append(f"порядок против перекрёстного умножения расходится "
                           f"на {other!r}")
        if not deep:
            return out
        if self + (-self) != Q(0):
            out.append("x + (−x) ≠ 0")
        if self ** 2 != self * self:
            out.append("x² ≠ x·x")
        if self.n != 0 and Q(1) / self != self ** -1:
            out.append("1/x ≠ x^−1")
        return out

    def record(self) -> str:
        return f"Φ(C=1 RATIONAL, value={self})"

    def __repr__(self):
        return str(self.n) if self.d == 1 else f"{self.n}/{self.d}"

    __str__ = __repr__


def nearest_int(q: "Q") -> int:
    """Ближайшее целое к паре; ровно половина — ОТ НУЛЯ.

    `[вычислено здесь]` Правило названо, а не унаследовано. Встроенный
    `round()` округляет половину к чётному (банковское), и это разумно для
    статистики и неразумно там, где округление входит в доказательство
    включения: правило должно быть одно и записанное. Здесь оно одно.

    Работает на целых: ни float, ни Fraction не участвуют.
    """
    n, d = q.n, q.d
    if n < 0:
        return -nearest_int(Q(-n, d))
    whole, rest = divmod(n, d)
    return whole + 1 if 2 * rest >= d else whole


def limit_denominator(q: "Q", maxden: int) -> "Q":
    """Наилучшее приближение пары при знаменателе ≤ maxden — цепной дробью.

    `[вычислено здесь]` Тот же ответ, что у `Fraction.limit_denominator`,
    но без захода в стандартную дробь: дверь для float у неё открыта
    (`Fraction(0.1)` проходит молча), а у `Q` закрыта построением, и
    открывать её ради одного алгоритма незачем.
    """
    if not isinstance(maxden, int) or isinstance(maxden, bool) or maxden < 1:
        raise ValueError(f"предел знаменателя — целое ≥ 1, получено {maxden!r}")
    if q.d <= maxden:
        return q
    p0, q0, p1, q1 = 0, 1, 1, 0
    n, d = q.n, q.d
    while True:
        a = n // d
        q2 = q0 + a * q1
        if q2 > maxden:
            break
        p0, q0, p1, q1 = p1, q1, p0 + a * p1, q2
        n, d = d, n - a * d
        if d == 0:
            break
    k = (maxden - q0) // q1
    lo, hi = Q(p0 + k * p1, q0 + k * q1), Q(p1, q1)
    return hi if abs(hi - q) <= abs(lo - q) else lo


def _parse(text: str, base: int = 10):
    """Разбор "a/b", целой и позиционной записи — целыми, без float.

    Позиционная запись есть СОКРАЩЁННАЯ ЗАПИСЬ СУММЫ КОЛЁС (слово автора
    09.09): цифра на месте k — мантисса колеса номер k, отсутствие цифры —
    отсутствие колеса. Потому её разбирает КОЛЕСО (`wheels.Record.parse`), а
    не две целочисленные строки с зашитой десяткой; основание — параметр.

    Ввоз ленивый по необходимости: `wheels` ввозит `rational`, и ввоз верхнего
    уровня замкнул бы круг загрузки.

    Замер 22.09 `[вычислено здесь]`: за наборами `test_phicore`,
    `test_phicore_rational` и `test_phicore_container` эта дверь звалась 6 раз,
    из них 2 позиционных, — путь не горячий, и цена колеса (в 57 раз дороже
    прежних двух действий на разбор) в прогоне не видна.
    """
    from .wheels import _DIGITS, Record as _WheelRecord   # круг загрузки: ленивый
    t = text.strip()
    if "/" in t:
        a, b = t.split("/", 1)
        return int(a.strip()), int(b.strip())
    if "." in t:
        тело = t.lstrip("+-").lower()
        whole, frac = тело.split(".", 1)
        if not (whole + frac) or not set(whole + frac) <= set(_DIGITS[:base]):
            raise ValueError(f"не пара: {text!r}")
        v = _WheelRecord.parse(t, base, sep=".").value
        return v.n, v.d
    return int(t), 1
