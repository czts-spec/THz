"""Значение в разложенном виде: ± ∏ p^(a_p) — эталон этапа 3 на Python.

`docs/PLAN.md`, этап 3: «значение — разреженный массив пар (простое,
показатель) с фиксированной ёмкостью. Умножение — слияние двух
отсортированных массивов; деление — то же со знаком; корень степени n —
деление показателей». Здесь та же величина написана на Python и служит
ЭТАЛОНОМ, с которым С-ядро обязано совпасть побитово.

Почему разложенный вид, а не пара целых. У пары (n, d) умножение стоит
умножения длинных целых, а корень не берётся вовсе. У разложенного вида
умножение есть сложение показателей при слиянии двух отсортированных
списков, деление — вычитание, корень степени n — деление показателей на
n, и все три не выходят за целые. Ровно это и синтезируемо (INF-3): в
кремний идёт слияние, а не арифметика длинных чисел.

Показатели рациональны и хранятся ОБЩИМ знаменателем: значение есть
(sign, den, [(p, num_p)]), показатель p равен num_p/den. Общий
знаменатель — не экономия, а условие: корень степени n умножает ровно
одно поле, а не каждый показатель по отдельности.

Канон: члены отсортированы по простому, нулевые выброшены, знаменатель
сокращён с НОД всех числителей. Отсюда равенство есть СОВПАДЕНИЕ ЗАПИСИ,
а не вызов упрощателя, — то самое требование этапа 1.

Отказы поимённо: чётный корень из отрицательного (вещественного нет —
пара (re, im) живёт на линейках 1 и 3), перевод в пару целых при
дробном показателе, разложение числа, не поддавшегося пробному делению.

Авторство: постановка (значение как вектор показателей, отказ вместо
подмены) — AIAΣA (автор проекта); разложение на простые — школьное
пробное деление.
"""
from __future__ import annotations

from math import gcd

from .errors import ExactRepresentationNotAvailable, Unresolved
from .rational import Q

TRIAL_LIMIT = 1 << 20          # до какого делителя идёт пробное деление
PRIME_MAX = (1 << 31) - 1      # простое обязано помещаться в int32 (этап 3)


def _factor(n: int) -> dict:
    """Разложение положительного целого; отказ, если не разложилось."""
    out: dict = {}
    if n == 1:
        return out
    d = 2
    while d * d <= n and d <= TRIAL_LIMIT:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        if n > PRIME_MAX:
            raise Unresolved("разложение", str(n),
                             hint=f"остаток {n} не разложен пробным делением до "
                                  f"{TRIAL_LIMIT} и не помещается в int32")
        out[n] = out.get(n, 0) + 1
    return out


class Factored:
    """± ∏ p^(num_p/den). Канон при построении, равенство — совпадением."""

    __slots__ = ("sign", "den", "terms")

    def __init__(self, sign: int = 1, den: int = 1, terms=None):
        if sign not in (1, -1):
            raise ValueError(f"знак — 1 или −1, получено {sign!r}")
        if not isinstance(den, int) or isinstance(den, bool) or den < 1:
            raise ValueError(f"общий знаменатель — целое ≥ 1, получено {den!r}")
        items = {}
        for p, a in dict(terms or {}).items():
            p, a = int(p), int(a)
            if a == 0:
                continue
            if p < 2 or p > PRIME_MAX:
                raise ValueError(f"простое вне int32: {p}")
            items[p] = items.get(p, 0) + a
        items = {p: a for p, a in items.items() if a != 0}
        g = den
        for a in items.values():
            g = gcd(g, abs(a))
        if g > 1:
            den //= g
            items = {p: a // g for p, a in items.items()}
        object.__setattr__(self, "sign", sign)
        object.__setattr__(self, "den", den)
        object.__setattr__(self, "terms", tuple(sorted(items.items())))

    def __setattr__(self, *_):
        raise AttributeError("Factored неизменяем")

    # --- построение ----------------------------------------------------
    @classmethod
    def from_rational(cls, q) -> "Factored":
        q = Q(q)
        if q.n == 0:
            raise Unresolved("разложение", "0",
                             hint="ноль показателями не выражается: он живёт "
                                  "состоянием с провенансом (zero_state)")
        sign = 1 if q.n > 0 else -1
        terms = _factor(abs(q.n))
        for p, a in _factor(q.d).items():
            terms[p] = terms.get(p, 0) - a
        return cls(sign, 1, terms)

    @classmethod
    def from_int(cls, k: int) -> "Factored":
        return cls.from_rational(Q(k))

    # --- канон и равенство ---------------------------------------------
    def key(self) -> tuple:
        """Каноническая запись целиком — то, что сравнивается побитово."""
        return (self.sign, self.den, self.terms)

    def __eq__(self, other) -> bool:
        return isinstance(other, Factored) and self.key() == other.key()

    def __hash__(self) -> int:
        return hash(self.key())

    @property
    def support(self) -> int:
        return len(self.terms)

    # --- действия -------------------------------------------------------
    def _align(self, other: "Factored") -> tuple:
        """Привести к общему знаменателю показателей — НОК двух."""
        d = self.den * other.den // gcd(self.den, other.den)
        ka, kb = d // self.den, d // other.den
        a = {p: n * ka for p, n in self.terms}
        b = {p: n * kb for p, n in other.terms}
        return d, a, b

    def __mul__(self, other: "Factored") -> "Factored":
        d, a, b = self._align(other)
        for p, n in b.items():
            a[p] = a.get(p, 0) + n
        return Factored(self.sign * other.sign, d, a)

    def __truediv__(self, other: "Factored") -> "Factored":
        d, a, b = self._align(other)
        for p, n in b.items():
            a[p] = a.get(p, 0) - n
        return Factored(self.sign * other.sign, d, a)

    def root(self, k: int) -> "Factored":
        """Корень степени k — деление показателей: знаменатель ×k."""
        if not isinstance(k, int) or isinstance(k, bool) or k < 1:
            raise ValueError(f"степень корня — целое ≥ 1, получено {k!r}")
        if self.sign < 0 and k % 2 == 0:
            raise Unresolved(f"√{k}", self.record(),
                             hint=f"чётная степень корня из отрицательного: "
                                  f"вещественного корня нет, пара (re, im) — "
                                  f"на линейках 1 и 3")
        return Factored(self.sign, self.den * k, dict(self.terms))

    def power(self, k: int) -> "Factored":
        """Целая степень — умножение показателей."""
        if not isinstance(k, int) or isinstance(k, bool):
            raise TypeError("степень — целое")
        sign = 1 if (self.sign > 0 or k % 2 == 0) else -1
        return Factored(sign, self.den, {p: n * k for p, n in self.terms})

    def inverse(self) -> "Factored":
        return Factored(self.sign, self.den, {p: -n for p, n in self.terms})

    # --- переводы -------------------------------------------------------
    def self_check(self, deep: bool = False) -> list:
        """Самопроверка разложенного вида ВСТРЕЧНОЙ дорогой.

        `[вычислено здесь]` Дороги: (1) канон записи пересчитывается здесь
        заново — простые строго возрастают, нулевых показателей нет, общий
        знаменатель ≥ 1 и взаимно прост со всеми показателями (НОД считается
        тут, а не берётся у конструктора), простые лежат в объявленном
        пределе; (2) ОБРАТНАЯ дорога к величине: при целых показателях
        `to_rational()` собирает пару перемножением, а `from_rational()`
        раскладывает её ЗАНОВО пробным делением (`_factor`), и запись
        обязана вернуться той же — то есть разложение и сборка сверяются
        друг против друга, а не каждое само с собой; (3) контракт eq/hash
        по канонической записи.

        Граница объявлена: при ДРОБНЫХ показателях (общий знаменатель > 1)
        пары целых не существует по правилу самого модуля, и дорога (2) не
        идёт — её молчание находкой не считается; корень из отрицательного
        и неразложенный остаток отбиты построением и сюда не доходят.
        `deep` добавляет тождества действий: двойное обращение и
        `power(k).root(k)` возвращают ту же запись — последнее ТОЛЬКО при
        чётном k и знаке +1: у отрицательной записи √(x²) = |x|, и
        расхождение там есть факт математики, а не поломка устройства
        (граница найдена этой же самопроверкой 12.09 на записи −45/8).
        Исключений не бросает.
        """
        out = []
        if self.sign not in (1, -1):
            out.append(f"знак не ±1: {self.sign!r}")
        if not isinstance(self.den, int) or self.den < 1:
            out.append(f"общий знаменатель не целое ≥ 1: {self.den!r}")
        if out:
            return out
        primes = [p for p, _ in self.terms]
        if primes != sorted(primes) or len(set(primes)) != len(primes):
            out.append(f"простые не строго возрастают: {primes}")
        for prime, a in self.terms:
            if a == 0:
                out.append(f"нулевой показатель у {prime} остался в записи")
            if prime < 2 or prime > PRIME_MAX:
                out.append(f"простое вне предела: {prime}")
        g = self.den
        for _, a in self.terms:
            g = gcd(g, abs(a))
        if self.terms and g != 1:
            out.append(f"дробь показателей не сокращена: общий делитель {g}")
        if not self.terms and self.den != 1:
            out.append(f"пустая запись с знаменателем {self.den}")
        if out:
            return out
        if self.den == 1:
            try:
                q = self.to_rational()
                back = Factored.from_rational(q)
                if back != self:
                    out.append(f"разложение и сборка расходятся: {self.record()} "
                               f"→ {q} → {back.record()}")
            except Unresolved:
                if self.terms:
                    out.append("величина есть, а обратное разложение отказало")
            except Exception as e:                              # noqa: BLE001
                out.append(f"обратная дорога отказала: {type(e).__name__}: {e}")
        twin = Factored(self.sign, self.den, dict(self.terms))
        if twin != self or hash(twin) != hash(self):
            out.append("контракт eq/hash нарушен: равные записи — разные ключи")
        if not deep:
            return out
        if self.inverse().inverse() != self:
            out.append("двойное обращение не вернуло запись")
        for k in (2, 3):
            if k % 2 == 0 and self.sign < 0:
                continue          # √(x²) = |x|: тождество неверно математически
            if self.power(k).root(k) != self:
                out.append(f"power({k}).root({k}) не вернул запись")
        return out

    def to_rational(self) -> Q:
        """Пара целых — только при целых показателях; иначе отказ."""
        if self.den != 1:
            raise ExactRepresentationNotAvailable(
                f"{self.record()}: показатели дробные (знаменатель {self.den}) — "
                f"пары целых нет; перевод по явному запросу, не молча")
        n, d = 1, 1
        for p, a in self.terms:
            if a > 0:
                n *= p ** a
            else:
                d *= p ** (-a)
        return Q(self.sign * n, d)

    def record(self) -> str:
        """Явная запись: знак и члены, показатель дробью, если он дробный."""
        if not self.terms:
            body = "1"
        else:
            body = "·".join(
                f"{p}^{n}" if self.den == 1 else f"{p}^({n}/{self.den})"
                for p, n in self.terms)
        return ("-" if self.sign < 0 else "") + body

    def __repr__(self) -> str:
        return f"Factored({self.record()})"


def fits(v: Factored, profile) -> bool:
    return not breaches(v, profile)


def breaches(v: Factored, profile) -> list:
    """Что выходит за профиль: (имя, замерено, предел). Пусто — влезает."""
    out = []
    if v.support > profile.support:
        out.append(("носитель", v.support, profile.support))
    if v.den > profile.exp_den:
        out.append(("знаменатель показателя", v.den, profile.exp_den))
    for p, n in v.terms:
        if abs(n).bit_length() > profile.mantissa_bits:
            out.append((f"показатель при {p}, бит", abs(n).bit_length(),
                        profile.mantissa_bits))
    return out
