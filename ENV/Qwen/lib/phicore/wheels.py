"""Колёса: позиционная запись как сумма Ф-контейнеров; колесо-ранг W(x; сетка).

Слово автора 10.09 `[в корпусе]` (коррекция 583e1146dd): «колесо на 1–9 —
НЕ ДОГМА Ф-чисел. Ф-число БЕЗРАЗЛИЧНО к размеру колеса, НЕТ предпочтений».
Слово автора 09.09 `[в корпусе]` (коррекции 439054471a … 06f1d4402e):
«десятичное число представлено колёсами 1–9; 0 — отсутствие колеса;
основная сетка — разряд единиц; класс 22 несёт номер целого колеса, класс
21 — дробного; каждое колесо — полноценный контейнер; знак + подразумевается:
541814 = 500000+40000+1000+800+10+4; W_b(x) — инвариант пары (число, сетка)».

Численная форма. Колесо — контейнер (m1 = m, m2 = 1, e1 = e, e2 = 1) над
α = b: значение m·b^e, m ∈ ℤ \ {0} ЛЮБОЙ ширины, e ∈ ℤ. Размер колеса —
не свойство числа, а параметр сетки: сетка (b, w) допускает |m| < b^w
(разряды основания b^w; при w = 1 — цифры 1…b−1 десятичного примера);
сетка без ширины (w = None) — колесо любой величины. Ф-число безразлично к
выбору: значение одно при всех (b, w). Номер колеса относительно
разделителя — по знаку экспоненты: e > 0 — целое колесо k = e + 1
(указатель 22); e = 0 — основная сетка (без указателя); e < 0 — дробное
колесо k = −e (указатель 21). Запись — конечное множество колёс с попарно
разными e; «+» подразумевается; пустая запись — класс 0 (ноль есть
отсутствие колеса, не число).

Алгебра `[вычислено здесь]`. Для сетки (b, w), R = b^w: ненормированные
записи — свободный коммутативный моноид ℕ^(ℤ) (кратности единиц R^e);
оценка ε(f) = Σ f(e)·R^e — гомоморфизм на ℤ[1/b]_{≥0}; перенос
f(e) ≥ R → (f(e) − R, f(e+1) + 1) сохраняет ε, завершается и конфлюэнтен;
нормальные формы f(e) ≤ R − 1 ↔ числа взаимно однозначно:
ℤ[1/b]_{≥0} ≅ ⟨u_e | R·u_e = u_{e+1}⟩ (коммутативный моноид). Без ширины
переноса НЕТ: колёса на одной позиции складываются мантиссами, и всякий
x ∈ ℤ[1/b] есть ОДНО колесо (m = N, e = −s); с контейнером (слот m2) одно
колесо есть всякое рациональное; со слотом e2 — и α^(p/q). Перенос, заём,
«схлопывание» 999 + 1 = 1000 — свойства ГРАНИЦЫ ширины колеса (регистра),
не числа и не колеса.

Закон ранга `[вычислено здесь]`. W(x; b, w) — число колёс нормальной
формы. W < ∞ ⟺ x ∈ ℤ[1/b] — не зависит от w; иначе запись периодична:
предпериод ℓ = min{s : q_b | b^s}, период t = ord_{q'}(b) в разрядах
основания b; в разрядах b^w — ⌈ℓ/w⌉ и t/gcd(t, w) (порядок b^w).
Монотонность по ширине: W(x; b, kw) ≤ W(x; b, w) (укрупнение колеса не
увеличивает счёт: ненулевой широкий разряд содержит ненулевой узкий);
W(x; b, None) = 1 для всех x ∈ ℤ[1/b] \ {0}. Субаддитивность
W(x + y) ≤ W(x) + W(y) и W(x·y) ≤ 2·W(x)·W(y) — при любой ширине
(доказательство инъекцией позиций; произведение двух разрядов < R²).
Число ненормированных записей целого n при ширине w — b^w-ичные
разбиения p_R(n) = p_R(n − 1) + [R | n]·p_R(n/R) (b = 2, w = 1 — OEIS
A018819; b = 3 — A062051; асимптотика — Малер 1940, де Брёйн 1948)
`[внешний источник]`.

Цена `[вычислено здесь]`. Колесо стоит bits(|m|) + бит индекса + бит знака
/зоны; запись — сумма по колёсам плюс счёт (uniform_cost). Раз размер колеса
свободен, нарезка числа на колёса — задача о цене, не о числе: optimal()
решает её точно динамикой по ненулевым разрядам (колесо покрывает отрезок
разрядов, внутренние нули стоят log₂ b бит каждый; разрыв нулей режется,
когда он дороже индекса). Слово автора 10.09 `[в корпусе]`: «Ф-число
никогда не хуже своего собственного частного случая» — здесь это теорема:
всякая сетка (b, w) есть одна из нарезок (колесо ширины w покрывает блок из
w разрядов), потому uniform_cost(optimal(x)) ≤ uniform_cost(запись x в
сетке (b, w)) при всех w, и W(x; free) = 1 ≤ W(x; b, w). Пример автора
5000037,000000800000000000000001: пять цифр 1–9 — 53 бита при полях
регистра (4 бита цифра); лучшая нарезка — 4 колеса (5·10⁶, 37, 8·10⁻⁷,
1·10⁻²⁴), 41 бит; одно колесо — 109. Байт колеса по правилам (WheelByte) —
одна из сеток: ширина 1, позиции ±K, (2K+1)·(b−1) ≤ 256.

Область применения — слово автора 10.09 `[в корпусе]` (коррекция после
поправки о ширине): «десятичный пример НЕ является правилом обязательного
разделения числа на колёса. Аппаратная база та же — самые производительные
блоки. Такое представление — лишь оптимизация там, где int8, float32, AVX,
GPU и другие стандартные блоки переполняются или отказывают. Там, где всё
работает нативно хорошо, действия НЕ ТРЕБУЕТСЯ». Значит: нарезка на колёса —
не обязательный шаг записи, а ответ на переполнение или отказ стандартного
блока (окно window.py показало образ с поглощением, регистр не вмещает
мантиссу); при нативном точном счёте число остаётся одним контейнером.

Граница. Float внутри нет; знаменатель с простым вне b — отказ по имени
(Unresolved) с предпериодом и периодом; байт вне досягаемости — отказ.
Замер 10.09 `[вычислено здесь]` (прогон «производительность с сокращением
нулевых», ENV/AIAΣA/notes/A): операции над разрядами — столбик S / S²,
колёса W_a + W_b / W_a·W_b (точные счётчики); время на этой машине при
S = 1024: умножение столбиком 185 мс, колёсами при W = 1 — 9 мкс, при
W = 1024 — 172 мс, int Python (стандартный блок) — 13 мкс: выигрыш
сокращения нулевых — (S/W)² против столбика; против стандартного блока,
пока он не переполнен, выигрыша нет — слово автора 10.09 подтверждено
замером. Такты железа — счёт по описаниям архитектур в proto/registers.py.
"""
from __future__ import annotations

from math import gcd
from typing import Dict, Iterable, List, Optional, Tuple

from .container import LockedInt, NumberClass, PhiNumberV6
from .errors import Unresolved
from .rational import Q

ZONE_INT = 22    #: целое колесо k ≥ 2 (e = k − 1 > 0) — указатель класса 22
ZONE_MAIN = 1    #: основная сетка, разряд единиц (e = 0) — контейнер класса 1 без указателя
ZONE_FRAC = 21   #: дробное колесо k ≥ 1 (e = −k < 0) — указатель класса 21

_DIGITS = "0123456789abcdefghijklmnopqrstuvwxyz"


class _Infinite:
    """Бесконечный ранг: колёс потребовалось бы неограниченно много (ω)."""
    __slots__ = ()
    is_finite = False

    def __repr__(self) -> str:
        return "ω"

    def __eq__(self, other) -> bool:
        return isinstance(other, _Infinite)

    def __hash__(self) -> int:
        return hash("ω")


INFINITE = _Infinite()


def _check_int(v, name: str) -> int:
    if isinstance(v, bool) or not isinstance(v, int):
        raise TypeError(f"{name} — только int, получено {type(v).__name__}")
    return v


def _check_base(base) -> int:
    base = _check_int(base, "base")
    if base < 2:
        raise ValueError(f"база {base} < 2")
    return base


def bits_for(n: int) -> int:
    """Бит, чтобы хранить значения 0…n (n ≥ 0): bit_length."""
    return _check_int(n, "n").bit_length()


# ---------------------------------------------------------------------------
# Колесо
# ---------------------------------------------------------------------------

class Wheel:
    """Колесо m·b^e: контейнер (m, 1, e, 1) над α = b; m ∈ ℤ \\ {0} любой ширины, e ∈ ℤ.
    Размер колеса — параметр сетки (Record.width), не догма колеса."""

    __slots__ = ("base", "m", "e")

    def __init__(self, base: int, m: int, e: int):
        base = _check_base(base)
        m = _check_int(m, "m")
        e = _check_int(e, "e")
        if m == 0:
            raise ValueError("колеса 0 не существует: 0 — отсутствие колеса, не мантисса")
        self.base, self.m, self.e = base, m, e

    @property
    def digit(self) -> int:
        """Мантисса колеса (имя «цифра» — от десятичного примера; ширина любая)."""
        return self.m

    @property
    def width(self) -> int:
        """Ширина колеса в разрядах основания b: число цифр |m|."""
        return len(digits_of(self.m, self.base))

    # --- значение и зона ---------------------------------------------------

    @property
    def value(self) -> Q:
        if self.e >= 0:
            return Q(self.m * self.base ** self.e)
        return Q(self.m, self.base ** (-self.e))

    @property
    def zone(self) -> int:
        """22 при e > 0, 21 при e < 0, 1 (основная сетка) при e = 0 — знак экспоненты, без float."""
        return ZONE_INT if self.e > 0 else (ZONE_FRAC if self.e < 0 else ZONE_MAIN)

    @property
    def index(self) -> int:
        """Номер колеса относительно разделителя: целое k = e + 1, дробное k = −e; колеса 0 нет."""
        return self.e + 1 if self.e >= 0 else -self.e

    def tag(self) -> tuple:
        """[класс 22, k, d] | [класс 21, k, d] | (d,) — основная сетка пишется голым контейнером."""
        z = self.zone
        if z == ZONE_MAIN:
            return (self.m,)
        return (z, self.index, self.m)

    @classmethod
    def from_tag(cls, base: int, tag) -> "Wheel":
        tag = tuple(tag)
        if len(tag) == 1:
            return cls(base, tag[0], 0)
        if len(tag) != 3:
            raise ValueError(f"тег колеса — (d,) или (класс, k, d), получено {tag!r}")
        zone, k, d = tag
        k = _check_int(k, "k")
        if k < 1:
            raise ValueError(f"номер колеса {k} < 1: колеса 0 не существует")
        if zone == ZONE_INT:
            return cls(base, d, k - 1) if k > 1 else cls(base, d, 0)
        if zone == ZONE_FRAC:
            return cls(base, d, -k)
        raise ValueError(f"зона {zone!r} — не 21 и не 22")

    # --- контейнер -----------------------------------------------------------

    def fields(self) -> Tuple[int, int, int, int]:
        """Регистры контейнера (m1, m2, e1, e2) над α = b."""
        return (self.m, 1, self.e, 1)

    def container(self) -> PhiNumberV6:
        """Полноценный контейнер класса 1: (d/1)·b^(e/1); при b = 2 — база неявная."""
        if self.base == 2:
            return PhiNumberV6(LockedInt(self.m), LockedInt(1), LockedInt(self.e), LockedInt(1),
                               NumberClass.NATURAL_POWER)
        return PhiNumberV6(LockedInt(self.m), LockedInt(1), LockedInt(self.e), LockedInt(1),
                           NumberClass.NATURAL_POWER, base_class=NumberClass.NATURAL_POWER,
                           base_params={"value": self.base})

    # --- самопроверка (ось «AIAΣA спрашивает AIAΣA», 04.09) --------------------

    def self_check(self, deep: bool = False) -> list:
        """Самопроверка колеса ВСТРЕЧНОЙ дорогой (не пересчётом тем же кодом).

        `[вычислено здесь]` Дороги: (1) величина — вторым счётом: здесь она
        собирается ветвлением по знаку показателя (`m·b^e` либо `m/b^(−e)`),
        встречная дорога умножает пару `Q(m)` на `Q(base) ** e` средствами
        линейки 1, где отрицательная степень — своё правило, и числа обязаны
        совпасть; (2) тег и обратно — `from_tag(base, tag())` обязан вернуть
        то же колесо (запись зоны 21/22 с НОМЕРОМ колеса против пары (m, e):
        номер целого колеса k = e + 1, дробного k = −e); (3) мост в СОСЕДНИЙ
        модуль — `container()`: контейнер (m, 1, e, 1) над α = b при b = 2
        обязан отдать ту же величину точной парой, при b ≠ 2 — назвать
        основание полем `base_params`; (4) зона по ЗНАКУ показателя (22 при
        e > 0, 21 при e < 0, основная сетка при e = 0) против записи.

        Граница названа: при b ≠ 2 точной пары у контейнера нет по его
        собственному правилу (`to_exact_int_ratio` отказывает при заданном
        основании) — там дорога (3) сверяет только основание и регистры, и
        её молчание о величине находкой НЕ считается. Ширина колеса есть
        свойство СЕТКИ (b, w), а не колеса, и здесь не проверяется — её
        держит `Record`. Возвращает список находок; исключений не бросает.
        """
        out = []
        for name in ("base", "m", "e"):
            v = getattr(self, name)
            if not isinstance(v, int) or isinstance(v, bool):
                out.append(f"поле {name} не целое: {v!r}")
        if out:
            return out
        if self.base < 2:
            out.append(f"основание {self.base} < 2")
        if self.m == 0:
            out.append("мантисса 0: колеса 0 не существует — 0 есть "
                       "отсутствие колеса")
        if out:
            return out
        want = Q(self.m) * (Q(self.base) ** self.e)
        if self.value != want:
            out.append(f"величина двумя дорогами расходится: {self.value} "
                       f"против {want}")
        try:
            back = Wheel.from_tag(self.base, self.tag())
            if back != self:
                out.append(f"тег и обратно дал иное колесо: {back!r}")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"разбор своего тега отказал: {type(e).__name__}: {e}")
        zone_want = ZONE_INT if self.e > 0 else (ZONE_FRAC if self.e < 0 else ZONE_MAIN)
        if self.zone != zone_want:
            out.append(f"зона {self.zone} не отвечает показателю {self.e}")
        if self.index < 1:
            out.append(f"номер колеса {self.index} < 1: колеса 0 не существует")
        try:
            c = self.container()
            if (c.m1.v, c.m2.v, c.e1.v, c.e2.v) != self.fields():
                out.append(f"регистры контейнера {(c.m1.v, c.m2.v, c.e1.v, c.e2.v)} "
                           f"не те же, что у колеса {self.fields()}")
            if self.base == 2:
                n, d = c.to_exact_int_ratio()
                if Q(n, d) != self.value:
                    out.append(f"величина через контейнер {Q(n, d)} против "
                               f"величины колеса {self.value}")
            elif c.base_params.get("value") != self.base:
                out.append(f"контейнер назвал основание "
                           f"{c.base_params.get('value')!r}, у колеса {self.base}")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"мост в контейнер отказал: {type(e).__name__}: {e}")
        if not deep:
            return out
        digits = digits_of(self.m, self.base)
        n = abs(self.m)
        count = 0
        while n:
            n //= self.base
            count += 1
        if len(digits) != count:
            out.append(f"ширина {len(digits)} не сходится со счётом делений "
                       f"{count}")
        return out

    # --- служебное -------------------------------------------------------------

    def __eq__(self, other) -> bool:
        return (isinstance(other, Wheel) and self.base == other.base
                and self.m == other.m and self.e == other.e)

    def __hash__(self) -> int:
        return hash((self.base, self.m, self.e))

    def __repr__(self) -> str:
        t = self.tag()
        if len(t) == 1:
            return f"{t[0]}"
        return f"[класс {t[0]}, {t[1]}, {t[2]}]"


# ---------------------------------------------------------------------------
# Знаменатель: часть из простых базы и остаток; предпериод и период
# ---------------------------------------------------------------------------

def top_wheel(x, base: int) -> int:
    """Номер СТАРШЕГО колеса величины x в сетке (base, 1) — закон, общая дорога.

    Для x ∈ ℤ[1/base] это `max(Record.of(x, base, width=1).counts())`; для
    прочих запись бесконечна, но старшее колесо всё равно конечно и находится
    сравнением степеней. Дорога общая и медленная — замер 22.09
    `[вычислено здесь]`: в 41 раз дороже битовой тропы `top_wheel2`.
    """
    x = Q(x)
    if x == 0:
        raise ValueError("у нуля нет старшего колеса: ноль есть отсутствие колеса")
    x = -x if x < 0 else x
    qb, rest, ell = split_denominator(x.d, base)
    if rest == 1:
        return max(Record.of(x, base, width=1).counts())
    e = 0
    сетка = Q(base)
    while сетка ** (e + 1) <= x:
        e += 1
    while сетка ** e > x:
        e -= 1
    return e


def top_wheel2(n: int, d: int) -> int:
    """Тот же номер при base = 2, целыми по длинам двоичной записи — ТРОПА.

    e с 2^e ≤ n/d < 2^(e+1) для n, d > 0. Это НЕ отдельный закон, а быстрый
    путь к `top_wheel(Q(n, d), 2)`; равносильность проверяется переборно
    приёмкой (`proto/test_phicore_wheels.py`, с подлогом сдвига).

    Тропа живёт здесь, в модуле закона, и ОДНА. Замер 22.09
    `[вычислено здесь]`: прежде она была переписана руками дважды —
    `ieee._binade` и `window._floor_log2`, — совпадала на 671 паре, и ничто
    этого не сверяло; обе двери горячие (32 580 и 359 473 вызова за свои
    наборы), потому колесом их не заменить, но и копиям места нет.
    """
    e = n.bit_length() - d.bit_length()
    if (n << max(0, -e)) < (d << max(0, e)):
        e -= 1
    return e


def split_denominator(den: int, base: int) -> Tuple[int, int, int]:
    """(q_b, q', ℓ): q_b — часть den из простых делителей base, q' — остаток
    (взаимно прост с base), ℓ = min{s : q_b | base^s} — предпериод."""
    den = _check_int(den, "den")
    base = _check_base(base)
    if den <= 0:
        raise ValueError("знаменатель должен быть > 0")
    rest = den
    g = gcd(rest, base)
    while g > 1:
        rest //= g
        g = gcd(rest, base)
    qb = den // rest
    ell, power = 0, 1
    while power % qb:
        power *= base
        ell += 1
    return qb, rest, ell


def mult_order(base: int, m: int) -> int:
    """Порядок base по модулю m (gcd = 1, m ≥ 2): наименьшее t ≥ 1 с base^t ≡ 1."""
    base, m = _check_base(base), _check_int(m, "m")
    if m < 2:
        raise ValueError("порядок определён при m ≥ 2")
    if gcd(base, m) != 1:
        raise ValueError(f"gcd({base}, {m}) ≠ 1 — порядка нет")
    t, cur = 1, base % m
    while cur != 1:
        cur = (cur * base) % m
        t += 1
    return t


def period(x, base: int, width: int = 1) -> Tuple[int, int]:
    """(предпериод ℓ, период t) записи x в разрядах основания b^w; t = 0 — запись конечна.
    В разрядах b^w: ℓ_w = ⌈ℓ/w⌉, t_w = ord_{q'}(b^w) = t / gcd(t, w)."""
    x = Q(x)
    base = _check_base(base)
    width = _check_int(width, "width")
    if width < 1:
        raise ValueError("ширина ≥ 1")
    if x == 0:
        return (0, 0)
    _, rest, ell = split_denominator(x.d, base)
    if rest == 1:
        return (-(-ell // width), 0)
    t = mult_order(base, rest)
    return (-(-ell // width), t // gcd(t, width))


def radix_of(base: int, width) -> Optional[int]:
    """Основание разряда сетки (b, w): b^w; None — колесо без ширины."""
    base = _check_base(base)
    if width is None:
        return None
    width = _check_int(width, "width")
    if width < 1:
        raise ValueError("ширина колеса ≥ 1 разряда либо None (без ширины)")
    return base ** width


def wheel_rank(x, base: int, width=None):
    """W(x; b, w): число колёс нормальной формы в сетке (b, w); INFINITE (ω), если
    x ∉ ℤ[1/b] (от ширины не зависит); без ширины — 1 для всякого x ≠ 0 из ℤ[1/b]."""
    x = Q(x)
    R = radix_of(base, width)
    if x == 0:
        return 0
    _, rest, ell = split_denominator(x.d, base)
    if rest != 1:
        return INFINITE
    if R is None:
        return 1
    # позиции считаются в разрядах b^w: сдвиг до кратного w
    shift = (-ell) % width
    n = abs(x.n) * base ** (ell + shift) // x.d
    count = 0
    while n:
        n, r = divmod(n, R)
        count += 1 if r else 0
    return count


def digits_of(n: int, base: int) -> List[int]:
    """Цифры |n| в базе, от младшей к старшей; [] для 0."""
    n = abs(_check_int(n, "n"))
    out: List[int] = []
    while n:
        n, r = divmod(n, base)
        out.append(r)
    return out


def expansion(x, base: int, limit: int = 10_000) -> Tuple[List[int], List[int], List[int]]:
    """Точное деление столбиком: (цифры целой части от старшей, предпериод, период).
    Период найден по повторению остатка — независимый путь к period()."""
    x = Q(x)
    base = _check_base(base)
    if x < 0:
        raise ValueError("expansion — для x ≥ 0; знак несёт запись")
    int_part, rem = divmod(x.n, x.d)
    int_digits = list(reversed(digits_of(int_part, base)))
    seen: Dict[int, int] = {}
    frac: List[int] = []
    while rem and rem not in seen:
        if len(frac) >= limit:
            raise Unresolved("expansion", f"{x.n}/{x.d}", str(base),
                             hint=f"{limit} цифр без повторения остатка — предел прогона, не факт о числе")
        seen[rem] = len(frac)
        rem *= base
        d, rem = divmod(rem, x.d)
        frac.append(d)
    if not rem:
        return int_digits, frac, []
    start = seen[rem]
    return int_digits, frac[:start], frac[start:]


def representations(n: int, base: int, width: int = 1) -> int:
    """b^w-ичные разбиения: число ненормированных записей целого n ≥ 0 в сетке (b, w)
    (мультимножества степеней R = b^w с суммой n): p(n) = p(n−1) + [R | n]·p(n/R)."""
    n = _check_int(n, "n")
    R = radix_of(base, width)
    if n < 0:
        raise ValueError("n ≥ 0")
    p = [1] * (n + 1)
    for k in range(1, n + 1):
        p[k] = p[k - 1] + (p[k // R] if k % R == 0 else 0)
    return p[n]


# ---------------------------------------------------------------------------
# Запись: конечное множество колёс с попарно разными e; ширина — параметр сетки
# ---------------------------------------------------------------------------

def normalize(counts: Dict[int, int], base: int, width=1) -> Tuple[Dict[int, int], int]:
    """Перенос как переписывание в сетке (b, w), R = b^w: f(e) ≥ R → f(e) − R, f(e+1) + 1.
    Без ширины (None) переноса нет: кратности просто складываются — свойство границы
    ширины, не числа. Возвращает (нормальная форма, число переносов); значение хранится."""
    R = radix_of(base, width)
    f = {e: c for e, c in counts.items() if c}
    for e, c in f.items():
        if c < 0:
            raise ValueError("кратности ≥ 0: вычитание — через значение (заём), не через ℕ^(ℤ)")
    if R is None:
        return dict(f), 0
    carries = 0
    while any(c >= R for c in f.values()):
        for e in sorted(f):
            c = f.get(e, 0)
            if c >= R:
                q, r = divmod(c, R)
                f[e] = r
                f[e + 1] = f.get(e + 1, 0) + q
                carries += q
    return {e: c for e, c in f.items() if c}, carries


class Record:
    """Запись числа колёсами: знак и колёса с попарно разными e. Ширина колеса —
    параметр сетки: width = w → |m| < b^w и позиции e кратны w (нормальная форма
    сетки (b, w)); width = None — колёса любой величины (Ф-число безразлично).
    Пустая запись — ноль как отсутствие колеса (класс 0)."""

    __slots__ = ("base", "sign", "width", "_w")

    def __init__(self, base: int, wheels: Iterable[Wheel] = (), sign: int = 1, width=None):
        base = _check_base(base)
        R = radix_of(base, width)
        if sign not in (1, -1):
            raise ValueError("знак записи — ±1")
        w: Dict[int, int] = {}
        for wh in wheels:
            if not isinstance(wh, Wheel):
                raise TypeError("запись собирается из Wheel")
            if wh.base != base:
                raise ValueError(f"колесо базы {wh.base} в записи базы {base}: смена сетки — rebase")
            if wh.e in w:
                raise ValueError(f"два колеса на позиции {wh.e}: не нормальная форма — Record.from_counts")
            if R is not None:
                if not 1 <= wh.m < R:
                    raise ValueError(f"мантисса {wh.m} вне 1…{R - 1} сетки ширины {width}: "
                                     "перенос или знак — не нормальная форма этой сетки")
                if wh.e % width:
                    raise ValueError(f"позиция {wh.e} не кратна ширине {width}")
            w[wh.e] = wh.m
        self.base, self.width = base, width
        self.sign = 1 if not w else sign
        self._w = w

    # --- построение ------------------------------------------------------------

    @classmethod
    def from_counts(cls, base: int, counts: Dict[int, int], sign: int = 1, width=1) -> "Record":
        f, _ = normalize(counts, base, width)
        return cls(base, [Wheel(base, c, e) for e, c in f.items()], sign, width)

    @classmethod
    def of(cls, x, base: int, width=None) -> "Record":
        """Нормальная форма x ∈ ℤ[1/b] в сетке (b, w); без ширины — одно колесо;
        иначе (x ∉ ℤ[1/b]) — отказ по имени с предпериодом и периодом."""
        x = Q(x)
        R = radix_of(base, width)
        if x == 0:
            return cls(base, (), 1, width)
        qb, rest, ell = split_denominator(x.d, base)
        if rest != 1:
            raise Unresolved("запись колёсами", f"{x.n}/{x.d}", f"база {base}",
                             hint=f"W = ω: знаменатель несёт простое вне базы (остаток {rest}); "
                                  f"запись периодична — предпериод {ell}, период {mult_order(base, rest)}")
        sign = 1 if x > 0 else -1
        if R is None:
            n = abs(x.n) * base ** ell // x.d
            return cls(base, [Wheel(base, n, -ell)], sign, None)
        shift = (-ell) % width
        n = abs(x.n) * base ** (ell + shift) // x.d
        wheels = [Wheel(base, d, (i - (ell + shift) // width) * width)
                  for i, d in enumerate(digits_of(n, R)) if d]
        return cls(base, wheels, sign, width)

    @classmethod
    def parse(cls, text: str, base: int = 10, sep: str = ",", width=1) -> "Record":
        """Позиционная запись с разделителем (',' или '.'; явные нули — традиция читаемости);
        по умолчанию — сетка ширины 1 (цифры), как в тексте; width=None — одно колесо."""
        base = _check_base(base)
        s = text.strip().replace(" ", "").replace("_", "").lower()
        sign = 1
        if s[:1] in "+-−":
            sign = -1 if s[0] in "-−" else 1
            s = s[1:]
        s = s.replace(".", sep).replace(",", sep)
        if s.count(sep) > 1:
            raise ValueError("больше одного разделителя")
        head, _, tail = s.partition(sep)
        if not head and not tail:
            raise ValueError("пустая запись")
        total = Q(0)
        for i, ch in enumerate(reversed(head)):
            d = _DIGITS.index(ch)
            if d >= base:
                raise ValueError(f"цифра {ch!r} вне базы {base}")
            total = total + Q(d * base ** i)
        for k, ch in enumerate(tail, start=1):
            d = _DIGITS.index(ch)
            if d >= base:
                raise ValueError(f"цифра {ch!r} вне базы {base}")
            total = total + Q(d, base ** k)
        return cls.of(total if sign > 0 else -total, base, width)

    @classmethod
    def from_tags(cls, base: int, tags, sign: int = 1, width=None) -> "Record":
        return cls(base, [Wheel.from_tag(base, t) for t in tags], sign, width)

    # --- чтение ----------------------------------------------------------------

    @property
    def wheels(self) -> Tuple[Wheel, ...]:
        return tuple(Wheel(self.base, self._w[e], e) for e in sorted(self._w, reverse=True))

    @property
    def positions(self) -> Tuple[int, ...]:
        return tuple(sorted(self._w, reverse=True))

    @property
    def rank(self) -> int:
        """W — число колёс."""
        return len(self._w)

    @property
    def span(self) -> int:
        """Разрядов основания b между старшим разрядом старшего колеса и позицией младшего; 0 у пустой."""
        if not self._w:
            return 0
        hi = max(e + len(digits_of(m, self.base)) - 1 for e, m in self._w.items())
        return hi - min(self._w) + 1

    @property
    def top_index(self) -> int:
        """Наибольший номер колеса K (для цены индекса); 0 у пустой записи."""
        return 0 if not self._w else max(w.index for w in self.wheels)

    @property
    def value(self) -> Q:
        total = Q(0)
        for w in self.wheels:
            total = total + w.value
        return total if self.sign > 0 else -total

    def tags(self) -> List[tuple]:
        return [w.tag() for w in self.wheels]

    def counts(self) -> Dict[int, int]:
        return dict(self._w)

    def is_zero(self) -> bool:
        return not self._w

    def container(self) -> PhiNumberV6:
        """Сумма колёс как один контейнер класса 1 (база 2 неявная); пустая запись — класс 0."""
        if not self._w:
            return PhiNumberV6(LockedInt(0), LockedInt(1), LockedInt(0), LockedInt(1), NumberClass.ZERO)
        v = self.value
        return PhiNumberV6.from_ratio(v.n, v.d)

    def containers(self) -> List[PhiNumberV6]:
        return [w.container() for w in self.wheels]

    def text(self, sep: str = ",") -> str:
        """Позиционная запись с явными нулями (традиция читаемости) — через значение."""
        if not self._w:
            return "0"
        v = abs(self.value)
        _, rest, ell = split_denominator(v.d, self.base)
        n = v.n * self.base ** ell // v.d
        ds = digits_of(n, self.base)
        while len(ds) <= ell:
            ds.append(0)
        ints = "".join(_DIGITS[d] for d in reversed(ds[ell:]))
        out = ("-" if self.sign < 0 else "") + ints
        if ell:
            out += sep + "".join(_DIGITS[d] for d in reversed(ds[:ell]))
        return out

    # --- алгебра: одометр — точное значение, затем нормальная форма той же сетки ----

    def _same_grid(self, other: "Record") -> None:
        if not isinstance(other, Record):
            raise TypeError("операнд — Record")
        if other.base != self.base or other.width != self.width:
            raise ValueError("разные сетки: сначала regrid")

    def __add__(self, other: "Record") -> "Record":
        self._same_grid(other)
        return Record.of(self.value + other.value, self.base, self.width)

    def __sub__(self, other: "Record") -> "Record":
        self._same_grid(other)
        return Record.of(self.value - other.value, self.base, self.width)

    def __mul__(self, other: "Record") -> "Record":
        self._same_grid(other)
        return Record.of(self.value * other.value, self.base, self.width)

    def __truediv__(self, other: "Record") -> "Record":
        self._same_grid(other)
        if other.is_zero():
            raise ZeroDivisionError("деление на отсутствие колеса")
        return Record.of(self.value / other.value, self.base, self.width)

    def __neg__(self) -> "Record":
        return Record(self.base, self.wheels, -self.sign, self.width)

    def union(self, other: "Record") -> Dict[int, int]:
        """Ненормированная запись суммы: мультимножество колёс (до переноса); знаки равные."""
        self._same_grid(other)
        if self.sign != other.sign and not (self.is_zero() or other.is_zero()):
            raise ValueError("объединение — для равных знаков; разность — через заём (значение)")
        out = dict(self._w)
        for e, c in other._w.items():
            out[e] = out.get(e, 0) + c
        return out

    def scaled(self, k: int) -> "Record":
        """Сдвиг по оси колёс: значение·b^k — те же колёса, W не меняется (k кратно ширине)."""
        k = _check_int(k, "k")
        if self.width is not None and k % self.width:
            raise ValueError(f"сдвиг {k} не кратен ширине {self.width}")
        return Record(self.base, [Wheel(self.base, w.m, w.e + k) for w in self.wheels], self.sign, self.width)

    def regrid(self, width) -> "Record":
        """Та же величина в сетке иной ширины (безразличие Ф-числа к размеру колеса).

        Одно понятие с `gridcore.Notch.regrid` и `affine.Affine.regrid` —
        смена сетки при той же величине. Различие — в исходе, когда величина
        в сетку не входит: здесь исхода два — точно либо отказ по имени
        (`Record.of`: `Unresolved` с предпериодом и периодом), округления нет
        никогда, потому что запись держит ОДНО число; засечка `Notch` на
        грубой сетке округляется НАРУЖУ (отрезок шире), форма `Affine` режет
        коэффициенты к нулю и отдаёт недостачу новому символу. Сетка здесь —
        ширина колеса w при том же основании; у соседей — шаг 2^−g при
        основании 2 (сдвиг по оси колёс).
        """
        return Record.of(self.value, self.base, width)

    def rebase(self, base: int, width=None) -> "Record":
        """Смена основания: та же величина, другая колёсная судьба (W — инвариант пары).

        Основание здесь ПОЗИЦИОННОЕ: b ∈ ℕ, b ≥ 2, колёса — m·b^e; смена
        основания меняет ранг (W₁₀(1/3) = ω, W₃(1/3) = 1) и может отказать по
        имени. Одноимённая дверь `zech.Zech.rebase` / `zech.ZechExponent.rebase`
        — то же понятие (та же величина при другом основании), но основание
        ТРАНСЦЕНДЕНТНОЕ: α = β^s, показатель масштабируется на s, колёс и
        отказов нет — там меняется база логарифма сложения, здесь — сетка
        цифр. Общее у обеих — величина; различное — что именно есть «основание».
        Дверь `Hereditary.rebase` — та же смена основания, но над деревом
        показателей: там сохраняется ФОРМА, а величина пересчитывается
        (шаг Гудстейна).
        """
        return Record.of(self.value, base, width)

    # --- цена ----------------------------------------------------------------------

    def cost(self, index_bits: Optional[int] = None, mantissa_bits: Optional[int] = None) -> Dict[str, int]:
        """Биты: колёса (мантиссы + индекс + зона/знак + счёт), позиционная пара, позиционные цифры.
        Поле мантиссы: в сетке ширины w — фиксированное bits(b^w − 1) на колесо (поле регистра);
        без ширины — bits(|m|) каждого колеса; mantissa_bits задаёт поле явно."""
        W, K, b = self.rank, self.top_index, self.base
        idx = bits_for(K) if index_bits is None else index_bits
        if mantissa_bits is not None:
            mant = W * mantissa_bits
        elif self.width is not None:
            mant = W * bits_for(b ** self.width - 1)
        else:
            mant = sum(bits_for(m) for m in self._w.values())
        wheels = mant + W * (idx + 1) + bits_for(W)
        v = self.value
        pair_num, pair_den = bits_for(abs(v.n)), bits_for(v.d)
        digits = self.span * bits_for(b - 1)
        return {"wheels": wheels, "mantissa_bits": mant, "index_bits": idx, "digit_bits": bits_for(b - 1),
                "count_bits": bits_for(W), "pair_num": pair_num, "pair_den": pair_den,
                "pair": pair_num + pair_den, "digits": digits, "span": self.span,
                "empty_positions": self.span - sum(len(digits_of(m, b)) for m in self._w.values())}

    # --- служебное ---------------------------------------------------------------------

    def __eq__(self, other) -> bool:
        return (isinstance(other, Record) and self.base == other.base and self.width == other.width
                and self.sign == other.sign and self._w == other._w)

    def __hash__(self) -> int:
        return hash((self.base, self.width, self.sign, tuple(sorted(self._w.items()))))

    def __repr__(self) -> str:
        if not self._w:
            return "∅ (класс 0: колеса нет)"
        body = " + ".join(repr(w) for w in self.wheels)
        return ("−(" + body + ")") if self.sign < 0 else body


# ---------------------------------------------------------------------------
# Оптимальная нарезка: размер колеса свободен — колёса выбирает цена, не догма
# ---------------------------------------------------------------------------

def span_overhead(x, base: int) -> int:
    """Биты индекса и знака/зоны на колесо для x: bits(K) + 1, K — наибольший номер разряда
    x относительно разделителя; одинаков для всех записей одного x (единая мера цены)."""
    x = Q(x)
    base = _check_base(base)
    if x == 0:
        return 1
    _, rest, ell = split_denominator(x.d, base)
    if rest != 1:
        raise Unresolved("цена колёс", f"{x.n}/{x.d}", f"база {base}", hint="W = ω")
    n = abs(x.n) * base ** ell // x.d
    ds = digits_of(n, base)
    K = max(abs(i - ell) + (1 if i - ell >= 0 else 0) for i, d in enumerate(ds) if d)
    return bits_for(K) + 1


def uniform_cost(rec: Record, overhead: Optional[int] = None) -> int:
    """Единая мера цены записи: Σ bits(|m_i|) + W·overhead + bits(W); overhead по умолчанию —
    span_overhead значения (один для всех записей того же x — сравнимы сетки и нарезки)."""
    ov = span_overhead(rec.value, rec.base) if overhead is None else _check_int(overhead, "overhead")
    return sum(bits_for(m) for m in rec._w.values()) + rec.rank * ov + bits_for(rec.rank)


def optimal(x, base: int, overhead: Optional[int] = None) -> Tuple[Record, int]:
    """Запись x ∈ ℤ[1/b] колёсами любой ширины с наименьшей ценой
    Σ (bits(|m_i|) + overhead) + bits(W): колесо покрывает отрезок разрядов
    основания b, внутренние нули стоят как разряды; overhead — биты индекса и
    знака/зоны на колесо (по умолчанию bits(K) + 1 при K — наибольшем номере
    разряда). Динамика по ненулевым разрядам, точно; возвращает (запись, цена)."""
    x = Q(x)
    base = _check_base(base)
    if x == 0:
        return Record(base), 0
    _, rest, ell = split_denominator(x.d, base)
    if rest != 1:
        raise Unresolved("нарезка колёсами", f"{x.n}/{x.d}", f"база {base}",
                         hint="W = ω: x ∉ ℤ[1/b] — колёс конечного числа нет")
    n = abs(x.n) * base ** ell // x.d
    ds = digits_of(n, base)                       # ds[i] — разряд позиции i − ell
    nz = [i for i, d in enumerate(ds) if d]       # позиции ненулевых, по возрастанию
    ov = span_overhead(x, base) if overhead is None else _check_int(overhead, "overhead")

    def mant(lo: int, hi: int) -> int:            # мантисса колеса, покрывающего ds[lo..hi]
        v = 0
        for i in range(hi, lo - 1, -1):
            v = v * base + ds[i]
        return v

    N = len(nz)
    best = [0] + [None] * N                       # best[j] — цена первых j ненулевых разрядов
    cut = [0] * (N + 1)
    for j in range(1, N + 1):
        for i in range(j):                        # колесо покрывает nz[i] … nz[j−1]
            c = best[i] + bits_for(mant(nz[i], nz[j - 1])) + ov
            if best[j] is None or c < best[j]:
                best[j], cut[j] = c, i
    wheels: List[Wheel] = []
    j = N
    while j > 0:
        i = cut[j]
        wheels.append(Wheel(base, mant(nz[i], nz[j - 1]), nz[i] - ell))
        j = i
    rec = Record(base, wheels, 1 if x > 0 else -1, None)
    return rec, best[N] + bits_for(rec.rank)


# ---------------------------------------------------------------------------
# Байт колеса по правилам — одна из сеток: ширина 1, позиции e ∈ [−K, K] × цифры 1…b−1 ≤ 256
# ---------------------------------------------------------------------------

class WheelByte:
    """Один байт на колесо, упаковка по правилам (не по биту на значение):
    код = (e + K)·(b − 1) + (d − 1); K = ⌊(⌊256/(b−1)⌋ − 1)/2⌋; вне [−K, K] — отказ.
    Это сетка ширины 1 (цифры 1…b−1) — выбор представления, не свойство Ф-числа."""

    __slots__ = ("base", "reach", "codes")

    def __init__(self, base: int):
        base = _check_base(base)
        positions = 256 // (base - 1)
        if positions < 1:
            raise ValueError(f"база {base}: цифр больше 256 — байта нет")
        K = (positions - 1) // 2
        self.base, self.reach, self.codes = base, K, (2 * K + 1) * (base - 1)

    def pack(self, wheel: Wheel) -> int:
        """Колесо сетки ширины 1 → один байт по правилам (позиция и цифра).

        Одноимённые двери соседей — НЕ то же: `ccore.pack` переводит партию
        пар в структуры C-ядра (перевод значения через границу библиотеки),
        `tag.pack` кладёт в байт адрес таблицы состояний (индекс состояния и
        код α, два поля фиксированной ширины). Здесь байт несёт колесо:
        код = (e + K)·(b − 1) + (d − 1), вне досягаемости — отказ по имени.
        """
        if wheel.base != self.base:
            raise ValueError("колесо другой базы")
        if not 1 <= wheel.m < self.base:
            raise Unresolved("байт колеса", repr(wheel), f"база {self.base}",
                             hint=f"мантисса {wheel.m} шире одного разряда — байт кодирует сетку ширины 1")
        if abs(wheel.e) > self.reach:
            raise Unresolved("байт колеса", repr(wheel), f"база {self.base}",
                             hint=f"позиция {wheel.e} вне досягаемости ±{self.reach}")
        return (wheel.e + self.reach) * (self.base - 1) + (wheel.m - 1)

    def unpack(self, byte: int) -> Wheel:
        """Байт → колесо сетки ширины 1; запасной код — отказ по имени.

        Не `ccore.unpack` (структуры C-ядра → значения эталона) и не
        `tag.unpack` (байт → состояние таблицы и код α): те двери не знают
        колеса, эта — только его.
        """
        byte = _check_int(byte, "byte")
        if not 0 <= byte < 256:
            raise ValueError("байт — 0…255")
        if byte >= self.codes:
            raise Unresolved("байт колеса", str(byte), f"база {self.base}",
                             hint=f"запасной код (кодов по правилам {self.codes})")
        e, d = divmod(byte, self.base - 1)
        return Wheel(self.base, d + 1, e - self.reach)

    def pack_record(self, rec: Record) -> bytes:
        if rec.base != self.base:
            raise ValueError("запись другой базы")
        return bytes(self.pack(w) for w in rec.wheels)

    def unpack_record(self, data: bytes, sign: int = 1) -> Record:
        return Record(self.base, [self.unpack(b) for b in data], sign, 1)


# ── Наследственная форма: рекурсивное самоприменение колеса ──────────────────
#
# Слово автора 24.09 `[в корпусе]` (коррекция d1eba5c2eb): «десятичное
# представление числа есть сокращённая форма записи сумм колёс 1-9. какие из
# этого выводы сделаны? … рекурсивное самоприменение ОБЯЗАТЕЛЬНО».
#
# Вывод, который отсюда следует буквально: у колеса m·b^e показатель e —
# тоже число, и его запись — тоже сумма колёс; так до тех пор, пока показатели
# не станут меньше b (одно колесо на основной сетке) либо нулём (отсутствие
# колеса). Это не новая конструкция дома, а классический предмет
# `[внешний источник]`: наследственное представление по основанию
# (hereditary base-b representation) — R. L. Goodstein, «On the restricted
# ordinal theorem», The Journal of Symbolic Logic, vol. 9, no. 2, 1944,
# pp. 33–41, https://doi.org/10.2307/2268019; независимость теоремы Гудстейна
# от арифметики Пеано — L. Kirby, J. Paris, «Accessible independence results
# for Peano arithmetic», Bulletin of the London Mathematical Society, vol. 14,
# no. 4, 1982, pp. 285–293, https://doi.org/10.1112/blms/14.4.285.
# Проверяемый след: последовательность Гудстейна от 4 — OEIS A056193
# (https://oeis.org/A056193): 4, 26, 41, 60, 83, 109, 139, 173, 211, 253, 299,
# 348, 401, 458, 519, 584, 653, 726, 803, 884, 969, 1058, 1151, 1222, 1295,
# 1370, 1447, 1526, 1607, 1690 — приёмка `proto/test_phicore_hereditary.py`
# считает их домом (Python и ядро C) и сверяет дословно.
#
# Самоприменение: `Hereditary` есть `Record.of`, применённый к собственному
# выходу — к показателям своих же колёс; неподвижные точки рекурсии — числа
# меньше основания (одно колесо, e = 0). Смена основания дерева (Гудстейн:
# «bump the base») — то же понятие, что `Record.rebase` у записи и
# `zech.ZechExponent.rebase` у показателя над трансцендентной базой: число
# остаётся, колёсная судьба меняется; здесь смена основания ЧИТАЕТ дерево в
# новой базе, потому значение меняется — это и есть шаг Гудстейна.

class Hereditary:
    """Наследственная (по Гудстейну) форма целого n ≥ 0 в базе b.

    Колёса `Record.of(n, b, 1)` (цифры 1…b−1 с показателями); показатель
    каждого колеса e ≥ 1 — сам `Hereditary(e, b)`; e = 0 — отсутствие
    (лист). Значение — Σ m·b^{value(e)}. `value(other)` читает то же дерево
    в другой базе — шаг Гудстейна без вычитания единицы.

    Инварианты `[вычислено здесь]`: value() == n; depth(n) ≤ ⌈log_b n⌉ и
    убывает как итерированный логарифм; count() ≥ W(n; b, 1) (колёса всех
    уровней не меньше колёс верхнего); Hereditary(n, b) при n < b — одно
    колесо без потомков (неподвижная точка).
    """

    __slots__ = ("base", "wheels")

    def __init__(self, n, base: int):
        base = _check_base(base)
        n = _check_int(n, "n")
        if n < 0:
            raise ValueError("наследственная форма — для n ≥ 0; знак несёт запись")
        self.base = base
        rec = Record.of(n, base, 1)
        ws = []
        for w in rec.wheels:
            m, _, e, _ = w.fields()
            ws.append((m, None if e == 0 else Hereditary(e, base)))
        self.wheels = tuple(ws)

    def value(self, base=None) -> int:
        """Значение дерева, прочитанного в базе `base` (по умолчанию — своей)."""
        b = self.base if base is None else _check_base(base)
        total = 0
        for m, sub in self.wheels:
            total += m * b ** (0 if sub is None else sub.value(b))
        return total

    def rebase(self, base: int) -> "Hereditary":
        """То же дерево, прочитанное в базе `base`, снова как форма в этой базе.

        Одно понятие с `Record.rebase` (запись), `zech.Zech.rebase` и
        `ZechExponent.rebase` (показатель над трансцендентной базой): смена
        основания. Различие —
        там величина сохраняется, здесь сохраняется ФОРМА, а величина
        пересчитывается (шаг Гудстейна)."""
        return Hereditary(self.value(base), base)

    def depth(self) -> int:
        if not self.wheels:
            return 0
        return 1 + max((0 if sub is None else sub.depth()) for _, sub in self.wheels)

    def count(self) -> int:
        """Колёс на всех уровнях."""
        return sum(1 + (0 if sub is None else sub.count()) for _, sub in self.wheels)

    def tree(self):
        """Вложенные кортежи (m, поддерево | None) — от младшего колеса."""
        return tuple((m, None if sub is None else sub.tree()) for m, sub in self.wheels)

    def containers(self) -> list:
        """Верхний уровень как Ф-контейнеры (m, 1, e, 1) над α = b; e — число."""
        return [Wheel(self.base, m, 0 if sub is None else sub.value()).container()
                for m, sub in self.wheels]

    def __eq__(self, other) -> bool:
        return isinstance(other, Hereditary) and self.base == other.base and self.tree() == other.tree()

    def __hash__(self) -> int:
        return hash((self.base, self.tree()))

    def __repr__(self) -> str:
        def show(t, b):
            parts = []
            for m, sub in t:
                if sub is None:
                    parts.append(f"{m}")
                else:
                    parts.append(f"{m}·{b}^({show(sub, b)})")
            return " + ".join(parts) if parts else "0"
        return f"Hereditary(base={self.base}: {show(self.tree(), self.base)})"


def goodstein(n: int, steps: int) -> List[int]:
    """Последовательность Гудстейна: G₀ = n; G_{k+1} = (наследственная форма G_k
    в базе k+2, прочитанная в базе k+3) − 1; останов на нуле. Возвращает не
    более `steps + 1` членов. `[внешний источник]` Goodstein 1944; OEIS A056193."""
    n = _check_int(n, "n")
    steps = _check_int(steps, "steps")
    if n < 0 or steps < 0:
        raise ValueError("n ≥ 0, steps ≥ 0")
    out = [n]
    base = 2
    for _ in range(steps):
        if out[-1] == 0:
            break
        out.append(Hereditary(out[-1], base).value(base + 1) - 1)
        base += 1
    return out
