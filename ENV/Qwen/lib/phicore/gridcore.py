"""gridcore — новый прототип целочисленного ядра: ПОЛНОЕ ДИСКРЕТНОЕ ПРЕДСТАВЛЕНИЕ.

Слово автора 03.09: «создавай НОВЫЙ прототип intcore и ставь в него всё,
что необходимо проекту; осмысленные эксперименты и проверки доказуемости
гипотез — рабочий инструмент»; «мы дискретные: засечка π на линейке
конечной длины имеет 2 значения в беззнаковой и 3 в знаковой линейке»;
«многие гипотезы становятся доказуемыми или тривиальными в ПОЛНОМ
дискретном представлении».

Устройство. Линейка конечной длины с шагом 2^−g. Величина на ней — ЗАСЕЧКА
`Notch`: два целых узла lo ≤ hi (беззнаковая линейка — два значения),
между которыми величина; на знаковой линейке — три значения: знак,
|lo|, |hi| (`signed()`); ноль — засечка, накрывающая узел 0. Вся
арифметика — на целых: сложение и вычитание узлов, произведение и
частное с округлением НАРУЖУ к сетке, корень целым счётом; ни одного
float, ни одной пары в горячем цикле. Носитель любой линейки входит
через свой сертификат сетки (`Notch.of`): ручка уточнения крутится,
пока засечка не станет не шире 2^−g, либо объявляется достигнутая ширина.

Что даёт дискретное представление, чего не даёт гладкое: ПОРЯДОК и ЗНАК
двух разных величин решаются за конечное число шагов (две разные
величины разводит любая достаточно тонкая сетка); равенство на сетке g
— совпадение засечек (определение, а не гипотеза); ноль на сетке —
засечка через узел 0. Что оно не даёт: равенство в пределе g → ∞ засечка
не доказывает (иначе — ложный класс, §316.4); это остаётся форме
(ratfunc) там, где форма есть, и не утверждается там, где её нет.

`Mixed` — СТРУКТУРНАЯ ЗАПИСЬ значения, у которого линейки-носителя нет
(e + π, π + √2, K·E, j₀,₁ + j₀,₂): дерево действий над носителями.
Запись точна как есть (ничего не округлено), её значение — предел
засечек по сетке (предел как способ получения — класс 19), засечка на
любой сетке — по явному запросу. Равенство записей — только структурное
или на названной сетке; «=» в пределе — не утверждается. Так снимается
отказ «разные базы» (слово автора 03.09: «отказ НЕ принят») без
костыля: ни Шануэль, ни общее основание не требуются — их требовала
гладкая постановка, а не дискретная.

Авторство (признаётся безоговорочно, слово автора 03.09: «Кесарю —
кесарево, слесарю — слесарево»): интервальная арифметика — Moore 1966;
стандарт интервалов — IEEE 1788-2015; точный аккумулятор скалярного
произведения — Kulisch; логарифмическая система счисления — Swartzlander
и Alexopoulos 1975, для нейросетей — Miyashita, Lee, Murmann 2016;
логарифм сложения — Leonelli 1803, Gauss 1812, Zech; целочисленный вывод с
двоичными весами — Jacob и др. 2018; unum/posit как иная дискретная
линейка — Gustafson 2017. Здесь своё — только состав: засечка как два
узла, три значения на знаковой линейке, структурная запись с засечкой
по запросу и место всего этого в онтологии классов.
"""
from __future__ import annotations

from .enclosure import Interval, int_root_floor
from .errors import Unresolved
from .rational import Q

_NUM = (int, Q)
_KNOBS = (30, 60, 120, 240, 480, 960)      # ручка уточнения сертификата носителя


def _floor_div(a: int, b: int) -> int:
    return a // b


def _ceil_div(a: int, b: int) -> int:
    return -((-a) // b)


class Notch:
    """Засечка на сетке 2^−g: целые узлы lo ≤ hi, величина в [lo; hi]·2^−g.

    ЗАКОН КОЛЕСА (слово автора 09.09: число есть сумма колёс, 0 — отсутствие
    колеса; `phicore/wheels.py`). Узел засечки `units·2^−g` есть контейнер
    (units, 1, −g, 1) над α = 2 — одно колесо базы 2 без ширины, а в сетке
    (2, 1) — сумма колёс по единичным битам `units`; сдвиг сетки g — сдвиг по
    оси колёс, не смена величины. Сдвиги и деления на 2^g в этом классе —
    БЫСТРАЯ ТРОПА к тому же числу (двери горячие: `regrid` зовётся из каждой
    пары действий через `_align`), а не отдельная истина: приёмка
    `proto/test_phicore_gridcore.py` сверяет узлы с `wheels.Record.of(Q(units,
    2**g), 2)` перебором и с подлогом, мост наружу — `wheel_records`.
    """

    __slots__ = ("g", "lo", "hi")

    def __init__(self, g: int, lo: int, hi: int):
        if not isinstance(g, int) or isinstance(g, bool) or g < 0:
            raise TypeError("сетка g — целое ≥ 0")
        if not (isinstance(lo, int) and isinstance(hi, int)) \
                or isinstance(lo, bool) or isinstance(hi, bool):
            raise TypeError("узлы засечки — целые")
        if lo > hi:
            raise ValueError(f"пустая засечка: {lo} > {hi}")
        object.__setattr__(self, "g", g)
        object.__setattr__(self, "lo", lo)
        object.__setattr__(self, "hi", hi)

    def __setattr__(self, *_):
        raise AttributeError("засечка неизменяема")

    # Удаление поля закрыто ТОЙ ЖЕ формой, что присваивание: один объект
    # функции на две двери, а не вторая копия строки отказа (§14 п. 11).
    # Копия расходится — разнобой §14 пп. 9, 10 был именно в ВИДЕ отказа;
    # здесь `del n.lo` физически не может ответить иначе, чем `n.lo = 5`.
    # Цена замерена и равна нулю — числа в докстринге неизменяемости `Affine`.
    __delattr__ = __setattr__

    # --- построение ------------------------------------------------------
    @classmethod
    def from_q(cls, q, g: int) -> "Notch":
        """Пара на сетке: узел точно, если делится; иначе два соседних."""
        q = Q(q)
        scale = 1 << g
        lo = _floor_div(q.numerator * scale, q.denominator)
        hi = _ceil_div(q.numerator * scale, q.denominator)
        return cls(g, lo, hi)

    @classmethod
    def from_interval(cls, iv: Interval, g: int) -> "Notch":
        """Отрезок ⊂ ℚ — засечкой: концы округляются наружу."""
        scale = 1 << g
        lo = _floor_div(iv.lo.numerator * scale, iv.lo.denominator)
        hi = _ceil_div(iv.hi.numerator * scale, iv.hi.denominator)
        return cls(g, lo, hi)

    @classmethod
    def of(cls, x, g: int) -> "Notch":
        """Любой носитель — засечкой на сетке g: пара точно; носитель с
        сертификатом — ручка уточнения до ширины ≤ 2^−g (иначе засечка
        шире, и это видно по width_units).

        Почему float и bool идут в дверь `from_q`, а не в ветку сертификата
        (правка §14 п. 10). До правки они падали в конец и получали Ф-отказ
        `Unresolved` «у носителя нет сертификата сетки», тогда как
        `Notch.from_q(1.5, 8)` на тот же вход отвечал TypeError «float
        запрещён построением; явная дверь — Q.from_float». Один вход давал
        два разных ВИДА отказа, а ветвление у вызывающего идёт по типу.
        Диагноз `Unresolved` при этом был неточен: `Unresolved` называет
        ВЕЛИЧИНУ и путь её уточнения, а у float пути уточнения нет — на него
        стоит запрет ТИПА, канон корпуса (`phicore/rational.py`: «float/bool
        на входе — TypeError; единственная дверь из мира IEEE — явный
        `Q.from_float`»). Поэтому здесь не заведён свой отказ с теми же
        словами (две копии канона разошлись бы снова), а вход просто пущен в
        ту же дверь, что и у соседа: `Q` отказывает сам, и обе двери модуля
        отвечают дословно одним сообщением. `bool` — подкласс `int`, так что
        снятая оговорка `not isinstance(x, bool)` отправляет его туда же, к
        отказу «bool не число». `Unresolved` остаётся за теми, у кого
        сертификата действительно нет (`'π'`, носитель без `enclosure`).
        """
        if isinstance(x, Notch):
            return x.regrid(g)
        if isinstance(x, _NUM) or isinstance(x, float):
            return cls.from_q(x, g)
        if isinstance(x, Interval):
            return cls.from_interval(x, g)
        enc = getattr(x, "enclosure", None)
        if not callable(enc):
            raise Unresolved("засечка", repr(x),
                             hint="у носителя нет сертификата сетки — на "
                                  "дискретную линейку не входит")
        best = None
        for k in _KNOBS:
            iv = enc(k)
            n = cls.from_interval(iv, g)
            if best is None or n.width_units < best.width_units:
                best = n
            if n.width_units <= 1:
                break
        return best

    # --- устройство ------------------------------------------------------
    @property
    def width_units(self) -> int:
        """Ширина в узлах: 0 — величина на узле, 1 — между соседями."""
        return self.hi - self.lo

    @property
    def sign(self) -> int:
        """+1 / −1, если засечка отделена от нуля; 0 — накрывает узел 0
        (ноль на этой сетке; не есть доказательство нуля величины)."""
        if self.lo > 0:
            return 1
        if self.hi < 0:
            return -1
        return 0

    def signed(self) -> tuple:
        """Три значения знаковой линейки: (знак, |lo|, |hi|) — при знаке
        0 модули берутся как есть (узел 0 внутри)."""
        s = self.sign
        if s > 0:
            return (1, self.lo, self.hi)
        if s < 0:
            return (-1, -self.hi, -self.lo)
        return (0, self.lo, self.hi)

    def regrid(self, g: int) -> "Notch":
        """Та же засечка на другой сетке: тоньше — точно (сдвиг), грубее —
        наружу.

        Одно понятие с `wheels.Record.regrid` (смена сетки при той же
        величине) и с `affine.Affine.regrid`; различие в том, что делается,
        когда величина в новую сетку не входит. Тоньше — у всех троих точно:
        узлы сдвигаются влево, и запись колёсами базы 2 у узла ТА ЖЕ
        (`Record.of(Q(lo, 2**g), 2, 1)` не зависит от g — безразличие
        Ф-числа к сетке). Грубее — здесь засечка округляется НАРУЖУ
        (`lo >> sh`, `⌈hi / 2^sh⌉`: включение соблюдено, величина ШИРЕ);
        `Record.regrid` наружу не округляет никогда — либо точен, либо
        отказ по имени (`Unresolved` с предпериодом и периодом), потому что
        запись держит ОДНО число, а засечка — отрезок; `Affine.regrid`
        режет коэффициенты к нулю и отдаёт отброшенное новому символу.
        Различие названо, а не сглажено: приёмка сверяет все три исхода.
        """
        if g == self.g:
            return self
        if g > self.g:
            sh = g - self.g
            return Notch(g, self.lo << sh, self.hi << sh)
        sh = self.g - g
        return Notch(g, self.lo >> sh, _ceil_div(self.hi, 1 << sh))

    def _canon_triple(self) -> tuple:
        """Тройка канона (g, lo, hi): та же засечка на самой ГРУБОЙ сетке, где
        оба узла ещё целые. Общие младшие нули узлов снимаются вместе с сеткой:
        ⟦16; 16⟧·2^−4 → ⟦1; 1⟧·2^−0. Узлы нулевые общих нулей не ограничивают,
        потому засечка узла 0 канонизуется до сетки 0 (⟦0; 0⟧·2^−4 → ⟦0; 0⟧·2^−0).
        """
        k = self.g
        for v in (self.lo, self.hi):
            if v:
                k = min(k, (v & -v).bit_length() - 1)
        if not k:
            return (self.g, self.lo, self.hi)
        return (self.g - k, self.lo >> k, self.hi >> k)

    def canon(self) -> "Notch":
        """Канон записи: равные засечки дают ОДНУ запись канона (§14 п. 5).

        Сдвиг вправо здесь точен, а не «наружу», как у `regrid` на грубую
        сетку: снимаются только общие нулевые младшие биты обоих узлов, то
        есть деление без остатка. Потому канон — та же величина, а не её
        расширение: `canon().regrid(g)` возвращает исходную засечку.
        """
        g, lo, hi = self._canon_triple()
        return self if g == self.g else Notch(g, lo, hi)

    def to_interval(self) -> Interval:
        scale = 1 << self.g
        return Interval(Q(self.lo, scale), Q(self.hi, scale))

    def wheel_records(self, width=1) -> tuple:
        """Узлы засечки как записи колёсами базы 2: (запись lo, запись hi).

        Мост к закону (`wheels.Record.of`), не тропа: узел `units·2^−g` лежит
        в ℤ[1/2] построением, потому запись конечна всегда и отказа здесь
        нет; `width=1` — сетка (2, 1), колесо на каждый единичный бит,
        `None` — одно колесо (units, e = −g после снятия общих нулей).
        Ввоз ленивый: `wheels` ввозит `container`, а засечка — ядро без него.
        """
        from .wheels import Record
        scale = 1 << self.g
        return (Record.of(Q(self.lo, scale), 2, width),
                Record.of(Q(self.hi, scale), 2, width))

    def self_check(self, deep: bool = False) -> list:
        """Самопроверка засечки независимой дорогой (ось «AIAΣA спрашивает AIAΣA»).

        `[вычислено здесь]` Устройство проверяет само себя не тем же счётом, каким
        строило: узлы сверяются умножением обратно (`lo·2^−g ≤ hi·2^−g`), ширина —
        разностью узлов против `width_units`, принадлежность — сравнением дробей,
        а не сравнением узлов. Возвращает список находок; пустой список — молчание.
        Очередь п. 18: расширение самоприменения с трёх модулей.
        """
        out = []
        lo, hi = Q(self.lo, 1 << self.g), Q(self.hi, 1 << self.g)
        if lo > hi:
            out.append(f"узлы переставлены: {self.lo} > {self.hi}")
        if self.hi - self.lo != self.width_units:
            out.append(f"ширина не сходится: узлов {self.hi - self.lo}, объявлено {self.width_units}")
        for q in (lo, hi, (lo + hi) / 2 if lo != hi else lo):
            if not self.contains(q):
                out.append(f"величина внутри засечки не признана своей: {q}")
        if deep:
            back = Notch.from_q(lo, self.g)
            if back.lo != self.lo:
                out.append(f"обратное построение не совпало: {back.lo} против {self.lo}")
            for outside in (lo - Q(1, 1 << self.g), hi + Q(1, 1 << self.g)):
                if self.contains(outside):
                    out.append(f"чужая величина признана своей: {outside}")
            twin = self.regrid(self.g + 3)          # та же величина, другая сетка
            if twin == self and hash(twin) != hash(self):
                out.append("хэш не согласован с равенством")
            back = self.canon().regrid(self.g)       # канон точен: назад узел в узел
            if (back.g, back.lo, back.hi) != (self.g, self.lo, self.hi):
                out.append(f"канон записи не та же величина: {self.canon().record()}")
        return out

    def contains(self, q) -> bool:
        q = Q(q)
        scale = 1 << self.g
        return (self.lo * q.denominator <= q.numerator * scale
                <= self.hi * q.denominator)

    def record(self) -> str:
        return f"⟦{self.lo}; {self.hi}⟧·2^−{self.g}"

    __repr__ = record

    def __eq__(self, other):
        """Совпадение засечек на одной сетке — равенство НА СЕТКЕ (определение)."""
        if not isinstance(other, Notch):
            return NotImplemented
        a, b = _align(self, other)
        return a.lo == b.lo and a.hi == b.hi

    def __hash__(self):
        """Контракт eq/hash: равные засечки — равный хэш (§14 п. 5).

        Считался по тройке (g, lo, hi) как она записана, тогда как равенство
        приводит засечки к общей сетке: `Notch(0,1,1) == Notch(4,16,16)`
        истинно, а хэши расходились, и засечка, положенная в словарь на одной
        сетке, по равной ей на другой НЕ находилась (замер 12.09:
        `{Notch(0,1,1): …}.get(Notch(4,16,16))` → `None`, `len({a, b})` → 2).
        Дыра той же природы, что п. 13 главы 06 у `ZechExponent`, и закрыта
        тем же приёмом корпуса (`zech.py:172`: «контракт eq/hash: E == q ⇒
        hash(E) == hash(q)») — КАНОНОМ записи, а не приведением пары к общей
        сетке: канон у каждой засечки свой и считается без второго операнда,
        поэтому одинаков у всех равных ей, сколько бы их ни было. У засечки,
        уже стоящей в каноне, хэш считается по той же тройке, что и до правки
        (`hash(Notch(8,85,86)) == hash(("Notch", 8, 85, 86))` — истинно).
        Замер после правки перебором 1071 засечки (g ≤ 6, узлы −8…8): равных
        пар 1446, хэш расходится у 0.
        """
        return hash(("Notch",) + self._canon_triple())

    # --- порядок: решается, когда засечки не пересекаются ---------------
    def __lt__(self, other):
        a, b = _align(self, _as_notch(other, self.g))
        return a.hi < b.lo

    def __gt__(self, other):
        a, b = _align(self, _as_notch(other, self.g))
        return a.lo > b.hi

    def separated(self, other) -> bool:
        """Засечки не пересекаются: величины различны — ДОКАЗАНО сеткой."""
        a, b = _align(self, _as_notch(other, self.g))
        return a.hi < b.lo or b.hi < a.lo

    def indistinguishable(self, other) -> bool:
        """Засечки пересекаются: на этой сетке величины не различимы."""
        return not self.separated(other)

    # --- арифметика на целых ---------------------------------------------
    def __add__(self, other):
        """Сумма засечек; ЦЕЛОЕ ЗДЕСЬ — ВЕЛИЧИНА, не узел.

        Оговорка 12.09 (§14 п. 8 главы 08: «целое значит разное в двух модулях
        одной главы, и ни один докстринг этого не оговаривает»). Замер, одна и
        та же единица у двух соседей:
          * засечка: `3 + Notch.from_q(1, 4)` → `⟦64; 64⟧·2^−4`, то есть 4 —
            целое вошло ВЕЛИЧИНОЙ через дверь `Notch.of` (`_as_notch`);
          * форма: `A + 1` при `A = Affine.from_notch(⟦256; 512⟧·2^−8)` →
            `(385 + 128·ε)·2^−8` — целое вошло УЗЛОМ сетки, то есть 2^−8
            (`affine.Affine.__add__`, там та же оговорка с другой стороны).
        Разница намеренная: засечка несёт ВЕЛИЧИНУ на линейке, и всякий вход
        проходит её единственной дверью `of`; аффинная форма живёт в целых
        единицах сетки, и её поля — узлы. Чтобы сложить с УЗЛОМ, узел
        называется засечкой явно: `n + Notch(n.g, 1, 1)`."""
        a, b = _align(self, _as_notch(other, self.g))
        return Notch(a.g, a.lo + b.lo, a.hi + b.hi)

    __radd__ = __add__

    def __neg__(self):
        return Notch(self.g, -self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-_as_notch(other, self.g))

    def __rsub__(self, other):
        return (-self) + _as_notch(other, self.g)

    def __mul__(self, other):
        a, b = _align(self, _as_notch(other, self.g))
        scale = 1 << a.g
        c = (a.lo * b.lo, a.lo * b.hi, a.hi * b.lo, a.hi * b.hi)
        return Notch(a.g, _floor_div(min(c), scale), _ceil_div(max(c), scale))

    __rmul__ = __mul__

    def inverse(self) -> "Notch":
        if self.sign == 0:
            raise Unresolved("1/x", self.record(),
                             hint="засечка накрывает узел 0 — обращение "
                                  "на этой сетке не сведено; уточнить сетку")
        big = 1 << (2 * self.g)
        return Notch(self.g, _floor_div(big, self.hi), _ceil_div(big, self.lo))

    def __truediv__(self, other):
        return self * _as_notch(other, self.g).inverse()

    def __rtruediv__(self, other):
        return _as_notch(other, self.g) * self.inverse()

    def __pow__(self, n: int):
        if not isinstance(n, int) or isinstance(n, bool):
            return NotImplemented
        if n == 0:
            return Notch.from_q(1, self.g)
        if n < 0:
            return (self ** (-n)).inverse()
        out = Notch.from_q(1, self.g)
        for _ in range(n):
            out = out * self
        return out

    def root(self, n: int) -> "Notch":
        """Корень n-й степени целым счётом (монотонность); чётный n при
        отрицательном узле — отказ.

        Степень отбивается ПОИМЁННО, как у соседей по модулю (`__init__`:
        «сетка g — целое ≥ 0»; `__pow__`: bool показателем не считается;
        `from_q`: «bool не число»). Отдельная проверка `isinstance(n, bool)`
        нужна потому, что одно сравнение `n < 1` для `True` ложно — как для
        1: `root(True)` тихо считал корень первой степени и возвращал
        значение вместо отказа (§14 п. 9). Тем же условием закрыт и всякий
        нецелый тип степени: раньше он доходил до сдвига `1 << (g·(n−1))` и
        давал TypeError интерпретатора, а не отказ этой двери.
        """
        if not isinstance(n, int) or isinstance(n, bool) or n < 1:
            raise TypeError("степень корня — целое ≥ 1")
        if self.lo < 0 and n % 2 == 0:
            raise Unresolved("root", self.record(),
                             hint="чётный корень из засечки с отрицательным "
                                  "узлом не веществен")
        if self.hi <= 0 and n % 2 == 1:
            return -((-self).root(n))
        scale_pow = 1 << (self.g * (n - 1))          # x/2^g = r ⇒ r^n = x·2^{g(n−1)}/2^{gn}
        lo = int_root_floor(max(self.lo, 0) * scale_pow, n)
        t = self.hi * scale_pow
        r = int_root_floor(t, n)
        hi = r if r ** n == t else r + 1
        return Notch(self.g, lo, hi)


def _as_notch(x, g: int) -> Notch:
    return x if isinstance(x, Notch) else Notch.of(x, g)


def _align(a: Notch, b: Notch) -> tuple:
    g = max(a.g, b.g)
    return a.regrid(g), b.regrid(g)


# ---------- проверки гипотез на сетке: инструмент, не приговор ----------

def separate(x, y, g_max: int = 512, g0: int = 32):
    """Первая сетка, на которой засечки x и y расходятся; None — не
    разошлись до g_max (различие не доказано, равенство не утверждается)."""
    g = g0
    while g <= g_max:
        if Notch.of(x, g).separated(Notch.of(y, g)):
            return g
        g *= 2
    return None


def sign_on_grid(x, g: int) -> int:
    """Знак величины на сетке g: ±1 доказан, 0 — узел 0 внутри засечки."""
    return Notch.of(x, g).sign


# ---------- структурная запись: значение без линейки-носителя ----------

class Mixed:
    """Дерево действий над носителями: (op, left, right) либо лист.

    Запись точна (ничего не округлено); засечка — по запросу на любой
    сетке; равенство — структурное или на сетке; класс — 19 (предел как
    способ получения величины: предел засечек по сетке)."""

    __slots__ = ("op", "left", "right")
    _OPS = ("+", "-", "*", "/")

    def __init__(self, op, left, right=None):
        """Сырой конструктор ОБОРАЧИВАЕТ листья (§14 п. 4).

        Было: `Mixed('+', 1, 2).record()` → AttributeError интерпретатора
        («'int' object has no attribute 'record'»), и то же у `.notch(8)`, —
        тогда как докстринг класса описывает «(op, left, right) либо лист» и
        единственной дверью `combine`/`leaf` не называет. Из двух выходов —
        оборачивать или отказать с указанием двери — взят первый: докстринг
        тогда становится верен как написан, а не сужается, и сырой
        конструктор даёт РОВНО то же, что дверь (`Mixed('+',1,2) ==
        Mixed.combine(1,2,'+')`, хэши равны). Обёртка идемпотентна:
        `leaf` возвращает `Mixed` как есть, так что `combine`, который уже
        оборачивал, не оборачивает дважды.

        Неполная и противоречивая запись отказывает ПОИМЁННО, а не молча:
        действие без правого операнда раньше доходило до `record()` и падало
        там же AttributeError интерпретатора, а лист с правым операндом
        (`Mixed(None, 1, 2)`) молча забывал второй аргумент.
        """
        if op is not None and op not in self._OPS:
            raise ValueError(f"действие {op!r} вне {self._OPS}")
        if op is None and right is not None:
            raise ValueError(f"лист не имеет правого операнда: {right!r} — "
                             f"дверь действия: combine")
        if op is not None:
            if right is None:
                raise ValueError(f"действие {op!r} без правого операнда — "
                                 f"дверь: combine(a, b, op)")
            left, right = Mixed.leaf(left), Mixed.leaf(right)
        object.__setattr__(self, "op", op)
        object.__setattr__(self, "left", left)
        object.__setattr__(self, "right", right)

    def __setattr__(self, *_):
        raise AttributeError("структурная запись неизменяема")

    # Удаление — той же формой (§14 п. 11): `del m.op` отказывает тем же
    # словом, что `m.op = '-'`. Без этого запись теряла узел молча и падала
    # потом в чужом месте: у дерева без `op` отказывает уже не она, а
    # `record()`/`notch()` — AttributeError интерпретатора, как в п. 4.
    # Цена замерена и равна нулю — числа в докстринге неизменяемости `Affine`.
    __delattr__ = __setattr__

    @classmethod
    def leaf(cls, x) -> "Mixed":
        return x if isinstance(x, Mixed) else cls(None, x)

    @classmethod
    def combine(cls, a, b, op: str) -> "Mixed":
        return cls(op, cls.leaf(a), cls.leaf(b))

    @property
    def is_leaf(self) -> bool:
        return self.op is None

    # --- запись и структурное тождество -----------------------------------
    def record(self) -> str:
        if self.is_leaf:
            x = self.left
            r = getattr(x, "record", None)
            return r() if callable(r) else f"Φ(C=1, value={Q(x)})"
        return f"Φ(C=19 STRUCT, ({self.left.record()}) {self.op} ({self.right.record()}))"

    __repr__ = record

    def __eq__(self, other):
        if not isinstance(other, Mixed):
            return NotImplemented
        if self.is_leaf and other.is_leaf:
            return self.left == other.left
        return (self.op == other.op and self.left == other.left
                and self.right == other.right)

    def __hash__(self):
        return hash(("Mixed", self.record()))

    # --- дискретная сторона -------------------------------------------------
    def notch(self, g: int) -> Notch:
        """Засечка записи на сетке g: засечки листьев, арифметика на целых."""
        if self.is_leaf:
            return Notch.of(self.left, g)
        a, b = self.left.notch(g), self.right.notch(g)
        if self.op == "+":
            return a + b
        if self.op == "-":
            return a - b
        if self.op == "*":
            return a * b
        return a / b

    def enclosure(self, terms: int = 64) -> Interval:
        """Сертификат сетки в общем интерфейсе носителей: засечка на сетке
        g = terms, отрезком ⊂ ℚ."""
        return self.notch(terms).to_interval()

    def equal_on_grid(self, other, g: int) -> bool:
        return self.notch(g) == Mixed.leaf(other).notch(g)

    def sign_on_grid(self, g: int) -> int:
        return self.notch(g).sign

    # --- арифметика: только строением ---------------------------------------
    def __add__(self, other):
        return Mixed.combine(self, other, "+")

    def __radd__(self, other):
        return Mixed.combine(other, self, "+")

    def __sub__(self, other):
        return Mixed.combine(self, other, "-")

    def __rsub__(self, other):
        return Mixed.combine(other, self, "-")

    def __mul__(self, other):
        return Mixed.combine(self, other, "*")

    def __rmul__(self, other):
        return Mixed.combine(other, self, "*")

    def __truediv__(self, other):
        return Mixed.combine(self, other, "/")

    def __rtruediv__(self, other):
        return Mixed.combine(other, self, "/")

    def __neg__(self):
        return Mixed.combine(Q(0), self, "-")

    # --- сокращённая запись суммы: части и слияние однородных ---------------
    # Слово автора 10.10: правило колёс (запись есть сокращённая сумма по
    # разрядам) распространяется на все классы — разнородный итог хранится
    # как есть, последующие действия идут только с соответствующей частью,
    # приведение к однородной записи — по требованию (засечка, образ IEEE).
    # Поле порядка так и устроено (ω + 5 + ε − 36/7 → 1·ω − 1/7 + 1·ε,
    # учебник 09 §5); здесь то же для дерева класса 19.
    def parts(self, sign: int = 1, out: list | None = None) -> list:
        """Плоский список (знак, часть) по дереву сумм и разностей;
        произведение и частное — одной частью, как есть."""
        out = [] if out is None else out
        if self.is_leaf:
            out.append((sign, self.left))
        elif self.op == "+":
            self.left.parts(sign, out)
            self.right.parts(sign, out)
        elif self.op == "-":
            self.left.parts(sign, out)
            self.right.parts(-sign, out)
        else:
            out.append((sign, self))
        return out

    @staticmethod
    def _merge(a, b):
        """Слияние двух частей ИХ СОБСТВЕННЫМ оператором (ℚ с ℚ, q·π^e с
        тем же e, c·e^y с тем же y, …); отказ носителя — None, части
        остаются порознь. Дерево в ответе — не слияние."""
        try:
            r = a + b
        except (Unresolved, ArithmeticError, TypeError, ValueError):
            return None
        if r is NotImplemented or isinstance(r, Mixed):
            return None
        return r

    def collect(self) -> "Mixed":
        """Сокращённая запись суммы: однородные части слиты своим
        оператором, нулевые части сняты, порядок — первое появление;
        результат — левая свёртка частей. Произведения и частные не
        раскрываются (они — части). Идемпотентна: collect(collect(x)) ==
        collect(x)."""
        acc: list = []
        for sign, x in self.parts():
            if isinstance(x, Mixed):
                acc.append([sign, x])
                continue
            for item in acc:
                if isinstance(item[1], Mixed):
                    continue
                y = x if item[0] == sign else -x
                merged = Mixed._merge(item[1], y)
                if merged is not None:
                    # часть, слившаяся в ноль, снимается сразу — иначе она
                    # поглотит следующую чужую часть (x − x давал −1/3 + 1/3)
                    if Mixed._is_zero(merged):
                        acc.remove(item)
                    else:
                        item[1] = merged
                    break
            else:
                acc.append([sign, x])
        live = [(s, x) for s, x in acc if not Mixed._is_zero(x)]
        if not live:
            return Mixed.leaf(Q(0))
        s0, x0 = live[0]
        res = Mixed.leaf(x0) if s0 > 0 else Mixed.combine(Q(0), x0, "-")
        for s, x in live[1:]:
            res = Mixed.combine(res, x, "+" if s > 0 else "-")
        return res

    @staticmethod
    def _is_zero(x) -> bool:
        if isinstance(x, Mixed):
            return False
        try:
            return bool(x == 0) or bool(x == Q(0))
        except (Unresolved, ArithmeticError, TypeError, ValueError):
            return False
