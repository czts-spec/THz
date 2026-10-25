"""Аффинные засечки на целочисленной сетке — лечение проблемы зависимости.

Дефект из очереди (CONTEXT §4; GLOSSARY §1): интервальная арифметика по
Муру не помнит, что два интервала — одна величина: x − x даёт не ноль, а
удвоенный радиус (на засечках gridcore: ⟦−2; 2⟧ при узле ±1 — замер
04.09). Лечение известно и старо — аффинная арифметика (Comba & Stolfi
1993; Stolfi & de Figueiredo 1997): величина хранится как
x₀ + Σ xᵢ·εᵢ, εᵢ ∈ [−1, 1] — символы шума ОБЩИЕ для всех величин, потому
линейные операции точны по коэффициентам (x − x = 0 в точности), а
нелинейные добавляют новый символ с радиусом остатка. Здесь всё в целых
узлах сетки 2^−g: коэффициенты — целые единицы, округление коэффициентов
при смене сетки — наружу, в отдельный символ (сумма отброшенного не
превосходит числа округлённых коэффициентов — так соблюдается
включение). Засечка формы: [x₀ − Σ|xᵢ|, x₀ + Σ|xᵢ|] — то, что видит Мур;
разность форм видит больше.

Векторная часть (массивы (n, S, H): n образцов, S символов, H выходов)
— для сетей: слой — точное умножение коэффициентов на веса (ядро
exact_gemm), ReLU — преобразователь зонотопа (Singh, Gehr, Mirman,
Püschel, Vechev 2018, DeepZ): на смешанном отрезке [l, u] — y = λx + μ +
μ·ε_new, λ = u/(u−l), μ = −λl/2; сертификат запаса — аффинная форма
разности логитов, а не разность границ.

Авторство: Comba–Stolfi 1993, Stolfi–de Figueiredo 1997 (аффинная
арифметика), Singh et al. 2018 (зонотопы для ReLU-сетей), Мур 1966
(интервалы), Rump 1999 (шаровая арифметика — то, что лечится); засечки
как узлы сетки — AIAΣA (автор проекта).

Обратный элемент и деление (04.09): нелинейная функция 1/x приближается
аффинной y = αx + ζ + δ·ε_new на засечке формы — два канонических выбора
(Stolfi–de Figueiredo 1997, §3.4). **Чебышёв** (минимакс): α = −1/(ab),
касание в u = √(ab), δ минимален из возможных. **min-range**: α = f′(b) =
−1/b², засечка итога не шире образа отрезка — ценой большего δ. Здесь оба
считаются точной рациональной арифметикой, √(ab) берётся парой целых
границ, а всё округление к сетке уходит НАРУЖУ, в тот же новый символ:
включение соблюдается по построению, а не по надежде. Засечка, накрывающая
ноль, — отказ поимённо (`Unresolved`), не бесконечность.

Корень √x взят тем же приёмом, но δ у него считается СТРОГО, а не по
формуле: α и ζ берутся рациональными приближениями (корень концов
иррационален), а затем e(x) = √x − αx − ζ ограничивается по концам и по
точке касания x* = 1/(4α²) — она рациональна, потому что рационально α.
Неточность приближения целиком уходит в δ и ничего не ломает. Засечка,
заходящая левее нуля, — отказ поимённо: вещественного корня там нет.

Ширина зонотопов замерена до D = 2048 (`proto/bench_zono_width.py`,
04.09): пик массива коэффициентов — тождество n·D²·8 байт (входной
зонотоп), символов D + 3, включение — 0 нарушений на каждом размере.

DeepPoly замерен 04.09 (`phicore/deeppoly.py`, `proto/bench_deeppoly.py`):
на сети с одним скрытым слоем он сертифицирует ровно столько же, сколько
зонотоп, при более широких границах логитов и цене в двести раз выше —
подстановка назад через один слой ReLU не даёт ничего сверх общих
символов шума.

NOT_MEASURED — разобран 12.09 на два разных утверждения (строка остаётся
записью о том, что было):

  * **пик памяти — ЗАКРЫТО точной дорогой; что осталось, есть ОТКАЗ СРЕДЫ.**
    Пик не оценивается, а СЧИТАЕТСЯ: массив коэффициентов входного зонотопа
    есть (n, D, D) целых по 8 байт, то есть ровно n·D²·8 байт (`zono_bytes`;
    сверено с `C.nbytes` самого `zono_input` — совпало на пяти парах (n, D)).
    При n = 32 и D = 4096 это 4 294 967 296 байт — РОВНО 4 ГиБ, без округления.
    Поправка замером `[вычислено здесь]`: «память процесса втрое больше» —
    число не отсюда; замеренная здесь лестница прироста RSS над базой даёт
    кратность 6,34 при (n=4, D=512), 5,20 при (4, 1024) и 4,80 при (8, 1024),
    то есть кратность УБЫВАЕТ с размером и втрое не бывает. По нижней
    замеренной кратности 4,80 пик процесса при (32, 4096) есть ≈ 19,2 ГиБ
    против MemTotal 15,7 ГиБ этой среды — не помещается, и это счёт, а не
    осторожность. Порог назван обоими числами и обязан перемеряться в другой
    среде, а не переноситься словом.

  * **приближения выше первой степени — ОТКАЗ ДОМА, названный ценой.**
    Предмет: Тейлор второго порядка для 1/x требует коэффициента при паре
    символов (εᵢ·εⱼ), а носитель формы несёт коэффициенты ЛИНЕЙНО —
    `Affine.__slots__` хранит словарь символ → коэффициент, места под пару в
    нём нет по построению. Причина: не арифметика не сходится, а носителя
    квадратичной части дом не несёт; произведение двух аффинных форм аффинным
    не бывает — ровно поэтому здесь и стоит приём нового символа.
    Чем мерилось бы: переход от D коэффициентов на форму к D·(D+1)/2, то есть
    от n·D к n·D²/2 в памяти и к тому же множителю в цене слоя; выигрыш —
    только на сильно нелинейных участках, где δ первой степени широка.
    Это «не посчитано здесь», а не «посчитать нельзя».
"""
from __future__ import annotations

import itertools
import math
from numbers import Integral
from types import MappingProxyType

import numpy as np

from .errors import Unresolved
from .exact_gemm import INF, exact_gemm
from .gridcore import Notch
from .rational import Q, nearest_int

_fresh = itertools.count(1)
_reserved = 0

# Одна запись отказа на все двери формы (§14 п. 7). Двух копий здесь нет
# намеренно: разошедшиеся копии одного канона — уже пройденная поломка
# (`Notch.of` против `Notch.from_q`, §14 п. 10), и разойтись им нельзя.
_NOT_AFFINE = "не аффинная форма и не целое: "
#: `bool` числом в доме не считается (канон `phicore/rational.py`; сосед по
#: главе отвечает тем же: `Notch.from_q(True, 8)` → «bool не число»,
#: `Notch.__pow__(True)` → NotImplemented). Замер 12.09: `A + True` давало
#: узел 385, `A * True` — множитель 1, а `A / True` уже отказывало — один вход
#: получал два ответа в одном модуле. Решено не переадресацией, а КАНОНОМ:
#: bool отбит во всех дверях формы одной записью.
_NOT_BOOL = "bool не число"
#: Поле записи формы — ЦЕЛАЯ единица сетки. Запись отдельная от `_NOT_AFFINE`
#: намеренно: та об ОПЕРАНДЕ действия («не аффинная форма и не целое»), эта — о
#: ПОЛЕ записи, и путать их у вызывающего нельзя. Замер 12.09 (глава 08 §14
#: п. 3, дыра, найденная прогоном главы): конструктор усекал всё, что
#: приводится к int, — `Affine(8, 3.9, {1: 2.9})` давал `(3 + 2·ε1)·2^−8`,
#: `Affine(8.7, 1).g` → 8, `Affine(8, True).c0` → 1, `Affine(8, '5').c0` → 5,
#: тогда как соседняя дверь того же модуля отказывала: `A / 1.5` → TypeError.
_NOT_UNIT = "не целая единица сетки: "


def _unit(v, name: str) -> int:
    """Поле записи формы — целое (в единицах сетки 2^−g); прочее отвергается
    ПОИМЁННО, а не усекается.

    Канон дома — `phicore/rational.py` («float/bool на входе — TypeError;
    единственная дверь из мира IEEE — явный `Q.from_float`») и соседи по
    главе: `Notch.__init__` («узлы засечки — целые»), `Notch.from_q` («float
    запрещён построением»), `exact_gemm._ints` («float запрещён построением
    (усечение молча меняет предмет)»). Целые numpy проходят (зонотоп считает
    ими), bool — нет: он подкласс int, и канон дома числом его не считает."""
    if isinstance(v, bool):
        raise TypeError(f"{name}: {_NOT_BOOL}")
    if isinstance(v, Integral):
        return int(v)
    raise TypeError(f"{name}: {_NOT_UNIT}{v!r}")


def fresh_symbol() -> int:
    """Новый символ шума. Пространство имён символов ОБЩЕЕ для всех форм:
    в этом весь смысл аффинной арифметики. Потому счётчик обязан обходить и
    те номера, что пришли в форму со стороны — иначе новый символ совпадёт
    с чужим, коэффициенты сложатся и включение молча сломается (замер 04.09:
    δ обратного элемента обнулялся, 1/x переставал накрывать истинное)."""
    while True:
        k = next(_fresh)
        if k > _reserved:
            return k


def _reserve(ids) -> None:
    global _reserved
    for k in ids:
        if k > _reserved:
            _reserved = k


class Affine:
    """Скалярная аффинная форма на сетке 2^−g: x₀ + Σ xᵢ εᵢ, целые единицы.

    ЗАКОН КОЛЕСА (слово автора 09.09; `phicore/wheels.py`): центр `c0·2^−g`
    и каждый коэффициент `xᵢ·2^−g` суть контейнеры (c, 1, −g, 1) над α = 2 —
    по одному колесу базы 2 без ширины, в сетке (2, 1) — сумма колёс по
    единичным битам; форма есть центр-колесо плюс колёса при символах шума.
    Сдвиги `<< (g − g')` и деления на 2^(g − g') здесь — БЫСТРАЯ ТРОПА к тому
    же числу (путь зонотопа и `regrid` горячие), закон сверяется приёмкой
    `proto/test_phicore_affine.py` с `wheels.Record.of(Q(c, 2**g), 2)`
    перебором и с подлогом; мост наружу — `wheel_records`. Сосед по закону —
    засечка `gridcore.Notch` (те же контейнеры над α = 2, два узла вместо
    центра и радиуса).
    """
    __slots__ = ("g", "c0", "_coef", "_view")

    def __init__(self, g: int, c0: int, coef: dict | None = None):
        """Форма в целых единицах сетки 2^−g; всякое поле — целое.

        Дверь закрыта поимённо с 12.09 (§14 п. 3 главы 08): прежде здесь
        стояло `int(...)` на каждом поле, и запись молча меняла предмет —
        `Affine(8, 3.9, {1: 2.9})` становилась `(3 + 2·ε1)·2^−8` (замер ДО),
        `Affine(8.7, 1).g` — 8, `Affine(8, True).c0` — 1, `Affine(8, '5').c0`
        — 5. Теперь каждое поле проходит `_unit` и отвечает той же записью,
        что соседи по дому; целое (в том числе целое numpy) проходит как
        прежде, и ни одна внутренняя сборка формы не изменилась.
        ЦЕЛОЕ ЗДЕСЬ ЕСТЬ УЗЕЛ, а не величина: `Affine(8, 1)` — это 2^−8,
        тогда как у соседа `Notch` целое входит ВЕЛИЧИНОЙ (`3 + Notch.from_q(1,4)`
        → `⟦64; 64⟧·2^−4`, то есть 4) — §14 п. 8, оговорено в `__add__`."""
        g = _unit(g, "Affine g")
        c0 = _unit(c0, "Affine c0")
        d = {}
        for k, v in (coef or {}).items():
            key = _unit(k, "Affine символ шума")
            val = _unit(v, f"Affine коэффициент ε{key}")
            if val != 0:
                d[key] = val
        object.__setattr__(self, "g", g)
        object.__setattr__(self, "c0", c0)
        object.__setattr__(self, "_coef", d)
        _reserve(d)

    @property
    def coef(self):
        """Коэффициенты наружу — неизменяемым видом, и ни одного лишнего хода.

        Вид строится ЛЕНИВО и запоминается, а внутри модуля счёт идёт по самому
        словарю (`_coef`): за неизменяемость словаря платит тот, кто словарь
        читает, а не вся арифметика. Замер встречным прогоном старого и нового
        модуля вперемешку, минимум из девяти ходов (12.09, машина под нагрузкой
        10 на 4 ядра): `radius` 20 000× — 8 783 против 8 804 мкс, отношение
        1,00; 3000×(+, ·, −, radius) — 16 069 против 20 928 мкс, 1,30 (надбавка
        целиком из конструктора, см. `__setattr__`). С видом, построенным в
        конструкторе, то же место стоило 23 211 мкс, 1,48 — ленивость сняла
        треть надбавки.
        """
        try:
            return self._view
        except AttributeError:
            object.__setattr__(self, "_view", MappingProxyType(self._coef))
            return self._view

    def __setattr__(self, *_):
        """Форма неизменяема — той же формой, что соседи по главе (§14 п. 6).

        Было: `mut.c0 = 5` проходило и меняло уже выданную наружу форму
        (замер 12.09: `(256 + 128·ε1)·2^−8` → `(5 + 128·ε1)·2^−8`), тогда как
        `Notch` отвечает «засечка неизменяема», а `Mixed` — «структурная
        запись неизменяема». Здесь изменяемость опаснее, чем у соседей: форма
        хранит номера ОБЩИХ символов шума, и подмена коэффициента задним
        числом рвёт включение у всех форм, с которыми она эти символы делит, —
        молча и в другом месте.

        Словарь коэффициентов закрыт вместе с полями, а не только ссылка на
        него: `mut.coef[99] = 7` тоже проходило и добавляло символ радиусом 7
        (замер 12.09: `(5 + 128·ε1)·2^−8` → `(5 + 128·ε1 + 7·ε99)·2^−8`), то
        есть объявить форму неизменяемой, оставив словарь открытым, значило бы
        объявить неверное. `coef` отдаётся неизменяемым видом
        (`MappingProxyType`) поверх того же словаря: чтение, `.items()`,
        `dict(...)`, `set(...)`, `max(...)`, `not coef` работают как раньше и
        ничего не копируют, а запись отказывает сам вид — `coef[99] = 7` →
        TypeError «'mappingproxy' object does not support item assignment»,
        `del coef[1]` → TypeError «doesn't support item deletion».

        Удаление поля закрыто ТОЙ ЖЕ формой, а не второй копией отказа:
        `__delattr__ = __setattr__`, один объект функции на две двери — и так
        же у обоих соседей по главе (§14 п. 11). Было: `del form.c0`
        проходило, хотя `form.c0 = 5` уже отбивалось, и оставляло выданную
        наружу форму полусобранной — `record()` после этого падал
        AttributeError интерпретатора («'Affine' object has no attribute
        'c0'»), то есть отказывала уже не форма и не в том месте, где её
        сломали. Опаснее присваивания здесь то же самое: `del form._coef`
        уносит номера ОБЩИХ символов шума. Копия строки отказа расходится —
        разнобой §14 пп. 9, 10 был именно в ВИДЕ отказа, — а при одном
        объекте функции `del` ответить иначе, чем присваивание, физически не
        может: у `Notch` — «засечка неизменяема», у `Mixed` — «структурная
        запись неизменяема», здесь — «аффинная форма неизменяема».

        Цена закрытия удаления замерена и равна нулю: форма не добавляет ни
        одной команды ни на один путь (одно имя в словаре класса при
        импорте), `__slots__` не тронуты, размер экземпляра тот же — 56, 56 и
        64 байта у `Notch`, `Mixed`, `Affine` до и после. Встречный прогон
        старого и нового модуля вперемешку, минимум из девяти ходов, два
        прогона 12.09 (машина под нагрузкой 2,9 и 5,1 на 4 ядра), 20 000×:
        конструктор `Notch` 14 435 → 14 766 и 14 488 → 14 689 мкс (1,02 и
        1,01); `Notch` «+» 23 742 → 24 071 и 23 721 → 23 995 (1,01 и 1,01);
        конструктор `Mixed` 42 472 → 42 618 и 41 905 → 42 219 (1,00 и 1,01);
        `Mixed` «==» 12 938 → 12 972 и 12 946 → 13 019 (1,00 и 1,01);
        конструктор `Affine` 24 072 → 23 874 и 23 990 → 23 865 (0,99 и 0,99);
        `Affine.radius` 7 573 → 7 514 и 7 498 → 7 635 (0,99 и 1,02). Знак
        отклонения меняется между прогонами — это шум машины, а не цена.

        Цена неизменяемости замерена и названа: вся она в конструкторе, где
        `object.__setattr__` идёт вместо прямого присваивания слота — 20 000
        построений формы 15 820 → 23 448 мкс, отношение 1,48, то есть 0,38 мкс
        на форму (замер 12.09 встречным прогоном старого и нового модуля,
        минимум из девяти ходов). Ту же цену той же строкой платят соседи по
        главе: `Notch.__init__` и `Mixed.__init__` пишут поля так же. На тех
        путях, где форма не строится в цикле, цены нет: путь зонотопа (numpy,
        `zono_input`→`affine`→`relu`→`margin`) — 2 474 против 2 505 мкс, 1,01
        при совпадении ответов узел в узел; нелинейное (40×(div, 1/x, √x)) —
        12 374 против 12 535 мкс, 1,01.
        """
        raise AttributeError("аффинная форма неизменяема")

    __delattr__ = __setattr__

    @classmethod
    def from_notch(cls, n: Notch) -> "Affine":
        """Новая величина из засечки — новый символ (радиус — половина ширины, наружу)."""
        mid2 = n.lo + n.hi
        if mid2 % 2 == 0:
            return cls(n.g, mid2 // 2, {fresh_symbol(): (n.hi - n.lo) // 2})
        # нечётная ширина: центр между узлами — расширяем наружу на единицу
        return cls(n.g, (mid2 - 1) // 2, {fresh_symbol(): (n.hi - n.lo + 1) // 2})

    @classmethod
    def exact(cls, g: int, units: int) -> "Affine":
        """Форма без символов шума: величина `units·2^−g` ТОЧНО (радиус 0).

        Одноимённые двери соседей — ДРУГИЕ понятия, и это не одно устройство:
        `ieee.exact(x)` берёт точную сторону ОДНОГО float (пару из
        `as_integer_ratio`, ±ω, состояние для nan), `window.Window.exact(x)`
        — показание одометра целиком (образ плюс поглощённое, `Q`). Здесь же
        `exact` — конструктор формы из узла сетки: целое есть УЗЕЛ, не
        величина (§14 п. 8, `__add__`).
        """
        return cls(g, units)

    def wheel_records(self, width=1) -> dict:
        """Центр и коэффициенты как записи колёсами базы 2: {0: центр, k: xₖ}.

        Мост к закону (`wheels.Record.of`), не тропа: узлы `c·2^−g` лежат в
        ℤ[1/2] построением, отказа нет; `width=1` — сетка (2, 1), `None` —
        одно колесо. Ввоз ленивый: `wheels` ввозит `container`, форма живёт
        без него.
        """
        from .wheels import Record
        scale = 1 << self.g
        out = {0: Record.of(Q(self.c0, scale), 2, width)}
        for k, v in self._coef.items():
            out[k] = Record.of(Q(v, scale), 2, width)
        return out

    @property
    def radius(self) -> int:
        return sum(abs(v) for v in self._coef.values())

    def notch(self) -> Notch:
        r = self.radius
        return Notch(self.g, self.c0 - r, self.c0 + r)

    def _same_grid(self, other: "Affine") -> None:
        """Общая дверь сложения и умножения: сначала ТИП, потом сетка (§14 п. 7).

        Было: тип не проверялся вовсе, и `A * 1.5` падало AttributeError
        интерпретатора («'float' object has no attribute 'g'») — отказом не
        этой двери и даже не этого модуля, — тогда как `A / 1.5` отвечало
        TypeError поимённо («не аффинная форма и не целое: 1.5»). Один вид
        отказа на оба действия получен здесь, а не в `__mul__`: именно тут
        обрывался путь обоих, и тем же одним местом закрыт `A + 1.5` (то же
        AttributeError интерпретатора) и правые формы `1.5 * A`, `1.5 + A`.
        Слова взяты не копией, а той же записью `_NOT_AFFINE`, что у
        `_as_affine`, — разойтись им нельзя.
        """
        if not isinstance(other, Affine):
            raise TypeError(_NOT_AFFINE + repr(other))
        if self.g != other.g:
            raise ValueError(f"разные сетки: 2^-{self.g} и 2^-{other.g}")

    def __add__(self, other):
        """Сумма форм на одной сетке; ЦЕЛОЕ ЗДЕСЬ — УЗЕЛ, не величина.

        Оговорка 12.09 (§14 п. 8 главы 08: «целое значит разное в двух модулях
        одной главы, и ни один докстринг этого не оговаривает»). Замер, одна
        и та же единица у двух соседей:
          * форма: `A + 1` при `A = Affine.from_notch(⟦256; 512⟧·2^−8)` даёт
            `(385 + 128·ε)·2^−8` — прибавлен ОДИН УЗЕЛ сетки, то есть 2^−8;
          * засечка: `3 + Notch.from_q(1, 4)` даёт `⟦64; 64⟧·2^−4`, то есть 4 —
            прибавлена ВЕЛИЧИНА 3 (дверь `Notch.of`).
        Разница намеренная и следует из предмета: форма живёт В ЦЕЛЫХ ЕДИНИЦАХ
        сетки (докстринг класса), и её поля — узлы; засечка несёт величину, и
        её дверь `of` принимает носителя. Чтобы прибавить к форме ВЕЛИЧИНУ,
        дверь названа явно: `A + Affine.exact(A.g, 1 << A.g)` (величина 1) либо
        `Affine.from_notch(Notch.of(x, A.g))` для произвольного носителя.
        То же правило действует у `__sub__`, `__mul__` (целое — множитель, не
        узел) и у `scale`."""
        if isinstance(other, bool):
            raise TypeError(_NOT_BOOL)
        if isinstance(other, int):
            return Affine(self.g, self.c0 + other, self._coef)
        self._same_grid(other)
        c = dict(self._coef)
        for k, v in other._coef.items():
            c[k] = c.get(k, 0) + v
        return Affine(self.g, self.c0 + other.c0, c)

    __radd__ = __add__

    def __neg__(self):
        return Affine(self.g, -self.c0, {k: -v for k, v in self._coef.items()})

    def __sub__(self, other):
        """Разность: та же дверь, что у сложения (§14 п. 7, попутно).

        Было `-int(other)` — третье, самое широкое определение «целого» в
        модуле: `A - 1.5` молча давало `A - 1`, `A - '5'` — `A - 5` (замеры
        12.09: `(255 + 128·ε1)·2^−8` и `(251 + 128·ε1)·2^−8`). Молчаливое
        усечение на публичном действии хуже отказа: ответ неверен и ничем не
        помечен. Целое здесь понимается ровно как в `__add__` (`bool` — как 1,
        чтобы `A - True` не разошлось с `A + True`), остальное уходит в общую
        дверь и отказывает поимённо. Усечение в КОНСТРУКТОРЕ (§14 п. 3) этой
        правкой не закрыто: другая дверь, другой пункт очереди.
        """
        if isinstance(other, Affine):
            return self + (-other)
        if isinstance(other, bool):
            raise TypeError(_NOT_BOOL)
        if isinstance(other, int):
            return self + (-other)
        raise TypeError(_NOT_AFFINE + repr(other))      # та же запись, не копия

    def __rsub__(self, other):
        return (-self) + other

    def scale(self, k: int) -> "Affine":
        """Умножение на целое — точно; `bool` числом не считается (канон дома)."""
        if isinstance(k, bool):
            raise TypeError(_NOT_BOOL)
        return Affine(self.g, self.c0 * k, {s: v * k for s, v in self._coef.items()})

    def __mul__(self, other):
        """Произведение форм: линейная часть точна, остаток Σxᵢεᵢ·Σyᵢεᵢ —
        новый символ радиуса rx·ry (Stolfi–de Figueiredo). Итог на сетке
        2^−(g+g): единицы перемножаются, сетка складывается."""
        if isinstance(other, bool):
            raise TypeError(_NOT_BOOL)
        if isinstance(other, int):
            return self.scale(other)
        self._same_grid(other)
        c: dict = {}
        for k, v in self._coef.items():
            c[k] = c.get(k, 0) + other.c0 * v
        for k, v in other._coef.items():
            c[k] = c.get(k, 0) + self.c0 * v
        rem = self.radius * other.radius
        if rem:
            c[fresh_symbol()] = rem
        return Affine(self.g + other.g, self.c0 * other.c0, c)

    __rmul__ = __mul__

    def regrid(self, g: int) -> "Affine":
        """Смена сетки с сохранением включения: коэффициенты к нулю, отброшенное
        — в новый символ (не больше числа округлённых членов плюс центр).

        Одно понятие с `wheels.Record.regrid` и `gridcore.Notch.regrid` —
        та же величина в другой сетке; расходятся трое там, где величина в
        новую сетку не входит. Тоньше (g ≥ self.g) — у всех точно: центр и
        коэффициенты умножаются на 2^(g − self.g), и запись колёсами базы 2
        у каждого узла ТА ЖЕ (`Record.of(Q(c, 2**g), 2, 1)` от g не зависит).
        Грубее — здесь каждый коэффициент режется К НУЛЮ и недостача уходит
        новому символу шума (форма ШИРЕ, включение соблюдено); `Notch.regrid`
        округляет оба узла НАРУЖУ без нового символа; `Record.regrid` не
        округляет никогда — либо точен, либо отказ по имени: запись держит
        одно число, а форма и засечка — множества. Приёмка сверяет все три.
        """
        if g >= self.g:
            f = 1 << (g - self.g)
            return Affine(g, self.c0 * f, {k: v * f for k, v in self._coef.items()})
        d = 1 << (self.g - g)
        c0 = self.c0 // d if self.c0 >= 0 else -((-self.c0) // d)
        coef = {}
        for k, v in self._coef.items():
            coef[k] = v // d if v >= 0 else -((-v) // d)
        budget = 1 + sum(1 for v in self._coef.values() if v % d != 0) + (1 if self.c0 % d else 0)
        coef[fresh_symbol()] = budget
        return Affine(g, c0, coef)

    # ── нелинейное: обратный элемент и деление ────────────────────────

    def _value(self, units: int) -> Q:
        return Q(units, 1 << self.g)

    def _affine_of(self, alpha: Q, zeta: Q, delta: Q, g_out: int) -> "Affine":
        """y = α·self + ζ ± δ на сетке 2^−g_out.

        Округление коэффициентов — к ближайшему узлу, а вся недостача
        (центр, коэффициенты, δ) — наружу, в новый символ. Тогда для любого
        набора εᵢ истинное значение представимо той же формой: расхождение
        не больше радиуса нового символа.
        """
        scale = 1 << g_out
        y0 = (alpha * self._value(self.c0) + zeta) * scale
        c0 = nearest_int(y0)
        slack = abs(y0 - c0)
        coef = {}
        for k, v in self._coef.items():
            yi = alpha * self._value(v) * scale
            u = nearest_int(yi)
            slack += abs(yi - u)
            if u:
                coef[k] = u
        need = delta * scale + slack
        d_units = -((-need.numerator) // need.denominator)      # ceil, наружу
        if d_units:
            coef[fresh_symbol()] = d_units
        return Affine(g_out, c0, coef)

    @staticmethod
    def _sqrt_bounds(x: Q, bits: int = 96) -> tuple:
        """Пара рациональных границ √x: lo ≤ √x ≤ hi, x > 0.

        Когда корень берётся точно (числитель и знаменатель — полные
        квадраты, в том числе при a = b), границы совпадают и приближение
        не платит за незнание корня ничего.
        """
        rn, rd = math.isqrt(x.numerator), math.isqrt(x.denominator)
        if rn * rn == x.numerator and rd * rd == x.denominator:
            exact = Q(rn, rd)
            return exact, exact
        scale = 1 << bits
        y = x * scale * scale
        num, den = y.numerator, y.denominator
        lo = Q(math.isqrt(num // den), scale)
        hi = Q(math.isqrt(-((-num) // den)) + 1, scale)
        return lo, hi

    def _default_out_grid(self) -> int:
        """Сетка обратного по умолчанию: столько же значащих бит, сколько у x.

        `[вычислено здесь]` Оставлять сетку прежней нельзя: у величины
        порядка 2^m обратная имеет порядок 2^−m, и на сетке 2^−g от неё
        остаётся g − m бит вместо g + m. Замер 04.09 на ⟦63; 65⟧, сетка
        2^−6: x/x давало ширину 2,06 против 0,063 у Мура — аффинная запись
        проигрывала интервальной вчистую, и не из-за зависимости, а из-за
        сетки. При g_out = g + 2m счёт значащих бит сохраняется.
        """
        m = max(abs(self.c0), 1).bit_length() - 1 - self.g
        return max(0, self.g + 2 * m)

    def reciprocal(self, g_out: int | None = None, mode: str = "chebyshev") -> "Affine":
        """1/x аффинной формой: Чебышёв (минимакс) или min-range.

        Отказ, а не бесконечность, если засечка накрывает ноль: обратного
        на ней нет, и правдоподобное число здесь было бы подменой.
        """
        if mode not in ("chebyshev", "minrange"):
            raise ValueError(f"неизвестное приближение: {mode!r}")
        g_out = self._default_out_grid() if g_out is None else int(g_out)
        n = self.notch()
        if n.lo <= 0 <= n.hi:
            raise Unresolved("1/", self.record(),
                             hint=f"засечка ⟦{n.lo}; {n.hi}⟧·2^−{n.g} накрывает 0: "
                                  f"обратного элемента на ней нет")
        if n.hi < 0:                       # 1/x = −(1/(−x)); слева от нуля
            return -((-self).reciprocal(g_out, mode))
        a, b = self._value(n.lo), self._value(n.hi)
        if mode == "minrange":
            alpha = -1 / (b * b)
            d_max = 1 / a + a / (b * b)    # разность в левом конце
            d_min = 2 / b                  # и в правом; f − αx убывает
            zeta = (d_max + d_min) / 2
            delta = (d_max - d_min) / 2
        else:
            alpha = -1 / (a * b)
            half = (1 / a + 1 / b) / 2
            u_lo, u_hi = self._sqrt_bounds(a * b)
            zeta = half + 1 / u_hi         # чуть ниже точного ζ* = half + 1/√(ab)
            delta = half - 1 / u_hi + (1 / u_lo - 1 / u_hi)   # δ* сверху + сдвиг ζ
        return self._affine_of(alpha, zeta, delta, g_out)

    def sqrt(self, g_out: int | None = None, mode: str = "chebyshev") -> "Affine":
        """√x аффинной формой: тот же выбор — минимакс или min-range.

        Здесь α и ζ берутся рациональными приближениями (корень концов
        иррационален), а δ считается СТРОГО: e(x) = √x − αx − ζ вогнута,
        её максимум — в x* = 1/(4α²), и точка эта рациональна; минимум —
        на конце отрезка. Потому неточность α и ζ ничего не ломает: она
        целиком уходит в δ, а δ доказан, а не оценён на глаз.
        """
        if mode not in ("chebyshev", "minrange"):
            raise ValueError(f"неизвестное приближение: {mode!r}")
        g_out = self.g if g_out is None else int(g_out)
        n = self.notch()
        if n.lo < 0:
            raise Unresolved("√", self.record(),
                             hint=f"засечка ⟦{n.lo}; {n.hi}⟧·2^−{n.g} заходит левее нуля: "
                                  f"вещественного корня на ней нет")
        a, b = self._value(n.lo), self._value(n.hi)
        sa_lo, sa_hi = (Q(0), Q(0)) if a == 0 else self._sqrt_bounds(a)
        sb_lo, sb_hi = self._sqrt_bounds(b)
        if b == 0:
            return self._affine_of(Q(0), Q(0), Q(0), g_out)
        if mode == "minrange":
            alpha = 1 / (2 * sb_hi)          # ≤ f′(b): g = √x − αx растёт
            g_a_lo, g_a_hi = sa_lo - alpha * a, sa_hi - alpha * a
            g_b_lo, g_b_hi = sb_lo - alpha * b, sb_hi - alpha * b
            zeta = (g_a_lo + g_b_hi) / 2
        else:
            alpha = 1 / (sa_hi + sb_hi)      # ≤ (√b−√a)/(b−a)
            top = 1 / (4 * alpha)            # f(u) − αu в точке касания u = 1/(4α²)
            zeta = (sa_lo - alpha * a + top) / 2
        # строгая оценка δ: вогнутая e(x) = √x − αx − ζ
        u = 1 / (4 * alpha * alpha)
        highs = [sa_hi - alpha * a - zeta, sb_hi - alpha * b - zeta]
        lows = [sa_lo - alpha * a - zeta, sb_lo - alpha * b - zeta]
        if a <= u <= b:
            highs.append(1 / (4 * alpha) - zeta)
        delta = max(max(highs), -min(lows), Q(0))
        return self._affine_of(alpha, zeta, delta, g_out)

    def _as_affine(self, other) -> "Affine":
        if isinstance(other, Affine):
            return other
        if isinstance(other, bool):
            raise TypeError(_NOT_BOOL)                  # та же запись, что у прочих дверей
        if isinstance(other, int):
            return Affine.exact(0, other).regrid(self.g) if self.g else Affine.exact(0, other)
        raise TypeError(_NOT_AFFINE + repr(other))

    def div(self, other, mode: str = "chebyshev") -> "Affine":
        """x/y = x·(1/y): символы шума общие, потому зависимость видна.

        Сетка итога — сумма сеток, как у умножения.
        """
        other = self._as_affine(other)
        r = other.reciprocal(mode=mode)
        g = max(self.g, r.g)          # вверх по сетке — точно, вниз бы теряло
        return self.regrid(g) * r.regrid(g)

    def __truediv__(self, other):
        return self.div(other)

    def __rtruediv__(self, other):
        return self._as_affine(other).div(self)

    # ── самопроверка формы (ось самоприменения, 09.10) ───────────────

    def self_check(self, deep: bool = False) -> list:
        """Самопроверка формы встречными дорогами (ось «AIAΣA спрашивает
        AIAΣA», слово автора 04.09; очередь строка 18; 09.10). Список
        находок, пусто — здорова; исключений не бросает.

        `[вычислено здесь]` (1) слоты: g ≥ 0, c0 и коэффициенты — целые,
        нулевых коэффициентов нет (конструктор их не хранит); (2) радиус
        есть сумма модулей коэффициентов, засечка формы — центр ± радиус;
        (3) ВСТРЕЧНЫЕ дороги: лечение зависимости — разность формы с собой
        есть точный ноль без символов (то, ради чего форма существует), а по
        Муру та же разность имела бы ширину 2·радиус; смена сетки — тоньше
        точна (узлы засечки умножаются на 2^k), обратно грубее включает
        исходную засечку; закон колеса — центр и коэффициенты как записи
        колёсами базы 2 совпадают с `wheels.Record.of` узла; обратный элемент
        (засечка вне нуля) — произведение формы на свою обратную форму несёт
        единицу на сетке произведения: умножение форм против приближения
        Чебышёва; корень (засечка правее нуля) — квадрат формы корня несёт
        засечку формы. При deep — ещё min-range у обоих."""
        out = []
        g, c0, coef = self.g, self.c0, self._coef
        if not isinstance(g, int) or isinstance(g, bool) or g < 0:
            return [f"слот g вне канона: {g!r}"]
        if not isinstance(c0, int) or isinstance(c0, bool):
            return [f"слот c0 не целое: {c0!r}"]
        if any(not isinstance(v, int) or isinstance(v, bool) or v == 0 for v in coef.values()):
            return ["коэффициенты формы: не целое либо нуль в словаре"]
        r = sum(abs(v) for v in coef.values())
        if self.radius != r:
            out.append(f"радиус {self.radius} ≠ сумма модулей коэффициентов {r}")
        n = self.notch()
        if (n.g, n.lo, n.hi) != (g, c0 - r, c0 + r):
            out.append(f"засечка формы ⟦{n.lo}; {n.hi}⟧·2^−{n.g} ≠ центр ± радиус")
        try:
            z = self - self
            if z.c0 != 0 or z.radius != 0:
                out.append(f"лечение зависимости: x − x = {z.record()}, а не точный ноль")
            if (self + self).radius != 2 * r:
                out.append("x + x: радиус не удвоен")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"разность с собой: {type(e).__name__}: {e}")
        try:
            fine = self.regrid(g + 4)
            nf = fine.notch()
            if (nf.lo, nf.hi) != (n.lo * 16, n.hi * 16) or fine.radius != 16 * r:
                out.append("смена сетки тоньше не точна")
            back = fine.regrid(g)
            nb = back.notch()
            if nb.lo > n.lo or nb.hi < n.hi:
                out.append(f"смена сетки грубее не включает исходную засечку: ⟦{nb.lo}; {nb.hi}⟧ ⊉ ⟦{n.lo}; {n.hi}⟧")
            if g >= 2:
                coarse = self.regrid(g - 2).notch()
                if coarse.lo * 4 > n.lo or coarse.hi * 4 < n.hi:
                    out.append("смена сетки на две грубее не включает исходную засечку")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"смена сетки: {type(e).__name__}: {e}")
        try:
            from .wheels import Record
            recs = self.wheel_records()
            if recs[0] != Record.of(Q(c0, 1 << g), 2, 1) or any(recs[k] != Record.of(Q(v, 1 << g), 2, 1) for k, v in coef.items()):
                out.append("закон колеса: записи центра или коэффициентов не совпали с Record.of")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"закон колеса: {type(e).__name__}: {e}")
        modes = ("chebyshev", "minrange") if deep else ("chebyshev",)
        if not (n.lo <= 0 <= n.hi):
            for mode in modes:
                try:
                    y = self.reciprocal(mode=mode)
                    common = max(self.g, y.g)                 # тоньше — точно, включение цело
                    prod = (self.regrid(common) * y.regrid(common)).notch()
                    one = 1 << (2 * common)
                    if not (prod.lo <= one <= prod.hi):
                        out.append(f"обратный элемент ({mode}): x·(1/x) = ⟦{prod.lo}; {prod.hi}⟧·2^−{prod.g} не несёт 1")
                except Exception as e:                          # noqa: BLE001
                    out.append(f"обратный элемент ({mode}): {type(e).__name__}: {e}")
        if n.lo > 0:
            for mode in modes:
                try:
                    s = self.sqrt(mode=mode)
                    sq = (s * s).notch()
                    shift = sq.g - g
                    lo_x, hi_x = n.lo << shift, n.hi << shift
                    if sq.hi < lo_x or sq.lo > hi_x:
                        out.append(f"корень ({mode}): (√x)² = ⟦{sq.lo}; {sq.hi}⟧·2^−{sq.g} не пересекает засечку x")
                except Exception as e:                          # noqa: BLE001
                    out.append(f"корень ({mode}): {type(e).__name__}: {e}")
        return out

    def record(self) -> str:
        terms = " ".join(f"{'+' if v >= 0 else '-'} {abs(v)}·ε{k}" for k, v in sorted(self._coef.items()))
        return f"({self.c0} {terms})·2^−{self.g}"

    def __repr__(self) -> str:
        return f"Affine{self.record()}"


# ── векторная часть: зонотопы для сетей ──────────────────────────────

def zono_bytes(n: int, D: int, itemsize: int = 8) -> int:
    """Пик массива коэффициентов входного зонотопа в байтах — тождеством.

    Массив есть (n, D, D) целых по `itemsize` байт, потому пик равен
    n·D²·itemsize ТОЧНО, целыми, без оценки и без округления. Дорога
    независимая от него — `C.nbytes` самого `zono_input`; расхождение двух
    есть отказ, а не повод округлить."""
    if min(n, D, itemsize) < 1:
        raise ValueError("n, D и itemsize — целые ≥ 1")
    return n * D * D * itemsize


def zono_input(x_units: np.ndarray, eps_units: int) -> tuple:
    """Входной ящик как зонотоп: по символу на координату. c0: (n, D); C: (n, D, D) диагональ."""
    n, D = x_units.shape
    C = np.zeros((n, D, D), np.int64)
    idx = np.arange(D)
    C[:, idx, idx] = eps_units
    return x_units.astype(np.int64).copy(), C


def zono_affine(c0: np.ndarray, C: np.ndarray, W: np.ndarray, b: np.ndarray, shift: int, bias_shift: int) -> tuple:
    """Слой: центр и каждый коэффициент — точно через ядро; сдвиг сетки наружу
    — в новый символ (радиус = число округлённых членов + 2)."""
    n, S, D = C.shape
    H = W.shape[1]
    v0, t0 = exact_gemm(c0, W)
    v1, t1 = exact_gemm(C.reshape(n * S, D), W)
    if (t0 == INF).any() or (t1 == INF).any():
        raise OverflowError("коэффициенты за int64 — зонотоп на этой ширине не замерен")
    v1 = v1.reshape(n, S, H)
    # округление к нулю одним проходом: для a < 0 это ceil, т. е. (a + 2^s − 1) >> s;
    # отброшенное считаем по остатку (для a ≥ 0 — младшие биты, для a < 0 — то же по модулю)
    m = np.int64((1 << shift) - 1)
    def tz(a):
        return (a + ((a < 0) * m)) >> shift
    c0n = tz(v0) + (b << bias_shift)
    Cn = tz(v1)
    dropped = ((v0 & m) != 0).astype(np.int64) + ((v1 & m) != 0).sum(axis=1) + 2
    Cn = np.concatenate([Cn, dropped[:, None, :]], axis=1)
    return c0n, Cn


def zono_bounds(c0: np.ndarray, C: np.ndarray) -> tuple:
    r = np.abs(C).sum(axis=1)
    return c0 - r, c0 + r


def zono_relu(c0: np.ndarray, C: np.ndarray) -> tuple:
    """Преобразователь DeepZ на целых: y = λx + μ + μ·ε_new, λ = u/(u−l),
    μ = −λl/2; коэффициенты к нулю, отброшенное — в тот же новый символ."""
    l, u = zono_bounds(c0, C)
    n, S, H = C.shape
    dead = u <= 0
    alive = l >= 0
    mixed = ~dead & ~alive
    c0n = c0.copy(); Cn = C.copy()
    c0n[dead] = 0; Cn[np.broadcast_to(dead[:, None, :], Cn.shape)] = 0
    new = np.zeros((n, H), np.int64)
    if mixed.any():
        if float(np.abs(u).max()) * float(np.abs(C).max()) >= 2.0 ** 62 or float(np.abs(u).max()) * float(np.abs(c0).max()) >= 2.0 ** 62:
            raise OverflowError("λ·c за int64 — преобразователь ReLU на этой ширине не считается")
        den = (u - l)                                    # > 0 на смешанных
        den_safe = np.where(mixed, den, 1)
        # λ·c = (u·c) // den к нулю, одним делением: для num < 0 это −((−num) // d) = (num + d − 1) // d;
        # отброшенное считаем по неточности деления (num − q·d ≠ 0)
        def lam(a):
            num = u[:, None, :] * a if a.ndim == 3 else u * a
            d = den_safe[:, None, :] if a.ndim == 3 else den_safe
            q = (num + ((num < 0) * (d - 1))) // d
            return q, (num != q * d)
        # μ = u·(−l) / (2·den): наружу (вверх)
        mu = -(-(u * (-l)) // (2 * den_safe))
        C_mixed, inexact = lam(C)
        c0_mixed, _ = lam(c0)
        c0_mixed = c0_mixed + mu
        dropped = inexact.sum(axis=1) + 2
        Cn = np.where(mixed[:, None, :], C_mixed, Cn)
        c0n = np.where(mixed, c0_mixed, c0n)
        new = np.where(mixed, mu + dropped, 0)
    Cn = np.concatenate([Cn, new[:, None, :]], axis=1)
    return c0n, Cn


def zono_margin(c0: np.ndarray, C: np.ndarray, cls: np.ndarray) -> tuple:
    """Запас z_cls − z_j как аффинная форма (общие символы вычитаются) —
    нижняя граница по всем j ≠ cls. Возвращает (lo, hi) запаса."""
    n, S, K = C.shape
    idx = np.arange(n)
    own0 = c0[idx, cls]; ownC = C[idx, :, cls]
    lo = np.full(n, np.iinfo(np.int64).max); hi = np.full(n, np.iinfo(np.int64).min)
    for j in range(K):
        m = cls != j
        d0 = own0 - c0[:, j]
        dC = ownC - C[:, :, j]
        r = np.abs(dC).sum(axis=1)
        lo = np.where(m, np.minimum(lo, d0 - r), lo)
        hi = np.where(m, np.maximum(hi, d0 + r), hi)
    return lo, hi


# ── самопроверка модуля (ось самоприменения, 09.10) ───────────────────

def self_check(deep: bool = False) -> list:
    """Самопроверка аффинных засечек встречными дорогами; список находок,
    пусто — здорова; исключений не бросает.

    `[вычислено здесь]` Лечение зависимости против Мура: у формы радиуса r
    разность x − x есть точный ноль, а разность засечек по Муру имеет ширину
    2r [внешний источник: Comba–Stolfi 1993]; общие символы — (x + y) − y
    возвращает x точно; `from_notch` при чётной и нечётной ширине включает
    засечку; `exact` без символов; обратный элемент на засечке с нулём —
    отказ поимённо; окно форм (две с символами, точная, отрицательная, из
    засечки нечётной ширины; deep — ещё широкая) проверяет себя
    (`Affine.self_check`). Граница — окно."""
    out = []
    try:
        x = Affine(8, 256, {1: 128})
        y = Affine(8, 40, {2: 8})
        if (x - x).radius != 0 or (x - x).c0 != 0:
            out.append("лечение зависимости: x − x не точный ноль")
        if (x.notch() - x.notch()).hi - (x.notch() - x.notch()).lo != 2 * (x.notch().hi - x.notch().lo):
            out.append("Мур: разность засечек с собой не удвоила ширину (дефект, который лечит форма, исчез из засечки)")
        back = (x + y) - y
        if back.c0 != x.c0 or dict(back.coef) != dict(x.coef):
            out.append(f"общие символы: (x + y) − y = {back.record()} ≠ x")
        for lo, hi in ((100, 104), (100, 105)):
            f = Affine.from_notch(Notch(8, lo, hi))
            nf = f.notch()
            if nf.lo > lo or nf.hi < hi:
                out.append(f"from_notch ⟦{lo}; {hi}⟧: форма не включает засечку")
        if Affine.exact(6, 100).radius != 0 or Affine.exact(6, 100).c0 != 100:
            out.append("exact: форма с символами или чужим центром")
        try:
            Affine(8, 10, {3: 20}).reciprocal()
            out.append("обратный элемент на засечке с нулём не отказал")
        except Unresolved:
            pass
    except Exception as e:                                      # noqa: BLE001
        out.append(f"лечение зависимости: {type(e).__name__}: {e}")
    window = [Affine(8, 256, {1: 128}), Affine(10, 3000, {2: 50, 3: -20}), Affine.exact(6, 100), Affine(8, -300, {4: 10}),
              Affine.from_notch(Notch(8, 100, 105))]
    if deep:
        window.append(Affine(12, 50000, {5: 4000, 6: 3000}))
    for f in window:
        try:
            out += [f"{f.record()}: {m}" for m in f.self_check(deep)]
        except Exception as e:                                  # noqa: BLE001
            out.append(f"{f.record()}: {type(e).__name__}: {e}")
    return out
