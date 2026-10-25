"""Аппарат Гольджи — явная выходная дверь дома со встроенной мини-ФIOS.

Слово автора 09.10: «Фios. от ДНК мы взяли саму структуру … однако аппарат
Гольджи, по аналогии, ни разу не задействован?» и «сделай явный Гольджи
аппарат со встроенным миниФIOS. это полезно?». Замер 09.10 (записка
«аппарат Гольджи в доме») `[вычислено здесь]`: выходных дверей у дома
двенадцать (record, IEEE, профиль кремния, C, Rust, BLAS, tag.pack, …),
контроль качества на пути стоит только у `profile.enforce`, на главном
пути route → record вызовов `self_check` ноль. То есть величина уходила из
дома без досмотра и без одного адреса.

Аналогия, взятая у клетки [внешний источник: Alberts et al., Molecular
Biology of the Cell, гл. 13 «Intracellular Vesicular Traffic», NCBI
Bookshelf NBK26941]: аппарат Гольджи принимает продукт, достраивает его,
СОРТИРУЕТ по назначению и упаковывает в везикулу с адресом. Здесь то же —
четырьмя ходами, ни один не новый для дома:

1. **контроль качества** — `self_check()` носителя обязан молчать; находка
   — отказ поимённо, величина наружу не идёт;
2. **модификация** — приведение к форме адресата (засечка на сетке,
   точная дробь, образ IEEE с названной потерей);
3. **сортировка** — мини-ФIOS: декларативная таблица (класс, адресат) →
   упаковщик, как у `Фios.TABLE` — ни одного ветвления по классам в
   коде; пустая ячейка — иммунный отказ с подсказкой;
4. **везикула** — `Vesicle`: значение в форме адресата, адресат, класс,
   метка провенанса и счёт досмотра.

Что не делается здесь: сведение действий (это ФIOS), хранение (это
контейнер), перевод между классами (это линейки). Гольджи — только выход.
"""
from __future__ import annotations

from dataclasses import dataclass

from . import Фios as fios
from . import ieee as _ieee
from .errors import Unresolved
from .gridcore import Notch
from .rational import Q

#: адресаты выхода
RECORD, FRACTION, NOTCH, IEEE = "record", "fraction", "notch", "ieee"
DESTINATIONS = (RECORD, FRACTION, NOTCH, IEEE)

#: сетка засечки по умолчанию (разрядов после запятой в двоичной записи)
DEFAULT_GRID = 40

PROVENANCE = "[вычислено здесь]"


@dataclass(frozen=True)
class Vesicle:
    """Везикула: упакованное значение с адресом и меткой."""
    value: object
    destination: str
    key: tuple
    packer: str
    checks: int            # сколько досмотров пройдено (0 — носитель без self_check)
    provenance: str = PROVENANCE
    note: str = ""         # цена/потеря, если адресат её несёт (IEEE)


# ---------- упаковщики: модификация под адресата ----------

def _rec(x: object) -> str:
    rec = getattr(x, "record", None)
    if callable(rec):
        return rec()
    if type(x) in (int, Q):
        return f"Φ(C=1, value={Q(x)})"
    return f"Φ(?={x!r})"


def pack_record(x, g: int):
    return _rec(x), ""


def pack_fraction(x, g: int):
    if type(x) in (int, Q):
        return Q(x), ""
    r = getattr(x, "to_exact_int_ratio", None)
    if callable(r):
        n, d = r()
        return Q(n, d), ""
    f = getattr(x, "to_fraction", None)
    if callable(f):
        return f(), ""
    raise Unresolved("golgi", _rec(x), hint="носитель не отдаёт точную дробь")


def pack_notch(x, g: int):
    n = Notch.of(x, g)
    return n, f"ширина {n.width_units} узлов сетки 2^−{g}"


def pack_ieee_exact(x, g: int):
    q = pack_fraction(x, g)[0]
    f = _ieee.round_to(q)
    return f, f"потеря {q - _ieee.dyadic(f)}"


def pack_ieee_from_notch(x, g: int):
    """Образ IEEE по вилке: вилка обязана не накрывать границу округления;
    если накрывает — сетка удваивается, до трёх раз (цена названа в
    разрядах), дальше отказ `image_from_enclosure` идёт как есть."""
    last = None
    for gg in (g, 2 * g, 4 * g, 8 * g):
        n = Notch.of(x, gg)
        try:
            f = _ieee.image_from_enclosure(n.to_interval())
        except Unresolved as e:
            last = e
            continue
        return f, f"образ по вилке шириной {n.width_units} узлов сетки 2^−{gg}"
    raise last


PACKERS = {
    "pack_record": pack_record,
    "pack_fraction": pack_fraction,
    "pack_notch": pack_notch,
    "pack_ieee_exact": pack_ieee_exact,
    "pack_ieee_from_notch": pack_ieee_from_notch,
}

# ---------- мини-ФIOS: строки → таблица (класс, адресат) → упаковщик ----------

#: строки таблицы: область классов — «live» (все живые), «gridded» (с
#: вещественным сертификатом сетки) либо точный ключ; порядок — точные
#: строки впереди общих, как в ФIOS
ROWS = (
    ("live", RECORD, "pack_record"),
    (fios.E1, FRACTION, "pack_fraction"),
    (fios.E21, FRACTION, "pack_fraction"),
    (fios.E22, FRACTION, "pack_fraction"),
    (fios.E1, IEEE, "pack_ieee_exact"),
    ("gridded", NOTCH, "pack_notch"),
    ("gridded", IEEE, "pack_ieee_from_notch"),
)


def _build(rows) -> dict:
    table: dict = {}
    for scope, dest, packer in rows:
        if packer not in PACKERS:
            raise ValueError(f"мини-ФIOS: упаковщик {packer!r} не существует")
        keys = (fios.LIVE if scope == "live" else fios.GRIDDED if scope == "gridded" else (scope,))
        for k in keys:
            table.setdefault((k, dest), packer)
    return table


TABLE = _build(ROWS)


def coverage() -> dict:
    """Заполненность мини-ФIOS: |LIVE|·|DESTINATIONS| — знаменатель, пустое
    объявляется числом (ноль инструмента не есть ноль предмета)."""
    denom = len(fios.LIVE) * len(DESTINATIONS)
    per = {d: sum(1 for (k, dd) in TABLE if dd == d) for d in DESTINATIONS}
    return {"знаменатель": denom, "заполнено": len(TABLE), "пусто": denom - len(TABLE), "по_адресатам": per}


# ---------- дверь ----------

def inspect(x: object) -> int:
    """Контроль качества: самопроверка носителя обязана молчать.

    Возвращает число пройденных досмотров (1 — носитель с self_check, 0 —
    без); находка — отказ поимённо, величина наружу не идёт."""
    sc = getattr(x, "self_check", None)
    if not callable(sc):
        return 0
    try:
        findings = sc()
    except Exception as e:                                       # noqa: BLE001
        # самопроверка носителя упала исключением — это тоже находка о
        # носителе (порченый слот), и наружу она идёт отказом поимённо,
        # а не исключением интерпретатора
        raise Unresolved("golgi", _rec(x), hint=f"контроль качества: самопроверка упала: {type(e).__name__}: {e}")
    if findings:
        raise Unresolved("golgi", _rec(x), hint=f"контроль качества: {findings[0]}")
    return 1


def ship(x: object, destination: str, g: int = DEFAULT_GRID) -> Vesicle:
    """Выпустить величину из дома везикулой по адресу.

    Порядок: контроль качества → сортировка по мини-ФIOS → модификация
    упаковщиком → везикула. Несведённое — отказ `Unresolved` с подсказкой;
    float на входе — отказ построением (класс 23 хранит его битово)."""
    if isinstance(x, float):
        raise Unresolved("golgi", _rec(x), hint="float запрещён построением: класс 23 IGNORANCE хранит его битово, не считает")
    if destination not in DESTINATIONS:
        raise Unresolved("golgi", _rec(x), hint=f"адресат {destination!r} вне {DESTINATIONS}")
    key = fios.classify(x)
    if key is None:
        raise Unresolved("golgi", _rec(x), hint="операнд вне онтологии 24 классов")
    checks = inspect(x)
    packer = TABLE.get((key, destination))
    if packer is None:
        raise Unresolved("golgi", _rec(x), hint=f"ячейка мини-ФIOS ({key[0]},{key[1]}) → «{destination}» пуста")
    value, note = PACKERS[packer](x, g)
    return Vesicle(value, destination, key, packer, checks, PROVENANCE, note)


# ---------- самоприменение ----------

def self_check(deep: bool = False) -> list:
    """Самопроверка Гольджи ВСТРЕЧНЫМИ дорогами (ось самоприменения).

    `[вычислено здесь]` (1) таблица как данные: каждая ячейка называет
    существующий упаковщик, у каждого живого класса есть адресат record,
    у каждого класса с сертификатом сетки — notch и ieee; (2) дверь против
    прямого счёта: дробь возвращается той же дробью, засечка содержит
    известную точку (Q — точно; π через PiMonomial — вилка содержит
    отрезок Архимеда 223/71 < π < 22/7 при g ≥ 8), образ IEEE точной
    двоичной дроби несёт нулевую потерю, недвоичной — ненулевую;
    (3) отказы по договору — float, чужой адресат, носитель с находкой
    самопроверки (подложный носитель) — суть `Unresolved` с подсказкой,
    а не тихая везикула. При `deep` засечка берётся на g = 80.
    Возвращает список находок; исключений не бросает."""
    out: list = []
    for (k, d), p in TABLE.items():
        if p not in PACKERS:
            out.append(f"ячейка {k} → {d} называет несуществующий упаковщик {p!r}")
    for k in fios.LIVE:
        if (k, RECORD) not in TABLE:
            out.append(f"у класса {k} нет адресата record")
    for k in fios.GRIDDED:
        for d in (NOTCH, IEEE):
            if (k, d) not in TABLE:
                out.append(f"у класса с сертификатом {k} нет адресата {d}")
    g = 80 if deep else DEFAULT_GRID
    try:
        v = ship(Q(3, 2), FRACTION)
        if v.value != Q(3, 2) or v.checks != 1:
            out.append(f"дробь вернулась не той же: {v}")
        if ship(Q(3, 2), IEEE).note != "потеря 0":
            out.append("образ IEEE двоичной дроби 3/2 несёт потерю")
        if "потеря 0" == ship(Q(1, 3), IEEE).note:
            out.append("образ IEEE дроби 1/3 без названной потери")
        n = ship(Q(1, 3), NOTCH, g).value
        if not n.contains(Q(1, 3)):
            out.append("засечка 1/3 не содержит 1/3")
        from .pi_value import PiMonomial
        n = ship(PiMonomial(1, 1), NOTCH, g).value
        iv = n.to_interval()
        if not (iv.lo <= Q(223, 71) and Q(22, 7) <= iv.hi or (iv.lo >= Q(223, 71) and iv.hi <= Q(22, 7))):
            out.append(f"засечка π на g={g} не согласна с отрезком Архимеда: {iv.lo}…{iv.hi}")
        if ship(PiMonomial(1, 1), RECORD).value != PiMonomial(1, 1).record():
            out.append("запись π через Гольджи не равна record() носителя")
    except Exception as e:                                       # noqa: BLE001
        out.append(f"дверь отказала на здоровом входе: {type(e).__name__}: {e}")
    # подложный носитель — не объявление класса, а сборка типом: приёмка
    # согласованности считает классы с методом self_check по ast, и
    # вспомогательный класс внутри самопроверки был бы сосчитан как носитель
    _Bad = type("_Bad", (), {"self_check": lambda self, deep=False: ["подложный носитель: слот порчен"],
                             "record": lambda self: "Φ(подлог)"})
    for args, word in (((1.5, RECORD), "float"), ((Q(1), "x64"), "адресат"), ((_Bad(), RECORD), "онтологии")):
        try:
            ship(*args)
            out.append(f"ship{args} не отказала, хотя обязана")
        except Unresolved as e:
            if word not in str(e):
                out.append(f"отказ ship{args} не называет «{word}»: {e}")
        except Exception as e:                                   # noqa: BLE001
            out.append(f"ship{args} отказала чужим исключением {type(e).__name__}")
    try:
        inspect(_Bad())
        out.append("контроль качества пропустил носитель с находкой")
    except Unresolved as e:
        if "контроль качества" not in str(e):
            out.append(f"отказ контроля качества без имени: {e}")
    return out
