"""Регистр-указатель класса — упаковка по правилам, не по битам.

Слово автора 05.09 `[в корпусе]`: «Int8 — указатель на класс числа; 256
значений — регистр должен иметь внутренний формат, чтобы упаковать флаги,
само число и α»; «мы не можем себе позволить битовую роскошь: 8 бит — это
256 значений, значит упаковывать по правилам; 1 бит на значение — норма
XIX века»; «нужны физические регистры: сколько, каких — самые быстрые,
кешируемые»; «α — мантисса классов 1–19, не Ф-число: циклическое
самоопределение недопустимо».

Численная форма. Поля контейнера v6 (класс, четыре бита блокировки, род
нуля, род ω/ε, род невежества) не независимы — их СОВМЕСТНЫЕ состояния
перечислены правилами контейнера: классы 1–19 с блокировкой «нет» либо
«все четыре» (38), нули класса 0 — ∅ absolutum / chaos ordinabile (2),
Center-0 класса 20 — абсолютный либо структурный рода ω / ε (3), указатели
21/22 — заперта пара порядка (2), невежество класса 23 — nan / inf /
finite (3): всего 48 состояний. Основание — код 0 (α = 2, умолчание и
«потерялось») либо слот α-таблицы вычисления 1…4. Байт = состояние·5 +
код α: 48·5 = 240 значений из 256, 16 запасных — против 5 + 8 + 3 = 16 бит
«по биту на значение». Число — четыре ЦЕЛЫХ РЕГИСТРА m1, m2, e1, e2 (одометр). Колесом они не
являются: колесо дома есть `wheels.Wheel(b, m, e)` — одно слагаемое суммы
в сетке (b, w) (`phicore/wheels.py`), а здесь — поля пары (M = m1/m2,
E = e1/e2). Замер 22.09 `[вычислено здесь]`: слово «колесо» стояло в этом
модуле трижды при нуле ввозов `phicore.wheels`; имя снято, потому что оно
обещало устройство, которого тут нет.

Алгебра. Упаковка — смешанная система счисления (состояние, α) → байт,
взаимно однозначная на 240 кодах; состояние вне правил (узор блокировки,
род без класса) — отказ по имени, не ложный код. α-таблица: элемент —
мантисса классов 1–19 над базой 2; Φα(Φα(…)) — отказ (самоопределение).
Геометрия. Два нуля (две микросхемы по 8 бит, у каждой свой) — в окне это
образ 0 с записью ±поглощённого; знак числа — в m1, не в теге.
Физические регистры — proto/registers.py и заметка 05.09. Ни одного float.
"""
from __future__ import annotations

from .container import NumberClass, PhiNumberV6
from .errors import Unresolved
from .rational import Q

ALPHA_RADIX = 5                      # код 0 — основание 2; 1…4 — слоты α-таблицы
ALPHA_SLOTS = ALPHA_RADIX - 1

# --- состояния по правилам контейнера -------------------------------------------

def _states():
    out = []
    for c in range(1, 20):                                   # обычные классы: блокировки нет | все
        out.append((NumberClass(c), "none", False, None, None))
        out.append((NumberClass(c), "all", False, None, None))
    out.append((NumberClass.ZERO, "all", False, None, None))                    # ∅ nihil absolutum
    out.append((NumberClass.ZERO, "all", True, None, None))                     # ∅ chaos ordinabile
    out.append((NumberClass.UNCOMPUTABLE, "all", False, None, None))            # Center-0 absolute
    out.append((NumberClass.UNCOMPUTABLE, "all", True, "infinite", None))       # ω−ω
    out.append((NumberClass.UNCOMPUTABLE, "all", True, "infinitesimal", None))  # ε−ε
    out.append((NumberClass.INFINITESIMAL, "e", False, None, None))
    out.append((NumberClass.INFINITY, "e", False, None, None))
    for kind in ("nan", "inf", "finite_translated"):
        out.append((NumberClass.IGNORANCE, "all", False, None, kind))
    return tuple(out)


#: (класс, узор блокировки, структурный ноль, род ω/ε, род невежества) — 48 состояний.
STATES = _states()
_INDEX = {s: i for i, s in enumerate(STATES)}
CODES = len(STATES) * ALPHA_RADIX    # 240


def _as_int(x, name: str) -> int:
    if isinstance(x, bool) or not isinstance(x, int):
        raise TypeError(f"{name} — только int")
    return x


def state_of(x: PhiNumberV6) -> tuple:
    """Состояние контейнера по правилам; узор вне правил — отказ по имени."""
    if not isinstance(x, PhiNumberV6):
        raise TypeError("state_of — контейнер v6")
    locks = (x.m1.b, x.m2.b, x.e1.b, x.e2.b)
    if locks == (False, False, False, False):
        pattern = "none"
    elif locks == (True, True, True, True):
        pattern = "all"
    elif locks == (False, False, True, True):
        pattern = "e"
    else:
        raise Unresolved("состояние", repr(x), hint=f"узор блокировки {locks} вне правил контейнера (нет | все | пара порядка)")
    kind = x.structural_kind if (x.number_class == NumberClass.UNCOMPUTABLE and x.zero_provenance) else None
    ign = x.ignorance_kind if x.number_class == NumberClass.IGNORANCE else None
    s = (x.number_class, pattern, bool(x.zero_provenance), kind, ign)
    if s not in _INDEX:
        raise Unresolved("состояние", repr(x), hint=f"состояние {s} вне 48 перечисленных правилами")
    return s


# --- α-таблица вычисления ----------------------------------------------------------

def _is_plain_base(alpha) -> bool:
    """α — мантисса классов 1–19 над базой 2: пара либо носитель без своего α."""
    if isinstance(alpha, (bool, float)):
        return False
    if isinstance(alpha, (int, Q)):
        return True
    if type(alpha).__name__ == "PhiAlpha":
        return False                                   # Φα(Φα(…)) — циклическое самоопределение
    if isinstance(alpha, PhiNumberV6):
        return alpha.base_class is None
    return callable(getattr(alpha, "record", None))


class AlphaTable:
    """Регистровый банк оснований вычисления: до четырёх α (коды 1…4); код 0 — 2."""

    __slots__ = ("_slots",)

    def __init__(self):
        self._slots: list = []

    def code_of(self, alpha) -> int:
        if alpha is None or (isinstance(alpha, (int, Q)) and not isinstance(alpha, bool) and Q(alpha) == 2):
            return 0
        if not _is_plain_base(alpha):
            raise Unresolved("α", repr(alpha),
                             hint="основание — мантисса классов 1–19 над базой 2, не Ф-число со своим α: "
                                  "циклическое самоопределение недопустимо (слово автора 05.09)")
        for i, a in enumerate(self._slots):
            if a == alpha:
                return i + 1
        if len(self._slots) >= ALPHA_SLOTS:
            raise Unresolved("α-таблица", repr(alpha),
                             hint=f"регистровый банк оснований полон: {ALPHA_SLOTS} слота на байт (48·5 = 240 ≤ 256); "
                                  f"пятое основание — второй байт, не молчаливая подмена")
        self._slots.append(alpha)
        return len(self._slots)

    def alpha(self, code: int):
        code = _as_int(code, "code")
        if code == 0:
            return Q(2)
        if not 1 <= code <= len(self._slots):
            raise ValueError(f"код α {code}: слот пуст")
        return self._slots[code - 1]

    def __len__(self) -> int:
        return len(self._slots)


# --- байт ---------------------------------------------------------------------------

def pack(state: tuple, alpha_code: int = 0) -> int:
    """Байт = индекс состояния · 5 + код α (0…239).

    Это АДРЕС в таблице — два поля фиксированной ширины (48 состояний × 5
    кодов α), а не позиционная запись величины: одноимённая дверь
    `wheels.WheelByte.pack` кладёт в байт КОЛЕСО (позиция и цифра сетки
    ширины 1), `ccore.pack` переводит партию пар в структуры C-ядра. Здесь
    величины нет вовсе — только состояние регистров и код основания.
    """
    if state not in _INDEX:
        raise ValueError(f"состояние {state} вне перечисленных")
    code = _as_int(alpha_code, "alpha_code")
    if not 0 <= code < ALPHA_RADIX:
        raise ValueError(f"код α {code} вне 0…{ALPHA_SLOTS}")
    return _INDEX[state] * ALPHA_RADIX + code


def unpack(byte: int) -> tuple:
    """(состояние, код α); байт ≥ 240 — запасной код, отказ по имени.

    Деление на `ALPHA_RADIX` здесь — разбор адреса таблицы на два поля, не
    разложение записи по основанию: колеса (`wheels`) в байте нет, и
    `wheels.WheelByte.unpack` (байт → колесо) с этой дверью не совпадает;
    `ccore.unpack` разбирает структуры C, а не байт.
    """
    byte = _as_int(byte, "byte")
    if not 0 <= byte < 256:
        raise ValueError("байт — 0…255")
    if byte >= CODES:
        raise Unresolved("байт", str(byte), hint=f"запасной код {byte} ≥ {CODES}: состояния за ним нет")
    return STATES[byte // ALPHA_RADIX], byte % ALPHA_RADIX


class _BaseKey:
    """Ключ основания контейнера (класс базы + параметры) — элемент α-таблицы."""
    __slots__ = ("key",)

    def __init__(self, key):
        self.key = key

    def __eq__(self, other):
        return isinstance(other, _BaseKey) and self.key == other.key

    def __hash__(self):
        return hash(self.key)

    def record(self) -> str:
        return f"α[{self.key[0].name}:{dict(self.key[1])}]"

    __repr__ = record


def self_check(deep: bool = False) -> list:
    """Самопроверка регистра-указателя ВСТРЕЧНОЙ дорогой (ось самоприменения).

    `[вычислено здесь]` Байт собирается умножением индекса состояния на
    основание α-кодов и сложением с кодом; встречная дорога — РАЗБОР
    обратно (`unpack`), то есть деление с остатком: для всех 240 живых
    байтов пара (состояние, код) обязана вернуться той же, а число кодов —
    сойтись с произведением 48 · 5, посчитанным здесь заново по самому
    перечню состояний. Вторая дорога внешняя: состояние берётся у
    СОСЕДНЕГО модуля — у настоящих контейнеров (`state_of`), и байт,
    собранный из него, обязан разобраться в то же состояние. Третья:
    запасные коды 240…255 обязаны отказывать поимённо (`Unresolved`), а не
    молча отдавать чужое состояние.

    Граница: перечень состояний закрыт правилами контейнера (48 совместных
    состояний, узор блокировки — нет | все | пара порядка), и самопроверка
    говорит о ПОЛНОТЕ разбора, не о полноте самого перечня: контейнер с
    узором вне правил сюда не входит — он отбит отказом в `state_of`, и это
    объявлено там. Вторая граница (12.09, найдена приёмкой первых принципов):
    живыми контейнерами класса 23 глубокая дорога НЕ пробует — проба
    потребовала бы построить float, а этот модуль мостом к железу не
    является и в список мостов (`bench_first_principles`: affine, container,
    exact_gemm, hyper, ieee) не входит; три состояния класса 23 покрыты
    основной дорогой — перебором всех 240 байтов.
    Возвращает список находок; исключений не бросает.
    """
    out = []
    if len(set(STATES)) != len(STATES):
        out.append(f"состояния не уникальны: {len(STATES)} записей, "
                   f"{len(set(STATES))} различных")
    if CODES != len(STATES) * ALPHA_RADIX:
        out.append(f"число кодов {CODES} ≠ {len(STATES)}·{ALPHA_RADIX}")
    for i, st in enumerate(STATES):
        for code in range(ALPHA_RADIX):
            try:
                byte = pack(st, code)
            except Exception as e:                              # noqa: BLE001
                out.append(f"сборка байта отказала на {st}/{code}: "
                           f"{type(e).__name__}: {e}")
                continue
            if byte != i * ALPHA_RADIX + code:
                out.append(f"байт {byte} не отвечает индексу {i} и коду {code}")
            try:
                back_state, back_code = unpack(byte)
            except Exception as e:                              # noqa: BLE001
                out.append(f"разбор байта {byte} отказал: {type(e).__name__}: {e}")
                continue
            if (back_state, back_code) != (st, code):
                out.append(f"разбор байта {byte} дал {back_state}/{back_code}, "
                           f"собран из {st}/{code}")
    for byte in range(CODES, 256):
        try:
            unpack(byte)
            out.append(f"запасной код {byte} разобран вместо отказа")
        except Unresolved:
            pass
        except Exception as e:                                  # noqa: BLE001
            out.append(f"запасной код {byte} отказал не Ф-отказом: "
                       f"{type(e).__name__}")
    if not deep:
        return out
    probes = (PhiNumberV6.from_int(5), PhiNumberV6.from_ratio(3, 2),
              PhiNumberV6.void(), PhiNumberV6.chaos(), PhiNumberV6.center_zero(),
              PhiNumberV6.center_zero(True, "infinite"),
              PhiNumberV6.center_zero(True, "infinitesimal"),
              PhiNumberV6.omega(), PhiNumberV6.epsilon())
    table = AlphaTable()
    for x in probes:
        try:
            st = state_of(x)
            byte, регистры = registers_of(x, table)
        except Exception as e:                                  # noqa: BLE001
            out.append(f"регистры отказали на {x.record()}: {type(e).__name__}: {e}")
            continue
        if unpack(byte)[0] != st:
            out.append(f"байт живого контейнера разобрался в чужое состояние: "
                       f"{x.record()}")
        if регистры != (x.m1.v, x.m2.v, x.e1.v, x.e2.v):
            out.append(f"регистры не те: {регистры} против записи "
                       f"{x.record()}")
    return out


def registers_of(x: PhiNumberV6, table: AlphaTable) -> tuple:
    """(байт, (m1, m2, e1, e2)) — регистры числа: один байт состояния и
    основания, четыре целых регистра пары (не колёса — см. шапку модуля)."""
    s = state_of(x)
    if x.base_class is None:
        code = 0
    else:
        key = (x.base_class, tuple(sorted((k, repr(v)) for k, v in x.base_params.items() if k != "absorbed")))
        code = table.code_of(_BaseKey(key))
    return pack(s, code), (x.m1.v, x.m2.v, x.e1.v, x.e2.v)
