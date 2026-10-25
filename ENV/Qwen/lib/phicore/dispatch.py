"""Фасад сведения бинарных действий — тонкая дверь над ФIOS.

Вся маршрутизация — в декларативной таблице `fios.TABLE` с ключом
(контракт, класс)²: ни одного ветвления по классам здесь нет (слово
автора 01.09: «ФIOS чтобы не погрязнуть в if/or»). Фасад держит только
входные ворота: запрет float построением (класс 23 хранит его битово,
не считает), типографские синонимы действий и переходную совместимость
порядка аргументов.
"""

from .rational import Q

from .errors import Unresolved
from . import Фios as fios

#: Типографские синонимы действий.
_ALIAS = {"×": "*", "−": "-", "÷": "/"}

_OPS = ("+", "-", "*", "/")


def _rec(x: object) -> str:
    rec = getattr(x, "record", None)
    if callable(rec):
        return rec()
    if type(x) in (int, Q):
        return f"Φ(C=1, value={Q(x)})"
    return f"Φ(?={x!r})"


def _is_op(x: object) -> bool:
    return isinstance(x, str) and (x in _OPS or x in _ALIAS)


def combine(a: object, b: object, op: str):
    """Свести действие `a <op> b` таблицей ФIOS.

    Несведённое — отказ `Unresolved` с явными Ф-записями; отказы
    носителей и контейнера проходят насквозь как есть. Переходная
    совместимость: порядок combine(a, op, b) тоже принят.
    """
    if _is_op(b) and not _is_op(op):
        b, op = op, b
    op = _ALIAS.get(op, op)
    if op not in _OPS:
        raise Unresolved(op, _rec(a), _rec(b), hint="действие вне таблицы")
    for x in (a, b):
        if isinstance(x, float):
            raise Unresolved(op, _rec(a), _rec(b),
                             hint="float запрещён построением: класс 23 "
                                  "IGNORANCE хранит его битово, не считает")
    return fios.route(a, b, op)


# ---------- самоприменение ----------

def self_check(deep: bool = False) -> list:
    """Самопроверка двери ВСТРЕЧНОЙ дорогой (ось самоприменения).

    `[вычислено здесь]` Дверь сама ничего не считает — она нормализует вход
    и зовёт таблицу; потому сверяется не счёт, а ворота: (1) синонимы —
    каждое типографское действие сводится в то же, что его ASCII-двойник
    (сравнение по точному отношению `to_exact_int_ratio`, не по записи);
    (2) переходный порядок combine(a, op, b) даёт то же, что combine(a, b,
    op); (3) итог на точных входах — против арифметики дробей `Q`: дверь
    над таблицей обязана давать ровно a ∘ b для четырёх действий и окна
    пар целых и дробей; (4) отказы по договору — float на любой стороне и
    действие вне таблицы суть `Unresolved` с названной подсказкой, а не
    иное исключение и не тихий результат; (5) `_rec` даёт Ф-запись для
    int, Q и носителя с `record`, и «?» для чужого. При `deep` окно пар
    расширяется отрицательными и большими знаменателями.

    Граница объявлена: ячейки таблицы и сами носители здесь не судятся —
    у них свои самопроверки (`Фios.self_check`, линейки классов); деление
    на нуль-состояние — предмет таблицы, не двери. Возвращает список
    находок; исключений не бросает.
    """
    out: list = []
    pairs = [(2, 3), (Q(1, 2), Q(1, 3)), (7, Q(2, 5)), (Q(3, 4), 1)]
    if deep:
        pairs += [(-5, 3), (Q(-7, 3), Q(2, 9)), (1000003, 999983), (Q(1, 1 << 40), 3)]

    def ratio(x):
        return tuple(x.to_exact_int_ratio())

    def q_of(x):
        n, d = ratio(x)
        return Q(n, d)

    # (3) итог против арифметики дробей
    for a, b in pairs:
        for op in _OPS:
            try:
                got = q_of(combine(a, b, op))
            except Exception as e:                               # noqa: BLE001
                out.append(f"дверь отказала на {a} {op} {b}: {type(e).__name__}: {e}")
                continue
            qa, qb = Q(a), Q(b)
            want = {"+": qa + qb, "-": qa - qb, "*": qa * qb, "/": qa / qb}[op]
            if got != want:
                out.append(f"итог двери {a} {op} {b} = {got} ≠ {want} по арифметике дробей")
    # (1) синонимы и (2) порядок
    a, b = 6, Q(4, 3)
    canon = {"×": "*", "−": "-", "÷": "/"}            # ожидание записано здесь, не читается из _ALIAS
    if dict(_ALIAS) != canon:
        out.append(f"таблица синонимов {_ALIAS} ≠ ожидаемой {canon}")
    for alias, op in canon.items():
        try:
            if ratio(combine(a, b, alias)) != ratio(combine(a, b, op)):
                out.append(f"синоним «{alias}» сводится не в «{op}»")
        except Exception as e:                                   # noqa: BLE001
            out.append(f"синоним «{alias}» отказал: {type(e).__name__}: {e}")
    for op in _OPS:
        try:
            if ratio(combine(a, op, b)) != ratio(combine(a, b, op)):
                out.append(f"переходный порядок combine(a, «{op}», b) расходится с combine(a, b, «{op}»)")
        except Exception as e:                                   # noqa: BLE001
            out.append(f"переходный порядок «{op}» отказал: {type(e).__name__}: {e}")
    # (4) отказы по договору
    for args, word in (((1.5, 2, "+"), "float"), ((2, 1.5, "*"), "float"), ((1, 2, "^"), "вне таблицы")):
        try:
            combine(*args)
            out.append(f"combine{args} не отказала, хотя обязана")
        except Unresolved as e:
            if word not in str(e):
                out.append(f"отказ combine{args} не называет «{word}»: {e}")
        except Exception as e:                                   # noqa: BLE001
            out.append(f"combine{args} отказала чужим исключением {type(e).__name__}, не Unresolved")
    # (5) запись операнда
    if _rec(3) != "Φ(C=1, value=3)" or "value=1/2" not in _rec(Q(1, 2)):
        out.append(f"_rec для точных даёт {_rec(3)!r}, {_rec(Q(1, 2))!r}")
    if not _rec(object()).startswith("Φ(?="):
        out.append("_rec для чужого объекта не ставит «?»")
    if not _rec(combine(2, 3, "+")).startswith("Φ"):
        out.append("_rec для носителя не берёт его record()")
    return out
