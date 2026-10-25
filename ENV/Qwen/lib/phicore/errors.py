"""Отказы ядра. Отказ — ответ, а не сбой (иммунитет AIAΣA v0.0.1)."""


class Unresolved(Exception):
    """Действие не сведено ни одним правилом.

    Несёт явные Ф-записи операндов: отказ обязан называть, ЧТО не
    сведено. Никогда не подменяется правдоподобным классом.
    """

    def __init__(self, op: str, left: str, right: str = "", hint: str = ""):
        self.op, self.left, self.right, self.hint = op, left, right, hint
        msg = f"не сведено: {left} {op} {right}".rstrip()
        if hint:
            msg += f"; {hint}"
        super().__init__(msg)


class ExactRepresentationNotAvailable(ArithmeticError):
    """Точного результата в запрошенной записи нет; смена точности —
    трансляция по явному запросу, не скрытое вычисление (решение автора
    09.07.2026, phi_number_v6). Единственный канонический класс пакета:
    контейнер v6 поднимает ЕГО ЖЕ (шов двух одноимённых классов найден
    мостом этапа E и закрыт 02.09)."""


# ---------- самоприменение ----------

def self_check(deep: bool = False) -> list:
    """Самопроверка отказов ядра ВСТРЕЧНОЙ дорогой (ось самоприменения).

    `[вычислено здесь]` Отказ — ответ, а не сбой, и у ответа есть форма:
    (1) запись отказа собирается из слотов — «не сведено: left op right;
    hint» — сверяется с прямой сборкой строки по тем же слотам, а не чтением
    самой записи; правый операнд пуст — хвост без лишнего пробела, подсказки
    нет — хвост без «;»; (2) слоты несутся как есть (op, left, right, hint) —
    отказ обязан называть, ЧТО не сведено; (3) два рода отказа различимы
    построением: `Unresolved` — не арифметический (не ArithmeticError),
    `ExactRepresentationNotAvailable` — арифметический и не `Unresolved`,
    так что ловушка одного рода не глотает другой; (4) отказ не подменяется
    правдоподобным классом: `Unresolved` не наследует ValueError/TypeError
    интерпретатора. При `deep` форма сверяется на окне из восьми наборов
    слотов с пустыми и непустыми полями. Возвращает список находок;
    исключений не бросает.
    """
    out: list = []
    window = [("+", "Φ(C=1, value=1)", "Φ(C=10 TRIGONOMETRIC, q=1, e=1)", "π независим: правила нет"),
              ("/", "Φ(C=1, value=1)", "Φ(C=0 ZERO)", ""),
              ("golgi", "Φ(?=1.5)", "", "float запрещён построением")]
    if deep:
        window += [("-", "a", "b", "h"), ("*", "a", "", ""), ("**", "", "b", "h; h2"), ("sqrt", "x", "", "граница"), ("∘", "f", "g", "")]
    for op, left, right, hint in window:
        try:
            e = Unresolved(op, left, right, hint=hint)
        except Exception as ex:                                   # noqa: BLE001
            out.append(f"построение отказа {op!r} упало: {type(ex).__name__}: {ex}")
            continue
        want = ("не сведено: " + left + " " + op + " " + right).rstrip() + (("; " + hint) if hint else "")
        if str(e) != want:
            out.append(f"запись отказа {str(e)!r} ≠ прямой сборке из слотов {want!r}")
        if (e.op, e.left, e.right, e.hint) != (op, left, right, hint):
            out.append(f"слоты отказа не несутся как есть: {(e.op, e.left, e.right, e.hint)}")
        if not hint and str(e).endswith(";"):
            out.append("пустая подсказка оставила висячую «;»")
    if issubclass(Unresolved, ArithmeticError):
        out.append("Unresolved — арифметический отказ: роды отказов слились")
    if not issubclass(ExactRepresentationNotAvailable, ArithmeticError) or issubclass(ExactRepresentationNotAvailable, Unresolved):
        out.append("ExactRepresentationNotAvailable — не арифметический отказ либо подменён Unresolved")
    for bad in (ValueError, TypeError, KeyError):
        if issubclass(Unresolved, bad):
            out.append(f"Unresolved наследует {bad.__name__} интерпретатора — отказ смешан со сбоем")
    try:
        raise Unresolved("+", "a", "b", hint="h")
    except ArithmeticError:
        out.append("Unresolved пойман ловушкой арифметического отказа")
    except Unresolved:
        pass
    return out
