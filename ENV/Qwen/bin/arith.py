#!/usr/bin/env python3
"""Замер изнутри: что арифметика среды делает на самом деле.

Промпт называет пустоту среды словами — «арифметику имеют 8 классов из
19». Изнутри это проверяется прогоном, и прогон говорит больше, чем
слова: важно не только **сколько** классов работает, но и **что
происходит там, где правила нет**.

Три разных исхода различаются здесь и не сваливаются в один:

* `точно`   — действие выполнено, класс результата верен;
* `ложный класс` — действие выполнено, класс результата **неверен**;
* `падение` — исключение.

Худший из трёх — средний, и именно он не назван в промпте. Отказ виден,
падение видно, ложный класс выглядит как ответ.

    python3 -B ENV/AIAΣA/bin/arith.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

BIN = Path(__file__).resolve().parent
sys.path.insert(0, str(BIN))

import phi  # noqa: E402

warnings.filterwarnings("ignore")

#: Представитель каждого достижимого класса: краткая запись фасада.
REPS = {1: "7", 2: "root:2:3", 4: "unity:8:1", 5: "e", 6: "ln:2",
        10: "cos:1/3", 11: "pi", 12: "cosh:1"}

#: Классы 1—19: вычислимые. Знаменатель всех долей ниже.
COMPUTABLE = tuple(range(1, 20))

#: Класс, в который библиотека складывает всё, чему правила не нашлось.
FALLBACK_CLASS = 3


def reachable() -> set[int]:
    """Классы, которые вообще можно построить. Замер, а не список."""
    return {c for c in phi.classes_reached() if c in COMPUTABLE}


def unreachable() -> list[int]:
    return sorted(set(COMPUTABLE) - reachable())


def _f(term) -> float | None:
    try:
        return float(term.evaluate("float"))
    except Exception:
        return None


def op_report() -> dict[int, dict[str, str]]:
    """Для каждого достижимого класса: что даёт x+x, x*x, x/x, x-x."""
    out: dict[int, dict[str, str]] = {}
    for c, spec in sorted(REPS.items()):
        row: dict[str, str] = {}
        for name, fn in (("+", lambda a: a + a), ("*", lambda a: a * a),
                         ("/", lambda a: a / a), ("-", lambda a: a - a)):
            try:
                t = phi.build_term(spec)
                r = fn(t)
                row[name] = phi.lib().NumberClass(int(r.number_class)).name
            except Exception as ex:
                row[name] = f"падение:{type(ex).__name__}"
        out[c] = row
    return out


def crashes() -> list[tuple[int, str]]:
    """Действия, кончающиеся исключением. Падение — не отказ."""
    return sorted((c, op) for c, row in op_report().items()
                  for op, v in row.items() if v.startswith("падение"))


#: Тождества, которые обязаны держаться при любой записи числа.
#: Проверяются численно: `x − x = 0` и `x / x = 1`.
def identities() -> dict[str, list[int]]:
    """Где `x − x` не ноль и где `x / x` не единица. Числами, не словами."""
    bad_sub: list[int] = []
    bad_div: list[int] = []
    for c, spec in sorted(REPS.items()):
        try:
            v = _f(phi.build_term(spec) - phi.build_term(spec))
            if v is None or abs(v) > 1e-15:
                bad_sub.append(c)
        except Exception:
            bad_sub.append(c)
        try:
            v = _f(phi.build_term(spec) / phi.build_term(spec))
            if v is None or abs(v - 1.0) > 1e-15:
                bad_div.append(c)
        except Exception:
            bad_div.append(c)
    return {"x−x≠0": bad_sub, "x/x≠1": bad_div}


def negation_is_lost() -> bool:
    """`−x` возвращает `x` там, где знак некуда положить.

    Отсюда `√3 − √3 = 2√3`: вычитание есть сложение с отрицанием, а
    отрицание у класса 2 не доходит ни до записи, ни до параметров.
    """
    t = phi.build_term("root:2:3")
    a, b = _f(t), _f(-t)
    return a is not None and b is not None and a == b and a > 0


def false_class_cases() -> list[dict]:
    """Случаи, где класс результата назван и неверен.

    Проверяется не мнением: `π + e` не есть корень многочлена с целыми
    коэффициентами (алгебраическая независимость π и e не доказана, но
    и алгебраичность суммы не доказана — а класс 3 утверждает её). Здесь
    ловится более грубое: символьный терм результата сохранён верно и
    противоречит выставленному классу.
    """
    out = []
    for left, right, why in (
        ("pi", "e", "сумма π и e объявлена корнем многочлена"),
        ("pi", "ln:2", "сумма π и ln 2 объявлена корнем многочлена"),
    ):
        try:
            r = phi.build_term(left) + phi.build_term(right)
        except Exception:
            continue
        cls = int(r.number_class)
        sym = str((r.class_params or {}).get("sympy_expr", ""))
        if cls == FALLBACK_CLASS and sym:
            out.append({"действие": f"{left} + {right}",
                        "класс": cls, "терм": sym, "почему": why})
    return out


def fallback_is_a_class_not_a_refusal() -> bool:
    """Главная находка: несработавшее правило даёт класс, а не отказ."""
    return len(false_class_cases()) > 0


def census() -> dict:
    ident = identities()
    return {
        "классов вычислимых": len(COMPUTABLE),
        "достижимо": len(reachable()),
        "недостижимо": len(unreachable()),
        "падений": len(crashes()),
        "x−x≠0": len(ident["x−x≠0"]),
        "x/x≠1": len(ident["x/x≠1"]),
        "ложных классов": len(false_class_cases()),
    }


NOT_MEASURED = (
    "Достижимость меряется через фасад `phi.py`: класс без фабрики в "
    "фасаде считается недостижимым. Если библиотека умеет строить его "
    "иначе, замер занизит охват — и это ноль инструмента, а не ноль "
    "предмета.",
    "Тождества проверяются численно с порогом 1e-15. Порог есть "
    "соглашение: он ловит потерю знака и не ловит расхождения ниже "
    "порога.",
    "«Ложный класс» ловится только там, где сохранён символьный терм. "
    "Случаи без терма не проверяются вовсе и в счёт не идут.",
    "Перебор идёт по одному представителю на класс. Один представитель "
    "не есть класс: другой параметр может дать другой исход.",
)


if __name__ == "__main__":
    print("достижимо классов:", sorted(reachable()))
    print("недостижимо:", unreachable())
    print()
    for c, row in op_report().items():
        name = phi.lib().NumberClass(c).name
        print(f"{c:>2} {name:<16} " +
              "  ".join(f"{k}→{v}" for k, v in row.items()))
    print()
    for k, v in census().items():
        print(f"{k}: {v}")
    print()
    for case in false_class_cases():
        print(f"ЛОЖНЫЙ КЛАСС: {case['действие']} → класс {case['класс']}, "
              f"терм сохранён верно: {case['терм']}")
    print("отрицание теряется:", negation_is_lost())
