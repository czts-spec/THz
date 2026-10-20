#!/usr/bin/env python3
"""Перевод по запросу: d знаков и объявленный остаток.

Правило среды: вычисление всегда, перевод — по запросу. Перевод есть
обыкновенное действие: запросить `d` знаков, получить `d` знаков и либо
нулевой остаток, либо **объявленное усечение**. `π` здесь не труднее
`1/3`.

Тонкость, ради которой этот файл существует, замерена изнутри. У
библиотеки есть свой перевод — `term.evaluate('mpmath', prec=N)`, — и он
идёт **не от класса, а от быстрого представления**: `PhiNumber` хранит
`float64`, и mpmath получает его же, дополненный до N бит. Запросив у
`π` двести бит, получают 200 бит, из которых верны **12 знаков**, и
усечение не объявляется никак.

Поэтому перевод здесь идёт от класса и параметров (символьный путь), а
быстрый оставлен для сравнения: `backend_gap` меряет расхождение.

    python3 -B ENV/AIAΣA/bin/translate.py pi 40
    python3 -B ENV/AIAΣA/bin/translate.py --gap
"""

from __future__ import annotations

import sys
from pathlib import Path

BIN = Path(__file__).resolve().parent
sys.path.insert(0, str(BIN))

import phi  # noqa: E402

#: Запас знаков, на котором проверяется обрыв развёртки.
GUARD = 12


def _sympy():
    import sympy
    return sympy


def digits(spec: str, d: int = 20) -> dict:
    """`d` знаков величины и честный ответ, точны ли они.

    Возвращает `остаток`: `нулевой` — развёртка оборвалась и перевод
    точен; `усечён` — знаки после d-го есть, и это объявлено.
    """
    sp = _sympy()
    t = phi.build_term(spec)
    exact = sp.nsimplify(t.symbolic)
    wide = sp.N(exact, d + GUARD)
    text = str(wide)
    head, _, tail = text.partition(".")
    body = (head + "." + tail[:d]) if tail else head
    rest = tail[d:d + GUARD]
    terminates = rest.strip("0") == ""
    return {"запись": spec, "знаков": d, "зна��ение": body,
            "остаток": "нулевой" if terminates else "усечён",
            "точно": terminates, "путь": "класс и параметры"}


def backend_gap(spec: str, d: int = 40) -> dict:
    """Сколько знаков даёт быстрый перевод библиотеки и сколько из них верны.

    Разрыв здесь — не ошибка счёта, а следствие пути: быстрый перевод
    идёт через `float64`, и точность класса до него не доходит.
    """
    sp = _sympy()
    t = phi.build_term(spec)
    try:
        import mpmath
        mpmath.mp.dps = d + GUARD
        fast = mpmath.nstr(t.evaluate("mpmath", prec=4 * d), d)
        exact = sp.nsimplify(t.symbolic)
        true = mpmath.nstr(mpmath.mpf(str(sp.N(exact, d + GUARD))), d)
    except Exception as ex:
        return {"запись": spec, "ошибка": type(ex).__name__, "довод": str(ex)}
    agree = 0
    for a, b in zip(fast, true):
        if a != b:
            break
        agree += 1
    return {"запись": spec, "запрошено знаков": d,
            "быстрый путь": fast, "точный путь": true,
            "совпало знаков": agree, "полное совпадение": fast == true,
            "усечение объявлено": False}


GAP_SPECS = ("pi", "e", "root:2:3", "ln:2", "cos:1/3")


def rational_translation_is_exact() -> bool:
    """`7/8` и `1/4` обрываются; `1/3` — нет. Это свойство основания."""
    return (digits("7/8", 12)["точно"] and digits("1/4", 12)["точно"]
            and not digits("1/3", 12)["точно"])


def the_fast_path_loses_digits() -> bool:
    """Быстрый перевод даёт меньше верных знаков, чем запрошено."""
    g = backend_gap("pi", 40)
    return "совпало знаков" in g and g["совпало знаков"] < 20


def truncation_is_declared() -> bool:
    """Здесь усечение объявляется полем, а не молчанием."""
    return digits("pi", 20)["остаток"] == "усечён"


def every_spec_translates() -> bool:
    for s in GAP_SPECS:
        try:
            digits(s, 12)
        except Exception:
            return False
    return True


NOT_MEASURED = (
    "Точность символьного пути обеспечивается SymPy: без него перевод "
    "невозможен вовсе, и дом это объявляет, а не обходит.",
    "«Остаток нулевой» проверяется на запасе в 12 знаков за пределом "
    "запроса. Развёртка, обрывающаяся позже запаса, будет названа "
    "усечённой — ошибка в сторону осторожности, но ошибка.",
    "Разрыв быстрого пути замерен на пяти величинах, а не на всех "
    "достижимых классах.",
)


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if "--gap" in sys.argv:
        for s in GAP_SPECS:
            g = backend_gap(s, 40)
            if "ошибка" in g:
                print(f"{s:<10} ошибка: {g['ошибка']}")
            elif g["полное совпадение"]:
                print(f"{s:<10} запрошено 40, разрыва нет: "
                      f"величина конечна ({g['быстрый путь']})")
            else:
                print(f"{s:<10} запрошено 40, верно {g['совпало знаков']:>3} "
                      f"знаков  быстрый: {g['быстрый путь'][:24]}…")
        raise SystemExit(0)
    spec = args[0] if args else "pi"
    d = int(args[1]) if len(args) > 1 else 20
    r = digits(spec, d)
    print(f"{r['запись']}: {r['значение']}")
    print(f"знаков {r['знаков']}, остаток {r['остаток']}, путь — {r['путь']}")
