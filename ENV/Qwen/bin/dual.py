#!/usr/bin/env python3
"""Двойственность: величина строится двумя путями, и пути сравниваются.

Правило среды: совпадение — здоровье, расхождение — сигнал, не требующий
внешней задачи. Один путь даёт правдоподобное число и не даёт повода ему
не верить; два пути дают повод.

Здесь пары исполнимые — те, для которых в собранной библиотеке есть оба
пути. Пар, требующих недостижимых классов (ζ, Γ, эллиптические), тут
нет, и это не умолчание: их считает `arith.py` как недостижимые.

    python3 -B ENV/AIAΣA/bin/dual.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

BIN = Path(__file__).resolve().parent
sys.path.insert(0, str(BIN))

import phi  # noqa: E402

warnings.filterwarnings("ignore")

#: Допуск сравнения. Соглашение, а не свойство величин: он объявлен
#: числом именно затем, чтобы не выдавать себя за точность.
TOL = 1e-12


def _v(t) -> float:
    return float(t.evaluate("float"))


def _pairs():
    """Пары путей. Каждая — две независимые дороги к одной величине."""
    b = phi.build_term
    return (
        ("cos(π/3)", "класс 10", lambda: _v(b("cos:1/3")),
         "рациональное, класс 1", lambda: _v(b("1/2"))),
        ("sin(π/6)", "класс 10", lambda: _v(b("sin:1/6")),
         "рациональное, класс 1", lambda: _v(b("1/2"))),
        ("√3·√3", "класс 2, умножение", lambda: _v(b("root:2:3") * b("root:2:3")),
         "целое, класс 1", lambda: _v(b("3"))),
        ("log₂8", "класс 6", lambda: _v(b("log:2:8")),
         "целое, класс 1", lambda: _v(b("3"))),
        ("sinh 1", "класс 12", lambda: _v(b("sinh:1")),
         "(e − 1/e)/2 через класс 5",
         lambda: (_v(b("e")) - 1.0 / _v(b("e"))) / 2),
        ("cosh 1", "класс 12", lambda: _v(b("cosh:1")),
         "(e + 1/e)/2 через класс 5",
         lambda: (_v(b("e")) + 1.0 / _v(b("e"))) / 2),
        ("x − x", "вычитание класса 2",
         lambda: _v(b("root:2:3") - b("root:2:3")),
         "ноль, класс 1", lambda: _v(b("0"))),
        ("x/x", "деление класса 11", lambda: _v(b("pi") / b("pi")),
         "единица, класс 1", lambda: _v(b("1"))),
    )


def report() -> list[dict]:
    out = []
    for name, wa, fa, wb, fb in _pairs():
        rec = {"величина": name, "путь A": wa, "путь B": wb}
        try:
            a = fa()
            rec["A"] = a
        except Exception as ex:
            rec["A"] = None
            rec["итог"] = f"путь A упал: {type(ex).__name__}"
            out.append(rec)
            continue
        try:
            b = fb()
            rec["B"] = b
        except Exception as ex:
            rec["B"] = None
            rec["итог"] = f"путь B упал: {type(ex).__name__}"
            out.append(rec)
            continue
        rec["итог"] = "сошлось" if abs(a - b) <= TOL else "РАСХОЖДЕНИЕ"
        rec["разность"] = abs(a - b)
        out.append(rec)
    return out


def agreed() -> int:
    return sum(1 for r in report() if r["итог"] == "сошлось")


def diverged() -> list[str]:
    return [r["величина"] for r in report() if r["итог"] != "сошлось"]


def pairs_count() -> int:
    return len(_pairs())


def divergence_is_found() -> bool:
    """Двойственность работает: она поймала расхождение, а не подтвердила всё.

    Пары, которые сходятся всегда, ничего не проверяют. Ценность приёма
    в том, что часть пар расходится, и расхождение указывает на место.
    """
    return len(diverged()) > 0


def the_tolerance_is_declared() -> bool:
    return TOL > 0


NOT_MEASURED = (
    "Пути сравниваются по значению `float`, а не символьно: совпадение "
    "значений не есть тождество величин (две разные величины совпадают "
    "как float64 побитно — замер STUDY-04 §1).",
    "Ошибка, лежащая в общем для обоих путей коде, проходит оба и не "
    "ловится. Двойственность ловит расхождение реализаций, а не общую "
    "неправду.",
    "Расхождение указывает место, но не виновного: который из двух путей "
    "неверен, пара не решает.",
    "Пар восемь, и все — на достижимых классах. Для 11 недостижимых "
    "классов пары нет вовсе.",
)


if __name__ == "__main__":
    for r in report():
        a = f"{r['A']:.12g}" if isinstance(r.get("A"), float) else "—"
        b = f"{r['B']:.12g}" if isinstance(r.get("B"), float) else "—"
        print(f"{r['величина']:<10} {r['итог']:<14} A={a:<18} B={b:<18} "
              f"({r['путь A']} / {r['путь B']})")
    print()
    print(f"пар: {pairs_count()}, сошлось: {agreed()}, "
          f"разошлось: {len(diverged())} — {', '.join(diverged()) or '—'}")
