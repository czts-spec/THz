#!/usr/bin/env python3
"""Каталог двадцати четырёх классов, сверенный с библиотекой, а не переписанный.

Основная система координат среды — не биты, а класс: природа величины.
Каталог здесь записан отдельно от библиотеки нарочно: записанное и
собранное сверяются, и расхождение называется, а не сглаживается.

    python3 -B ENV/AIAΣA/bin/classes.py
"""

from __future__ import annotations

import sys
from pathlib import Path

BIN = Path(__file__).resolve().parent
sys.path.insert(0, str(BIN))

import phi  # noqa: E402

#: Записанный каталог: номер → (имя в библиотеке, чем задаётся).
CATALOG = {
    0:  ("ZERO", "онтологическая пустота, предшествует числу"),
    1:  ("NATURAL_POWER", "x = a^n, a, n ∈ ℕ"),
    2:  ("ALGEBRAIC_ROOT", "x^n = a, a ∈ ℚ, n ∈ ℕ"),
    3:  ("POLY_ROOT", "P(x) = 0, P ∈ ℤ[x]"),
    4:  ("UNITY_ROOT", "x^n = 1"),
    5:  ("EXPONENTIAL", "x = e^y"),
    6:  ("LOGARITHMIC", "x = ln y"),
    7:  ("POWER_SERIES", "x = Σ aₙ rⁿ сходится"),
    8:  ("INFINITE_PRODUCT", "x = Π pₙ^aₙ"),
    9:  ("INTEGRAL", "x = ∫ f(t) dt, f рациональна"),
    10: ("TRIGONOMETRIC", "sin, cos, tan"),
    11: ("ARCTRIG", "arcsin, arccos, arctan"),
    12: ("HYPERBOLIC", "sinh, cosh, tanh"),
    13: ("ZETA", "ζ(s) = Σ n^(−s)"),
    14: ("GAMMA", "Γ(z)"),
    15: ("ELLIPTIC", "K(k), E(k)"),
    16: ("BESSEL_ROOT", "J_n(x) = 0"),
    17: ("LAMBERT_W", "W(y): x·e^x = y"),
    18: ("DE_SOLUTION", "решения дифференциальных уравнений"),
    19: ("ITERATION_LIMIT", "lim f^(n)(x₀)"),
    20: ("UNCOMPUTABLE", "познанный предел, Ω Хайтина"),
    21: ("INFINITESIMAL", "a·ε"),
    22: ("INFINITY", "a·ω"),
    23: ("IGNORANCE", "навык: float/nan/inf как есть — 23 на входе, 23 на выходе, не мешать"),
}

#: Разбиение каталога на области. Границы — не украшение: у областей
#: разные правила, и смешение их даёт ложные утверждения.
BANDS = {
    "пустота": (0, 0),
    "вычислимые": (1, 19),
    "предел": (20, 20),
    "нестандартные": (21, 22),
    "нативный IEEE 754": (23, 23),
}


def library_names() -> dict[int, str]:
    NC = phi.lib().NumberClass
    return {int(c): c.name for c in NC}


def mismatches() -> list[str]:
    """Расхождения записанного каталога и собранной библиотеки."""
    lib = library_names()
    out = []
    for num, (name, _) in sorted(CATALOG.items()):
        if num not in lib:
            out.append(f"{num}: в библиотеке нет")
        elif lib[num] != name:
            out.append(f"{num}: записано {name}, в библиотеке {lib[num]}")
    for num, name in sorted(lib.items()):
        if num not in CATALOG:
            out.append(f"{num}: {name} есть в библиотеке и не записан")
    return out


def catalogue_matches_the_library() -> bool:
    return not mismatches()


def bands_cover_the_catalogue() -> bool:
    """Области покрывают каталог целиком и не пересекаются."""
    seen: set[int] = set()
    for lo, hi in BANDS.values():
        band = set(range(lo, hi + 1))
        if band & seen:
            return False
        seen |= band
    return seen == set(CATALOG)


def computable_count() -> int:
    lo, hi = BANDS["вычислимые"]
    return hi - lo + 1


def every_class_says_how_it_is_given() -> bool:
    """У класса записано, чем он задаётся. Имя без этого — ярлык."""
    return all(len(how) >= 10 for _, how in CATALOG.values())


def zero_is_not_uncomputable() -> bool:
    """Класс 0 и класс 20 арифметически родня, но происхождением разные.

    Смешивать их нельзя: `ZERO` — «ещё не познано», `UNCOMPUTABLE` —
    «познано как непознаваемое». Проверяется тем, что оба есть и это
    разные номера с разными именами.
    """
    return (CATALOG[0][0] != CATALOG[20][0]
            and "предшеств" in CATALOG[0][1] and "предел" in CATALOG[20][1])


NOT_MEASURED = (
    "Сверяются номер и имя. Что библиотека вкладывает в класс по существу, "
    "отсюда не проверяется: совпадение имени не есть совпадение смысла.",
    "Класс 23 в библиотеке объявлен отдельным слоем (PhiIgnorance), а не "
    "ядром перечисления; каталог этого различения не несёт.",
)


if __name__ == "__main__":
    lib = library_names()
    for num, (name, how) in sorted(CATALOG.items()):
        mark = " " if lib.get(num) == name else "×"
        band = next(b for b, (lo, hi) in BANDS.items() if lo <= num <= hi)
        print(f"{mark} {num:>2} {name:<17} {band:<14} {how}")
    print()
    print("каталог сходится с библиотекой:", catalogue_matches_the_library())
    for m in mismatches():
        print("  расхождение:", m)
