#!/usr/bin/env python3
"""Фасад среды: величина строится по краткой записи, а не по памяти о API.

Изнутри среды первое действие — не «посчитать», а **назвать класс и
рациональные параметры**. Фасад держит краткий язык записи, чтобы это
dействие стоило одну строку:

    7            класс 1, натуральная степень
    1/3          класс 1 из рационального
    root:2:3     класс 2, корень степени 2 из 3   (√3)
    unity:8:1    класс 4, корень ��з единицы       (ζ₈)
    e            класс 5
    ln:2         класс 6
    cos:1/3      класс 10, аргумент в долях π
    sin:1/6      класс 10
    pi           класс 11 (4·arctan 1)
    cosh:1       класс 12
    sinh:1       класс 12

Того, чего в библиотеке нет, здесь нет тоже: язык не изображает
полноты. Что именно отсутствует — считает `arith.py`, а не этот файл.

    python3 -B ENV/AIAΣA/bin/phi.py pi e root:2:3
"""

from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

BIN = Path(__file__).resolve().parent
sys.path.insert(0, str(BIN))

import assemble  # noqa: E402

_LIB = None


def lib():
    """Собранная библиотека. Собирается при первом обращении."""
    global _LIB
    if _LIB is None:
        d = assemble.build_dir()
        if str(d) not in sys.path:
            sys.path.insert(0, str(d))
        import libphi  # noqa: E402
        _LIB = libphi
    return _LIB


def _frac(s: str) -> Fraction:
    return Fraction(s)


#: Краткая запись → постр��ение. Ключ без двоеточия — величина без
#: параметров; с двоеточием — параметры через двоеточие.
def build_term(spec: str):
    """Построить инвариант по краткой записи. Неизвестное — ValueError."""
    L = lib()
    IT = L.InvariantTerm
    head, _, rest = spec.partition(":")
    args = rest.split(":") if rest else []
    if head == "pi":
        return IT.pi()
    if head == "e":
        return IT.e()
    if head == "root":
        return IT.algebraic_root(int(args[0]), _frac(args[1]))
    if head == "unity":
        return IT.unity_root(int(args[0]), int(args[1]) if len(args) > 1 else 1)
    if head == "ln":
        return IT.ln(_frac(args[0]))
    if head == "log":
        return IT.log(_frac(args[0]), _frac(args[1]))
    if head == "cos":
        return IT.cos(_frac(args[0]))
    if head == "sin":
        return IT.sin(_frac(args[0]))
    if head == "cosh":
        return IT.cosh(_frac(args[0]))
    if head == "sinh":
        return IT.sinh(_frac(args[0]))
    try:
        return IT.from_rational(Fraction(spec))
    except (ValueError, ZeroDivisionError):
        raise ValueError(f"нет такой записи: {spec!r}")


#: Записи, покрытые фасадом. Пополняется только вместе с библиотекой.
SPECS = ("7", "1/3", "root:2:3", "unity:8:1", "e", "ln:2", "cos:1/3",
         "sin:1/6", "pi", "cosh:1", "sinh:1", "log:2:8")


def describe(term) -> dict:
    """Что известно об инварианте: класс, параметры, представление."""
    L = lib()
    cls = L.NumberClass(int(term.number_class))
    return {"класс": int(cls), "имя": cls.name,
            "параметры": {k: str(v) for k, v in
                          (term.class_params or {}).items()
                          if k != "sympy_obj"},
            "запись": repr(term)}


def every_spec_builds() -> bool:
    """Каждая объявленная запись строится. Объявить и не построить нельзя."""
    for s in SPECS:
        try:
            build_term(s)
        except Exception:
            return False
    return True


def unknown_spec_is_refused() -> bool:
    """Неизвестная запись отвергается, а не подставляется ближайшей."""
    try:
        build_term("bessel:1")
    except ValueError:
        return True
    return False


def classes_reached() -> set[int]:
    """Какие классы вообще достижимы через фасад. Замер, не обещание."""
    out = set()
    for s in SPECS:
        try:
            out.add(int(build_term(s).number_class))
        except Exception:
            pass
    return out


if __name__ == "__main__":
    specs = sys.argv[1:] or list(SPECS)
    for s in specs:
        try:
            d = describe(build_term(s))
            print(f"{s:<12} класс {d['класс']:>2} {d['имя']:<16} {d['параметры']}")
        except ValueError as ex:
            print(f"{s:<12} ОТКАЗ: {ex}")
