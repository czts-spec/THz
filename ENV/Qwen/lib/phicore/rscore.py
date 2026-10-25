"""Связь с ядром на Rust (этап 4 плана) — ctypes, тот же приём, что у C.

Структура значения и коды возврата — общие с `ccore` (там же и разбор
кодов в отказы). Отличие этапа не в скорости: она та же. Отличие в том,
что отсутствие небезопасного кода в арифметике здесь ПРОВЕРЕНО
компилятором — модуль `math` в `rscore.rs` объявлен
`#![forbid(unsafe_code)]`, и `unsafe` остаётся только в тонком слое
связи, который ничего не считает.

Rust может быть не установлен. Тогда это отказ поимённо, а не тихий
переход на C или на Python.
"""
from __future__ import annotations

import ctypes
import hashlib
import os
import pathlib
import shutil
import subprocess
import tempfile

from .ccore import CODE_NAMES, PHI_OK, PhiT, from_c, to_c  # noqa: F401
from .ccore import _raise as _raise_code
from .errors import Unresolved
from .factored import Factored

HERE = pathlib.Path(__file__).resolve().parent
SOURCE = HERE / "rscore.rs"
FLAGS = ("--crate-type=cdylib", "-O", "--edition", "2021")

_LIB = None


def library():
    """Собранное ядро на Rust; отказ поимённо, если собрать нечем."""
    global _LIB
    if _LIB is not None:
        return _LIB
    src = SOURCE.read_bytes()
    tag = hashlib.sha256(src).hexdigest()[:16]
    so = pathlib.Path(tempfile.gettempdir()) / f"phicore-rscore-{tag}.so"
    if not so.is_file():
        rustc = shutil.which("rustc")
        if rustc is None:
            raise Unresolved("сборка ядра", "rscore.rs",
                             hint="rustc не найден: ядро на Rust не собрано; "
                                  "молча на C или Python не переходим")
        tmp = so.with_name(f"{so.name}.{os.getpid()}.tmp")
        r = subprocess.run([rustc, *FLAGS, "-o", str(tmp), str(SOURCE)],
                           capture_output=True, text=True)
        if r.returncode != 0:
            pathlib.Path(tmp).unlink(missing_ok=True)
            raise Unresolved("сборка ядра", "rscore.rs",
                             hint=f"rustc вернул {r.returncode}: {r.stderr.strip()[:300]}")
        os.replace(tmp, so)
    lib = ctypes.CDLL(str(so))
    p = ctypes.POINTER(PhiT)
    lib.phir_capacity.restype = ctypes.c_int
    lib.phir_sizeof.restype = ctypes.c_int
    lib.phir_set.argtypes = [p, ctypes.c_int32, ctypes.c_int32,
                             ctypes.POINTER(ctypes.c_int32),
                             ctypes.POINTER(ctypes.c_int32), ctypes.c_int32]
    for name in ("phir_mul", "phir_div"):
        getattr(lib, name).argtypes = [p, p, p]
    lib.phir_root.argtypes = [p, ctypes.c_int32, p]
    lib.phir_pow.argtypes = [p, ctypes.c_int32, p]
    lib.phir_inverse.argtypes = [p, p]
    lib.phir_eq.argtypes = [p, p]
    lib.phir_bench_chain.argtypes = [p, p, p, ctypes.c_int32]
    lib.phir_merge_pairs.argtypes = [p, p, p, ctypes.POINTER(ctypes.c_int32),
                                     ctypes.c_int32, ctypes.c_int32]
    if lib.phir_sizeof() != ctypes.sizeof(PhiT):
        raise Unresolved("связь с ядром", "rscore.rs",
                         hint=f"размер структуры {lib.phir_sizeof()} против "
                              f"{ctypes.sizeof(PhiT)}: описания разошлись")
    _LIB = lib
    return lib


def available() -> bool:
    try:
        library()
        return True
    except Unresolved:
        return False


def to_rs(v: Factored) -> PhiT:
    """Эталон → структура ядра Rust (тот же вид, своя проверка ширины)."""
    out = PhiT()
    n = len(v.terms)
    primes = (ctypes.c_int32 * max(n, 1))(*[p for p, _ in v.terms])
    nums = (ctypes.c_int32 * max(n, 1))(*[a for _, a in v.terms])
    code = library().phir_set(ctypes.byref(out), v.sign, v.den, primes, nums, n)
    if code != PHI_OK:
        _raise_code(code, v.record())
    return out


def _binary(name: str, a: Factored, b: Factored) -> Factored:
    ca, cb, out = to_rs(a), to_rs(b), PhiT()
    code = getattr(library(), name)(ctypes.byref(ca), ctypes.byref(cb), ctypes.byref(out))
    if code != PHI_OK:
        _raise_code(code, f"{a.record()} и {b.record()}", op=name)
    return from_c(out)


def mul(a: Factored, b: Factored) -> Factored:
    return _binary("phir_mul", a, b)


def div(a: Factored, b: Factored) -> Factored:
    return _binary("phir_div", a, b)


def root(a: Factored, k: int) -> Factored:
    ca, out = to_rs(a), PhiT()
    code = library().phir_root(ctypes.byref(ca), k, ctypes.byref(out))
    if code != PHI_OK:
        _raise_code(code, a.record(), op=f"√{k}")
    return from_c(out)


def power(a: Factored, k: int) -> Factored:
    ca, out = to_rs(a), PhiT()
    code = library().phir_pow(ctypes.byref(ca), k, ctypes.byref(out))
    if code != PHI_OK:
        _raise_code(code, a.record(), op=f"^{k}")
    return from_c(out)


def inverse(a: Factored) -> Factored:
    ca, out = to_rs(a), PhiT()
    code = library().phir_inverse(ctypes.byref(ca), ctypes.byref(out))
    if code != PHI_OK:
        _raise_code(code, a.record(), op="^-1")
    return from_c(out)


def equal(a: Factored, b: Factored) -> bool:
    ca, cb = to_rs(a), to_rs(b)
    return bool(library().phir_eq(ctypes.byref(ca), ctypes.byref(cb)))


def bench_chain(a: Factored, b: Factored, reps: int) -> Factored:
    ca, cb, out = to_rs(a), to_rs(b), PhiT()
    code = library().phir_bench_chain(ctypes.byref(ca), ctypes.byref(cb),
                                      ctypes.byref(out), reps)
    if code != PHI_OK:
        _raise_code(code, f"{a.record()} и {b.record()}", op="phir_bench_chain")
    return from_c(out)


def merge_pairs(pairs, sgn: int = 1):
    """Партия слияний за один переход границы — то же, что у C.

    Раскладка структуры и разбор кодов общие, потому берётся тот же
    сборщик партии: своя библиотека, своё имя функции, свой перевод.
    """
    from .ccore import merge_pairs as _pairs
    return _pairs(pairs, sgn, lib=library(), fn="phir_merge_pairs", conv=to_rs)


def bits(v: Factored) -> bytes:
    """Байты структуры ядра Rust — то, что сравнивается ПОБИТОВО."""
    return bytes(memoryview(to_rs(v)).cast("B"))


# ── Колёса и наследственная форма на ядре Rust (обратное преобразование) ─────
# Та же раскладка вызовов, что у ccore.wheels_of / hereditary / goodstein;
# счёт — в модуле `wheels_math` под `#![forbid(unsafe_code)]`.

def _wheel_lib():
    lib = library()
    if not getattr(lib, "_wheels_bound", False):
        i32p, i64p = ctypes.POINTER(ctypes.c_int32), ctypes.POINTER(ctypes.c_int64)
        lib.phir_wheels.argtypes = [ctypes.c_int64, ctypes.c_int32, i32p, i32p, ctypes.c_int32]
        lib.phir_wheels.restype = ctypes.c_int32
        lib.phir_wheels_batch.argtypes = [i64p, ctypes.c_int32, ctypes.c_int32, i32p, i32p, i32p, ctypes.c_int64]
        lib.phir_wheels_batch.restype = ctypes.c_int64
        lib.phir_hereditary.argtypes = [ctypes.c_int64, ctypes.c_int32, ctypes.c_int32, i64p, i32p, i64p]
        lib.phir_hereditary.restype = ctypes.c_int32
        lib.phir_goodstein.argtypes = [ctypes.c_int64, ctypes.c_int32, i64p]
        lib.phir_goodstein.restype = ctypes.c_int32
        lib._wheels_bound = True
    return lib


def wheels_of(n: int, base: int, cap: int = 64) -> list:
    from .ccore import _wheel_code
    lib = _wheel_lib()
    m = (ctypes.c_int32 * cap)()
    e = (ctypes.c_int32 * cap)()
    k = lib.phir_wheels(int(n), int(base), m, e, cap)
    if k < 0:
        _wheel_code(k, "phir_wheels")
    return [(m[i], e[i]) for i in range(k)]


def wheels_batch(ns, base: int) -> list:
    from .ccore import _wheel_code
    lib = _wheel_lib()
    ns = [int(x) for x in ns]
    count = len(ns)
    cap = 64 * max(1, count)
    n_arr = (ctypes.c_int64 * max(1, count))(*ns)
    m = (ctypes.c_int32 * cap)()
    e = (ctypes.c_int32 * cap)()
    cnt = (ctypes.c_int32 * max(1, count))()
    total = lib.phir_wheels_batch(n_arr, count, int(base), m, e, cnt, cap)
    if total < 0:
        _wheel_code(int(total), "phir_wheels_batch")
    out, off = [], 0
    for i in range(count):
        out.append([(m[off + j], e[off + j]) for j in range(cnt[i])])
        off += cnt[i]
    return out


def hereditary(n: int, base: int, newbase: int | None = None) -> tuple:
    from .ccore import _wheel_code
    lib = _wheel_lib()
    v, d, w = ctypes.c_int64(), ctypes.c_int32(), ctypes.c_int64()
    rc = lib.phir_hereditary(int(n), int(base), int(base if newbase is None else newbase),
                             ctypes.byref(v), ctypes.byref(d), ctypes.byref(w))
    if rc:
        _wheel_code(rc, "phir_hereditary")
    return v.value, d.value, w.value


def goodstein(n: int, steps: int) -> list:
    from .ccore import _wheel_code
    lib = _wheel_lib()
    out = (ctypes.c_int64 * (int(steps) + 1))()
    k = lib.phir_goodstein(int(n), int(steps), out)
    if k < 0:
        _wheel_code(k, "phir_goodstein")
    return [out[i] for i in range(k)]
