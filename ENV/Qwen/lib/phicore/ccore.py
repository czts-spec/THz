"""Связь с ядром на C (этап 3 плана) — ctypes, без заголовков Python.

Ядро `ccore.c` собирается по требованию и кладётся в кэш по хешу
исходника: при правке исходника собирается заново, при неизменном —
берётся готовое. Заголовков Python оно не знает нарочно: так тот же
алгоритм годится этапам 4 (Rust) и 5 (FPGA), а связь остаётся тонкой.

Компилятора может не быть. Тогда это ОТКАЗ поимённо (`Unresolved`), а не
тихий переход на Python: подмена быстрого ядра медленным без ведома
вызывающего — то самое молчание, которого в устройстве нет.

Коды возврата C переводятся в те же отказы, что поднимает эталон
(`factored.py`): переполнение носителя и показателя — событие ширины,
чётный корень из отрицательного — `Unresolved`. Совпадение эталона и
ядра проверяется побитово, по канону, а не «с точностью до».
"""
from __future__ import annotations

import ctypes
import hashlib
import pathlib
import subprocess
import tempfile

from .errors import ExactRepresentationNotAvailable, Unresolved
from .factored import Factored

HERE = pathlib.Path(__file__).resolve().parent
SOURCE = HERE / "ccore.c"
CC = ("cc", "gcc", "clang")
FLAGS = ("-O2", "-std=c11", "-Wall", "-Wextra", "-fPIC", "-shared")

PHI_OK, PHI_E_SUPPORT, PHI_E_DEN, PHI_E_EXP, PHI_E_EVEN_ROOT, PHI_E_INPUT = range(6)
CODE_NAMES = {
    PHI_E_SUPPORT: "носитель шире ёмкости ядра",
    PHI_E_DEN: "общий знаменатель показателей за int32",
    PHI_E_EXP: "показатель за int32",
    PHI_E_EVEN_ROOT: "чётная степень корня из отрицательного",
    PHI_E_INPUT: "негодный вход",
}
CAP = 16                                   # обязано совпасть с PHI_CAP в C


class PhiT(ctypes.Structure):
    _fields_ = [("sign", ctypes.c_int32), ("den", ctypes.c_int32),
                ("n", ctypes.c_int32),
                ("prime", ctypes.c_int32 * CAP), ("num", ctypes.c_int32 * CAP)]


_LIB = None


def _which_cc() -> str | None:
    import shutil
    for c in CC:
        if shutil.which(c):
            return c
    return None


def library():
    """Собранное ядро; отказ поимённо, если собрать нечем."""
    global _LIB
    if _LIB is not None:
        return _LIB
    src = SOURCE.read_bytes()
    tag = hashlib.sha256(src).hexdigest()[:16]
    so = pathlib.Path(tempfile.gettempdir()) / f"phicore-ccore-{tag}.so"
    if not so.is_file():
        cc = _which_cc()
        if cc is None:
            raise Unresolved("сборка ядра", "ccore.c",
                             hint=f"компилятора нет ({', '.join(CC)}): ядро на C "
                                  f"не собрано; молча на Python не переходим")
        # Наборы гоняются в несколько процессов: сборка идёт во временное
        # имя и переставляется на место одним os.replace — так двое
        # собирающих не портят файл друг другу.
        import os
        tmp = so.with_name(f"{so.name}.{os.getpid()}.tmp")
        r = subprocess.run([cc, *FLAGS, "-o", str(tmp), str(SOURCE)],
                           capture_output=True, text=True)
        if r.returncode != 0:
            tmp.unlink(missing_ok=True)
            raise Unresolved("сборка ядра", "ccore.c",
                             hint=f"{cc} вернул {r.returncode}: {r.stderr.strip()[:300]}")
        os.replace(tmp, so)
    lib = ctypes.CDLL(str(so))
    p = ctypes.POINTER(PhiT)
    lib.phi_capacity.restype = ctypes.c_int
    lib.phi_sizeof.restype = ctypes.c_int
    lib.phi_set.argtypes = [p, ctypes.c_int32, ctypes.c_int32,
                            ctypes.POINTER(ctypes.c_int32),
                            ctypes.POINTER(ctypes.c_int32), ctypes.c_int32]
    for name in ("phi_mul", "phi_div"):
        getattr(lib, name).argtypes = [p, p, p]
    lib.phi_root.argtypes = [p, ctypes.c_int32, p]
    lib.phi_pow.argtypes = [p, ctypes.c_int32, p]
    lib.phi_inverse.argtypes = [p, p]
    lib.phi_bench_chain.argtypes = [p, p, p, ctypes.c_int32]
    lib.phi_merge_pairs.argtypes = [p, p, p, ctypes.POINTER(ctypes.c_int32),
                                    ctypes.c_int32, ctypes.c_int32]
    lib.phi_eq.argtypes = [p, p]
    if lib.phi_capacity() != CAP or lib.phi_sizeof() != ctypes.sizeof(PhiT):
        raise Unresolved("связь с ядром", "ccore.c",
                         hint=f"ёмкость {lib.phi_capacity()} против {CAP} или размер "
                              f"{lib.phi_sizeof()} против {ctypes.sizeof(PhiT)}: "
                              f"описания разошлись")
    _LIB = lib
    return lib


def available() -> bool:
    try:
        library()
        return True
    except Unresolved:
        return False


def _raise(code: int, what: str, op: str = "") -> None:
    name = CODE_NAMES.get(code, f"код {code}")
    if code == PHI_E_EVEN_ROOT:
        raise Unresolved(op or "√", what, hint=name)
    raise ExactRepresentationNotAvailable(f"{what} {op}".strip() + f": {name}")


def to_c(v: Factored) -> PhiT:
    """Эталон → структура ядра. Ширина сверяется здесь, а не внутри."""
    out = PhiT()
    n = len(v.terms)
    primes = (ctypes.c_int32 * max(n, 1))(*[p for p, _ in v.terms])
    nums = (ctypes.c_int32 * max(n, 1))(*[a for _, a in v.terms])
    code = library().phi_set(ctypes.byref(out), v.sign, v.den, primes, nums, n)
    if code != PHI_OK:
        _raise(code, v.record())
    return out


def from_c(s: PhiT) -> Factored:
    return Factored(s.sign, s.den, {s.prime[i]: s.num[i] for i in range(s.n)})


def _binary(name: str, a: Factored, b: Factored) -> Factored:
    ca, cb, out = to_c(a), to_c(b), PhiT()
    code = getattr(library(), name)(ctypes.byref(ca), ctypes.byref(cb), ctypes.byref(out))
    if code != PHI_OK:
        _raise(code, f"{a.record()} и {b.record()}", op=name)
    return from_c(out)


def mul(a: Factored, b: Factored) -> Factored:
    return _binary("phi_mul", a, b)


def div(a: Factored, b: Factored) -> Factored:
    return _binary("phi_div", a, b)


def root(a: Factored, k: int) -> Factored:
    ca, out = to_c(a), PhiT()
    code = library().phi_root(ctypes.byref(ca), k, ctypes.byref(out))
    if code != PHI_OK:
        _raise(code, a.record(), op=f"√{k}")
    return from_c(out)


def power(a: Factored, k: int) -> Factored:
    ca, out = to_c(a), PhiT()
    code = library().phi_pow(ctypes.byref(ca), k, ctypes.byref(out))
    if code != PHI_OK:
        _raise(code, a.record(), op=f"^{k}")
    return from_c(out)


def inverse(a: Factored) -> Factored:
    ca, out = to_c(a), PhiT()
    code = library().phi_inverse(ctypes.byref(ca), ctypes.byref(out))
    if code != PHI_OK:
        _raise(code, a.record(), op="^-1")
    return from_c(out)


def pack(pairs, conv=None):
    """Перевод партии в структуры ядра — отдельно от вызова.

    Разделение не косметическое: замер 04.09 показал, что граница — это
    ПЕРЕВОД значения, а не переход в библиотеку. Пока перевод сидел
    внутри вызова, партия не давала выигрыша и причина была не видна.

    Одноимённые двери соседей — ДРУГИЕ понятия: `tag.pack` кладёт в байт
    адрес таблицы состояний (индекс и код α), `wheels.WheelByte.pack` — колесо
    сетки ширины 1 в байт по правилам. Здесь ни байта, ни колеса: массивы
    ctypes для перехода границы C.
    """
    conv = to_c if conv is None else conv
    n = len(pairs)
    arr = PhiT * max(n, 1)
    ca, cb, out = arr(), arr(), arr()
    codes = (ctypes.c_int32 * max(n, 1))()
    for i, (a, b) in enumerate(pairs):
        ca[i] = conv(a)
        cb[i] = conv(b)
    return ca, cb, out, codes, n


def call_packed(packed, sgn: int = 1, lib=None, fn: str = "phi_merge_pairs") -> None:
    """Только переход границы: партия считается, итог остаётся в структурах."""
    ca, cb, out, codes, n = packed
    lib = library() if lib is None else lib
    rc = getattr(lib, fn)(ca, cb, out, codes, n, sgn)
    if rc != PHI_OK:
        _raise(rc, f"партия из {n} пар", op=fn)


def unpack(packed) -> list:
    """Только обратный перевод: структуры ядра → значения эталона.

    Не `tag.unpack` (байт → состояние и код α) и не `wheels.WheelByte.unpack`
    (байт → колесо): вход здесь — структуры C, а не байт.
    """
    _, _, out, codes, n = packed
    return [from_c(out[i]) if codes[i] == PHI_OK else int(codes[i]) for i in range(n)]


def merge_packed(packed, sgn: int = 1, lib=None, fn: str = "phi_merge_pairs"):
    """Партия переведённых пар: переход границы и обратный перевод."""
    call_packed(packed, sgn, lib, fn)
    return unpack(packed)


def merge_pairs(pairs, sgn: int = 1, lib=None, fn: str = "phi_merge_pairs",
                conv=None):
    """Партия слияний: [(a, b), …] → [Factored | код возврата].

    Значение, если пара сведена, и код (int), если у неё событие ширины.
    Партия НЕ обрывается на первом событии — иначе пропали бы итоги всех
    последующих пар.
    """
    return merge_packed(pack(pairs, conv), sgn, lib, fn)


def bench_chain(a: Factored, b: Factored, reps: int) -> Factored:
    """Цепочка acc = (acc·b)/b, reps раз, одним вызовом — 2·reps действий.

    Итерации зависимы нарочно: независимый повтор компилятор вправе
    свернуть в одну итерацию, и замер солгал бы. Итог обязан вернуться к
    a — верность замера проверяется тем же прогоном.
    """
    ca, cb, out = to_c(a), to_c(b), PhiT()
    code = library().phi_bench_chain(ctypes.byref(ca), ctypes.byref(cb),
                                     ctypes.byref(out), reps)
    if code != PHI_OK:
        _raise(code, f"{a.record()} и {b.record()}", op="phi_bench_chain")
    return from_c(out)


def equal(a: Factored, b: Factored) -> bool:
    ca, cb = to_c(a), to_c(b)
    return bool(library().phi_eq(ctypes.byref(ca), ctypes.byref(cb)))


def bits(v: Factored) -> bytes:
    """Байты структуры ядра — то, что сравнивается ПОБИТОВО."""
    return bytes(memoryview(to_c(v)).cast("B"))


# ── Колёса и наследственная форма: обратное преобразование на ядре ──────────
# Слово автора 24.09 (коррекция d1eba5c2eb): обратное преобразование — из
# позиционной сокращённой записи (машинного целого) в сумму колёс — выполняется
# В ЯДРЕ; эталон — phicore/wheels.py; сверка побитовая.

def _wheel_lib():
    lib = library()
    if not getattr(lib, "_wheels_bound", False):
        i32p, i64p = ctypes.POINTER(ctypes.c_int32), ctypes.POINTER(ctypes.c_int64)
        lib.phi_wheels.argtypes = [ctypes.c_int64, ctypes.c_int32, i32p, i32p, ctypes.c_int32]
        lib.phi_wheels.restype = ctypes.c_int32
        lib.phi_wheels_batch.argtypes = [i64p, ctypes.c_int32, ctypes.c_int32, i32p, i32p, i32p, ctypes.c_int64]
        lib.phi_wheels_batch.restype = ctypes.c_int64
        lib.phi_hereditary.argtypes = [ctypes.c_int64, ctypes.c_int32, ctypes.c_int32, i64p, i32p, i64p]
        lib.phi_hereditary.restype = ctypes.c_int32
        lib.phi_goodstein.argtypes = [ctypes.c_int64, ctypes.c_int32, i64p]
        lib.phi_goodstein.restype = ctypes.c_int32
        lib._wheels_bound = True
    return lib


def _wheel_code(rc: int, what: str):
    if rc == -1:
        raise ValueError(f"{what}: негодный вход")
    if rc == -2:
        raise ExactRepresentationNotAvailable(f"{what}: буфер колёс мал")
    if rc == -3:
        raise ExactRepresentationNotAvailable(f"{what}: переполнение int64 при чтении дерева в новой базе")
    raise ExactRepresentationNotAvailable(f"{what}: код ядра {rc}")


def wheels_of(n: int, base: int, cap: int = 64) -> list:
    """Колёса |n| базы b на ядре: [(m, e), …] от младшего; [] для 0."""
    lib = _wheel_lib()
    m = (ctypes.c_int32 * cap)()
    e = (ctypes.c_int32 * cap)()
    k = lib.phi_wheels(int(n), int(base), m, e, cap)
    if k < 0:
        _wheel_code(k, "phi_wheels")
    return [(m[i], e[i]) for i in range(k)]


def wheels_batch(ns, base: int) -> list:
    """Партия: список колёс для каждого n одним вызовом ядра."""
    lib = _wheel_lib()
    ns = [int(x) for x in ns]
    count = len(ns)
    cap = 64 * max(1, count)
    n_arr = (ctypes.c_int64 * max(1, count))(*ns)
    m = (ctypes.c_int32 * cap)()
    e = (ctypes.c_int32 * cap)()
    cnt = (ctypes.c_int32 * max(1, count))()
    total = lib.phi_wheels_batch(n_arr, count, int(base), m, e, cnt, cap)
    if total < 0:
        _wheel_code(int(total), "phi_wheels_batch")
    out, off = [], 0
    for i in range(count):
        out.append([(m[off + j], e[off + j]) for j in range(cnt[i])])
        off += cnt[i]
    return out


def hereditary(n: int, base: int, newbase: int | None = None) -> tuple:
    """(значение дерева в newbase, глубина, колёс всех уровней) — на ядре."""
    lib = _wheel_lib()
    v, d, w = ctypes.c_int64(), ctypes.c_int32(), ctypes.c_int64()
    rc = lib.phi_hereditary(int(n), int(base), int(base if newbase is None else newbase),
                            ctypes.byref(v), ctypes.byref(d), ctypes.byref(w))
    if rc:
        _wheel_code(rc, "phi_hereditary")
    return v.value, d.value, w.value


def goodstein(n: int, steps: int) -> list:
    """Последовательность Гудстейна на ядре (останов на нуле; ≤ steps + 1 членов)."""
    lib = _wheel_lib()
    out = (ctypes.c_int64 * (int(steps) + 1))()
    k = lib.phi_goodstein(int(n), int(steps), out)
    if k < 0:
        _wheel_code(k, "phi_goodstein")
    return [out[i] for i in range(k)]
