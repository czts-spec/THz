"""Точное целочисленное матричное умножение на быстром ядре (BLAS) с классом на значении.

Слово автора 04.09: «реально высокопроизводительное ядро float32 даже не
начиналось». На этой машине ядро — BLAS (float32 152 GOPS, float64
60–100 GOPS); целые в NumPy идут без ядра (1,8 GOPS) и через object
(0,05 GOPS). Здесь целые считаются НА ядре, точно, при любой ширине:

  преобразование без ошибки (Ozaki, Ogita, Oishi, Rump 2012): каждый
  операнд расщепляется на куски по s бит так, что произведение кусков и
  сумма D таких произведений не превосходят 2^53 — тогда float64-GEMM
  над кусками точен (IEEE 754: целые до 2^53 представимы, сложение и
  умножение таких целых без округления). Куски перемножаются попарно
  (p_x · p_w вызовов GEMM), результаты — целые int64 без заворота; итог
  V = Σ_k parts[k]·2^(s·k).

  класс на значении: контракт — ширина c; |V| ≥ 2^(c−1) даёт класс 22
  (INF) в метке значения, а не заворот молча и не общий флаг. Решение
  принимается оценкой в float64 с запасом и точной сборкой целыми Python
  только для значений вблизи границы (их доля мала); внутри контракта
  сборка идёт в int64 по модулю 2^64 — она точна, когда истинное V в
  диапазоне.

  интервалы (засечки): [lo, hi]·W через центр и радиус — два GEMM вместо
  четырёх; всё в целых, деление на 2 точное (числитель чётен).

Авторство: Ozaki–Ogita–Oishi–Rump 2012 (error-free GEMM на быстрых
ядрах); Ootomo–Ozaki–Yokota 2024 (то же на целочисленных тензорных
ядрах — перспектива, здесь не замерена); Goto–van de Geijn 2008 (GEMM);
Harris 2020 (NumPy); Мур 1966 (интервалы); контракт и класс на
значении — AIAΣA (автор проекта).

NOT_MEASURED — разобран 12.09 на два разных утверждения (строка остаётся
записью о том, что было):

  * **ширина куска s — ЗАКРЫТО точной дорогой.** s не «не более семи», а
    ровно семь при D = 1024 на ядре float32, и это МАКСИМУМ: граница
    расщепления есть целочисленное неравенство D·(2^s − 1)² ≤ 2^b − 1
    (`split_bound`, считается bigint'ами Python, без единого float), оно
    держится при s и ломается при s + 1. Замер по таблице глубин
    `[вычислено здесь]`: D = 1024 даёт s = 21 на float64 (b = 53) и s = 7 на
    float32 (b = 24); D = 4096 — 20 и 6; D = 65536 — 18 и 4.

    ПОПРАВКА замером: «шесть вызовов против одного» — не одно число, а
    функция ширины операнда, и запись, называвшая одно число, названа
    дорогой счёта. Замер самим `split` при D = 1024: на операндах 14 бит
    float64 делает 1 вызов, float32 — 9; на 21 бите 1 против 16; на 31 бите
    4 против 25. Шесть против одного приходится на НЕСИММЕТРИЧНУЮ пару
    (8 бит × 15…20 бит: куски 2 × 3). Число кусков берётся у самого `split`,
    а не у формулы ⌈w/s⌉: знаковый перенос старшего куска добавляет один.

  * **GPU и целочисленные тензорные ядра — ОТКАЗ СРЕДЫ, названный замером.**
    Не «величины не существует», а «посчитать здесь нечем», и это измерено, а
    не заявлено `[вычислено здесь]`: `nvidia-smi` в PATH нет, устройств
    `/dev/nvidia*` нет, `/proc/driver/nvidia` нет, `cupy` нет, torch стоит
    сборкой `+cpu` и сообщает 0 устройств CUDA. Чем мерилось бы: пропускная
    способность int8-тензорных ядер (Ootomo–Ozaki–Yokota 2024) против
    здешних 152 GOPS float32-BLAS, на той же тройке (D, ширина, доля
    значений у границы контракта); дорога расщепления при этом не меняется —
    меняется только b в `split_bound` (int8-ядро копит в int32, то есть
    b = 31 вместо 53).
"""
from __future__ import annotations

import numpy as np

from .errors import ExactRepresentationNotAvailable

INF = 22           # класс: за пределом контракта (переполнение ширины)
RATIONAL = 1       # класс: конечное точное целое в контракте
_F64_EXACT = 53    # бит мантиссы float64 + 1: целые до 2^53 точны
_INT64 = 1 << 63   # ширина ответа: |v| < 2^63, иначе значение в него не входит


def split_bound(depth: int, s: int, exact_bits: int = _F64_EXACT) -> tuple:
    """Граница расщепления ЦЕЛЫМИ: (наибольшая сумма, помещается ли она).

    Сумма depth произведений двух кусков по s бит не превосходит
    depth·(2^s − 1)², и ядро точно ровно тогда, когда она не превосходит
    2^exact_bits − 1. Считается bigint'ами Python — ни одного float, ни
    одного логарифма: именно это утверждение `piece_bits` и обязан
    обеспечивать, и здесь оно проверяется прямо, а не через ⌈log2⌉."""
    if min(depth, s, exact_bits) < 1:
        raise ValueError("depth, s и exact_bits — целые ≥ 1")
    верх = depth * (2 ** s - 1) ** 2
    return верх, верх <= 2 ** exact_bits - 1


def piece_bits(depth: int, exact_bits: int = _F64_EXACT) -> int:
    """s: 2s + ⌈log2 depth⌉ ≤ exact_bits — произведение кусков и их сумма по D точны."""
    d = max(1, int(np.ceil(np.log2(max(depth, 1)))))
    return max(1, (exact_bits - d) // 2)


def split(a: np.ndarray, s: int) -> list:
    """a = Σ_i pieces[i]·2^(s·i); младшие куски в [0, 2^s), старший — знаковый, |·| < 2^(s−1).
    Число кусков считается по наибольшему модулю, без сканирования на каждом шаге."""
    a = np.asarray(a, dtype=np.int64)
    amax = int(np.abs(a).max()) if a.size else 0
    half = 1 << (s - 1)
    if amax < half:
        return [a]
    mask = np.int64((1 << s) - 1)
    pieces = []
    r = a
    while amax >= half:
        lo = r & mask
        pieces.append(lo)
        r = (r - lo) >> s
        amax = (amax >> s) + 1          # граница |r| после переноса
    pieces.append(r)
    return pieces


def _parts(x: np.ndarray, w: np.ndarray, s: int) -> dict:
    """parts[k] = Σ_{i+j=k} x_i @ w_j — каждый GEMM точен в float64, сумма в int64."""
    xs = [p.astype(np.float64) for p in split(x, s)]
    ws = [p.astype(np.float64) for p in split(w, s)]
    parts: dict = {}
    for i, a in enumerate(xs):
        for j, b in enumerate(ws):
            prod = np.rint(a @ b).astype(np.int64)
            parts[i + j] = parts[i + j] + prod if (i + j) in parts else prod
    return parts


def _estimate(parts: dict, s: int) -> np.ndarray:
    est = np.zeros(next(iter(parts.values())).shape, np.float64)
    for k, p in parts.items():
        est += p.astype(np.float64) * float(2.0 ** (s * k))
    return est


def _exact_at(parts: dict, s: int, idx) -> list:
    """Точные значения целыми Python в указанных позициях."""
    out = []
    for pos in zip(*idx):
        out.append(sum(int(p[pos]) << (s * k) for k, p in parts.items()))
    return out


def _wrap_combine(parts: dict, s: int) -> np.ndarray:
    """Σ parts[k]·2^(s·k) по модулю 2^64 (uint64 — определено), как int64."""
    acc = np.zeros(next(iter(parts.values())).shape, np.uint64)
    for k, p in parts.items():
        sh = s * k
        if sh >= 64:
            continue
        acc += (p.astype(np.uint64) << np.uint64(sh))
    return acc.view(np.int64)



def _ints(a, name: str) -> "np.ndarray":
    """Дверь ядра: целые принимаются, прочее отвергается ПОИМЁННО, не усекается.

    `[вычислено здесь]` Замер 11.09 (глава 08 учебника, расхождение 2): вход
    `[[1.5, 2.5]]` молча становился `[[3]]` — `np.asarray(..., dtype=np.int64)`
    усекает и меняет предмет, тогда как соседняя дверь того же дома отказывает
    прямо: `Notch.from_q(1.5, 8)` → «float запрещён построением; явная дверь —
    Q.from_float». Двери дома держатся одинаково.
    """
    arr = np.asarray(a)
    kind = arr.dtype.kind
    if kind == "b":
        raise TypeError(f"{name}: bool не число — дверь принимает целые")
    if kind == "f":
        raise TypeError(f"{name}: float запрещён построением (усечение молча меняет предмет); "
                        "явная дверь — целые либо Q.from_float с объявленной сеткой")
    if kind == "c":
        raise TypeError(f"{name}: комплексное не сведено — дверь принимает целые")
    if kind in ("U", "S"):
        raise TypeError(f"{name}: строка не число — дверь принимает целые")
    if kind == "O":
        flat = arr.reshape(-1)
        for v in flat[:1024]:
            if isinstance(v, bool) or not isinstance(v, int):
                raise TypeError(f"{name}: в массиве не целое ({type(v).__name__}) — дверь принимает целые")
        return arr.astype(np.int64)
    return arr.astype(np.int64)

def exact_gemm(x: np.ndarray, w: np.ndarray, contract_bits: int = 64):
    """Точное x @ w для целых любой ширины (int64 на входе).

    Возвращает (v, tag): v — int64, точное там, где tag == RATIONAL;
    tag == INF там, где |точное| ≥ 2^(contract_bits−1) — за пределом
    контракта (в v там заворот по модулю 2^64, значением не считается).

    ШИРИНА ОТВЕТА — 64 бита, и это граница устройства, названная ОТКАЗОМ
    (12.09, глава 08 §14 п. 1 — дыра устройства, найденная прогоном главы).
    Замер ДО: при `contract_bits=80` и истинном 2^76 = 75557863725914323419136
    ответ был `(0, метка 1)` — заворот по модулю 2^64 (`_wrap_combine`) выдан
    за ТОЧНОЕ значение; при `contract_bits=64` то же место помечалось классом
    22 верно. Теперь: значение, которое лежит ВНУТРИ контракта, но не
    помещается в int64 (|V| ≥ 2^63), не возвращается вовсе — поднимается
    `ExactRepresentationNotAvailable`, называющий позицию, точное значение и
    дверь без ограничения ширины (`exact_gemm_bigint`, целые Python). Класс 22
    остаётся тем, чем был: «за пределом КОНТРАКТА», и молчаливым заворотом
    больше не подменяется. Всё, что помещается в int64, считается как прежде:
    при `contract_bits=80` малое произведение отвечает `(7, метка 1)`.
    """
    x = _ints(x, "exact_gemm x"); w = _ints(w, "exact_gemm w")
    depth = x.shape[-1]
    xmax = int(np.abs(x).max()) if x.size else 0
    wmax = int(np.abs(w).max()) if w.size else 0
    if xmax * wmax * depth < (1 << _F64_EXACT):
        # быстрый путь: один GEMM float64 точен, |V| < 2^53 — в контракте ≥ 54 бит всегда
        v = np.rint(x.astype(np.float64) @ w.astype(np.float64)).astype(np.int64)
        tag = np.full(v.shape, RATIONAL, np.uint8)
        if contract_bits <= _F64_EXACT:
            tag[np.abs(v) >= (1 << (contract_bits - 1))] = INF
        return v, tag
    s = piece_bits(depth)
    parts = _parts(x, w, s)
    est = _estimate(parts, s)
    limit = 2.0 ** (contract_bits - 1)
    n_terms = len(parts)
    # относительная погрешность оценки ≤ n_terms·2^-52 (каждое слагаемое точно, сложение округляет)
    margin = limit * n_terms * 2.0 ** -50
    inside = np.abs(est) < limit - margin
    outside = np.abs(est) > limit + margin
    v = _wrap_combine(parts, s)
    tag = np.where(inside, RATIONAL, INF).astype(np.uint8)
    border = ~inside & ~outside
    if border.any():
        idx = np.nonzero(border)
        for pos, val in zip(zip(*idx), _exact_at(parts, s, idx)):
            over = abs(val) >= (1 << (contract_bits - 1))
            tag[pos] = INF if over else RATIONAL
            if not over:
                v[pos] = val
    if contract_bits > 64:
        # Внутри контракта, но шире ответа: int64 такое значение не несёт, и
        # прежде на его месте молча стоял заворот с меткой «точно» (§14 п. 1).
        # Подозрение берётся по оценке с тем же запасом, что и граница
        # контракта, а решается ТОЧНОЙ сборкой целыми Python — как и на краю.
        suspect = (tag == RATIONAL) & (np.abs(est) >= float(_INT64) - margin)
        if suspect.any():
            idx = np.nonzero(suspect)
            for pos, val in zip(zip(*idx), _exact_at(parts, s, idx)):
                if abs(val) >= _INT64:
                    raise ExactRepresentationNotAvailable(
                        f"exact_gemm: контракт {contract_bits} бит шире ответа "
                        f"(int64, |v| < 2^63), а точное значение в позиции "
                        f"{pos} равно {val} — вернуть его нечем, и заворот по "
                        f"модулю 2^64 значением не считается; дверь без "
                        f"ограничения ширины — exact_gemm_bigint (целые Python)")
    return v, tag


def exact_gemm_bigint(x: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Точный результат целыми Python (object) — эталон и путь класса 1 без ограничения."""
    x = _ints(x, "exact_gemm_bigint x"); w = _ints(w, "exact_gemm_bigint w")
    s = piece_bits(x.shape[-1])
    parts = _parts(x, w, s)
    out = np.zeros(next(iter(parts.values())).shape, dtype=object)
    for k, p in parts.items():
        out = out + p.astype(object) * (1 << (s * k))
    return out


def interval_gemm(lo: np.ndarray, hi: np.ndarray, w: np.ndarray, contract_bits: int = 64):
    """Точные границы [lo, hi] @ w для целых интервалов: lo' = Σ min, hi' = Σ max по
    знаку w — через центр и радиус (два GEMM): 2c = lo+hi, 2r = hi−lo,
    P = 2c @ w, R = 2r @ |w|; lo' = (P − R)/2, hi' = (P + R)/2 — числители чётны.
    Метка INF там, где P или R за пределом контракта."""
    lo = _ints(lo, "interval_gemm lo"); hi = _ints(hi, "interval_gemm hi"); w = _ints(w, "interval_gemm w")
    c2 = lo + hi
    r2 = hi - lo
    P, tp = exact_gemm(c2, w, contract_bits)
    R, tr = exact_gemm(r2, np.abs(w), contract_bits)
    tag = np.where((tp == INF) | (tr == INF), INF, RATIONAL).astype(np.uint8)
    # (P − R) и (P + R) могут переполнить int64 при |P|,|R| ~ 2^62 — метим как INF по оценке
    big = (np.abs(P.astype(np.float64)) + np.abs(R.astype(np.float64))) >= 2.0 ** 63
    tag[big] = INF
    lo_out = (P - R) >> 1
    hi_out = (P + R) >> 1
    return lo_out, hi_out, tag


# ---------- самоприменение ----------

def self_check(deep: bool = False) -> list:
    """Самопроверка точного GEMM ВСТРЕЧНЫМИ дорогами (ось самоприменения).

    `[вычислено здесь]` (1) граница расщепления: `split_bound` сверяется с
    прямым целым неравенством D·(2^s − 1)² ≤ 2^b − 1, посчитанным здесь
    заново, а `piece_bits(D)` обязан давать МАКСИМУМ — граница держится при
    s и ломается при s + 1 (окно глубин 1, 2, 3, 100, 1024, 4096); (2)
    расщепление: Σ pieces[i]·2^(s·i), собранная целыми Python, равна
    исходному массиву, младшие куски в [0, 2^s), старший по модулю меньше
    2^(s−1); (3) итог на ядре против `exact_gemm_bigint` (целые Python,
    без ядра) на матрице с произведениями и внутри, и за контрактом 64:
    значения равны там, где метка RATIONAL, метка INF стоит ровно там, где
    |V| ≥ 2^63; (4) интервалы: [x − 1, x + 1]·W содержит точное x·W, нижняя
    не выше верхней; (5) договор: контракт шире ответа (80 бит) при 2^63 ≤ |V| < 2^79 — отказ
    `ExactRepresentationNotAvailable`, не заворот молча. При `deep` матрицы
    шире (8 × 16 на 16 × 4). Возвращает список находок; исключений не бросает.
    """
    out: list = []
    b = _F64_EXACT
    for depth in (1, 2, 3, 100, 1024, 4096):
        s = piece_bits(depth)
        for ss in (s, s + 1):
            direct = depth * (2 ** ss - 1) ** 2 <= 2 ** b - 1
            try:
                _, fits = split_bound(depth, ss)
            except Exception as e:                                # noqa: BLE001
                out.append(f"split_bound({depth}, {ss}) упал: {type(e).__name__}: {e}")
                continue
            if fits != direct:
                out.append(f"split_bound({depth}, {ss}) = {fits} ≠ прямому неравенству {direct}")
        if not (depth * (2 ** s - 1) ** 2 <= 2 ** b - 1):
            out.append(f"piece_bits({depth}) = {s} — граница не держится")
        if depth * (2 ** (s + 1) - 1) ** 2 <= 2 ** b - 1:
            out.append(f"piece_bits({depth}) = {s} — не максимум: держится и при {s + 1}")
    a = np.array([0, 1, -1, (1 << 40) + 12345, -(1 << 45) + 7, (1 << 62) - 1, -(1 << 62)], dtype=np.int64)
    for s in (7, 21):
        pieces = split(a, s)
        back = [sum(int(p[i]) * (1 << (s * k)) for k, p in enumerate(pieces)) for i in range(a.size)]
        if back != [int(v) for v in a]:
            out.append(f"split(s={s}) не собирается обратно: {back[:3]}… ≠ {a[:3]}…")
        for p in pieces[:-1]:
            if int(p.min()) < 0 or int(p.max()) >= (1 << s):
                out.append(f"split(s={s}): младший кусок вне [0, 2^s)")
        if int(np.abs(pieces[-1]).max()) >= (1 << (s - 1)):
            out.append(f"split(s={s}): старший кусок не меньше 2^(s−1) по модулю")
    rows, inner, cols = (8, 16, 4) if deep else (3, 4, 2)
    x = np.array([[(i * 7 + j * 3) % 11 - 5 for j in range(inner)] for i in range(rows)], dtype=np.int64)
    w = np.array([[(i * 5 + j) % 9 - 4 for j in range(cols)] for i in range(inner)], dtype=np.int64)
    x[0, 0] = 1 << 40
    w[0, 0] = 1 << 40                                            # произведение 2^80 — за контрактом 64
    try:
        v, tag = exact_gemm(x, w)
        big = exact_gemm_bigint(x, w)
        for i in range(rows):
            for j in range(cols):
                exact = int(big[i, j])
                over = abs(exact) >= (1 << 63)
                if (tag[i, j] == INF) != over:
                    out.append(f"метка в [{i},{j}]: {int(tag[i, j])} при |V| {'≥' if over else '<'} 2^63")
                elif not over and int(v[i, j]) != exact:
                    out.append(f"ядро в [{i},{j}] = {int(v[i, j])} ≠ целым Python {exact}")
    except Exception as e:                                        # noqa: BLE001
        out.append(f"exact_gemm/bigint упали на здоровом входе: {type(e).__name__}: {e}")
    try:
        lo, hi, tag = interval_gemm(x - 1, x + 1, w)
        big = exact_gemm_bigint(x, w)
        for i in range(rows):
            for j in range(cols):
                if tag[i, j] == INF:
                    continue
                if not (int(lo[i, j]) <= int(big[i, j]) <= int(hi[i, j])):
                    out.append(f"интервал в [{i},{j}] не содержит точного: {int(lo[i, j])}…{int(hi[i, j])} ∌ {int(big[i, j])}")
    except Exception as e:                                        # noqa: BLE001
        out.append(f"interval_gemm упал на здоровом входе: {type(e).__name__}: {e}")
    try:
        # |V| = 2^70: внутри контракта 80, но за шириной ответа int64 — отказ по договору
        exact_gemm(np.array([[1 << 35]], dtype=np.int64), np.array([[1 << 35]], dtype=np.int64), contract_bits=80)
        out.append("контракт шире ответа при 2^63 ≤ |V| < 2^79 не отказал — заворот молча")
    except ExactRepresentationNotAvailable:
        pass
    except Exception as e:                                        # noqa: BLE001
        out.append(f"контракт шире ответа отказал чужим исключением {type(e).__name__}")
    return out
