"""Спецметка нуля — who is who.

Слово автора 04.09 (чат): «у нас в терминологии нулей разных много; понять,
какой конкретный ноль имеется в виду, теперь не может никто, и ИИ тоже:
интерференция смысла полная. Требуется спецметка, отделяющая разные нули:
у них разные функции и нагрузки, обязательства и права».

Перепись нулей репозитория (замер 04.09 по коду и корпусу) и их права —
таблица ниже; каждому — спецметка вида 0{код} из КОНЕЧНОГО списка, и
функция mark_of(x), дающая метку любому нулевому объекту phicore. Права —
не декларация: они меряются пробами (proto/test_phicore_zeros.py):
ноль ∘ проба для проб 5, ω, ε под + и ×; строка прав = наблюдённый исход.
Два имени с одной строкой прав — один ноль под двумя именами (омоним),
одно имя с двумя строками — интерференция; и то и другое печатается.

| метка | где живёт | что это | права под + | права под × |
|---|---|---|---|---|
| 0{A1} | zero_state.ABSOLUTE; ФIOS E0 | мгновение определения (аксиома A1) | нейтрален (уступает) | аннигилирует: 0·x = 0, 0·ω = 0 |
| 0{void} | container class 0, void | nihil absolutum — пустота, предшествует числу | уступает | поглощает (итог ZERO) |
| 0{chaos} | container class 0, chaos | chaos ordinabile — та же пустота с провенансом «структурный» | уступает | поглощает (итог ZERO) |
| 0{Ω} | container class 20, center_zero(absolute); ФIOS U20 | познанный предел (Хайтин) | участвует как 0 | поглощает (итог UNCOMPUTABLE) |
| 0{ω−ω} | zero_state.STRUCTURAL_INFINITE; container center_zero(structural,'infinite') | неограниченная неопределённость | ПОГЛОЩАЕТ | поглощает, род хранит |
| 0{ε−ε} | zero_state.STRUCTURAL_INFINITESIMAL; container (structural,'infinitesimal') | ограниченная неопределённость (идеал бесконечно малых) | уступает | поглощает, род хранит |
| 0{∅} | zero_state.NOT_COMPUTED | «вычисление не производилось», 0∉S — состояние | ОТКАЗ | ОТКАЗ |
| 0{1} | rational.Q(0), PhiAlpha m=0, container from_int(0), int 0 | арифметический ноль класса 1 | нейтрален | аннигилирует; обратный — отказ |
| 0{ω−ω}= 0{ε−ε}= 0{a−a}= | hyper.Hyper.zero(origin) | арифметический ноль с записью пути (поле порядка) | нейтрален | аннигилирует; запись пути хранит |
| 0{g} | gridcore.Notch, sign == 0 | ноль СЕТКИ: засечка накрывает узел 0 при разрешении 2^−g — не доказательство нуля величины | как интервал | как интервал |
| 0{форма} | affine.Affine, c0 = 0 и коэффициентов нет | точное сокращение формы (x − x) — ноль с сертификатом ширины 0 | нейтрален | аннигилирует |
| 0{метка} | bench_classnet.ZERO (тег 5 бит) | значение сети ровно 0 — класс на значении | по правилу склейки меток | по правилу склейки меток |
| 0{IEEE±0} | класс 23, float ±0.0 | знаковый ноль IEEE 754: 1/−0 = −inf | как у float | как у float |
| 0{ОП+0}, 0{ОП−0}, 0{ОП0} | op_axis.Seam — шов дискретной оси-тора | пересечение шва вперёд / назад / шов без направления (выведено из первых принципов, 04.09) | уступает; +0 + (−0) = 0{ОП0} | знак произведения; x/(±0) = ±ω; 1/0{ОП0} — отказ |
| 0{множество} | GLOSSARY §1 «структурный ноль» (Энигма) | множество, которое не возникает никогда; принадлежность разрешима за O(1) | не арифметика | не арифметика |

Омонимы, замеренные по корпусу: «структурный ноль» = 0{chaos} (класс 0),
0{ω−ω} и 0{ε−ε} (класс 20), 0{множество} (комбинаторика) — четыре
разных объекта под одним словом; «абсолютный» = 0{A1}, 0{void}, 0{Ω} —
три. Спецметка их разводит; слово без метки в записи не употребляется.

Авторство: роды нуля, ноль-состояние, пустота/хаос — AIAΣA (автор
проекта; v6, OP-ш Theory); IEEE 754-2019 (±0); Мур 1966 (нуль сетки).
NOT_MEASURED — ЗАКРЫТО 12.09 замером (строка остаётся записью о том, что
было). Замер `carriers_of_marks()` двумя независимыми дорогами: прогон по
носителям дома и разбор исходника самого `mark_of`. Обе дают одно и то же —
из 19 меток носителем дома достижимы 17, недостижимы 2, и это РАЗНЫЕ
утверждения, каждое со своим доказательством:

  * `0{множество}` — ВЕЛИЧИНЫ НЕТ: это множество из GLOSSARY §1 (Энигма),
    предмет комбинаторики, а не носитель `phicore`; принадлежность
    разрешима за O(1), арифметики у него нет вовсе. Мерить нечего — не
    «не умею посчитать», а «объекта такого рода дом не несёт»;
  * `0{метка}` — НОСИТЕЛЯ НЕТ, хотя строка таблицы есть: `bench_classnet.ZERO`
    (`proto/bench_classnet.py`) есть КОД класса — пятибитный тег со значением
    int 0, а не объект-ноль; сам код дом метит как `0{1}` (замер:
    `mark_of(0)` даёт ARITH). Это класс на значении сети, и объекта под
    строкой нет.

Знаковые нули сетки: `Notch.signed()` даёт ровно ТРИ значения знака
(−1, 0, +1 — замер по засечкам), ноль-метку `0{g}` получает ровно одно из
трёх (знак 0), а сама тройка `signed()` ноль-меткой не метится: это
значения линейки, не нули.
"""
from __future__ import annotations

import math
from enum import Enum


class ZeroMark(Enum):
    A1 = "0{A1}"
    VOID = "0{void}"
    CHAOS = "0{chaos}"
    OMEGA_LIMIT = "0{Ω}"
    STRUCT_INF = "0{ω−ω}"
    STRUCT_EPS = "0{ε−ε}"
    NOT_COMPUTED = "0{∅}"
    ARITH = "0{1}"
    HYPER_INF = "0{ω−ω}="
    HYPER_EPS = "0{ε−ε}="
    HYPER_FIN = "0{a−a}="
    GRID = "0{g}"
    FORM = "0{форма}"
    TAG = "0{метка}"
    IEEE = "0{IEEE±0}"
    OP_PLUS = "0{ОП+0}"
    OP_MINUS = "0{ОП−0}"
    OP_NEUTRAL = "0{ОП0}"
    SET = "0{множество}"


# слова корпуса → метки (омонимы названы явно)
HOMONYMS = {
    "структурный ноль": (ZeroMark.CHAOS, ZeroMark.STRUCT_INF, ZeroMark.STRUCT_EPS, ZeroMark.SET),
    "абсолютный ноль": (ZeroMark.A1, ZeroMark.VOID, ZeroMark.OMEGA_LIMIT),
    "ноль": tuple(ZeroMark),
}


def carriers_of_marks() -> dict:
    """Метка → пути построения её носителя в доме (замер ПРОГОНОМ).

    Строится каждый носитель, какой дом несёт, и метится `mark_of`. Метка с
    пустым списком есть заявление о ПРЕДМЕТЕ: носителя под ней в `phicore`
    нет, и причина названа в шапке модуля поимённо.

    Среда, у которой нет numpy, не делает метку недостижимой: путь формы
    (`affine.Affine`) в этом случае помечается словом «среда без numpy», а
    не молчанием — ноль инструмента не есть ноль предмета."""
    from .container import PhiNumberV6 as P
    from .gridcore import Notch
    from .hyper import Hyper
    from .op_axis import Seam
    from .rational import Q
    from .zero_state import Zero
    w, e = Hyper.omega(), Hyper.epsilon()
    пути = {
        "Zero.absolute()": Zero.absolute(),
        "Zero.structural_infinite()": Zero.structural_infinite(),
        "Zero.structural_infinitesimal()": Zero.structural_infinitesimal(),
        "Zero.not_computed()": Zero.not_computed(),
        "PhiNumberV6.void()": P.void(),
        "PhiNumberV6.chaos()": P.chaos(),
        "PhiNumberV6.center_zero(structural=False)": P.center_zero(structural=False),
        "PhiNumberV6.from_int(0)": P.from_int(0),
        "Q(0)": Q(0),
        "Hyper.omega() − Hyper.omega()": w - w,
        "Hyper.epsilon() − Hyper.epsilon()": e - e,
        "Hyper.const(5) − Hyper.const(5)": Hyper.const(5) - Hyper.const(5),
        "Notch(20, −1, 1)": Notch(20, -1, 1),
        "Notch(20, 0, 0)": Notch(20, 0, 0),
        "float −0.0 (класс 23)": -0.0,
        "Seam.plus()": Seam.plus(),
        "Seam.minus()": Seam.minus(),
        "Seam.neutral()": Seam.neutral(),
    }
    out = {m: [] for m in ZeroMark}
    try:
        from .affine import Affine
        a = Affine.from_notch(Notch(20, 5, 7))
        пути["Affine.from_notch(…) − сама себя"] = a - a
    except ImportError:                       # среда без numpy — названо, не скрыто
        out[ZeroMark.FORM].append("среда без numpy: phicore/affine.py не ввозится")
    for путь, x in пути.items():
        m = mark_of(x)
        if m is not None:
            out[m].append(путь)
    return out


def unreachable_marks() -> tuple:
    """Метки без носителя в доме — кортеж в порядке объявления.

    Утверждение о предмете, не об усилии: причина каждой названа в шапке
    модуля (одна — величины нет, другая — носителя нет, есть код класса)."""
    носители = carriers_of_marks()
    return tuple(m for m in ZeroMark if not носители[m])


def mark_of(x) -> ZeroMark | None:
    """Спецметка нулевого объекта; None — объект не ноль (или не phicore)."""
    from .zero_state import Zero, ZeroKind
    if isinstance(x, Zero):
        return {ZeroKind.ABSOLUTE: ZeroMark.A1, ZeroKind.STRUCTURAL_INFINITE: ZeroMark.STRUCT_INF,
                ZeroKind.STRUCTURAL_INFINITESIMAL: ZeroMark.STRUCT_EPS, ZeroKind.NOT_COMPUTED: ZeroMark.NOT_COMPUTED}[x.kind]
    from .container import NumberClass, PhiNumberV6
    if isinstance(x, PhiNumberV6):
        nc = x.number_class
        if nc == NumberClass.ZERO:
            return ZeroMark.CHAOS if x.zero_provenance else ZeroMark.VOID
        if nc == NumberClass.UNCOMPUTABLE:
            if not x.zero_provenance:
                return ZeroMark.OMEGA_LIMIT
            return ZeroMark.STRUCT_EPS if x.structural_kind == "infinitesimal" else ZeroMark.STRUCT_INF
        if x.m1.v == 0 and nc not in (NumberClass.INFINITY, NumberClass.INFINITESIMAL, NumberClass.IGNORANCE):
            return ZeroMark.ARITH
        if nc == NumberClass.IGNORANCE and x.base_params.get("raw_float") == 0.0:
            return ZeroMark.IEEE
        return None
    from .hyper import Hyper
    if isinstance(x, Hyper):
        if x.terms:
            return None
        return {"ω−ω": ZeroMark.HYPER_INF, "ε−ε": ZeroMark.HYPER_EPS, "a−a": ZeroMark.HYPER_FIN}.get(x.origin, ZeroMark.ARITH)
    from .gridcore import Notch
    if isinstance(x, Notch):
        return ZeroMark.GRID if x.sign == 0 else None
    from .affine import Affine
    if isinstance(x, Affine):
        return ZeroMark.FORM if (x.c0 == 0 and not x.coef) else (ZeroMark.GRID if x.notch().sign == 0 else None)
    from .op_axis import Seam
    if isinstance(x, Seam):
        return {1: ZeroMark.OP_PLUS, -1: ZeroMark.OP_MINUS, 0: ZeroMark.OP_NEUTRAL}[x.sign]
    from .rational import Q
    if isinstance(x, Q):
        return ZeroMark.ARITH if x.n == 0 else None
    try:
        from .alpha import PhiAlpha
        if isinstance(x, PhiAlpha):
            return ZeroMark.ARITH if x.m.n == 0 else None
    except ImportError:
        pass
    if isinstance(x, bool):
        return None
    if isinstance(x, int):
        return ZeroMark.ARITH if x == 0 else None
    if isinstance(x, float):
        return ZeroMark.IEEE if x == 0.0 and not math.isnan(x) else None
    return None


def record(x) -> str:
    """Запись объекта со спецметкой: «0{код}» для нуля, иначе — его собственная запись."""
    m = mark_of(x)
    if m is not None:
        return m.value
    rec = getattr(x, "record", None)
    return rec() if callable(rec) else repr(x)


def rights_row(zero, probes: dict, ops: dict) -> tuple:
    """Наблюдённая строка прав: для каждой пробы и операции — исход:
    'уступает' (итог равен пробе), 'поглощает' (итог — ноль той же метки),
    'ноль:<метка>' (ноль другой метки), 'отказ' (исключение), 'иное'."""
    row = []
    my = mark_of(zero)
    for pname, p in probes.items():
        for oname, op in ops.items():
            for order in ("z∘p", "p∘z"):
                try:
                    r = op(zero, p) if order == "z∘p" else op(p, zero)
                except Exception as e:  # noqa: BLE001 — отказ есть исход
                    row.append((pname, oname, order, "отказ:" + type(e).__name__)); continue
                mr = mark_of(r)
                if mr is not None and mr is my:
                    out = "поглощает"
                elif mr is not None:
                    out = "ноль:" + mr.value
                elif _same(r, p):
                    out = "уступает"
                else:
                    out = "иное"
                row.append((pname, oname, order, out))
    return tuple(row)


def _same(a, b) -> bool:
    try:
        if a is b or a == b:
            return True
    except Exception:  # noqa: BLE001
        pass
    ra, rb = getattr(a, "record", None), getattr(b, "record", None)
    if callable(ra) and callable(rb):
        return ra() == rb()
    return repr(a) == repr(b)


# ---------- самопроверка переписи нулей (ось самоприменения, 09.10) ----------

def self_check(deep: bool = False) -> list:
    """Самопроверка спецметок нуля встречными дорогами (ось «AIAΣA спрашивает
    AIAΣA», слово автора 04.09; очередь строка 18; 09.10). Список находок,
    пусто — здорова; исключений не бросает.

    `[вычислено здесь]` (1) список меток конечен и значения различны;
    омонимы корпуса — подмножества списка, «ноль» покрывает все; (2) две
    НЕЗАВИСИМЫЕ дороги переписи: прогон носителей дома (`carriers_of_marks`)
    против разбора исходника самого `mark_of` (ast: какие метки он называет)
    — множества достижимых меток обязаны совпасть, недостижимы ровно
    0{метка} и 0{множество}, причины названы в шапке; (3) каждый носитель окна
    метится ожидаемой меткой, а не-нули (единица, засечка без нуля, половина,
    bool) не метятся вовсе; `record` нуля есть значение метки; (4) права
    замером проб: 0{A1} и 0{1} уступают под сложением и поглощают под
    умножением, 0{∅} отказывает под обоими, 0{ω−ω} поглощает под умножением
    [в корпусе: таблица прав модуля]. При deep — права ещё под пробами ω и ε."""
    import ast as _ast
    import operator as _op
    import pathlib as _pl
    out = []
    try:
        values = [m.value for m in ZeroMark]
        if len(set(values)) != len(values):
            out.append("значения меток не различны")
        for word, marks in HOMONYMS.items():
            if not set(marks) <= set(ZeroMark):
                out.append(f"омоним «{word}» называет метку вне списка")
        if set(HOMONYMS["ноль"]) != set(ZeroMark):
            out.append("омоним «ноль» покрывает не все метки")
    except Exception as e:                                      # noqa: BLE001
        out.append(f"список меток: {type(e).__name__}: {e}")
    try:
        carriers = carriers_of_marks()
        reachable = {m for m, paths in carriers.items() if paths}
        src = _ast.parse(_pl.Path(__file__).read_text(encoding="utf-8"))
        fn = next(n for n in _ast.walk(src) if isinstance(n, _ast.FunctionDef) and n.name == "mark_of")
        named = {n.attr for n in _ast.walk(fn)
                 if isinstance(n, _ast.Attribute) and isinstance(n.value, _ast.Name) and n.value.id == "ZeroMark"}
        if {m.name for m in reachable} != named:
            out.append(f"дороги переписи разошлись: прогон {sorted(m.name for m in reachable)} против исходника {sorted(named)}")
        if set(unreachable_marks()) != {ZeroMark.TAG, ZeroMark.SET}:
            out.append(f"недостижимые метки не те: {[m.value for m in unreachable_marks()]}")
        if "ВЕЛИЧИНЫ НЕТ" not in __doc__ or "НОСИТЕЛЯ НЕТ" not in __doc__:
            out.append("причины недостижимости не названы в шапке модуля")
    except Exception as e:                                      # noqa: BLE001
        out.append(f"перепись двумя дорогами: {type(e).__name__}: {e}")
    try:
        from .container import PhiNumberV6 as P
        from .gridcore import Notch
        from .hyper import Hyper
        from .op_axis import Seam
        from .rational import Q
        from .zero_state import Zero
        from .affine import Affine
        a = Affine.from_notch(Notch(20, 5, 7))
        expect = [(Zero.absolute(), ZeroMark.A1), (Zero.not_computed(), ZeroMark.NOT_COMPUTED), (P.void(), ZeroMark.VOID),
                  (P.chaos(), ZeroMark.CHAOS), (Q(0), ZeroMark.ARITH), (0, ZeroMark.ARITH), (Notch(20, -1, 1), ZeroMark.GRID),
                  (a - a, ZeroMark.FORM), (Seam.plus(), ZeroMark.OP_PLUS), (Seam.neutral(), ZeroMark.OP_NEUTRAL),
                  (Hyper.omega() - Hyper.omega(), ZeroMark.HYPER_INF), (-0.0, ZeroMark.IEEE)]
        for x, m in expect:
            if mark_of(x) is not m:
                out.append(f"метка {record(x)!r} ≠ {m.value} у носителя {x!r}")
            elif record(x) != m.value:
                out.append(f"record нуля {m.value} не есть значение метки")
        for x in (Q(1), Notch(20, 1, 2), 1, 0.5, True, False, Seam):
            if mark_of(x) is not None:
                out.append(f"не-ноль {x!r} получил метку {mark_of(x)}")
        probes = {"5": Q(5)}
        if deep:
            probes["ω"], probes["ε"] = Hyper.omega(), Hyper.epsilon()
        ops = {"+": _op.add, "×": _op.mul}
        def outcomes(z, op):
            return {r[3] for r in rights_row(z, {"5": Q(5)}, {op: ops[op]})}
        for z, name in ((Zero.absolute(), "0{A1}"), (Q(0), "0{1}")):
            if outcomes(z, "+") != {"уступает"} or outcomes(z, "×") != {"поглощает"}:
                out.append(f"права {name}: ожидались «уступает» под + и «поглощает» под ×, замерено {outcomes(z, '+')}, {outcomes(z, '×')}")
        if any(not r[3].startswith("отказ") for r in rights_row(Zero.not_computed(), {"5": Q(5)}, ops)):
            out.append("права 0{∅}: не отказывает под + или ×")
        if outcomes(Hyper.omega() - Hyper.omega(), "×") != {"поглощает"}:
            out.append("права 0{ω−ω}: не поглощает под ×")
        if deep:
            row = rights_row(Zero.absolute(), probes, ops)
            if len(row) != 12:
                out.append(f"права 0{{A1}} под тремя пробами: строк {len(row)} вместо 12")
    except Exception as e:                                      # noqa: BLE001
        out.append(f"носители и права: {type(e).__name__}: {e}")
    return out
