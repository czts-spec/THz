"""Комплиментарная таблица (ДНК): обязательные соответствия классов.

Постановка автора (консенсус корпуса, п.2.7): рядом с правилами живут
ОБЯЗАТЕЛЬНЫЕ СООТВЕТСТВИЯ — двойственности; их проверка есть репарация
без внешней задачи. До этого модуля двойственности жили только во
внешних тестах; здесь они — данные самой библиотеки: строка = (шаг §3,
классы, имя соответствия, проверка), прогонщик `repair()` исполняет все
и называет числа. Библиотека может проверить саму себя изнутри.

Каждая проверка ТОЧНАЯ (тождество объектов/Q) либо
СЕРТИФИКАТНАЯ (вложение на сетке); род назван в строке. Несведённое
соответствие — не тихий пропуск: строка с упавшей проверкой попадает в
отчёт репарации поимённо.

Самопроверка модуля `self_check(deep)` (ось самоприменения, 09.10) судит
таблицу как данные и прогонщик `repair()` как счётчик — встречными
дорогами, а не самой таблицей: шаги по карте §3 `STEPS_3`, классы по
`Фios.LIVE`, счёт против независимого прогона строк и против
синтетической таблицы с заранее известным ответом. Математику строк
судят их приёмки.
"""
from __future__ import annotations

import re

from .rational import Q

from .algebraic import AlgebraicNumber
from .arctrig import arctrig
from .bessel import bessel_root
from .de_solution import DESolutionValue
from .elliptic import complementary, elliptic_k, legendre_enclosure
from .enclosure import exp_enclosure, pi_enclosure
from .exp_log import LogValue
from .gamma_value import GammaValue
from .hyperbolic import pythagorean_check
from .intpoly import IntPoly
from .iteration import iteration_limit
from .lambert import lambert_w
from .pi_value import PI, PiMonomial
from .product import InfiniteProductValue
from .series import PowerSeriesValue
from .unity_root import UnityRoot
from .zeta import zeta


def _d3_annihilator():
    s = AlgebraicNumber.from_root(2, 2) + AlgebraicNumber.from_root(2, 3)
    return s.minimal_annihilator() == IntPoly([1, 0, -10, 0, 1])


def _d3_minimal():
    """3↔3 с факторизацией над ℤ (04.09): узел × строит аннулятор степени 6,
    минимальный многочлен — его неприводимый множитель с корнем на отрезке."""
    s = AlgebraicNumber.from_root(2, 2) + AlgebraicNumber.from_root(2, 3)
    sq = s * s
    return (sq.annihilator().degree == 6
            and sq.minimal_polynomial() == IntPoly([1, -10, 1]))


def _d4_3_cyclotomic():
    """4↔3 (04.09): ζ₅ + ζ₅⁻¹ = 2cos(2π/5) = φ⁻¹ — класс 3 с минимальным
    многочленом x²+x−1; re(ζ₅) = cos 72° той же дорогой."""
    z = UnityRoot(5, 1)
    s = z + UnityRoot(5, 4)
    return (isinstance(s, AlgebraicNumber) and s.p == IntPoly([-1, 1, 1])
            and z.re() * 2 == s and z - UnityRoot(10, 2) == Q(0))


def _d3_self_check():
    """Самоприменение линейки 3 (04.09): self_check(deep=True) пуст на
    узлах +/× семейства; факторизация сверяет произведение сама."""
    r2, r3, c2 = (AlgebraicNumber.from_root(2, 2), AlgebraicNumber.from_root(2, 3),
                  AlgebraicNumber.from_root(3, 2))
    fam = (r2 + r3, (r2 + r3) * (r2 + r3), c2 * r2, (r2 + Q(1)) * r2)
    return all(x.self_check(deep=True) == [] for x in fam)


def _d4_3_pair():
    """4↔3 через пару (re, im) (04.09): Σ_{k<5} ζ₅^k = 0 точно;
    |ζ₅ + ζ₅²|² = 2 + 2cos(2π/5) = φ² — минимальный многочлен x²−3x+1."""
    tot = Q(0)
    for k in range(5):
        tot = tot + UnityRoot.of(5, k)
    z = UnityRoot(5, 1) + UnityRoot(5, 2)
    return tot == Q(0) and z.re == Q(-1, 2) and z.magnitude_squared().minimal_polynomial() == IntPoly([1, -3, 1])


def _d3_mirror():
    """Зазеркалье (04.09): палиндромный аннулятор → θ + 1/θ половинной степени:
    (√2+√3) + 1/(√2+√3) = 2√3 точно."""
    s = AlgebraicNumber.from_root(2, 2) + AlgebraicNumber.from_root(2, 3)
    return s.p.is_palindromic() and s.mirror_sum() == AlgebraicNumber.from_root(2, 3) * 2


def _d3_sign():
    """Знак двумя дорогами (04.09, Гильберт 17 в дискретной форме): отрезок
    и Штурм согласны на √2−√3 (−), 2cos(2π/7) (+), 2cos(4π/7) (−), ∛2 (+)."""
    from .cyclotomic import two_cos
    r2, r3 = AlgebraicNumber.from_root(2, 2), AlgebraicNumber.from_root(2, 3)
    fam = ((r2 - r3, -1), (two_cos(1, 7), 1), (two_cos(2, 7), -1),
           (AlgebraicNumber(IntPoly([-2, 0, 0, 1]), Q(-3), Q(3)), 1))
    return all(v.sign_certificate() == {"отрезок": w, "Штурм": w, "знак": w}
               for v, w in fam)


def _d4_odometer():
    w = UnityRoot(12, 1)
    return w ** 12 == Q(1) and w ** 6 == Q(-1)


def _d5_12_pythagorean():
    return pythagorean_check(Q(3, 7)) == Q(1)


def _d6_log_canon():
    return LogValue.make(1, 8) == LogValue.make(3, 2)


def _d7_series_exp():
    ряд = PowerSeriesValue(Q(1), IntPoly([1]), IntPoly([1, 1]),
                           Q(1))
    return (ряд.partial_sum(5) == Q(163, 60)
            and ряд.enclosure(40).intersect(exp_enclosure(1, 30)) is not None)


def _d8_wallis_log():
    p = InfiniteProductValue(IntPoly([0, 0, 4]), IntPoly([-1, 0, 4]), 1)
    return p.log_partial(6) == LogValue.make(1, p.partial_product(6))


def _d8_7_series():
    """8↔7 (§3 шаг 3): ln f₁ множителя Валлиса (класс 6 от 8) и ряд
    −ln(1−x) = x·Σ xᵏ/(k+1) (класс 7) при x = 1/4 — сертификаты пересекаются."""
    wallis = InfiniteProductValue(IntPoly([0, 0, 4]), IntPoly([-1, 0, 4]), 1)
    lg = LogValue.make(Q(1), wallis.factor(1)).enclosure(60)
    sr = PowerSeriesValue(Q(1), IntPoly([1, 1]), IntPoly([2, 1]), Q(1, 4)).enclosure(40)
    lo, hi = sr.lo * Q(1, 4), sr.hi * Q(1, 4)
    return lg.lo <= hi and lo <= lg.hi


def _d9_integral():
    from .integral import integral
    return (integral(IntPoly([0, 0, 1]), IntPoly([1]), Q(0),
                     Q(1)) == Q(1, 3)
            and integral(IntPoly([1]), IntPoly([1, 0, 1]), Q(0),
                         Q(1)) == PiMonomial(Q(1, 4), 1))


def _d10_machin():
    a = arctrig("arctan", Q(1, 5)).enclosure(40)
    b = arctrig("arctan", Q(1, 239)).enclosure(40)
    return (a * 16 - b * 4).intersect(pi_enclosure(30)) is not None


def _d11_complement():
    return (arctrig("arctan", Q(2, 3))
            + arctrig("arctan", Q(3, 2))
            == PiMonomial(Q(1, 2), 1))


def _d13_functional_eq():
    lhs = zeta(Q(-1))
    from .gamma_value import sin_pi
    rhs = Q(1, 2)
    rhs = rhs * PiMonomial(Q(1), -2)
    rhs = rhs * sin_pi(Q(-1, 2))
    rhs = rhs * GammaValue.of(Q(2))
    rhs = rhs * zeta(Q(2))
    return lhs == Q(-1, 12) and rhs == Q(-1, 12)


def _d14_reflection():
    return (GammaValue.of(Q(1, 2))
            * GammaValue.of(Q(1, 2)) == PI)


def _d15_elliptic():
    ok_inv = all(complementary(complementary(Q(n, 7)))
                 == Q(n, 7) for n in range(1, 7))
    leg = legendre_enclosure(Q(1, 2), 48)
    return ok_inv and (elliptic_k(Q(0)) == PiMonomial(Q(1, 2), 1)
                       and leg.intersect(pi_enclosure(30) / 2) is not None)


def _d16_bessel_half():
    return (bessel_root(Q(1, 2), 1) == PI
            and bessel_root(Q(1, 2), 3) == PiMonomial(Q(3), 1))


def _d17_lambert():
    w = lambert_w(Q(1)).enclosure(50)
    e_lo = exp_enclosure(w.lo, 30)
    e_hi = exp_enclosure(w.hi, 30)
    prod_lo = min(w.lo * e_lo.lo, w.lo * e_hi.lo)
    prod_hi = max(w.hi * e_hi.hi, w.hi * e_lo.hi)
    return prod_lo <= 1 <= prod_hi


def _d18_substitution():
    y = DESolutionValue([Q(1)], [], Q(1))
    return y.check_substitution(15) and y.coeffs(6)[5] == Q(1, 120)


def _d19_fixed_point():
    lim = iteration_limit(IntPoly([2, 0, 1]), IntPoly([0, 2]), Q(2))
    return lim.equals(AlgebraicNumber.from_root(2, 2))


def _d19_17_lambert():
    """19↔17 (§3 шаг 11): неподвижная точка x = e^{−x} интервальной итерацией
    (f(I) ⊂ I — Брауэр на отрезке) накрывает сертификат W(1) класса 17."""
    from .enclosure import Interval
    from .exp_log import ExpValue

    def f(iv):
        lo = ExpValue.make(Q(1), -iv.hi).enclosure(40)
        hi = ExpValue.make(Q(1), -iv.lo).enclosure(40)
        return Interval(lo.lo, hi.hi)

    I = Interval(Q(1, 2), Q(13, 20))
    J = f(I)
    if not (I.lo <= J.lo and J.hi <= I.hi):
        return False
    for _ in range(30):
        J = f(I)
        I = Interval(max(I.lo, J.lo), min(I.hi, J.hi))
    w = lambert_w(Q(1)).enclosure(40)
    return w.lo <= I.hi and I.lo <= w.hi


#: Таблица соответствий — данные: (шаг §3, «классы», имя, род, проверка).
CORRESPONDENCES = (
    (1, "3↔3", "√2+√3 ↔ x⁴−10x²+1", "точное", _d3_annihilator),
    (1, "3↔3", "(√2+√3)² ↔ x²−10x+1: минимальный из аннулятора степени 6 (факторизация над ℤ)", "точное", _d3_minimal),
    (1, "3", "самопроверка линейки 3: self_check(deep) пуст на √2+√3, (√2+√3)², ∛2·√2, (1+√2)·√2", "точное", _d3_self_check),
    (1, "3", "зазеркалье θ ↔ 1/θ: палиндром x⁴−10x²+1 → (√2+√3) + 1/(√2+√3) = 2√3 половинной степенью", "точное", _d3_mirror),
    (1, "3", "знак двумя дорогами: отрезок ∥ Штурм на √2−√3, 2cos(2π/7), 2cos(4π/7), ∛2 (Гильберт 17 в дискретной форме)", "точное", _d3_sign),
    (0, "4", "одометр: ζ₁₂¹²=1, ζ₁₂⁶=−1", "точное", _d4_odometer),
    (0, "4↔3", "ζ₅ + ζ₅⁻¹ = 2cos(2π/5) = φ⁻¹ ↔ x²+x−1 (круговая величина на линейке 3)", "точное", _d4_3_cyclotomic),
    (0, "4↔3", "пара (re, im): Σ_{k<5} ζ₅^k = 0; |ζ₅ + ζ₅²|² = φ² ↔ x²−3x+1", "точное", _d4_3_pair),
    (0, "12↔5", "cosh²−sinh²=1 сведено классом 5", "точное",
     _d5_12_pythagorean),
    (0, "6", "канон ln 8 = 3·ln 2", "точное", _d6_log_canon),
    (2, "7↔5", "частные суммы ↔ e (S₅=163/60; сертификаты согласны)",
     "смешанное", _d7_series_exp),
    (3, "8↔6", "Σln fₙ сворачивается в ln Π (Валлис N=6)", "точное",
     _d8_wallis_log),
    (3, "8↔7", "ln f₁ Валлиса (класс 6 от 8) ∩ ряд −ln(1−x) класса 7", "сертификат", _d8_7_series),
    (9, "9↔1,10", "∫x²=1/3; ∫dx/(1+x²)=π/4", "точное", _d9_integral),
    (0, "11↔10", "Мэчин: 16·arctan⅕−4·arctan(1/239) ∋ π", "сертификат",
     _d10_machin),
    (0, "11", "дополнение arctan x + arctan 1/x = π/2", "точное",
     _d11_complement),
    (4, "13↔14,10", "функциональное уравнение ζ при s=−1 → −1/12",
     "точное", _d13_functional_eq),
    (5, "14↔10", "Γ(1/2)² = π (отражение Эйлера)", "точное", _d14_reflection),
    (8, "15", "инволюция дополнения; K(0)=π/2; Лежандр ∋ π/2",
     "смешанное", _d15_elliptic),
    (7, "16↔10", "корни J_{1/2} = kπ", "точное", _d16_bessel_half),
    (6, "17↔5", "W·e^W ∋ 1 (сертификат)", "сертификат", _d17_lambert),
    (10, "18", "проверка ДУ подстановкой; a₅=1/120", "точное",
     _d18_substitution),
    (11, "19↔3", "предел Ньютона x²−2 ≡ √2", "точное", _d19_fixed_point),
    (11, "19↔17", "неподвижная точка e^{−x} (Брауэр на отрезке) ∋ W(1)", "сертификат", _d19_17_lambert),
)

#: Строки, добавленные 04.09 по замеру покрытия пар «класс ↔ двойственный»
#: таблицы 00-ЛИНИЯ §3 (proto/bench_dna_coverage.py): до них 9 из 11 пар.
ADDED_04_09 = frozenset({_d3_minimal, _d8_7_series, _d19_17_lambert})

#: Карта §3 как данные (таблица 00-ЛИНИЯ §3; приёмка `proto/test_duality_grid.py`:
#: «11 строк в порядке §3: 3, 7, 8, 13, 14, 17, 16, 15, 9, 18, 19»): шаг k
#: строки ведёт класс STEPS_3[k − 1]; шаг 0 — строка вне сетки §3 (одометр
#: класса 4, канон логарифма 6, дополнение 11, Пифагор 12 ↔ 5).
STEPS_3 = (3, 7, 8, 13, 14, 17, 16, 15, 9, 18, 19)

#: Роды проверки, один из которых строка обязана назвать: точное тождество
#: объектов/Q, сертификат (вложение на сетке), смешанное (обе части в строке).
KINDS = ("точное", "сертификат", "смешанное")


def repair(rows=None) -> dict:
    """Прогнать все обязательные соответствия; репарация без внешней
    задачи. Возвращает числа и поимённый список упавших/отказавших.
    `rows` — таблица для прогона; по умолчанию CORRESPONDENCES (имя
    берётся из модуля при вызове, так подлог таблицы в приёмке виден)."""
    rows = CORRESPONDENCES if rows is None else rows
    held, broken = 0, []
    for step, classes, name, kind, chk in rows:
        try:
            ok = bool(chk())
        except Exception as e:  # noqa: BLE001 — репарация ловит любой сбой
            ok = False
            name = f"{name} [{type(e).__name__}: {e}]"
        if ok:
            held += 1
        else:
            broken.append((classes, name, kind))
    return {"всего": len(rows), "держится": held,
            "нарушено": broken}


def _direct(rows) -> list:
    """Прямой прогон строк мимо repair(): истинность каждой проверки
    независимым вызовом, сбой считается ложью. Встречная дорога к счёту."""
    out = []
    for row in rows:
        try:
            out.append(bool(row[4]()))
        except Exception:  # noqa: BLE001 — та же свобода, что у репарации
            out.append(False)
    return out


def _base_name(name: str, row_names) -> str:
    """Имя строки из отчёта repair без хвоста «[Исключение: …]»."""
    for rn in row_names:
        if name == rn or name.startswith(rn + " ["):
            return rn
    return name


def _synthetic_failure():
    raise ValueError("нарочный сбой синтетической строки")


def self_check(deep: bool = False) -> list:
    """Самопроверка комплиментарной таблицы ВСТРЕЧНЫМИ дорогами (ось самоприменения).

    `[вычислено здесь]` (1) таблица как данные: у каждой строки шаг §3
    существует (0 — вне сетки, 1…11 — по карте STEPS_3) и согласован с
    ведущим классом строки; все названные классы живые по `Фios.LIVE`;
    имя соответствия уникально; род проверки назван (KINDS); проверка
    вызываема; реестр ADDED_04_09 указывает на строки, которые в таблице
    есть. (2) `repair()` против прямого прогона строк: «всего» равно
    числу строк, «держится» — числу строк с истинной проверкой при
    независимом вызове каждой, «нарушено» поимённо совпадает с упавшими в
    прямом прогоне, сумма сходится; и `repair(rows)` на синтетической
    таблице из трёх строк (истина, ложь, сбой ValueError) обязана назвать
    всего 3, держится 1, нарушено 2 с именем исключения у сбойной строки.
    (3) двойственность самой таблицы: у каждого шага §3 есть хотя бы одна
    строка. При `deep` — детерминизм и независимость от порядка: второй
    прямой прогон и второй `repair()` совпадают с первыми, обратный
    порядок строк даёт те же числа и тот же состав упавших.
    Граница: содержимое проверок — математику тождеств — судят их приёмки
    (`proto/test_phicore_dna.py`, `proto/test_dna_coverage.py`), а здесь
    судятся таблица как данные и прогонщик как счётчик; ложная проверка в
    таблице — строка отчёта repair, здесь она лишь сверяется по имени.
    Возвращает список находок; исключений не бросает."""
    out: list = []
    rows = tuple(CORRESPONDENCES)
    try:
        from . import Фios as _fios  # лениво: Фios ввозит носители всех классов
        live = {k[1] for k in _fios.LIVE}
    except Exception as e:  # noqa: BLE001
        live = set()
        out.append(f"Фios.LIVE недоступен: {type(e).__name__}: {e}")
    # (1) таблица как данные
    seen: dict = {}
    for row in rows:
        if not isinstance(row, tuple) or len(row) != 5:
            out.append(f"строка {row!r} не есть пятёрка (шаг, классы, имя, род, проверка)")
            continue
        step, classes, name, kind, chk = row
        cls = [int(n) for n in re.findall(r"\d+", str(classes))]
        if not isinstance(step, int) or not 0 <= step <= len(STEPS_3):
            out.append(f"строка «{name}»: шага §3 {step!r} нет (шаги 0…{len(STEPS_3)})")
        elif step and cls and cls[0] != STEPS_3[step - 1]:
            out.append(f"строка «{name}»: стоит на шаге {step} (класс {STEPS_3[step - 1]}), а ведёт класс {cls[0]}")
        elif not step and cls and cls[0] in STEPS_3:
            out.append(f"строка «{name}»: ведёт класс {cls[0]} шага {STEPS_3.index(cls[0]) + 1} §3, а стоит вне сетки (0)")
        if not cls:
            out.append(f"строка «{name}»: классы в записи {classes!r} отсутствуют")
        for c in cls:
            if live and c not in live:
                out.append(f"строка «{name}»: класс {c} не живой по Фios.LIVE")
        if kind not in KINDS:
            out.append(f"строка «{name}»: род проверки {kind!r} вне списка {KINDS}")
        if not callable(chk):
            out.append(f"строка «{name}»: проверка {chk!r} не вызываема")
        seen[name] = seen.get(name, 0) + 1
    for name, n in seen.items():
        if n > 1:
            out.append(f"имя соответствия «{name}» повторено {n} раз")
    for f in ADDED_04_09:
        if not any(isinstance(r, tuple) and len(r) == 5 and r[4] is f for r in rows):
            out.append(f"реестр строк 04.09 называет проверку {getattr(f, '__name__', f)!r}, которой в таблице нет")
    # (2) repair() против прямого прогона и против синтетической таблицы
    good = [r for r in rows if isinstance(r, tuple) and len(r) == 5]
    direct = _direct(good)
    row_names = [r[2] for r in good]
    rep = None
    try:
        rep = repair()
        if rep["всего"] != len(rows):
            out.append(f"repair называет всего {rep['всего']}, строк в таблице {len(rows)}")
        if rep["держится"] != sum(direct):
            out.append(f"repair называет держится {rep['держится']}, прямой прогон строк — {sum(direct)}")
        if rep["держится"] + len(rep["нарушено"]) != rep["всего"]:
            out.append(f"repair: держится {rep['держится']} + нарушено {len(rep['нарушено'])} ≠ всего {rep['всего']}")
        broken_direct = sorted(n for n, ok in zip(row_names, direct) if not ok)
        broken_rep = sorted(_base_name(n, row_names) for _, n, _ in rep["нарушено"])
        if broken_direct != broken_rep:
            out.append(f"состав упавших расходится: прямой прогон {broken_direct}, repair {broken_rep}")
    except Exception as e:  # noqa: BLE001
        out.append(f"repair() отказала на здоровой таблице: {type(e).__name__}: {e}")
    synth = ((1, "3", "держится", "точное", lambda: True),
             (1, "3", "ложь", "точное", lambda: False),
             (1, "3", "сбой", "сертификат", _synthetic_failure))
    try:
        sr = repair(synth)
        if sr["всего"] != 3:
            out.append(f"синтетическая таблица: repair называет всего {sr['всего']} вместо 3")
        if sr["держится"] != 1:
            out.append(f"синтетическая таблица: repair называет держится {sr['держится']} вместо 1")
        if sr["держится"] + len(sr["нарушено"]) != sr["всего"]:
            out.append("синтетическая таблица: держится + нарушено ≠ всего")
        названы = sorted(_base_name(n, ("держится", "ложь", "сбой")) for _, n, _ in sr["нарушено"])
        if названы != ["ложь", "сбой"]:
            out.append(f"синтетическая таблица: упавшие названы как {названы} вместо ['ложь', 'сбой']")
        else:
            сбой = [n for _, n, _ in sr["нарушено"] if n.startswith("сбой")][0]
            if "ValueError" not in сбой:
                out.append(f"синтетическая таблица: у сбойной строки нет имени исключения: {сбой!r}")
            if [k for _, _, k in sr["нарушено"]] != ["точное", "сертификат"]:
                out.append("синтетическая таблица: род проверки не донесён до отчёта")
    except Exception as e:  # noqa: BLE001
        out.append(f"синтетическая таблица: repair(rows) отказала: {type(e).__name__}: {e}")
    # (3) двойственность самой таблицы — карта §3 покрыта
    present = {r[0] for r in good}
    for k, c in enumerate(STEPS_3, 1):
        if k not in present:
            out.append(f"шаг {k} §3 (класс {c}) без единой строки таблицы")
    if deep and rep is not None:
        if _direct(good) != direct:
            out.append("прямой прогон недетерминирован: второй проход расходится с первым")
        try:
            if repair() != rep:
                out.append("repair() недетерминирована: второй прогон расходится с первым")
            rev = repair(tuple(reversed(good)))
            if (rev["всего"], rev["держится"]) != (len(good), sum(direct)):
                out.append(f"обратный порядок строк меняет счёт: {rev['всего']} / {rev['держится']}")
            if sorted(_base_name(n, row_names) for _, n, _ in rev["нарушено"]) != sorted(
                    _base_name(n, row_names) for _, n, _ in rep["нарушено"]):
                out.append("обратный порядок строк меняет состав упавших")
        except Exception as e:  # noqa: BLE001
            out.append(f"глубокая дорога: repair отказала: {type(e).__name__}: {e}")
    return out
