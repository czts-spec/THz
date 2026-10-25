"""ФIOS — декларативный индекс правил O(1) по ключу (контракт, класс)².

Слово автора (01.09): «ФIOS чтобы НЕ ПОГРЯЗНУТЬ В IF/OR — декларативная
таблица. выискивает ошибки на раз-два». Устройство:

* ТАБЛИЦА — данные, не ветвления: строки `ROWS` разворачиваются в
  `TABLE[(ключ_a, ключ_b, op)] = (имя_правила, подсказка)`; сборщик
  падает на противоречии (двух разных правилах в одной ячейке) —
  таблица сама ловит ошибку составителя.
* Индекс НЕ вычисляет: правила — имена в реестре `RULE_IMPL`; правило
  вида `post_check:*` — воркер, которому разрешено смотреть значение
  (00-ЛИНИЯ п.2.4: значение-зависимые случаи не живут в таблице классов).
* Пустая ячейка — ИММУННЫЙ ОТКАЗ `Unresolved` с явными Ф-записями:
  несведённое не получает правдоподобный класс.
* Аудит: `coverage()` называет числа (заполнено/отказных/пусто из
  знаменателя), `sweep()` прогоняет каждую ячейку на представителях
  классов — несостыковка интерфейсов всплывает одним прогоном.

Замер-мотив из корпуса: в `phi_transition_table_24x24.json` 526 пар из
576 в `add` сплющены в C=1 — контейнер знал одну линейку и стирал
остальные. Здесь сплющивание невозможно построением: ячейка либо несёт
названное правило, либо отказывает.
"""
from __future__ import annotations

from .rational import Q

from .errors import Unresolved
from .algebraic import AlgebraicNumber
from .bessel import BesselRootValue
from .container import PhiNumberV6
from .de_solution import DESolutionValue
from .elliptic import EValue, KValue
from .arctrig import ArcTrigValue
from .exp_log import ExpValue, LogValue
from .gamma_value import GammaValue
from .hyperbolic import HyperbolicValue
from .unity_root import UnityRoot
from .integral import IntegralValue
from .lambert import WValue
from .pi_value import PiMonomial
from .pointers import ScalePointer
from .product import InfiniteProductValue
from .series import PowerSeriesValue
from .dfinite import DFiniteSeries
from .zero_state import Zero, ZeroKind
from .zeta import ZetaValue
from .alpha import PhiAlpha, _is_transcendental
from .gridcore import Mixed

EXACT, UNCOMP = "EXACT", "UNCOMP"
#: контракт записи Φα над трансцендентной базой (план 03.09, шаг 5):
#: запись несёт то, чего линейки не несут — суммы над одной базой
#: (e + e², 9·(ln 2)², π + π²), и в онтологию входит СВОИМ ключом
#: (ALPHA, класс базы), а не мостом-переводом в носителя, которого нет
ALPHA = "ALPHA"

# обозначения ключей (контракт, класс) — как в корпусной паре осей
E0, E1, E3, E4 = (EXACT, 0), (EXACT, 1), (EXACT, 3), (EXACT, 4)
E5, E6, E7 = (EXACT, 5), (EXACT, 6), (EXACT, 7)
E8, E9, E10 = (EXACT, 8), (EXACT, 9), (EXACT, 10)
E11, E12 = (EXACT, 11), (EXACT, 12)
E13, E14 = (EXACT, 13), (EXACT, 14)
E15, E16, E17, E18 = (EXACT, 15), (EXACT, 16), (EXACT, 17), (EXACT, 18)
E21, E22 = (EXACT, 21), (EXACT, 22)
U20, U23 = (UNCOMP, 20), (UNCOMP, 23)
A5, A6, A10 = (ALPHA, 5), (ALPHA, 6), (ALPHA, 10)
#: структурная запись (gridcore.Mixed): значение без линейки-носителя —
#: дерево действий над носителями, засечка на сетке по запросу; класс 19
#: (предел засечек как способ получения величины). Слово автора 03.09:
#: отказ для разных баз НЕ принят — дискретное представление даёт значение
STRUCT = "STRUCT"
M19 = (STRUCT, 19)

#: живые ключи этой версии; класс 19 не носит собственного типа —
#: результат iteration_limit живёт в классах 1/3 (перевод в фасад §3);
#: класс 2 покрыт from_root линейки 3
LIVE = (E0, E1, E3, E4, E5, E6, E7, E8, E9, E10, E11, E12, E13, E14,
        E15, E16, E17, E18, E21, E22, U20, U23, A5, A6, A10, M19)

#: ключи величин с ВЕЩЕСТВЕННЫМ сертификатом сетки: между любыми двумя
#: действие сводится структурной записью (E4 — не веществен, E18 —
#: сертификат по точке t; они вне этого круга)
GRIDDED = (E1, E3, E5, E6, E7, E8, E9, E10, E11, E12, E13, E14, E15,
           E16, E17, A5, A6, A10, M19)

_OPS = "+-*/"


# ---------- классификация носителей ----------

def classify(x: object):
    """(контракт, класс) операнда; None — вне онтологии (отказ)."""
    t = type(x)
    if t is int or t is Q:
        return E1
    if t is AlgebraicNumber:
        return E3
    if t is ExpValue:
        return E5
    if t is LogValue:
        return E6
    if t is PiMonomial:
        return E10
    if t is UnityRoot:
        return E4
    if t is ArcTrigValue:
        return E11
    if t is HyperbolicValue:
        return E12
    if t is PowerSeriesValue or t is DFiniteSeries:
        return E7
    if t is InfiniteProductValue:
        return E8
    if t is IntegralValue:
        return E9
    if t is ZetaValue:
        return E13
    if t is GammaValue:
        return E14
    if t is KValue or t is EValue:
        return E15
    if t is BesselRootValue:
        return E16
    if t is WValue:
        return E17
    if t is DESolutionValue:
        return E18
    if t is ScalePointer:
        return (EXACT, x.phi_class())
    if t is Zero:
        return E0 if x.kind is ZeroKind.ABSOLUTE else U20
    if t is PhiNumberV6:
        nc = int(x.number_class)
        return (UNCOMP, nc) if nc in (20, 23) else (EXACT, nc)
    if t is PhiAlpha and _is_transcendental(x.alpha):
        base = classify(x.alpha)
        return (ALPHA, base[1]) if base is not None else None
    if t is Mixed:
        return M19
    return None


def _rec(x: object) -> str:
    r = getattr(x, "record", None)
    if callable(r):
        return r()
    if type(x) in (int, Q):
        return f"Φ(C=1, value={Q(x)})"
    return f"Φ(?={x!r})"


# ---------- реестр правил (имя → воркер) ----------

def _try_ops(a, b, op, order):
    """Бинарный протокол в заданном порядке сторон; NotImplemented с обеих
    сторон или чужой TypeError — None (решает вызывающий)."""
    dunder = {"+": "__add__", "-": "__sub__", "*": "__mul__",
              "/": "__truediv__"}[op]
    rdunder = "__r" + dunder[2:]
    for who, name in ((a, dunder), (b, rdunder)) if order == "lr" \
            else ((b, rdunder), (a, dunder)):
        fn = getattr(type(who), name, None)
        if fn is None:
            continue
        try:
            r = fn(who, b if who is a else a)
        except TypeError:
            continue
        if r is not NotImplemented:
            return r
    return None


def _delegate(a, b, op, hint):
    r = _try_ops(a, b, op, "lr")
    if r is None:
        raise Unresolved(op, _rec(a), _rec(b), hint=hint or
                         "операторы носителей не сводят пару")
    return r


def _delegate_r(a, b, op, hint):
    r = _try_ops(a, b, op, "rl")
    if r is None:
        raise Unresolved(op, _rec(a), _rec(b), hint=hint or
                         "операторы носителей не сводят пару")
    return r


def _to_container(x):
    if type(x) is PhiNumberV6:
        return x
    if type(x) is int:
        return PhiNumberV6.from_int(x)
    return PhiNumberV6.from_ratio(x.numerator, x.denominator)


def _scalar_exact(a, b, op, hint):
    """Класс 1 — контейнер v6 как рабочий скаляр (миграция, этап A)."""
    return _delegate(_to_container(a), _to_container(b), op, hint)


def _inv_scale_r(a, b, op, hint):
    """число / величина = inverse(величина) · число."""
    return b.inverse() * a


def _zero_left_div(a, b, op, hint):
    """ноль-состояние / x (x — не ноль): род сохраняется; невычисленный —
    отказ; 0/0 отказывает таблицей (строка refuse, 05.09)."""
    if a.kind is ZeroKind.NOT_COMPUTED:
        raise Unresolved(op, _rec(a), _rec(b),
                         hint="вычисление не производилось")
    return a


def _refuse(a, b, op, hint):
    raise Unresolved(op, _rec(a), _rec(b), hint=hint or "правила нет")


def _norm_hyper(x):
    """Мост носителей 21/22: HYPER контейнера v6 читается как указатель
    k = ±порядок градуированной записи — с 04.09 контейнер несёт поле
    порядка (решение 5.1 вычислением), мост переводит его один в один."""
    if type(x) is PhiNumberV6 and (x.is_infinity or x.is_infinitesimal):
        return ScalePointer(Q(x.m1.v, x.m2.v), x.hyper_order)      # поле порядка контейнера (04.09): ±k
    return x


def _hyper_bridge(a, b, op, hint):
    return _delegate(_norm_hyper(a), _norm_hyper(b), op, hint)


def _carrier_or_refuse(x, op, other):
    """Запись Φα — носителем; нет носителя — отказ, называющий запись."""
    if type(x) is PhiAlpha:
        try:
            return x.as_carrier()
        except (ArithmeticError, Unresolved) as ex:
            raise Unresolved(op, _rec(x), _rec(other),
                             hint=f"запись над своей базой носителя не "
                                  f"имеет: {ex}") from ex
    return x


def _alpha_bridge(a, b, op, hint):
    """Две записи Φα: одна база — ось E (Гаусс, сложение показателей);
    разные базы одного класса — переводом в носителей (нет носителя —
    отказ поимённо, разные базы общего основания не имеют)."""
    if type(a) is PhiAlpha and type(b) is PhiAlpha and not a._same_base(b):
        return route(_carrier_or_refuse(a, op, b),
                     _carrier_or_refuse(b, op, a), op)
    return _delegate(a, b, op, hint)


def _lift(c):
    """Носитель классов 5/6/10 с коэффициентом-парой — записью Φα над
    первообразной базой своего класса (канон записи сведёт e^y к (e^1)^y,
    ln 8 к 3·ln 2)."""
    if type(c) is ExpValue:
        return PhiAlpha(c.coeff, ExpValue.make(1, c.y), 1)
    if type(c) is PiMonomial and type(c.coeff) is Q:
        return PhiAlpha(c.coeff, PiMonomial(1, c.e), 1)
    if type(c) is LogValue:
        return PhiAlpha(c.coeff, LogValue.make(1, c.arg), 1)
    return None


def _alpha_lift(a, b, op, hint):
    """Запись Φα и носитель её класса: носитель поднимается в запись над
    той же базой, дальше — ось E. Иная база (ln 3 при базе ln 2, 2e при
    дробном показателе) — отказ."""
    ra = a if type(a) is PhiAlpha else _lift(a)
    rb = b if type(b) is PhiAlpha else _lift(b)
    if ra is None or rb is None:
        raise Unresolved(op, _rec(a), _rec(b),
                         hint="носитель не поднимается в запись: "
                              "коэффициент — не пара")
    if not ra._same_base(rb):
        raise Unresolved(op, _rec(a), _rec(b),
                         hint="разные базы-носители: общего основания нет")
    return _delegate(ra, rb, op, hint)


def _mixed(a, b, op, hint):
    """Структурная запись: дерево действий, засечка по запросу. Ничего не
    вычисляется и не округляется здесь — значение есть строение."""
    return Mixed.combine(a, b, op)


RULE_IMPL = {
    "delegate": _delegate,
    "delegate_r": _delegate_r,
    "scalar_exact": _scalar_exact,
    "inv_scale_r": _inv_scale_r,
    "zero_left_div": _zero_left_div,
    "refuse": _refuse,
    "hyper_bridge": _hyper_bridge,
    # воркеры со взглядом на значение (00-ЛИНИЯ п.2.4) — делегирование
    # носителю, который внутри решает по значению и честно отказывает:
    "post_check:hyper": _hyper_bridge,
    "post_check:unity": _delegate,
    "post_check:unity_sum": _delegate,
    "post_check:arctan_sum": _delegate,
    "post_check:zero": _delegate,
    "post_check:exp_add": _delegate,
    "post_check:log_add": _delegate,
    "post_check:gamma_mul": _delegate,
    "post_check:alpha": _alpha_bridge,
    "post_check:alpha_lift": _alpha_lift,
    "post_check:series": _delegate,
    "post_check:product": _delegate,
    "post_check:integral": _delegate,
    "post_check:elliptic": _delegate,
    "post_check:mixed": _mixed,
}

#: правила, которым разрешено отвечать отказом ПО ЗНАЧЕНИЮ операндов
_VALUE_DEPENDENT = tuple(n for n in RULE_IMPL if n.startswith("post_check:"))


# ---------- сама таблица: строки-данные ----------
# строка: (ключ_a, ключ_b, строка_операций, имя_правила, подсказка)

ROWS = [
    # класс 1 — скаляры на контейнере v6
    (E1, E1, "+-*/", "scalar_exact", ""),

    # линейка 3 POLY_ROOT
    (E3, E3, "+-*/", "delegate", ""),
    (E3, E1, "+-*/", "delegate", ""),
    (E1, E3, "+-*", "delegate", ""),
    (E1, E3, "/", "inv_scale_r", ""),

    # мономы π (класс 10)
    (E10, E10, "*/", "delegate", ""),
    (E10, E10, "+-", "post_check:hyper", "сведение только одноимённой степени π"),
    (E10, E1, "*/", "delegate", ""),
    (E1, E10, "*/", "delegate_r", ""),
    (E10, E1, "+-", "refuse", "π независим: правила нет"),
    (E1, E10, "+-", "refuse", "π независим: правила нет"),
    (E10, E3, "*/", "delegate", ""),
    (E3, E10, "*/", "delegate_r", ""),
    (E10, E3, "+-", "refuse", "π независим: правила нет"),
    (E3, E10, "+-", "refuse", "π независим: правила нет"),

    # класс 4 UNITY_ROOT — одометр: умножение всегда сведено (цикл)
    (E4, E4, "*/", "delegate", ""),
    (E4, E1, "*/", "post_check:unity", "сведены только скаляры 0, ±1"),
    (E1, E4, "*", "post_check:unity", "сведены только скаляры 0, ±1"),
    (E4, E4, "+-", "post_check:unity_sum",
     "вещественный итог — 2cos(2πa/n) (класс 1/3) или 0; комплексный — пара (re, im) над линейками 1/3 (структурная запись)"),
    (E4, E1, "+-", "post_check:unity_sum", "итог не веществен — пара (re, im): re = cos + q, im = sin"),
    (E1, E4, "+-", "post_check:unity_sum", "итог не веществен — пара (re, im)"),

    # класс 11 ARCTRIG — единственное направленное сложение (дополнение)
    (E11, E11, "+", "post_check:arctan_sum",
     "сведено только дополнение arctan(x)+arctan(1/x), x>0"),
    (E11, E11, "-*/", "refuse", "правила нет"),
    (E11, E1, "+-*/", "refuse", "правила нет"),
    (E1, E11, "+-*/", "refuse", "правила нет"),

    # класс 12 HYPERBOLIC — структурная пара над классом 5; арифметика
    # носителей не сведена (тождество — методом pythagorean_check)
    (E12, E12, "+-*/", "refuse", "правила нет; тождество — pythagorean_check"),
    (E12, E1, "+-*/", "refuse", "правила нет"),
    (E1, E12, "+-*/", "refuse", "правила нет"),

    # класс 5 EXPONENTIAL
    (E5, E5, "*/", "delegate", ""),
    (E5, E5, "+-", "post_check:exp_add", "сведение только одноимённого показателя"),
    (E5, E1, "*/", "delegate", ""),
    (E1, E5, "*/", "delegate_r", ""),
    (E5, E1, "+-", "refuse", "e^y независим: правила нет"),
    (E1, E5, "+-", "refuse", "e^y независим: правила нет"),

    # класс 6 LOGARITHMIC
    (E6, E6, "+-", "post_check:log_add", "ln складывается умножением аргументов"),
    (E6, E1, "*/", "delegate", ""),
    (E1, E6, "*", "delegate_r", ""),
    (E6, E1, "+-", "refuse", "ln независим: правила нет"),
    (E1, E6, "+-", "refuse", "ln независим: правила нет"),

    # класс 14 GAMMA
    (E14, E14, "*/", "post_check:gamma_mul",
     "отражение Эйлера / одинаковый base; прочее — отказ"),
    (E14, E1, "*/", "delegate", ""),
    (E1, E14, "*", "delegate_r", ""),
    (E14, E14, "+-", "refuse", "сложение Γ не сведено"),
    (E14, E1, "+-", "refuse", "Γ независима: правила нет"),
    (E1, E14, "+-", "refuse", "Γ независима: правила нет"),

    # класс 7 POWER_SERIES: одна рекуррента — по a0; масштаб парой в a0
    (E7, E7, "+-", "post_check:series", "одна рекуррента и один r — по a0; иначе замыкание D-конечных рядов (05.09)"),
    (E7, E7, "*", "post_check:series", "свёртка Коши — замыкание D-конечных рядов (Стэнли 1980), рекуррента результата точна"),
    (E7, E7, "/", "refuse", "деление рядов: общего замыкания D-конечных нет (1/f — лишь частные случаи)"),
    (E7, E1, "*/", "delegate", ""),
    (E1, E7, "*", "delegate_r", ""),
    (E7, E1, "+-", "post_check:series", "число — ряд c + 0 + 0 + …; сумма замыканием D-конечных (05.09)"),
    (E1, E7, "+-", "post_check:series", "число — ряд c + 0 + 0 + …; сумма замыканием D-конечных (05.09)"),

    # класс 8 INFINITE_PRODUCT: слияние строений при одном старте; слота
    # коэффициента нет — масштаб парой только единицей
    (E8, E8, "*", "post_check:product", "слияние строений при одном старте n0; иначе — свёртка Коши рядов телескопа (7↔8, 05.09)"),
    (E8, E8, "/", "post_check:product", "слияние строений при одном старте n0; иначе отказ по теореме (частное рядов)"),
    (E8, E8, "+-", "post_check:product", "Π = p_{n0} + Σ (p_N − p_{N−1}) — сумма замыканием D-конечных рядов (7↔8 телескопом, 05.09)"),
    (E8, E1, "*/", "post_check:product", "1 — само строение; иное число — масштаб ряда телескопа"),
    (E1, E8, "*/", "post_check:product", "1/Π — обращение строением; q·Π, q/Π — ряд телескопа"),
    (E8, E1, "+-", "post_check:product", "число — ряд c + 0 + 0 + …; сумма замыканием D-конечных (7↔8, 05.09)"),
    (E1, E8, "+-", "post_check:product", "число — ряд c + 0 + 0 + …; сумма замыканием D-конечных (7↔8, 05.09)"),
    (E7, E8, "+-*", "post_check:series", "ряд и произведение — оба D-конечны (телескоп): сумма и свёртка замыканием (05.09)"),
    (E8, E7, "+-*", "post_check:product", "произведение и ряд — оба D-конечны (телескоп): сумма и свёртка замыканием (05.09)"),
    (E8, E9, "+-*", "post_check:product", "произведение и интеграл — оба D-конечны: сумма и свёртка замыканием (05.09)"),
    (E9, E8, "+-*", "post_check:integral", "интеграл и произведение — оба D-конечны: сумма и свёртка замыканием (05.09)"),

    # класс 9 INTEGRAL: линейность над общим отрезком, смежные отрезки
    (E9, E9, "+-", "post_check:integral", "общий отрезок либо смежные отрезки одной подынтегральной; иначе замыкание D-конечных рядов (05.09)"),
    (E9, E9, "*", "post_check:integral", "произведение интегралов — свёртка Коши D-конечных рядов (05.09)"),
    (E9, E9, "/", "refuse", "частное интегралов: общего замыкания деления D-конечных рядов нет"),
    (E9, E1, "*/", "delegate", ""),
    (E1, E9, "*", "delegate_r", ""),
    (E9, E1, "+-", "post_check:integral", "число — ряд c + 0 + 0 + …; сумма замыканием D-конечных (05.09)"),
    (E1, E9, "+-", "post_check:integral", "число — ряд c + 0 + 0 + …; сумма замыканием D-конечных (05.09)"),
    (E7, E9, "+-*", "post_check:series", "ряд и интеграл — оба D-конечны: сумма и свёртка замыканием (05.09)"),
    (E9, E7, "+-*", "post_check:integral", "интеграл и ряд — оба D-конечны: сумма и свёртка замыканием (05.09)"),

    # класс 15 ELLIPTIC: K = (π/2)·S_K, E = (π/2)·S_E — гипергеометрические ряды класса 7 (05.09)
    (E15, E15, "+-*", "post_check:elliptic", "π-моном × D-конечный ряд: сумма и произведение гипергеометрических частей замыканием (05.09)"),
    (E15, E1, "*/", "post_check:elliptic", "масштаб — в π-мономе: c·K = (cπ/2)·S_K"),
    (E1, E15, "*", "post_check:elliptic", "масштаб — в π-мономе"),
    (E15, E7, "*", "post_check:elliptic", "K·ряд = (π/2)·(S_K·ряд) — свёртка Коши замыканием"),
    (E7, E15, "*", "post_check:elliptic", "ряд·K = (π/2)·(ряд·S_K) — свёртка Коши замыканием"),
    (E15, E9, "*", "post_check:elliptic", "K·∫ = (π/2)·(S_K·∫) — оба D-конечны"),
    (E9, E15, "*", "post_check:elliptic", "∫·K — оба D-конечны"),
    (E15, E8, "*", "post_check:elliptic", "K·Π = (π/2)·(S_K·Π-телескоп)"),
    (E8, E15, "*", "post_check:elliptic", "Π·K — оба D-конечны"),

    # Указатели масштаба 21/22 (значение-зависимые случаи — внутри носителя).
    # Счётчик k этих классов и экспонента — одно и то же `[в корпусе]`
    # (`inbox/SKILL_op_geometry_phi_numbers.md:264`: «тавтология колеса»), и в
    # позиционной сетке тот же k есть НОМЕР КОЛЕСА (`wheels.Wheel.index`):
    # целое колесо k = e + 1 — класс 22, дробное k = −e — класс 21. Строки
    # (E21, E1, «+-») ниже потому и refuse: сумма КОЛЁС живёт в классе 1, а не
    # в указателе, и тег колеса не есть указатель ω/ε.
    (E21, E21, "+-*/", "post_check:hyper", ""),
    (E21, E22, "+-*/", "post_check:hyper", ""),
    (E22, E21, "+-*/", "post_check:hyper", ""),
    (E22, E22, "+-*/", "post_check:hyper", ""),
    (E21, E1, "*/", "hyper_bridge", ""),
    (E22, E1, "*/", "hyper_bridge", ""),
    (E1, E21, "*/", "hyper_bridge", ""),
    (E1, E22, "*/", "hyper_bridge", ""),
    (E21, E1, "+-", "refuse", "сумма указателя с числом живёт в поле порядка: hyper.Hyper"),
    (E22, E1, "+-", "refuse", "сумма указателя с числом живёт в поле порядка: hyper.Hyper"),
    (E1, E21, "+-", "refuse", "сумма указателя с числом живёт в поле порядка: hyper.Hyper"),
    (E1, E22, "+-", "refuse", "сумма указателя с числом живёт в поле порядка: hyper.Hyper"),

    # IGNORANCE (класс 23, липкий) — на контейнере
    (U23, U23, "+-*/", "scalar_exact", ""),
    (U23, E1, "+-*/", "scalar_exact", ""),
    (E1, U23, "+-*/", "scalar_exact", ""),
]

#: Отказ разных баз называет ТРЕБОВАНИЕ, а не факт (план 03.09, шаг 6 —
#: замер, не костыль): логарифму сложения нужно общее основание — обе
#: базы как рациональные степени одного α (модуль Котса, 1714).
_NO_COMMON_BASE = ("требуется общее основание (Котс): обе базы — "
                   "рациональные степени одного α; между e, π и ln a "
                   "оно не установлено (ln π ∈ ℚ? — открыто; независимость — "
                   "гипотеза Шануэля, не аппарат)")
_NO_COMMON_BASE_ALG = ("требуется общее основание (Котс): алгебраическое "
                       "число не есть рациональная степень трансцендентной "
                       "базы (Линдеман 1882) — записи одной базы нет; "
                       "сумма жила бы в поле с алгебраическими "
                       "коэффициентами, которого здесь нет")

# записи Φα над трансцендентной базой (план 03.09, шаг 5): свои строки
for _A, _E in ((A5, E5), (A6, E6), (A10, E10)):
    ROWS += [
        (_A, _A, "+-*/", "post_check:alpha",
         "одна база — ось E (Гаусс); разные базы — переводом либо отказ"),
        (_A, E1, "*/", "delegate", ""),
        (E1, _A, "*", "delegate_r", ""),
        (E1, _A, "/", "inv_scale_r", ""),
        (_A, E1, "+-", "post_check:alpha",
         "логарифм сложения; знак вычитания — сертификатом базы"),
        (E1, _A, "+-", "post_check:alpha",
         "логарифм сложения; знак вычитания — сертификатом базы"),
        (_A, _E, "+-*/", "post_check:alpha_lift",
         "носитель поднимается в запись над той же базой; иная база — отказ"),
        (_E, _A, "+-*/", "post_check:alpha_lift",
         "носитель поднимается в запись над той же базой; иная база — отказ"),
        (_A, E3, "*/", "refuse",
         "мантисса записи — пара: алгебраический множитель не несётся"),
        (E3, _A, "*/", "refuse",
         "мантисса записи — пара: алгебраический множитель не несётся"),
        (_A, E3, "+-", "refuse", _NO_COMMON_BASE_ALG),
        (E3, _A, "+-", "refuse", _NO_COMMON_BASE_ALG),
    ]
    for _B in (A5, A6, A10):
        if _B != _A:
            ROWS.append((_A, _B, "+-*/", "refuse", _NO_COMMON_BASE))
del _A, _E, _B

# нули (класс 0 и Center-0 класса 20): род решает носитель — post-check
for _Z in (E0, U20):
    for _X in LIVE:
        ROWS.append((_Z, _X, "+-*", "post_check:zero", ""))
        if _X in (E0, U20):
            # 0/0 — отказ таблицей, не вычисляющим правилом (свежие глаза
            # 05.09: 0/0 давал абсолютный ноль; как 0⁰ — ни нуля, ни единицы)
            ROWS.append((_Z, _X, "/", "refuse", "0/0 не сведено: ни нуля, ни единицы (как 0⁰)"))
        else:
            ROWS.append((_Z, _X, "/", "zero_left_div", ""))
        if _X not in (E0, U20):
            ROWS.append((_X, _Z, "+-*", "post_check:zero", ""))
            ROWS.append((_X, _Z, "/", "refuse", "x/0 здесь не сведено"))
del _Z, _X


#: законность пары (контракт, класс) — правило корпуса (STUDY-08 §3: ZERO
#: только EXACT, UNCOMPUTABLE только UNCOMP, IGNORANCE — контракт незнания;
#: здесь APPROX снят, класс 23 живёт под UNCOMP); структурная запись —
#: только STRUCT 19; ALPHA — только классы базы 5/6/10 этой версии.
#: Незаконная пара НЕ нормализуется молча (дефект v6 project∘unproject,
#: CONTEXT §4) — она отвергается сборкой таблицы и ловится самопроверкой.
_ALPHA_CLASSES = (5, 6, 10)


def legal_key(key) -> bool:
    if not (isinstance(key, tuple) and len(key) == 2):
        return False
    contract, cls = key
    if not isinstance(cls, int) or not 0 <= cls <= 23:
        return False
    if cls == 0:
        return contract == EXACT
    if cls in (20, 23):
        return contract == UNCOMP
    if cls == 19:
        return contract == STRUCT
    if contract == ALPHA:
        return cls in _ALPHA_CLASSES
    return contract == EXACT


def _build(rows):
    table = {}
    for ka, kb, ops, rule, hint in rows:
        if rule not in RULE_IMPL:
            raise ValueError(f"неизвестное правило {rule!r} в строке {ka}×{kb}")
        for k in (ka, kb):
            if not legal_key(k):
                raise ValueError(f"незаконная пара (контракт, класс) {k} в строке {ka}×{kb}: "
                                 "ZERO — только EXACT, 20/23 — только UNCOMP, 19 — только STRUCT, ALPHA — 5/6/10")
        for op in ops:
            key = (ka, kb, op)
            if key in table and table[key][0] != rule:
                raise ValueError(
                    f"противоречие таблицы: {key} → {table[key][0]} и {rule}")
            table[key] = (rule, hint)
    return table


TABLE = _build(ROWS)

#: Дискретное замыкание (слово автора 03.09, «отказ НЕ принят»): всякая
#: пара величин с вещественным сертификатом сводится структурной записью —
#: пустые и отказные ячейки среди GRIDDED получают правило mixed. Правила,
#: сводящие точно (линейки, ось E, мост), остаются впереди.
_MIXED_HINT = ("структурная запись: линейки-носителя нет, значение — строение "
               "и засечка на сетке по запросу; равенство в пределе не утверждается")
for _ka in GRIDDED:
    for _kb in GRIDDED:
        for _op in _OPS:
            _cell = TABLE.get((_ka, _kb, _op))
            if _cell is None or _cell[0] == "refuse":
                TABLE[(_ka, _kb, _op)] = ("post_check:mixed", _MIXED_HINT)
del _ka, _kb, _op, _cell


# ---------- маршрутизация ----------

def _norm_alpha(x):
    """Мост записи с произвольным α (фундамент 02.09): Φ=m·α^e над
    рациональной/алгебраической/особой базой переводится в носителя её
    класса — то же направление, что у _norm_hyper: запись в линейку, не
    наоборот; несведённый перевод отказывает сам, поимённо. Запись над
    ТРАНСЦЕНДЕНТНОЙ базой остаётся записью — у неё свой ключ (ALPHA,
    класс базы) и свои строки (шаг 5): линейка её суммы не несёт."""
    if type(x) is PhiAlpha and not _is_transcendental(x.alpha):
        return x.as_carrier()
    return x


def route(a: object, b: object, op: str):
    """Свести действие по таблице. Пустая ячейка — иммунный отказ."""
    a, b = _norm_alpha(a), _norm_alpha(b)
    ka, kb = classify(a), classify(b)
    if ka is None or kb is None:
        raise Unresolved(op, _rec(a), _rec(b),
                         hint="операнд вне онтологии 24 классов")
    cell = TABLE.get((ka, kb, op))
    if cell is None:
        raise Unresolved(op, _rec(a), _rec(b),
                         hint=f"ячейка ФIOS ({ka[0]},{ka[1]})×({kb[0]},{kb[1]})"
                              f" «{op}» пуста")
    rule, hint = cell
    return RULE_IMPL[rule](a, b, op, hint)


# ---------- аудит ----------

def coverage():
    """Числа по осям: заполнено / из них отказных / пусто; знаменатель —
    |LIVE|²·4 (зеро-инструмента: неполнота объявляется числом)."""
    denom = len(LIVE) * len(LIVE) * len(_OPS)
    ruled = len(TABLE)
    refusals = sum(1 for r, _ in TABLE.values() if r == "refuse")
    per_op = {}
    for op in _OPS:
        cells = [k for k in TABLE if k[2] == op]
        per_op[op] = {
            "заполнено": len(cells),
            "отказных": sum(1 for k in cells if TABLE[k][0] == "refuse"),
        }
    return {"знаменатель": denom, "заполнено": ruled,
            "отказных": refusals, "пусто": denom - ruled, "по_действиям": per_op}


def self_check(table=None):
    """ФIOS ищет ошибки В СЕБЕ — инварианты таблицы как данных.

    Ловит: (1) асимметрию коммутативных действий — ячейка (a,b) считает,
    а зеркальная (b,a) пуста или отказывает; (2) дыры нулевой оси —
    живой класс без правил взаимодействия с нулём-состоянием; (3) отказ
    без подсказки — немой отказ запрещён. Плюс два инварианта, которые
    ловятся ещё СБОРКОЙ (_build): противоречие ячейки и неизвестное имя
    правила. Возвращает список находок; пусто — таблица согласована.

    Предел честности: самопроверка ловит НЕСОГЛАСОВАННОСТЬ слоёв, но не
    ошибку, одинаково внесённую в таблицу и в носитель; такие ловятся
    только независимым путём — двойственностью-приёмкой линейки и
    внешней сверкой (поэтому они в проекте обязательны).
    """
    t = TABLE if table is None else table
    findings = []
    for (ka, kb, op), (rule, hint) in t.items():
        if op in "+*" and rule != "refuse":
            mirror = t.get((kb, ka, op))
            if mirror is None:
                findings.append((f"{ka}×{kb} «{op}»",
                                 "асимметрия: зеркальная ячейка пуста"))
            elif mirror[0] == "refuse":
                findings.append((f"{ka}×{kb} «{op}»",
                                 "асимметрия: зеркало отказывает, ячейка считает"))
        if rule == "refuse" and not hint:
            findings.append((f"{ka}×{kb} «{op}»", "немой отказ: нет подсказки"))
    for x in LIVE:
        for z in (E0, U20):
            if not any((z, x, op) in t for op in _OPS):
                findings.append((f"{z}×{x}", "дыра нулевой оси"))
    for (ka, kb, op) in t:
        for k in (ka, kb):
            if not legal_key(k):
                findings.append((f"{ka}×{kb} «{op}»", f"незаконная пара (контракт, класс) {k}"))
    # псевдоним `phicore.fios` (латинское имя, 03.09) ничего не несёт сам:
    # его единственный инвариант — тождество с этим модулем и с этой таблицей;
    # проверяется здесь, потому что в самом fios.py определение жить не может
    # (sys.modules отдаёт под его именем ЭТОТ модуль) — ось самоприменения 09.10
    try:
        import importlib as _il
        import sys as _s
        alias = _il.import_module("phicore.fios")
        if alias is not _s.modules.get(__name__):
            findings.append(("fios", "псевдоним fios не тождествен модулю Фios"))
        elif getattr(alias, "TABLE", None) is not TABLE:
            findings.append(("fios", "псевдоним fios несёт другую таблицу"))
    except Exception as e:                                       # noqa: BLE001
        findings.append(("fios", f"псевдоним fios не ввозится: {type(e).__name__}: {e}"))
    return findings


def sweep(samples: dict):
    """Прогон каждой ячейки на представителях классов.

    samples: ключ → фабрика представителя. Возвращает список несостыковок:
    ячейка обещает вычисление, а прогон падает чужим исключением или
    отказом (для не-post-check правил); ячейка «refuse» — а прогон
    вычисляет. Пустой список — таблица согласована с носителями.
    """
    problems = []
    for (ka, kb, op), (rule, hint) in sorted(TABLE.items(),
                                             key=lambda kv: repr(kv[0])):
        if ka not in samples or kb not in samples:
            continue
        a, b = samples[ka](), samples[kb]()
        try:
            RULE_IMPL[rule](a, b, op, hint)
            outcome = "вычислено"
        except Unresolved:
            outcome = "отказ"
        except ArithmeticError:
            outcome = "отказ"  # честные отказы контейнера — тоже ответ
        except Exception as e:  # noqa: BLE001 — аудит ловит любой сбой
            problems.append(((ka, kb, op), rule,
                             f"сбой {type(e).__name__}: {e}"))
            continue
        if rule == "refuse" and outcome != "отказ":
            problems.append(((ka, kb, op), rule, "обещан отказ, вычислено"))
        elif rule not in _VALUE_DEPENDENT and rule != "refuse" \
                and outcome == "отказ":
            problems.append(((ka, kb, op), rule, "обещано вычисление, отказ"))
    return problems
