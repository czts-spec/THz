"""
container.py — контейнер Φ-числа v6.0 «очищенная структура» в phicore.

ПЕРЕНОС из корпуса: inbox/raw/phi_number_v6 (2).py (версия 1.1.0 со
structural_kind от свипа ассоциативности 2026-07-26) — байт-в-байт по
семантике; изменены только эта шапка и глосса класса 23 (навык, не
невежество — слово автора 02.09). Слово автора 01.09.2026: рабочий
скаляр библиотеки — этот контейнер, не чужой универсальный тип
(тот выселен тотально 02.09; собственная пара целых — rational.Q,
контейнер — эталон семантики поверх неё). Мантисса —
пара целых (m1, m2): любое рациональное точно; ось экспоненты E несёт
масштаб и корни (√2 = 2^(1/2)); биты блокировки на компонентах; отказ
ExactRepresentationNotAvailable / CrossBaseOperationNotSupported вместо
тихого приближения. approx_value() — явная трансляция наружу, не
вычисление внутри.

Исходная шапка перенесённого файла:
phi_number_v6.py — Φ-число v6.0 "очищенная структура"
Версия: 1.1.0 — РЕКОНСТРУИРОВАНО 2026-07-25 (см. RECOVERY_LOG.md)
Автор спецификации: Cu2ZnSnS4@gmail.com (проект OP-геометрия / Φ-числа)

ВАЖНО — ПРОИСХОЖДЕНИЕ ЭТОЙ ВЕРСИИ:
Оригинальная v1.1.0 (887 строк, 95/95 тестов) была утрачена при сбросе
песочницы 2026-07-25. /mnt/project/ также оказался недоступен в этой сессии,
поэтому байтовое восстановление невозможно. Данный файл — ПЕРЕСБОРКА по
точной архитектурной спецификации, зафиксированной в резюме разговора
(версия v1.0.0 была дословно извлечена из /mnt/transcripts/, остальное
дописано заново по описанию решений, которые пользователь уже подтвердил
2026-07-09). Функциональная эквивалентность оригиналу ПРОВЕРЯЕТСЯ, а не
предполагается: см. test_phi_number_v6.py (95 проверок) — это объективный
критерий приёмки согласно принципу "не гадать, писать тест".

Реализация: по трём источникам, загруженным пользователем 2026-07-09:
  [A] "2xIEEE754_Φ-число_v6_0_очищенная_структура_(без_знака_и_основания).pdf"
  [B] "Дополнение_к_теории_Φ-чисел_мультипликативное_представление_(Уровень_III).pdf"
  [C] "Приложение_3__ОП_геометрия_как_универсальная_модель.pdf"

ПОЧЕМУ НЕ contaminated/APPROX (архитектурное решение, 2026-07-09):
APPROX-Contract как отдельная координата ФIOS полностью дублируется общей
лестницей Mixed-Precision Contract; всё, что APPROX-Contract "может" — это
членство в Class 23 (IGNORANCE). PhiNumberV6 хранит ТОЛЬКО точные (или явно
заблокированные/UNCOMPUTABLE/HYPER) значения. Операция, которая не может
быть выполнена точно, поднимает ExactRepresentationNotAvailable либо (по
явному запросу через .approx_value()) отдаёт голый Python float — другой
тип, не претендующий быть точным. Смена точности — ТРАНСЛЯЦИЯ, а не скрытое
вычисление.

ОСНОВАНИЕ ЭКСПОНЕНТЫ α: НЕ ВСЕГДА 2 (уточнение, подтверждено 2026-07-09)
base_class/base_params описывают α рекурсивно (та же схема, что уже
проверена в InvariantTerm: NumberClass + class_params). base_class=None ⇒
неявная база 2, ноль накладных расходов (подавляющее большинство чисел).

ИСТОРИЯ ВЕРСИЙ:
  v1.0.0 (2026-07-09 05:45 UTC) — LockedInt, PhiNumberV6, mul/div/inverse/
    pow(int)/add(целая экспонента), base_class/base_params для α≠2,
    center_zero() для класса 20, check_boundary_signal().
  v1.1.0 (2026-07-09, пересобрано 2026-07-25) — ДОБАВЛЕНО:
    - ZERO класс (void/chaos: nihil absolutum vs chaos ordinabile),
      с правилами поглощения (*,/) и уступки (+,-);
    - HYPER класс (ω/ε — INFINITY/INFINITESIMAL): полная 8-образная таблица
      арифметики, включая критичный фикс "ε+0=ε" (эпсилон не должен
      схлопываться в 0 при сложении с точным нулём — независимо найден
      как БАГ v3.10.3 в libphi, воспроизведён и исправлен здесь тоже);
    - IGNORANCE класс (from_float: nan/inf/конечное) — липкий/поглощающий
      под всеми операциями, битово точное хранение конечных float;
    - ИСПРАВЛЕНА семантика UNCOMPUTABLE в сложении: участвует КАК 0
      (не поглощает), возрождается в Center0 ТОЛЬКО если сумма точно 0
      (критичная регрессия — была БАГ-10 в libphi, где Center0-in-add
      ошибочно поглощал);
    - common-exponent alignment для сложения: раньше требовалось, чтобы
      КАЖДЫЙ операнд по отдельности имел целую экспоненту; теперь
      достаточно, чтобы РАЗНОСТЬ экспонент E_a-E_b была целой (тогда
      мантисса операнда с БОЛЬШЕЙ экспонентой домножается на 2^diff,
      результат хранит МЕНЬШУЮ экспоненту) — работает для базы 2, ГДЕ
      2^целое всегда целое; для экзотических оснований (α=e,π,...) этот
      трюк НЕ применяется (α^diff иррационально для целого diff≠0) —
      это архитектурное отличие, не баг: см. класс CrossBase-логику ниже.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from enum import IntEnum
from typing import Optional, Tuple, Dict, Any, List


# =============================================================================
# ИСКЛЮЧЕНИЯ — явные сигналы границы точности (вместо скрытого APPROX-флага)
# =============================================================================

# Канонический класс отказа — общий для всего пакета (шов двух
# одноимённых классов найден мостом этапа E и закрыт 02.09).
from .errors import ExactRepresentationNotAvailable, Unresolved  # noqa: F401


class CrossBaseOperationNotSupported(ArithmeticError):
    """Операция над двумя Φ-числами с РАЗНЫМИ основаниями α не сводится к
    покомпонентной операции без потери точности — нужна явная трансляция."""
    pass


# =============================================================================
# NumberClass — 24 онтологических класса (значения совпадают с libphi.NumberClass)
# =============================================================================

class NumberClass(IntEnum):
    ZERO = 0
    NATURAL_POWER = 1
    ALGEBRAIC_ROOT = 2
    POLY_ROOT = 3
    UNITY_ROOT = 4
    EXPONENTIAL = 5
    LOGARITHMIC = 6
    POWER_SERIES = 7
    INFINITE_PRODUCT = 8
    INTEGRAL = 9
    TRIGONOMETRIC = 10
    ARCTRIG = 11
    HYPERBOLIC = 12
    ZETA = 13
    GAMMA = 14
    ELLIPTIC = 15
    BESSEL_ROOT = 16
    LAMBERT_W = 17
    DE_SOLUTION = 18
    ITERATION_LIMIT = 19
    UNCOMPUTABLE = 20
    INFINITESIMAL = 21
    INFINITY = 22
    IGNORANCE = 23


# =============================================================================
# Вспомогательные целочисленные функции (без Q) — n-й корень
# =============================================================================

def _int_nth_root(n: int, k: int) -> int:
    """Наибольшее целое r такое, что r^n <= k. Чистая int-арифметика."""
    if k <= 0:
        return 0
    if n == 1:
        return k
    if n == 2:
        return math.isqrt(k)
    r = int(round(k ** (1.0 / n)))
    while True:
        r1 = ((n - 1) * r + k // max(r ** (n - 1), 1)) // n
        if r1 >= r:
            break
        r = r1
    return r


def _perfect_nth_root(n: int, k: int) -> Tuple[bool, int]:
    """(True, root) если k == root**n точно, иначе (False, 0). k<0 для нечётных n."""
    if n <= 0:
        raise ValueError(f"n должно быть > 0, получено {n}")
    if k == 0:
        return True, 0
    if k == 1:
        return True, 1
    if k < 0:
        if n % 2 == 0:
            return False, 0
        ok, r = _perfect_nth_root(n, -k)
        return (True, -r) if ok else (False, 0)
    r = _int_nth_root(n, k)
    for candidate in (r, r + 1):
        if candidate > 0 and candidate ** n == k:
            return True, candidate
    return False, 0


def _two_adic(n: int) -> Tuple[int, int]:
    """(k, u): n = 2^k · u, u нечётно — 2-адическое разложение целого n ≥ 0.

    Единственность пары (k, u) — единственность разложения на простые
    [внешний источник: Hardy–Wright, An Introduction to the Theory of
    Numbers, теорема 2]. n = 0 отдаёт (0, 0): у нуля разложения нет, и
    корень из него берёт `_perfect_nth_root(q, 0)` = (True, 0) как прежде."""
    if n < 0:
        raise ValueError(f"n должно быть ≥ 0, получено {n}")
    if n == 0:
        return 0, 0
    k = (n & -n).bit_length() - 1
    return k, n >> k


def _ieee_lift(x):
    """Точная сторона операнда для класса 23 (ieee.py, слово автора 05.09):
    класс 23 → его exact; указатель → одночлен Hyper своего порядка; число с
    точным отношением → константа; иное (α ≠ 2, дробная экспонента) → None
    («не сведено»); Hyper/Zero/int/Q — как есть."""
    from .ieee import _lift_exact
    if not isinstance(x, PhiNumberV6):
        return _lift_exact(x)
    from .hyper import Hyper
    from .rational import Q
    if x.is_ignorance:
        return x.base_params.get('exact')
    if x.is_infinity or x.is_infinitesimal:
        return Hyper({x.hyper_order: Q(x.m1.v, x.m2.v)})
    if x.is_zero_class or x.is_uncomputable:
        return None
    try:
        n, d = x.to_exact_int_ratio()
    except ExactRepresentationNotAvailable:
        return None
    return Hyper.const(Q(n, d))


def _ieee_params(val: float, a, b, op: str) -> dict:
    """base_params класса 23 с точной стороной: raw_float — как считает IEEE
    (не мешать), exact — точно (кольцо элементов), loss — что потеряно."""
    from .ieee import _exact_op, loss
    fa = _full_value(a) if isinstance(a, PhiNumberV6) else _ieee_lift(a)     # с поглощённым: канцеллер указателя
    fb = _full_value(b) if isinstance(b, PhiNumberV6) else _ieee_lift(b)
    ex = _exact_op(fa, fb, op)
    return {'raw_float': val, 'exact': ex, 'loss': loss(ex, val)}


def _ieee_kind(val: float) -> str:
    return 'nan' if math.isnan(val) else ('inf' if math.isinf(val) else 'finite_translated')


def _full_value(x):
    """Полная величина в чтении элемента: подъём + поглощённое (Hyper); None — не сведено."""
    from .hyper import Hyper
    from .zero_state import Zero as Zero_
    base = _ieee_lift(x)
    if isinstance(x, PhiNumberV6) and (x.is_zero_class or x.is_uncomputable):
        base = Hyper.zero()                         # состояние по величине — ноль (метка — запись пути)
    if isinstance(x, PhiNumberV6) and x.is_ignorance and isinstance(base, Zero_):
        return base                                 # nan-состояние: величины нет
    absorbed = x.base_params.get('absorbed') if isinstance(x, PhiNumberV6) else None
    if base is None or absorbed is None and 'absorbed' in (x.base_params if isinstance(x, PhiNumberV6) else {}):
        return None
    return base if absorbed is None else base + absorbed


def _carry_absorbed(new: "PhiNumberV6", a: "PhiNumberV6", b: "PhiNumberV6") -> "PhiNumberV6":
    """Состояние из двух состояний наследует поглощённое обоими (иначе запись
    потери исчезала бы при c0 + c0)."""
    if 'absorbed' not in a.base_params and 'absorbed' not in b.base_params:
        return new
    from dataclasses import replace
    from .hyper import Hyper
    xa, xb = a.base_params.get('absorbed'), b.base_params.get('absorbed')
    if ('absorbed' in a.base_params and xa is None) or ('absorbed' in b.base_params and xb is None):
        total = None
    else:
        total = (xa if xa is not None else Hyper.zero()) + (xb if xb is not None else Hyper.zero())
    return replace(new, base_params={**new.base_params, 'absorbed': total})


def _settle(new, a, b, op: str):
    """Единый канцеллер указателя: если полные величины операндов сведены,
    absorbed(new) = full(a) ∘ full(b) − lift(new) точно (0 — записи нет);
    иначе остаётся запись места поглощения (_absorb: «не сведено»)."""
    if not isinstance(new, PhiNumberV6) or new.is_ignorance:
        return new                                   # класс 23 несёт полноту в exact
    from dataclasses import replace
    from .hyper import Hyper
    from .zero_state import Zero as Zero_
    fa, fb = _full_value(a), (_full_value(b) if isinstance(b, PhiNumberV6) else _ieee_lift(b))
    if fa is None or fb is None or isinstance(fa, Zero_) or isinstance(fb, Zero_):
        return new
    try:
        if op == '+':
            exact = fa + fb
        elif op == '*':
            exact = fa * fb
        else:
            exact = fa / fb                          # обратный к многочлену — отказ Hyper → запись не сведена
    except Exception:
        return new
    if not isinstance(exact, Hyper):
        return new
    base = Hyper.zero() if (new.is_zero_class or new.is_uncomputable) else _ieee_lift(new)
    if base is None:
        return new
    absorbed = exact - base
    params = dict(new.base_params)
    if absorbed.terms:
        params['absorbed'] = absorbed
    else:
        params.pop('absorbed', None)
    return replace(new, base_params=params)


def _absorb(keeper: "PhiNumberV6", dropped) -> "PhiNumberV6":
    """Канцеллер указателя (05.09, письмо «IEEE error canceller»): чтение
    указателя поглощает младший масштаб (ω + 5 = ω, 5 + ε = 5, ω−ω + 5 —
    состояние); поглощённое не проглатывается, а записывается:
    base_params['absorbed'] — сумма поглощённого в чтении элемента (Hyper),
    keeper + absorbed = точная сумма. Как fl + r у IEEE."""
    from dataclasses import replace
    from .hyper import Hyper
    d = _full_value(dropped)
    prev = keeper.base_params.get('absorbed', Hyper.zero()) if 'absorbed' in keeper.base_params else Hyper.zero()
    total = None if (d is None or prev is None) else prev + d
    return replace(keeper, base_params={**keeper.base_params, 'absorbed': total})


def _is_hyper(x) -> bool:
    """Многочлен Лорана phicore.hyper.Hyper — без ввоза модуля (hyper не ввозит контейнер)."""
    return type(x).__name__ == "Hyper" and hasattr(x, "terms") and hasattr(x, "origin")


_RECORD_CLASS = re.compile(r"C=(\d+)")


def _exponent_class(text: str):
    """Класс носителя из его записи «Φ(C=<класс> …»; None — записи без класса."""
    m = _RECORD_CLASS.search(text)
    return int(m.group(1)) if m else None


def _pow_hint(power, text: str) -> str:
    """Подсказка отказа `**` — по КЛАССУ показателя, не одна на всех (строка 79, п. 2)."""
    cls = _exponent_class(text)
    if _is_hyper(power) or cls in (21, 22):
        return (f"показатель — величина поля порядка (класс {cls}); степень контейнера "
                "в указателе не определена (правило: отказ вместо ложного класса)")
    if cls == 1:
        return ("показатель — носитель класса 1 (пара целых); дверь рациональной степени "
                "p/q — .pow_rational(p, q)")
    if cls is None:
        return ("показатель — носитель без класса в записи; степень контейнера в нём "
                "не определена (правило: отказ вместо ложного класса)")
    return (f"показатель — носитель класса {cls}; степень контейнера в величине класса "
            f"{cls} не определена: показатель контейнера — пара целых e1/e2 "
            "(правило: отказ вместо ложного класса)")


def _reduce_ratio(num: int, den: int) -> Tuple[int, int]:
    """GCD-сокращение (num,den), den всегда положителен на выходе.

    НУЛЕВОЙ ЧИСЛИТЕЛЬ СОКРАЩАЕТСЯ ТОЖЕ: gcd(0, d) = d, значит 0/d → 0/1.
    До 12.09 здесь стоял сторож `if num != 0`, и разность показателей 0/4
    не становилась 0/1 — от этого отказывало сложение записей с РАВНЫМ
    дробным показателем (√2 + 8·√2, дыра устройства п. 1 главы 09:
    сведено 4 из 12 против 12 из 12 у Φα) и ломалось равенство
    (from_ratio(4,9).pow_rational(1,2) против from_ratio(2,3) — п. 2).
    Каким каноном: единственность записи объявлена конструктором
    from_mantissa_exponent (reduce=True — каноническая форма) и
    инвариантом __post_init__ (e2.v > 0, «каноническая форма»); gcd(0, d)
    = d есть свойство самого НОД, а не выбор — сторож отменял его ровно
    в том случае, где он единственно и нужен. den = 0 остаётся (0, 0):
    это вход вне области, его ловят конструкторы."""
    if den < 0:
        num, den = -num, -den
    g = math.gcd(abs(num), den)
    if g > 1:
        num //= g
        den //= g
    return num, den


# =============================================================================
# Род структурной неопределённости — перенос через все пути арифметики
# =============================================================================
#
# Канон (12.09, дыра устройства п. 7 главы 09): zero_state.py:34–36 о роде
# нуля — «STRUCTURAL_INFINITESIMAL, от ε−ε: идеал бесконечно малых. УСТУПАЕТ
# в сложении (нейтрален), В УМНОЖЕНИИ СОХРАНЯЕТ СВОЙ РОД»; докстринги
# structural_kind (:359) и _add_core (:1097) о том же: бесконечно малые
# образуют идеал в кольце конечных гипердействительных чисел. Идеал замкнут
# относительно + и ·, поэтому род переносится через сложение, умножение,
# деление на число и натуральную степень. Обращения в этом каноне НЕТ: 1/I
# из кольца конечных выходит вовсе, и род поднимается до неограниченного
# 'infinite' (_kind_after_inverse) — это решение о смысле, принятое по
# канону замкнутости идеала, а не по удобству. Неуточнённый род (None)
# читается консервативно как 'infinite' (докстринг center_zero) и остаётся
# None: обращение не добавляет к нему знания.


def _structural_kind(x: "PhiNumberV6") -> Optional[str]:
    """Род структурной неопределённости операнда; None — операнд не
    структурный Центр 0 (у абсолютного Центра 0 и у числа рода нет)."""
    return x.structural_kind if (x.is_uncomputable and x.zero_provenance) else None


def _merge_kind(ka: Optional[str], kb: Optional[str]) -> Optional[str]:
    """Старшинство родов при встрече двух операндов: 'infinite'
    (неограниченная неопределённость ω±ω) старше 'infinitesimal' (идеал
    ε±ε); None старшинства не даёт и не отбирает."""
    if ka == 'infinite' or kb == 'infinite':
        return 'infinite'
    if ka == 'infinitesimal' or kb == 'infinitesimal':
        return 'infinitesimal'
    return None


def _kind_after_inverse(kind: Optional[str]) -> Optional[str]:
    """Род после ОБРАЩЕНИЯ: идеал бесконечно малых замкнут относительно +
    и ·, но не относительно 1/x — обратное к заведомо бесконечно малому
    бесконечно мало́ уже не будет, оно неограниченно. Поэтому
    'infinitesimal' → 'infinite'; 'infinite' и None остаются собой."""
    return 'infinite' if kind == 'infinitesimal' else kind


#: Одна запись отказа на все восемь дверей арифметики контейнера (строка 79
#: очереди, п. 1): предмет, дверь, причина — причина берётся у канона
#: `rational.Q`, а не переписывается (две копии разошлись бы, как разошлись
#: `Notch.of` и `Notch.from_q` до 12.09 — §14 п. 10 главы 08).
_IEEE_DOOR = "контейнер, дверь {door}: {kind} {x!r} не принят — {why}"


def _lift_pair(x, door: str = "арифметика"):
    """Точная пара целых дома (int, Q) — контейнером; None — не пара.

    Канон: шапка модуля («рабочий скаляр библиотеки — этот контейнер…
    собственная пара целых — rational.Q, контейнер — эталон семантики
    поверх неё»), та же дверь у оси основания (alpha._as_container) и у
    фасада (dispatch.combine транслирует Q на входе — приёмка
    test_phicore_container.py «Q транслируется в контейнер на входе»), и
    уже действующий подъём пары в _ieee_lift/_settle. До 12.09 дверь
    самого контейнера пару не принимала, и определённое действие не
    выполнялось: `Φα(ω-носитель) + Φα(пара)` падало голым TypeError
    интерпретатора (дыра устройства п. 9 главы 09), хотя _special_sum
    переводит в контейнер соседнюю ветвь int/Q.
    float и bool НЕ поднимаются, а отбиваются ПОИМЁННО — одной записью на
    все восемь дверей арифметики (строка 79 очереди, п. 1; 24.09). До того
    дверь отвечала им NotImplemented, и замер 24.09 `[вычислено здесь]` по
    восьми дверям (+, r+, −, r−, ·, r·, /, r/) дал: названо 0 из 8 при float
    и 0 из 8 при bool — все шестнадцать входов уходили интерпретатору голым
    «unsupported operand type(s)», и лишь `**` называл предмет сам. Соседи
    закрыты 12.09 тем же приёмом: `gridcore.Notch.of` пускает float и bool
    в дверь `from_q`, где отказывает сам `Q`; `affine._unit` — одной записью
    на все поля. Текст ПРИЧИНЫ здесь не копируется — его говорит канон дома,
    конструктор пары `rational.Q` («float запрещён построением; явная дверь
    — Q.from_float», «bool не число»); эта дверь добавляет ПРЕДМЕТ (тип и
    repr) и ИМЯ ДВЕРИ. Ворота против float в dispatch.py:35 остаются — они
    о фасаде, эта запись — о самом контейнере. Прочие чужие типы по-прежнему
    получают None → NotImplemented: у них может быть своя правая дверь
    (Φα, носители), у float и bool — нет."""
    if isinstance(x, (bool, float)):
        from .rational import Q
        try:
            Q(x)                                   # канон отказывает сам — одна копия текста
        except TypeError as why:
            raise TypeError(_IEEE_DOOR.format(door=door, kind=type(x).__name__,
                                              x=x, why=why)) from None
        raise TypeError(_IEEE_DOOR.format(door=door, kind=type(x).__name__, x=x,
                                          why="канон rational.Q промолчал"))
    if isinstance(x, int):
        return PhiNumberV6.from_int(x)
    from .rational import Q
    if isinstance(x, Q):
        return PhiNumberV6.from_ratio(int(x.n), int(x.d))
    return None


# =============================================================================
# Запись и хэш — протокол дома (12.09: пп. 3, 4, 5, 6 главы 09)
# =============================================================================
#
# Канон записи — README phicore («всякий record()») и соседи по пакету:
# у линейки со своим устройством запись идёт её старшим знаком и совпадает
# с печатью — `Φα(M=…, α=…, E=…)` у оси основания (alpha.py),
# «⟦lo; hi⟧·2^−g» у засечки (gridcore.Notch), «log_α(1 + k·α^d)» у
# логарифма сложения (zech.py), «Φ(C=23 IEEE raw=…, точно=…, потеря=…)» у
# класса 23 (ieee.py) — и у всех `__repr__ = record`; носители одного
# класса печатают «Φ(C=<класс>, …)» (zero_state.Zero, pointers, rational.Q).
# Контейнер держит ВСЕ 24 класса и четыре компонента с битами блокировки,
# поэтому его знак — «Φv6(…)», класс назван внутри числом и именем, а
# полями печатается ровно то, что класс несёт: у нуля и Центра 0 — род, у
# указателей — коэффициент и порядок, у класса 23 — raw/точно/потеря
# (форма ieee.py, где raw есть сам образ IEEE, хранимый битово), у обычных
# чисел — пара мантиссы, основание, пара показателя.
# ФЛОАТ В ЗАПИСЬ НЕ ВХОДИТ: вместо приближения — точная пара, когда она
# существует (to_exact_int_ratio); approx_value() остаётся ЯВНОЙ дверью
# наружу, как объявлено шапкой модуля.
#
# Хэш — каноном Notch.__hash__ (12.09, gridcore.py): равные записи дают
# равный хэш, и считается он по СВОЕЙ записи, без второго операнда.

#: Служебные поля base_params: печатаются своими ветвями записи, не как
#: параметры основания (канцеллер указателя и точная сторона класса 23).
_SERVICE_PARAMS = ('absorbed', 'raw_float', 'exact', 'loss')

#: Метка «поля нет» для различителя оснований: отличима от None,
#: которое бывает законным ЗНАЧЕНием поля (absorbed=None — «не сведено»).
_NO_FIELD = object()


def _freeze(v):
    """Значение поля-словаря — хэшируемым ключом ТОГО ЖЕ различения.

    Равенство контейнера структурное (поле-в-поле), и base_params входит в
    него как словарь: равные словари обязаны дать равный ключ. Словарь —
    отсортированными парами, список/кортеж — покомпонентно, прочее — самим
    значением, когда оно хэшируемо (Q хэшируется совместимо с int, Hyper —
    по своим членам, ровно как сравнивается), иначе его Ф-записью. Порядок
    ключей словаря в ключ не входит — потому один словарь даёт один ключ,
    как бы он ни был собран."""
    if isinstance(v, dict):
        return ('{}',) + tuple(sorted((str(k), _freeze(x)) for k, x in v.items()))
    if isinstance(v, (list, tuple)):
        return ('[]',) + tuple(_freeze(x) for x in v)
    try:
        hash(v)
    except TypeError:
        rec = getattr(v, 'record', None)
        return ('rec', rec() if callable(rec) else repr(v))
    return v


def _params_record(p: Dict[str, Any]) -> str:
    """Параметры основания записью: ключи по имени, значения — record()
    носителя, иначе repr. Служебные поля (_SERVICE_PARAMS) не параметры
    основания и печатаются своими ветвями."""
    items = []
    for k in sorted(p):
        if k in _SERVICE_PARAMS:
            continue
        v = p[k]
        rec = getattr(v, 'record', None)
        items.append(f"{k}={rec() if callable(rec) else v!r}")
    return ", ".join(items)


def _base_text(x: "PhiNumberV6") -> str:
    """Основание α записью: «2» — неявная база; иначе «[C=<номер> <ИМЯ>:
    <параметры>]». Печать НЕ зовёт _resolve_base_float и потому работает на
    всех 24 классах основания: п. 6 главы 09 замерил 5 из 24 (отказывали
    19, в том числе ζ = 13 и fn='arcctg' у ARCTRIG — через approx_value)."""
    if x.base_class is None:
        return "2"
    return f"[C={int(x.base_class)} {x.base_class.name}: {_params_record(x.base_params)}]"


def _base_mismatch(a: "PhiNumberV6", b: "PhiNumberV6") -> str:
    """Чем ИМЕННО расходятся основания — п. 8 главы 09: прежде отказ
    печатал «self.base_class=1, other.base_class=1» (числа в сообщении
    РАВНЫ, а различие лежало в base_params и не печаталось вовсе), у
    сложения — ни слова о том, какие базы столкнулись. Называются: класс
    (номером и именем), затем ключи параметров с обоими значениями.

    Служебные поля записи (канцеллер указателя и точная сторона класса 23)
    называются ПОИМЁННО и с пометкой, что сами основания совпали: случай
    достижим (запись с поглощённым против такой же записи без него даёт
    равные base_class и разные base_params — ровно то место, где старое
    сообщение печатало «self.base_class=5, other.base_class=5» и больше
    ничего). Сторожа «различий не найдено» здесь нет: различитель зовут
    только когда класс либо параметры разошлись, и пустым его ответ быть не
    может — недостижимая ветвь есть украшение, а не диагноз."""
    out = []
    if a.base_class != b.base_class:
        an = "нет (неявная база 2)" if a.base_class is None else f"{int(a.base_class)} {a.base_class.name}"
        bn = "нет (неявная база 2)" if b.base_class is None else f"{int(b.base_class)} {b.base_class.name}"
        out.append(f"класс основания {an} против {bn}")
    for k in sorted(set(a.base_params) | set(b.base_params)):
        va, vb = a.base_params.get(k, _NO_FIELD), b.base_params.get(k, _NO_FIELD)
        if va is _NO_FIELD and vb is _NO_FIELD or va == vb:
            continue
        what = "служебное поле записи" if k in _SERVICE_PARAMS else "параметр"
        tail = " (сами основания совпали)" if k in _SERVICE_PARAMS else ""
        la = "нет" if va is _NO_FIELD else repr(va)
        lb = "нет" if vb is _NO_FIELD else repr(vb)
        out.append(f"{what} {k}: {la} против {lb}{tail}")
    return "; ".join(out)


# =============================================================================
# LockedInt — один "рациональный компонент" (m1, m2, e1 или e2)
# =============================================================================

@dataclass(frozen=True, slots=True)
class LockedInt:
    """Один компонент Φ-числа: целое произвольной точности v + независимый
    бит блокировки b. b=False: активен, участвует в операциях. b=True:
    заблокирован — при операциях передаётся в результат без изменений
    (см. правила PhiNumberV6.__mul__/__truediv__), другой операнд
    игнорируется."""
    v: int
    b: bool = False

    def __repr__(self) -> str:
        lock = "🔒" if self.b else ""
        return f"{self.v}{lock}"


# =============================================================================
# PhiNumberV6 — контейнер Φ-числа v6.0
# =============================================================================

@dataclass(frozen=True, slots=True)
class PhiNumberV6:
    """
    Φ = (M, E), M=((m1,bm1),(m2,bm2)), E=((e1,be1),(e2,be2)).
    val(Φ) = (m1/m2) · α^(e1/e2), где α=2 по умолчанию (base_class=None)
    либо задаётся через (base_class, base_params) для e, π и т.п.

    Знак числа — в знаке m1.v. Инварианты (проверяются в __post_init__):
      - m2.v > 0, e2.v > 0 (каноническая форма);
      - NumberClass.ZERO / UNCOMPUTABLE требуют m1.v == 0.

    zero_provenance — ТОЛЬКО для ZERO/UNCOMPUTABLE: False="абсолютный"
    (nihil absolutum / Center0 обычный), True="структурный"
    (chaos ordinabile / Center0 от сокращения ω-ω).

    ignorance_kind — ТОЛЬКО для IGNORANCE: 'nan'|'inf'|'finite_translated'.
    Для IGNORANCE фактическое float-значение хранится в
    base_params['raw_float'] (переиспользование поля — у IGNORANCE нет
    понятия "основание степени", поле свободно для другой роли).

    РАВЕНСТВО ЗДЕСЬ СТРУКТУРНОЕ — О ЗАПИСИ, НЕ О ВЕЛИЧИНЕ (п. 3 главы 09;
    прежде об этом молчали и класс, и модуль). `==` есть равенство
    @dataclass поле-в-поле: совпадают мантисса, показатель, класс,
    основание, провенанс, род и БИТЫ БЛОКИРОВКИ. Отсюда 2 и 1·2^1 — две
    записи одной величины и НЕ равны. Замер (свип главы `sweep_container_eq`,
    шесть записей двойки: from_int(2), 1·2^(1/1), 4/2, 4·2^(−1/1),
    2·2^(0/2), 2/1): 15 пар, равных по ВЕЛИЧИНЕ 15, равных по `==` 6 (было
    3 — три записи свело закрытие п. 2, канонизация нулевого показателя).
    Дверь к равенству по величине — явная и уже стоит: `to_exact_int_ratio()`
    (точная пара) либо `PhiRecord.from_container` (старшая запись Φ(s, α, M,
    E)); у соседа по главе `Φα` равенство ПО ЗНАЧЕНИЮ, у соседа по главе 08
    та же двойственность оговорена (`Notch.__eq__` о величине на сетке,
    `Mixed.__eq__` о записи). Смысл предмета этим не меняется: контейнер
    есть ЗАПИСЬ Φ-числа, и её различение — его работа.

    ХЭШ СОГЛАСОВАН С ЭТИМ РАВЕНСТВОМ (`__hash__`, п. 4 главы 09): считается
    по тем же полям, словарь основания — каноническим ключом `_freeze`.

    ЗАПИСЬ: `record()` — явная Ф-запись «Φv6(…, C=<класс>)» (п. 5 главы 09),
    и `__repr__` — ТОТ ЖЕ текст (канон соседей Φα, Notch, Zech, Ieee: одна
    запись на предмет, `__repr__ = record`); float не печатается (п. 6),
    трансляция наружу остаётся явной дверью `approx_value()`.
    """
    m1: LockedInt
    m2: LockedInt
    e1: LockedInt
    e2: LockedInt
    number_class: NumberClass = NumberClass.NATURAL_POWER
    base_class: Optional[NumberClass] = None
    base_params: Dict[str, Any] = field(default_factory=dict)
    zero_provenance: bool = False
    ignorance_kind: Optional[str] = None
    structural_kind: Optional[str] = None  # ТОЛЬКО при zero_provenance=True на
    # UNCOMPUTABLE: 'infinite' (от ω±ω — неограниченная неопределённость,
    # результат мог бы быть ЛЮБЫМ: конечным/бесконечным/бесконечно малым —
    # поглощает дальнейшую арифметику) | 'infinitesimal' (от ε±ε —
    # ОГРАНИЧЕННАЯ неопределённость: сумма/разность двух бесконечно малых
    # математически гарантированно остаётся бесконечно малой величиной или
    # нулём — бесконечно малые образуют идеал в кольце конечных
    # гипердействительных чисел, замкнутый относительно +/- — поэтому
    # уступает конечному/бесконечному числу так же, как обычный ε).

    def __post_init__(self):
        if self.m2.v <= 0:
            raise ValueError(f"m2.v должен быть > 0 (знак — в m1); получено m2.v={self.m2.v}")
        if self.e2.v <= 0:
            raise ValueError(f"e2.v должен быть > 0 (каноническая форма); получено e2.v={self.e2.v}")
        if self.number_class == NumberClass.ZERO and self.m1.v != 0:
            raise ValueError("NumberClass.ZERO требует m1.v == 0")
        if self.number_class == NumberClass.UNCOMPUTABLE and self.m1.v != 0:
            raise ValueError("NumberClass.UNCOMPUTABLE (Center 0) требует m1.v == 0")

    def __hash__(self) -> int:
        """Контракт eq/hash: равные записи — равный хэш (п. 4 главы 09).

        До 12.09 хэша не было вовсе при живом `==`: @dataclass(frozen=True)
        собирал его по кортежу полей, а поле base_params — словарь, и
        `hash(from_int(3))` отвечал TypeError «unhashable type: 'dict'»
        (замер: хэш брался у 0 записей из 36 пробных; эталон семантики в
        словарь не положить, в множество не собрать). Соседи хэшируемы все:
        `hash(A(3,2,1)) == hash(A(3,4,Q(1,2)))`, Notch, Mixed, Zero,
        ScalePointer, Ieee, Q.
        Каким каноном: `Notch.__hash__` (gridcore.py, 12.09) и `zech.py:172`
        — КАНОНОМ ЗАПИСИ, то есть по своим полям, без второго операнда,
        поэтому одинаков у всех равных, сколько бы их ни было. Равенство
        контейнера структурное (докстринг класса), значит в хэш входят ровно
        его поля; словарь основания — каноническим ключом `_freeze`
        (порядок ключей в ключ не входит). Замер после правки на пробном
        наборе 36 записей (все классы, состояния, блокировки, поглощённое,
        экзотические основания): 630 пар, равных 7, хэш расходится у 0; хэш
        берётся у 36 из 36 (было 0 из 36), словарь находит запись по равной
        ей, множество сводит равные в одну."""
        return hash(("Φv6", self.m1, self.m2, self.e1, self.e2,
                     self.number_class, self.base_class,
                     _freeze(self.base_params), self.zero_provenance,
                     self.ignorance_kind, self.structural_kind))

    # ── Классификация ──────────────────────────────────────────────────

    @property
    def sign(self) -> int:
        return 1 if self.m1.v < 0 else 0

    @property
    def exponent_is_integer(self) -> bool:
        return self.e1.v % self.e2.v == 0

    @property
    def m_pair_locked(self) -> bool:
        return self.m1.b or self.m2.b

    @property
    def e_pair_locked(self) -> bool:
        return self.e1.b or self.e2.b

    @property
    def is_zero_class(self) -> bool:
        return self.number_class == NumberClass.ZERO

    @property
    def is_void(self) -> bool:
        return self.number_class == NumberClass.ZERO and not self.zero_provenance

    @property
    def is_chaos_ordinabile(self) -> bool:
        return self.number_class == NumberClass.ZERO and self.zero_provenance

    @property
    def is_uncomputable(self) -> bool:
        return self.number_class == NumberClass.UNCOMPUTABLE

    @property
    def is_absolute_zero(self) -> bool:
        return self.number_class == NumberClass.UNCOMPUTABLE and not self.zero_provenance

    @property
    def is_structural_zero(self) -> bool:
        return self.number_class == NumberClass.UNCOMPUTABLE and self.zero_provenance

    @property
    def is_infinity(self) -> bool:
        return self.number_class == NumberClass.INFINITY

    @property
    def is_infinitesimal(self) -> bool:
        return self.number_class == NumberClass.INFINITESIMAL

    @property
    def is_ignorance(self) -> bool:
        return self.number_class == NumberClass.IGNORANCE

    # ── Значение (только для отображения/отладки — НЕ для точных вычислений) ──

    def approx_value(self) -> float:
        """Явная ТРАНСЛЯЦИЯ в обычный float. Возвращаемый тип — float, НЕ
        PhiNumberV6: это осознанный выход из точного контейнера, а не
        скрытая APPROX-контаминация внутри него."""
        if self.number_class == NumberClass.UNCOMPUTABLE:
            return 0.0
        if self.number_class == NumberClass.ZERO:
            return 0.0
        if self.number_class == NumberClass.IGNORANCE:
            return self.base_params.get('raw_float', float('nan'))
        if self.is_infinity:                       # 05.09: указатель масштаба — не 2^order (было ω ≈ 2.0)
            return math.copysign(float('inf'), self.m1.v)
        if self.is_infinitesimal:
            return math.copysign(0.0, self.m1.v)
        base = 2.0 if self.base_class is None else self._resolve_base_float()
        return (self.m1.v / self.m2.v) * (base ** (self.e1.v / self.e2.v))

    def _resolve_base_float(self) -> float:
        bc = self.base_class
        p = self.base_params
        if bc is None:
            return 2.0
        if bc == NumberClass.NATURAL_POWER:
            return float(p.get('value', 2))
        if bc == NumberClass.ALGEBRAIC_ROOT:
            q = p.get('root_degree', 2)
            rad = p.get('radicand', 2)
            return float(rad) ** (1.0 / q)
        if bc == NumberClass.EXPONENTIAL:
            return math.exp(float(p.get('arg', 1)))
        if bc == NumberClass.LOGARITHMIC:
            return math.log(float(p.get('arg', math.e)))
        if bc == NumberClass.ARCTRIG:
            fn = p.get('fn', 'arctan')
            arg = float(p.get('arg', 1))
            coeff = float(p.get('coeff', 1))
            fns = {'arctan': math.atan, 'arcsin': math.asin, 'arccos': math.acos}
            if fn not in fns:
                raise NotImplementedError(f"_resolve_base_float: ARCTRIG fn={fn!r} не поддержан")
            return coeff * fns[fn](arg)
        raise NotImplementedError(f"_resolve_base_float: base_class={bc.name} ещё не реализован")

    def to_exact_int_ratio(self) -> Tuple[int, int]:
        """Точное (num, den) как ЦЕЛЫЕ (не Q) — только когда основание
        неявное (2) и экспонента целая. Иначе — явное исключение."""
        if self.number_class in (NumberClass.UNCOMPUTABLE, NumberClass.ZERO,
                                  NumberClass.INFINITY, NumberClass.INFINITESIMAL,
                                  NumberClass.IGNORANCE):
            raise ExactRepresentationNotAvailable(
                f"класс {self.number_class.name} не имеет обычного рационального значения"
            )
        if self.base_class is not None:
            raise ExactRepresentationNotAvailable(
                f"base_class={self.base_class.name} задан (основание не 2) — точное (num,den) "
                f"в этой форме не определено; используйте .approx_value()"
            )
        if not self.exponent_is_integer:
            raise ExactRepresentationNotAvailable(
                f"экспонента e1/e2={self.e1.v}/{self.e2.v} не целая — значение иррационально "
                f"относительно базы 2; точного (num,den) не существует. Используйте .approx_value()"
            )
        n = self.e1.v // self.e2.v
        num, den = self.m1.v, self.m2.v
        if n >= 0:
            num *= 2 ** n
        else:
            den *= 2 ** (-n)
        if not (self.m1.b or self.m2.b):
            num, den = _reduce_ratio(num, den)
        return num, den

    # ── Конструкторы ─────────────────────────────────────────────────

    @classmethod
    def from_int(cls, n: int, number_class: NumberClass = NumberClass.NATURAL_POWER) -> "PhiNumberV6":
        return cls(LockedInt(n), LockedInt(1), LockedInt(0), LockedInt(1), number_class)

    @classmethod
    def from_ratio(cls, num: int, den: int, number_class: NumberClass = NumberClass.NATURAL_POWER) -> "PhiNumberV6":
        if den == 0:
            raise ValueError("den (m2) не может быть 0")
        num, den = _reduce_ratio(num, den)
        return cls(LockedInt(num), LockedInt(den), LockedInt(0), LockedInt(1), number_class)

    @classmethod
    def from_mantissa_exponent(cls, m_num: int, m_den: int, e_num: int, e_den: int,
                                number_class: NumberClass = NumberClass.NATURAL_POWER,
                                reduce: bool = True) -> "PhiNumberV6":
        if m_den == 0 or e_den == 0:
            raise ValueError("m_den и e_den не могут быть 0")
        if reduce:
            m_num, m_den = _reduce_ratio(m_num, m_den)
            e_num, e_den = _reduce_ratio(e_num, e_den)
        else:
            if m_den < 0:
                m_num, m_den = -m_num, -m_den
            if e_den < 0:
                e_num, e_den = -e_num, -e_den
        return cls(LockedInt(m_num), LockedInt(m_den), LockedInt(e_num), LockedInt(e_den), number_class)

    @classmethod
    def exotic_base(cls, m_num: int, m_den: int, e_num: int, e_den: int,
                     base_class: NumberClass, base_params: Dict[str, Any],
                     number_class: NumberClass) -> "PhiNumberV6":
        """Конструктор для α≠2. Точное (num,den)-представление НЕ определено
        для таких чисел (см. to_exact_int_ratio) — доступно только .approx_value()."""
        if m_den <= 0 or e_den <= 0:
            raise ValueError("m_den и e_den должны быть > 0 для exotic_base")
        return cls(LockedInt(m_num), LockedInt(m_den), LockedInt(e_num), LockedInt(e_den),
                   number_class=number_class, base_class=base_class, base_params=dict(base_params))

    @classmethod
    def void(cls) -> "PhiNumberV6":
        """Class 0 (ZERO), nihil absolutum — "чистая" пустота: 0+x=x, 0*x=0."""
        return cls(LockedInt(0, True), LockedInt(1, True), LockedInt(0, True), LockedInt(1, True),
                   number_class=NumberClass.ZERO, zero_provenance=False)

    @classmethod
    def chaos(cls) -> "PhiNumberV6":
        """Class 0 (ZERO), chaos ordinabile — "структурный" нуль (например,
        от неопределённого сокращения), отличается провенансом от void()."""
        return cls(LockedInt(0, True), LockedInt(1, True), LockedInt(0, True), LockedInt(1, True),
                   number_class=NumberClass.ZERO, zero_provenance=True)

    @classmethod
    def center_zero(cls, structural: bool = False, kind: Optional[str] = None) -> "PhiNumberV6":
        """Class 20 (UNCOMPUTABLE) — Центр 0. structural=False: абсолютный
        (напр. константа Хайтина). structural=True: возник из сокращения —
        kind='infinite' (ω±ω, неограниченная неопределённость, поглощает
        дальнейшую арифметику) | kind='infinitesimal' (ε±ε, ограниченная
        неопределённость — уступает конечному/бесконечному числу) |
        kind=None (структурный, но род неопределённости не уточнён —
        трактуется консервативно как 'infinite', см. _add_impl)."""
        return cls(LockedInt(0, True), LockedInt(1, True), LockedInt(0, True), LockedInt(1, True),
                   number_class=NumberClass.UNCOMPUTABLE, zero_provenance=structural,
                   structural_kind=kind if structural else None)

    @classmethod
    def omega(cls, coeff_num: int = 1, coeff_den: int = 1, order: int = 1) -> "PhiNumberV6":
        """Class 22 (INFINITY) — coeff·ω^order, масштабный указатель бесконечности.
        Поле порядка (решение 5.1 вычислением, 04.09): order ≥ 1 хранится в e1;
        ω·ω = ω², ω·ε = 1 — умножение ассоциативно (свип 1728 троек: 32 → 0)."""
        if not isinstance(order, int) or isinstance(order, bool) or order < 1:
            raise ValueError("порядок ω — целое ≥ 1 (ω^0 — обычное число, ω^−k — ε^k)")
        coeff_num, coeff_den = _reduce_ratio(coeff_num, coeff_den)
        return cls(LockedInt(coeff_num), LockedInt(coeff_den), LockedInt(order, True), LockedInt(1, True),
                   number_class=NumberClass.INFINITY)

    @classmethod
    def epsilon(cls, coeff_num: int = 1, coeff_den: int = 1, order: int = 1) -> "PhiNumberV6":
        """Class 21 (INFINITESIMAL) — coeff·ε^order = coeff·ω^(−order); order ≥ 1 в e1."""
        if not isinstance(order, int) or isinstance(order, bool) or order < 1:
            raise ValueError("порядок ε — целое ≥ 1")
        coeff_num, coeff_den = _reduce_ratio(coeff_num, coeff_den)
        return cls(LockedInt(coeff_num), LockedInt(coeff_den), LockedInt(order, True), LockedInt(1, True),
                   number_class=NumberClass.INFINITESIMAL)

    @property
    def hyper_order(self) -> int:
        """Знаковый порядок указателя: +k для coeff·ω^k, −k для coeff·ε^k, 0 — не указатель."""
        if self.is_infinity:
            return self.e1.v
        if self.is_infinitesimal:
            return -self.e1.v
        return 0

    @classmethod
    def from_hyper(cls, coeff_num: int, coeff_den: int, order: int) -> "PhiNumberV6":
        """coeff·ω^order по знаковому порядку: > 0 — ω^k, < 0 — ε^k, 0 — обычное число."""
        if order > 0:
            return cls.omega(coeff_num, coeff_den, order)
        if order < 0:
            return cls.epsilon(coeff_num, coeff_den, -order)
        return cls.from_ratio(coeff_num, coeff_den, NumberClass.NATURAL_POWER)

    @classmethod
    def from_float(cls, x: float) -> "PhiNumberV6":
        """Class 23 (IGNORANCE) — явная трансляция float→Φ-число-обёртка.
        Класс 23 — НАВЫК, не невежество (слово автора 02.09: нативный
        IEEE 754; 23 на входе — 23 на выходе, не мешать, ничего не
        делать). Ярлык «невежество» — инверсия смысла одной из моделей;
        фраза «невежда может не учиться, для него есть класс 23» — о
        пользователе, не о классе. raw_float хранит БИТОВО ТОЧНОЕ
        исходное значение (не приближение приближения)."""
        if math.isnan(x):
            kind = 'nan'
        elif math.isinf(x):
            kind = 'inf'
        else:
            kind = 'finite_translated'
        from .ieee import exact, loss
        ex = exact(x)
        return cls(LockedInt(0, True), LockedInt(1, True), LockedInt(0, True), LockedInt(1, True),
                   number_class=NumberClass.IGNORANCE, ignorance_kind=kind,
                   base_params={'raw_float': x, 'exact': ex, 'loss': loss(ex, x)})

    # ── Диагностика границы (переполнение/антипереполнение) ──────────

    def check_boundary_signal(self, min_v: int, max_v: int) -> List[str]:
        """Имена компонентов, которые ОДНОВРЕМЕННО заблокированы И лежат на
        границе [min_v, max_v]. Это СИГНАЛ переполнения/антипереполнения
        (источник [A]), НЕ смена number_class."""
        hits = []
        for name, comp in (('m1', self.m1), ('m2', self.m2), ('e1', self.e1), ('e2', self.e2)):
            if comp.b and comp.v in (min_v, max_v):
                hits.append(name)
        return hits

    # ── Отрицание (нужно для унификации __sub__ через __add__) ───────

    def __neg__(self) -> "PhiNumberV6":
        if self.number_class in (NumberClass.UNCOMPUTABLE, NumberClass.ZERO):
            return self  # -0 = 0 (обеих разновидностей), -Center0 = Center0
        if self.number_class == NumberClass.IGNORANCE:
            raw = self.base_params.get('raw_float')
            new_raw = -raw if raw is not None else raw
            from .ieee import loss
            from .hyper import Hyper
            ex = self.base_params.get('exact')
            ex = -ex if isinstance(ex, Hyper) else ex
            return PhiNumberV6(self.m1, self.m2, self.e1, self.e2,
                               number_class=self.number_class, base_class=self.base_class,
                               base_params={**self.base_params, 'raw_float': new_raw, 'exact': ex, 'loss': loss(ex, new_raw)},
                               zero_provenance=self.zero_provenance, ignorance_kind=self.ignorance_kind)
        params = dict(self.base_params)
        if params.get('absorbed') is not None:
            params['absorbed'] = -params['absorbed']          # канцеллер: −(x + поглощённое)
        return PhiNumberV6(LockedInt(-self.m1.v, self.m1.b), self.m2, self.e1, self.e2,
                           number_class=self.number_class, base_class=self.base_class,
                           base_params=params,
                           zero_provenance=self.zero_provenance, ignorance_kind=self.ignorance_kind)

    # ── Арифметика: умножение / деление ───────────────────────────────
    #
    # Приоритет (источник [A] + расширения HYPER/IGNORANCE, подтверждено
    # 2026-07-09): ZERO (поглощает безусловно) > UNCOMPUTABLE (поглощает)
    # > IGNORANCE (липкий) > HYPER (особые правила) > обычные числа.
    #
    # Для обычных чисел: заблокированная ПАРА (M целиком или E целиком —
    # см. обоснование в v1.0.0 docstring) передаётся в результат без
    # изменений; другой операнд игнорируется. Активные пары комбинируются
    # как обычно.

    def __mul__(self, other: "PhiNumberV6") -> "PhiNumberV6":
        return self._combine(other, op='mul')

    def __truediv__(self, other: "PhiNumberV6") -> "PhiNumberV6":
        return self._combine(other, op='div')

    def __rmul__(self, other):
        if _is_hyper(other):
            return self._mul_hyper(other, 'mul')
        lifted = _lift_pair(other, "__rmul__")         # пара дома — дверью (12.09)
        return NotImplemented if lifted is None else lifted._combine(self, 'mul')

    def __rtruediv__(self, other):
        if _is_hyper(other):
            return self._mul_hyper(other, 'rdiv')
        lifted = _lift_pair(other, "__rtruediv__")
        return NotImplemented if lifted is None else lifted._combine(self, 'div')

    def _mul_hyper(self, h, op: str):
        """Контейнер с Hyper: 'mul' — self·h, 'div' — self/h, 'rdiv' — h/self.
        Состояния — по правилам _combine (ноль и Центр 0 поглощают, класс 23
        липкий), остальное — в поле порядка, результат остаётся Hyper;
        обратный к многочлену — отказ Hyper.inverse (ряд Лорана)."""
        if self.is_zero_class:
            return PhiNumberV6.chaos() if self.zero_provenance else PhiNumberV6.void()
        if self.is_uncomputable:
            # род переносится; при 'rdiv' Центр 0 стоит ДЕЛИТЕЛЕМ — обращение
            # выводит из идеала бесконечно малых (_kind_after_inverse)
            kind = self.structural_kind if op != 'rdiv' else _kind_after_inverse(self.structural_kind)
            return PhiNumberV6.center_zero(structural=self.zero_provenance, kind=kind)
        if self.is_ignorance:
            a, b = self.approx_value(), h.approx()
            if op == 'mul':
                val = a * b
                params = _ieee_params(val, self, h, '*')
            elif op == 'div':
                val = a / b if b != 0 else float('nan')
                params = _ieee_params(val, self, h, '/')
            else:
                val = b / a if a != 0 else float('nan')
                params = _ieee_params(val, h, self, '/')
            return PhiNumberV6(LockedInt(0, True), LockedInt(1, True), LockedInt(0, True), LockedInt(1, True),
                               number_class=NumberClass.IGNORANCE, ignorance_kind=_ieee_kind(val),
                               base_params=params)
        if op == 'rdiv' and self.m1.v == 0 and not (self.is_infinity or self.is_infinitesimal):
            return PhiNumberV6.center_zero()          # x / 0 — как в _combine_ordinary
        s = self._to_hyper()
        if op == 'mul':
            return s * h
        # деление идёт через `/` поля порядка, а не через `h.inverse()`: у `inverse`
        # отказ называет один операнд («1/x»), и делимое терялось — замер 11.09
        # (сверка очереди с деревом, строка 39): `контейнер / Hyper.zero()` давал
        # «не сведено: 0 1/x» без пятёрки. У `/` отказ называет обе Ф-записи.
        if op in ('div', 'rdiv'):
            # отказ называет операнды КАК ДАНЫ (контейнер — своей записью, не после
            # подъёма в поле порядка): замер 11.09 — прежде здесь стояло
            # `h.inverse()`, и в отказе не было ни делимого, ни исходной записи
            try:
                return (s / h) if op == 'div' else (h / s)
            except Unresolved as e:
                left, right = (str(self), h.record()) if op == 'div' else (h.record(), str(self))
                raise Unresolved(getattr(e, "op", "/"), left, right,
                                 hint=getattr(e, "hint", None) or "x/0 здесь не сведено") from None
        return h * s.inverse()

    def _combine(self, other: "PhiNumberV6", op: str) -> "PhiNumberV6":
        if _is_hyper(other):
            return self._mul_hyper(other, op)
        if not isinstance(other, PhiNumberV6):
            lifted = _lift_pair(other, "__mul__" if op == 'mul' else "__truediv__")
            if lifted is None:
                return NotImplemented
            other = lifted
        return _settle(self._combine_core(other, op), self, other, '*' if op == 'mul' else '/')

    def _combine_core(self, other: "PhiNumberV6", op: str) -> "PhiNumberV6":

        # 1. ZERO — поглощает безусловно (даже над UNCOMPUTABLE/HYPER)
        if self.is_zero_class or other.is_zero_class:
            prov = (self.zero_provenance if self.is_zero_class else False) or \
                   (other.zero_provenance if other.is_zero_class else False)
            return PhiNumberV6.chaos() if prov else PhiNumberV6.void()

        # 2. UNCOMPUTABLE — поглощает (Center0 * x = Center0, x / Center0 = Center0)
        if self.is_uncomputable or other.is_uncomputable:
            prov = (self.zero_provenance if self.is_uncomputable else False) or \
                   (other.zero_provenance if other.is_uncomputable else False)
            # РОД неопределённости переносится (12.09): zero_state.py:34 —
            # «в умножении сохраняет свой род». Иначе уступающий Центр 0 от
            # ε−ε после умножения на ЕДИНИЦУ становился поглощающим
            # (kind=None читается как 'infinite'), и (ε−ε)·1 + 5 расходилось
            # с (ε−ε) + 5 — дыра устройства п. 7 главы 09. Делитель входит
            # обратным, и для него род поднимается: 2/(ε−ε) неограниченно.
            ks = _structural_kind(self)
            ko = _structural_kind(other)
            if op == 'div':
                ko = _kind_after_inverse(ko)
            return PhiNumberV6.center_zero(structural=prov, kind=_merge_kind(ks, ko))

        # 3. IGNORANCE — липкий/поглощающий
        if self.is_ignorance or other.is_ignorance:
            a, b = self.approx_value(), other.approx_value()
            val = a * b if op == 'mul' else (a / b if b != 0 else float('nan'))
            # 05.09: рядом с IEEE — точная сторона и потеря (ieee.py, слово автора)
            return PhiNumberV6(LockedInt(0, True), LockedInt(1, True), LockedInt(0, True), LockedInt(1, True),
                               number_class=NumberClass.IGNORANCE, ignorance_kind=_ieee_kind(val),
                               base_params=_ieee_params(val, self, other, '*' if op == 'mul' else '/'))

        # 4. HYPER (ω/ε) — особые правила
        if self.is_infinity or self.is_infinitesimal or other.is_infinity or other.is_infinitesimal:
            return self._combine_hyper(other, op)

        # 5. Обычные числа (включая экзотическое основание)
        return self._combine_ordinary(other, op)

    def _combine_hyper(self, other: "PhiNumberV6", op: str) -> "PhiNumberV6":
        s_inf, s_eps = self.is_infinity, self.is_infinitesimal
        o_inf, o_eps = other.is_infinity, other.is_infinitesimal

        def coeff_mul(a: "PhiNumberV6", b: "PhiNumberV6") -> Tuple[int, int]:
            return _reduce_ratio(a.m1.v * b.m1.v, a.m2.v * b.m2.v)

        def coeff_div(a: "PhiNumberV6", b: "PhiNumberV6") -> Tuple[int, int]:
            if b.m1.v == 0:
                raise ZeroDivisionError("HYPER: деление на нулевой коэффициент")
            return _reduce_ratio(a.m1.v * b.m2.v, a.m2.v * b.m1.v)

        s_h, o_h = (s_inf or s_eps), (o_inf or o_eps)
        if op == 'mul':
            if s_h and o_h:
                n, d = coeff_mul(self, other)
                return PhiNumberV6.from_hyper(n, d, self.hyper_order + other.hyper_order)
            # HYPER * ordinary_n (или наоборот): порядок указателя сохраняется
            hyper_side, ord_side = (self, other) if s_h else (other, self)
            if ord_side.m1.v == 0:
                return PhiNumberV6.from_int(0)  # 0 * HYPER = обычный незапертый 0
            n, d = coeff_mul(hyper_side, ord_side)
            return PhiNumberV6.from_hyper(n, d, hyper_side.hyper_order)
        else:  # div
            if s_h and o_h:
                n, d = coeff_div(self, other)
                return PhiNumberV6.from_hyper(n, d, self.hyper_order - other.hyper_order)
            if s_h:
                # HYPER / ordinary_n
                if other.m1.v == 0:
                    raise ZeroDivisionError("HYPER / 0")
                n, d = coeff_div(self, other)
                return PhiNumberV6.from_hyper(n, d, self.hyper_order)
            # ordinary_n / HYPER: n/ω^k → ε^k·(n/coeff); n/ε^k → ω^k·(n/coeff)
            if self.m1.v == 0:
                return PhiNumberV6.from_int(0)
            n, d = coeff_div(self, other)
            return PhiNumberV6.from_hyper(n, d, -other.hyper_order)

    def _combine_ordinary(self, other: "PhiNumberV6", op: str) -> "PhiNumberV6":
        if self.base_class is not None or other.base_class is not None:
            if self.base_class != other.base_class or self.base_params != other.base_params:
                # 12.09, п. 8 главы 09: отказ называет ОБЕ базы и ПУТЬ. Прежде
                # в сообщении стояли «self.base_class=1, other.base_class=1» —
                # числа РАВНЫ, а расхождение лежало в base_params ({'value': 3}
                # против {'value': 5}) и не печаталось; класс шёл числом без
                # имени. Отказ был, уточнить по нему было нельзя.
                raise CrossBaseOperationNotSupported(
                    f"{op}: разные основания α — путь _combine_ordinary(op={op!r}); "
                    f"слева α={_base_text(self)}, справа α={_base_text(other)}; "
                    f"расхождение: {_base_mismatch(self, other)}; покомпонентная "
                    f"операция потеряла бы точность; нужна явная трансляция (смена "
                    f"системы координат = трансляция, не вычисление)"
                )
            result_base_class = self.base_class
            result_base_params = dict(self.base_params)
        else:
            result_base_class = None
            result_base_params = {}

        # ---- Мантисса M ----
        a_locked, b_locked = self.m_pair_locked, other.m_pair_locked
        if a_locked:
            m1_res, m2_res = self.m1, self.m2
        elif b_locked:
            m1_res, m2_res = other.m1, other.m2
        else:
            if op == 'mul':
                num, den = self.m1.v * other.m1.v, self.m2.v * other.m2.v
            else:
                if other.m1.v == 0:
                    return PhiNumberV6.center_zero()
                num, den = self.m1.v * other.m2.v, self.m2.v * other.m1.v
            num, den = _reduce_ratio(num, den)
            m1_res, m2_res = LockedInt(num), LockedInt(den)

        # ---- Экспонента E ----
        ae_locked, be_locked = self.e_pair_locked, other.e_pair_locked
        if ae_locked:
            e1_res, e2_res = self.e1, self.e2
        elif be_locked:
            e1_res, e2_res = other.e1, other.e2
        else:
            if op == 'mul':
                e1n = self.e1.v * other.e2.v + other.e1.v * self.e2.v
            else:
                e1n = self.e1.v * other.e2.v - other.e1.v * self.e2.v
            e2n = self.e2.v * other.e2.v
            e1n, e2n = _reduce_ratio(e1n, e2n)
            e1_res, e2_res = LockedInt(e1n), LockedInt(e2n)

        result_class = self.number_class if self.number_class == other.number_class else NumberClass.NATURAL_POWER
        return PhiNumberV6(m1_res, m2_res, e1_res, e2_res,
                           number_class=result_class,
                           base_class=result_base_class, base_params=result_base_params)

    def inverse(self) -> "PhiNumberV6":
        """Структурная инверсия: swap(m1,m2) с переносом знака (m2 остаётся
        >0 — инвариант), негация e1. Для HYPER — дополнительно меняет класс
        (INFINITY↔INFINITESIMAL)."""
        if self.is_uncomputable:
            # род сохраняется как ФАКТ пути, но обращение выводит из идеала
            # бесконечно малых: 1/(ε−ε) уже не ограничено эпсилоном
            return PhiNumberV6.center_zero(structural=self.zero_provenance,
                                           kind=_kind_after_inverse(self.structural_kind))
        if self.is_zero_class:
            return self  # 1/0(void) не определено арифметически как обычное число — по спецификации ZERO не поглощается инверсией
        if self.is_ignorance:
            raw = self.base_params.get('raw_float')
            new_raw = (1.0 / raw) if raw not in (None, 0.0) else float('nan')
            from .hyper import Hyper
            from .rational import Q
            return PhiNumberV6(self.m1, self.m2, self.e1, self.e2,
                               number_class=NumberClass.IGNORANCE,
                               ignorance_kind=_ieee_kind(new_raw),
                               base_params=_ieee_params(new_raw, Hyper.const(Q(1)), self, '/'))
        if self.m1.v == 0:
            return PhiNumberV6.center_zero()
        sign = 1 if self.m1.v > 0 else -1
        new_m1 = LockedInt(sign * self.m2.v, self.m2.b)
        new_m2 = LockedInt(abs(self.m1.v), self.m1.b)
        new_e1 = LockedInt(-self.e1.v, self.e1.b)
        new_e2 = LockedInt(self.e2.v, self.e2.b)
        new_class = self.number_class
        if self.is_infinity or self.is_infinitesimal:
            new_e1 = LockedInt(self.e1.v, self.e1.b)          # порядок сохраняется, знак — классом
        if self.is_infinity:
            new_class = NumberClass.INFINITESIMAL
        elif self.is_infinitesimal:
            new_class = NumberClass.INFINITY
        return PhiNumberV6(new_m1, new_m2, new_e1, new_e2,
                           number_class=new_class,
                           base_class=self.base_class, base_params=dict(self.base_params))

    # ── Арифметика: степень ───────────────────────────────────────────
    #
    # Подсказка отказа называет КЛАСС показателя (строка 79 очереди, п. 2).
    # Замер 24.09 `[вычислено здесь]`: один и тот же текст «показатель —
    # величина поля порядка» стоял у показателей классов 1 (пара и контейнер),
    # 3, 5, 6, 10, 21 и 22 — верен он лишь для 21/22 (ω, ε) и Hyper; у пары
    # была объявлена дверь `.pow_rational`, а подсказка её не называла.
    # Класс читается из ЗАПИСИ показателя («Φ(C=<класс>…», протокол дома), а
    # не из типа: типов носителей контейнер не ввозит.

    def __pow__(self, power) -> "PhiNumberV6":
        if not isinstance(power, int) or isinstance(power, bool):
            # величина в показателе — не ошибка двери, а названный отказ: степень
            # контейнера в указателе (ω, ε) поля порядка не определена, и отказ
            # обязан нести обе Ф-записи (замер 11.09, строка 39 очереди: здесь
            # стоял голый TypeError без операндов)
            rec = getattr(power, "record", None)
            if callable(rec):
                text = rec()
                raise Unresolved("**", str(self), text, hint=_pow_hint(power, text))
            raise TypeError(f"PhiNumberV6 ** {power!r} ({type(power).__name__}): степень принимает int; "
                            "для рациональной степени p/q — .pow_rational(p, q)")
        return self.pow_rational(power, 1)

    def pow_rational(self, p: int, q: int = 1) -> "PhiNumberV6":
        if q <= 0:
            raise ValueError(f"q должно быть > 0, получено {q}")
        if self.is_uncomputable:
            # натуральная степень идеал не покидает (I·I ⊆ I), отрицательная
            # покидает через обращение — род поднимается
            kind = self.structural_kind if p >= 0 else _kind_after_inverse(self.structural_kind)
            return PhiNumberV6.center_zero(structural=self.zero_provenance, kind=kind)
        if self.is_zero_class:
            return self  # ZERO ** n = ZERO (пустота не становится единицей даже при n=0)
        if self.is_ignorance:
            # Класс 23 — навык: 23 на входе, 23 на выходе; степень считает
            # сам IEEE 754 над битово точным raw_float, вид (nan/inf/finite)
            # определяется заново по итогу. Прежде ветка отсутствовала и
            # степень шла общим путём мантиссы, теряя провенанс
            # (IGNORANCE:inf ** 1 → IGNORANCE:None → nan); поправка 03.09
            # по слову автора «выгода и цена — расчётом».
            raw = self.base_params.get('raw_float', float('nan'))
            try:
                val = raw ** (p / q) if q != 1 else raw ** p
            except (OverflowError, ZeroDivisionError):
                val = float('inf') if raw != 0 else float('nan')
            if isinstance(val, complex):
                val = float('nan')
            val = float(val)
            from .ieee import loss
            from .hyper import Hyper
            ex = self.base_params.get('exact')
            if q == 1 and isinstance(ex, Hyper):
                try:
                    ex = ex ** p
                except Exception:
                    ex = None
            elif q != 1:
                ex = None                      # дробная степень: точная сторона не сведена
            return PhiNumberV6(LockedInt(0, True), LockedInt(1, True), LockedInt(0, True), LockedInt(1, True),
                               number_class=NumberClass.IGNORANCE, ignorance_kind=_ieee_kind(val),
                               base_params={'raw_float': val, 'exact': ex, 'loss': loss(ex, val)})
        if p == 0:
            return PhiNumberV6.from_int(1)
        if self.is_infinity or self.is_infinitesimal:
            g = math.gcd(abs(p), q)
            p, q = p // g, q // g
            if q != 1:
                raise ExactRepresentationNotAvailable("pow_rational: дробная степень HYPER-числа не поддержана")
            if p < 0:
                return self.inverse().pow_rational(-p, 1)
            n, d = _reduce_ratio(self.m1.v ** p, self.m2.v ** p)
            return PhiNumberV6.from_hyper(n, d, self.hyper_order * p)
        if self.base_class is not None:
            raise ExactRepresentationNotAvailable(
                "pow_rational: точная степень для base_class != None (α≠2) не реализована — "
                "используйте .approx_value()"
            )

        g_pq = math.gcd(abs(p), q)
        p, q = p // g_pq, q // g_pq

        if q == 1:
            m1n = self.m1.v ** abs(p)
            m2n = self.m2.v ** abs(p)
            if p < 0:
                if m1n == 0:
                    raise ZeroDivisionError("pow_rational: 0 в отрицательной степени")
                sign = 1 if m1n > 0 else -1
                m1n, m2n = sign * m2n, abs(m1n)
            e1n = self.e1.v * p
            e2n = self.e2.v
            return PhiNumberV6.from_mantissa_exponent(m1n, m2n, e1n, e2n, self.number_class)

        m1_pow = self.m1.v ** abs(p)
        m2_pow = self.m2.v ** abs(p)
        if p < 0:
            if m1_pow == 0:
                raise ZeroDivisionError("pow_rational: 0 в отрицательной степени")
            sign = 1 if m1_pow > 0 else -1
            m1_pow, m2_pow = sign * m2_pow, abs(m1_pow)
        sign_final = 1 if m1_pow >= 0 else -1
        # Корень берётся с НЕЧЁТНОЙ части мантиссы, а не с мантиссы целиком
        # (расхождение D1 записки «стены и известные законы», 09.10). Мантисса
        # несёт множитель базы: |m1_pow| = 2^k1·u1, m2_pow = 2^k2·u2 с нечётными
        # u1, u2 — единственно [внешний источник: Hardy–Wright, теорема 2].
        # Тогда (m1_pow/m2_pow)^(1/q) = (u1/u2)^(1/q) · 2^((k1−k2)/q): в форме
        # (m1/m2)·2^(e1/e2) корень точен тогда и только тогда, когда u1 и u2 —
        # точные q-е степени; множитель 2^((k1−k2)/q) ложится в показатель,
        # который по форме контейнера дробным быть вправе (инвариант класса —
        # лишь e2 > 0; эта же ветвь и прежде строила e2·q). До правки
        # [вычислено здесь, 09.10]: 2^(1/2) и 18^(1/2) — отказ при точном ответе
        # 1·2^(1/2) и 3·2^(1/2); после: оба сведены, 12^(1/2) — отказ остаётся
        # (нечётная часть 3 не квадрат). Целая часть (k1−k2) div q
        # возвращается в мантиссу: там, где мантисса была точной q-й степенью
        # (4^(1/2) = 2·2^0, 8^(1/3) = 2·2^0), запись не меняется — в показатель
        # уходит лишь остаток b/q, 0 ≤ b < q.
        k1, u1 = _two_adic(abs(m1_pow))
        k2, u2 = _two_adic(m2_pow)
        num_ok, num_root = _perfect_nth_root(q, u1)
        den_ok, den_root = _perfect_nth_root(q, u2)
        if not (num_ok and den_ok):
            raise ExactRepresentationNotAvailable(
                f"pow_rational({p}/{q}): мантисса {self.m1.v}/{self.m2.v} в степени {p}/{q} "
                f"не имеет точного представления (не точный q-й корень); используйте .approx_value()"
            )
        a, b = divmod(k1 - k2, q)
        new_m1 = sign_final * num_root
        new_m2 = den_root
        if a >= 0:
            new_m1 *= 2 ** a
        else:
            new_m2 *= 2 ** (-a)
        # показатель: e1·p/(e2·q) + b/q = (e1·p + b·e2) / (e2·q)
        e1n = self.e1.v * p + b * self.e2.v
        e2n = self.e2.v * q
        return PhiNumberV6.from_mantissa_exponent(new_m1, new_m2, e1n, e2n, self.number_class)

    # ── Арифметика: сложение / вычитание ──────────────────────────────
    #
    # __sub__ реализован через __add__(-other) — это позволяет иметь ОДНУ
    # ветвь приоритетной логики (_add_impl) без дублирования знаковых
    # случаев. Приоритет: ZERO (уступает) > UNCOMPUTABLE (участвует как 0,
    # НЕ поглощает) > IGNORANCE (липкий) > HYPER (особые правила,
    # включая критичный фикс ε+0=ε) > обычные числа (common-exponent
    # alignment).

    def __add__(self, other: "PhiNumberV6"):
        if _is_hyper(other):
            return self._add_hyper(other)
        if not isinstance(other, PhiNumberV6):
            lifted = _lift_pair(other, "__add__ (она же __radd__)")   # пара дома — дверью (12.09)
            if lifted is None:
                return NotImplemented
            other = lifted
        return self._add_impl(other)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other: "PhiNumberV6"):
        if _is_hyper(other):
            return self._add_hyper(-other)
        if not isinstance(other, PhiNumberV6):
            lifted = _lift_pair(other, "__sub__")
            if lifted is None:
                return NotImplemented
            other = lifted
        return self._add_impl(-other)

    def __rsub__(self, other):
        if _is_hyper(other):
            return (-self)._add_hyper(other)
        lifted = _lift_pair(other, "__rsub__")
        return NotImplemented if lifted is None else lifted._add_impl(-self)

    # ── Мост в поле порядка (05.09, очередь п. 39(c); свежие глаза, глава 03) ──
    #
    # Два чтения указателей 21/22 (00-ЛИНИЯ §2.6, A9 — класс против элемента):
    #   контейнер — чтение указателя (консенсус §2.2 «указатели масштаба, не
    #   числа»): ω + 5 = ω, 5 + ε = 5, ω − ω — ноль-состояние, поглощающее
    #   (слово автора, zero_state №2); phicore.hyper.Hyper — чтение элемента:
    #   кольцо многочленов Лорана над ℚ, ω + 5 = 1·ω + 5, ω − ω = 0{ω−ω}.
    # Замер 05.09 (свип 1728 троек, proto/test_phicore_container_order.py):
    #   чтение указателя — 16 нарушений (+), чтение элемента — 4 (все класс 23,
    #   IEEE inf − inf), смешение (сумма в кольце, ноль поглощает) — 24.
    #   Смешение хуже обоих чистых — на контейнере чтение не смешивается.
    # Мост: контейнер ⊕ Hyper считается в чтении элемента и ОСТАЁТСЯ Hyper
    #   (чтение липкое, как класс 23); состояния контейнера решают сами
    #   (ноль уступает, Центр 0 рода ω±ω поглощает, класс 23 липкий).
    #   Конечное без точного отношения (α ≠ 2, дробная экспонента) в поле
    #   порядка над ℚ не поднимается — отказ поимённо, не поглощение.

    def _to_hyper(self):
        """Подъём: указатель → одночлен своего порядка, обычное число →
        константа по точному отношению; иначе ExactRepresentationNotAvailable."""
        from .hyper import Hyper
        from .rational import Q
        if self.is_infinity or self.is_infinitesimal:
            return Hyper({self.hyper_order: Q(self.m1.v, self.m2.v)})
        try:
            num, den = self.to_exact_int_ratio()
        except ExactRepresentationNotAvailable as ex:
            raise ExactRepresentationNotAvailable(
                f"контейнер ⊕ Hyper считается в поле порядка над ℚ; "
                f"{self!r} точного отношения не имеет — отказ поимённо, не поглощение: {ex}") from ex
        return Hyper.const(Q(num, den))

    def _add_hyper(self, h):
        """Контейнер + Hyper: состояния решают сами (ноль уступает, Центр 0
        рода ω±ω поглощает, класс 23 липкий), остальное — в поле порядка,
        результат остаётся Hyper (чтение элемента липкое)."""
        if self.is_zero_class:
            return h
        if self.is_uncomputable:
            return self if (self.zero_provenance and self.structural_kind != 'infinitesimal') else h
        if self.is_ignorance:
            val = self.approx_value() + h.approx()
            return PhiNumberV6(LockedInt(0, True), LockedInt(1, True), LockedInt(0, True), LockedInt(1, True),
                               number_class=NumberClass.IGNORANCE, ignorance_kind=_ieee_kind(val),
                               base_params=_ieee_params(val, self, h, '+'))
        return self._to_hyper() + h

    def _add_impl(self, other: "PhiNumberV6") -> "PhiNumberV6":
        return _settle(self._add_core(other), self, other, '+')

    def _add_core(self, other: "PhiNumberV6") -> "PhiNumberV6":
        # 1. ZERO — уступает место бытию
        self_zero, other_zero = self.is_zero_class, other.is_zero_class
        if self_zero and other_zero:
            prov = self.zero_provenance or other.zero_provenance
            return _carry_absorbed(PhiNumberV6.chaos() if prov else PhiNumberV6.void(), self, other)
        if self_zero:
            return other
        if other_zero:
            return self

        # 2. UNCOMPUTABLE — участвует КАК 0 (НЕ поглощает) для АБСОЛЮТНОГО
        #    Center0 (напр. Ω Хайтина) И для СТРУКТУРНОГО с kind='infinitesimal'
        #    (ε±ε — результат ГАРАНТИРОВАННО бесконечно мал или нулевой,
        #    бесконечно малые образуют идеал, замкнутый относительно +/- —
        #    значит "неизвестно, но заведомо не больше эпсилон" — уступает
        #    конечному/бесконечному числу так же, как обычный ε).
        #    ПОГЛОЩАЕТ только structural с kind='infinite' (ω±ω — результат
        #    МОГ БЫ быть любым: конечным/бесконечным/бесконечно малым —
        #    настоящая неограниченная индетерминация, аналог ∞-∞).
        #    Найдено систематическим тестом test_associativity_sweep.py,
        #    2026-07-26 — ДВА прохода: первый фикс ("любой structural
        #    поглощает") починил (ω+(-ω))+5, но СЛОМАЛ (ε+(-ε))+5 (раньше
        #    корректно дававшее 5) — сумма/разность двух бесконечно малых
        #    не может стать неограниченной индетерминацией, это доказанный
        #    факт (замкнутость идеала бесконечно малых), поэтому потребовалось
        #    различать РОД структурной неопределённости, а не только сам факт
        #    её структурности.
        def _absorbs(x: "PhiNumberV6") -> bool:
            return x.zero_provenance and x.structural_kind != 'infinitesimal'

        # род при встрече двух Центров 0 — общим правилом старшинства
        # (_merge_kind, выше в модуле): с 12.09 тем же, что у умножения,
        # деления, степени и обращения — прежде правило жило только здесь
        def _combine_kind(a: "PhiNumberV6", b: "PhiNumberV6") -> Optional[str]:
            return _merge_kind(_structural_kind(a), _structural_kind(b))

        self_unc, other_unc = self.is_uncomputable, other.is_uncomputable
        if self_unc and other_unc:
            prov = self.zero_provenance or other.zero_provenance
            return _carry_absorbed(PhiNumberV6.center_zero(structural=prov, kind=_combine_kind(self, other)), self, other)
        if self_unc:
            if _absorbs(self):
                return _absorb(self, other)          # поглощено — записано (канцеллер указателя)
            return other
        if other_unc:
            if _absorbs(other):
                return _absorb(other, self)
            return self

        # 3. IGNORANCE — липкий/поглощающий
        if self.is_ignorance or other.is_ignorance:
            val = self.approx_value() + other.approx_value()
            # 05.09: рядом с IEEE — точная сторона и потеря (ieee.py, слово автора)
            return PhiNumberV6(LockedInt(0, True), LockedInt(1, True), LockedInt(0, True), LockedInt(1, True),
                               number_class=NumberClass.IGNORANCE, ignorance_kind=_ieee_kind(val),
                               base_params=_ieee_params(val, self, other, '+'))

        # 4. HYPER (ω/ε)
        s_inf, s_eps = self.is_infinity, self.is_infinitesimal
        o_inf, o_eps = other.is_infinity, other.is_infinitesimal
        if s_inf or s_eps or o_inf or o_eps:
            if (s_inf and o_inf) or (s_eps and o_eps):
                if self.hyper_order != other.hyper_order:
                    # старший порядок доминирует; младший — поглощён и записан (канцеллер указателя, 05.09)
                    return _absorb(self, other) if self.hyper_order > other.hyper_order else _absorb(other, self)
                n = self.m1.v * other.m2.v + other.m1.v * self.m2.v
                d = self.m2.v * other.m2.v
                if n == 0:
                    return _carry_absorbed(PhiNumberV6.center_zero(structural=True, kind='infinite' if s_inf else 'infinitesimal'), self, other)
                n, d = _reduce_ratio(n, d)
                return _carry_absorbed(PhiNumberV6.from_hyper(n, d, self.hyper_order), self, other)
            if s_inf or o_inf:
                # ∞ доминирует над ε и над обычным числом — поглощённое записано
                return _absorb(self, other) if s_inf else _absorb(other, self)
            # ровно один операнд — INFINITESIMAL, другой — обычное число
            eps_side, ord_side = (self, other) if s_eps else (other, self)
            if ord_side.m1.v == 0:
                return eps_side  # КРИТИЧНЫЙ ФИКС: ε + 0 = ε (не 0!)
            return _absorb(ord_side, eps_side)  # ε пренебрежимо мало относительно ненулевого числа — записано
            # чтение указателя (консенсус 00-ЛИНИЯ §2.2: «указатели масштаба, не числа»);
            # чтение элемента (ω + 5 = 1·ω + 5) — в phicore.hyper, мост _add_hyper; смешение
            # чтений на контейнере замерено 05.09: 24 нарушения против 20 (указатель) и 4 (элемент)

        # 5. Обычные числа / экзотическое основание
        return self._add_ordinary(other)

    def _add_ordinary(self, other: "PhiNumberV6") -> "PhiNumberV6":
        any_locked = (self.m1.b or self.m2.b or self.e1.b or self.e2.b or
                      other.m1.b or other.m2.b or other.e1.b or other.e2.b)
        if any_locked:
            # "результат блокируется (точное значение не определено)":
            # приоритет self, если у self ЕСТЬ хоть один установленный бит;
            # иначе (лок только у other) — приоритет other.
            src = self if (self.m1.b or self.m2.b or self.e1.b or self.e2.b) else other
            return PhiNumberV6(
                LockedInt(src.m1.v, True), LockedInt(src.m2.v, True),
                LockedInt(src.e1.v, True), LockedInt(src.e2.v, True),
                number_class=src.number_class, base_class=src.base_class,
                base_params=dict(src.base_params)
            )

        if self.base_class is not None or other.base_class is not None:
            if self.base_class != other.base_class or self.base_params != other.base_params:
                # 12.09, п. 8 главы 09: здесь сообщение не содержало вовсе
                # ничего о том, какие базы столкнулись, — теперь обе базы и путь
                raise CrossBaseOperationNotSupported(
                    f"+: разные основания α — путь _add_ordinary(op='+'); "
                    f"слева α={_base_text(self)}, справа α={_base_text(other)}; "
                    f"расхождение: {_base_mismatch(self, other)}; сложение чисел с "
                    f"разными основаниями α не имеет общей точной формы в этом "
                    f"контейнере (нужна явная трансляция либо символьная сумма термов)"
                )
            # Для экзотического основания трюк "сдвиг мантиссы на α^diff" РАБОТАЕТ
            # только когда diff=0 (совпадающие показатели) — для diff≠0 α^diff
            # иррационально (α трансцендентно/алгебраически "неудобно"), в отличие
            # от базы 2, где 2^целое всегда целое. Поэтому здесь — ТОЛЬКО точное
            # совпадение показателей, иначе явное исключение.
            if self.e1.v * other.e2.v != other.e1.v * self.e2.v:
                raise ExactRepresentationNotAvailable(
                    "сложение с экзотическим основанием при РАЗНЫХ показателях степени "
                    "не имеет представления вида (m1/m2)·α^(e1/e2) в этом контейнере; "
                    "нужна символьная сумма термов (см. InvariantTerm), не этот контейнер"
                )
            new_m_num = self.m1.v * other.m2.v + other.m1.v * self.m2.v
            new_m_den = self.m2.v * other.m2.v
            new_m_num, new_m_den = _reduce_ratio(new_m_num, new_m_den)
            return PhiNumberV6.exotic_base(new_m_num, new_m_den, self.e1.v, self.e2.v,
                                           self.base_class, self.base_params, self.number_class)

        # База 2: common-exponent alignment. diff = E_self - E_other. Точное
        # сложение возможно, если diff ЦЕЛОЕ (не обязательно ноль!) — тогда
        # операнд с БОЛЬШЕЙ экспонентой домножается на 2^|diff|, результат
        # хранит МЕНЬШУЮ экспоненту (см. проверку на sqrt(2)+8·sqrt(2)=9·sqrt(2)
        # в test_phi_number_v6.py).
        # 12.09: РАВНЫЙ дробный показатель (diff = 0/d) здесь отказывал вопреки
        # этому докстрингу — причина была не тут, а в стороже нулевого
        # числителя в _reduce_ratio (0/4 не становилось 0/1). Замер свипом
        # 12 пар: сведено 4 → 12, отказ 8 → 0.
        diff_num = self.e1.v * other.e2.v - other.e1.v * self.e2.v
        diff_den = self.e2.v * other.e2.v
        diff_num, diff_den = _reduce_ratio(diff_num, diff_den)
        if diff_den != 1:
            raise ExactRepresentationNotAvailable(
                f"сложение: показатели степени различаются на нецелую величину "
                f"{diff_num}/{diff_den} — точного представления в базе 2 не существует "
                f"(используйте .approx_value() на каждом слагаемом, либо переход к "
                f"символьному уровню)"
            )
        diff = diff_num  # целое: E_self - E_other

        if diff > 0:
            # self.E больше — сдвигаем мантиссу self, результат хранит other.E (меньшую)
            shifted_num = self.m1.v * (2 ** diff)
            shifted_den = self.m2.v
            new_num = shifted_num * other.m2.v + other.m1.v * shifted_den
            new_den = shifted_den * other.m2.v
            result_e1, result_e2 = other.e1.v, other.e2.v
        elif diff < 0:
            # other.E больше — сдвигаем мантиссу other, результат хранит self.E (меньшую)
            shifted_num = other.m1.v * (2 ** (-diff))
            shifted_den = other.m2.v
            new_num = self.m1.v * shifted_den + shifted_num * self.m2.v
            new_den = self.m2.v * shifted_den
            result_e1, result_e2 = self.e1.v, self.e2.v
        else:
            new_num = self.m1.v * other.m2.v + other.m1.v * self.m2.v
            new_den = self.m2.v * other.m2.v
            result_e1, result_e2 = self.e1.v, self.e2.v

        new_num, new_den = _reduce_ratio(new_num, new_den)
        if new_num == 0:
            return PhiNumberV6.from_int(0)
        return PhiNumberV6.from_mantissa_exponent(new_num, new_den, result_e1, result_e2,
                                                  NumberClass.NATURAL_POWER)

    # ── Самопроверка (ось «AIAΣA спрашивает AIAΣA», 04.09) ─────────────

    def self_check(self, deep: bool = False) -> list:
        """Самопроверка записи ВСТРЕЧНОЙ дорогой, а не тем же счётом.

        `[вычислено здесь]` Дороги: (1) инварианты записи (m2 > 0, e2 > 0;
        классы ZERO и UNCOMPUTABLE требуют m1 = 0) — те, что ставит
        `__post_init__`, здесь пересчитываются на УЖЕ собранном объекте, в
        том числе на собранном в обход конструктора; (2) канон сокращения —
        НОД пары считается тут заново (`math.gcd`), а не доверяется
        `_reduce_ratio`; (3) величина через СОСЕДНИЙ модуль: старшая запись
        `record.PhiRecord.from_container(self).to_fraction()` (у неё знак —
        отдельное поле s, а не знак m1) обязана дать ту же пару, что
        `to_exact_int_ratio()`; (4) контракт eq/hash — пересборка тех же
        полей даёт равный объект и равный ключ; (5) запись и печать — один
        текст (`__repr__ = record`), и класс в ней назван.

        Граница объявлена: точной пары нет у состояний (ZERO, Центр 0,
        указатели 21/22, класс 23), при α ≠ 2 и при дробном показателе —
        там дорога (3) не идёт вовсе, и её молчание находкой НЕ считается;
        мост через старшую запись отказывает на тех же входах поимённо.
        Вторая граница названа самой этой самопроверкой (12.09): старшая
        запись несёт M ПАРОЙ Q и битов блокировки не хранит, поэтому
        заблокированная мантисса 4🔒/2🔒 приходит оттуда сведённой к 2/1.
        Сравнение поэтому идёт по ВЕЛИЧИНЕ (пара против Q), а не по паре
        целых: блокировка есть метаданное записи, в значение не входит
        (канон «lock is metadata, not value»), и расхождение записей здесь
        законно, расхождение величин — находка.
        `deep` добавляет дорогу канцеллера: полная величина `_full_value`
        (удержанное плюс поглощённое) обязана совпасть с подъёмом самой
        записи, когда поглощённого нет. Возвращает список находок; пустой
        список — молчание, исключений не бросает.
        """
        out = []
        for name, comp in (("m1", self.m1), ("m2", self.m2),
                           ("e1", self.e1), ("e2", self.e2)):
            if not isinstance(comp, LockedInt) or not isinstance(comp.v, int) \
                    or isinstance(comp.v, bool):
                out.append(f"компонент {name} не целое в LockedInt: {comp!r}")
        if out:
            return out
        if self.m2.v <= 0:
            out.append(f"m2.v не положителен: {self.m2.v} (знак — в m1)")
        if self.e2.v <= 0:
            out.append(f"e2.v не положителен: {self.e2.v} (каноническая форма)")
        if self.number_class in (NumberClass.ZERO, NumberClass.UNCOMPUTABLE) \
                and self.m1.v != 0:
            out.append(f"класс {self.number_class.name} требует m1.v == 0, "
                       f"записано {self.m1.v}")
        if out:
            return out
        if not (self.m1.b or self.m2.b):
            g = math.gcd(abs(self.m1.v), self.m2.v)
            if g > 1 and not (self.is_infinity or self.is_infinitesimal):
                out.append(f"мантисса не сокращена: gcd({self.m1.v}, "
                           f"{self.m2.v}) = {g}")
        if not (self.e1.b or self.e2.b):
            g = math.gcd(abs(self.e1.v), self.e2.v)
            if g > 1:
                out.append(f"показатель не сокращён: gcd({self.e1.v}, "
                           f"{self.e2.v}) = {g}")
        if self.structural_kind is not None and not (
                self.is_uncomputable and self.zero_provenance):
            out.append(f"род {self.structural_kind!r} стоит не на структурном "
                       f"Центре 0, а на классе {self.number_class.name}")
        # встречная дорога: величина через старшую запись соседнего модуля
        try:
            pair = self.to_exact_int_ratio()
        except ExactRepresentationNotAvailable:
            pair = None                       # объявленная граница, не находка
        if pair is not None:
            try:
                from .record import PhiRecord
                from .rational import Q as _Q
                f = PhiRecord.from_container(self).to_fraction()
                if f != _Q(pair[0], pair[1]):
                    out.append(f"величина двумя дорогами расходится: точная "
                               f"пара {pair}, старшая запись "
                               f"{(f.numerator, f.denominator)}")
            except ExactRepresentationNotAvailable:
                if self.m1.v != 0:
                    out.append("точная пара есть, а старшая запись Φ(s, α, M, E) "
                               "её не принимает")
            except Exception as e:                              # noqa: BLE001
                out.append(f"мост в старшую запись отказал не по границе: "
                           f"{type(e).__name__}: {e}")
        twin = PhiNumberV6(self.m1, self.m2, self.e1, self.e2,
                           number_class=self.number_class,
                           base_class=self.base_class,
                           base_params=dict(self.base_params),
                           zero_provenance=self.zero_provenance,
                           ignorance_kind=self.ignorance_kind,
                           structural_kind=self.structural_kind)
        if twin != self:
            out.append("пересборка тех же полей не равна себе")
        elif hash(twin) != hash(self):
            out.append("хэш не согласован с равенством: у равных записей "
                       "разные ключи")
        try:
            if self.record() != repr(self):
                out.append("запись и печать разошлись (канон: __repr__ = record)")
            if f"C={int(self.number_class)}" not in self.record() \
                    and self.number_class not in (NumberClass.ZERO,
                                                  NumberClass.UNCOMPUTABLE):
                out.append(f"класс не назван в записи: {self.record()}")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"запись отказала: {type(e).__name__}: {e}")
        if not deep:
            return out
        if 'absorbed' not in self.base_params and pair is not None:
            full = _full_value(self)
            lift = _ieee_lift(self)
            if full is not None and lift is not None and full != lift:
                out.append("канцеллер: полная величина без поглощённого не "
                           "совпала с подъёмом записи")
        return out

    # ── Отображение ────────────────────────────────────────────────────

    def _absorbed_record(self) -> str:
        """Поглощённое канцеллером указателя — хвостом записи; «не сведено»,
        когда поглощённая величина в поле порядка не поднялась."""
        if 'absorbed' not in self.base_params:
            return ""
        a = self.base_params['absorbed']
        return "; поглощено=" + ("не сведено" if a is None else a.describe())

    def _exact_side_record(self) -> str:
        """Точная сторона класса 23 записью (форма ieee.Ieee.exact_record)."""
        ex = self.base_params.get('exact')
        if ex is None:
            return "не сведено"
        return ex.describe() if hasattr(ex, "describe") else ex.record()

    def _exact_pair_record(self) -> str:
        """Точная пара целых — хвостом записи, РОВНО когда она существует
        (неявная база 2 и целый показатель). Это замена приближения из
        п. 6 главы 09: печаталось «≈ 1.5» через approx_value(), тогда как
        трансляция наружу объявлена явной дверью. Пара не печатается при
        нулевом показателе (она и есть мантисса) и когда её нет — отсутствие
        точной пары видно по самой записи показателя (2^(1/2))."""
        if self.e1.v == 0 or self.base_class is not None:
            return ""
        try:
            n, d = self.to_exact_int_ratio()
        except ExactRepresentationNotAvailable:
            return ""
        return f", точная пара={n}/{d}"

    def record(self) -> str:
        """Явная Ф-запись контейнера — протокол дома (п. 5 главы 09).

        До 12.09 метода не было (`hasattr(from_int(3), 'record')` ложно), и
        следствий было два, оба замерены главой: фасад ФIOS печатал
        контейнер ЗАПАСНОЙ ветвью `Φ(?=…)` (`dispatch.py:22`), затаскивая в
        текст Ф-отказа приближение из печати, а `from_ratio(3,2) **
        from_int(2)` давал голый TypeError вместо названного отказа, который
        обещает комментарий `__pow__` (у носителя с `record` — PI,
        Hyper.omega(1) — отказ именно названный: сторож зовёт `record`
        показателя).

        ФОРМА — СВОЯ, КАК У СОСЕДЕЙ, И ОДНА НА ПРЕДМЕТ. У дома два вида
        записи: носители одного класса печатают `Φ(C=<класс>, …)`
        (`zero_state.Zero`, `pointers.ScalePointer`, `ieee.Ieee`,
        `rational.Q`), а линейки со своим устройством — своим старшим
        знаком: `Φα(M=…, α=…, E=…)` у оси основания (`alpha.py`),
        `⟦lo; hi⟧·2^−g` у засечки (`gridcore.Notch`), `log_α(1 + k·α^d)` у
        логарифма сложения (`zech.py`) — и у ВСЕХ ТРЁХ `__repr__ = record`:
        одна запись, а не две. Контейнер держит все 24 класса и четыре
        компонента с битами блокировки, поэтому его старший знак — `Φv6(…)`,
        а класс назван внутри числом и именем, как у носителей: `Φv6(3/2 ·
        2^(0/1), C=1 NATURAL_POWER)`. Печать и запись — ОДИН текст
        (`__repr__ = record`), потому что отказы дома требуют именно этого:
        приёмка `proto/test_hyper_bridge.py` мерит, что в отказе стоит
        `record()` операнда (замер 12.09: пока record() и repr расходились,
        78 исходов зоны контейнера читались как «отказ без имени операнда»,
        хотя обе записи были верны, — расхождение форм и есть разнобой
        п. 12 главы 09 в чистом виде).
        Float не печатается никогда (п. 6): вместо приближения — точная
        пара, когда она существует. Работает на всех 24 классах числа и на
        всех 24 классах основания — запись не зовёт `_resolve_base_float`.
        Запись НЕ есть равенство: две разные записи одной величины дают два
        разных `record()` (докстринг класса о структурном равенстве).
        12.09 закрыта и косметика §21: род структурного Центра 0 печатается
        (`structural_kind` решает арифметику, а в записи его не было)."""
        ab = self._absorbed_record()
        if self.is_void:
            return f"Φv6(∅ nihil absolutum, C=0 ZERO{ab})"
        if self.is_chaos_ordinabile:
            return f"Φv6(∅ chaos ordinabile, C=0 ZERO{ab})"
        if self.is_absolute_zero:
            return f"Φv6(Center-0 absolute, C=20 UNCOMPUTABLE{ab})"
        if self.is_structural_zero:
            kind = "род не уточнён" if self.structural_kind is None else self.structural_kind
            return f"Φv6(Center-0 structural[{kind}], C=20 UNCOMPUTABLE{ab})"
        if self.is_ignorance:
            # форма ieee.Ieee.record: образ IEEE стоит в записи САМ — он
            # хранится битово и не есть приближение вычисления дома
            return (f"Φv6(C=23 IGNORANCE raw={self.base_params.get('raw_float')!r}; "
                    f"точно={self._exact_side_record()}; "
                    f"потеря={self.base_params.get('loss')})")
        if self.is_infinity:
            k = "" if self.e1.v == 1 else f"^{self.e1.v}"
            return f"Φv6(ω{k}·{self.m1.v}/{self.m2.v}, C=22 INFINITY{ab})"
        if self.is_infinitesimal:
            k = "" if self.e1.v == 1 else f"^{self.e1.v}"
            return f"Φv6(ε{k}·{self.m1.v}/{self.m2.v}, C=21 INFINITESIMAL{ab})"
        base_str = "2" if self.base_class is None else f"α{_base_text(self)}"
        return (f"Φv6({self.m1!r}/{self.m2!r} · {base_str}^({self.e1!r}/{self.e2!r})"
                f"{self._exact_pair_record()}, "
                f"C={int(self.number_class)} {self.number_class.name}{ab})")

    #: Печать и запись — один текст (канон соседей: Φα, Notch, Zech, Ieee).
    #: Прежняя печать несла float («≈ 1.5» через approx_value) и отказывала
    #: на 19 классах основания из 24 — п. 6 главы 09; теперь основание идёт
    #: именем класса и параметрами, значение — точной парой, когда она есть.
    __repr__ = record


if __name__ == "__main__":
    a = PhiNumberV6.from_ratio(3, 2)
    b = PhiNumberV6.from_ratio(4, 1)
    print("a =", a)
    print("b =", b)
    print("a*b =", a * b)
    print("a/b =", a / b)
    print("inverse(a) =", a.inverse())
    print("a**2 =", a ** 2)
    print("a + b =", a + b)
    root2 = PhiNumberV6.from_mantissa_exponent(1, 1, 1, 2)
    print("sqrt(2) =", root2, "approx=", root2.approx_value())
    w = PhiNumberV6.omega(1, 1)
    e = PhiNumberV6.epsilon(1, 1)
    print("w*e =", w * e, "| w+5 =", w + PhiNumberV6.from_int(5), "| e+0 =", e + PhiNumberV6.from_int(0))
