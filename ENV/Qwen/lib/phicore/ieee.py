"""Точная сторона IEEE 754: что теряет float, где и на сколько — числом.

Слово автора 05.09: «мы точно знаем, что IEEE теряет, в каком месте,
численно — предсказуемо точно; целочисленная арифметика без акселератора;
IEEE-коррекция — в приоритет». Основание — три факта, все точные:

1. Всякий конечный float есть двоичная дробь m·2^e (Q) — битово точно
   (`float.as_integer_ratio()`, целые без потерь). Класс 23 хранит
   raw_float битово; здесь он же читается как точная величина.
2. Потеря одной операции IEEE — точная величина: loss = exact(a ∘ b) −
   Q(fl(a ∘ b)), считается целыми (числители и степени двойки), без
   акселератора. Классическая дорога к тому же числу — безошибочные
   преобразования (Møller 1965, Knuth 1969, TwoSum; Dekker 1971,
   TwoProduct): ошибка суммы/произведения двух float сама есть float. Обе
   дороги обязаны совпадать — приёмка сверяет.
3. Бесконечности и nan — не значения, а указатели и состояния:
   ±inf ↔ ±ω (указатель масштаба, поле порядка), nan на входе ↔ состояние
   NOT_COMPUTED («путь потерян»: nan не помнит, из чего возник);
   inf − inf ↔ ω − ω (структурный ноль-состояние), 0·inf ↔ 0·ω = 0,
   inf/inf ↔ ω/ω = 1. Точная сторона считается в чтении элемента
   (phicore.hyper: кольцо многочленов Лорана над ℚ) — там ω + 0.1 − ω = 0.1
   точно, тогда как IEEE даёт nan: потеря названа, не проглочена.

`Ieee` — число класса 23 с точной стороной: raw (float, как посчитал бы
IEEE: «на входе 23 — на выходе 23, не мешать»), exact (Hyper | Zero),
loss (накопленная потеря: exact − Q(raw), либо именованная —
«переполнение», «nan: путь», «не сведено»). Ни один результат не
подменяется: IEEE остаётся IEEE, рядом лежит точная величина и разность.

Род нуля и субнормальный край замерены 11.09 `[вычислено здесь]`. Знак нуля —
ровно тот род, который теряет модуль (регистр §512, «скалярная амнезия»): для
величины его нет (−0 = +0 = 0 точно), для IEEE он есть и решает знак
бесконечности при делении. Потому он не подменяется и не выдумывается: точная
сторона читает 0, а род берётся из raw (`zero_sign`, `Ieee.zero_kind`) и
объявляется отдельно — один бит, живущий только в представлении.

NOT_MEASURED (запись 05.09) — РАЗОБРАН 17.09 на два РАЗНЫХ утверждения; строка
остаётся записью о том, что здесь было объявлено, а каждое закрыто своим родом
доказательства.

1. Режимы округления кроме ближайшего-чётного — ЗАКРЫТО точной дорогой. Оговорка
   «Python использует только его» верна о ПЛОЩАДКЕ и неверна о предмете:
   направление округления есть правило выбора одного из двух соседей на линейке,
   и оно считается ЦЕЛЫМИ, без float и без сетки. Все пять направлений §4.3
   стоят в `round_to` (`MODES`), соседи — в `neighbours`; переполнение разведено
   по §7.4 (ближайшие дают ±∞, направленное — ∞ лишь в свою сторону, против неё
   наибольший конечный; замер 17.09 `[вычислено здесь]`: первый ход отдавал ±∞
   всем пятерым, и это было неверно для трёх из пяти). Род нуля сохраняется
   знаком величины. Две независимые дороги приёмки: КРЕМНИЙ — образ Python для
   + × / на 2000 парах обязан совпасть с целым счётом при RNE (совпал 2000 из
   2000 по каждой операции); LIBM — `math.nextafter` обязан делать пару
   `neighbours` смежной, а вилка — зажимать величину (2000 из 2000). Ничья
   разводит RNE и RNA: на 4000 ровных серединах они разошлись ровно 2000 раз.

2. Операции сверх + − × / и целой степени — ЗАКРЫТО точной дорогой для той, что
   стандарт ТРЕБУЕТ, и названным отказом с ценой для тех, что он лишь советует.
   Это два разных утверждения, и смешивать их нельзя:

   * √ — пятая обязательная операция §5.4.1, и она считается точно: разряд по
     разряду подкоренного, засечка через `isqrt`, выбор соседа сравнением 4A с
     (2k+1)²·B — целыми, без `math.sqrt` (`sqrt_image`). Замер: 3000 значений,
     целая дорога совпала с libm 3000 из 3000; направленные режимы зажимают
     истину 500 из 500; √(−0) = −0, √(отрицательного) — nan по §7.2. Точная
     сторона разведена честно: √v рационален ⟺ A = k²·B, и тогда потеря ровно 0;
     иначе величина СУЩЕСТВУЕТ в ℝ, но конечной ℚ-записи у неё нет, и потеря
     названа ВИЛКОЙ, а не выдана числом, которого не существует.

   * exp, log, sin, cos, tan, atan, нецелая степень (`RECOMMENDED`) — операции
     §9.2, РЕКОМЕНДУЕМЫЕ: верного округления стандарт для них не требует вовсе.
     ОТКАЗ ДОМА, названный целиком. ПРЕДМЕТ: точная сторона этих операций как
     ЧИСЛО. ПРИЧИНА: величина трансцендентна (Линдеман 1882 — exp при
     рациональном ненулевом доводе), и конечной ℚ-записи у неё нет; это НЕ «не
     умею посчитать» и НЕ «величины не существует». АДРЕС: `image_from_enclosure`
     здесь, `phicore/enclosure.py` (exp/ln/sin/cos — целые вилки). ЧЕМ МЕРИЛОСЬ
     БЫ и цена: образ решается ТОГДА И ТОЛЬКО ТОГДА, когда вилка не накрывает
     границу округления, и цена платится разрядами — замер 17.09
     `[вычислено здесь]`: exp на 41 доводе решён 29 при 24 членах и все 41 при
     32; ln — 22 при 40 членах, 36 при 64, все 41 при 96. Не решённых на верхнем
     разряде — 0 из 82: отказ ВРЕМЕНЕН и оплачиваем, и потому назван ценой, а не
     бессилием. Где вилка решила, образ сверен с libm на одном и том же точном
     доводе — 180 из 180.

   Попутный замер, важный для чтения чисел: сверять exp(k/7) с `math.exp` от
   ПОДНЯТОГО довода нельзя — это две разные величины (разошлись 15 из 41). Довод
   обеим дорогам обязан даваться один и тот же, точно.
"""
from __future__ import annotations

import math

from .enclosure import sqrt_enclosure
from .errors import Unresolved
from .hyper import Hyper
from .rational import Q
from .wheels import top_wheel2
from .zero_state import Zero

_NUM = (int, Q)


# ---------------------------------------------------------------------------
# точная сторона одного float
# ---------------------------------------------------------------------------

def exact(x: float):
    """Точная сторона float: Q-константа Hyper для конечного (двоичная дробь
    битово точно), ±ω для ±inf, состояние NOT_COMPUTED для nan.

    Одноимённые двери соседей — ДРУГИЕ понятия: `window.Window.exact(x)` —
    показание одометра окна целиком (образ плюс поглощённое, для элемента
    контейнера), `affine.Affine.exact(g, units)` — конструктор формы без
    символов шума из узла сетки. Здесь вход — один float, и только он.
    """
    if isinstance(x, bool) or not isinstance(x, (float, int)):
        raise TypeError("exact — только float (или int)")
    x = float(x)
    if math.isnan(x):
        return Zero.not_computed()
    if math.isinf(x):
        return Hyper.omega(1 if x > 0 else -1)
    n, d = x.as_integer_ratio()
    return Hyper.const(Q(n, d))


def dyadic(x: float) -> Q:
    """Конечный float как точная двоичная дробь; inf/nan — TypeError."""
    if isinstance(x, bool) or not isinstance(x, float) or not math.isfinite(x):
        raise TypeError("dyadic — только конечный float")
    n, d = x.as_integer_ratio()
    return Q(n, d)


def is_dyadic(q) -> bool:
    """Представима ли величина float-дробью точно — по закону РАНГА колеса.

    «Знаменатель — степень двойки» есть частный случай закона
    `W_b(x) < ∞ ⟺ x ∈ ℤ[1/b]` при b = 2 (`phicore/wheels.py`): конечная запись
    колёсами сетки (2, 1) и есть точная float-дробь. Дверь холодная — замер
    22.09 `[вычислено здесь]`: 2 вызова за набор `test_phicore_ieee`, — потому
    здесь стоит закон, а не битовый приём.
    """
    from .wheels import INFINITE, wheel_rank          # круг загрузки: ленивый
    return wheel_rank(Q(q), 2) is not INFINITE


#: наименьшее нормальное число double; ниже него — субнормальный край
NORMAL_MIN = 2.0 ** -1022


def zero_sign(x: float) -> int | None:
    """Род нуля: +1 у `+0.0`, −1 у `−0.0`, None у всего остального.

    `[вычислено здесь]` Для величины рода нет: `−0.0 == 0.0` истинно, и точная
    сторона читает оба как точный 0. Род живёт только в представлении — один
    бит, который IEEE хранит и пускает в дело при делении (`1/(−0) = −inf`).
    Тот же случай, что и выброшенный знак в теореме Пифагора: модуль рода не
    несёт, а деление и обратная величина — несут.
    """
    if isinstance(x, bool) or not isinstance(x, (float, int)):
        raise TypeError("zero_sign — только float (или int)")
    x = float(x)
    if x != 0 or math.isnan(x):
        return None
    return int(math.copysign(1.0, x))


def subnormal(x: float) -> bool:
    """Лежит ли конечное число в субнормальном крае: 0 < |x| < 2^−1022.

    `[вычислено здесь]` Замер 11.09 на 3000 парах субнормальных: сложение
    остаётся точным 3000 из 3000 (сумма не покидает сетку 2^−1074), умножение
    точно 0 из 3000 (произведение уходит под засечку целиком). Потому потеря
    края — не «мелочь округления», а полное исчезновение произведения.
    """
    if isinstance(x, bool) or not isinstance(x, (float, int)):
        raise TypeError("subnormal — только float (или int)")
    x = float(x)
    return math.isfinite(x) and x != 0 and abs(x) < NORMAL_MIN


def loss(exact_value, raw: float):
    """Потеря IEEE: точная величина минус то, что дал float.

    Конечный raw и конечная точная сторона — Q (0 означает «без потерь»);
    raw = ±inf при конечной точной — «переполнение: точно <величина>»;
    точная сторона бесконечна (несёт ω) и raw = ±inf — потеря есть
    конечная часть, которую IEEE отбросил (Q); raw = nan — «nan: точно
    <запись>»; точная сторона неизвестна — «не сведено»."""
    if exact_value is None:
        return "не сведено"
    if isinstance(exact_value, Zero):
        return f"состояние: {exact_value.record()}"
    if not isinstance(exact_value, Hyper):
        return "не сведено"
    if math.isnan(raw):
        return f"nan: точно {exact_value.describe()}"
    order = exact_value.order
    finite_part = exact_value.terms.get(0, Q(0))
    if math.isinf(raw):
        if order is not None and order > 0:
            # IEEE удержал масштаб (inf при ω), отбросил всё ниже старшего члена
            dropped = Hyper({k: c for k, c in exact_value.terms.items() if k < order})
            return dropped.terms.get(0, Q(0)) if dropped.is_monomial or not dropped.terms else f"отброшено: {dropped.describe()}"
        return f"переполнение: точно {exact_value.describe()}"
    if order is not None and order != 0:
        # точная сторона несёт ω или ε, а IEEE дал конечное: масштаб потерян
        return f"масштаб потерян: точно {exact_value.describe()}, IEEE {raw!r}"
    return finite_part - dyadic(raw)


# ---------------------------------------------------------------------------
# линейка IEEE как целое (письмо автора 05.09 «IEEE error canceller»)
# ---------------------------------------------------------------------------

#: наименьшая засечка линейки double — наименьший субнормальный 2^−1074
TICK = Q(1, 2 ** 1074)
#: наибольший конечный double = 2^1024 − 2^971; в засечках — целое K
RULER_K = 2 ** 2098 - 2 ** 2045
#: засечек от −max до +max включительно (0 — одна засечка): целое
RULER_TICKS = 2 * RULER_K + 1


def tick_index(x: float) -> int:
    """Номер засечки конечного double на целочисленной линейке: x = k·2^−1074,
    k ∈ ℤ, |k| ≤ RULER_K. Всякий конечный double — целое число засечек."""
    k = dyadic(x) / TICK
    if k.denominator != 1:
        raise ArithmeticError("double вне целочисленной линейки — невозможно по построению IEEE 754")
    return int(k.numerator)


def finest_unit(q) -> int:
    """Показатель e, при котором q = j·2^e с целым j (знаменатель q — степень
    двойки); для q вне двоичных дробей — ValueError."""
    q = Q(q)
    d = q.denominator
    if d & (d - 1):
        raise ValueError("не двоичная дробь: единицы засечки нет")
    return -(d.bit_length() - 1)


def representable(q) -> bool:
    """Лежит ли двоичная дробь q на линейке double точно (после округления в
    double значение совпадает с q)."""
    q = Q(q)
    try:
        f = q.numerator / q.denominator
    except OverflowError:
        return False
    return math.isfinite(f) and dyadic(f) == q


def remainder(a: float, b: float, op: str) -> dict:
    """Канцеллер: a ∘ b = fl + r точно. fl — образ IEEE; r — остаток целыми
    (для «/»: r = a − fl·b, так что a = fl·b + r); «на линейке» — лежит ли r
    на линейке double (тогда засечка, потерянная с одной стороны, приобретена
    с другой — на той же линейке); иначе r — целое число засечек более мелкой
    линейки с единицей 2^unit (конец линейки: субнормальный край или INF)."""
    if op not in "+-*/":
        raise ValueError("op ∈ + − × /")
    if not (math.isfinite(a) and math.isfinite(b)):
        raise ValueError("канцеллер — для конечных операндов; inf/nan — указатели и состояния (exact)")
    da, db = dyadic(a), dyadic(b)
    if op == "+":
        fl, ex = a + b, da + db
    elif op == "-":
        fl, ex = a - b, da - db
    elif op == "*":
        fl, ex = a * b, da * db
    else:
        if b == 0:
            raise ZeroDivisionError("a/0: остаток не определён (x/0 — не число)")
        fl, ex = a / b, None
    if math.isinf(fl):
        # конец линейки INF: образ — указатель, точная величина конечна (ex) — вся она «потеряна» образом;
        # для «/» точное частное — не двоичная дробь вообще, остаток к inf не определён
        return {"fl": fl, "r": ex, "на линейке": False, "unit": None, "конец": "INF"}
    if op == "/":
        r = da - dyadic(fl) * db
    else:
        r = ex - dyadic(fl)
    on = representable(r)
    return {"fl": fl, "r": r, "на линейке": on, "unit": finest_unit(r) if r != 0 else 0,
            "конец": None if on else ("субнормальный край" if r != 0 and abs(r) < TICK * 2 ** 53 else "за линейкой")}

# ---------------------------------------------------------------------------
# направления округления IEEE 754 §4.3: устройство, которым ЗАКРЫТО утверждение 1
# пометки NOT_MEASURED (шапка модуля, разбор 17.09) — счёт целыми, ни одного float
# ---------------------------------------------------------------------------

#: Пять направлений округления стандарта. Python считает ТОЛЬКО первым, но это
#: свойство площадки, а не предмета: направление есть правило выбора одного из
#: двух соседей на линейке, и оно считается целыми — без float и без сетки.
RNE = "ближайший-чётный"        # roundTiesToEven — умолчание IEEE 754 и всё, что даёт Python
RNA = "ближайший-от-нуля"       # roundTiesToAway
RTZ = "к нулю"                  # roundTowardZero
RTP = "к +∞"                    # roundTowardPositive
RTN = "к −∞"                    # roundTowardNegative
MODES = (RNE, RNA, RTZ, RTP, RTN)

#: наибольший конечный double в засечках старшего разряда: (2^53 − 1)·2^971
_MAX_K, _MAX_U = 2 ** 53 - 1, 971


def _scale(k: int, u: int) -> float:
    """k·2^u как float ТОЧНО: |k| ≤ 2^53 представимо, ldexp — сдвиг порядка.

    `[вычислено здесь]` Это построение, а не арифметика: во float здесь ничего
    не складывается и не умножается (правило дома — float-арифметики в phicore
    нет). Точность построения не обещается словом, а проверяется приёмкой
    обратным чтением `dyadic()` в ℚ.
    """
    if abs(k) > 2 ** 53:
        raise ValueError("засечка шире 53 бит: построение перестало быть точным")
    return math.ldexp(float(k), u)


def _binade(n: int, d: int) -> int:
    """e с 2^e ≤ n/d < 2^(e+1) для n, d > 0 — целыми, по длинам записи.

    ЗАКОН здесь — старший номер колеса записи в сетке (2, 1)
    (`wheels.Record.of(x, 2, width=1)`): бинада есть позиция старшего колеса,
    и ничего сверх этого. Длины двоичной записи — БЫСТРАЯ ТРОПА к тому же
    числу, а не отдельная истина.

    Тропа оставлена замером `[вычислено здесь]` 22.09: за набор
    `test_phicore_ieee` эта дверь зовётся 32 580 раз, а колесо на ней дороже в
    41 раз; равносильность проверяется переборно приёмкой
    (`proto/test_phicore_ieee.py`, раздел «бинада против колеса», с подлогом).
    Если тропа и закон разойдутся хоть на одной паре — краснеет приёмка, а не
    молчит библиотека.

    Замер 22.09 `[вычислено здесь]`: та же тропа была переписана руками и в
    `window._floor_log2` — две копии одного закона, ничем не сверяемые между
    собой. Копия убрана: тропа живёт ОДНА в модуле закона (`wheels.top_wheel2`),
    здесь — только вызов.
    """
    return top_wheel2(n, d)


def neighbours(q) -> tuple:
    """Два соседних double вокруг точной величины q: (lo, hi, точно ли).

    `[вычислено здесь]` Целыми: разряд по длинам записи, засечка k = ⌊|q|/2^u⌋
    делением целых. lo — наибольший double ≤ q, hi — наименьший double ≥ q; при
    точном q оба равны ему. Ни одного float в счёте: float появляется только
    построением `_scale` и обязан читаться обратно той же дробью.
    """
    q = Q(q)
    n, d = int(q.numerator), int(q.denominator)
    if n == 0:
        return 0.0, 0.0, True
    sign, n = (-1, -n) if n < 0 else (1, n)
    e = _binade(n, d)
    u = max(e - 52, -1074)
    if u <= 0:
        num, den = n << (-u), d
    else:
        num, den = n, d << u
    k, rem = divmod(num, den)
    if rem == 0:
        v = _over(sign, k, u, RTZ if abs(k) else RNE)
        return v, v, True
    lo_k, hi_k = (k, k + 1) if sign > 0 else (k + 1, k)
    # нижний сосед — как округление к −∞, верхний — как к +∞ (на концах линейки
    # это и есть разница между наибольшим конечным и ±∞)
    return _over(sign, lo_k, u, RTN), _over(sign, hi_k, u, RTP), False


def _over(sign: int, k: int, u: int, mode: str = RNE):
    """Засечка k·2^u со знаком — либо float, либо переполнение по §7.4.

    Переполнение стандарт разводит ПО НАПРАВЛЕНИЮ, а не одинаково: ближайшие
    дают ±∞; направленное даёт ∞ только в свою сторону, а против неё —
    наибольший конечный. Замер 17.09 `[вычислено здесь]`: первый ход отдавал
    ±∞ всем пятерым, и это было неверно для трёх из пяти.
    """
    # сравнение целыми: k·2^u против (2^53 − 1)·2^971 — без float и без логарифма
    if (k << max(0, u - _MAX_U)) <= (_MAX_K << max(0, _MAX_U - u)):
        return _scale(sign * k, u)
    big = _scale(sign * _MAX_K, _MAX_U)                  # наибольший конечный того же знака
    inf = math.inf if sign > 0 else -math.inf
    if mode in (RNE, RNA):
        return inf
    if mode == RTZ:
        return big
    if mode == RTP:
        return inf if sign > 0 else big
    return inf if sign < 0 else big                      # RTN


def round_to(q, mode: str = RNE) -> float:
    """Образ точной величины q в double по направлению `mode` — счёт целыми.

    `[вычислено здесь]` Здесь пометка NOT_MEASURED ЗАКРЫТА 17.09, утверждение 1: «режимы
    округления кроме ближайшего-чётного» не были ни невычислимы, ни
    несуществующи — их не считали. Направление есть правило выбора одного из
    двух соседей (`neighbours`), и оно решается сравнением ЦЕЛЫХ:

      * к нулю / к +∞ / к −∞ — сосед назван направлением;
      * ближайший — сравнением 2·остаток с засечкой; при ровной середине
        RNE берёт чётную засечку, RNA — дальнюю от нуля.

    Переполнение по §7.4 стандарта: ближайшие дают ±∞, направленные — ±∞ лишь
    в свою сторону, а против своего направления возвращают наибольший конечный.
    Род нуля сохраняется знаком величины (−0 при отрицательном q, идущем в 0).
    """
    if mode not in MODES:
        raise ValueError(f"направление округления — одно из {MODES}")
    q = Q(q)
    n, d = int(q.numerator), int(q.denominator)
    if n == 0:
        return 0.0
    sign, n = (-1, -n) if n < 0 else (1, n)
    e = _binade(n, d)
    u = max(e - 52, -1074)
    if u <= 0:
        num, den = n << (-u), d
    else:
        num, den = n, d << u
    k, rem = divmod(num, den)
    if rem == 0:
        return _over(sign, k, u, mode)
    if mode == RTZ:
        near = k
    elif mode == RTP:
        near = k if sign < 0 else k + 1
    elif mode == RTN:
        near = k + 1 if sign < 0 else k
    else:
        twice = 2 * rem
        if twice > den:
            near = k + 1
        elif twice < den:
            near = k
        else:                                    # ровная середина — правило ничьей
            near = k + 1 if (mode == RNA or k % 2) else k
    if near == 0:
        return -0.0 if sign < 0 else 0.0         # род нуля живёт в представлении
    return _over(sign, near, u, mode)


def rounding_table(q) -> dict:
    """Все пять образов одной величины разом — для сверки и для записи."""
    return {m: round_to(q, m) for m in MODES}


# ---------------------------------------------------------------------------
# операции сверх + − × / и целой степени: здесь ЗАКРЫТО утверждение 2 пометки
# NOT_MEASURED (шапка модуля, разбор 17.09) — обязательная √ точной дорогой,
# рекомендуемые §9.2 — названным отказом с оплачиваемой ценой
# ---------------------------------------------------------------------------

def sqrt_image(x, mode: str = RNE) -> dict:
    """√x с ВЕРНЫМ округлением — пятая обязательная операция IEEE 754 §5.4.1.

    `[вычислено здесь]` Считается целыми, без float и без `math.sqrt`: разряд
    итога по разряду подкоренного, засечка k = ⌊√v·2^(−u)⌋ через `isqrt`
    (k = isqrt(A·B)//B при v = A/B), выбор соседа — сравнением 4A с (2k+1)²·B,
    то есть сравнением ЦЕЛЫХ, а не сравнением приближений.

    Точная сторона разведена честно. √v рационален ТОГДА И ТОЛЬКО ТОГДА, когда
    A = k²·B; тогда потеря ровно 0 и записана дробью. Иначе величина
    СУЩЕСТВУЕТ в ℝ, но конечной ℚ-записи у неё нет — и это третье утверждение,
    не «не умею» и не «величины нет»: потеря тогда названа ВИЛКОЙ на объявленной
    сетке (`sqrt_enclosure`), а не выдана числом, которого не существует.
    """
    if mode not in MODES:
        raise ValueError(f"направление округления — одно из {MODES}")
    if isinstance(x, Ieee):
        v, raw_in = x.exact, x.raw
    else:
        v, raw_in = exact(x), float(x)
    if isinstance(v, Zero):
        return {"образ": math.nan, "точно": v, "потеря": "nan: путь потерян", "вилка": None}
    if not isinstance(v, Hyper) or v.order not in (0, None):
        return {"образ": math.nan, "точно": v, "потеря": "масштаб вне линейки: √ω не число IEEE", "вилка": None}
    q = v.terms.get(0, Q(0))
    if q < 0:
        return {"образ": math.nan, "точно": Zero.not_computed(),
                "потеря": "√(отрицательного) — не вещественно: nan по §7.2", "вилка": None}
    if q == 0:
        return {"образ": math.copysign(0.0, raw_in), "точно": Hyper.const(Q(0)),
                "потеря": Q(0), "вилка": None}
    n, d = int(q.numerator), int(q.denominator)
    e_r = _binade(n, d) // 2
    u = max(e_r - 52, -1074)
    A, B = n << max(0, -2 * u), d << max(0, 2 * u)
    k = math.isqrt(A * B) // B
    точно = (A == k * k * B)
    if точно:
        near = k
    else:
        cmp = 4 * A - (2 * k + 1) ** 2 * B          # v против середины (k + ½)·2^u — целыми
        if mode == RTZ or mode == RTN:
            near = k
        elif mode == RTP:
            near = k + 1
        elif cmp > 0:
            near = k + 1
        elif cmp < 0:
            near = k
        else:                                        # ровная середина: у √ невозможна, но правило названо
            near = k + 1 if (mode == RNA or k % 2) else k
    образ = _over(1, near, u, mode)
    if точно:
        val = Q(near) * Q(2) ** u if u >= 0 else Q(near, 2 ** -u)
        return {"образ": образ, "точно": Hyper.const(val),
                "потеря": (val - dyadic(образ)) if math.isfinite(образ) else "переполнение",
                "вилка": None}
    iv = sqrt_enclosure(q, bits=80)
    return {"образ": образ, "точно": None,
            "потеря": f"иррационален: точной ℚ-записи нет; вилка [{iv.lo}; {iv.hi}] шириной ≤ 2^−79",
            "вилка": iv}


#: Операции §9.2 стандарта — РЕКОМЕНДУЕМЫЕ, а не обязательные: верного
#: округления стандарт для них не требует вовсе, и потому «образ» у них не один.
RECOMMENDED = ("exp", "log", "sin", "cos", "tan", "atan", "pow(нецелая степень)")

_PRICE = ("величина существует в ℝ, но конечной ℚ-записи у неё нет (трансцендентность: "
          "Линдеман 1882 для exp при рациональном ненулевом доводе), потому точная "
          "сторона здесь не число, а ВИЛКА. Образ решается тогда и только тогда, когда "
          "вилка не накрывает границу округления; цена — разряды сетки вилки")


def image_from_enclosure(iv, mode: str = RNE) -> float:
    """Образ вещественной величины, данной ВИЛКОЙ [lo, hi], — либо отказ с ценой.

    `[вычислено здесь]` Это и есть честная граница утверждения 2 пометки. Если
    оба конца вилки округляются в один и тот же double, образ ОПРЕДЕЛЁН и верен
    — какова бы ни была величина внутри. Если концы расходятся, вилка не решает,
    и отказ называет цену в разрядах, а не говорит «не умею»: ПРЕДМЕТ —
    верно округлённый образ; ПРИЧИНА — вилка шире границы округления (дилемма
    составителя таблиц), не отсутствие дороги; АДРЕС — здесь и
    `phicore/enclosure.py`; ЧЕМ МЕРИЛОСЬ БЫ — уже́нием вилки: каждый разряд
    сетки делит ширину надвое, и достаточная ширина известна замером.
    """
    if mode not in MODES:
        raise ValueError(f"направление округления — одно из {MODES}")
    lo, hi = round_to(iv.lo, mode), round_to(iv.hi, mode)
    if lo == hi and not (math.isnan(lo) or math.isnan(hi)):
        return lo
    raise Unresolved("округление", f"вилка [{iv.lo}; {iv.hi}]", f"направление {mode}",
                     hint=f"вилка не решает образ: концы дают {lo!r} и {hi!r}. " + _PRICE)


# ---------------------------------------------------------------------------
# две дороги к одной потере: безошибочные преобразования
# ---------------------------------------------------------------------------

def two_sum(a: float, b: float) -> tuple[float, float]:
    """Knuth 1969: s = fl(a + b), e — точная ошибка, сама float; a + b = s + e."""
    s = a + b
    bb = s - a
    e = (a - (s - bb)) + (b - bb)
    return s, e


def two_prod(a: float, b: float) -> tuple[float, float]:
    """Dekker 1971 (разбиение Veltkamp): p = fl(a·b), e — точная ошибка; a·b = p + e."""
    p = a * b
    split = 134217729.0  # 2^27 + 1
    def _split(x):
        c = split * x
        hi = c - (c - x)
        return hi, x - hi
    ah, al = _split(a)
    bh, bl = _split(b)
    e = ((ah * bh - p) + ah * bl + al * bh) + al * bl
    return p, e


# ---------------------------------------------------------------------------
# число класса 23 с точной стороной
# ---------------------------------------------------------------------------

def _lift_exact(x):
    """Точная сторона операнда: Ieee → exact; float → exact(float); int/Q →
    Hyper-константа; Hyper/Zero — как есть; иное — None (не сведено)."""
    if isinstance(x, Ieee):
        return x.exact
    if isinstance(x, bool):
        return None
    if isinstance(x, float):
        return exact(x)
    if isinstance(x, _NUM):
        return Hyper.const(Q(x))
    if isinstance(x, (Hyper, Zero)):
        return x
    return None


def _raw_of(x) -> float:
    if isinstance(x, Ieee):
        return x.raw
    if isinstance(x, float):
        return x
    if isinstance(x, _NUM):
        return int(Q(x).numerator) / int(Q(x).denominator) if Q(x).denominator != 1 else float(int(Q(x).numerator))
    if isinstance(x, Hyper):
        return x.approx()
    if isinstance(x, Zero):
        return 0.0
    raise TypeError(f"IEEE-образа нет: {x!r}")


def _fdiv(a: float, b: float) -> float:
    """Деление как в IEEE 754: x/0 = ±inf (знак — произведение знаков), 0/0 и nan — nan."""
    if b == 0:
        if a == 0 or math.isnan(a):
            return math.nan
        return math.copysign(math.inf, a) * math.copysign(1.0, b)
    return a / b


def _exact_op(a, b, op: str):
    """Точная сторона операции; состояния решают сами (NOT_COMPUTED — отказ →
    остаётся состоянием), деление на многочлен — не сведено (None)."""
    if a is None or b is None:
        return None
    try:
        if op == "+":
            return a + b
        if op == "-":
            return a - b
        if op == "*":
            return a * b
        if op == "/":
            if isinstance(b, Zero):
                return Zero.not_computed()
            if isinstance(b, Hyper) and not b.terms:
                return Zero.not_computed()       # x/0 — не число
            return a / b
    except Exception:
        return Zero.not_computed()
    return None


class Ieee:
    """Класс 23 с точной стороной: raw — как считает IEEE, exact — точно,
    loss — что потеряно (число либо имя)."""

    __slots__ = ("raw", "exact", "loss", "path")

    def __init__(self, raw: float, exact_value=None, loss_value=None, path: str = ""):
        if isinstance(raw, bool) or not isinstance(raw, float):
            raise TypeError("raw — только float")
        self.raw = raw
        self.exact = exact(raw) if exact_value is None else exact_value
        self.loss = loss(self.exact, raw) if loss_value is None else loss_value
        self.path = path

    @classmethod
    def of(cls, x) -> "Ieee":
        """float → класс 23 с точной стороной без потерь; int/Q — образ IEEE и
        точная величина (потеря — округление ввода, если дробь не двоичная)."""
        if isinstance(x, Ieee):
            return x
        if isinstance(x, bool):
            raise TypeError("bool — не число")
        if isinstance(x, float):
            return cls(x)
        if isinstance(x, _NUM):
            q = Q(x)
            raw = int(q.numerator) / int(q.denominator)
            return cls(raw, Hyper.const(q), loss(Hyper.const(q), raw), "ввод")
        raise TypeError(f"Ieee.of: {x!r}")

    def _binary(self, other, op: str, swap: bool = False):
        o = other if isinstance(other, Ieee) else None
        if o is None:
            if isinstance(other, bool) or not isinstance(other, (float, Hyper, Zero, *_NUM)):
                return NotImplemented
            o = Ieee(_raw_of(other), _lift_exact(other), Q(0) if not isinstance(other, float) else Q(0), "операнд")
        a, b = (o, self) if swap else (self, o)
        ra, rb = a.raw, b.raw
        try:
            raw = {"+": ra + rb, "-": ra - rb, "*": ra * rb, "/": _fdiv(ra, rb)}[op]
        except OverflowError:
            raw = math.inf
        ex = _exact_op(a.exact, b.exact, op)
        return Ieee(raw, ex, loss(ex, raw), f"({a.path or 'x'} {op} {b.path or 'y'})")

    def __add__(self, other):
        return self._binary(other, "+")

    def __radd__(self, other):
        return self._binary(other, "+", swap=True)

    def __sub__(self, other):
        return self._binary(other, "-")

    def __rsub__(self, other):
        return self._binary(other, "-", swap=True)

    def __mul__(self, other):
        return self._binary(other, "*")

    def __rmul__(self, other):
        return self._binary(other, "*", swap=True)

    def __truediv__(self, other):
        return self._binary(other, "/")

    def __rtruediv__(self, other):
        return self._binary(other, "/", swap=True)

    def __neg__(self):
        ex = None if self.exact is None else (-self.exact if isinstance(self.exact, Hyper) else self.exact)
        return Ieee(-self.raw, ex, loss(ex, -self.raw), f"−({self.path or 'x'})")

    def __pow__(self, n: int):
        if isinstance(n, bool) or not isinstance(n, int):
            return NotImplemented
        try:
            raw = self.raw ** n
        except OverflowError:
            raw = math.inf
        except ZeroDivisionError:
            raw = math.inf
        ex = None
        if isinstance(self.exact, Hyper):
            try:
                ex = self.exact ** n
            except Exception:
                ex = None
        elif isinstance(self.exact, Zero):
            ex = self.exact
        return Ieee(float(raw), ex, loss(ex, float(raw)), f"({self.path or 'x'})^{n}")

    # --- запись --------------------------------------------------------------

    @property
    def kind(self) -> str:
        return "nan" if math.isnan(self.raw) else ("inf" if math.isinf(self.raw) else "finite")

    def zero_kind(self) -> str | None:
        """Род нуля значения: «+0», «−0» либо None. Точная сторона его не несёт."""
        s = zero_sign(self.raw)
        return None if s is None else ("+0" if s > 0 else "−0")

    def subnormal(self) -> bool:
        """Лежит ли значение в субнормальном крае (край назван, а не проглочен)."""
        return subnormal(self.raw)

    def self_check(self, deep: bool = False) -> list:
        """Самопроверка величины класса 23 независимой дорогой (очередь п. 18).

        `[вычислено здесь]` Проверяется не тем счётом, каким строилось:
        потеря пересчитывается из пары (точная сторона, raw) заново и
        сверяется с записанной; у конечного значения обязано выполняться
        тождество линейки `dyadic(raw) + потеря = точная величина`; род
        значения (`finite/inf/nan`), род нуля и субнормальный край обязаны
        отвечать самому raw. При `deep` — двойной проход: значение лежит на
        целой засечке линейки (`tick_index(raw)·2^−1074 = dyadic(raw)`), оно
        представимо, и пересборка `Ieee.of(raw)` даёт ту же величину без
        потери.

        Граница: путь (`path`) здесь не переигрывается — запись о дороге, а
        не величина; сверка пути есть отдельная работа `remainder`.
        """
        out = []
        if isinstance(self.raw, bool) or not isinstance(self.raw, float):
            out.append(f"raw не float: {type(self.raw).__name__}")
            return out
        again = loss(self.exact, self.raw)
        if type(again) is not type(self.loss) or again != self.loss:
            out.append(f"потеря пересчётом {again!r} не сходится с записанной {self.loss!r}")
        want = "nan" if math.isnan(self.raw) else ("inf" if math.isinf(self.raw) else "finite")
        if self.kind != want:
            out.append(f"род значения {self.kind} не отвечает raw={self.raw!r}")
        if (self.zero_kind() is None) != (zero_sign(self.raw) is None):
            out.append(f"род нуля {self.zero_kind()} не отвечает raw={self.raw!r}")
        if self.subnormal() != subnormal(self.raw):
            out.append("субнормальный край назван не по raw")
        if math.isfinite(self.raw) and isinstance(self.exact, Hyper) and self.exact.order in (0, None) \
                and isinstance(self.loss, Q):
            finite_part = self.exact.terms.get(0, Q(0))
            if dyadic(self.raw) + self.loss != finite_part:
                out.append("тождество линейки нарушено: dyadic(raw) + потеря ≠ точная величина")
        if not deep:
            return out
        if math.isfinite(self.raw):
            q = dyadic(self.raw)
            if not representable(q):
                out.append(f"величина raw не представима на линейке: {q}")
            if tick_index(self.raw) * TICK != q:
                out.append("номер засечки не сходится с двоичной дробью raw")
            back = Ieee.of(self.raw)
            if back.raw != self.raw or back.loss != Q(0):
                out.append(f"пересборка Ieee.of(raw) дала иное: {back.record()}")
        return out

    def exact_record(self) -> str:
        if self.exact is None:
            return "не сведено"
        return self.exact.describe() if isinstance(self.exact, Hyper) else self.exact.record()

    def record(self) -> str:
        l = self.loss if isinstance(self.loss, str) else f"{self.loss}"
        return f"Φ(C=23 IEEE raw={self.raw!r}, точно={self.exact_record()}, потеря={l})"

    __repr__ = record

    def __eq__(self, other):
        if not isinstance(other, Ieee):
            return NotImplemented
        return self.exact == other.exact and (self.raw == other.raw or (math.isnan(self.raw) and math.isnan(other.raw)))

    def __hash__(self):
        return hash(("ieee", self.raw if not math.isnan(self.raw) else "nan"))
