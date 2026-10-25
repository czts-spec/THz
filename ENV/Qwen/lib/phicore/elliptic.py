"""Линейка класса 15 ELLIPTIC: полные интегралы K(m) и E(m), m = k² ∈ ℚ.

Величина хранится структурно: единственный слот — параметр m (квадрат
модуля), канонически m ∈ (0, 1). Края области приземляются фабриками:
K(0) = E(0) = π/2 — моном класса 10 (двойственность, уже объявленная в
pi_value); E(1) = 1 — вырождение в класс 1; K(1) — расходимость, честный
отказ. Вне [0, 1] сведение (мнимый и обратный модуль) не реализовано —
отказ, не выдумка.

Ось двойственности — дополнительный модуль m′ = 1 − m (complementary):
инволюция m′′ = m держится тождеством 1 − (1 − m) = m, самодуальная
точка m = 1/2. Связь сторон — соотношение Лежандра
E(m)·K(m′) + E(m′)·K(m) − K(m)·K(m′) = π/2 [внешний источник: Лежандр,
Абрамовиц—Стиган 17.3.13]; направленного правила-приземления оно здесь
не даёт (произведения линейкой не несутся) и потому несётся
сертификатом сетки (legendre_enclosure).

Десятичная сторона — только сертификатом сетки (enclosure) по явному
запросу: AGM в интервальной арифметике,

    K(m) = π / (2·AGM(1, √(1−m))),
    E(m) = K(m)·(1 − Σ_{n≥0} 2^{n−1}·c_n²),  c_0² = m,
    a_{n+1} = (a_n+b_n)/2,  b_{n+1} = √(a_n·b_n),  c_{n+1} = (a_n−b_n)/2

[внешний источник: метод AGM, Абрамовиц—Стиган 17.6]. Корни берутся
sqrt_enclosure от концов отрезка (монотонность √), хвост суммы
мажорируется геометрической прогрессией из тождества
c_{n+1} = c_n²/(4·a_{n+1}) — теоремы анализа, применённые в точной
арифметике Q; ни одного float.

Двойственность с классом 7 (05.09, очередь п. 41 «K·E»): K(m) = (π/2)·S_K,
E(m) = (π/2)·S_E, где S_K = ₂F₁(½, ½; 1; m) и S_E = ₂F₁(−½, ½; 1; m) —
гипергеометрические ряды с рациональным отношением членов
(2n+1)²/(4(n+1)²)·m и (4n²−1)/(4(n+1)²)·m [внешний источник: Абрамовиц—
Стиган 17.3.9–17.3.12], то есть D-конечные последовательности порядка 1.
Метод series() отдаёт S как ряд класса 7 с сертификатом от самой линейки
(AGM, делённый на π/2 — работает при всяком m ∈ (0, 1), тогда как
геометрическая мажоранта ряда держит лишь m < 4/9). Отсюда K·E, K·K′,
K ± E, c·K, K/c, Kⁿ, −K, K·ряд/интеграл/произведение — структурная запись
π-моном × D-конечный ряд (произведение и сумма рядов — замыканием
dfinite.py, рекуррента точна).

Единица и ноль множителем записи НЕ строят: класс ответа зависит здесь от
ЗНАЧЕНИЯ довода, а не только от его рода (замер 12.09, гл. 07 §19 п. 9).
c = 1 отдаёт САМ НОСИТЕЛЬ — `K(1/3)·1`, `1·K(1/3)`, `K(1/3)/1` и
`K(1/3)**1` дают тот же `KValue` (и тот же объект), то же у `EValue`;
c = 0 — точный ноль класса 1 (`K(1/3)·0` → `Q(0)`); деление на 0 — отказ;
всякое иное c даёт класс 19 (`·3` → π-моном q=3/2, `/3` → q=1/6,
`·(−1)` и `/(−1)` → q=−1/2, `**2` → q=1/4, e=2). Тождество на единице
есть решение устройства, а не сокращение счёта: слота масштаба у линейки
нет, и записи «1·K» она не несёт — потому единственный способ остаться в
классе 15 есть вернуть носитель. У соседей по главе класс от значения
довода НЕ зависит, потому что масштаб у них лежит в слоте: у класса 14 —
в `mult` (`Γ(1/3)·1` и `Γ(1/3)·3` оба `GammaValue`), у класса 9 — в
числителе (`∫₀¹dt/(1+t³)·1` и `·3` оба `IntegralValue`).

Преобразование Ландена (12.09, очередь стр. 77) — вторая ось связи
дополнительных модулей, и в отличие от Лежандра она НАПРАВЛЕННАЯ:
нисходящий шаг m → ((1−s)/(1+s))², s = √(1−m) = k′, и восходящий
m → 4t/(1+t)², t = √m [внешний источник: DLMF 19.8.11–19.8.12 и 19.8.16,
https://dlmf.nist.gov/19.8; Абрамовиц—Стиган 17.5.1, 17.5.3, 17.5.11,
стр. 597—598]. Корень входит и в слот образа, и в множители, потому
ТОЧНАЯ дорога (ℚ → ℚ) есть ровно там, где он рационален: вниз — когда
1 − m квадрат в ℚ, вверх — когда квадрат сам m. Вне этих точек
landen_descending/landen_ascending ОТКАЗЫВАЮТ поимённо, называя, какой
корень нерационален, — а величина образа несётся отрезком
(landen_slot_enclosure) и замеряется сертификатом (landen_enclosure)
монотонностью K и E по m. Замок двойственности: дополнение
перестановочно с Ланденом (1 − Ланден↓(m) = Ланден↑(1 − m) точно в ℚ), и
оттого K(m₁′)·K(m) = 2·K(m′)·K(m₁) — тождество с РАЦИОНАЛЬНЫМ
коэффициентом, где иррациональный множитель (1+s) гасится между
сторонами (landen_duality_enclosure). Один шаг АГМ и есть один
нисходящий шаг Ландена [внешний источник: Борвейн—Борвейн, SIAM Review
26(3):351—366, 1984, DOI 10.1137/1026073; «Pi and the AGM», Wiley 1987,
ISBN 9780471831389] — потому _agm_certificate ниже уже несёт Ландена
внутри, а эта дверь делает его ВИДИМЫМ и точным.

Соотношение Лежандра как ряд:
S_E(m)·S_K(m′) + S_E(m′)·S_K(m) − S_K(m)·S_K(m′) = 2/π — legendre_series.
Частное K/E и суммы с числом — вне направленных правил (частное рядов не
замкнуто; π-моном + число — не моном): отказ Unresolved с Ф-записями,
структурная запись — по ФIOS.
"""
from __future__ import annotations

from math import isqrt

from .rational import Q

from .enclosure import Interval, pi_enclosure, sqrt_enclosure
from .errors import Unresolved
from .intpoly import IntPoly
from .pi_value import PiMonomial
from .zero_state import Zero

#: Допустимые скаляры. float исключён построением.
_NUM = (int, Q)

#: `[вычислено здесь]` 10.10, очередь строка 105, шаг восьмой: дверь оператора
#: носителя = дверь ФIOS (как у π, e^y/ln, arctrig, ζ, Γ, корня J_ν). Для ЧУЖОГО
#: класса с записью операторы зовут `dispatch.combine`; класс 1 (int, Q) и
#: одноимённая линейка — прежние правила и собственные отказы носителя.
#: Ввоз ленивый, внутри функций.
_IN_DOOR: set = set()
_PINNED: frozenset = frozenset()


def _pinned(other, op: str) -> bool:
    from .Фios import classify, EXACT      # ленивый ввоз
    k = classify(other)
    return k is not None and k[0] == EXACT and (k[1], op) in _PINNED


def _door(a, b, op: str):
    from .dispatch import combine          # ленивый ввоз
    key = (id(a), id(b), op)
    if key in _IN_DOOR:
        return NotImplemented
    _IN_DOOR.add(key)
    try:
        return combine(a, b, op)
    finally:
        _IN_DOOR.discard(key)


def _foreign_record(x) -> bool:
    """Чужой класс с записью: не класс 1 и не одноимённая линейка."""
    return (not isinstance(x, (_NUM, _Elliptic))) and callable(getattr(x, "record", None))
#: Одна запись отказа на bool во всех дверях модуля: bool числом в доме не
#: считается и называется bool'ом (канон `phicore/rational.py`; сосед по
#: главе отвечает тем же: `Notch.from_q(True, 8)` → «bool не число»).
#: Двух копий одного канона здесь нет намеренно: разошедшиеся копии — уже
#: пройденная поломка (`affine.py:75`).
_NOT_BOOL = "bool не число"

#: Предел итераций AGM. Сходимость квадратичная — предел лишь страховка;
#: его исчерпание получает честный отказ, не молчание.
_MAX_AGM = 300

_MULT_HINT = "слот линейки — только m: множитель линейкой не несётся"
_ADD_HINT = ("эллиптический интеграл независим: сложение сведено "
             "только с нолём-состоянием")


def _rec(x: object) -> str:
    """Ф-запись операнда для сообщений отказа."""
    r = getattr(x, "record", None)
    if callable(r):
        return r()
    if isinstance(x, _NUM):
        return f"Φ(C=1 RATIONAL, value={Q(x)})"
    return f"Φ(?={x!r})"


def _param(m) -> Q:
    """Параметр m = k²: только int/Q; float — TypeError."""
    if isinstance(m, float):
        raise TypeError("m — только int/Q, float запрещён")
    if not isinstance(m, _NUM):
        raise TypeError("m — только int/Q")
    return Q(m)


def _check_bits(bits) -> None:
    """Запрос сетки: целое число двоичных разрядов ширины."""
    if isinstance(bits, bool):
        raise TypeError(f"bits — {_NOT_BOOL}; дверь принимает int")
    if not isinstance(bits, int):
        raise TypeError("bits — только int, float запрещён")
    if bits < 1:
        raise ValueError("bits ≥ 1: запрошенная ширина сетки есть 2^−bits")


def complementary(m) -> Q:
    """Дополнительный модуль m′ = 1 − m — ось двойственности линейки.

    Самодуальность в устройстве: инволюция m′′ = m держится тождеством
    1 − (1 − m) = m, точной арифметикой Q; неподвижная
    (самодуальная) точка — m = 1/2.
    """
    return 1 - _param(m)


# ---------- сертификаты сетки: AGM в интервальной арифметике ----------

def _sqrt_iv(iv: Interval, bits: int) -> Interval:
    """Вложение √ на отрезке: √ монотонен, корни концов обрамляют.

    Нижний конец при положительном входе отделяется от ноля уже́нием
    сетки — точный шаг, не округление.
    """
    lo_e = sqrt_enclosure(iv.lo, bits)
    while iv.lo > 0 and lo_e.lo <= 0:
        bits *= 2
        lo_e = sqrt_enclosure(iv.lo, bits)
    hi_e = sqrt_enclosure(iv.hi, bits)
    return Interval(lo_e.lo, hi_e.hi)


def _agm_certificate(m: Q, bits: int) -> tuple[Interval, Interval]:
    """Совместный ход AGM для m ∈ (0, 1): вложения AGM(1, √(1−m)) и
    суммы Σ_{n≥0} 2^{n−1}·c_n² с мажорантой хвоста.

    Инвариант вложения: истинные a_n, b_n лежат в своих отрезках, и
    b_n ≤ AGM ≤ a_n, потому AGM ∈ [b.lo, a.hi]. Итерация — до сходимости
    ширины: зазор a.hi − b.lo ≤ 2^−(bits+16) (квадратичная сходимость
    среднего). Хвост суммы: из тождества c_{n+1} = c_n²/(4·a_{n+1})
    [внешний источник: метод AGM, Абрамовиц—Стиган 17.6] и нижней грани
    B ≤ b_0 ≤ a_n следует c_{n+k} ≤ r^k·C при C ≥ |c_n|, r = C/(4B), и

        Σ_{k≥1} 2^{n+k−1}·c_{n+k}² ≤ 2^{n−1}·C²·2r²/(1 − 2r²),

    геометрическая мажоранта при 2r² < 1 — целиком в Q.
    """
    wbits = bits + 48                 # запас на накопление ширины
    sbits = wbits + 16                # сетка корней мельче рабочей
    eps = Q(1, 2 ** (bits + 16))
    b = sqrt_enclosure(1 - m, sbits)  # b_0 = √(1−m)
    while b.lo <= 0:                  # отделить b_0 от ноля уже́нием
        sbits *= 2
        b = sqrt_enclosure(1 - m, sbits)
    B = b.lo                          # нижняя грань всех a_n, b_n
    a = Interval.point(1)
    s = Interval.point(m / 2)         # член n = 0: 2^{−1}·c_0², c_0² = m точен
    pow2 = 1                          # 2^{n−1} на шаге n = 1
    for _ in range(_MAX_AGM):
        c = (a - b) * Q(1, 2)          # c_n = (a_{n−1} − b_{n−1})/2
        s = s + c * c * pow2                  # член 2^{n−1}·c_n²
        pow2 *= 2
        a, b = (a + b) * Q(1, 2), _sqrt_iv(a * b, sbits)
        gap = a.hi - b.lo
        big = max(abs(c.lo), abs(c.hi))       # верхняя грань |c_n|
        r = big / (4 * B)
        if gap <= eps and 2 * r * r < 1:
            # 2^{n−1} = pow2/2 после удвоения; мажоранта хвоста ≥ 0
            tail = pow2 * big * big * r * r / (1 - 2 * r * r)
            if tail <= eps:
                return Interval(b.lo, a.hi), s + Interval(0, tail)
    raise Unresolved("AGM", f"Φ(C=1 RATIONAL, value={m})",
                     hint=f"сертификат сетки не сошёлся за {_MAX_AGM} "
                          f"итераций — уточнение по явному запросу")


def _ke_certificates(m: Q, bits: int) -> tuple[Interval, Interval]:
    """Вложения K(m) и E(m) одним ходом AGM, m ∈ (0, 1); при m < 0 — через
    ПРЕОБРАЗОВАНИЕ МНИМОГО МОДУЛЯ (02.10) [внешний источник: Абрамовиц—Стиган
    17.4.17–17.4.18; DLMF 19.7.5, https://dlmf.nist.gov/19.7.E5]:
        K(−μ) = K(μ/(1+μ)) / √(1+μ),   E(−μ) = √(1+μ) · E(μ/(1+μ)),  μ > 0,
    где μ/(1+μ) ∈ (0, 1) — область доказанного AGM-сертификата; √(1+μ) —
    сертификатом корня. Прямой ход AGM при m < 0 здесь НЕ берётся: его
    мажоранта хвоста опирается на B = b₀ как нижнюю грань всех средних, а
    при b₀ = √(1−m) > 1 = a₀ эта грань ложна (a₁ < b₀) — замер 02.10."""
    if m < 0:
        mu = -m
        kp, ep = _ke_certificates(mu / (1 + mu), bits)
        root = sqrt_enclosure(1 + mu, bits + 16)
        return kp * root.inverse(), ep * root
    agm, s = _agm_certificate(m, bits)
    k = pi_enclosure(max(30, bits)) * agm.inverse() * Q(1, 2)
    return k, k * (1 - s)


def agm_enclosure(m, bits: int = 64) -> Interval:
    """Публичная AGM-дверь: сертификат M(1, √(1−m)) ⊂ ℚ, m ∈ (0, 1).

    Та же величина, что внутри `_agm_certificate` (первый элемент пары), —
    до 01.10 дверь была непубличной, и приёмка двойственности брала её
    через модуль (сдача 25.09, очередь строка 104). Связь с линейкой —
    K(m)·M(1, √(1−m)) = π/2 [внешний источник: DLMF 19.8.5]; проверяется
    приёмкой пересечением с `pi_enclosure`. На краях сведения нет: m = 1
    даёт b₀ = 0 (M = 0, K — полюс), m = 0 — M = 1 тривиально; оба — отказ
    по имени, как у `legendre_enclosure`.
    """
    m = _param(m)
    _check_bits(bits)
    if not 0 < m < 1:
        raise Unresolved("AGM", _rec(m),
                         hint="M(1, √(1−m)) сведено на m ∈ (0, 1): m = 1 — "
                              "b₀ = 0 и полюс K, m = 0 — M = 1 тривиально")
    return _agm_certificate(m, bits)[0]


def ke_enclosures(m, bits: int = 64) -> tuple[Interval, Interval]:
    """Публичная дверь: сертификаты K(m) и E(m) ОДНИМ ходом AGM, m ∈ (0, 1) —
    те же отрезки, что `KValue(m).enclosure(bits)` и `EValue(m).enclosure(bits)`
    порознь (приёмка двойственности C); до 01.10 — только через модуль."""
    m = _param(m)
    _check_bits(bits)
    if m == 0 or m >= 1:
        raise Unresolved("K, E", _rec(m),
                         hint="сертификаты одним ходом AGM сведены на m ∈ (0, 1) и, "
                              "преобразованием мнимого модуля, на m < 0; m = 0 — π/2 "
                              "точно, m = 1 — полюс K, m > 1 — K и E комплексны")
    return _ke_certificates(m, bits)


# ---------- сама линейка ----------

class _SeriesView:
    """Гипергеометрическая часть S = X/(π/2) линейки X ∈ {K, E}: лист ряда
    класса 7, чей сертификат — от самой линейки (AGM) делением на π/2."""

    __slots__ = ("x",)

    def __init__(self, x):
        self.x = x

    def enclosure(self, N: int = 64) -> Interval:
        bits = max(64, int(N))
        half_pi = pi_enclosure(bits) * Interval.point(Q(1, 2))
        return self.x.enclosure(bits) / half_pi

    def record(self) -> str:
        return f"{self.x.record()} / (π/2)"

    def __eq__(self, other):
        return isinstance(other, _SeriesView) and self.x == other.x

    def __hash__(self):
        return hash(("_SeriesView", self.x))


def _lift_series(other):
    """Поднятие второго операнда к D-конечному ряду (ряд, интеграл, произведение); None — не поднимается."""
    from .dfinite import DFiniteSeries
    if isinstance(other, (bool, float)) or isinstance(other, _NUM):
        return None
    return DFiniteSeries._lift(other)


class _Elliptic:
    """Общий носитель слота m ∈ (0, 1); буква линейки — у наследников.

    Края и внешность области сведены фабриками elliptic_k/elliptic_e,
    не конструктором: канон держится построением.
    """

    __slots__ = ("m",)

    #: Буква линейки в Ф-записи: "K" | "E" (задаёт наследник).
    _KIND = "?"

    def __init__(self, m):
        m = _param(m)
        if m == 0 or m >= 1:
            raise ValueError("канонический слот требует m < 1, m ≠ 0 (m < 0 — мнимый модуль, 02.10); края "
                             "и внешность области — фабриками "
                             "elliptic_k/elliptic_e")
        self.m = m

    # --- двойственность с классом 7: гипергеометрическая часть -----------

    def series(self):
        """S = X/(π/2) как D-конечный ряд класса 7 (порядок 1); сертификат — от линейки."""
        from .dfinite import DFiniteSeries
        from .series import PowerSeriesValue
        num = IntPoly([1, 4, 4]) if self._KIND == "K" else IntPoly([-1, 0, 4])
        ps = PowerSeriesValue(Q(1), num, IntPoly([4, 8, 4]), self.m)
        return DFiniteSeries(DFiniteSeries.from_power(ps).seq, ("leaf", _SeriesView(self)))

    # --- самопроверка носителя (ось самоприменения, 08.10) ----------------

    def self_check(self, deep: bool = False) -> list:
        """Самопроверка носителя встречными дорогами (ось «AIAΣA спрашивает
        AIAΣA», слово автора 04.09; очередь строка 18; 08.10). Список
        находок, пусто — здоров; исключений не бросает.

        `[вычислено здесь]` (1) слот: m — ℚ, m < 1, m ≠ 0 (края и m > 1 живут
        у фабрик); (2) сертификат enclosure(64): концы упорядочены, ширина
        ≤ 2^−32, при m ∈ (0, 1) K ≥ π/2 и 1 ≤ E ≤ π/2 [внешний источник:
        Абрамовиц—Стиган 17.3]; (3) ВСТРЕЧНЫЕ дороги при m ∈ (0, 1), все
        независимые от AGM по тождеству: ряд класса 7 (`series`) — у K члены
        положительны, потому (π/2)·S_N снизу под сертификатом, у E все члены
        после первого отрицательны, потому (π/2)·S_N сверху над ним
        [А—С 17.3.9–17.3.12]; соотношение Лежандра по дополнительному
        модулю несёт π/2 [А—С 17.3.13]; невязки шага Ландена вниз несут ноль
        [DLMF 19.8.11–19.8.12]. При m < 0 (мнимый модуль) дороги (3) не
        строятся — граница объявлена. При deep — ряд до N = 48 и Ланден
        вверх."""
        out = []
        m = self.m
        if not isinstance(m, Q):
            return [f"слот m чужого носителя: {type(m).__name__}"]
        if m == 0 or m >= 1:
            return [f"слот m = {m} вне канона (0 и m ≥ 1 — у фабрик)"]
        try:
            iv = self.enclosure(64)
        except Exception as e:                                  # noqa: BLE001
            return out + [f"сертификат не построен: {type(e).__name__}: {e}"]
        if iv.lo > iv.hi:
            out.append(f"концы сертификата переставлены: {iv.record()}")
        if iv.width() > Q(1, 2 ** 32):
            out.append(f"ширина сертификата {iv.width()} > 2^−32")
        if m < 0:
            return out
        try:
            half_pi = pi_enclosure(40) / 2
            if self._KIND == "K" and iv.hi < half_pi.lo:
                out.append(f"K({m}) ниже π/2: {iv.record()}")
            if self._KIND == "E" and (iv.lo > half_pi.hi or iv.hi < 1):
                out.append(f"E({m}) вне [1, π/2]: {iv.record()}")
            N = 48 if deep else 24
            s_n = self.series().partial_sum(N)
            if self._KIND == "K" and half_pi.lo * s_n > iv.hi:
                out.append(f"ряд класса 7: (π/2)·S_{N} = {half_pi.lo * s_n} выше сертификата K({m}) {iv.record()}")
            if self._KIND == "E" and half_pi.hi * s_n < iv.lo:
                out.append(f"ряд класса 7: (π/2)·S_{N} = {half_pi.hi * s_n} ниже сертификата E({m}) {iv.record()}")
            leg = legendre_enclosure(m, 64)
            if leg.intersect(half_pi) is None:
                out.append(f"соотношение Лежандра при m = {m} не несёт π/2: {leg.record()}")
            kinds = (LANDEN_DOWN, LANDEN_UP) if deep else (LANDEN_DOWN,)
            for kind in kinds:
                rk, re_ = landen_enclosure(m, 64, kind)
                if not rk.contains(Q(0)) or not re_.contains(Q(0)):
                    out.append(f"невязка Ландена {kind} при m = {m} не несёт ноль: K {rk.record()}, E {re_.record()}")
        except Unresolved as e:
            out.append(f"встречная дорога отказала при m = {m}: {e}")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"встречная дорога не построена при m = {m}: {type(e).__name__}: {e}")
        return out

    @staticmethod
    def _pi_times(coeff, e, s):
        """π-моном × D-конечный ряд — структурная запись (класс 19): значение — строение, засечка по запросу."""
        from .gridcore import Mixed
        return Mixed.combine(PiMonomial(Q(coeff), Q(e)), s, "*")

    # --- сложение / вычитание ------------------------------------------

    def _refuse(self, op: str, other: object, self_left: bool, hint: str):
        left, right = ((self.record(), _rec(other)) if self_left
                       else (_rec(other), self.record()))
        raise Unresolved(op, left, right, hint=hint)

    def _add_like(self, other, sign: int, self_left: bool):
        if isinstance(other, _Elliptic):
            s = self.series() + other.series() if sign > 0 else (self.series() - other.series() if self_left else other.series() - self.series())
            return self._pi_times(Q(1, 2), 1, s) if not isinstance(s, Q) else PiMonomial(Q(1, 2), 1) * s if s != 0 else Q(0)
        op = "+" if sign > 0 else "-"
        if _foreign_record(other) and not _pinned(other, op):
            a, b = (self, other) if self_left else (other, self)
            return _door(a, b, op)              # чужой класс — дверь ФIOS
        self._refuse(op, other, self_left, _ADD_HINT)

    def __add__(self, other):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__radd__
        return self._add_like(other, 1, True)

    def __radd__(self, other):
        return self._add_like(other, 1, False)

    def __sub__(self, other):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__rsub__
        return self._add_like(other, -1, True)

    def __rsub__(self, other):
        return self._add_like(other, -1, False)

    def __neg__(self):
        return self._pi_times(Q(-1, 2), 1, self.series())

    # --- умножение / деление -------------------------------------------

    def __mul__(self, other):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__rmul__
        if isinstance(other, _NUM) and not isinstance(other, bool):
            q = Q(other)
            if q == 0:
                return Q(0)  # масштаб нолём — точный ноль класса 1
            if q == 1:
                return self         # единица — тождество
            return self._pi_times(q / 2, 1, self.series())
        if isinstance(other, _Elliptic):
            return self._pi_times(Q(1, 4), 2, self.series() * other.series())
        lifted = _lift_series(other)
        if lifted is not None:
            return self._pi_times(Q(1, 2), 1, self.series() * lifted)
        if _foreign_record(other) and not _pinned(other, "*"):
            return _door(self, other, "*")      # чужой класс — дверь ФIOS
        r = getattr(other, "record", None)
        if callable(r):
            raise Unresolved("*", self.record(), r(),
                             hint="правила сведения с этой линейкой нет")
        return NotImplemented

    def __rmul__(self, other):
        """Скаляр, ряд и одноимённая линейка — коммутативно, тем же правилом;
        чужой класс — дверь с порядком операндов как в выражении."""
        if _foreign_record(other) and _lift_series(other) is None and not _pinned(other, "*"):
            return _door(other, self, "*")
        return self.__mul__(other)

    def __truediv__(self, other):
        if isinstance(other, _NUM) and not isinstance(other, bool):
            q = Q(other)
            if q == 0:
                raise Unresolved("/", self.record(), _rec(other),
                                 hint="x/0 здесь не сведено")
            if q == 1:
                return self
            return self._pi_times(1 / (2 * q), 1, self.series())
        if isinstance(other, Zero) or (_foreign_record(other) and not _pinned(other, "/")):
            return _door(self, other, "/")      # структурная запись — по ФIOS, как и сказано
        r = getattr(other, "record", None)
        if callable(r):
            raise Unresolved("/", self.record(), r(),
                             hint="частное рядов не замкнуто (теорема): структурная запись — по ФIOS")
        return NotImplemented

    def __rtruediv__(self, other):
        if isinstance(other, _NUM):
            raise Unresolved("/", _rec(other), self.record(),
                             hint="обратная величина этой линейкой "
                                  "не несётся (обращение ряда не замкнуто)")
        if _foreign_record(other) and not _pinned(other, "/"):
            return _door(other, self, "/")
        return NotImplemented

    def __pow__(self, n):
        if not isinstance(n, int) or isinstance(n, bool):
            return NotImplemented
        if n == 0:
            return Q(1)
        if n == 1:
            return self
        if n >= 2:
            return self._pi_times(Q(1, 2 ** n), n, self.series() ** n)
        raise Unresolved("^", self.record(), _rec(n),
                         hint="отрицательная степень — обращение ряда, не замкнуто")

    # --- запись и сравнение --------------------------------------------

    def record(self) -> str:
        return f"Φ(C=15 ELLIPTIC {self._KIND}, m={self.m})"

    __repr__ = record

    def __eq__(self, other):
        if isinstance(other, _Elliptic):
            # K строго возрастает, E строго убывает по m, а их области
            # значений на (0, 1) разделены точкой π/2: равенство величин
            # есть совпадение буквы и слота
            return type(self) is type(other) and self.m == other.m
        # равенство с рациональным/алгебраическим этой линейкой не
        # решается — ответ не выдумывается
        return NotImplemented

    def __hash__(self):
        return hash((type(self), self.m))


class KValue(_Elliptic):
    """K(m) — полный интеграл первого рода в каноническом слоте m."""

    __slots__ = ()
    _KIND = "K"

    def enclosure(self, bits: int = 64) -> Interval:
        """Сертификат сетки [lo, hi] ⊂ ℚ: K(m) = π/(2·AGM(1, √(1−m)));
        уточнение — явным запросом (рост bits)."""
        _check_bits(bits)
        return _ke_certificates(self.m, bits)[0]


class EValue(_Elliptic):
    """E(m) — полный интеграл второго рода в каноническом слоте m."""

    __slots__ = ()
    _KIND = "E"

    def enclosure(self, bits: int = 64) -> Interval:
        """Сертификат сетки [lo, hi] ⊂ ℚ: E(m) = K(m)·(1 − Σ 2^{n−1}·c_n²)
        с мажорантой хвоста; уточнение — явным запросом (рост bits)."""
        _check_bits(bits)
        return _ke_certificates(self.m, bits)[1]


# ---------- фабрики: приземление краёв, канон внутри ----------

def elliptic_k(m) -> PiMonomial | KValue:
    """K(m) по рациональному параметру m = k².

    m = 0 → π/2 — приземление в класс 10; m = 1 — расходимость
    (K(m) → ∞), честный отказ; m < 0 — носитель KValue, сертификат
    преобразованием мнимого модуля (02.10; до того отказ «не
    реализовано»); m > 1 — отказ: k² > 1, K(m) комплексен [внешний
    источник: Абрамовиц—Стиган 17.4.15; DLMF 19.7.3], вещественный слот
    линейки его не несёт; иначе — носитель KValue.
    """
    m = _param(m)
    if m == 0:
        return PiMonomial(Q(1, 2), Q(1))
    if m == 1:
        raise Unresolved("K", _rec(m),
                         hint="расходимость: K(m) → ∞ при m → 1, "
                              "значения нет")
    if m > 1:
        raise Unresolved("K", _rec(m),
                         hint="m > 1 — вне [0, 1]: k² > 1, K(m) комплексен "
                              "(обратный модуль, А—С 17.4.15), вещественный "
                              "слот его не несёт")
    return KValue(m)


def elliptic_e(m) -> PiMonomial | Q | EValue:
    """E(m) по рациональному параметру m = k².

    m = 0 → π/2 — приземление в класс 10; m = 1 → точная единица
    (вырождение в класс 1); m < 0 — носитель EValue (мнимый модуль,
    02.10); m > 1 — отказ: E(m) комплексен; иначе — носитель EValue.
    """
    m = _param(m)
    if m == 0:
        return PiMonomial(Q(1, 2), Q(1))
    if m == 1:
        return Q(1)
    if m > 1:
        raise Unresolved("E", _rec(m),
                         hint="m > 1 — вне [0, 1]: k² > 1, E(m) комплексен "
                              "(обратный модуль, А—С 17.4.15), вещественный "
                              "слот его не несёт")
    return EValue(m)


# ---------- двойственность 7↔15: соотношение Лежандра рядом ----------

def legendre_series(m):
    """Соотношение Лежандра рядами: S_E(m)·S_K(m′) + S_E(m′)·S_K(m) − S_K(m)·S_K(m′) = 2/π.

    Три произведения — D-конечные ряды (свёртка Коши замыканием, порядок 6–7
    каждый, доли секунды); их сумма держится СТРУКТУРНОЙ записью (класс 19):
    замыкание суммы через ОДУ даёт рекуррентту порядка 22 за минуты (замер
    05.09: 170 с при m = 1/3) — минимальной рекурренты dfinite.py не ищет,
    и это названо в очереди, не спрятано. Сертификат — по частям (от AGM
    линеек), частные суммы — legendre_partial. m ∈ (0, 1)."""
    from .gridcore import Mixed
    p1, p2, p3 = _legendre_products(m)
    return Mixed.combine(Mixed.combine(p1, p2, "+"), p3, "-")


#: Кэш произведений ряда Лежандра по слоту m. Замер 25.09 (сдача эллиптики,
#: очередь строка 104): `legendre_partial` пересобирал замыкание трёх свёрток
#: Коши при КАЖДОМ вызове — ≈1,3 с при m = 1/3, ≈0,7 с при 1/2; частные суммы
#: S_N по одному и тому же замыканию считаются из него, замыкание — один раз.
_LEGENDRE_CACHE: dict = {}


def _legendre_products(m):
    m = _param(m)
    if not 0 < m < 1:
        raise Unresolved("legendre_series", _rec(m), hint="требует m ∈ (0, 1): на краях K — полюс")
    got = _LEGENDRE_CACHE.get(m)
    if got is not None:
        return got
    mc = complementary(m)
    sk, se = KValue(m).series(), EValue(m).series()
    skc, sec = KValue(mc).series(), EValue(mc).series()
    got = (se * skc, sec * sk, sk * skc)
    _LEGENDRE_CACHE[m] = got
    return got


def legendre_partial(m, N: int) -> Q:
    """Частная сумма ряда Лежандра до N: S_N(p1) + S_N(p2) − S_N(p3) → 2/π."""
    p1, p2, p3 = _legendre_products(m)
    return p1.partial_sum(N) + p2.partial_sum(N) - p3.partial_sum(N)


# ---------- двойственность: соотношение Лежандра сертификатом ----------

def legendre_enclosure(m, bits: int = 64) -> Interval:
    """Сертификат левой части соотношения Лежандра по дополнительному
    модулю m′ = 1 − m:

        E(m)·K(m′) + E(m′)·K(m) − K(m)·K(m′)

    [внешний источник: Лежандр, Абрамовиц—Стиган 17.3.13]; отрезок
    обязан нести π/2 — сверка с pi_enclosure и есть замер двойственности.
    В самодуальной точке m = 1/2 выражение стягивается в
    2·E(1/2)·K(1/2) − K(1/2)². Требует m ∈ (0, 1): на краях K — полюс.
    """
    m = _param(m)
    _check_bits(bits)
    if not 0 < m < 1:
        raise Unresolved("Лежандр", _rec(m),
                         hint="соотношение сведено на m ∈ (0, 1): "
                              "на краю K(m) или K(m′) — полюс")
    k1, e1 = _ke_certificates(m, bits)
    k2, e2 = _ke_certificates(1 - m, bits)
    return e1 * k2 + e2 * k1 - k1 * k2


# ---------- преобразование Ландена: точная дорога и сертификатная ----------

#: Род шага Ландена. Вниз — модуль убывает (k₁ < k), вверх — растёт.
LANDEN_DOWN = "вниз"
LANDEN_UP = "вверх"

#: Сколько раз сетка удваивается, пока отрезок образа отделяется от краёв
#: 0 и 1. Исчерпание получает честный отказ, не молчание.
_MAX_NARROW = 8


def _exact_root(q: Q):
    """Точный рациональный √q или None. Ни одного float: целый isqrt.

    Канон Q (d > 0, gcd(n, d) = 1) делает проверку ПОЛНОЙ: у несократимой
    пары корень рационален тогда и только тогда, когда квадратны обе
    части (иначе в n/d остался бы неквадратный множитель, а сократить
    его не с чем). Это теорема, а не перебор.
    """
    if q < 0:
        return None
    n, d = q.numerator, q.denominator
    rn, rd = isqrt(n), isqrt(d)
    if rn * rn != n or rd * rd != d:
        return None
    return Q(rn, rd)


def _check_kind(kind) -> str:
    """Род шага: только два слова, третьего нет."""
    if not isinstance(kind, str):
        raise TypeError("род шага Ландена — строка "
                        f"{LANDEN_DOWN!r} либо {LANDEN_UP!r}")
    if kind not in (LANDEN_DOWN, LANDEN_UP):
        raise ValueError(f"род шага Ландена — {LANDEN_DOWN!r} (модуль "
                         f"убывает) либо {LANDEN_UP!r} (модуль растёт); "
                         f"получено {kind!r}")
    return kind


class LandenStep:
    """Точный шаг Ландена: слот → слот и три коэффициента, все в ℚ.

        K(m₁) = k_mult·K(m),   E(m₁) = e_self·E(m) + e_k·K(m).

    Слот образа и коэффициенты рациональны ПОСТРОЕНИЕМ: шаг строится лишь
    там, где нужный корень (√(1−m) вниз, √m вверх) точно рационален; иначе
    фабрика отказывает поимённо. float здесь невозможен.
    """

    __slots__ = ("kind", "m", "image", "root", "k_mult", "e_self", "e_k")

    def __init__(self, kind, m, image, root, k_mult, e_self, e_k):
        self.kind, self.m, self.image, self.root = kind, m, image, root
        self.k_mult, self.e_self, self.e_k = k_mult, e_self, e_k

    def record(self) -> str:
        знак = "−" if self.e_k < 0 else "+"
        модуль = -self.e_k if self.e_k < 0 else self.e_k
        return (f"Φ(Ланден {self.kind}, m={self.m} → m₁={self.image}, "
                f"корень={self.root}, K×{self.k_mult}, "
                f"E={self.e_self}·E {знак} {модуль}·K)")

    __repr__ = record

    def __eq__(self, other):
        if not isinstance(other, LandenStep):
            return NotImplemented
        return (self.kind, self.m, self.image, self.root, self.k_mult,
                self.e_self, self.e_k) == (other.kind, other.m, other.image,
                                           other.root, other.k_mult,
                                           other.e_self, other.e_k)

    def __hash__(self):
        return hash((LandenStep, self.kind, self.m, self.image, self.root,
                     self.k_mult, self.e_self, self.e_k))


def landen_descending(m) -> LandenStep:
    """Нисходящий Ланден: m → m₁ = ((1−s)/(1+s))², s = √(1−m) = k′.

        K(m₁) = ((1+s)/2)·K(m),   E(m₁) = (E(m) + s·K(m))/(1+s)

    [внешний источник: DLMF 19.8.11–19.8.12, https://dlmf.nist.gov/19.8 —
    k₁ = (1−k′)/(1+k′), K(k) = (1+k₁)K(k₁), E(k) = (1+k′)E(k₁) − k′K(k);
    то же в модульных углах — Абрамовиц—Стиган 17.5.1, 17.5.3, стр. 597—598].
    Здесь обе строки решены относительно ОБРАЗА: 1 + k₁ = 2/(1+s).

    Точная дорога существует ровно тогда, когда s рационален, то есть
    m′ = 1 − m есть квадрат в ℚ: m = (q²−p²)/q² даёт m₁ = ((q−p)/(q+p))²,
    k_mult = (q+p)/(2q), E(m₁) = (q·E(m) + p·K(m))/(q+p). Иначе — отказ
    поимённо; десятичной дороги эта дверь не открывает.
    """
    m = _param(m)
    if not 0 < m < 1:
        raise Unresolved("Ланден↓", _rec(m),
                         hint="шаг сведён на m ∈ (0, 1): на краю K — полюс "
                              "либо вырождение")
    s = _exact_root(1 - m)
    if s is None:
        raise Unresolved("Ланден↓", _rec(m),
                         hint=f"точной дороги нет: √(1−m) = √({1 - m}) не "
                              f"рационален (1−m — не квадрат в ℚ), потому "
                              f"образ слота ((1−√(1−m))/(1+√(1−m)))² "
                              f"иррационален и слотом линейки 15 не "
                              f"несётся; сертификатная дорога — "
                              f"landen_enclosure(m, bits, 'вниз')")
    up = 1 + s
    r = (1 - s) / up
    return LandenStep(LANDEN_DOWN, m, r * r, s, up / 2, 1 / up, s / up)


def landen_ascending(m) -> LandenStep:
    """Восходящий Ланден: m → m₁ = 4t/(1+t)², t = √m = k.

        K(m₁) = (1+t)·K(m),   E(m₁) = (2·E(m) − (1−m)·K(m))/(1+t)

    [внешний источник: DLMF 19.8.16, https://dlmf.nist.gov/19.8 —
    k₂ = 2√k/(1+k); Абрамовиц—Стиган 17.5.11, стр. 598 —
    (1+sin αₙ)(1+cos αₙ₊₁) = 2]. ВАЖНО и замерено: пара DLMF 19.8.16—19.8.17
    записана для НЕПОЛНОГО интеграла и полноту не сохраняет — при φ = π/2
    выходит 2φ₂ = π/2 + arcsin k ≠ π, то есть прямой строки «K через K»
    там нет. Полные формулы выше получены ОБРАЩЕНИЕМ 19.8.11–19.8.12:
    нисходящий шаг от k₂ возвращает ровно k ((1−k₂′)/(1+k₂′) = k при
    k₂′ = (1−k)/(1+k)), потому 19.8.12 читается справа налево. Вторая
    запись того же E-коэффициента: E(m₁) = 2E(m)/(1+t) − (1−t)·K(m)
    (равенство (1−m)/(1+t) = 1−t точно в ℚ).

    Точная дорога существует ровно тогда, когда t рационален, то есть сам
    m есть квадрат в ℚ: m = p²/q² даёт m₁ = 4pq/(q+p)², k_mult = (q+p)/q.
    Иначе — отказ поимённо.
    """
    m = _param(m)
    if not 0 < m < 1:
        raise Unresolved("Ланден↑", _rec(m),
                         hint="шаг сведён на m ∈ (0, 1): на краю K — полюс "
                              "либо вырождение")
    t = _exact_root(m)
    if t is None:
        raise Unresolved("Ланден↑", _rec(m),
                         hint=f"точной дороги нет: √m = √({m}) не рационален "
                              f"(m — не квадрат в ℚ), потому образ слота "
                              f"4√m/(1+√m)² иррационален и слотом линейки 15 "
                              f"не несётся; сертификатная дорога — "
                              f"landen_enclosure(m, bits, 'вверх')")
    up = 1 + t
    return LandenStep(LANDEN_UP, m, 4 * t / (up * up), t,
                      up, 2 / up, -(1 - t))


def landen(m, kind: str = LANDEN_DOWN) -> LandenStep:
    """Точный шаг Ландена названного рода — одна дверь на две фабрики."""
    _check_kind(kind)
    return (landen_descending(m) if kind == LANDEN_DOWN
            else landen_ascending(m))


# ---------- сертификатная дорога: корень иррационален ----------

def _root_iv(x: Q, bits: int) -> Interval:
    """Вложение √x при x > 0, отделённое от ноля уже́нием сетки."""
    e = sqrt_enclosure(x, bits)
    while e.lo <= 0:
        bits *= 2
        e = sqrt_enclosure(x, bits)
    return e


def landen_slot_enclosure(m, bits: int = 64,
                          kind: str = LANDEN_DOWN) -> Interval:
    """Отрезок [lo, hi] ⊂ ℚ ∩ (0, 1), несущий образ слота Ландена.

    Там, где корень иррационален, ТОЧНОЙ дороги нет (landen_descending
    отказывает), но величина образа есть, и её несёт отрезок с
    рациональными концами: √ берётся sqrt_enclosure, дальше — интервальная
    арифметика. Концы приводятся к сетке 2^−bits наружу, чтобы знаменатели
    не росли от шага к шагу. Отрезок обязан отделиться от обоих краёв
    (0 и 1): иначе K на конце не построить, и дверь отказывает поимённо
    после _MAX_NARROW удвоений сетки.
    """
    m = _param(m)
    _check_bits(bits)
    _check_kind(kind)
    if not 0 < m < 1:
        raise Unresolved("Ланден-слот", _rec(m),
                         hint="шаг сведён на m ∈ (0, 1)")
    b = bits
    for _ in range(_MAX_NARROW):
        if kind == LANDEN_DOWN:
            s = _root_iv(1 - m, b)
            if s.hi < 1:
                r = (1 - s) / (1 + s)
                iv = r * r
            else:
                iv = None
        else:
            t = _root_iv(m, b)
            iv = (t * 4) / ((1 + t) * (1 + t))
        if iv is not None:
            iv = Interval.on_grid(iv.lo, iv.hi, b)
            if iv.lo > 0 and iv.hi < 1:
                return iv
        b *= 2
    raise Unresolved("Ланден-слот", _rec(m),
                     hint=f"отрезок образа не отделён от краёв (0, 1) за "
                          f"{_MAX_NARROW} удвоений сетки — уточнение по "
                          f"явному запросу (рост bits)")


def _k_on_slot(iv: Interval, bits: int) -> Interval:
    """K на отрезке слота: K строго возрастает по m, концы обрамляют."""
    return Interval(KValue(iv.lo).enclosure(bits).lo,
                    KValue(iv.hi).enclosure(bits).hi)


def _e_on_slot(iv: Interval, bits: int) -> Interval:
    """E на отрезке слота: E строго убывает по m, концы меняются местами."""
    return Interval(EValue(iv.hi).enclosure(bits).lo,
                    EValue(iv.lo).enclosure(bits).hi)


def landen_enclosure(m, bits: int = 64,
                     kind: str = LANDEN_DOWN) -> tuple[Interval, Interval]:
    """Невязки шага Ландена сертификатами сетки: (K, E).

    Каждый отрезок ОБЯЗАН нести ноль — это и есть замер тождества числами:

        K(m₁) − k_mult·K(m) ∋ 0,   E(m₁) − (e_self·E(m) + e_k·K(m)) ∋ 0,

    где m₁ несётся отрезком (landen_slot_enclosure), а K и E на нём —
    монотонностью по m. Работает при ВСЯКОМ рациональном m ∈ (0, 1), в том
    числе там, где точной дороги нет: иррациональны здесь только концы
    отрезков, ни одного float. Ширина невязки называет цену замера.
    """
    m = _param(m)
    _check_bits(bits)
    _check_kind(kind)
    slot = landen_slot_enclosure(m, bits, kind)
    km, em = _ke_certificates(m, bits)
    if kind == LANDEN_DOWN:
        s = _root_iv(1 - m, bits)
        pred_k = (s + 1) * km * Q(1, 2)
        pred_e = (em + km * s) / (s + 1)
    else:
        t = _root_iv(m, bits)
        pred_k = (t + 1) * km
        pred_e = (em * 2 - km * (1 - m)) / (t + 1)
    return _k_on_slot(slot, bits) - pred_k, _e_on_slot(slot, bits) - pred_e


def landen_duality_enclosure(m, bits: int = 64) -> tuple[Interval, Interval]:
    """Двойственность дополнительного модуля ЧИСЛАМИ, вне самодуальной точки.

    Нисходящий шаг от m и ВОСХОДЯЩИЙ от m′ = 1 − m берут ОДИН И ТОТ ЖЕ
    корень s = √m′ и дают дополнительные друг другу слоты:

        m₁ = ((1−s)/(1+s))²,   1 − m₁ = 4s/(1+s)² = Ланден↑(m′),

    то есть дополнение перестановочно с Ланденом (точное тождество в ℚ).
    Отсюда величинами:

        K(m₁′) = (1+s)·K(m′),   K(m₁) = ((1+s)/2)·K(m)
        ⟹  K(m₁′)·K(m) = 2·K(m′)·K(m₁).

    Иррациональный множитель (1+s) гасится между сторонами: коэффициент
    тождества РАЦИОНАЛЕН (двойка), и сравнивать нечего, кроме самих
    величин. Возвращаются две стороны; они обязаны пересечься.

    Что именно даёт выход из самодуальной точки (замер 12.09): при m = 1/2
    множители K(m) и K(m′) суть ОДНО число и сокращаются, тождество
    стягивается в K(m₁′) = 2·K(m₁) — утверждение верное и нетривиальное
    (второе сингулярное значение), но связи ДВУХ РАЗНЫХ K в нём уже нет.
    Тавтологией оно не становится — в отличие от равенства
    K(complementary(1/2)) == K(1/2), которое есть утверждение об инволюции
    слота и величины не меряет вовсе. При m ≠ 1/2 в тождестве стоят четыре
    попарно различные величины, и дополнительный модуль несёт вес.
    """
    m = _param(m)
    _check_bits(bits)
    if not 0 < m < 1:
        raise Unresolved("Ланден-двойственность", _rec(m),
                         hint="сведено на m ∈ (0, 1): на краю K — полюс")
    mc = 1 - m
    k1 = _k_on_slot(landen_slot_enclosure(m, bits, LANDEN_DOWN), bits)
    k1c = _k_on_slot(landen_slot_enclosure(mc, bits, LANDEN_UP), bits)
    km = _ke_certificates(m, bits)[0]
    kmc = _ke_certificates(mc, bits)[0]
    return k1c * km, k1 * kmc * 2


# ---------- самопроверка линейки (ось самоприменения, 08.10) ----------

def self_check(deep: bool = False) -> list:
    """Самопроверка линейки класса 15 встречными дорогами; список находок,
    пусто — здорова; исключений не бросает.

    `[вычислено здесь]` Окно: m ∈ {1/2, 1/3, 3/4, 1/10} (deep — ещё 9/10 и
    мнимый модуль −1) для K и E — носители проверяют себя (`self_check`:
    ряд класса 7, Лежандр, Ланден); края фабрик: K(0) = E(0) = π/2 мономом
    класса 10, E(1) = 1, K(1) — отказ; инволюция дополнительного модуля
    m′′ = m; самодуальная точка — complementary(1/2) = 1/2. Граница — окно."""
    out = []
    pts = [Q(1, 2), Q(1, 3), Q(3, 4), Q(1, 10)]
    if deep:
        pts += [Q(9, 10), Q(-1)]
    for m in pts:
        for fab, letter in ((elliptic_k, "K"), (elliptic_e, "E")):
            try:
                v = fab(m)
                if isinstance(v, _Elliptic):
                    out += [f"{letter}({m}): {f}" for f in v.self_check(deep)]
                else:
                    out.append(f"{letter}({m}) приземлилось в {v!r} без носителя")
            except Exception as e:                              # noqa: BLE001
                out.append(f"{letter}({m}): {type(e).__name__}: {e}")
    try:
        hp = PiMonomial(Q(1, 2), Q(1))
        if elliptic_k(0) != hp or elliptic_e(0) != hp:
            out.append(f"край m = 0: K = {elliptic_k(0)!r}, E = {elliptic_e(0)!r} ≠ π/2")
        if elliptic_e(1) != Q(1):
            out.append(f"край m = 1: E = {elliptic_e(1)!r} ≠ 1")
        try:
            elliptic_k(1)
            out.append("край m = 1: K не отказал на расходимости")
        except Unresolved:
            pass
        for m in (Q(1, 3), Q(2, 7)):
            if complementary(complementary(m)) != m:
                out.append(f"инволюция m′′ ≠ m при m = {m}")
        if complementary(Q(1, 2)) != Q(1, 2):
            out.append("самодуальная точка 1/2 не неподвижна")
    except Exception as e:                                      # noqa: BLE001
        out.append(f"края и инволюция: {type(e).__name__}: {e}")
    return out
