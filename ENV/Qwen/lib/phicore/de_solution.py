"""Линейка класса 18 DE_SOLUTION: решение линейного ОДУ y′ = P(x)·y + Q(x).

Величина хранится структурно: многочлены P, Q ∈ ℚ[x] (коэффициенты
млад→ст) и начальное условие y(0) = y0 ∈ ℚ. Такая запись определяет
решение однозначно (задача Коши линейного уравнения), и первые N
коэффициентов Тейлора восстанавливаются точной рекуррентой

    a_0 = y0,   (n+1)·a_{n+1} = Σ_i P_i·a_{n−i} + Q_n

(свёртка по уже найденным a; Q_n = 0 вне степени Q) — счётом в Q,
без единого округления.

ДВОЙСТВЕННОСТЬ линейки — проверка подстановкой (структурная, точная):
check_substitution(N) собирает правую часть P·y + Q независимо от
рекурренты — полным произведением многочленов по уже полученным a — и
сверяет с левой частью (n+1)·a_{n+1} при всех n < N; равенство обязано
держаться до последнего бита.

Правил сведения решения с независимыми величинами нет — сумма с числом,
произведение, частное, степень вне {0, 1} получают отказ Unresolved с
явными Ф-записями обоих операндов, не ложный класс и не молчание.
Нулевая степень — единица по канону пустого произведения, одному на дом
(разбор — в докстринге `__pow__`; правка 12.09 по строке 75 очереди).
Единственная сведённая унарная операция — смена знака: по линейности
уравнения −y решает y′ = P·y − Q, y(0) = −y0.

Десятичная сторона — только сертификатом сетки (enclosure), по явному
запросу, и только там, где точный ответ известен самой линейке
(01.10, очередь строка 104; до того — лишь P = [c], Q = []):
  * Q = [] при любом P: y = y0·e^{F(t)}, F = ∫₀ᵗ P — многочлен в ℚ,
    значение F(t) точно, сертификат — y0·exp_enclosure(F(t));
  * P = [] при любом Q: y = y0 + ∫₀ᵗ Q — точная точка класса 1;
  * P = [c], Q = [q]: y = (y0 + q/c)·e^{c·t} − q/c — композиция дверей.
Общий случай (P непостоянен и Q ≠ [], либо deg Q ≥ 1 при P ≠ []) —
с 01.10 (вечер) сертификат ОЦЕНКОЙ КОШИ (`cauchy_enclosure`): решение
y — целая функция; на круге |z| ≤ ρ вдоль луча |y|′ ≤ P̃(ρ)·|y| + Q̃(ρ),
P̃(ρ) = Σ|P_i|·ρ^i, и по Грёнуоллу |y(z)| ≤ M(ρ) = (|y0| + ρ·Q̃(ρ))·e^{ρ·P̃(ρ)};
по Коши |a_n| ≤ M(ρ)/ρⁿ, и хвост Σ_{n≥N} a_n·tⁿ при |t| < ρ не превосходит
M(ρ)·(|t|/ρ)^N/(1 − |t|/ρ). Частная сумма — точно в ℚ; радиус ρ берётся
лучшим из восьми кратных |t| (1.5|t| … 5|t|), каждый — честная граница.
Замкнутые формы остаются первыми: они у́же. Ни одного float: float на
входе — TypeError построением.
"""
from __future__ import annotations

from .rational import Q

from .enclosure import Interval, exp_enclosure
from .errors import Unresolved
from .zero_state import Zero

#: Допустимые скаляры. float исключён построением.
_NUM = (int, Q)

#: `[вычислено здесь]` 10.10, очередь строка 105, шаг десятый: дверь оператора
#: носителя = дверь ФIOS (как у восьми прежних носителей). Для ЧУЖОГО класса с
#: записью единый `_refuse` ведёт в `dispatch.combine` с порядком операндов
#: как в выражении; класс 1 (int, Q) и одноимённая линейка — прежний отказ
#: носителя. Ввоз ленивый, внутри функций.
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
    return (not isinstance(x, (_NUM, DESolutionValue))) and callable(getattr(x, "record", None))
#: Одна запись отказа на bool во всех дверях модуля: bool числом в доме не
#: считается и называется bool'ом (канон `phicore/rational.py`; сосед по
#: главе отвечает тем же: `Notch.from_q(True, 8)` → «bool не число»).
#: Двух копий одного канона здесь нет намеренно: разошедшиеся копии — уже
#: пройденная поломка (`affine.py:75`).
_NOT_BOOL = "bool не число"


def _rec(x: object) -> str:
    """Ф-запись операнда для сообщений отказа."""
    r = getattr(x, "record", None)
    if callable(r):
        return r()
    if isinstance(x, _NUM):
        return f"Φ(C=1 RATIONAL, value={Q(x)})"
    return f"Φ(?={x!r})"


def _as_fraction(x, name: str) -> Q:
    """int/Q → Q; float и прочее — TypeError построением."""
    if isinstance(x, bool):
        raise TypeError(f"{name} — {_NOT_BOOL}; дверь принимает int/Q")
    if not isinstance(x, _NUM):
        raise TypeError(f"{name} — только int/Q, float запрещён")
    return Q(x)


def _as_poly(coeffs, name: str) -> "tuple[Q, ...]":
    """Список коэффициентов (млад→ст) → кортеж Q.

    Канон — без хвостовых нулей: старшие нулевые степени записи не несут
    (P = [1, 0] и P = [1] — одно уравнение). float — TypeError."""
    if not isinstance(coeffs, (list, tuple)):
        raise TypeError(f"{name} — список коэффициентов int/Q")
    out = []
    for c in coeffs:
        if isinstance(c, bool):
            raise TypeError(f"{name}: коэффициенты — {_NOT_BOOL}; "
                            "дверь принимает int/Q")
        if not isinstance(c, _NUM):
            raise TypeError(f"{name}: коэффициенты — только int/Q, "
                            "float запрещён")
        out.append(Q(c))
    while out and out[-1] == 0:
        out.pop()
    return tuple(out)


class DESolutionValue:
    """Решение задачи Коши y′ = P(x)·y + Q(x), y(0) = y0 — носитель.

    p, q — кортежи Q (многочлены млад→ст, канон без хвостовых
    нулей), y0 — Q."""

    __slots__ = ("p", "q", "y0")

    def __init__(self, p, q, y0):
        self.p = _as_poly(p, "P")
        self.q = _as_poly(q, "Q")
        self.y0 = _as_fraction(y0, "y0")

    # --- ряд Тейлора: точная рекуррента --------------------------------

    def coeffs(self, n: int) -> "list[Q]":
        """Первые n коэффициентов Тейлора a_0 … a_{n−1}.

        Точная рекуррента (m+1)·a_{m+1} = Σ_i P_i·a_{m−i} + Q_m:
        свёртка по уже найденным a, деление на (m+1) — точное
        (Q). Таблиц значений нет."""
        if isinstance(n, bool):
            raise TypeError(f"n — {_NOT_BOOL}; дверь принимает int")
        if not isinstance(n, int):
            raise TypeError("n — только int")
        if n < 0:
            raise ValueError("n ≥ 0")
        if n == 0:
            return []
        a = [self.y0]
        for m in range(n - 1):
            s = self.q[m] if m < len(self.q) else Q(0)
            for i, pi in enumerate(self.p):
                if i > m:
                    break
                s += pi * a[m - i]
            a.append(s / (m + 1))
        return a

    # --- двойственность: проверка подстановкой -------------------------

    def check_substitution(self, n: int) -> bool:
        """Подстановка ряда в уравнение — структурная точная сверка.

        Правая часть P·y + Q собирается НЕЗАВИСИМО от рекурренты —
        полным произведением многочлена P на отрезок ряда a
        (накопление по парам (i, j)), левая — (m+1)·a_{m+1};
        равенство обязано сойтись до последнего бита при всех m < n."""
        if isinstance(n, bool):
            raise TypeError(f"n — {_NOT_BOOL}; дверь принимает int")
        if not isinstance(n, int):
            raise TypeError("n — только int")
        if n < 0:
            raise ValueError("n ≥ 0")
        a = self.coeffs(n + 1)  # a_0 … a_n: слева нужен a_{m+1}, m < n
        rhs = [self.q[m] if m < len(self.q) else Q(0)
               for m in range(n)]
        for i, pi in enumerate(self.p):
            if pi == 0:
                continue
            for j, aj in enumerate(a):
                if i + j >= n:
                    break
                rhs[i + j] += pi * aj
        return all((m + 1) * a[m + 1] == rhs[m] for m in range(n))

    # --- арифметика: направленного правила нет, кроме ноля-состояния ---

    def _refuse(self, op: str, other: object, self_left: bool):
        """Чужой класс с записью — дверь ФIOS (порядок операндов как в
        выражении); класс 1, одноимённая линейка и прибитые ячейки —
        прежний отказ носителя."""
        if _foreign_record(other) and not _pinned(other, {"×": "*"}.get(op, op)):
            a, b = (self, other) if self_left else (other, self)
            return _door(a, b, {"×": "*"}.get(op, op))
        left, right = ((self.record(), _rec(other)) if self_left
                       else (_rec(other), self.record()))
        raise Unresolved(op, left, right,
                         hint="решение ОДУ независимо: направленного "
                              "правила арифметики нет")

    def __add__(self, other):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__radd__
        return self._refuse("+", other, True)

    def __radd__(self, other):
        if isinstance(other, Zero):
            return NotImplemented
        return self._refuse("+", other, False)

    def __sub__(self, other):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__rsub__
        return self._refuse("-", other, True)

    def __rsub__(self, other):
        if isinstance(other, Zero):
            return NotImplemented
        return self._refuse("-", other, False)

    def __mul__(self, other):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__rmul__
        return self._refuse("×", other, True)

    def __rmul__(self, other):
        if isinstance(other, Zero):
            return NotImplemented
        return self._refuse("×", other, False)

    def __truediv__(self, other):
        return self._refuse("/", other, True)

    def __rtruediv__(self, other):
        return self._refuse("/", other, False)

    def __pow__(self, n: int):
        """Целая степень: n = 0 — точная единица, n = 1 — тот же носитель.

        КАНОН НУЛЕВОЙ СТЕПЕНИ (один на дом, замерен 12.09 по строке 75
        очереди; прежде на один вопрос стояло два ответа — `IV ** 0`
        давало 1, здесь был отказ). Нулевая степень есть ПУСТОЕ
        ПРОИЗВЕДЕНИЕ: факт о строении операции, не о величине основания.
        Потому она равна 1 у ВСЯКОГО носителя, который несёт величину, —
        включая носитель, чья величина есть точный ноль:

          `Q(0) ** 0` → 1 (оговорено прямо, `rational.py:196—197`:
            «0^Q(0) остаётся отказом, тогда как 0**0 дало бы 1»),
          `PhiNumberV6.from_int(0) ** 0` → 1 (класс 1),
          `Hyper.const(0) ** 0` → 1·1 (класс 7, пустое произведение
            выписано телом цикла: `hyper.py:420`),
          `IntegralValue ** 0` → 1, `ArcTrigValue ** 0` → 1.

        Отказ канон держит НЕ за нулевой величиной, а за носителем,
        который величины не несёт вовсе: класс 0 (пустота) и класс 20
        (не вычислено) — там степень возвращает сам носитель
        (`container.py:1028`: «пустота не становится единицей даже при
        n=0»), а не единицу и не отказ. Класс 18 величину НЕСЁТ: задача
        Коши линейного уравнения определяет решение однозначно, и
        y⁰ есть пустое произведение этой величины. Прежний довод «даже
        y⁰ не 1 там, где решение несёт ноль» есть довод о ПРЕДЕЛЕ 0^0 —
        дом считает точно, а не пределом, и для класса 1 уже ответил 1.

        Степень выше первой не сведена: произведения решений линейка не
        несёт (правила арифметики класса 18 нет).
        """
        if not isinstance(n, int) or isinstance(n, bool):
            return NotImplemented
        if n == 0:
            return Q(1)
        if n == 1:
            return DESolutionValue(self.p, self.q, self.y0)
        raise Unresolved("^", self.record(), _rec(n),
                         hint="степень решения вне {0, 1} не сведена: "
                              "произведения решений линейка не несёт")

    def __neg__(self) -> "DESolutionValue":
        """Линейность уравнения: −y решает y′ = P·y − Q, y(0) = −y0."""
        return DESolutionValue(self.p, tuple(-c for c in self.q), -self.y0)

    # --- запись, сравнение, сертификат ---------------------------------

    def record(self) -> str:
        ps = ",".join(str(c) for c in self.p)
        qs = ",".join(str(c) for c in self.q)
        return f"Φ(C=18 DE_SOLUTION, P=[{ps}], Q=[{qs}], y0={self.y0})"

    __repr__ = record

    def __eq__(self, other):
        """Тождество НОСИТЕЛЯ (записи), не решений: разные (P, Q, y0)
        могут нести одно решение (y′=y и y′=2y при y0=0 оба несут
        тождественный ноль) — сравнение величин здесь не сведено."""
        if isinstance(other, DESolutionValue):
            return (self.p == other.p and self.q == other.q
                    and self.y0 == other.y0)
        return NotImplemented

    def __hash__(self):
        return hash((DESolutionValue, self.p, self.q, self.y0))

    # --- самопроверка носителя (ось самоприменения, 08.10) ----------------

    def self_check(self, deep: bool = False) -> list:
        """Самопроверка носителя встречными дорогами (ось «AIAΣA спрашивает
        AIAΣA», слово автора 04.09; очередь строка 18; 08.10). Список
        находок, пусто — здоров; исключений не бросает.

        `[вычислено здесь]` (1) слоты: p, q — кортежи Q без хвостовых нулей
        (канон `_as_poly`), y0 — Q; (2) двойственность подстановкой до
        n = 24 (`check_substitution`: правая часть собирается произведением
        многочленов, не рекуррентой); (3) ВСТРЕЧНАЯ дорога сертификата в
        точке t = 1/2: там, где есть замкнутая форма (`enclosure` дверью
        exp), — оценка Коши (`cauchy_enclosure`: частная сумма ряда плюс
        мажоранта хвоста по Грёнуоллу) обязана её пересекать — две дороги
        независимы по построению; в общем случае — оценка Коши на 30 и на
        60 членах пересекается, и длинная не шире; (4) линейность: у −y
        коэффициенты ряда противоположны. При deep — ещё t ∈ {1, −1/2} и
        подстановка до n = 48."""
        out = []
        p, q, y0 = self.p, self.q, self.y0
        if not isinstance(p, tuple) or not isinstance(q, tuple) or not all(isinstance(c, Q) for c in p + q):
            return ["слоты P/Q чужого носителя: ожидались кортежи Q"]
        if not isinstance(y0, Q):
            return [f"слот y0 чужого носителя: {type(y0).__name__}"]
        if (p and p[-1] == 0) or (q and q[-1] == 0):
            out.append("слоты вне канона: хвостовой нуль в P или Q")
        try:
            n = 48 if deep else 24
            if not self.check_substitution(n):
                out.append(f"подстановка ряда в уравнение разошлась до n = {n}")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"подстановка не проверена: {type(e).__name__}: {e}")
        try:
            neg = -self
            if neg.coeffs(12) != [-c for c in self.coeffs(12)]:
                out.append("линейность: коэффициенты −y не противоположны коэффициентам y")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"линейность не проверена: {type(e).__name__}: {e}")
        closed = len(q) == 0 or len(p) == 0 or (len(p) == 1 and len(q) == 1)
        for t in ((Q(1, 2), Q(1), Q(-1, 2)) if deep else (Q(1, 2),)):
            try:
                c1 = self.cauchy_enclosure(t, 30)
                if closed:
                    iv = self.enclosure(t, 30)
                    if iv.intersect(c1) is None:
                        out.append(f"замкнутая форма и оценка Коши разошлись при t = {t}: {iv.record()} против {c1.record()}")
                else:
                    c2 = self.cauchy_enclosure(t, 60)
                    if c1.intersect(c2) is None:
                        out.append(f"оценки Коши на 30 и 60 членах не пересекаются при t = {t}: {c1.record()} против {c2.record()}")
                    if c2.width() > c1.width():
                        out.append(f"оценка Коши на 60 членах шире, чем на 30, при t = {t}")
            except Unresolved as e:
                out.append(f"сертификат при t = {t} отказал: {e}")
            except Exception as e:                              # noqa: BLE001
                out.append(f"сертификат при t = {t}: {type(e).__name__}: {e}")
        return out

    def enclosure(self, t, terms: int = 30) -> Interval:
        """Сертификат сетки [lo, hi] ⊂ ℚ значения y(t), по явному
        запросу.

        Три замкнутых случая (01.10, очередь строка 104; прежде — только
        P = [c], Q = []): Q = [] при любом P — y0·e^{F(t)}, F = ∫₀ᵗ P точно
        в ℚ; P = [] при любом Q — точка y0 + ∫₀ᵗ Q; P = [c], Q = [q] —
        (y0 + q/c)·e^{c·t} − q/c. Остальное — `cauchy_enclosure`: частная
        сумма ряда и мажоранта хвоста оценкой Коши (01.10, вечер)."""
        t = _as_fraction(t, "t")
        if isinstance(terms, bool):
            raise TypeError(f"terms — {_NOT_BOOL}; дверь принимает int")
        if not isinstance(terms, int):
            raise TypeError("terms — только int")
        if terms < 1:
            raise ValueError("terms ≥ 1")
        p, q = self.p, self.q
        if len(q) == 0:
            # y = y0·e^{F(t)}, F(t) = Σ P_i·t^{i+1}/(i+1) — точно в ℚ
            F = sum((pi * t ** (i + 1) / (i + 1) for i, pi in enumerate(p)), Q(0))
            return exp_enclosure(F, terms) * self.y0
        if len(p) == 0:
            # y = y0 + ∫₀ᵗ Q — многочлен, значение точно: точка класса 1
            G = sum((qi * t ** (i + 1) / (i + 1) for i, qi in enumerate(q)), Q(0))
            return Interval.point(self.y0 + G)
        if len(p) == 1 and len(q) == 1:
            c, q0 = p[0], q[0]          # c ≠ 0: хвостовые нули канонизированы
            return exp_enclosure(c * t, terms) * (self.y0 + q0 / c) - q0 / c
        return self.cauchy_enclosure(t, terms)

    def cauchy_enclosure(self, t, terms: int = 30) -> Interval:
        """Сертификат общего случая: частная сумма + мажоранта хвоста
        оценкой Коши; работает при любых P, Q (и там, где есть замкнутая
        форма, — но шире её).

        `[вычислено здесь]` Решение y′ = P·y + Q с многочленными P, Q —
        целая функция. На луче z = r·e^{iθ}, 0 ≤ r ≤ ρ:
        |y|′ ≤ |y′| ≤ |P(z)|·|y| + |Q(z)| ≤ P̃(ρ)·|y| + Q̃(ρ), где
        P̃(ρ) = Σ|P_i|·ρ^i; интегральная форма Грёнуолла даёт
        |y(z)| ≤ (|y0| + ρ·Q̃(ρ))·e^{ρ·P̃(ρ)} =: M(ρ) на всём круге.
        Неравенство Коши для коэффициентов: |a_n| ≤ M(ρ)/ρⁿ; при |t| < ρ
        хвост |Σ_{n≥N} a_n·tⁿ| ≤ M(ρ)·r^N/(1 − r), r = |t|/ρ.
        e^{ρ·P̃(ρ)} сверху — верх `exp_enclosure` (сам с мажорантой).
        Радиус: восемь кандидатов ρ = |t|·(1 + k/2), k = 1..8, берётся
        наименьший хвост — каждый кандидат есть честная граница, минимум
        честных границ честен. При t = 0 — точка y0."""
        t = _as_fraction(t, "t")
        if isinstance(terms, bool):
            raise TypeError(f"terms — {_NOT_BOOL}; дверь принимает int")
        if not isinstance(terms, int):
            raise TypeError("terms — только int")
        if terms < 1:
            raise ValueError("terms ≥ 1")
        if t == 0:
            return Interval.point(self.y0)
        a = self.coeffs(terms)
        partial = sum((an * t ** n for n, an in enumerate(a)), Q(0))
        tail = None
        for k in range(1, 9):
            rho = abs(t) * (1 + Q(k, 2))
            p_abs = sum((abs(c) * rho ** i for i, c in enumerate(self.p)), Q(0))
            q_abs = sum((abs(c) * rho ** i for i, c in enumerate(self.q)), Q(0))
            bound = (abs(self.y0) + rho * q_abs) * exp_enclosure(rho * p_abs, terms).hi
            r = abs(t) / rho
            cand = bound * r ** terms / (1 - r)
            if tail is None or cand < tail:
                tail = cand
        return Interval(partial - tail, partial + tail)


# ---------- самопроверка линейки (ось самоприменения, 08.10) ----------

def self_check(deep: bool = False) -> list:
    """Самопроверка линейки класса 18 встречными дорогами; список находок,
    пусто — здорова; исключений не бросает.

    `[вычислено здесь]` Окно носителей: y′ = y (e^t), y′ = 2x·y (e^{x²}),
    y′ = y + 1 при y0 = 0 (e^t − 1), общие y′ = x·y + 1 и y′ = (1 + x)·y + x
    (deep — ещё y′ = (x² − 1)·y + 2 при y0 = 1/3) проверяют себя
    (`DESolutionValue.self_check`); внешние факты: коэффициенты e^t суть
    1/n!, у e^{x²} — 1/k! на чётных местах и нули на нечётных, у e^t − 1 —
    те же 1/n! при n ≥ 1 и нуль в начале [внешний источник: ряд
    экспоненты]; оценка Коши для e^t в t = 1 обязана нести e дверью
    `exp_enclosure` — дорога ряда против двери. Граница — окно."""
    out = []
    window = [([1], [], 1), ([0, 2], [], 1), ([1], [1], 0), ([0, 1], [1], 1), ([1, 1], [0, 1], 1)]
    if deep:
        window.append(([-1, 0, 1], [2], Q(1, 3)))
    for p, q, y0 in window:
        try:
            v = DESolutionValue(p, q, y0)
            out += [f"{v.record()}: {f}" for f in v.self_check(deep)]
        except Exception as e:                                  # noqa: BLE001
            out.append(f"носитель P={p}, Q={q}, y0={y0}: {type(e).__name__}: {e}")
    try:
        fact = [Q(1)]
        for n in range(1, 12):
            fact.append(fact[-1] / n)
        if DESolutionValue([1], [], 1).coeffs(12) != fact:
            out.append("y′ = y: коэффициенты не 1/n!")
        sq = [fact[k // 2] if k % 2 == 0 else Q(0) for k in range(12)]
        if DESolutionValue([0, 2], [], 1).coeffs(12) != sq:
            out.append("y′ = 2x·y: коэффициенты не 1/k! на чётных местах")
        if DESolutionValue([1], [1], 0).coeffs(12) != [Q(0)] + fact[1:]:
            out.append("y′ = y + 1, y0 = 0: коэффициенты не 1/n! при n ≥ 1")
        e_c = DESolutionValue([1], [], 1).cauchy_enclosure(Q(1), 40)
        if e_c.intersect(exp_enclosure(Q(1), 40)) is None:
            out.append(f"оценка Коши для e^t в t = 1 не несёт e: {e_c.record()}")
    except Exception as e:                                      # noqa: BLE001
        out.append(f"внешние факты: {type(e).__name__}: {e}")
    return out
