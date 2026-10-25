"""Носитель класса 7 POWER_SERIES: x = Σ aₙ·rⁿ с рекуррентой отношения.

Линейка держит ряд СТРОЕНИЕМ, а не готовой суммой: a0 ∈ ℚ, отношение
соседних коэффициентов aₙ₊₁/aₙ = ratio_num(n)/ratio_den(n) —
рациональная функция номера (ratio_num, ratio_den ∈ ℤ[n], intpoly.IntPoly),
точка r ∈ ℚ. Каждый член и каждая частная сумма точны: Q,
рекуррента с кэшем, ни одного float — float на входе даёт TypeError
построением.

Горячий цикл рекурренты — на парах целых с отложенной нормализацией
(этап F выселения Q): множитель ratio(k)·r — целая пара, числитель
члена и числитель частной суммы растут целыми произведениями на общем
знаменателе, gcd берётся один раз при выдаче Q наружу — форма
контейнера v6 (сокращение по запросу, не на каждом действии), результат
бит-в-бит тот же. Страховка от разбухания: при росте пары за _PAIR_BITS
бит — разовое сокращение общим gcd, и счёт продолжается (нормализация
отложенная, не отменённая).

Десятичная сторона — только сертификатом сетки (enclosure), по явному
запросу: частная сумма до n* ≥ N плюс строгий хвост геометрической
мажорантой. Отношение членов |term(n+1)/term(n)| = |ratio(n)·r|
накрывается монотонно невозрастающей рациональной гранью, и хвост
сводится к |term(n*+1)|/(1−ρ) при ρ < 1; грань монотонна лишь при
deg(ratio_num) ≤ deg(ratio_den) — иначе честный отказ Unresolved:
расходимость мажоранты называется, не прячется за приближением.

Двойственность с классом 5 EXPONENTIAL: ряд aₙ = 1/n! (a0=1,
ratio 1/(n+1)) в точке r=1 несёт ту же величину, что e¹; сертификаты
двух линеек обязаны пересекаться — согласие проверяется приёмкой.

Арифметика между рядами и с числом — с 05.09 замыканием D-конечных
рядов (dfinite.py, Стэнли 1980): сумма, разность, свёртка Коши и целая
степень ≥ 2 дают ряд класса 7 с точной рекуррентой (одна линейка, не
структурная запись); одна рекуррента и один r по-прежнему складываются
по a0 (точный ноль при взаимном уничтожении). Деление ряда на ряд и
обращение — отказ: общего замыкания нет (теорема). +/− с Zero —
нейтральность нуля-состояния, смена знака — по a0.
"""
from __future__ import annotations

from .rational import Q
from math import gcd

from .enclosure import Interval
from .errors import ExactRepresentationNotAvailable, Unresolved
from .intpoly import IntPoly
from .zero_state import Zero

#: Допустимые скаляры. float исключён построением.
_NUM = (int, Q)

#: Порог страховки: пара целых длиннее этого числа бит сокращается
#: разовым gcd прямо в цикле — отложенная нормализация, не отменённая.
_PAIR_BITS = 1 << 20

#: Предел поиска n* с ρ < 1 сверх запрошенного N.
_SEARCH_LIMIT = 4096
#: Ширина рационального замера отношения на диапазоне за n*.
_CHECK_RANGE = 64


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
    if isinstance(x, bool) or not isinstance(x, _NUM):
        raise TypeError(f"{name} — только int/Q, float запрещён")
    return Q(x)


def hole_hint(holes, at=None) -> str:
    """ОДИН текст о дыре ratio_den на оба слоя — линейку и замыкание.

    Прежде линейка называла дыру при встрече («ratio_den(3)=0»), а замыкание
    — сразу и полным списком («n ∈ [3]»): один факт, два разных диагноза
    (глава 11 §20 п. 5, замер 12.09). Текст теперь собирается здесь, и обе
    стороны берут его отсюда.

    `holes` — ПОЛНЫЙ список дыр, когда целые корни `ratio_den` отделены;
    `None` — когда отделение отказало, и тогда называется дыра `at`, на
    которой стоим, а неполнота объявляется числом. Пустой список есть
    утверждение «дыр нет», и говорить его о неотделённом многочлене нельзя:
    замер 17.09 `[вычислено здесь]` — при `ratio_den = (n−3)(n+p·q)`,
    p = 1 000 003, q = 1 000 033, отделение отказывает (свободный член не
    разложен пробным делением до 1 000 000), и прежний текст печатал
    «при n ∈ []», стоя на дыре n = 3. Ноль инструмента не есть ноль предмета.
    """
    if holes is None:
        return (f"ratio_den обращается в ноль при n = {at} — рекуррента "
                f"коэффициентов не определена, величина ряда там обрывается; "
                f"полный список дыр НЕ ОТДЕЛЁН (целые корни ratio_den не "
                f"сведены), названа та дыра, на которой счёт остановился")
    return (f"ratio_den обращается в ноль при n ∈ {sorted(holes)} — рекуррента "
            f"коэффициентов не определена, величина ряда там обрывается")


def _as_index(n, name: str) -> int:
    """Номер члена: целое ≥ 0; bool и float исключены построением."""
    if isinstance(n, bool) or not isinstance(n, int):
        raise TypeError(f"{name} — только int")
    if n < 0:
        raise ValueError(f"{name} должен быть ≥ 0, получено {n}")
    return n


class PowerSeriesValue:
    """Σ aₙ·rⁿ: a0 ≠ 0 (Q), ratio_num/ratio_den ∈ ℤ[n], r ∈ ℚ.

    term(n) = aₙ·rⁿ считается рекуррентой term(n+1) = term(n)·ratio(n)·r
    и кэшируется; канона записи у линейки нет — разные строения могут
    нести одну сумму, поэтому тождество решается только совпадением
    строения, иное сравнение — отказ."""

    __slots__ = ("a0", "ratio_num", "ratio_den", "r", "_terms")

    def __init__(self, a0, ratio_num, ratio_den, r):
        a0 = _as_fraction(a0, "a0")
        if a0 == 0:
            raise ValueError("a0=0 гасит рекуррентой весь ряд — ноль "
                             "живёт в zero_state")
        if not isinstance(ratio_num, IntPoly) \
                or not isinstance(ratio_den, IntPoly):
            raise TypeError("ratio_num и ratio_den — только intpoly.IntPoly")
        if ratio_den.is_zero():
            raise ValueError("ratio_den — нулевой многочлен: отношение "
                             "aₙ₊₁/aₙ не определено")
        self.a0 = a0
        self.ratio_num, self.ratio_den = ratio_num, ratio_den
        self.r = _as_fraction(r, "r")
        #: кэш троек целых: _terms[n] = (числитель aₙ·rⁿ, числитель
        #: Σ_{k≤n}, общий знаменатель > 0); term(0) = a0·r⁰ = a0.
        self._terms = [(a0.numerator, a0.numerator, a0.denominator)]

    # --- члены и частные суммы: точная рекуррента на парах целых --------

    def _extend(self, n: int) -> None:
        """Дотягивает кэш троек до номера n целой рекуррентой.

        Шаг k → k+1: множитель ratio(k)·r — целая пара (m_num, m_den),
        m_den > 0 (знак уходит в числитель); числитель члена умножается
        на m_num, общий знаменатель — на m_den, числитель частной суммы
        сносится на новый знаменатель и принимает член. gcd на шаге не
        берётся — нормализация отложена до выдачи Q наружу;
        страховка: тройка длиннее _PAIR_BITS бит сокращается разовым
        общим gcd, и счёт продолжается."""
        rp, rq = self.r.numerator, self.r.denominator
        while len(self._terms) <= n:
            k = len(self._terms) - 1
            pd = self.ratio_den.eval(k)
            if pd == 0:
                try:
                    holes = self.holes()
                except Unresolved:
                    holes = None          # список не отделён — назовём дыру под ногами
                raise Unresolved("term", self.record(), hint=hole_hint(holes, at=k))
            pn = self.ratio_num.eval(k)
            m_num, m_den = pn * rp, pd * rq
            if m_den < 0:
                m_num, m_den = -m_num, -m_den
            t, s, d = self._terms[-1]
            t *= m_num
            s = s * m_den + t
            d *= m_den
            if max(t.bit_length(), s.bit_length(),
                   d.bit_length()) > _PAIR_BITS:
                g = gcd(t, s, d)
                if g > 1:
                    t, s, d = t // g, s // g, d // g
            self._terms.append((t, s, d))

    def holes(self) -> list:
        """Неотрицательные целые n, где ratio_den(n) = 0 [вычислено здесь].

        Там рекуррента коэффициентов не определена, и величина ряда обрывается:
        линейка считает члены ДО первой дыры и отказывает на ней, а подъём в
        замыкание (`dfinite.from_power`) отказывает СРАЗУ — ему нужна
        рекуррента, определённая при всех n. Асимметрия объявлена нарочно:
        частная сумма до дыры существует, а суммы ряда не существует.

        Отделение корней может отказать само (предел пробного деления
        `_TRIAL`) — тогда отсюда уходит `Unresolved`, А НЕ ПУСТОЙ СПИСОК:
        пустой список есть утверждение «дыр нет», и прежний `return []`
        превращал отказ измерителя в ложное утверждение о предмете (замер
        17.09 — текст «при n ∈ []», сказанный стоя на дыре n = 3).
        """
        from .dfinite import _int_roots
        return [k for k in _int_roots(self.ratio_den) if k >= 0]

    def term(self, n: int) -> Q:
        """Точный член aₙ·rⁿ: term(n+1) = term(n)·ratio(n)·r, кэш.

        Внутри — пара целых на общем знаменателе; gcd один раз здесь,
        при выдаче Q наружу (бит-в-бит то же значение)."""
        n = _as_index(n, "n")
        self._extend(n)
        t, _, d = self._terms[n]
        return Q(t, d)

    def partial_sum(self, N: int) -> Q:
        """Точная частная сумма Σ_{n=0}^{N} aₙ·rⁿ (Q).

        Числитель суммы копится в тройке кэша на общем знаменателе
        члена; gcd один раз здесь, при выдаче Q наружу."""
        N = _as_index(N, "N")
        self._extend(N)
        _, s, d = self._terms[N]
        return Q(s, d)

    # --- сертификат сетки -----------------------------------------------

    def _dfinite_road(self, N: int, why: str) -> Interval:
        """Вторая дорога сертификата (02.10): тот же ряд как ГОЛЫЙ D-конечный
        (`DFiniteSeries.from_power` — рекуррента порядка 1 точно в ℤ[n]), и
        мажоранта хвоста доказывается сдвигом Тейлора с подбором (n*, ρ) по
        сетке (`dfinite._auto_tail`). Берётся, когда грубая мажоранта
        |r|·A/B не опускается ниже 1 либо ρ < 1 не найдена: у −ln(1−x) при
        x = 1/2 грубая даёт предел 1 (отказ до 02.10), точная — ρ = 3/4 с n* = 2.
        Отказ второй дороги склеивается с причиной первой — по имени."""
        from .dfinite import DFiniteSeries
        try:
            raw = DFiniteSeries(DFiniteSeries.from_power(self).seq)
            return raw.enclosure(N)
        except Unresolved as e:
            raise Unresolved("enclosure", self.record(),
                             hint=why + "; D-конечная дорога (точная рекуррента, "
                                  "знаки сдвигом Тейлора, сетка (n*, ρ)) тоже не "
                                  "доказала: хвост геометрией не сведён — "
                                  + str(e).split("; ")[-1][:120]) from None

    def enclosure(self, N: int = 40) -> Interval:
        """Сертификат сетки [lo, hi] ⊂ ℚ: частная сумма + строгий хвост.

        Мажоранта отношения: при n ≥ 1 верно |ratio_num(n)| ≤ A·n^dp и
        |ratio_den(n)| ≥ B·n^dq − C·n^{dq−1}, где A — сумма модулей
        коэффициентов числителя, B — модуль старшего коэффициента
        знаменателя, C — сумма модулей остальных; значит при B·n > C

            |ratio(n)·r| ≤ M(n) = |r|·A·n^{dp+1} / (n^{dq}·(B·n − C)),

        и M(n) монотонно не возрастает по n при dp ≤ dq (проверка
        отношением M(n+1)/M(n) ≤ 1). Берётся наименьшее n* ≥ N с
        ρ = M(n*+1) < 1; тогда для всех n ≥ n*+1 отношение ≤ ρ, ряд
        мажорируется геометрией и хвост ≤ |term(n*+1)|/(1−ρ). Грань
        дополнительно замеряется рационально на диапазоне за n*.

        При dp > dq мажоранта растёт — честный отказ Unresolved
        (расходимость мажоранты), не приближение."""
        N = _as_index(N, "N")
        p, q = self.ratio_num, self.ratio_den
        dp, dq = p.degree, q.degree
        if not p.is_zero() and dp > dq:
            raise Unresolved("enclosure", self.record(),
                             hint="deg(ratio_num) > deg(ratio_den): "
                                  "мажоранта отношения |ratio(n)·r| растёт "
                                  "с n — расходимость мажоранты, хвост "
                                  "геометрией не сведён")
        A = sum(abs(c) for c in p.c)
        B = abs(q.c[-1])
        C = sum(abs(c) for c in q.c[:-1])
        r_abs = abs(self.r)
        if dp == dq and r_abs * A >= B:
            # 02.10: грубая мажоранта (сумма модулей против старшего) —
            # предел ИНСТРУМЕНТА (замер 12.09, сетка двойственности): у
            # −ln(1−x) истинный радиус 1, а она даёт 1/2. Второй дорогой идёт
            # точная рекуррента класса 7 (D-конечная): знаки коэффициентов и
            # неравенство ρ·|p₁(n)| − |p₀(n)| ≥ 0 сдвигом Тейлора, (n*, ρ) по
            # сетке (`dfinite._auto_tail`); расходящийся ряд там честно
            # отказывает (ρ < 1 недоказуема).
            return self._dfinite_road(N, f"мажоранта |ratio(n)·r| не опускается "
                                         f"ниже 1 (её предел {r_abs * Q(A, B)} ≥ 1)")

        def major(m: int) -> Q:
            # Верхняя грань |ratio(n)·r| для всех n ≥ m ≥ 1 при B·m > C.
            return r_abs * A * Q(m ** (dp + 1), m ** dq * (B * m - C))

        n_star = N
        while True:
            m = n_star + 1
            if B * m > C:
                rho = major(m)
                if rho < 1:
                    break
            n_star += 1
            if n_star - N > _SEARCH_LIMIT:
                return self._dfinite_road(N, f"ρ < 1 не найдено при n* ∈ "
                                             f"[{N}; {N + _SEARCH_LIMIT}]")
        # Рациональный замер на диапазоне: отношение не превышает ρ.
        for n in range(m, m + _CHECK_RANGE):
            dv = q.eval(Q(n))
            if dv == 0 or abs(self.r * p.eval(Q(n)) / dv) > rho:
                raise Unresolved("enclosure", self.record(),
                                 hint=f"замер отношения на диапазоне "
                                      f"опроверг мажоранту при n={n}")
        s = self.partial_sum(n_star)
        tail = abs(self.term(n_star + 1)) / (1 - rho)
        return Interval(s - tail, s + tail)

    # --- арифметика: правил между рядами нет — честный отказ ------------

    def self_check(self, deep: bool = False) -> list:
        """Самопроверка ряда ВСТРЕЧНОЙ дорогой (ось самоприменения).

        `[вычислено здесь]` Члены и частные суммы считаются ОДНОЙ целой
        рекуррентой на общем знаменателе (кэш троек), поэтому встречные
        дороги идут мимо неё: (1) частная сумма собирается здесь заново —
        прямым сложением `term(i)` парами `Q`, без кэша сумм, и обязана
        совпасть с `partial_sum(N)`; (2) отношение соседних членов
        пересчитывается прямо из многочленов рекурренты:
        term(n+1)/term(n) обязано равняться ratio_num(n)/ratio_den(n)·r —
        то есть рекуррента сверяется со СВОИМ объявлением, а не сама с
        собой; (3) сертификат обязан нести свою частную сумму, а уточнение
        (больше членов) — вкладываться в грубый.

        Граница объявлена: сходимость здесь НЕ доказывается — мажоранта
        хвоста живёт в `enclosure` и честно отказывает, когда ряд ею не
        мажорируется; такой отказ записывается находкой только если он
        пришёл не из объявленной двери. Нулевой знаменатель рекурренты в
        проверяемом окне — находка, а вне окна (n > N) молчание не есть
        утверждение. Возвращает список находок; исключений не бросает.
        """
        out = []
        if not isinstance(self.a0, Q) or self.a0 == 0:
            out.append(f"a0 не ненулевая пара: {self.a0!r}")
            return out
        if not isinstance(self.r, Q):
            out.append(f"r не пара: {self.r!r}")
            return out
        N = 12 if not deep else 24
        try:
            acc = Q(0)
            for i in range(N + 1):
                acc = acc + self.term(i)
            if acc != self.partial_sum(N):
                out.append(f"частная сумма двумя дорогами расходится: прямым "
                           f"сложением {acc}, рекуррентой {self.partial_sum(N)}")
        except Unresolved as e:
            out.append(f"члены не сведены в окне n ≤ {N}: {e}")
            return out
        except Exception as e:                                  # noqa: BLE001
            out.append(f"счёт членов отказал: {type(e).__name__}: {e}")
            return out
        for k in range(N):
            t_k, t_k1 = self.term(k), self.term(k + 1)
            pd = self.ratio_den.eval(k)
            if pd == 0:
                out.append(f"знаменатель рекурренты нулевой при n = {k}")
                break
            want = Q(self.ratio_num.eval(k), pd) * self.r
            if t_k == 0:
                continue                      # член погашен: отношение не определено
            if t_k1 != t_k * want:
                out.append(f"отношение членов при n = {k} не отвечает "
                           f"объявленной рекурренте: {t_k1} против {t_k * want}")
                break
        try:
            iv = self.enclosure(N)
        except (Unresolved, ExactRepresentationNotAvailable):
            return out                        # объявленная дверь мажоранты
        except Exception as e:                                  # noqa: BLE001
            out.append(f"сертификат отказал не по объявленной двери: "
                       f"{type(e).__name__}: {e}")
            return out
        if not iv.contains(self.partial_sum(N)):
            out.append(f"сертификат не несёт свою частную сумму: "
                       f"{iv.record()} против {self.partial_sum(N)}")
        if not deep:
            return out
        try:
            fine = self.enclosure(N + 8)
        except Exception:                                       # noqa: BLE001
            return out
        if not iv.encloses(fine):
            out.append(f"уточнение не вложено в грубый сертификат: "
                       f"{fine.record()} ⊄ {iv.record()}")
        return out

    def _refuse(self, op: str, other, hint: str):
        if isinstance(other, _NUM) or callable(getattr(other, "record", None)):
            raise Unresolved(op, self.record(), _rec(other), hint=hint)
        return NotImplemented

    def _same_shape(self, other) -> bool:
        return (isinstance(other, PowerSeriesValue)
                and other.ratio_num == self.ratio_num
                and other.ratio_den == self.ratio_den and other.r == self.r)

    def _addsub(self, other, op: str, sign: int):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__radd__/__rsub__
        if self._same_shape(other):
            # одна рекуррента и один r: aₙ = a0·Π ratio, потому
            # Σ(aₙ ± bₙ)rⁿ = (a0 ± b0)·Σ(Π ratio)rⁿ — тот же ряд по a0;
            # взаимное уничтожение — точный ноль
            a0 = self.a0 + sign * other.a0
            if a0 == 0:
                return Q(0)
            return PowerSeriesValue(a0, self.ratio_num, self.ratio_den, self.r)
        from .dfinite import DFiniteSeries
        if not isinstance(other, bool) and DFiniteSeries._lift(other) is not None:
            # 05.09: иная пара, число, интеграл — замыкание D-конечных рядов
            # (dfinite.py): сумма одной линейкой с точной рекуррентой, сверено прогоном
            a, b = DFiniteSeries.from_power(self), DFiniteSeries._lift(other)
            return a + b if sign > 0 else a - b
        return self._refuse(op, other, "у операнда нет D-конечного ряда: правила нет")

    def __add__(self, other):
        return self._addsub(other, "+", +1)

    def __sub__(self, other):
        return self._addsub(other, "-", -1)

    def __radd__(self, other):
        if isinstance(other, _NUM) and not isinstance(other, bool):
            from .dfinite import DFiniteSeries
            return DFiniteSeries.const(other) + DFiniteSeries.from_power(self)
        return NotImplemented

    def __rsub__(self, other):
        if isinstance(other, _NUM) and not isinstance(other, bool):
            from .dfinite import DFiniteSeries
            return DFiniteSeries.const(other) - DFiniteSeries.from_power(self)
        return NotImplemented

    def __mul__(self, other):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__rmul__
        if isinstance(other, _NUM):
            q = Q(other)
            if q == 0:
                return Q(0)
            # масштаб парой — в a0: рекуррента множительна
            return PowerSeriesValue(self.a0 * q, self.ratio_num,
                                    self.ratio_den, self.r)
        from .dfinite import DFiniteSeries
        if not isinstance(other, bool) and DFiniteSeries._lift(other) is not None:
            # 05.09: свёртка Коши — замыкание D-конечных рядов (Стэнли 1980); интеграл — тем же мостом
            return DFiniteSeries.from_power(self) * DFiniteSeries._lift(other)
        return self._refuse("×", other, "произведение с этим операндом не сведено")

    __rmul__ = __mul__

    def __truediv__(self, other):
        if isinstance(other, _NUM):
            q = Q(other)
            if q == 0:
                raise Unresolved("/", self.record(), _rec(other),
                                 hint="x/0 здесь не сведено")
            return self * (1 / q)
        return self._refuse("/", other, "деление рядов не сведено")

    def __rtruediv__(self, other):
        if isinstance(other, _NUM):
            raise Unresolved("/", _rec(other), self.record(),
                             hint="обращение ряда не сведено")
        return NotImplemented

    def __pow__(self, n):
        if isinstance(n, bool) or not isinstance(n, int):
            return NotImplemented
        if n == 1:
            return self
        if n >= 2:
            from .dfinite import DFiniteSeries
            return DFiniteSeries.from_power(self) ** n
        # n=0 не сводится: сумма ряда может быть нулевой, x⁰ не решается
        # строением записи; отрицательная — обращение, замыкания нет.
        raise Unresolved("^", self.record(), _rec(n),
                         hint="степень ряда сведена для целых ≥ 1; 0 и отрицательные — нет")

    def __neg__(self) -> "PowerSeriesValue":
        """−Σ aₙ·rⁿ = Σ (−aₙ)·rⁿ: рекуррента множительна, знак несёт a0."""
        return PowerSeriesValue(-self.a0, self.ratio_num, self.ratio_den,
                                self.r)

    # --- запись и тождество ---------------------------------------------

    def record(self) -> str:
        num = ",".join(map(str, self.ratio_num.c))
        den = ",".join(map(str, self.ratio_den.c))
        return (f"Φ(C=7 POWER_SERIES, a0={self.a0}, "
                f"ratio=[{num}]/[{den}], r={self.r})")

    __repr__ = record

    def __eq__(self, other):
        if isinstance(other, PowerSeriesValue):
            if (self.a0 == other.a0 and self.ratio_num == other.ratio_num
                    and self.ratio_den == other.ratio_den
                    and self.r == other.r):
                return True
            # Канона записи нет: различие строений не решает различие
            # сумм — тождество величин не сведено, отказ вместо ложного
            # «не равно».
            raise Unresolved("=", self.record(), other.record(),
                             hint="тождество сумм разных строений не сведено")
        if isinstance(other, _NUM):
            raise Unresolved("=", self.record(), _rec(other),
                             hint="сумма ряда может быть рациональной — "
                                  "тождество с числом не решается строением")
        return NotImplemented

    def __hash__(self):
        return hash((PowerSeriesValue, self.a0, tuple(self.ratio_num.c),
                     tuple(self.ratio_den.c), self.r))
