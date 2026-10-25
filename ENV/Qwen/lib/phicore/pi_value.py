"""Носитель класса 10: мономы q·π^e — приземление двойственностей.

Сюда сводятся направленные правила других линеек: ζ(2n) (класс 13),
отражение Эйлера Γ(z)Γ(1−z) (класс 14), корни J_{1/2} (класс 16),
K(0)=E(0)=π/2 (класс 15), ∫dx/(1+x²) (класс 9). Величина хранится
структурно: коэффициент q (Q или AlgebraicNumber — sin πz
рационального аргумента алгебраичен) и показатель e ∈ ℚ (Γ(1/2)=π^{1/2}).

Правил сведения π с независимыми величинами нет: точные правила монома
(одноимённая степень, рациональный и алгебраический множитель) стоят
впереди, а для ЧУЖОГО класса оператор зовёт дверь ФIOS (`dispatch.combine`)
и возвращает её ответ — запись класса 19 там, где у таблицы правило mixed
(π·e, π + ln 2), и отказ двери насквозь там, где ячейка пуста (π · корень
единицы). Ложного класса нет ни на одной дороге. Ячейки, в которых приёмки
прибили отказ оператора поимённо (π + 1, π + √2, π · Γ, π ∘ корень Бесселя),
дверь не зовут — см. `_PINNED`. Десятичная сторона — только сертификатом
сетки (enclosure), по явному запросу.
"""
from __future__ import annotations

from .rational import Q

from .errors import Unresolved
from .algebraic import AlgebraicNumber
from .enclosure import Interval, pi_enclosure, root_enclosure
from .zero_state import Zero

_NUM = (int, Q)


def _is_zero_num(x) -> bool:
    """Нуль класса 1 — int 0 или Q(0), не bool: x ± 0 = x в любом классе.
    Это поглощение нейтрального элемента, не сведение разных классов;
    до 01.10 дверь отказывала «независим: правила нет» и на нуле (наблюдение
    сдачи 25.09, очередь строка 104; замер 01.10 `[вычислено здесь]`:
    π + 0, 0 + π, π − 0, 0 − π — четыре отказа из четырёх)."""
    return isinstance(x, _NUM) and not isinstance(x, bool) and Q(x) == 0


def _rec(x: object) -> str:
    r = getattr(x, "record", None)
    if callable(r):
        return r()
    if isinstance(x, _NUM):
        return f"Φ(C=1 RATIONAL, value={Q(x)})"
    return f"Φ(?={x!r})"


#: `[вычислено здесь]` Правило сложения ОДНО, потому и подсказка одна.
#: Прежде их было две — «сложение сведено только для одноимённой степени»
#: с одной стороны и «правила нет» с другой, — и вызывающий из `1 + π`
#: узнавал меньше, чем из `π + 1`, хотя отказ тот же (очередь п. 36).
_ADD_HINT = "π независим: сложение сведено только для одноимённой степени π^e"

#: `[вычислено здесь]` 09.10, очередь строка 105, шаг первый: дверь оператора
#: носителя = дверь ФIOS. Замер до правки (`proto/bench_operator_vs_route.py
#: --ten`, 400 ячеек): в 55 ячейках с классом 10 оператор отказывал, а
#: `Фios.route` той же пары давал запись класса 19; в 9 ячейках оператор падал
#: TypeError интерпретатора (NotImplemented с обеих сторон: 10 × 5, 10 × 6 в
#: × и /, 10 × 14 в /). Для чужого класса оператор теперь не отказывает сам и
#: не отвечает NotImplemented — он зовёт `dispatch.combine` и возвращает её
#: ответ; отказ двери (пустая ячейка) проходит насквозь как есть. Ввоз
#: ленивый, внутри функции: pi_value — модуль ядра, Фios ввозит его сам, и
#: ввоз на уровне модуля дал бы кольцо.
#:
#: Кольцо ВЫЗОВОВ: правила `delegate`/`delegate_r` двери зовут оператор
#: носителя обратно с ТЕМИ ЖЕ операндами — так в ячейках (E10, E1, «×/»), а
#: класс 1 несёт не только int/Q, но и контейнер PhiNumberV6, которого _NUM
#: не покрывает (замер 09.10: `PI * PhiNumberV6.from_int(2)` — TypeError,
#: `route` — отказ «операторы носителей не сводят пару»). Второй вход с той
#: же парой (id, id, действие) отвечает NotImplemented, и дверь отказывает
#: сама, поимённо: глубина — две двери, не бесконечность. Прочие правила
#: таблицы для класса 10 × X оператор монома не зовут (post_check:mixed —
#: `Mixed.combine`, alpha_lift — операторы Φα, zero — оператор нуля; grep
#: RULE_IMPL в Фios.py). Один поток: множество не защищено замком.
_IN_DOOR: set = set()

#: Ячейки (класс другого операнда, действие), в которых дверь НЕ зовётся:
#: приёмки прибили в них отказ оператора поимённо, и по правилу главной сессии
#: 09.10 («краснеет существующая приёмка — снять правку до зелёного, приёмку
#: не трогать») правка снята ровно в них. Дверь в тех же ячейках даёт запись
#: класса 19 — расхождение названо, не спрятано; снятие — вместе с пересмотром
#: приёмок, не раньше. Места `[вычислено здесь]` 09.10:
#:   (1, ±)  — proto/test_doors.py:126–130 (одна подсказка «одноимённой степени»
#:             с обеих сторон), proto/test_phicore_carriers.py:99,
#:             proto/test_phicore_zeta.py:398–400 (записи C=10 и C=1 в отказе);
#:   (3, ±)  — proto/test_phicore_gamma_dual.py:409 (проверка G: записи C=3 и
#:             C=10 в отказе обеими сторонами);
#:   (14, ×) — proto/test_phicore_gamma.py:197 (π·Γ(1/2) — отказ);
#:   (16, все четыре) — proto/test_phicore_bessel_dual.py:219–230 (проверка B:
#:             восемь форм, текст «не сведено» с записями C=16 и C=10).
_PINNED = frozenset({
    (1, "+"), (1, "-"),
    (3, "+"), (3, "-"),
    (14, "*"),
    (16, "+"), (16, "-"), (16, "*"), (16, "/"),
})

#: метка «чужой класс: точного правила монома нет» — отлична от
#: NotImplemented, который для нуля-состояния есть законный ответ
_FOREIGN = object()


def _pinned(other, op: str) -> bool:
    """Ячейка (класс другого операнда, действие) прибита приёмкой — см. _PINNED."""
    from .Фios import classify, EXACT      # ленивый ввоз: Фios ввозит pi_value
    k = classify(other)
    return k is not None and k[0] == EXACT and (k[1], op) in _PINNED


def _door(a, b, op: str):
    """Дверь ФIOS для чужого класса; второй вход с той же парой — NotImplemented
    (кольцо delegate → оператор → дверь разорвано, см. _IN_DOOR)."""
    from .dispatch import combine          # ленивый ввоз: Фios ввозит pi_value
    key = (id(a), id(b), op)
    if key in _IN_DOOR:
        return NotImplemented
    _IN_DOOR.add(key)
    try:
        return combine(a, b, op)
    finally:
        _IN_DOOR.discard(key)


class PiMonomial:
    """q·π^e: coeff ≠ 0 (Q | AlgebraicNumber), e ≠ 0 (Q)."""

    __slots__ = ("coeff", "e")

    def __init__(self, coeff, e):
        if isinstance(coeff, float) or isinstance(e, float):
            raise TypeError("float запрещён построением")
        if isinstance(coeff, _NUM):
            coeff = Q(coeff)
            if coeff == 0:
                raise ValueError("coeff=0 не есть моном — ноль живёт в "
                                 "zero_state")
        elif not isinstance(coeff, AlgebraicNumber):
            raise TypeError("coeff — Q или AlgebraicNumber")
        e = Q(e)
        if e == 0:
            raise ValueError("e=0 — обычное число, не моном от π")
        self.coeff, self.e = coeff, e

    # --- сборка результата: вырождение e=0 в сам коэффициент ---
    @classmethod
    def _make(cls, coeff, e: Q):
        if isinstance(coeff, Q) and coeff == 0:
            return Q(0)
        if e == 0:
            return coeff
        return cls(coeff, e)

    # --- умножение / деление: точные правила впереди, чужой класс — дверь ---
    def _mul_exact(self, other):
        """Точные правила произведения (коммутативны, потому одни на обе
        стороны); _FOREIGN — чужой класс, решает дверь."""
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__rmul__ / Zero.__mul__
        if isinstance(other, PiMonomial):
            return self._make(self.coeff * other.coeff, self.e + other.e)
        if isinstance(other, _NUM):
            q = Q(other)
            if q == 0:
                return Q(0)
            return self._make(self.coeff * q, self.e)
        if isinstance(other, AlgebraicNumber):
            return self._make(self.coeff * other, self.e)
        return _FOREIGN

    def __mul__(self, other):
        r = self._mul_exact(other)
        if r is not _FOREIGN:
            return r
        if _pinned(other, "*"):
            return NotImplemented   # слово другого носителя (его отражённый оператор)
        return _door(self, other, "*")

    def __rmul__(self, other):
        r = self._mul_exact(other)
        if r is not _FOREIGN:
            return r
        if _pinned(other, "*"):
            return NotImplemented
        return _door(other, self, "*")  # порядок операндов — как в выражении

    def inverse(self) -> "PiMonomial":
        c = (self.coeff.inverse() if isinstance(self.coeff, AlgebraicNumber)
             else 1 / self.coeff)
        return PiMonomial(c, -self.e)

    def __truediv__(self, other):
        if isinstance(other, PiMonomial):
            return self * other.inverse()
        if isinstance(other, _NUM):
            q = Q(other)
            if q == 0:
                raise Unresolved("/", self.record(), _rec(other),
                                 hint="x/0 здесь не сведено")
            return self._make(self.coeff / q, self.e)
        if isinstance(other, AlgebraicNumber):
            return self._make(self.coeff * other.inverse(), self.e)
        if _pinned(other, "/"):
            return NotImplemented
        # ноль-состояние справа тоже идёт дверью: у Zero отражённого деления
        # нет, и прежде π / 0-состояние падал TypeError; таблица отказывает
        # поимённо («x/0 здесь не сведено»)
        return _door(self, other, "/")

    def __rtruediv__(self, other):
        if isinstance(other, (_NUM + (AlgebraicNumber,))):
            inv = self.inverse()
            return inv * other
        if _pinned(other, "/"):
            return NotImplemented
        return _door(other, self, "/")

    def __pow__(self, n: int):
        if not isinstance(n, int) or isinstance(n, bool):
            return NotImplemented
        if n == 0:
            return Q(1)
        c = self.coeff
        out = c
        for _ in range(abs(n) - 1):
            out = out * c
        if n < 0:
            out = (out.inverse() if isinstance(out, AlgebraicNumber)
                   else 1 / out)
        return self._make(out, self.e * n)

    def __neg__(self) -> "PiMonomial":
        return PiMonomial(-self.coeff, self.e)

    # --- сложение: точно — только одноимённая степень; чужой класс — дверь ---
    def _addsub(self, other, op: str, sign: int):
        if isinstance(other, Zero):
            return NotImplemented  # сведёт Zero.__radd__/__rsub__
        if _is_zero_num(other):
            return self             # x ± 0 = x: нуль класса 1 поглощается
        if isinstance(other, PiMonomial) and other.e == self.e:
            s = self.coeff + sign * other.coeff
            # алгебраическая сумма, равная нулю, уже распознана как Q
            if isinstance(s, Q) and s == 0:
                return Q(0)
            return self._make(s, self.e)
        # одноимённый класс с разной степенью — отказ самого носителя: ячейка
        # (E10, E10, «±») таблицы есть post_check:hyper → delegate, и дверь
        # звала бы этот же оператор обратно (кольцо); прибитые ячейки — _PINNED
        if isinstance(other, PiMonomial) or _pinned(other, op):
            raise Unresolved(op, self.record(), _rec(other), hint=_ADD_HINT)
        return _door(self, other, op)

    def __add__(self, other):
        return self._addsub(other, "+", +1)

    def __sub__(self, other):
        return self._addsub(other, "-", -1)

    def __radd__(self, other):
        if _is_zero_num(other):
            return self             # 0 + x = x
        if _pinned(other, "+"):
            raise Unresolved("+", _rec(other), self.record(), hint=_ADD_HINT)
        return _door(other, self, "+")

    def __rsub__(self, other):
        if _is_zero_num(other):
            return -self            # 0 − x = −x
        if _pinned(other, "-"):
            raise Unresolved("-", _rec(other), self.record(), hint=_ADD_HINT)
        return _door(other, self, "-")

    # --- запись, сравнение, сертификат ---
    def record(self) -> str:
        c = (self.coeff.record() if isinstance(self.coeff, AlgebraicNumber)
             else str(self.coeff))
        return f"Φ(C=10 TRIGONOMETRIC, вид=q·π^e, q={c}, e={self.e})"

    __repr__ = record

    def __eq__(self, other):
        if isinstance(other, PiMonomial):
            return self.e == other.e and self.coeff == other.coeff
        if isinstance(other, (_NUM + (AlgebraicNumber,))):
            return False  # π трансцендентно: моном не равен алгебраическому
        return NotImplemented

    def __hash__(self):
        return hash((PiMonomial, self.e))

    def self_check(self, deep: bool = False) -> list:
        """Самопроверка монома от π ВСТРЕЧНОЙ дорогой (ось самоприменения).

        `[вычислено здесь]` Сертификат сетки нельзя проверять тем же рядом,
        которым он построен, поэтому встречная дорога здесь ВНЕШНЯЯ:
        `[внешний источник]` вилка Архимеда—Метия 333/106 < π < 22/7
        (Метий, 1585; нижняя — подходящая дробь цепной дроби π). Она
        возводится в целый показатель целыми средствами линейки 1 (при
        e > 0 степень монотонно растёт, при e < 0 — убывает, потому концы
        меняются местами), умножается на коэффициент со сменой концов при
        отрицательном коэффициенте, и СЕРТИФИКАТ ОБЯЗАН ЛЕЖАТЬ ВНУТРИ —
        ряд Мэчина об этой вилке ничего не знает. Вторая дорога: уточнение
        вкладывается — `enclosure(terms + 8) ⊆ enclosure(terms)` и не шире.
        Третья: контракт eq/hash и запись, называющая класс 10 и вид q·π^e.

        Граница объявлена: при ДРОБНОМ показателе и при АЛГЕБРАИЧЕСКОМ
        коэффициенте внешняя вилка прямо не возводится (нужен корень из
        отрезка и уточнение носителя), и дорога вилки не идёт — остаются
        вложенность и запись; молчание вилки там находкой НЕ считается.
        Возвращает список находок; исключений не бросает.
        """
        out = []
        if isinstance(self.coeff, Q) and self.coeff == 0:
            out.append("coeff = 0: ноль записан мономом — ноль живёт "
                       "состоянием (zero_state)")
        if not isinstance(self.e, Q):
            out.append(f"показатель не пара: {self.e!r}")
            return out
        if self.e == 0:
            out.append("e = 0: обычное число записано мономом от π")
        rec = self.record()
        if "C=10 TRIGONOMETRIC" not in rec or "q·π^e" not in rec:
            out.append(f"запись не называет класс и вид: {rec}")
        if hash(self) != hash(PiMonomial(self.coeff, self.e)):
            out.append("контракт eq/hash нарушен: равные мономы — разные ключи")
        if out:
            return out
        try:
            iv = self.enclosure(24)
        except Exception as e:                                  # noqa: BLE001
            out.append(f"сертификат отказал: {type(e).__name__}: {e}")
            return out
        if iv.lo > iv.hi:
            out.append(f"сертификат пуст: {iv.record()}")
        if isinstance(self.coeff, Q) and self.e.denominator == 1:
            lo_pi, hi_pi = Q(333, 106), Q(22, 7)              # [внешний источник]
            k = self.e.numerator
            if k > 0:
                v_lo, v_hi = lo_pi ** k, hi_pi ** k
            else:
                v_lo, v_hi = Q(1) / hi_pi ** (-k), Q(1) / lo_pi ** (-k)
            if self.coeff > 0:
                b_lo, b_hi = self.coeff * v_lo, self.coeff * v_hi
            else:
                b_lo, b_hi = self.coeff * v_hi, self.coeff * v_lo
            if not (b_lo <= iv.lo and iv.hi <= b_hi):
                out.append(f"сертификат вышел за внешнюю вилку "
                           f"⟦{b_lo}; {b_hi}⟧ (333/106 < π < 22/7): {iv.record()}")
        if not deep:
            return out
        try:
            fine = self.enclosure(32)
        except Exception as e:                                  # noqa: BLE001
            out.append(f"уточнение отказало: {type(e).__name__}: {e}")
            return out
        if not iv.encloses(fine):
            out.append(f"уточнение не вложено в грубый сертификат: "
                       f"{fine.record()} ⊄ {iv.record()}")
        if fine.width() > iv.width():
            out.append("уточнение шире грубого сертификата")
        return out

    def enclosure(self, terms: int = 30) -> Interval:
        """Сертификат сетки [lo, hi] ⊂ ℚ; уточнение — явным запросом."""
        pi_iv = pi_enclosure(terms)
        if self.e.denominator != 1:
            # дробная степень: корень из отрезка (монотонность), сетка не
            # грубее сертификата π (Мэчин ≈ 4,6 бита на член)
            pi_iv = root_enclosure(pi_iv, self.e.denominator,
                                   bits=max(64, 5 * terms))
        iv = pi_iv ** self.e.numerator
        if isinstance(self.coeff, Q):
            return iv * self.coeff
        # алгебраический коэффициент: его изолирующий отрезок — вложение
        a = self.coeff
        for _ in range(terms):
            a.refine()
        return iv * Interval(a.lo, a.hi)


#: π — единица носителя.
PI = PiMonomial(Q(1), Q(1))
