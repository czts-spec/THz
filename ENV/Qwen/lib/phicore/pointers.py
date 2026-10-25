"""Указатели масштаба сетки: классы 21/22 с полем порядка k (ω^k).

HYPER — указатели масштаба, не числа (консенсус корпуса, 00-ЛИНИЯ п.2.2):
coeff — множитель шага сетки, k — поле порядка. Прямых вычислений в зоне
нет: сложение с чужим масштабом или с обычным числом — отказ Unresolved;
перемасштабирование — по явному запросу.

Ассоциативность умножения держится ТОЖДЕСТВЕННО, по построению:
произведение сведено к сложению целых k и умножению coeff (Q),
и то и другое ассоциативно. Трилемма (чтение A9, §17.2):
идемпотентность ЭЛЕМЕНТА отпущена (ω·ω ≠ ω), идемпотентность КЛАССА
сохранена — ω² остаётся классом 22 (k=2 > 0), сторона масштаба при
умножении внутри стороны не меняется.

Ограничения:
  * никаких float; coeff — только int/Q. float в любой из 10 дверей (+ − · /
    ** в обоих порядках) — TypeError с запретом построением, теми же словами,
    что у поля порядка: 9 дверей из 10 — «float запрещён построением»
    (hyper._lift), дверь `float ** указатель` — «основание степени — не bool
    и не float (запрет построением)» (hyper.__rpow__). Замер 18.09
    `[вычислено здесь]` на 5 указателях × 10 дверей = 50: прежде запрет
    называли 10 (обе двери `**`), 40 отдавали float интерпретатору без слова
    («unsupported operand»); теперь 50 из 50 (45 + 5). Ф-записи у float нет
    (класс 23 хранит его битово), потому запрет назван словом, а не отказом
    Unresolved с записью;
  * k == 0 в публичном конструкторе запрещён: k=0 — обычное число;
  * все значения в сообщениях/record() — явной Ф-записью.

Свип 11.09 (очередь п. 39d, приёмка test_hyper_bridge, 2630 случаев): где
арифметика указателя не определена, ответ — Unresolved с обеими Ф-записями,
не TypeError интерпретатора: деление с ноль-состоянием (у Zero деления нет)
и степень с показателем-НЕвеличиной (b^ω, ω^ω — показательная функция, её в
кольце нет). Показатель — целое в любой одежде (int, Q со знаменателем 1,
контейнер-число с целым значением); целая константа Hyper в показателе — за
полем порядка (Hyper.__rpow__, результат остаётся Hyper). Дробный показатель
— тоже в любой одежде (Q, контейнер-число, константа Hyper) — уходит в поле
порядка, и о нём решает теорема кольца, а не дверь (hyper.rational_power:
h^(p/q) в ℚ[ω, ω⁻¹] ⟺ h есть q-я степень; для одночлена c·ω^k — q | k и c
есть q-я степень в ℚ); результат — Hyper, отказ называет указатель КАК ДАН.
Замер 18.09 `[вычислено здесь]`: прежде одежда Q отказывала, а константа
Hyper считала — (4ω²)^(1/2) давало отказ при `** Q(1, 2)` и 2·ω при
`** Hyper.const(Q(1, 2))`, две двери одного действия; на 42 парах
(указатель, p/q) двери сходились 0 раз, теперь 42 из 42, и обе сходятся с
явным rational_power (42 из 42).

Указатель ⊕ класс 23 (ieee.Ieee) — с 18.09 `[вычислено здесь]`: у класса 23
правило есть для Hyper (точная сторона — величина кольца), для указателя не
было — 200 из 200 случаев (5 указателей × 5 float-значений × 4 действия × 2
порядка) кончались TypeError интерпретатора при двух гражданах пакета.
Теперь указатель поднимается в Hyper (_to_hyper), и считает класс 23 по
своему правилу (липкий: результат Ieee, точная сторона несёт 1·ω + x), оба
порядка одно; дорога указателя сходится с дорогой Hyper 200 из 200. Это то
же чтение элемента, что у OMEGA + Hyper.const(5) = 1·ω + 5: решает тот, у
кого правило. Ввоза нет: ieee ввозит hyper, класс 23 узнаётся по имени.

Остаток 11.09 (очередь п. 39, замечание (c) «три семантики ω + 5»): третья
пара — указатель ⊕ контейнер v6 — прежде TypeError интерпретатора (160 из
160 случаев: 5 указателей × 4 контейнера × 4 действия × 2 порядка). Контейнер
держит то же чтение указателя (00-ЛИНИЯ §2.2), потому указатель поднимается в
него (to_container: coeff·ω^k → from_hyper(n, d, k)) и считает контейнер по
своим правилам — ω + 5 = ω с поглощённым 5, ω − ω = Центр 0, void уступает,
класс 23 липкий; результат — контейнер, оба порядка одно. Сам указатель с
числом по-прежнему отказывает (прямых вычислений в зоне нет) — правило у
контейнера, не у указателя. x/0 через контейнер — тот же отказ поимённо, что
ScalePointer / 0 (у контейнера — ZeroDivisionError). Показатель-контейнер с
целым значением — степень (OMEGA ** контейнер 5 = ω^5), иначе отказ поимённо.
Приёмка: test_hyper_bridge (свип 200 + 25) и test_phicore (набор указателей).
"""

import operator

from .rational import Q

from .errors import Unresolved
from .zero_state import Zero

#: Допустимые скаляры. float исключён построением.
_NUM = (int, Q)

#: `[вычислено здесь]` Прежде здесь стояло «смена масштаба по явному
#: запросу» — обещание действия, которого не было ни в `ScalePointer`, ни
#: в пакете (очередь п. 42). Действие с 05.09 есть, и подсказка называет
#: ЕГО, а не намерение: указатель несёт ОДИН масштаб и суммы нести не
#: может, сумма живёт в поле порядка, куда указатель поднимается мостом.
_HINT = ("в HYPER-зоне прямых вычислений нет: указатель несёт один масштаб, "
         "сумма живёт в поле порядка — hyper.Hyper.const(0) + x + y (ω + 5 = 1·ω + 5)")
#: Те же слова, что у поля порядка (hyper.py): ω/0 и ω^ω отвечают одинаково
#: в обоих чтениях.
_X0_HINT = "x/0 здесь не сведено"
_ZERO_DIV_HINT = "деление с ноль-состоянием не сведено: правила деления у Zero нет"
#: `[вычислено здесь]` 18.09: прежний текст называл ω^(1/2) «дробной степенью
#: вне кольца ℚ[ω, ω⁻¹]» — о кольце вообще, а это верно НЕ ВСЕГДА (hyper.py:
#: (4ω²)^(1/2) = 2ω внутри, решает теорема); дробный показатель теперь уходит
#: в поле порядка, и здесь остаётся только показатель-НЕвеличина.
_POW_HINT = ("показатель степени — целое (int, Q со знаменателем 1, целая константа "
             "Hyper, контейнер-число) либо рациональное (дробный показатель уходит в "
             "поле порядка: о нём решает теорема кольца — hyper.rational_power, "
             "результат Hyper); здесь же показателя-величины нет вовсе — b^ω и ω^ω "
             "суть показательная функция, а её в кольце ℚ[ω, ω⁻¹] нет")


def _citizen(x: object) -> bool:
    """Гражданин пакета (объект из phicore или целое): отказ ему — поимённо."""
    if isinstance(x, bool):
        return False
    return isinstance(x, int) or type(x).__module__.startswith("phicore")


def _is_container(x: object) -> bool:
    """Контейнер v6 (эталон семантики, чтение указателя) — по имени и фабрике
    from_hyper, без ввоза: контейнер ввозит указатели, не наоборот."""
    return type(x).__name__ == "PhiNumberV6" and callable(getattr(type(x), "from_hyper", None))


def _is_ieee(x: object) -> bool:
    """Класс 23 (`ieee.Ieee`) — по имени и паре полей raw/exact, без ввоза:
    ieee ввозит hyper, указатели никого не ввозят (hyper — лениво, в методах)."""
    return type(x).__name__ == "Ieee" and hasattr(x, "raw") and hasattr(x, "exact")


def _reject_float(x: object) -> None:
    """float — запрет построением, и запрет НАЗВАН той же фразой, что у поля
    порядка (hyper._lift): NotImplemented отдал бы его интерпретатору без слова
    (замер 18.09: 40 дверей из 50 молчали, см. шапку)."""
    if isinstance(x, float):
        raise TypeError("float запрещён построением")


#: Действия для подъёма в поле порядка (указатель ⊕ класс 23) — по имени действия.
_OPS = {"+": operator.add, "-": operator.sub, "*": operator.mul, "/": operator.truediv}


def _phi(x: object) -> str:
    """Ф-запись операнда для сообщений отказа."""
    rec = getattr(x, "record", None)
    if callable(rec):
        return rec()
    if isinstance(x, _NUM):
        return f"Φ(C=1 RATIONAL, value={Q(x)})"
    return str(x)


class ScalePointer:
    """Указатель масштаба: coeff ≠ 0 (Q), k ≠ 0 (int).

    Счётчик `k` и экспонента — одно и то же: слово автора `[в корпусе]`
    (`inbox/SKILL_op_geometry_phi_numbers.md:264`, «Ключевые инсайты», п. 1):
    «Экспонента e и счётчик k — одно и то же (полярная структура, тавтология
    колеса)». Тот же счётчик в ПОЗИЦИОННОЙ сетке есть НОМЕР КОЛЕСА
    (`wheels.Wheel.index`): целое колесо k = e + 1 — указатель 22, дробное
    k = −e — указатель 21. Одно чтение при разных масштабах: ω в поле
    порядка, основание b в позиционной сетке. Тег колеса при этом НЕ есть
    указатель ω/ε: сумма колёс живёт в классе 1, и строка (E21, E1, «+»)
    таблицы ФIOS — refuse (приёмка `proto/test_phicore_wheels.py`).
    """

    __slots__ = ("coeff", "k")

    def __init__(self, coeff, k: int):
        if not isinstance(coeff, _NUM):
            raise TypeError("coeff — только int/Q, float запрещён")
        if not isinstance(k, int) or isinstance(k, bool):
            raise TypeError("k — только int")
        coeff = Q(coeff)
        if coeff == 0:
            raise ValueError("coeff=0 не есть указатель масштаба")
        if k == 0:
            raise ValueError("k=0 — обычное число, не указатель")
        self.coeff = coeff
        self.k = k

    # --- мост в контейнер v6 (очередь п. 39, остаток (c)/(d), 11.09) -------

    def to_container(self, cls) -> object:
        """Подъём в контейнер того же чтения: coeff·ω^k → cls.from_hyper(n, d, k).
        Контейнер держит чтение указателя (00-ЛИНИЯ §2.2), потому указатель ⊕
        контейнер считает контейнер по своим правилам (ω + 5 = ω, поглощение;
        состояния и класс 23 — свои), и результат остаётся контейнером — в
        обоих порядках одно. Чтение элемента (1·ω + 5) — мост в Hyper."""
        return cls.from_hyper(int(self.coeff.n), int(self.coeff.d), self.k)

    def _via_container(self, other: object, op: str, self_left: bool) -> object:
        """Указатель ⊕ контейнер: поднять себя в контейнер и отдать ему счёт.
        x/0 у контейнера — ZeroDivisionError интерпретатора; здесь тот же
        отказ, что ScalePointer / 0: поимённо, с обеими записями."""
        me = self.to_container(type(other))
        a, b = (me, other) if self_left else (other, me)
        if op == "+":
            return a + b
        if op == "-":
            return a - b
        if op == "*":
            return a * b
        try:
            return a / b
        except ZeroDivisionError:
            left, right = (self.record(), _phi(other)) if self_left else (_phi(other), self.record())
            raise Unresolved("/", left, right, hint=_X0_HINT) from None   # операнды как даны

    def _to_hyper(self) -> object:
        """Подъём в поле порядка (чтение элемента): coeff·ω^k → Hyper({k: coeff}).
        Ввоз ленивый: hyper указатели не ввозит, указатели hyper — только здесь."""
        from .hyper import Hyper
        return Hyper({self.k: self.coeff})

    def _via_hyper(self, other: object, op: str, self_left: bool) -> object:
        """Указатель ⊕ класс 23: поднять себя в Hyper и отдать счёт классу 23 —
        у него правило есть для Hyper (точная сторона — величина кольца), для
        указателя нет; результат Ieee (класс 23 липкий), оба порядка одно."""
        me = self._to_hyper()
        a, b = (me, other) if self_left else (other, me)
        return _OPS[op](a, b)

    # --- умножение / деление -------------------------------------------

    def __mul__(self, other: object) -> object:
        _reject_float(other)
        if _is_container(other):
            return self._via_container(other, "*", True)
        if _is_ieee(other):
            return self._via_hyper(other, "*", True)
        if isinstance(other, ScalePointer):
            c = self.coeff * other.coeff
            kk = self.k + other.k
            if kk == 0:
                return c  # ω·ε = 1 — аксиома A6: указатель вырождается в число
            return ScalePointer(c, kk)
        if isinstance(other, _NUM):
            q = Q(other)
            if q == 0:
                return Q(0)  # ε×0 = обычный ноль (правило корпуса)
            return ScalePointer(self.coeff * q, self.k)
        return NotImplemented

    def __rmul__(self, other: object) -> object:
        if _is_container(other):
            return self._via_container(other, "*", False)
        if _is_ieee(other):
            return self._via_hyper(other, "*", False)
        return self.__mul__(other)                      # умножение на скаляр коммутативно

    def __truediv__(self, other: object) -> object:
        _reject_float(other)
        if _is_container(other):
            return self._via_container(other, "/", True)
        if _is_ieee(other):
            return self._via_hyper(other, "/", True)
        if isinstance(other, ScalePointer):
            c = self.coeff / other.coeff
            kk = self.k - other.k
            if kk == 0:
                return c
            return ScalePointer(c, kk)
        if isinstance(other, _NUM):
            q = Q(other)
            if q == 0:
                # x/0 — значение-зависимое правило post-check (00-ЛИНИЯ п.2.4)
                raise Unresolved("/", self.record(), _phi(other), hint=_X0_HINT)
            return ScalePointer(self.coeff / q, self.k)
        if isinstance(other, Zero):
            # у Zero нет __rtruediv__ — решать некому; отказ поимённо, не TypeError
            raise Unresolved("/", self.record(), other.record(), hint=_X0_HINT)
        return NotImplemented

    def __rtruediv__(self, other: object) -> object:
        """число / указатель = число · inverse(); 0/указатель = обычный ноль."""
        _reject_float(other)
        if _is_container(other):
            return self._via_container(other, "/", False)
        if _is_ieee(other):
            return self._via_hyper(other, "/", False)
        if isinstance(other, _NUM):
            q = Q(other)
            if q == 0:
                return Q(0)
            return ScalePointer(q / self.coeff, -self.k)
        if isinstance(other, Zero):
            raise Unresolved("/", other.record(), self.record(), hint=_ZERO_DIV_HINT)
        return NotImplemented

    def inverse(self) -> "ScalePointer":
        """Обратный указатель. −k ≠ 0 раз k ≠ 0 — всегда указатель,
        вырождения в число нет."""
        return ScalePointer(1 / self.coeff, -self.k)

    def __pow__(self, n: object) -> object:
        if isinstance(n, (bool, float)):
            raise TypeError("степень — целое (bool не число, float запрещён построением)")
        if isinstance(n, Q) and n.d == 1:
            n = int(n.n)
        if not isinstance(n, int):
            if type(n).__name__ == "Hyper":
                return NotImplemented   # поле порядка решает: целая константа — степень, дробная — теорема кольца, иначе отказ поимённо
            # одна одежда показателя на пакет (hyper указатели не ввозит — ввоз ленивый):
            # Q, контейнер-число с точным отношением — величина; иное — показателя-величины нет
            from .hyper import _rational_exponent, rational_power
            e = _rational_exponent(n)
            if e is not None:
                if e.d == 1:
                    return self ** int(e.n)             # целое в одежде контейнера: ω ** контейнер 5 = ω^5
                # дробный показатель есть ВЕЛИЧИНА: решает теорема кольца, не дверь — та же
                # дорога, что у `** Hyper.const(Q(1, 2))` (две двери одного действия свёрнуты
                # в одну, 18.09); результат — Hyper, отказ называет указатель КАК ДАН
                return rational_power(self._to_hyper(), e, given=n, given_base=self)
            if _citizen(n):
                raise Unresolved("**", self.record(), _phi(n), hint=_POW_HINT)
            return NotImplemented
        if n == 0:
            return Q(1)  # k·0 = 0 — обычное число
        return ScalePointer(self.coeff ** n, self.k * n)

    def __rpow__(self, base: object) -> object:
        """base ** указатель — показательная функция, вне зоны: отказ поимённо.
        Hyper-основание сюда не доходит (его __pow__ отвечает первым)."""
        if isinstance(base, (bool, float)):
            raise TypeError("основание степени — не bool и не float (запрет построением)")
        if _citizen(base):
            raise Unresolved("**", _phi(base), self.record(), hint=_POW_HINT)
        return NotImplemented

    # --- сложение / вычитание ------------------------------------------

    def _combine(self, other: object, op: str, sign: int) -> object:
        _reject_float(other)
        if isinstance(other, Zero):
            return NotImplemented  # обработает Zero.__radd__/__rsub__
        if _is_container(other):
            return self._via_container(other, op, True)
        if _is_ieee(other):
            return self._via_hyper(other, op, True)
        if isinstance(other, ScalePointer):
            if other.k == self.k:
                s = self.coeff + sign * other.coeff
                if s == 0:
                    # род нуля по происхождению: ω−ω / ε−ε
                    if self.k > 0:
                        return Zero.structural_infinite()
                    return Zero.structural_infinitesimal()
                return ScalePointer(s, self.k)
            raise Unresolved(op, self.record(), other.record(), hint=_HINT)
        if isinstance(other, _NUM):
            raise Unresolved(op, self.record(), _phi(other), hint=_HINT)
        return NotImplemented

    def __add__(self, other: object) -> object:
        return self._combine(other, "+", +1)

    def __sub__(self, other: object) -> object:
        return self._combine(other, "-", -1)

    def __radd__(self, other: object) -> object:
        _reject_float(other)
        if _is_container(other):
            return self._via_container(other, "+", False)
        if _is_ieee(other):
            return self._via_hyper(other, "+", False)
        if isinstance(other, _NUM):
            raise Unresolved("+", _phi(other), self.record(), hint=_HINT)
        return NotImplemented

    def __rsub__(self, other: object) -> object:
        _reject_float(other)
        if _is_container(other):
            return self._via_container(other, "-", False)
        if _is_ieee(other):
            return self._via_hyper(other, "-", False)
        if isinstance(other, _NUM):
            raise Unresolved("-", _phi(other), self.record(), hint=_HINT)
        return NotImplemented

    def __neg__(self) -> "ScalePointer":
        return ScalePointer(-self.coeff, self.k)

    # --- класс, запись, сравнение --------------------------------------

    def phi_class(self) -> int:
        """22 — сторона ω (k>0), 21 — сторона ε (k<0)."""
        return 22 if self.k > 0 else 21

    def self_check(self, deep: bool = False) -> list:
        """Самопроверка указателя ВСТРЕЧНОЙ дорогой (ось самоприменения).

        `[вычислено здесь]` Дороги: (1) класс по знаку порядка (22 ⟺ k > 0,
        21 ⟺ k < 0) — против того, что печатает запись; (2) мост в
        СОСЕДНИЙ модуль: подъём `to_container(PhiNumberV6)` обязан дать тот
        же знаковый порядок и тот же коэффициент, а контейнер держит порядок
        ДРУГИМ устройством — классом 21/22 и парой e1/e2, не полем k; (3)
        мост в ТРЕТЬЕ чтение — поле порядка `hyper.Hyper` (кольцо
        многочленов Лорана): одночлен {k: coeff} несёт тот же порядок, и
        произведение двух указателей, посчитанное там, обязано совпасть с
        посчитанным здесь; (4) контракт eq/hash.

        Граница: указатель есть ЗАПИСЬ масштаба, не число (coeff ≠ 0 и
        k ≠ 0 — построением), поэтому дорога (3) сверяет только порядок и
        коэффициент; сложение указателя с числом — чтение указателя
        (поглощение) и сверяется не здесь, а приёмкой моста. Возвращает
        список находок; исключений не бросает.
        """
        out = []
        if not isinstance(self.k, int) or isinstance(self.k, bool):
            out.append(f"порядок k не целое: {self.k!r}")
            return out
        if self.k == 0:
            out.append("k = 0: обычное число записано указателем масштаба")
        if not isinstance(self.coeff, Q):
            out.append(f"коэффициент не пара: {self.coeff!r}")
            return out
        if self.coeff == 0:
            out.append("coeff = 0: ноль записан указателем — ноль живёт "
                       "состоянием (zero_state)")
        if self.k and self.phi_class() != (22 if self.k > 0 else 21):
            out.append(f"класс {self.phi_class()} не отвечает порядку {self.k}")
        if out:
            return out
        try:
            from .container import PhiNumberV6
            c = self.to_container(PhiNumberV6)
            if c.hyper_order != self.k:
                out.append(f"мост в контейнер дал порядок {c.hyper_order}, "
                           f"записано {self.k}")
            if Q(c.m1.v, c.m2.v) != self.coeff:
                out.append(f"мост в контейнер дал коэффициент "
                           f"{Q(c.m1.v, c.m2.v)}, записано {self.coeff}")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"мост в контейнер отказал: {type(e).__name__}: {e}")
        try:
            from .hyper import Hyper
            h = Hyper({self.k: self.coeff})
            if h.order != self.k or h.terms.get(self.k) != self.coeff:
                out.append(f"мост в поле порядка дал {h.describe()}, "
                           f"записано {self.record()}")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"мост в поле порядка отказал: {type(e).__name__}: {e}")
        twin = ScalePointer(self.coeff, self.k)
        if twin != self or hash(twin) != hash(self):
            out.append("контракт eq/hash нарушен: равные указатели — разные ключи")
        if not deep:
            return out
        if self.inverse().inverse() != self:
            out.append("обращение дважды не вернуло указатель")
        try:
            from .hyper import Hyper
            other = ScalePointer(Q(3, 2), -self.k if self.k else 1)
            mine = self * other
            theirs = Hyper({self.k: self.coeff}) * Hyper({other.k: other.coeff})
            k_mine = mine.k if isinstance(mine, ScalePointer) else 0
            c_mine = mine.coeff if isinstance(mine, ScalePointer) else Q(mine)
            if theirs.order != k_mine or theirs.terms.get(k_mine) != c_mine:
                out.append(f"произведение двумя дорогами расходится: здесь "
                           f"{mine!r}, в поле порядка {theirs.describe()}")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"сверка произведения отказала: {type(e).__name__}: {e}")
        return out

    def record(self) -> str:
        return f"Φ(C={self.phi_class()}, coeff={self.coeff}, k={self.k})"

    def __repr__(self) -> str:
        return self.record()

    def __eq__(self, other: object) -> bool:
        if isinstance(other, ScalePointer):
            return self.coeff == other.coeff and self.k == other.k
        return NotImplemented

    def __hash__(self) -> int:
        return hash((ScalePointer, self.coeff, self.k))


#: ω — единичный шаг вверх по масштабу.
OMEGA = ScalePointer(Q(1), 1)
#: ε — единичный шаг вниз; ω·ε = 1 (A6).
EPSILON = ScalePointer(Q(1), -1)
