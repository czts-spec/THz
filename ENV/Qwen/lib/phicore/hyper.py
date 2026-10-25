"""Указатели масштаба с полем порядка: конечные суммы c_k·ω^k (ε = ω^−1).

Дефект корпуса (CONTEXT §4, STUDY-01 §4; воспроизведён на phicore.container
04.09): свип ассоциативности по 12 представителям классов давал 20 нарушений
для + и 32 для * — от двух свойств контейнера v6: у ω/ε не было порядка
(ω·ω = ω, ω·ε = 1 ⇒ (ω·ω)·ε ≠ ω·(ω·ε); снято 04.09 — порядок в e1, 32 → 0
вне класса 23) и одночленность (ε поглощается конечным: (5−5)+ε ≠ 5+(−5+ε);
остаётся — это чтение указателя, консенсус 00-ЛИНИЯ §2.2). Консенсус линии (00-ЛИНИЯ §2.6): «указателям
21/22 нужно то же поле порядка ω^k, что классу 4; выбор чтения A9 (класс
против элемента) за автором». Здесь — элемент с порядком, чтобы выбор был
замерен, а не угадан: Hyper = конечная сумма c_k·ω^k, k ∈ ℤ, c_k ∈ Q
(кольцо многочленов Лорана над ℚ — строка «Робинсон» таблицы аксиом:
сложение и умножение ассоциативны и коммутативны по построению).

Класс по старшему члену: k_max > 0 — INFINITY (22), k_max < 0 —
INFINITESIMAL (21), k_max = 0 — конечное, пусто — ноль. Ноль из сокращения
помнит происхождение (чтение v6), пока остаётся нулём: origin — запись о
пути, не величина, и в равенство величин не входит (замер 04.09: если
входит, ассоциативность ломается на тройках с нулями: 6 для + и 28 для · — путь
разный, величина одна); в сумме с ненулём он обычный 0. Обратный элемент — только у одночлена: обратный к сумме есть
бесконечный ряд, и он не представим конечно — отказ
ExactRepresentationNotAvailable, не приближение.

Контейнер v6 (эталон семантики) держит чтение указателя (00-ЛИНИЯ §2.2;
ω − ω — поглощающее ноль-состояние, слово автора); этот модуль — чтение
элемента, замеренная альтернатива. Таблица 2×2 (05.09, свип 1728 троек на
контейнере): указатель/поглощает 20, указатель/уступает 36, элемент/
поглощает 28, элемент/уступает 12 — у каждого чтения свой ноль, смешение
хуже обоих; на контейнере чтения не смешиваются. Мост 05.09: контейнер ⊕
Hyper считается здесь (чтение элемента липкое — результат остаётся Hyper;
_lift контейнер не поднимает: его состояния решают сами через отражённые
методы контейнера, число и указатель поднимаются PhiNumberV6._to_hyper).

Авторство: Робинсон 1966 (нестандартный анализ, поле порядка бесконечно
больших), Конвей 1976 (сюрреальные: ω^k как масштаб); классы 21/22 как
указатели масштаба и ноль-состояние — AIAΣA (автор проекта).

Свип 11.09 (очередь п. 39d, «Hyper и ScalePointer не считаются друг с
другом»): исчерпывающий перебор малых случаев в test_hyper_bridge — 9 Hyper ×
23 соседа (5 указателей, 5 чисел, 4 ноль-состояния, 9 Hyper) × 5 операций ×
2 порядка плюс зона указателей, 2630 случаев. До правки 668 исходов были
TypeError/ZeroDivisionError интерпретатора без Ф-записи (степень с чужим
показателем 486, деление с ноль-состоянием 112, x/0 70); после — 0: где
арифметика не определена, стоит Unresolved, называющий обе Ф-записи (record()
у Hyper — Φ(C=22 HYPER, 1·ω + 5)). Правило моста: гражданин пакета с
отражённым методом решает сам (NotImplemented); гражданин без него — отказ
поимённо; чужак (str, list) — TypeError интерпретатора; bool и float —
TypeError построением. x/0 — тот же отказ, что у указателя («x/0 здесь не
сведено»: правило post-check за контейнером, 00-ЛИНИЯ п.2.4). Степень:
показатель — целое, в любой одежде (int, Q со знаменателем 1, целая
константа Hyper); b^ω, ω^ω, ω^(1/2) — Unresolved: показательная функция и
дробная степень вне кольца ℚ[ω, ω⁻¹]. Уточнено замером 17.09 (см. п. 1 ниже):
«дробная степень вне кольца» верно НЕ ВСЕГДА — ω^(1/2) вне (порядок 1 не делится
на 2), а (4ω²)^(1/2) = 2ω внутри, и решает о том теорема, а не дверь.

Остаток 11.09 (повтор по очереди п. 39, приёмка test_hyper_bridge 43/0, свип
теперь в наборе: зона Hyper 2070, зона указателей 700, указатель ⊕ контейнер
200 + 25): отказ называет операнд КАК ДАН — 5 / 0{ω−ω} говорит Φ(C=1 RATIONAL,
value=5), OMEGA / 0{1} — Φ(C=22, coeff=1, k=1), не запись после подъёма (было
31 таких, 2 из них 0^(−k) — «1/x» без показателя; стало 0). Показатель —
ещё и контейнер-число с целым значением (ω ** контейнер 2 = ω²; √2 — отказ).
Третья пара трёх семантик ω + 5 (замечание (c)) — указатель ⊕ контейнер —
считается в pointers.py подъёмом указателя в контейнер (то же чтение
указателя, результат — контейнер, оба порядка одно; было 160 TypeError).
За контейнером (container.py, не правился): контейнер ** (Hyper | указатель)
— TypeError его двери «принимает int» (70 из 70), контейнер / Hyper-ноль —
отказ «1/x» из _mul_hyper без записи делимого (8).

Деление на многочлен замерено 11.09 `[вычислено здесь]` и вынесено отдельным
устройством `laurent_inverse` / `laurent_div`: обратный к сумме есть бесконечный
ряд Лорана, и конечной записи у него нет — потому `/` по-прежнему отказывает
(ложного класса не возникает), а разложение считается только по явному запросу,
с ОБЪЯВЛЕННЫМ порядком усечения и названным остатком. Итог всегда со
знаменателем: «частичная сумма до ω^m, остаток порядка ω^(m−1) и ниже».

NOT_MEASURED (запись 11.09) — РАЗОБРАН 17.09 на два РАЗНЫХ утверждения; строка
остаётся записью о том, что здесь было объявлено, а каждое утверждение закрыто
своим доказательством.

1. Степень с рациональным показателем — ЗАКРЫТО точной дорогой. Запись
   говорила «отказ поимённо, не вычисление», и замер 17.09 `[вычислено здесь]`
   её опроверг: устройство корня (`laurent_root`, `rational_power`) уже стояло,
   а ДВЕРЬ `**` его не звала — отказывала в 5 случаях из 6, где величина кольца
   ЕСТЬ ((4ω²)^(1/2) = 2ω, 9^(1/2) = 3, (4ε²)^(1/2) = 2ε, (4ω²)^(3/2) = 8ω³,
   (ω² + 2 + ε²)^(1/2) = ω + ε). Теперь дорогу берёт сама дверь и одинаково с
   обеих сторон (`__pow__`, `__rpow__`; 4 ** Hyper(1/2) = 2 — тот же счёт).
   Граница названа теоремой: h^(p/q) при несократимой p/q лежит в кольце ТОГДА
   И ТОЛЬКО ТОГДА, когда h есть q-я степень в нём. Где корня нет, отказ
   остаётся, но он О ПРЕДМЕТЕ (`_ROOT_HINT`): «не умею посчитать» — неверно
   (тем же q дом считает точно, когда h есть q-я степень), «величины не
   существует» — тоже неверно (√2 зажата вилкой в ℝ на сетке 2^−64, ω^(1/2)
   живёт у Леви-Чивита); верно третье — величина есть, конечной записи ЭТОГО
   кольца у неё нет. Две дороги приёмки: построение g^q арифметикой кольца
   против возврата рекуррентой — 32 пары (g, q); тождество (h^(p/q))^q = h^p
   целыми степенями — 96 троек.

2. Взаимодействие с классами 23 (IGNORANCE) и ZERO — ЗАКРЫТО замером, и запись
   «они остаются у контейнера» ПОПРАВЛЕНА: контейнера на пути нет ни в одном
   случае. Свип 17.09 `[вычислено здесь]` двумя независимыми дорогами (исход
   оператора против решающего модуля по протоколу Python) — 360 случаев, дороги
   сошлись 360 из 360. С ZERO 160: решает `zero_state` 120, сам Hyper 40;
   величин 90, отказов поимённо 70. С классом 23 — 200, и все 200 решает
   `ieee`, точная сторона несёт запись кольца (ω + 23(0.1) — точно 1·ω + 0.1).
   TypeError интерпретатора — 0 из 360. Правило моста прежнее: сосед решает сам
   отражённым методом, а не контейнер за него.

   ОТКАЗ ДОМА с ценой (та часть записи, что осталась верной по сути, но названа
   неточно). ПРЕДМЕТ: арифметика класса 23 на стороне hyper. ПРИЧИНА: порядок
   ввоза, не теорема — `ieee.py` ввозит `hyper`, потому назвать здесь класс 23
   значит завести круг ввоза. АДРЕС: `ieee.py` (`from .hyper import Hyper`),
   здешние `_lift` и `_rhs`. ЧЕМ МЕРИЛОСЬ БЫ: ценой слияния двух модулей в один
   либо ленивого ввоза внутри метода — оба меняют дом, а не счёт, и выигрыша в
   точности не дают вовсе, ибо счёт уже точен на стороне соседа (замер выше:
   200 из 200).
"""
from __future__ import annotations

from .enclosure import int_root_floor
from .errors import ExactRepresentationNotAvailable, Unresolved
from .rational import Q
from .zero_state import Zero

INFINITY = 22
INFINITESIMAL = 21
FINITE = 1
ZERO = 0


def _q(v) -> Q:
    return v if isinstance(v, Q) else Q(v)


#: Подсказки отказов — тот же текст, что у указателя (pointers.py), чтобы
#: ω/0 в обоих чтениях отвечал одним словом.
_X0_HINT = "x/0 здесь не сведено"
#: Показателя-ВЕЛИЧИНЫ здесь нет вовсе (ω, указатель, ноль-состояние): дробный
#: показатель сюда не попадает — его с 17.09 берёт сама дверь (`rational_power`).
_POW_HINT = ("показатель степени — целое (int, Q со знаменателем 1, целая константа "
             "Hyper) либо рациональное (его дверь считает сама: `rational_power`); "
             "здесь же показателя-величины нет вовсе — b^ω и ω^ω суть "
             "показательная функция, а её в кольце ℚ[ω, ω⁻¹] нет")
#: Дробный показатель: отказ О ПРЕДМЕТЕ с теоремой, а не «не умею». Замер 12.09
#: `[вычислено здесь]`: прежде ОДИН текст накрывал и показательную функцию, и
#: корень — и потому отказывал даже там, где корень точен ((4·ω²)^(1/2) = 2·ω).
_ROOT_HINT = ("корня степени {q} у этой записи в ℚ[ω, ω⁻¹] нет. Утверждение точное, "
              "не оценка: при gcd(p, q) = 1 величина h^(p/q) лежит в кольце ТОГДА И "
              "ТОЛЬКО ТОГДА, когда h есть q-я степень в нём, и извлечение "
              "(`laurent_root(h, {q})`) записи не дало. Величина при этом "
              "СУЩЕСТВУЕТ — в расширении (Леви-Чивита, сюрреальные Конвея) либо в ℝ "
              "для коэффициента; конечной записи этого кольца у неё нет")
_NO_RULE_HINT = ("у моста поля порядка правила для этого операнда нет, и отражённым "
                 "методом он не отвечает")
_ZERO_DIV_HINT = "деление с ноль-состоянием не сведено: правила деления у Zero нет"


def _phi(x) -> str:
    """Ф-запись операнда для отказа: record() гражданина, число — линейка 1."""
    rec = getattr(x, "record", None)
    if callable(rec):
        return rec()
    if isinstance(x, (int, Q)) and not isinstance(x, bool):
        return f"Φ(C=1 RATIONAL, value={Q(x)})"
    return repr(x)


def _citizen(x) -> bool:
    """Гражданин пакета: объект из phicore (Q, указатель, ноль-состояние,
    контейнер) или целое. Для него отказ обязан быть поимённым — Unresolved с
    обеими Ф-записями, не TypeError интерпретатора; чужак (str, list …)
    получает NotImplemented, и дверь закрывает интерпретатор."""
    if isinstance(x, bool):
        return False
    return isinstance(x, int) or type(x).__module__.startswith("phicore")


def _int_exponent(n):
    """Целый показатель в любой одежде: int, Q со знаменателем 1, целая
    константа Hyper (нулевой Hyper — 0). Иначе None."""
    if isinstance(n, bool):
        return None
    if isinstance(n, int):
        return n
    if isinstance(n, Q):
        return int(n.n) if n.d == 1 else None
    if isinstance(n, Hyper):
        if not n.terms:
            return 0
        if set(n.terms) == {0}:
            c = n.terms[0]
            return int(c.n) if c.d == 1 else None
    if _is_container(n) and n.hyper_order == 0 and not n.is_zero_class and not n.is_uncomputable:
        try:
            return _int_exponent(n._to_hyper())         # контейнер-число: подъём по точному отношению
        except ExactRepresentationNotAvailable:
            return None                                 # точного отношения нет (√2) — не целое
    return None


def _rational_exponent(n):
    """Рациональный показатель в любой одежде: int, Q, константа Hyper,
    контейнер-число. Иначе None — там показателя-величины нет вовсе (ω,
    указатель, ноль-состояние), и отвечает `_POW_HINT`.

    Отличие от `_int_exponent` ровно одно — знаменатель не обязан быть 1:
    дробный показатель есть ВЕЛИЧИНА, и решать о нём обязана теорема кольца
    (`rational_power`), а не дверь.
    """
    if isinstance(n, (bool, float)):
        return None
    if isinstance(n, int):
        return Q(n)
    if isinstance(n, Q):
        return n
    if isinstance(n, Hyper):
        if not n.terms:
            return Q(0)
        if set(n.terms) == {0}:
            return n.terms[0]
        return None
    if _is_container(n) and n.hyper_order == 0 and not n.is_zero_class and not n.is_uncomputable:
        try:
            return _rational_exponent(n._to_hyper())
        except ExactRepresentationNotAvailable:
            return None                     # точного отношения нет (√2) — не рационален
    return None


def _is_container(x) -> bool:
    """Контейнер v6 — по имени и полю порядка, без ввоза (hyper не ввозит контейнер)."""
    return type(x).__name__ == "PhiNumberV6" and hasattr(x, "hyper_order")


def _lift(x):
    """Мост в поле порядка (свежие глаза 05.09, глава 03, замечания 4, 15, 19:
    «Hyper и ScalePointer друг с другом не считаются», «три семантики ω + 5»).
    Число (int/Q) — константа; указатель масштаба (coeff, k) — одночлен
    coeff·ω^k; ноль-состояние НЕ поднимается: его правило (поглощение,
    нейтральность, отказ NOT_COMPUTED) решает сам Zero через отражённый
    метод, чтобы оба порядка операндов давали одно. float — запрет
    построением. Прочее — None: кому решать, смотрит Hyper._rhs."""
    if isinstance(x, Hyper):
        return x
    if isinstance(x, bool):
        return None
    if isinstance(x, float):
        raise TypeError("float запрещён построением")
    if isinstance(x, (int, Q)):
        return Hyper.const(x)
    if _is_container(x):
        return None        # контейнер v6 решает сам (отражённые методы): состояния — по своим правилам,
                           # число и указатель — подъём _to_hyper, результат-одночлен спускается в контейнер
    coeff, k = getattr(x, "coeff", None), getattr(x, "k", None)
    if coeff is not None and isinstance(k, int) and type(x).__name__ == "ScalePointer":
        return Hyper({k: coeff})
    return None


class Hyper:
    __slots__ = ("terms", "origin")

    def __init__(self, terms=None, origin: str | None = None):
        t = {}
        for k, c in (terms or {}).items():
            c = _q(c)
            if c != 0:
                t[int(k)] = c
        self.terms = t
        self.origin = origin if not t else None    # происхождение живёт только у нуля

    # --- фабрики
    @classmethod
    def omega(cls, coeff=1, order: int = 1) -> "Hyper":
        return cls({order: coeff})

    @classmethod
    def epsilon(cls, coeff=1, order: int = 1) -> "Hyper":
        return cls({-order: coeff})

    @classmethod
    def const(cls, value) -> "Hyper":
        return cls({0: value})

    @classmethod
    def zero(cls, origin: str | None = None) -> "Hyper":
        return cls({}, origin)

    # --- строение
    @property
    def order(self) -> int | None:
        return max(self.terms) if self.terms else None

    @property
    def leading(self) -> Q | None:
        return self.terms[self.order] if self.terms else None

    @property
    def number_class(self) -> int:
        k = self.order
        if k is None:
            return ZERO
        return INFINITY if k > 0 else (INFINITESIMAL if k < 0 else FINITE)

    @property
    def is_monomial(self) -> bool:
        return len(self.terms) == 1

    def __eq__(self, other) -> bool:
        """Равенство величин; происхождение нуля — запись о пути, в равенство не входит
        (иначе (0ₛ+ω)+(−ω) ≠ 0ₛ+(ω−ω): путь разный, величина одна)."""
        return isinstance(other, Hyper) and self.terms == other.terms

    def same_record(self, other) -> bool:
        return self == other and self.origin == other.origin

    def __hash__(self) -> int:
        return hash(tuple(sorted(self.terms.items())))

    def record(self) -> str:
        """Ф-запись: класс по старшему члену и печать describe()."""
        return f"Φ(C={self.number_class} HYPER, {self.describe()})"

    def _rhs(self, other, op: str, reflected: str | None, hint: str = _NO_RULE_HINT):
        """Второй операнд мостом. reflected — имя отражённого метода операнда,
        когда self слева (`__radd__` …); None — когда self справа и спрашивать
        больше некого. None в ответе — NotImplemented: гражданин с отражённым
        методом решает сам (ноль-состояние, контейнер), чужак — TypeError
        интерпретатора; гражданин без правила — отказ поимённо."""
        h = _lift(other)
        if h is not None:
            return h
        if _citizen(other) and not (reflected and hasattr(type(other), reflected)):
            left, right = (self.record(), _phi(other)) if reflected else (_phi(other), self.record())
            raise Unresolved(op, left, right, hint=hint)
        return None

    # --- арифметика (кольцо многочленов Лорана)
    def __add__(self, other: "Hyper") -> "Hyper":
        other = self._rhs(other, "+", "__radd__")
        if other is None:
            return NotImplemented                          # чужой операнд решает сам (ноль-состояние: __radd__)
        t = dict(self.terms)
        for k, c in other.terms.items():
            t[k] = t.get(k, Q(0)) + c
        r = Hyper(t)
        if not r.terms:
            if self.terms or other.terms:
                k = max(self.terms) if self.terms else max(other.terms)      # старший порядок сокращённого
                r.origin = "ω−ω" if k > 0 else ("ε−ε" if k < 0 else "a−a")   # спецметка нуля: код из конечного списка
            else:
                r.origin = self.origin or other.origin
        return r

    def __neg__(self) -> "Hyper":
        return Hyper({k: -c for k, c in self.terms.items()}, self.origin)

    def __sub__(self, other: "Hyper") -> "Hyper":
        other = self._rhs(other, "-", "__rsub__")
        if other is None:
            return NotImplemented
        return self + (-other)

    def __radd__(self, other) -> "Hyper":
        other = self._rhs(other, "+", None)
        if other is None:
            return NotImplemented
        return other + self                                # кольцо коммутативно: тот же __add__

    def __rsub__(self, other) -> "Hyper":
        other = self._rhs(other, "-", None)
        if other is None:
            return NotImplemented
        return other + (-self)

    def __mul__(self, other: "Hyper") -> "Hyper":
        other = self._rhs(other, "*", "__rmul__")
        if other is None:
            return NotImplemented
        t: dict = {}
        for k1, c1 in self.terms.items():
            for k2, c2 in other.terms.items():
                t[k1 + k2] = t.get(k1 + k2, Q(0)) + c1 * c2
        r = Hyper(t)
        if not r.terms:
            r.origin = self.origin or other.origin or "0·x"
        return r

    def __rmul__(self, other) -> "Hyper":
        other = self._rhs(other, "*", None)
        if other is None:
            return NotImplemented
        return other * self                                # кольцо коммутативно: тот же __mul__

    def inverse(self) -> "Hyper":
        if not self.terms:
            raise Unresolved("1/x", self.record(), hint=_X0_HINT)   # тот же отказ, что OMEGA / 0
        if not self.is_monomial:
            raise ExactRepresentationNotAvailable(
                f"обратный к сумме {self.record()} — бесконечный ряд Лорана, конечно не представим")
        k, c = next(iter(self.terms.items()))
        return Hyper({-k: Q(1) / c})

    def truncate(self, floor: int) -> "Hyper":
        """Обрезка снизу: члены порядка ниже `floor` отброшены (обрезка объявляется)."""
        return Hyper({k: c for k, c in self.terms.items() if k >= floor})

    def self_check(self, deep: bool = False) -> list:
        """Самопроверка величины Лорана независимой дорогой (очередь п. 18).

        `[вычислено здесь]` Проверяется не тем счётом, каким строилось:
        запись нормальна (нулевых коэффициентов нет, показатели целы),
        происхождение живёт только у нуля, класс числа отвечает старшему
        члену, обрезка идемпотентна, хеш согласован с равенством. При
        `deep` — обратный ряд встречной дорогой: `h · laurent_inverse(h, m)`
        обязан давать единицу, и всё расхождение лежать не выше объявленной
        границы усечения `floor = −k − m`; у монома — точную единицу.

        Граница: сумма ряда ниже `floor` здесь не проверяется — её нет в
        конечной записи, там стоит названный остаток (`laurent_inverse`).
        """
        out = []
        for k, c in self.terms.items():
            if not isinstance(k, int):
                out.append(f"показатель не целый: {k!r}")
            if c == 0:
                out.append(f"нулевой коэффициент при ω^{k}: член не член")
        if self.terms and self.origin is not None:
            out.append(f"происхождение «{self.origin}» живёт только у нуля, а запись не нулевая")
        k = self.order
        by_order = ZERO if k is None else (INFINITY if k > 0 else (INFINITESIMAL if k < 0 else FINITE))
        if self.number_class != by_order:
            out.append(f"класс числа {self.number_class} не отвечает старшему члену ω^{k}")
        if self.terms and self.leading != self.terms.get(k):
            out.append("старший коэффициент не тот, что стоит при старшем члене")
        if hash(self) != hash(Hyper(dict(self.terms))):
            out.append("хеш не согласован с равенством")
        if not deep:
            return out
        if self.terms:
            lo = min(self.terms)
            if self.truncate(lo) != self or self.truncate(lo).truncate(lo) != self.truncate(lo):
                out.append(f"обрезка по младшему члену ω^{lo} не тождественна")
            for m in (4, 8):
                try:
                    inv = laurent_inverse(self, m)
                except Exception as e:                      # noqa: BLE001
                    out.append(f"обратный ряд отказал на порядке {m}: {type(e).__name__}")
                    continue
                prod = self * inv["частичная сумма"]
                if prod.terms.get(0) != Q(1):
                    out.append(f"h·h⁻¹ не даёт единицы на порядке {m}: {prod.describe()}")
                floor = -k - m
                high = [j for j in prod.terms if j != 0 and j > floor]
                if high:
                    out.append(f"расхождение h·h⁻¹ выше объявленной границы ω^{floor}: порядки {sorted(high)}")
                if self.is_monomial and prod != Hyper({0: Q(1)}):
                    out.append(f"моном: h·h⁻¹ не точная единица, а {prod.describe()}")
        return out

    def __truediv__(self, other: "Hyper") -> "Hyper":
        # делитель — ноль-состояние: у Zero деления нет, отказ x/0 поимённо
        hint = _X0_HINT if isinstance(other, Zero) else _NO_RULE_HINT
        given = other                                   # отказ называет операнд как дан, не после подъёма
        other = self._rhs(other, "/", "__rtruediv__", hint)
        if other is None:
            return NotImplemented
        if not other.terms:
            raise Unresolved("/", self.record(), _phi(given), hint=_X0_HINT)
        return self * other.inverse()

    def __rtruediv__(self, other) -> "Hyper":
        hint = _X0_HINT if not self.terms else (_ZERO_DIV_HINT if isinstance(other, Zero) else _NO_RULE_HINT)
        given = other
        other = self._rhs(other, "/", None, hint)
        if other is None:
            return NotImplemented
        if not self.terms:
            raise Unresolved("/", _phi(given), self.record(), hint=_X0_HINT)
        return other * self.inverse()

    def __pow__(self, n) -> "Hyper":
        # bool отвергается наравне с прочим не-целым: `+`, `−`, `·`, `/`
        # его уже отвергали (мост `_lift`, `pointers.py`, `Q`), а степень
        # принимала — `Hyper.omega() ** True` давало `1·ω`. Замер 05.09
        # (глава 03, расхождения): дверь одна на пакет, и держать её надо
        # одинаково со всех сторон. Свип 11.09: показатель-гражданин без
        # целого значения (Q(1, 2), указатель, Hyper, ноль-состояние) —
        # отказ поимённо, не TypeError интерпретатора.
        if isinstance(n, (bool, float)):
            raise TypeError("степень — целое (bool не число, float запрещён построением)")
        k = _int_exponent(n)
        if k is None:
            e = _rational_exponent(n)
            if e is not None:
                # дробный показатель есть ВЕЛИЧИНА: решает теорема кольца, не дверь
                return rational_power(self, e, given=n)
            if _citizen(n):
                raise Unresolved("**", self.record(), _phi(n), hint=_POW_HINT)
            raise TypeError("степень — целое либо рациональное (bool не число, "
                            "float запрещён построением)")
        if k < 0:
            if not self.terms:
                raise Unresolved("**", self.record(), _phi(n), hint=_X0_HINT)   # 0^(−k) = 1/0: оба операнда названы
            return self.inverse() ** (-k)
        r = Hyper.const(1)
        for _ in range(k):
            r = r * self
        return r

    def __rpow__(self, base) -> "Hyper":
        """base ** self: целая константа в показателе — степень поднятого
        основания (результат остаётся Hyper); иначе b^ω — отказ поимённо."""
        if isinstance(base, (bool, float)):
            raise TypeError("основание степени — не bool и не float (запрет построением)")
        k = _int_exponent(self)
        if k is None:
            e = _rational_exponent(self)
            b = _lift(base) if e is not None else None
            if b is not None:
                return rational_power(b, e, given=self, given_base=base)  # дверь одна: тот же счёт справа
            if _citizen(base):
                raise Unresolved("**", _phi(base), self.record(), hint=_POW_HINT)
            return NotImplemented
        b = self._rhs(base, "**", None)
        if b is None:
            return NotImplemented
        return b ** k

    # --- представление
    def approx(self) -> float:
        """Значение как float для сверки с классом 23: ω → inf, ε → 0."""
        k = self.order
        if k is None:
            return 0.0
        if k > 0:
            return float("inf") if self.leading > 0 else float("-inf")
        c = self.terms.get(0, Q(0))
        return int(c.n) / int(c.d)

    def describe(self) -> str:
        if not self.terms:
            return "0{" + (self.origin or "1") + "}"                     # спецметка: 0{ω−ω}, 0{ε−ε}, 0{a−a}, 0{1}
        parts = []
        for k in sorted(self.terms, reverse=True):
            c = self.terms[k]
            if k == 0:
                parts.append(f"{c}")
            elif k == 1:
                parts.append(f"{c}·ω")
            elif k == -1:
                parts.append(f"{c}·ε")
            elif k > 0:
                parts.append(f"{c}·ω^{k}")
            else:
                parts.append(f"{c}·ε^{-k}")
        return " + ".join(parts)

    def __repr__(self) -> str:
        return f"Hyper({self.describe()})"


def laurent_inverse(h: "Hyper", order: int = 8) -> dict:
    """Обратный к многочлену Лорана как РЯД, с объявленным порядком усечения.

    `[вычислено здесь]` `1/h` для суммы — бесконечный ряд, конечной записи нет
    (потому `Hyper.inverse` отказывает и отказ остаётся). Здесь считается
    разложение в окрестности `ω → ∞`: `h = c·ω^k·(1 + u)`, где у `u` все порядки
    строго ниже нуля, и `1/h = (1/c)·ω^(−k)·Σ_{j≥0} (−u)^j`. Каждый следующий
    член опускает порядок хотя бы на единицу, потому сумма до `order` членов
    точна во всех порядках до объявленной границы, а ниже — названный остаток.

    Возвращается со знаменателем: частичная сумма, порядок старшего и младшего
    члена, граница остатка и сам порядок усечения — чтобы ни одно число не
    читалось как точное там, где оно частичное.
    """
    if isinstance(h, Zero) or not isinstance(h, Hyper):
        raise Unresolved("1/x", _phi(h), hint=_NO_RULE_HINT)
    if not h.terms:
        raise Unresolved("1/x", h.record(), hint=_X0_HINT)
    if order < 0:
        raise ValueError("порядок усечения не бывает отрицательным")
    k = max(h.terms)
    c = h.terms[k]
    u = Hyper({j - k: Q(v) / c for j, v in h.terms.items() if j != k})
    floor = -k - order                     # ниже этой границы члены не считаются
    total = Hyper({0: Q(1)})
    term = Hyper({0: Q(1)})
    minus_u = Hyper({j: -v for j, v in u.terms.items()})
    for _ in range(order):
        term = (term * minus_u).truncate(floor + k)
        if not term.terms:
            break
        total = total + term
    partial = (Hyper({-k: Q(1) / c}) * total).truncate(floor)
    return {"частичная сумма": partial,
            "порядок усечения": order,
            "старший член": partial.order,
            "младший член": (min(partial.terms) if partial.terms else None),
            "остаток": f"порядка ω^{floor - 1} и ниже",
            "точно": h.is_monomial}


def laurent_div(a: "Hyper", b: "Hyper", order: int = 8) -> dict:
    """`a / b` рядом Лорана с объявленным порядком: `a · laurent_inverse(b)`."""
    inv = laurent_inverse(b, order)
    k = max(b.terms)
    out = dict(inv)
    out["частичная сумма"] = (a * inv["частичная сумма"]).truncate(-k - order + (a.order or 0))
    out["точно"] = b.is_monomial
    return out


# ---------------------------------------------------------------------------
# корень и рациональная степень в кольце ℚ[ω, ω⁻¹]: устройство, которым ЗАКРЫТО
# утверждение 1 пометки NOT_MEASURED (шапка модуля, разбор 17.09); дверь `**`
# зовёт его сама, счёт целыми, ни одного float
# ---------------------------------------------------------------------------

def _rational_root(c: Q, q: int):
    """Точный корень степени q из рационального c, либо None.

    Счёт целыми: числитель и знаменатель по отдельности через `int_root_floor`
    (та же дорога, какой дом извлекает целый корень в сертификатах —
    `phicore/enclosure.py`; второй копии здесь нет намеренно), проверка —
    возведением обратно. Чётный корень из отрицательного вещественного
    значения не имеет, и это названо, а не обойдено.
    """
    n, d = int(c.numerator), int(c.denominator)
    sign = 1
    if n < 0:
        if q % 2 == 0:
            return None                      # чётный корень из отрицательного не веществен
        sign, n = -1, -n
    rn, rd = int_root_floor(n, q), int_root_floor(d, q)
    if rn ** q != n or rd ** q != d:
        return None                          # не точная степень — корень иррационален
    return Q(sign * rn, rd)


def _poly_pow_trunc(g: list, q: int, upto: int) -> list:
    """(Σ g_j·x^j)^q, обрезанный по степени upto — счёт в ℚ, без float."""
    r = [Q(1)] + [Q(0)] * upto
    for _ in range(q):
        nxt = [Q(0)] * (upto + 1)
        for i, a in enumerate(r):
            if a == 0:
                continue
            for j, b in enumerate(g):
                if i + j > upto or b == 0:
                    continue
                nxt[i + j] += a * b
        r = nxt
    return r


def laurent_root(h: "Hyper", q: int):
    """Точный корень степени q ≥ 1 в кольце ℚ[ω, ω⁻¹]: g с g^q == h, либо None.

    `[вычислено здесь]` Корень ЕСТЬ в кольце тогда и только тогда, когда h есть
    q-я степень в нём. Извлечение конечно и точно:

      * h = ω^lo·H(ω), H(0) ≠ 0; у корня g = ω^m·G(ω) обязано быть q·m = lo и
        q·deg G = deg H — иначе корня нет и счёт не начинается;
      * g_0 — точный корень степени q из младшего коэффициента (`_rational_root`);
      * дальше коэффициенты определяются ОДНОЗНАЧНО снизу вверх: у (Σ g)^q
        коэффициент при x^t есть q·g_0^(q−1)·g_t плюс члены, где участвуют лишь
        g_1 … g_{t−1} — потому g_t берётся делением, а не поиском.

    Итог проверяется ВСТРЕЧНОЙ дорогой — возведением g^q арифметикой самого
    кольца (`Hyper.__mul__`), которая рекурренту коэффициентов не знает: где
    дороги разошлись, корня нет и возвращается None, а не «почти корень».
    """
    if not isinstance(h, Hyper):
        raise TypeError("laurent_root — только Hyper (величина кольца ℚ[ω, ω⁻¹])")
    if isinstance(q, bool) or not isinstance(q, int):
        raise TypeError("степень корня — только int (bool не число)")
    if q < 1:
        raise ValueError("степень корня — целое ≥ 1")
    if q == 1:
        return Hyper(dict(h.terms), h.origin)
    if not h.terms:
        return Hyper.zero(h.origin)          # 0 = 0^q при любом q ≥ 1
    lo, hi = min(h.terms), max(h.terms)
    span = hi - lo
    if lo % q or span % q:
        return None                          # порядок или протяжённость не делятся — корня в кольце нет
    m, d = lo // q, span // q
    H = [h.terms.get(lo + j, Q(0)) for j in range(span + 1)]
    g0 = _rational_root(H[0], q)
    if g0 is None:
        return None                          # младший коэффициент не q-я степень в ℚ
    G = [g0] + [Q(0)] * d
    denom = q * g0 ** (q - 1)
    for t in range(1, d + 1):
        G[t] = Q(0)
        part = _poly_pow_trunc(G, q, t)[t]   # вклад уже известных g_1 … g_{t−1}
        G[t] = (H[t] - part) / denom
    g = Hyper({m + j: G[j] for j in range(d + 1)})
    return g if g ** q == h else None        # встречная дорога: арифметика кольца


def rational_power(h: "Hyper", e, given=None, given_base=None) -> "Hyper":
    """h^(p/q) ТОЧНО в ℚ[ω, ω⁻¹], где величина там есть; иначе отказ поимённо.

    `[вычислено здесь]` Граница названа теоремой, а не осторожностью. Пусть
    e = p/q несократима, q ≥ 1. Тогда h^e лежит в кольце ⟺ h есть q-я степень
    в нём. Необходимость: если f = h^(p/q) ∈ R, то f^q = h^p; взяв a, b с
    a·p + b·q = 1, получаем h = (f^a·h^b)^q, то есть h — q-я степень (мономы
    кольца обратимы, потому отрицательное b законно). Достаточность очевидна:
    h = g^q ⇒ h^(p/q) = g^p.

    Потому отказ здесь есть отказ О ПРЕДМЕТЕ, и он отличает себя от двух
    соседних утверждений: «не умею посчитать» — неверно (извлечение конечно и
    точно, `laurent_root`), «величины не существует» — тоже неверно (ω^(1/2)
    живёт в поле Леви-Чивита, √5 — в ℝ). Верно третье: величина есть, а
    конечной записи ЭТОГО кольца у неё нет.

    Граница фасада `[вычислено здесь]` (18.09, очередь п. 88): `laurent_root` и
    `rational_power` в `phicore.__all__` не выведены — 0 из 2 при 2 из 2 у
    `laurent_inverse` / `laurent_div`. Различие не в чине, а в двери: у ряда
    Лорана оператор `/` на сумму ОТКАЗЫВАЕТ (и обязан), потому явный запрос —
    единственная дорога, и она стоит в фасаде; сюда же дорогу берёт сама дверь
    `**` с обеих сторон (приёмка: дверь и явный вызов дают одно на 96 из 96
    троек), и второй дверью того же действия фасад не удваивается. Явный вызов
    — `from phicore.hyper import rational_power, laurent_root`.
    """
    if not isinstance(h, Hyper):
        raise TypeError("rational_power — только Hyper")
    e = _q(e)
    given = e if given is None else given
    p, q = int(e.numerator), int(e.denominator)
    if q == 1:
        return h ** p
    g = laurent_root(h, q)
    if g is None:
        # операнды называются КАК ДАНЫ, не после подъёма (правило моста 11.09):
        # 5 ** Hyper(1/2) обязан сказать Φ(C=1 RATIONAL, value=5), не Φ(C=1 HYPER, 5)
        left = h.record() if given_base is None else _phi(given_base)
        raise Unresolved("**", left, _phi(given), hint=_ROOT_HINT.format(q=q))
    if p < 0 and not g.terms:
        # 0^(−p/q) = 1/0: тот же отказ x/0, и операнды названы КАК ДАНЫ,
        # а не как целая степень поднятого корня (иначе −1/2 печаталось бы как −1)
        left = h.record() if given_base is None else _phi(given_base)
        raise Unresolved("**", left, _phi(given), hint=_X0_HINT)
    return g ** p
