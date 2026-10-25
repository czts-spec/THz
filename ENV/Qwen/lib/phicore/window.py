"""Окно представления: двоичный формат IEEE как окно контейнера v6.

Слово автора 05.09 `[в корпусе]`: «меняется не само число IEEE, а его
представление, интерпретация; весь матаппарат Ф-чисел к этому
представлению был применим, однако не применён»; «одометр не помнит, как
быстро или медленно накапливался полный путь — одометр показывает точный
результат». Здесь аппарат применён к представлению; ни одного float внутри.

Численная форма. Число — элемент контейнера Φ = (m1/m2)·2^(e1/e2); при
m2 = e2 = 1 это засечка k·2^e целой линейки ℤ[1/2] (кольцо: сложение,
вычитание, умножение замкнуты; деление выводит в ℚ — слот m2). Формат
(p, emin, emax) — не другое число, а ОКНО на эту линейку: показаны
мантиссы |k| < 2^p на колёсах 2^(e−p+1), emin ≤ e ≤ emax, и субнормальные
засечки младшего колеса 2^(emin−p+1). Образ числа в окне — тот же элемент
контейнера (класс 1, m2 = e2 = 1); за верхним краем — указатель ω
(класс 22), состояния (не вычислено) — zero_state.

Геометрия (ОП). Окно — одометр: мантисса — лимб на 2^p засечек, показатель —
счётчик оборотов (поле порядка); перенос между колёсами точен («засечка,
теряя с одной стороны, приобретает с другой»). Что мельче младшего колеса или
за верхним краем — не потеряно, а не показано: оно записано в
base_params['absorbed'] тем же канцеллером указателя, что поглощает
младший масштаб у ω и ε (container._absorb). Сертификат — тождество:
образ + поглощённое = число, всегда, по построению.

Алгебра. Действие формата: a ∘ b := project(exact(a) ∘ exact(b)) —
одометр: состояние после любого пути равно точному элементу кольца (для
÷ — поля ℚ); образ зависит от пути только потому, что окно его не показывает.
Тот же формат с забыванием (remember=False): действие над образами,
поглощённое не переносится — это и есть арифметика IEEE, интерпретация
окна, а не свойство числа. Отказ — Unresolved с Ф-записями, не ложный класс.

Указатели и точность (слово автора 05.09, ~16:2x). Засечки одной линейки
равны по определению линейки: внутри неё нет меры мельче засечки, чтобы их
различить; равномерность в идеальном смысле не гарантирована и не
утверждается. Указатели ε (класс 21) и ω (класс 22) — не числа, а смена
линейки: колесо мельче или крупнее; потому они не EXACT — точность на новой
линейке наследуется от самой грубой линейки расчёта (главной базы).
Численная запись истинной точности — долями базы: младшее колесо / база =
2^−(p·n_in) при n_in колёсах внутрь (`base_fraction`), размах — единицами
крупнейшего 2^(p·n_out) (`reach`); как ppm кварца при базе герц/секунда.
"""
from __future__ import annotations

from dataclasses import replace

from .container import NumberClass, PhiNumberV6, _absorb
from .errors import Unresolved
from .rational import Q
from .wheels import top_wheel2
from .zero_state import Zero

#: Допустимые скаляры. float исключён построением.
_NUM = (int, Q)


def _as_int(x, name: str) -> int:
    if isinstance(x, bool) or not isinstance(x, int):
        raise TypeError(f"{name} — только int")
    return x


def _rec(x) -> str:
    r = getattr(x, "record", None)
    if callable(r):
        return r()
    if isinstance(x, PhiNumberV6):
        return repr(x)
    if isinstance(x, _NUM):
        return f"Φ(C=1 RATIONAL, value={Q(x)})"
    return f"Φ(?={x!r})"


def _pow2(e: int) -> Q:
    return Q(1 << e) if e >= 0 else Q(1, 1 << (-e))


def _floor_log2(q: Q) -> int:
    """⌊log₂|q|⌋ целыми: q ≠ 0 — он же СТАРШИЙ НОМЕР КОЛЕСА сетки (2, 1).

    Замер 22.09 `[вычислено здесь]`: здесь стояла вторая рукописная копия того
    же закона, что и в `ieee._binade`; копии совпадали на 671 паре, и ничто их
    не сверяло. Тропа теперь ОДНА, в модуле закона (`wheels.top_wheel2`), её
    равносильность колесу проверяется переборно приёмкой. Дверь горячая —
    359 473 вызова за набор, — потому именно тропа, а не общая дорога колеса.
    """
    return top_wheel2(abs(q.numerator), q.denominator)


class Window:
    """Окно (p, emin, emax): p — засечек лимба (бит мантиссы), emin ≤ e ≤ emax —
    обороты (показатель величины, IEEE-соглашение: 2^(p−1) ≤ |k| < 2^p при
    значении k·2^(e−p+1)). binary64 — Window(53, −1022, 1023); binary32 —
    Window(24, −126, 127); приёмка ведётся на малых p (стратегия снижения
    разрядности, слово автора 05.09): окно перечислимо целиком."""

    __slots__ = ("p", "emin", "emax")

    def __init__(self, p: int, emin: int, emax: int):
        p, emin, emax = _as_int(p, "p"), _as_int(emin, "emin"), _as_int(emax, "emax")
        if p < 1:
            raise ValueError("p ≥ 1")
        if emin > emax:
            raise ValueError("emin ≤ emax")
        self.p, self.emin, self.emax = p, emin, emax

    # --- геометрия окна ---------------------------------------------------

    @property
    def tick(self) -> Q:
        """Младшее колесо: единица субнормальных засечек 2^(emin−p+1)."""
        return _pow2(self.emin - self.p + 1)

    @property
    def top(self) -> Q:
        """Старшая засечка окна (2^p − 1)·2^(emax−p+1)."""
        return Q((1 << self.p) - 1) * _pow2(self.emax - self.p + 1)

    def base_fraction(self, scales_in: int = 0) -> Q:
        """Истинная точность долями базы: единица колеса на scales_in масштабов
        внутрь / единица базового колеса = 2^−(p·(scales_in + 1))."""
        scales_in = _as_int(scales_in, "scales_in")
        if scales_in < 0:
            raise ValueError("scales_in ≥ 0")
        return Q(1, 1 << (self.p * (scales_in + 1)))

    def reach(self, scales_out: int = 0) -> int:
        """Размах единицами крупнейшего колеса: 2^(p·scales_out) базовых единиц."""
        scales_out = _as_int(scales_out, "scales_out")
        if scales_out < 0:
            raise ValueError("scales_out ≥ 0")
        return 1 << (self.p * scales_out)

    def count(self) -> int:
        """Число конечных элементов окна: 0 и ± по оборотам и субнормальным."""
        normals = (1 << (self.p - 1)) * (self.emax - self.emin + 1)
        sub = (1 << (self.p - 1)) - 1
        return 1 + 2 * (normals + sub)

    def elements(self):
        """Все конечные элементы окна (Q), для исчерпывающей приёмки при малом p."""
        yield Q(0)
        half = 1 << (self.p - 1)
        for m in range(1, half):                       # субнормальные
            v = Q(m) * self.tick
            yield v
            yield -v
        for e in range(self.emin, self.emax + 1):
            unit = _pow2(e - self.p + 1)
            for m in range(half, 1 << self.p):
                v = Q(m) * unit
                yield v
                yield -v

    def unit_at(self, q: Q) -> Q:
        """Колесо, на которое ложится величина |q|: 2^(e−p+1), не мельче tick."""
        if q == 0:
            return self.tick
        e = _floor_log2(q)
        return _pow2(max(e, self.emin) - self.p + 1)

    def contains(self, q) -> bool:
        """Есть ли засечка окна ровно в q."""
        q = Q(q)
        if abs(q) > self.top:
            return False
        r = q / self.unit_at(q)
        return r.denominator == 1

    def round(self, q):
        """Ближайшая засечка окна (при равенстве — чётная); None — за верхним краем."""
        q = Q(q)
        unit = self.unit_at(q)
        r = q / unit
        k, rem = divmod(r.numerator, r.denominator)
        if 2 * rem > r.denominator or (2 * rem == r.denominator and k & 1):
            k += 1
        v = Q(k) * unit
        if abs(v) > self.top:
            return None
        return v

    # --- элементы контейнера ---------------------------------------------

    @staticmethod
    def _element(q: Q) -> PhiNumberV6:
        """Засечка q ∈ ℤ[1/2] как элемент контейнера класса 1 (m2 = e2 = 1)."""
        if q == 0:
            return PhiNumberV6.from_int(0)
        n, d = q.numerator, q.denominator
        if d & (d - 1):
            raise ValueError("засечка окна — двоичная дробь")
        e = -(d.bit_length() - 1)
        while n % 2 == 0:
            n //= 2
            e += 1
        return PhiNumberV6.from_mantissa_exponent(n, 1, e, 1)

    @staticmethod
    def _absorbed_q(x: PhiNumberV6):
        """Поглощённое как Q; Q(0) — записи нет; None — не сведено."""
        if "absorbed" not in x.base_params:
            return Q(0)
        h = x.base_params["absorbed"]
        if h is None:
            return None
        if any(k != 0 for k in h.terms):
            return None
        return h.terms.get(0, Q(0))

    def exact(self, x):
        """Показание одометра целиком: образ + поглощённое (Q); None — состояние.

        Не `ieee.exact` (точная сторона одного float: пара, ±ω, состояние) и
        не `affine.Affine.exact` (форма без шума из узла сетки): здесь вход —
        элемент контейнера, int или Q, а ответ — величина ДО поглощения окном.
        """
        if isinstance(x, bool):
            raise TypeError("bool не число")
        if isinstance(x, _NUM):
            return Q(x)
        if isinstance(x, Zero):
            return None
        if not isinstance(x, PhiNumberV6):
            raise TypeError("exact — элемент контейнера, int или Q")
        if x.is_zero_class or x.is_uncomputable or x.is_ignorance:
            return None
        a = self._absorbed_q(x)
        if a is None:
            return None
        if x.is_infinity or x.is_infinitesimal:
            return a                                   # указатель — край окна, величина в записи
        n, d = x.to_exact_int_ratio()
        return Q(n, d) + a

    @staticmethod
    def image(x):
        """Что показывает окно: элемент без записи поглощённого (взгляд IEEE)."""
        if isinstance(x, PhiNumberV6) and "absorbed" in x.base_params:
            return replace(x, base_params={k: v for k, v in x.base_params.items() if k != "absorbed"})
        return x

    def project(self, x):
        """Проекция на окно: образ с записью поглощённого; образ + поглощённое = x."""
        if isinstance(x, Zero):
            return x
        if isinstance(x, PhiNumberV6) and (x.is_zero_class or x.is_uncomputable or x.is_ignorance):
            return x
        q = self.exact(x)
        if q is None:
            raise Unresolved("project", _rec(x), hint="полная величина не сведена — окну нечего показывать")
        v = self.round(q)
        if v is None:
            img = PhiNumberV6.omega(1 if q > 0 else -1, 1, 1)
            return _absorb(img, PhiNumberV6.from_ratio(q.numerator, q.denominator))
        img = self._element(v)
        r = q - v
        if r == 0:
            return img
        return _absorb(img, PhiNumberV6.from_ratio(r.numerator, r.denominator))

    # --- действие формата ------------------------------------------------

    def op(self, a, b, op: str, remember: bool = True):
        """a ∘ b в окне. remember=True — одометр: действие над полными
        величинами, проекция; remember=False — IEEE: действие над образами
        арифметикой контейнера (указатели по его правилам), проекция."""
        if op not in "+-*/" or len(op) != 1:
            raise ValueError("op ∈ {+, -, *, /}")
        if not remember:
            ia, ib = self.image(a), self.image(b)
            if isinstance(ia, Zero) or isinstance(ib, Zero):
                return Zero.not_computed()             # состояние распространяется (nan у IEEE)
            if isinstance(ia, _NUM):
                ia = self.project(ia)
            if isinstance(ib, _NUM):
                ib = self.project(ib)
            if op == "+":
                r = ia + ib
            elif op == "-":
                r = ia - ib
            elif op == "*":
                r = ia * ib
            else:
                r = ia / ib
            if isinstance(r, PhiNumberV6) and r.number_class == NumberClass.NATURAL_POWER and r.base_class is None:
                return self.project(self.image(r))
            return r                                   # указатель или состояние — по правилам контейнера
        qa, qb = self.exact(a), self.exact(b)
        if qa is None or qb is None:
            return Zero.not_computed()
        if op == "+":
            q = qa + qb
        elif op == "-":
            q = qa - qb
        elif op == "*":
            q = qa * qb
        else:
            if qb == 0:
                return Zero.not_computed()
            q = qa / qb
        return self.project(q)

    def self_check(self, deep: bool = False) -> list:
        """Самопроверка окна ВСТРЕЧНОЙ дорогой (ось самоприменения).

        `[вычислено здесь]` Дороги: (1) ТОЖДЕСТВО ОДОМЕТРА, на котором стоит
        весь модуль (слово автора 05.09): образ плюс поглощённое есть само
        число — `exact(project(q)) == q` для проб, причём `project` округляет
        и записывает остаток, а `exact` собирает его обратно другой дорогой
        (чтение `base_params['absorbed']`), так что сверяются две разные
        стороны одного действия; (2) число засечек ДВУМЯ дорогами — формулой
        `count()` (нормальные плюс субнормальные, оба знака и ноль) и прямым
        ПЕРЕЧИСЛЕНИЕМ `elements()`; там же младшая и старшая засечки
        сверяются с `tick` и `top`; (3) всякая засечка перечисления обязана
        опознаваться `contains`, а всякое округление — попадать в окно и
        отстоять от величины не дальше половины своего колеса.

        Граница объявлена числом: перечисление идёт ТОЛЬКО при малом окне
        (p ≤ 4 и не более 4096 засечек) — у binary64 их 2^64, и дорога (2)
        там не идёт вовсе; её молчание находкой не считается, а тождество
        одометра (1) проверяется на любом окне пробами. Возвращает список
        находок; исключений не бросает.
        """
        out = []
        for name in ("p", "emin", "emax"):
            if not isinstance(getattr(self, name), int):
                out.append(f"поле {name} не целое: {getattr(self, name)!r}")
        if out:
            return out
        if self.p < 1:
            out.append(f"p = {self.p} < 1")
        if self.emin > self.emax:
            out.append(f"emin {self.emin} > emax {self.emax}")
        if out:
            return out
        probes = [Q(0), self.tick, -self.tick, self.tick / Q(3), self.top,
                  self.top * Q(2, 3), -self.top / Q(7), Q(1, 3), Q(5, 8)]
        for q in probes:
            try:
                v = self.project(q)
            except Unresolved:
                continue                    # за верхним краем — край окна, объявлено
            except Exception as e:                              # noqa: BLE001
                out.append(f"проекция отказала на {q}: {type(e).__name__}: {e}")
                continue
            back = self.exact(v)
            if back is None:
                out.append(f"одометр не собрался обратно на {q}: величина "
                           f"не сведена")
            elif back != q:
                out.append(f"тождество одометра нарушено: образ + поглощённое "
                           f"= {back}, а было {q}")
            img = self.image(v)
            if isinstance(img, PhiNumberV6) and "absorbed" in img.base_params:
                out.append(f"образ окна несёт запись поглощённого: {img.record()}")
            r = self.round(q)
            if r is not None:
                if not self.contains(r):
                    out.append(f"округление дало величину вне окна: {r}")
                if abs(q - r) * 2 > self.unit_at(q):
                    out.append(f"округление дальше половины колеса: {q} → {r}")
        if self.p <= 4 and self.count() <= 4096:
            seen = list(self.elements())
            if len(seen) != self.count():
                out.append(f"число засечек двумя дорогами расходится: формула "
                           f"{self.count()}, перечисление {len(seen)}")
            if len(set(seen)) != len(seen):
                out.append(f"перечисление повторяется: {len(seen)} штук, "
                           f"{len(set(seen))} различных")
            pos = [v for v in seen if v > 0]
            if pos:
                if min(pos) != self.tick:
                    out.append(f"младшая засечка {min(pos)} не равна tick "
                               f"{self.tick}")
                if max(pos) != self.top:
                    out.append(f"старшая засечка {max(pos)} не равна top "
                               f"{self.top}")
            for v in seen:
                if not self.contains(v):
                    out.append(f"своя засечка не опознана: {v}")
                    break
        if not deep:
            return out
        a_q, b_q = self.top / Q(3), self.tick * Q(5)
        for op in ("+", "-", "*"):
            r = self.op(a_q, b_q, op, remember=True)
            want = {"+": a_q + b_q, "-": a_q - b_q, "*": a_q * b_q}[op]
            got = self.exact(r)
            if got is None:
                if abs(want) <= self.top:
                    out.append(f"одометр не сведён на {op}, хотя величина в окне")
            elif got != want:
                out.append(f"одометр на {op} дал {got}, точная величина {want}")
        return out

    def record(self) -> str:
        return f"Φ-окно(p={self.p}, emin={self.emin}, emax={self.emax}, засечек={self.count()})"

    __repr__ = record
