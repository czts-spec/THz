"""Произвольное основание α — член семейства контейнеров.

Находка автора (письмо 02.09): «"α=2" нахожу, что произвольное
основание экспоненты фиктивно… как туда положить число класса 1–19…
это фундамент контейнера». Здесь основание положено носителем:

    Φ = m · α^e,   m и e — пары (rational.Q), α — носитель классов 1–19
                   (пара, линейка 3, моном π, экспонента, …) или пара.

Контейнер v6 — эталон семантики, не размещения (слово автора 02.09):
формула та же, но ось α открыта. Фиктивность основания при этом
ИЗМЕРИМА: одно значение в записях с α=2, 4, 16, 1/2 равно себе —
основание есть свойство ЗАПИСИ, не значения.

Граничные условия (первые принципы, каноном построения):

* m=0 — нулевая пара: контейнер несёт точный ноль (α и e гаснут);
* α=1 — ТРИВИАЛЬНОЕ СОСТОЯНИЕ (автор: «1/Ф=Ф/1 допускает α=Ф=1/Ф=1»):
  1^e=1 при любом e, ось экспоненты гаснет — e:=0; при m=±1 контейнер
  равен собственной инверсии (Ф=1/Ф);
* e=0 — основание фиктивно целиком: α:=1 (значение равно m);
* α=0 — нулевая ось: e>0 даёт точный ноль, e=0 — отказ (0⁰ не
  сведено), e<0 — ZeroDivisionError;
* α<0 — знак выносится в мантиссу чётностью числителя e (нечётный
  знаменатель), чётный знаменатель — не веществен, отказ построением;
* 0<α<1 — двойственность 1/Ф=Ф/1: канон α:=1/α, e:=−e (α>1 всегда);
* рациональное α сводится к ПЕРВООБРАЗНОМУ основанию: α=β^k → α:=β,
  e:=k·e (8^(1/3) и 2^1 — одна запись); соизмеримые базы после канона
  совпадают буквально.

Кросс-база в контейнере — честный отказ CrossBaseOperationNotSupported;
сведение — переводом `as_carrier()` (ФIOS зовёт его мостом сам:
контейнерная запись переводится в линейку, не наоборот).

Логарифм сложения (Гаусс 1812, Цех; слово автора 02.09 «известное —
аппаратом»): сумма над ОДНОЙ базой при разных показателях остаётся одной
записью над той же базой — показатель несёт величину класса 6:
m₁α^e₁ + m₂α^e₂ = m₁·α^(e₁ + log_α(1 + (m₂/m₁)·α^(e₂−e₁))). Носитель
показателя — `zech.ZechExponent`; умножение на оси E — по-прежнему
сложение. Над рациональной базой запись переводится в носителя точно
(пара / линейка 3); над основанием-носителем (e, π) запись точна как
есть, а перевода нет — носителя для суммы над трансцендентной базой в
библиотеке нет, и это названо, а не спрятано.

Нормальная форма над трансцендентной базой (план 03.09, шаг 1): при
β = α^(1/N) запись есть рациональная функция от β; β трансцендентна
(Эрмит 1873, Линдеман 1882), потому ℚ(β) ≅ ℚ(x) — равенство чисел есть
равенство несократимых дробей `ratfunc.RatFunc`, а не структуры записи.
Канон: если форма — моном c·x^n, запись сворачивается к m·α^(n/N)
(история сложений гаснет); ноль распознаётся точно (P = 0); знак
1 + k·α^d — интервальным Горнером по сертификату β (корень из отрезка,
шаг 2: дробная степень π сертифицирована).
"""
from __future__ import annotations

from math import gcd

from .container import (CrossBaseOperationNotSupported, NumberClass,
                        PhiNumberV6)
from .errors import ExactRepresentationNotAvailable, Unresolved
from .rational import Q
from .zech import (Zech, ZechExponent, exp_add, exp_neg, exp_rebase,
                   exp_scale)


def _prime_powers(n: int) -> dict:
    """Показатели простых; пробное деление до 10⁴, остаток атомом
    (предел — разрядность, не величина: правило 02.09)."""
    fac = {}
    d, m = 2, n
    while d * d <= m and d <= 10 ** 4:
        while m % d == 0:
            fac[d] = fac.get(d, 0) + 1
            m //= d
        d += 1 if d == 2 else 2
    if m > 1:
        fac[m] = fac.get(m, 0) + 1
    return fac


def _int_root_exact(t: int, k: int):
    """Целый корень k-й степени точно или None."""
    from .algebraic import _int_root_floor
    r = _int_root_floor(t, k)
    return r if r ** k == t else None


_PRIM_CACHE: dict = {}      # (n, d) → (β, k): факторизация одной базы — один раз
_PRIM_CACHE_CAP = 4096


def _primitive_base(a: Q) -> tuple[Q, int]:
    """(β, k): a = β^k с наибольшим целым k; β — первообразное основание.

    a > 1 каноном; best-effort на атомах факторизации — несведённый
    атом даёт k=1, запись остаётся верной (β=a). Итог запоминается по
    паре (n, d): канон зовётся при каждом построении записи."""
    key = (a.n, a.d)
    hit = _PRIM_CACHE.get(key)
    if hit is not None:
        return hit
    r = _primitive_base_uncached(a)
    if len(_PRIM_CACHE) >= _PRIM_CACHE_CAP:
        _PRIM_CACHE.clear()
    _PRIM_CACHE[key] = r
    return r


def _primitive_base_uncached(a: Q) -> tuple[Q, int]:
    n, d = a.numerator, a.denominator
    exps = list(_prime_powers(n).values())
    if d > 1:
        exps += list(_prime_powers(d).values())
    k0 = 0
    for e in exps:
        k0 = gcd(k0, e)
    for k in range(k0, 1, -1):
        if k0 % k:
            continue
        rn = _int_root_exact(n, k)
        rd = _int_root_exact(d, k) if d > 1 else 1
        if rn is not None and rd is not None:
            return Q(rn, rd), k
    return a, 1




# ---------- радикальное поле ℚ(β), β = α^(1/N): суммы над рациональной базой ----------

def _exp_lcm(e, acc: int) -> int:
    """НОК знаменателей всех показателей записи (рекурсивно по членам)."""
    from math import lcm
    if isinstance(e, ZechExponent):
        acc = lcm(acc, e.q.denominator)
        for _, z in e.terms:
            acc = _exp_lcm(z.d, acc)
        return acc
    return lcm(acc, Q(e).denominator)


_CARRIERS: tuple = ()      # (ExpValue, PiMonomial, LogValue) — один импорт на процесс


def _carriers() -> tuple:
    global _CARRIERS
    if not _CARRIERS:
        from .exp_log import ExpValue, LogValue
        from .pi_value import PiMonomial
        _CARRIERS = (ExpValue, PiMonomial, LogValue)
    return _CARRIERS


def _is_transcendental(a) -> bool:
    """Базы, чья трансцендентность — теорема: coeff·e^y (y∈ℚ, y≠0; Эрмит–
    Линдеман), q·π^e (e∈ℚ, e≠0, q алгебраическое; Линдеман) и c·ln(a)
    (a∈ℚ₊, a≠1; следствие Линдемана–Вейерштрасса: e^{ln a} = a). Прочие
    носители сюда не записываются: недоказанное не есть аппарат."""
    return type(a) in _carriers()


# ---------- соизмеримость баз-носителей: модуль Котса (план 03.09, шаг 3) ----------

def _carrier_parts(a):
    """(c, y, rebuild) для носителя вида c·X^y: X = e (класс 5), X = π
    (класс 10); для c·ln(a) (класс 6) y = None — только скалярная ось.
    None — носитель другого устройства (алгебраический коэффициент и т.п.)."""
    ExpValue, PiMonomial, LogValue = _carriers()
    t = type(a)
    if t is ExpValue:
        return a.coeff, a.y, ExpValue.make
    if t is PiMonomial and type(a.coeff) is Q:
        return a.coeff, a.e, PiMonomial._make
    if t is LogValue:
        return a.coeff, None, (lambda c, _y: LogValue.make(c, a.arg))
    return None


def _exponent_is_integral(e) -> bool:
    """Все показатели записи целые (q и каждое d — пары со знаменателем 1,
    без вложений): тогда c^E и c^d рациональны."""
    if isinstance(e, ZechExponent):
        return (e.q.d == 1 and all(isinstance(z.d, Q) and z.d.d == 1
                                   for _, z in e.terms))
    return Q(e).d == 1


def _scale_exponent(e, c: Q):
    """(c·β)^E = c^q·β^(q + Σ cᵢ·log_β(1 + kᵢ·c^{dᵢ}·β^{dᵢ})) при целых
    показателях: возвращает (множитель в мантиссу, показатель над β)."""
    if isinstance(e, ZechExponent):
        terms = tuple((ci, Zech(z.k * c ** int(z.d), z.d)) for ci, z in e.terms)
        return c ** int(e.q), ZechExponent(e.q, terms)
    return c ** int(e), e


def _canon_carrier(m: Q, a, e):
    """Канон основания-носителя (Котс 1714: модуль системы; канон
    первообразного основания, как для пары):

    1. степенная ось: c·X^y = (c^(1/k)·X^(y/k))^k с наибольшим k (знак k —
       по знаку y, чтобы у β показатель X был положителен: 1/Ф=Ф/1);
       при c=1 первообразное — сам X: e^y ≡ (e^1)^y, π^e ≡ (π^1)^e;
    2. скалярная ось при ЦЕЛЫХ показателях: (c·β)^E → множитель c^E в
       мантиссу, основание — β с коэффициентом 1 (ln 8 = 3·ln 2 →
       Φα(1, ln 8, 1) ≡ Φα(3, ln 2, 1)); при дробном показателе основание
       остаётся своим (√2·√e над базой 2e — мантисса пары его не несёт),
       а отрицательный коэффициент — отказ построением.
    Возвращает (m, α, e)."""
    for _ in range(4):                        # оси чередуются до неподвижности
        parts = _carrier_parts(a)
        if parts is None:
            break
        c, y, rebuild = parts
        if c.n == 1 and c.d == 1 and (y is None or (y.n == 1 and y.d == 1)):
            break                             # уже первообразное: X^1, коэффициент 1
        if y is not None:
            if c == 1:
                k = Q(y)
            elif c > 0:
                base, p = _primitive_base(c)
                k = Q(p if y > 0 else -p)
                c = base if k > 0 else 1 / base
            else:
                k = Q(1)
            if k != 1:
                a = rebuild(c, Q(y) / k)
                e = exp_rebase(e, k)
                continue
        if c != 1:
            if _exponent_is_integral(e):
                factor, e = _scale_exponent(e, c)
                m = m * factor
                a = rebuild(Q(1), y)
                continue
            if c < 0:
                raise ExactRepresentationNotAvailable(
                    f"знак основания-носителя {a.record()} при дробном "
                    f"показателе {e} не отделён")
        break
    return m, a, e
    c, y, rebuild = parts
    if y is not None:
        if c == 1:
            k = y
        elif c > 0:
            base, p = _primitive_base(c)
            p = p if y > 0 else -p
            k = Q(p)
            c = base if p > 0 else 1 / base
        else:
            k = Q(1)
        if k != 1:
            a = rebuild(c, y / k)
            e = exp_rebase(e, k)
            c = a.coeff if hasattr(a, "coeff") else Q(1)
    if c != 1:
        if _exponent_is_integral(e):
            factor, e = _scale_exponent(e, c)
            m = m * factor
            a = rebuild(Q(1), y)
        elif c < 0:
            raise ExactRepresentationNotAvailable(
                f"знак основания-носителя {a.record()} при дробном "
                f"показателе {e} не отделён")
    return m, a, e


def _beta_enclosure(a, N: int, terms: int):
    """Сертификат сетки β = α^(1/N) — корень из отрезка носителя."""
    from .enclosure import root_enclosure
    iv = a.enclosure(terms)
    return root_enclosure(iv, N, bits=max(64, 5 * terms)) if N > 1 else iv


class _Radical:
    """Элементы ℚ[x]/(x^N − α): списки из N пар, x = β = α^(1/N)."""

    __slots__ = ("N", "alpha")

    def __init__(self, N: int, alpha: Q):
        self.N, self.alpha = N, alpha

    def one(self):
        return [Q(1)] + [Q(0)] * (self.N - 1)

    def mono(self, n: int):
        """β^n = α^(n div N)·β^(n mod N)."""
        q, r = divmod(n, self.N)
        out = [Q(0)] * self.N
        out[r] = self.alpha ** q
        return out

    def add(self, a, b):
        return [x + y for x, y in zip(a, b)]

    def scale(self, a, k: Q):
        return [x * k for x in a]

    def mul(self, a, b):
        N = self.N
        prod = [Q(0)] * (2 * N - 1)
        for i, x in enumerate(a):
            if x == 0:
                continue
            for j, y in enumerate(b):
                if y != 0:
                    prod[i + j] += x * y
        out = prod[:N]
        for k in range(N, 2 * N - 1):
            out[k - N] += prod[k] * self.alpha      # β^N = α
        return out

    def inv(self, a):
        """Обратный расширенным Евклидом в ℚ[x] по модулю x^N − α."""
        def strip(p):
            p = list(p)
            while p and p[-1] == 0:
                p.pop()
            return p

        def divmod_(u, v):
            u, v = strip(u), strip(v)
            q = [Q(0)] * max(1, len(u) - len(v) + 1)
            while len(u) >= len(v) and u:
                c = u[-1] / v[-1]
                d = len(u) - len(v)
                q[d] = c
                for i, y in enumerate(v):
                    u[i + d] -= c * y
                u = strip(u)
            return q, u

        m = [Q(0)] * (self.N + 1)
        m[0], m[self.N] = -self.alpha, Q(1)
        r0, r1 = m, strip(a)
        s0, s1 = [Q(0)], [Q(1)]
        while r1:
            q, r = divmod_(r0, r1)
            # s = s0 − q·s1
            qs = [Q(0)] * (len(q) + len(s1) - 1)
            for i, x in enumerate(q):
                for j, y in enumerate(s1):
                    qs[i + j] += x * y
            s = [(s0[i] if i < len(s0) else Q(0)) - (qs[i] if i < len(qs) else Q(0))
                 for i in range(max(len(s0), len(qs)))]
            r0, r1, s0, s1 = r1, r, s1, s
        r0 = strip(r0)
        if len(r0) != 1:
            raise ExactRepresentationNotAvailable(
                "радикальное поле: элемент не обратим по модулю x^N − α")
        inv = [x / r0[0] for x in s0]
        inv = (inv + [Q(0)] * self.N)[: self.N]
        return inv

    def pow(self, a, c: Q):
        if c.d != 1:
            raise ExactRepresentationNotAvailable(
                "радикальное поле: дробная степень множителя не сведена")
        n = c.n
        base = self.inv(a) if n < 0 else a
        n = abs(n)
        out, acc = self.one(), base
        while n:
            if n & 1:
                out = self.mul(out, acc)
            n >>= 1
            if n:
                acc = self.mul(acc, acc)
        return out

    def of_exponent(self, e):
        """α^E как элемент поля."""
        if isinstance(e, ZechExponent):
            out = self.mono(int(e.q * self.N))
            for c, z in e.terms:
                x = self.add(self.one(), self.scale(self.of_exponent(z.d), z.k))
                out = self.mul(out, self.pow(x, c))
            return out
        return self.mono(int(Q(e) * self.N))


#: Классы, у которых значение не есть величина на линейке: пустота, познанный
#: предел, бесконечно малое и большое, измеренное невежество. Их семантику
#: несёт контейнер v6 — эталон (слово автора 02.09: «семантика, не
#: размещение»), и Φα её НЕ ПОВТОРЯЕТ, а делегирует: значение принадлежит
#: классу своего носителя, и здесь это доведено до арифметики.
SPECIAL_CLASSES = (NumberClass.ZERO, NumberClass.UNCOMPUTABLE,
                   NumberClass.INFINITESIMAL, NumberClass.INFINITY,
                   NumberClass.IGNORANCE)


def _class3_exponent(m: Q, alpha, e):
    """Показатель с линейки 3 (04.09, дверь Гильберта 7). Рациональный на
    линейке — обычная пара (2^θ при θ = √4 есть 2²). Иррациональный над
    алгебраическим основанием — Гельфонд–Шнайдер (1934): значение
    трансцендентно, отказ теоремой (1^θ = 1 и 0^θ = 0 при θ > 0 сведены).
    Над трансцендентным основанием — отказ поимённо: нормальная форма
    ratfunc над β = α^(1/N) требует рационального показателя. Особое
    основание (ω, ε) — как и дробная степень: не определено."""
    from .algebraic import AlgebraicNumber, _power_by_theorem
    r = e.as_fraction()
    if r is not None:
        return m, alpha, r
    if is_special(alpha):
        raise ExactRepresentationNotAvailable(
            f"степень класса 3 основания класса "
            f"{alpha.number_class.name} не определена")
    base = Q(alpha) if isinstance(alpha, int) else alpha
    if isinstance(base, (Q, AlgebraicNumber)):
        v = _power_by_theorem(base, e)      # 1 или 0 — иначе отказ теоремой
        return m * v, Q(1), Q(0)
    rec = alpha.record() if hasattr(alpha, "record") else repr(alpha)
    raise Unresolved("^", rec, e.record(),
                     hint="показатель класса 3 над трансцендентным основанием: "
                          "нормальная форма (ratfunc) требует рационального "
                          "показателя, носителя нет")


def is_special(v) -> bool:
    """Носитель особого класса: его арифметику знает контейнер, не запись."""
    return (isinstance(v, PhiNumberV6)
            and v.number_class in SPECIAL_CLASSES)


def _as_container(q):
    """Пара — контейнером, чтобы умножать особого носителя на неё им же."""
    return PhiNumberV6.from_ratio(int(q.n), int(q.d))


def _normal_form(alpha, e, N: int | None = None):
    """α^E над трансцендентной базой как элемент ℚ(x), x = α^(1/N) —
    (RatFunc, N). Дробный коэффициент члена Цеха в ℚ(x) не живёт — отказ."""
    from .ratfunc import RatFunc
    N = _exp_lcm(e, 1) if N is None else N

    def of(e):
        if isinstance(e, ZechExponent):
            out = RatFunc.mono(int(e.q * N))
            for c, z in e.terms:
                if c.d != 1:
                    raise ExactRepresentationNotAvailable(
                        "нормальная форма: дробная степень множителя "
                        f"(1 + k·α^d)^({c}) в ℚ(x) не сведена")
                out = out * (RatFunc.const(1) + of(z.d) * z.k) ** c.n
            return out
        return RatFunc.mono(int(Q(e) * N))

    return of(e), N


def _liftable(x) -> bool:
    """Носители, которые ФIOS поднимает в запись над их базой (строки
    alpha_lift): классы 5, 6, 10."""
    return type(x) in _carriers()


# ---------- двери Φα против мира IEEE: одна копия правила отказа ----------

#: Причина отказа — по ПРЕДМЕТУ, не по типу; названа один раз и берётся
#: отсюда каждой дверью. Шесть копий правила в дверях были шестью
#: молчаниями: копия, которой нет, отказать не может.
_IEEE_WHY = {
    "float": "float запрещён построением; явная дверь — Q.from_float "
             "(rational.py:39): float несёт разрядность, которой у величины "
             "нет, и принять его молча значит впустить приближение в точную "
             "величину — потеря при явной двери записана классом 23",
    "bool": "bool не число (rational.py:37): bool есть подкласс int, потому "
            "isinstance(x, int) пропускает True как 1 и False как 0 — "
            "величина подменяется истинностным значением",
}


def _ieee_name(x):
    """Предмет из мира IEEE или None. bool ПРОВЕРЯЕТСЯ ПЕРВЫМ: он подкласс
    int (и не подкласс float), и порядок проверок здесь и есть та ловушка,
    о которой говорит `_IEEE_WHY['bool']`."""
    if isinstance(x, bool):
        return "bool"
    if isinstance(x, float):
        return "float"
    return None


def _refuse_ieee(x, door: str, takes: str, depth: int = 1) -> None:
    """Отказ двери Φα предмету из мира IEEE: ПРЕДМЕТ, ПРИЧИНА, АДРЕС.

    Правило дома: отказ обязан быть НАЗВАН — что пришло, почему не принято
    и откуда отказ (§316.4: несработавшее правило даёт отказ, а не класс).
    Одна копия правила на все двери; `door` и `takes` несут только имя
    места и что оно принимает — причина здесь не повторяется.

    Адрес берётся КАДРОМ вызывающей двери (`sys._getframe(depth)`), а не
    строковой пометкой: пометка врёт после всякой правки файла, кадр не
    врёт никогда. `depth` — сколько кадров до двери: 1 — дверь зовёт канон
    сама, 2 — через `_refuse_ieee_build` (иначе адресом встал бы сам канон,
    и отказ называл бы место ПРАВИЛА вместо места ДВЕРИ: замер 12.09 —
    все три слота построения показывали одну строку 546).

    Замеры (все — [вычислено здесь], 12.09, строка 79 очереди):

    * двери Φα, принимающие чужой предмет: 13 МЕСТ ВХОДА — три слота
      построения (`m`, `α`, `E`), девять действий и равенство; связанных
      имён 11 (двенадцатое, `__neg__`, операнда не принимает), а различных
      функций 9 (`__mul__` = `__rmul__`, `__add__` = `__radd__`). Входов
      «место × предмет» — 26. Прежний замер знал восемь входов (четыре
      действия × два предмета) и не считал ни отражённых имён, ни слотов
      построения, ни равенства;
    * ДО правки дом отказывал САМ на 8 входах из 26 (построение 6: `m` и
      `e` через пару, `alpha` — своей проверкой; степень 2), ПРЕДМЕТ
      называл на 6 (степень не называла: «показатель — int или пара» —
      о том, что дверь принимает, а не о том, что пришло), АДРЕС не давал
      ни на одном; 16 входов отдавались интерпретатору голым «unsupported
      operand type(s)»; 2 входа (равенство) молча отвечали False;
    * ПОСЛЕ дом отказывает сам на 24 входах из 26 — предмет, причина и
      адрес в каждом; интерпретатору не отдан ни один вход. Счёт обеих
      сторон — `proto/test_phicore_alpha_doors.py`;
    * корень ловушки bool в доме: `rational.py:37` (пара), `rational.py:173`
      и `:241` (степень пары), `container.py:439` (bool и float не поднимаются в
      контейнер);
    * канон соседа по степени: `container.py:1245-1257` — отказ несёт обе
      Ф-записи (строка 39 очереди, замер 11.09);
    * прежний замер главы 09 (адреса ДО этой правки): §17 «Двери» и §21
      п. 11 — `docs/УЧЕБНИК/09-контейнер-и-ось-основания.md:710-716` и
      `:998-1002`, где двери названы как `alpha.py:757, :894, :910, :934`;
      счётчик разнобоя — `docs/УЧЕБНИК/прогоны/09-вердикты.py:740`
      (контейнер: названо 2, интерпретатором 6; у Φα было то же 2/6).

    Равенство (`__eq__`, ниже) СЮДА НЕ ЗОВЁТ, и это решение, не пропуск:
    канон соседей — определённый ответ False, а не отказ (`Q(3) == 3.0` →
    False, `rational.py:305`; контейнер то же), и отказ в одном Φα развёл
    бы соседей по одному и тому же вопросу. Названо остатком, не закрыто
    молча.
    """
    if not isinstance(x, (bool, float)):
        return                      # быстрый путь двери: чужого предмета нет
    kind = _ieee_name(x)
    import sys
    raise TypeError(f"alpha.py:{sys._getframe(depth).f_lineno} дверь Φα.{door}: "
                    f"{takes}; получен {kind} {x!r} — {_IEEE_WHY[kind]}")


#: Что принимает каждый слот построения — данные двери, не правило.
_BUILD_SLOTS = (("__init__(m=…)", "мантисса — пара (Q), int или строка пары"),
                ("__init__(α=…)", "основание — пара или носитель, не float/bool"),
                ("__init__(E=…)", "показатель — пара, int или показатель Цеха"))


def _refuse_ieee_build(m, alpha, e) -> None:
    """Тот же канон для трёх слотов построения: отказывает ПЕРВЫЙ слот, в
    который пришёл предмет из мира IEEE, и называет, какой это слот.

    Цена двери замерена (12.09, машина под воркерами, min из пяти прогонов
    по 20000 вызовов): перебор всех трёх слотов вызовами стоил построению
    Φα 1.6 → 2.1 мкс. Быстрый путь ниже — ТА ЖЕ одна копия правила, но
    чистое построение за перебор слотов не платит."""
    if not (isinstance(m, (bool, float)) or isinstance(alpha, (bool, float))
            or isinstance(e, (bool, float))):
        return                      # быстрый путь построения
    for value, (door, takes) in zip((m, alpha, e), _BUILD_SLOTS):
        if _ieee_name(value) is not None:
            _refuse_ieee(value, door, takes, depth=2)   # адрес — место ДВЕРИ, не правила


class PhiAlpha:
    """Φ = m·α^e с произвольным основанием-носителем. Запись, не новый
    класс: значение принадлежит классу своего носителя, в онтологию
    ФIOS входит переводом `as_carrier()`."""

    __slots__ = ("m", "alpha", "e")

    def __init__(self, m, alpha, e):
        _refuse_ieee_build(m, alpha, e)         # канон дверей: предмет, причина, адрес
        m = m if type(m) is Q else Q(m)
        if isinstance(e, ZechExponent) and e.is_rational:
            e = e.q
        zech = isinstance(e, ZechExponent)
        if not zech and hasattr(e, "as_fraction") and hasattr(e, "p"):
            m, alpha, e = _class3_exponent(m, alpha, e)   # линейка 3 в показателе
        if not zech and type(e) is not Q:
            e = Q(e)
        if alpha is None:
            raise TypeError("основание не задано: None — не носитель")
        if isinstance(alpha, PhiAlpha) or (isinstance(alpha, PhiNumberV6) and alpha.base_class is not None):
            # слово автора 05.09 («α ≠ 2»): основание — мантисса классов 1–19 над базой 2,
            # не Ф-число со своим основанием; циклическая ссылка самоопределения недопустима
            raise ValueError("основание — мантисса классов 1–19 над базой 2, не Ф-число со своим α: "
                             "циклическое самоопределение недопустимо")
        if is_special(alpha):
            # Особый носитель основанием допустим только при целом показателе:
            # ω^½ и ε^⅓ библиотекой не определены, и это отказ поимённо, а не
            # молчаливая подстановка (§316.4 — несработавшее правило обязано
            # давать отказ, а не класс).
            if zech:
                raise Unresolved(
                    "log по особому основанию",
                    f"Φα(α={alpha.number_class.name}, E=…)",
                    hint="логарифм сложения по основанию класса "
                         f"{alpha.number_class.name} не сведён")
            if Q(e).d != 1:
                raise ExactRepresentationNotAvailable(
                    f"дробная степень {e} основания класса "
                    f"{alpha.number_class.name} не определена")
        if isinstance(alpha, int):
            alpha = Q(alpha)

        if type(alpha) is Q:
            if alpha.n == 0:
                if zech:
                    raise Unresolved("0^E", f"Φα(α=0, E={e.record()})",
                                     hint="нулевое основание при "
                                          "логарифмическом показателе не сведено")
                if e < 0:
                    raise ZeroDivisionError("0 в отрицательной степени")
                if e == 0:
                    raise Unresolved("0^0", "Φα(α=0, E=0)",
                                     hint="0⁰ не сведено: ни нуля, ни единицы")
                m = Q(0)                      # 0^e = 0, e>0
            elif alpha.n < 0:
                if zech:
                    raise ExactRepresentationNotAvailable(
                        "знак отрицательного основания при логарифмическом "
                        "показателе не отделён")
                if e.d % 2 == 0:
                    raise ExactRepresentationNotAvailable(
                        f"чётный корень из отрицательного основания "
                        f"{alpha} не веществен; пара (re, im) — на линейках 1/3 "
                        f"(complex_pair), над основанием Φα не строится")
                if e.n % 2:
                    m = -m                    # знак — чётностью числителя
                alpha = -alpha
        if m.n == 0:
            alpha, e, zech = Q(1), Q(0), False
        if type(alpha) is Q:
            if alpha.n == alpha.d or (not zech and e.n == 0):
                alpha, e = Q(1), Q(0)         # тривиальное состояние
            else:
                if alpha.n < alpha.d:
                    alpha, e = 1 / alpha, exp_rebase(e, -1)   # 1/Ф=Ф/1
                beta, k = _primitive_base(alpha)
                if k != 1:
                    alpha, e = beta, exp_rebase(e, k)
        elif not zech and e.n == 0:
            alpha = Q(1)                      # фиктивность основания
        elif is_special(alpha):
            pass                              # особый носитель: как есть
        else:
            m, alpha, e = _canon_carrier(m, alpha, e)
            if isinstance(e, ZechExponent) and e.is_rational:
                e = e.q
            zech = isinstance(e, ZechExponent)
            if isinstance(alpha, Q):          # носитель выродился в пару
                if m != 0 and not (alpha == 1 or (not zech and e == 0)):
                    alpha, e = (1 / alpha, exp_rebase(e, -1)) if alpha < 1 \
                        else (alpha, e)
                    beta, k = _primitive_base(alpha)
                    alpha, e = beta, exp_rebase(e, k)
                else:
                    alpha, e = Q(1), Q(0)
        if not isinstance(alpha, Q) and not is_special(alpha) \
                and isinstance(e, ZechExponent):
            # канон: нормальная форма — моном ⇒ показатель рационален
            try:
                nf, N = _normal_form(alpha, e)
            except ExactRepresentationNotAvailable:
                nf = None
            mono = nf.is_monomial if nf is not None else None
            if mono is not None:
                c, n = mono
                m, e, zech = m * c, Q(n, N), False
                if e == 0:
                    alpha = Q(1)
        object.__setattr__(self, "m", m)
        object.__setattr__(self, "alpha", alpha)
        object.__setattr__(self, "e", e)

    def __setattr__(self, *_):
        raise AttributeError("Φα неизменяем")

    # --- запись ---------------------------------------------------------
    def record(self) -> str:
        a = repr(self.alpha) if isinstance(self.alpha, Q) \
            else self.alpha.record() if hasattr(self.alpha, "record") \
            else repr(self.alpha)
        e = self.e.record() if isinstance(self.e, ZechExponent) else repr(self.e)
        return f"Φα(M={self.m}, α={a}, E={e})"

    __repr__ = record

    # --- перевод в носителя: дверь в онтологию --------------------------
    def as_carrier(self):
        """Точное значение носителем его класса: m·(α^e).

        Несведённый перевод — честный отказ построением, поимённо."""
        a, e = self.alpha, self.e
        if isinstance(e, ZechExponent):
            from .algebraic import AlgebraicNumber
            if isinstance(a, AlgebraicNumber):
                # алгебраическая база: форма ℚ(x) как выражение (не канон —
                # у β есть соотношение), значение — на линейке 3 в точке
                # β = α^(1/N): P(β)/Q(β) точно, результантами eval_poly
                nf, N = _normal_form(a, e)
                beta = a ** Q(1, N) if N > 1 else a
                if not isinstance(beta, (Q, AlgebraicNumber)):
                    # ОТРИЦАТЕЛЬНОЕ АЛГЕБРАИЧЕСКОЕ основание: чётный корень
                    # из него вещественным не бывает, и β выходит парой
                    # (re, im). До 12.09 этот вход дверь проходил и падал
                    # здесь сырым TypeError («unsupported operand type(s)
                    # for **: 'ComplexPair' and 'int'»), тогда как
                    # отрицательная ПАРА в основании отказывает поимённо
                    # ещё построением (ветка alpha.n < 0 выше). Теперь
                    # форма отказа одна и та же на оба входа; считать
                    # значение парой нельзя, не правя complex_pair.
                    raise ExactRepresentationNotAvailable(
                        f"чётный корень (N={N}) из отрицательного основания "
                        f"{a.record()} не веществен; пара (re, im) — на "
                        f"линейках 1/3 (complex_pair), над основанием Φα "
                        f"не строится")
                num = beta.eval_poly(nf.num) if isinstance(beta, AlgebraicNumber) \
                    else sum((c * beta ** i for i, c in enumerate(nf.num)), Q(0))
                den = beta.eval_poly(nf.den) if isinstance(beta, AlgebraicNumber) \
                    else sum((c * beta ** i for i, c in enumerate(nf.den)), Q(0))
                if den == 0:
                    raise ExactRepresentationNotAvailable(
                        "перевод: знаменатель формы обращается в ноль в β — "
                        "запись не определена")
                inv = den.inverse() if isinstance(den, AlgebraicNumber) else 1 / den
                return num * inv * self.m
            if not isinstance(a, Q):
                raise ExactRepresentationNotAvailable(
                    f"перевод: сумма над основанием-носителем "
                    f"{type(a).__name__} носителя не имеет; запись точна как "
                    f"есть (нормальная форма {self.normal_form()[0].record()})"
                    if _is_transcendental(a) else
                    f"перевод: сумма над основанием-носителем "
                    f"{type(a).__name__} носителя не имеет; запись точна как есть")
            try:
                poly, N = self._radical()
            except ExactRepresentationNotAvailable:
                # дробная степень множителя — медленный путь результантами
                val = self.m * (a ** e.q)
                for c, z in e.terms:
                    x = 1 + z.k * PhiAlpha(1, a, z.d).as_carrier()
                    val = val * (x ** c)
                return val
            if all(c == 0 for c in poly[1:]):
                return poly[0]
            from .algebraic import AlgebraicNumber
            beta = AlgebraicNumber.from_root(N, a)
            if isinstance(beta, Q):                  # α — полная N-я степень
                return sum((c * beta ** i for i, c in enumerate(poly)), Q(0))
            return beta.eval_poly(poly)
        if is_special(a):
            # Значение принадлежит классу носителя: степень и множитель
            # считает сам контейнер. Умножать носителя на пару нельзя —
            # пара переводится в контейнер (иначе TypeError, а не отказ).
            #
            # `[вычислено здесь]` Здесь стоял обход: E=1 не возводилось,
            # потому что у класса 23 `x ** 1` терял провенанс
            # (`IGNORANCE:inf` → `IGNORANCE:None` → `nan`). Дефект эталона
            # снят 03.09 ветвью степени в контейнере; замер 05.09 на всех
            # особых классах (0, 20, 21, 22, 23) даёт x¹ = x тождественно,
            # потому обход убран, а не оставлен «на всякий случай»:
            # мёртвая ветвь прячет, снят дефект или нет.
            return (a ** int(e)) * _as_container(self.m)
        if isinstance(a, Q):
            return self.m * (a ** e)          # пара/линейка 3 — точно
        from .algebraic import AlgebraicNumber
        from .exp_log import ExpValue
        from .pi_value import PiMonomial
        if isinstance(a, AlgebraicNumber):
            return (a ** e) * self.m
        if isinstance(a, PiMonomial):
            c = a.coeff ** e if isinstance(a.coeff, Q) else None
            if c is None:
                raise ExactRepresentationNotAvailable(
                    "перевод: алгебраический коэффициент монома π в "
                    "дробной степени не сведён")
            return PiMonomial._make(c * self.m, a.e * e)
        if isinstance(a, ExpValue):
            c = a.coeff ** e
            if not isinstance(c, Q):
                raise ExactRepresentationNotAvailable(
                    f"перевод: коэффициент {a.coeff} экспоненты в степени "
                    f"{e} не остаётся парой")
            return ExpValue.make(c * self.m, a.y * e)
        if e.d == 1:
            return (a ** int(e)) * self.m     # целая степень носителя
        raise ExactRepresentationNotAvailable(
            f"перевод: дробная степень носителя "
            f"{type(a).__name__} не сведена")

    def _radical(self, N: int | None = None):
        """m·α^E над рациональной базой как элемент ℚ(α^(1/N)) — (poly, N).
        Один общий N для всех показателей; сравнение двух записей над
        одной базой — сравнение списков, без результантов."""
        if not isinstance(self.alpha, Q):
            raise ExactRepresentationNotAvailable("радикальное поле — только "
                                                  "над рациональной базой")
        N = _exp_lcm(self.e, 1) if N is None else N
        F = _Radical(N, self.alpha)
        return F.scale(F.of_exponent(self.e), self.m), N

    def enclosure(self, terms: int = 30):
        """Сертификат сетки записи в общем интерфейсе носителей: носитель
        есть — его сертификат; над трансцендентной базой с показателем
        Цеха — форма ℚ(x) в засечке β = α^(1/N) (интервальный Горнер)."""
        from .enclosure import Interval
        try:
            v = self.as_carrier()
        except (ArithmeticError, Unresolved):
            v = None
        if v is not None:
            if isinstance(v, Q):
                return Interval.point(v)
            enc = getattr(v, "enclosure", None)
            if callable(enc):
                return enc(terms)
            raise Unresolved("enclosure", self.record(),
                             hint="у носителя записи нет сертификата сетки")
        if _is_transcendental(self.alpha):
            nf, N = self.normal_form()
            return nf.value_at(_beta_enclosure(self.alpha, N, terms)) * self.m
        raise Unresolved("enclosure", self.record(),
                         hint="сертификат записи над этой базой не сведён")

    def normal_form(self, N: int | None = None):
        """m·α^E над трансцендентной базой как несократимая дробь ℚ(x),
        x = α^(1/N) — (RatFunc, N). Равенство чисел = равенство форм."""
        if not _is_transcendental(self.alpha):
            raise ExactRepresentationNotAvailable(
                "нормальная форма ℚ(x) — только над трансцендентной базой "
                "(над рациональной — радикальное поле _radical)")
        nf, N = _normal_form(self.alpha, self.e, N)
        return nf * self.m, N

    # --- арифметика: одноимённая база в контейнере ----------------------
    def _same_base(self, other: "PhiAlpha") -> bool:
        a, b = self.alpha, other.alpha
        if a is b:
            return True
        if type(a) is Q and type(b) is Q:
            return a.n == b.n and a.d == b.d
        return a == b

    def __mul__(self, other):
        _refuse_ieee(other, "__mul__ (она же __rmul__)",
                     "множитель — Φα, int, пара или носитель классов 5/6/10")
        if isinstance(other, PhiAlpha):
            if type(self.e) is Q and self.e.n == 0:
                return PhiAlpha(self.m * other.m, other.alpha, other.e)
            if type(other.e) is Q and other.e.n == 0:
                return PhiAlpha(self.m * other.m, self.alpha, self.e)
            if self._same_base(other):
                return PhiAlpha(self.m * other.m, self.alpha,
                                exp_add(self.e, other.e))
            raise CrossBaseOperationNotSupported(
                f"базы {self.record()} и {other.record()} различны; "
                f"сведение — переводом as_carrier (ФIOS зовёт его сам)")
        if isinstance(other, (int, Q)):
            return PhiAlpha(self.m * Q(other), self.alpha, self.e)
        if _liftable(other):
            from .Фios import route
            return route(self, other, "*")
        return NotImplemented

    __rmul__ = __mul__

    # --- сложение: одна база — по мантиссе или логарифмом сложения ------
    def _sign_of_one_plus(self, k, d) -> int:
        """Знак X = 1 + k·α^d — вторая таблица Гаусса (логарифм вычитания).

        Рациональная и алгебраическая база — точно (пара / линейка 3).
        Основание-носитель (e, π): сертификатом сетки носителя α^d с
        уточнением по явному запросу; ноль сертификат не доказывает, но
        для трансцендентной базы и рационального k, d≠0 нуля и нет
        (α^d рационален только при d=0). Не отделилось — честный отказ."""
        a = self.alpha
        if _is_transcendental(a):
            return self._sign_by_form(k, d)
        try:
            v = PhiAlpha(1, a, d).as_carrier()
        except ArithmeticError as ex:
            raise Unresolved("−", self.record(), f"α^({d})",
                             hint=f"логарифм вычитания: α^d не переводится "
                                  f"({type(ex).__name__})") from ex
        from .algebraic import AlgebraicNumber
        if isinstance(v, Q):
            x = 1 + k * v
            return (x > 0) - (x < 0)
        if isinstance(v, AlgebraicNumber):
            return (1 + k * v)._separate_sign()
        if not hasattr(v, "enclosure"):
            raise Unresolved("−", self.record(), v.record()
                             if hasattr(v, "record") else repr(v),
                             hint="логарифм вычитания: у носителя нет "
                                  "сертификата сетки")
        for terms in (30, 60, 120, 240, 480):
            iv = 1 + v.enclosure(terms) * k
            if iv.lo > 0:
                return 1
            if iv.hi < 0:
                return -1
        raise Unresolved("−", self.record(), v.record(),
                         hint="логарифм вычитания: знак 1 + k·α^d не "
                              "отделился сертификатом до 480 членов")

    def _sign_by_form(self, k, d) -> int:
        """Знак 1 + k·α^d над трансцендентной базой: нормальная форма в
        ℚ(x) — ноль точно (P = 0), иначе интервальный Горнер по сертификату
        β = α^(1/N) с уточнением 30→480 членов."""
        from .ratfunc import RatFunc
        try:
            f, N = _normal_form(self.alpha, d)
        except ExactRepresentationNotAvailable as ex:
            raise Unresolved("−", self.record(), f"α^({d})",
                             hint=f"логарифм вычитания: {ex}") from ex
        g = RatFunc.const(1) + f * k
        if g.is_zero:
            return 0
        for terms in (30, 60, 120, 240, 480):
            sgn = g.sign_at(_beta_enclosure(self.alpha, N, terms))
            if sgn is not None:
                return sgn
        raise Unresolved("−", self.record(), g.record(),
                         hint="логарифм вычитания: знак 1 + k·α^d не "
                              "отделился сертификатом до 480 членов")

    def _gauss(self, m2, e2):
        """Логарифм сложения (Гаусс 1812, Цех): m₁α^e₁ + m₂α^e₂ =
        m₁·α^(e₁ + log_α(1 + (m₂/m₁)·α^(e₂−e₁))). Ведущий — меньший
        показатель, когда оба рациональны (канон коммутативности).
        Рациональная база и целая разность — множитель рационален и
        уходит в мантиссу. Вычитание (k<0) — вторая таблица Гаусса: знак
        1 + k·α^d над рациональной/алгебраической базой точно, над
        основанием-носителем — сертификатом сетки (см. _sign_of_one_plus)."""
        a = self.alpha
        m1, e1 = self.m, self.e
        m2 = Q(m2)
        if m2 == 0:
            return self
        if isinstance(e1, Q) and isinstance(e2, Q) and e2 < e1:
            m1, e1, m2, e2 = m2, e2, m1, e1
        k = m2 / m1
        d = exp_add(e2, exp_neg(e1))
        if isinstance(a, Q) and isinstance(d, Q) and d.d == 1:
            return PhiAlpha(m1 * (1 + k * a ** d), a, e1)
        if k < 0:
            sgn = self._sign_of_one_plus(k, d)
            if sgn == 0:
                return Q(0)
            if sgn < 0:
                m1, e1, m2, e2 = m2, e2, m1, e1
                k = m2 / m1
                d = exp_add(e2, exp_neg(e1))
        E = exp_add(e1, ZechExponent(0, ((Q(1), Zech(k, d)),)))
        return PhiAlpha(m1, a, E)

    def _special_sum(self, other, minus: bool):
        """Сумма с участием особого носителя — контейнером, эталоном.

        Возвращает НОСИТЕЛЯ, а не Φα: значение принадлежит классу своего
        носителя, и запись над ω или ε записью уже не является.
        """
        a = self.as_carrier() if isinstance(self, PhiAlpha) else self
        if isinstance(other, PhiAlpha):
            b = other.as_carrier()
        elif isinstance(other, PhiNumberV6):
            b = other
        else:
            b = _as_container(Q(other))
        return a - b if minus else a + b

    def _has_special(self, other) -> bool:
        return (is_special(self.alpha) or is_special(other)
                or (isinstance(other, PhiAlpha) and is_special(other.alpha)))

    def __add__(self, other):
        """m·α^e + m'·α^e — по мантиссе; одна база, разные показатели —
        логарифм сложения; пара при e≠0 — тем же Гауссом (q = q·α^0);
        разные базы — переводом в ФIOS (носитель сведёт либо честный
        отказ по строке таблицы)."""
        _refuse_ieee(other, "__add__ (она же __radd__)",
                     "слагаемое — Φα, int, пара или носитель классов 5/6/10")
        if self._has_special(other):
            return self._special_sum(other, minus=False)
        if isinstance(other, PhiAlpha):
            if other.e == 0:
                return self + other.m
            if self.e == 0:
                return other + self.m
            if self._same_base(other):
                if self.e == other.e:
                    return PhiAlpha(self.m + other.m, self.alpha, self.e)
                return self._gauss(other.m, other.e)
        elif isinstance(other, (int, Q)):
            q = Q(other)
            if q == 0:
                return self
            if self.e == 0:
                return PhiAlpha(self.m + q, Q(1), 0)
            return self._gauss(q, Q(0))
        elif not _liftable(other):
            return NotImplemented
        from .Фios import route
        return route(self, other, "+")

    __radd__ = __add__

    def __sub__(self, other):
        _refuse_ieee(other, "__sub__", "вычитаемое — Φα, int, пара или "
                     "носитель классов 5/6/10")
        if self._has_special(other):
            return self._special_sum(other, minus=True)
        if isinstance(other, PhiAlpha):
            return self + (-other)
        if isinstance(other, (int, Q)):
            return self + (-Q(other))
        if _liftable(other):
            from .Фios import route
            return route(self, other, "-")
        return NotImplemented

    def __rsub__(self, other):
        _refuse_ieee(other, "__rsub__", "уменьшаемое — int, пара или "
                     "носитель классов 5/6/10")
        if isinstance(other, (int, Q)):
            return (-self) + Q(other)
        if _liftable(other):
            from .Фios import route
            return route(other, self, "-")
        return NotImplemented

    def inverse(self) -> "PhiAlpha":
        """1/Φ = (1/m)·α^(−e); при α=1, m=±1 — Ф=1/Ф (тривиальное)."""
        if self.m == 0:
            raise ZeroDivisionError("ноль не обратим")
        return PhiAlpha(1 / self.m, self.alpha, exp_neg(self.e))

    def __truediv__(self, other):
        _refuse_ieee(other, "__truediv__", "делитель — Φα, int, пара или "
                     "носитель классов 5/6/10")
        if isinstance(other, PhiAlpha):
            return self * other.inverse()
        if isinstance(other, (int, Q)):
            return PhiAlpha(self.m / Q(other), self.alpha, self.e)
        if _liftable(other):
            from .Фios import route
            return route(self, other, "/")
        return NotImplemented

    def __rtruediv__(self, other):
        _refuse_ieee(other, "__rtruediv__", "делимое — int, пара или "
                     "носитель классов 5/6/10")
        if isinstance(other, (int, Q)):
            return self.inverse() * Q(other)
        if _liftable(other):
            from .Фios import route
            return route(other, self, "/")
        return NotImplemented

    def __neg__(self):
        return PhiAlpha(-self.m, self.alpha, self.e)

    def __pow__(self, p):
        """Степень точная всегда: ось e несёт дробь; мантисса без
        полного корня уходит значением на линейку (запись меняется,
        значение точно)."""
        _refuse_ieee(p, "__pow__", "показатель — int или пара")
        if isinstance(p, int):
            p = Q(p)
        if not isinstance(p, Q):
            return NotImplemented
        mp = self.m ** p                      # Q либо переход на линейку 3
        ep = exp_scale(self.e, p)
        if isinstance(mp, Q):
            return PhiAlpha(mp, self.alpha, ep)
        rest = PhiAlpha(1, self.alpha, ep).as_carrier()
        return mp * rest if rest != 1 else mp

    # --- равенство: по значению, основание фиктивно ---------------------
    def __eq__(self, other):
        if isinstance(other, PhiAlpha):
            if (self.m == other.m and self.e == other.e
                    and self._same_base(other)):
                return True
            if (isinstance(self.alpha, Q) and isinstance(other.alpha, Q)
                    and self.alpha == other.alpha):
                from math import lcm
                N = lcm(_exp_lcm(self.e, 1), _exp_lcm(other.e, 1))
                try:
                    return self._radical(N)[0] == other._radical(N)[0]
                except ExactRepresentationNotAvailable:
                    pass
            if _is_transcendental(self.alpha) and self._same_base(other):
                from math import lcm
                N = lcm(_exp_lcm(self.e, 1), _exp_lcm(other.e, 1))
                try:
                    return self.normal_form(N)[0] == other.normal_form(N)[0]
                except ExactRepresentationNotAvailable:
                    pass
            try:
                other = other.as_carrier()
            except (ArithmeticError, Unresolved):
                return NotImplemented
        try:
            mine = self.as_carrier()
        except (ArithmeticError, Unresolved):
            return NotImplemented
        return mine == other

    # --- самопроверка записи (ось самоприменения, 08.10) ------------------

    def self_check(self, deep: bool = False) -> list:
        """Самопроверка записи встречными дорогами (ось «AIAΣA спрашивает
        AIAΣA», слово автора 04.09; очередь строка 18; 08.10). Список
        находок, пусто — здорова; исключений не бросает.

        `[вычислено здесь]` (1) слоты: m — Q, e — Q либо показатель Цеха, α —
        носитель (не Φα и не контейнер со своим α); (2) канон построения над
        рациональной базой: α ≥ 1, α = 1 ⟺ e = 0 (при m ≠ 0), α первообразна
        (`_primitive_base` даёт показатель 1); (3) ВСТРЕЧНЫЕ дороги:
        фиктивность основания — та же величина над α² с показателем e/2 равна
        записи (равенство решает радикальное поле либо нормальная форма, не
        совпадение слотов); инверсия — запись, умноженная на свою инверсию,
        есть единица, и инверсия инверсии есть запись; сертификат —
        квадрат отрезка записи пересекает сертификат квадрата записи
        (умножение на оси E против интервальной арифметики); отказ перевода
        в носителя (несведённая база) — объявленная граница. При deep —
        ещё куб и 60 членов."""
        out = []
        m, a, e = self.m, self.alpha, self.e
        if not isinstance(m, Q):
            return [f"слот m чужого носителя: {type(m).__name__}"]
        if not isinstance(e, (Q, ZechExponent)):
            return [f"слот e чужого носителя: {type(e).__name__}"]
        if a is None or isinstance(a, PhiAlpha):
            return [f"слот α вне границы входа: {a!r}"]
        if m == 0:
            return out
        if isinstance(a, Q) and isinstance(e, Q):
            if a < 1:
                out.append(f"канон: α = {a} < 1 (двойственность 1/Φ обязана была обратить базу)")
            if (a == 1) != (e == 0):
                out.append(f"канон: α = {a} и e = {e} — тривиальное состояние требует α = 1 ⟺ e = 0")
            if a > 1 and _primitive_base(a)[1] != 1:
                out.append(f"канон: α = {a} не первообразна (α = β^{_primitive_base(a)[1]})")
        try:
            if isinstance(e, Q) and not isinstance(a, Q) or (isinstance(a, Q) and isinstance(e, Q) and a > 1):
                twin = PhiAlpha(m, a * a if isinstance(a, Q) else a ** 2, e / 2)
                if not (twin == self):
                    out.append(f"фиктивность основания: запись над α² с e/2 ({twin.record()}) не равна {self.record()}")
        except (Unresolved, ExactRepresentationNotAvailable, ArithmeticError):
            pass                                                # граница: равенство над этой базой не сведено
        except Exception as ex:                                 # noqa: BLE001
            out.append(f"фиктивность основания не проверена: {type(ex).__name__}: {ex}")
        try:
            inv = self.inverse()
            if not (self * inv == Q(1)):
                out.append(f"инверсия: Φ·(1/Φ) = {(self * inv)!r} ≠ 1")
            if not (inv.inverse() == self):
                out.append("инверсия инверсии не равна записи")
        except (Unresolved, ExactRepresentationNotAvailable, ArithmeticError, CrossBaseOperationNotSupported):
            pass
        except Exception as ex:                                 # noqa: BLE001
            out.append(f"инверсия не проверена: {type(ex).__name__}: {ex}")
        try:
            terms = 60 if deep else 30
            iv = self.enclosure(terms)
            sq = (self * self).enclosure(terms)
            if (iv * iv).intersect(sq) is None:
                out.append(f"сертификат: квадрат отрезка {(iv * iv).record()} не пересекает сертификат квадрата {sq.record()}")
            if deep:
                cube = (self * self * self).enclosure(terms)
                if (iv * iv * iv).intersect(cube) is None:
                    out.append("сертификат: куб отрезка не пересекает сертификат куба")
        except (Unresolved, ExactRepresentationNotAvailable, ArithmeticError, CrossBaseOperationNotSupported):
            pass                                                # граница: сертификат над этой базой не сведён
        except Exception as ex:                                 # noqa: BLE001
            out.append(f"сертификат не проверен: {type(ex).__name__}: {ex}")
        return out

    def __hash__(self):
        try:
            return hash(self.as_carrier())
        except (ArithmeticError, Unresolved):
            pass
        if _is_transcendental(self.alpha):
            try:
                nf, N = self.normal_form()
                return hash(nf.reduce_index(N))
            except ExactRepresentationNotAvailable:
                pass
        return hash((self.m, self.e))


# ---------- самопроверка записи α (ось самоприменения, 08.10) ----------

def self_check(deep: bool = False) -> list:
    """Самопроверка контейнера с произвольным основанием встречными
    дорогами; список находок, пусто — здоров; исключений не бросает.

    `[вычислено здесь]` Фиктивность основания — одно значение в пяти записях
    (α = 2, 4, 16, 1/2, 256) равно себе и двойке; граничные условия
    построения: α = 1 гасит показатель, e = 0 гасит базу, α = 0 при e > 0 —
    точный ноль, m = 0 — ноль; первообразная база: 8^(1/3) и 2^1 — одна
    запись; окно записей — 3·2^(1/2), 2^(−3/2), 5·3^(2/3), 7^(1/4) и над
    трансцендентной базой π² и e^(1/2) (deep — ещё 2^(5/6) и π^(3/2))
    проверяют себя (`PhiAlpha.self_check`). Граница — окно."""
    out = []
    try:
        forms = (PhiAlpha(1, 2, 1), PhiAlpha(1, 4, Q(1, 2)), PhiAlpha(1, 16, Q(1, 4)), PhiAlpha(1, Q(1, 2), -1), PhiAlpha(1, 256, Q(1, 8)))
        if not all(f == forms[0] and f == Q(2) for f in forms):
            out.append("фиктивность основания: пять записей двойки не равны между собой")
        if PhiAlpha(5, 1, Q(99)) != Q(5) or PhiAlpha(3, 7, 0) != Q(3) or PhiAlpha(3, 7, 0).alpha != Q(1):
            out.append("граничные условия: α = 1 или e = 0 не гасят ось")
        if PhiAlpha(1, 0, Q(1, 2)) != Q(0) or PhiAlpha(0, 5, Q(3)) != Q(0):
            out.append("граничные условия: α = 0 при e > 0 либо m = 0 не дают точный ноль")
        if PhiAlpha(1, 8, Q(1, 3)) != PhiAlpha(1, 2, 1) or PhiAlpha(1, 8, Q(1, 3)).alpha != Q(2):
            out.append("первообразная база: 8^(1/3) не сведена к 2^1")
    except Exception as e:                                      # noqa: BLE001
        out.append(f"граничные условия: {type(e).__name__}: {e}")
    window = [(3, 2, Q(1, 2)), (1, 2, Q(-3, 2)), (5, 3, Q(2, 3)), (1, 7, Q(1, 4))]
    if deep:
        window.append((1, 2, Q(5, 6)))
    try:
        from .pi_value import PiMonomial
        from .exp_log import ExpValue
        window += [(1, PiMonomial(Q(1), Q(1)), Q(2)), (1, ExpValue.make(Q(1), Q(1)), Q(1, 2))]
        if deep:
            window.append((1, PiMonomial(Q(1), Q(1)), Q(3, 2)))
    except Exception as e:                                      # noqa: BLE001
        out.append(f"трансцендентные базы не построены: {type(e).__name__}: {e}")
    for m, a, e in window:
        try:
            v = PhiAlpha(m, a, e)
            out += [f"{v.record()}: {f}" for f in v.self_check(deep)]
        except Exception as ex:                                 # noqa: BLE001
            out.append(f"запись m={m}, α={a!r}, e={e}: {type(ex).__name__}: {ex}")
    return out
