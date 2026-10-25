"""Логарифм сложения — точный носитель показателя над основанием α.

Известное (Леонелли 1803, Гаусс 1812 — «аддитивные логарифмы»; Цех —
то же над конечным полем; LNS 1975 — то же в железе):

    α^a + α^b = α^(a + Z(b−a)),   Z(x) = log_α(1 + α^x).

Слово автора 02.09: «делаем ту часть, которую можем и знаем как… нельзя
идти дальше, пока известное не является твоим аппаратом». Здесь
известное — аппарат: показатель записи Φα = m·α^E несёт величины
класса 6 над своим основанием,

    E = q + Σ cᵢ · log_α(1 + kᵢ·α^{dᵢ}),   q, cᵢ, kᵢ — пары; dᵢ — показатель.

Z хранится СИМВОЛОМ (в LNS-железе его приближают таблицей — здесь ни
одного приближения): сумма над одной базой остаётся одной записью над
той же базой, а умножение на оси E — по-прежнему сложение. Перевод в
носителя — по явному запросу и только там, где носитель есть.
"""
from __future__ import annotations

from .rational import Q


def _exp_scale(d, s: Q):
    """Показатель при смене основания α = β^s: α^d = β^(s·d), а log_β = s·log_α,
    потому показатель d над α становится s·d над β (для ZechExponent — rebase:
    рациональная часть и показатели внутри логарифмов масштабируются на s)."""
    return d.rebase(s) if isinstance(d, ZechExponent) else Q(d) * s


class Zech:
    """Z = log_α(1 + k·α^d), k ≠ 0. Основание α принадлежит записи."""

    __slots__ = ("k", "d")

    def __init__(self, k, d):
        k = Q(k)
        if k == 0:
            raise ValueError("логарифм сложения с k=0 — это ноль, не член")
        # инвариант 1 + k·α^d > 0 (очередь п. 48, 10.09): при d = 0 он не зависит от α —
        # 1 + k ≤ 0 значит аргумент логарифма ≤ 0; при d ≠ 0 знак проверяется там, где
        # известна база (alpha.py, _sign_of_one_plus при присоединении показателя)
        if not isinstance(d, ZechExponent) and Q(d) == 0 and k <= -1:
            raise ValueError(f"log_α(1 + {k}·α^0): аргумент логарифма 1 + {k} ≤ 0 — не число")
        if isinstance(d, ZechExponent) and d.is_rational and d.q == 0 and k <= -1:
            raise ValueError(f"log_α(1 + {k}·α^0): аргумент логарифма 1 + {k} ≤ 0 — не число")
        object.__setattr__(self, "k", k)
        object.__setattr__(self, "d", d if isinstance(d, ZechExponent)
                           else Q(d))

    def __setattr__(self, *_):
        raise AttributeError("Zech неизменяем")

    def self_check(self, deep: bool = False) -> list:
        """Самопроверка показателя Цеха независимой дорогой (очередь п. 18).

        `[вычислено здесь]` Проверяется не тем счётом, каким строилось:
        коэффициент не нулевой, аргумент логарифма положителен там, где он от
        базы не зависит (d = 0), перенос базы обратим (`rebase(s).rebase(1/s)`
        возвращает исходное), запись читается, хеш согласован с равенством.

        Граница: при d ≠ 0 знак `1 + k·α^d` зависит от базы и здесь не
        проверяется — это место знает `alpha.py`, и там оно проверено.
        """
        out = []
        if self.k == 0:
            out.append("коэффициент нулевой: член не член")
        d0 = (not isinstance(self.d, ZechExponent) and Q(self.d) == 0) or \
             (isinstance(self.d, ZechExponent) and self.d.is_rational and self.d.q == 0)
        if d0 and self.k <= -1:
            out.append(f"аргумент логарифма 1 + {self.k} ≤ 0 при d = 0")
        if deep:
            for s in (Q(2), Q(1, 3)):
                try:
                    back = self.rebase(s).rebase(Q(1) / s)
                except Exception as e:                      # noqa: BLE001
                    out.append(f"перенос базы отказал на s = {s}: {type(e).__name__}")
                    continue
                if back != self:
                    out.append(f"перенос базы необратим при s = {s}: {back.record()}")
            if hash(self) != hash(Zech(self.k, self.d)):
                out.append("хеш не согласован с равенством")
            if "log_α(1 + " not in self.record():
                out.append("запись не читается как логарифм суммы")
        return out

    def rebase(self, s: Q) -> "Zech":
        """Смена основания α = β^s: log_α(1 + k·α^d) = log_β(1 + k·β^{s·d}) / s
        — показатель d масштабируется, k не меняется.

        Основание здесь ТРАНСЦЕНДЕНТНОЕ (база логарифма сложения), и смена —
        точное масштабирование без отказов. Одноимённая дверь
        `wheels.Record.rebase` — то же понятие (та же величина при другом
        основании), но основание ПОЗИЦИОННОЕ b ∈ ℕ: там меняется сетка цифр,
        ранг колёс может стать ω, и дверь отказывает по имени. Родство —
        величина; различие — что есть «основание».
        """
        return Zech(self.k, _exp_scale(self.d, s))

    def __eq__(self, other):
        return (isinstance(other, Zech) and self.k == other.k
                and self.d == other.d)

    def __hash__(self):
        return hash(("Z", self.k, self.d))

    def record(self) -> str:
        d = self.d.record() if isinstance(self.d, ZechExponent) \
            else repr(self.d)
        return f"log_α(1 + {self.k}·α^({d}))"

    __repr__ = record


class ZechExponent:
    """E = q + Σ cᵢ·Zᵢ — показатель с логарифмами сложения."""

    __slots__ = ("q", "terms")

    def __init__(self, q=0, terms=()):
        merged: dict = {}
        order = []
        for c, z in terms:
            c = Q(c)
            if z not in merged:
                merged[z] = Q(0)
                order.append(z)
            merged[z] = merged[z] + c
        kept = tuple((merged[z], z) for z in order if merged[z] != 0)
        object.__setattr__(self, "q", Q(q))
        object.__setattr__(self, "terms", kept)

    def __setattr__(self, *_):
        raise AttributeError("ZechExponent неизменяем")

    @property
    def is_rational(self) -> bool:
        return not self.terms

    # --- арифметика показателей: сложение, знак, масштаб, смена базы ------
    def __add__(self, other):
        if isinstance(other, ZechExponent):
            return ZechExponent(self.q + other.q, self.terms + other.terms)
        if isinstance(other, (int, Q)) and not isinstance(other, bool):
            return ZechExponent(self.q + Q(other), self.terms)
        return NotImplemented

    __radd__ = __add__

    def __neg__(self):
        return ZechExponent(-self.q, tuple((-c, z) for c, z in self.terms))

    def __sub__(self, other):
        r = self + (-other if isinstance(other, ZechExponent)
                    else -Q(other))
        return r

    def __rsub__(self, other):
        return (-self) + Q(other)

    def scale(self, s) -> "ZechExponent":
        """(α^E)^s = α^(s·E): чистое масштабирование показателя."""
        s = Q(s)
        return ZechExponent(self.q * s, tuple((c * s, z) for c, z in self.terms))

    def rebase(self, s) -> "ZechExponent":
        """Смена основания α = β^s: E_β = s·q + Σ cᵢ·log_β(1 + kᵢ·β^{s·dᵢ}).

        Делегирует каждому члену `Zech.rebase` (одно понятие, одна дверь);
        от `wheels.Record.rebase` отличается тем же, чем и та: основание
        трансцендентное, не позиционное.
        """
        s = Q(s)
        return ZechExponent(self.q * s,
                            tuple((c, z.rebase(s)) for c, z in self.terms))

    # --- сравнение: структурное; с парой — только рациональный случай ----
    def __eq__(self, other):
        if isinstance(other, ZechExponent):
            return self.q == other.q and set(self.terms) == set(other.terms)
        if isinstance(other, (int, Q)) and not isinstance(other, bool):
            return self.is_rational and self.q == Q(other)
        return NotImplemented

    def __hash__(self):
        if not self.terms:                      # контракт eq/hash: E == q ⇒ hash(E) == hash(q)
            return hash(self.q)
        return hash(("E", self.q, frozenset(self.terms)))

    def record(self) -> str:
        parts = [repr(self.q)] if (self.q != 0 or not self.terms) else []
        for c, z in self.terms:
            parts.append(f"{c}·{z.record()}" if c != 1 else z.record())
        return " + ".join(parts)

    __repr__ = record


def exp_add(a, b):
    """Сумма показателей (пара или ZechExponent) — пара, если рациональна."""
    if isinstance(a, ZechExponent) or isinstance(b, ZechExponent):
        r = (a if isinstance(a, ZechExponent) else ZechExponent(a)) + b
        return r.q if r.is_rational else r
    return Q(a) + Q(b)


def exp_neg(a):
    return -a if isinstance(a, ZechExponent) else -Q(a)


def exp_scale(a, s):
    if isinstance(a, ZechExponent):
        r = a.scale(s)
        return r.q if r.is_rational else r
    return Q(a) * Q(s)


def exp_rebase(a, s):
    if isinstance(a, ZechExponent):
        r = a.rebase(s)
        return r.q if r.is_rational else r
    return Q(a) * Q(s)
