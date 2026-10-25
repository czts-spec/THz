"""Профили ограниченной ширины — этап 2 плана.

Неограниченная разрядность несинтезируема: в кремний идёт не «вся
библиотека», а ограниченное ядро значения (`docs/PLAN.md`, INF-3). Профиль
называет три предела на Ф-запись Φ(s, α, M, E):

* `mantissa_bits` — разрядность числителя и знаменателя мантиссы M;
* `exp_den` — предел на знаменатель показателя E (дробные степени);
* `support` — число различных простых в мантиссе (носитель вектора
  показателей: то, что в этапе 3 станет разреженным массивом пар).

Два действия над значением и профилем, и их нельзя путать:

**`enforce` — рабочее.** Выход за профиль есть СОБЫТИЕ КОНТРАКТА
(`ProfileOverflow`), а не молчаливый перенос. Отказ называет поимённо,
что именно превышено и во сколько раз. Значение внутри профиля проходит
через `enforce` неизменным — тождественно, а не «почти».

**`project` — замерочное.** Ближайшее представимое в профиле значение.
Оно нужно ровно затем, чтобы посчитать, сколько различений профиль
теряет, и НЕ используется в рабочем пути. Округление знаменателя —
наилучшее рациональное приближение при ограниченном знаменателе
(цепные дроби), то есть подстановка не произвольная, а каноническая.

Отсюда критерий приёмки этапа 2: для каждого профиля публикуется
**таблица измеренного невежества** — сколько пар значений профиль
сливает, сколько ложно различает (у детерминированной подстановки — ноль)
и какая доля образцов из профиля выпадает событием контракта. Выбор
профиля становится решением «сколько различений мы согласны потерять», а
не «сколько бит взять».

Носитель считается разложением на простые пробным делением до `SUPPORT_LIMIT`;
если остаётся сомножитель крупнее, носитель не определён и это ОТКАЗ, а не
догадка: `Unresolved`.

Авторство: постановка (профиль, событие контракта, мера невежества) —
AIAΣA (автор проекта); наилучшее приближение при ограниченном знаменателе
— классические цепные дроби.
"""
from __future__ import annotations

from dataclasses import dataclass
from .errors import ExactRepresentationNotAvailable, Unresolved
from .rational import Q, limit_denominator
from .record import PhiRecord

SUPPORT_LIMIT = 1 << 20        # до какого простого идёт пробное деление


class ProfileOverflow(ExactRepresentationNotAvailable):
    """Событие контракта: значение не помещается в профиль.

    Несёт список превышений поимённо — что, сколько и каков предел.
    Молчаливого переноса нет: либо значение в профиле, либо событие.
    """

    def __init__(self, profile: "Profile", value: PhiRecord, breaches: list):
        self.profile, self.value, self.breaches = profile, value, list(breaches)
        detail = "; ".join(f"{name}: {got} > {limit}" for name, got, limit in breaches)
        super().__init__(f"профиль «{profile.name}» не вмещает {value.format()} — {detail}")


@dataclass(frozen=True)
class Profile:
    """Три предела ширины. Имя — не украшение: оно идёт в событие контракта."""

    name: str
    mantissa_bits: int
    exp_den: int
    support: int

    def __post_init__(self):
        for field in ("mantissa_bits", "exp_den", "support"):
            v = getattr(self, field)
            if isinstance(v, bool) or not isinstance(v, int) or v < 1:
                raise ValueError(f"{field} — целое ≥ 1, получено {v!r}")

    @property
    def max_mantissa(self) -> int:
        """Наибольшее целое, помещающееся в разрядность мантиссы."""
        return (1 << self.mantissa_bits) - 1

    def limits(self) -> dict:
        return {"мантисса, бит": self.mantissa_bits, "знаменатель E": self.exp_den,
                "носитель": self.support}


# Лестница профилей: от широкого к тому, что и вправду идёт в кремний.
LADDER = (
    Profile("широкий", mantissa_bits=64, exp_den=64, support=8),
    Profile("средний", mantissa_bits=32, exp_den=16, support=6),
    Profile("узкий", mantissa_bits=16, exp_den=8, support=4),
    Profile("кремний", mantissa_bits=8, exp_den=4, support=3),
)


def _prime_support(n: int) -> set:
    """Различные простые числа n; отказ, если разложение не дошло до конца."""
    n = abs(int(n))
    if n in (0, 1):
        return set()
    out, d = set(), 2
    while d * d <= n and d <= SUPPORT_LIMIT:
        if n % d == 0:
            out.add(d)
            while n % d == 0:
                n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        if n <= SUPPORT_LIMIT * SUPPORT_LIMIT:
            out.add(n)                       # остаток простой: d² > n
        else:
            raise Unresolved("носитель", str(n),
                             hint=f"пробное деление до {SUPPORT_LIMIT} не разложило "
                                  f"остаток {n}: носитель не определён")
    return out


def support_of(m: Q) -> int:
    """Носитель мантиссы: сколько различных простых в числителе и знаменателе."""
    return len(_prime_support(m.n) | _prime_support(m.d))


def breaches(rec: PhiRecord, profile: Profile) -> list:
    """Что именно выходит за профиль: (имя, замерено, предел). Пусто — влезает."""
    out = []
    for name, got in (("числитель M, бит", abs(rec.M.n).bit_length()),
                      ("знаменатель M, бит", rec.M.d.bit_length()),
                      ("числитель E, бит", abs(rec.E.n).bit_length())):
        if got > profile.mantissa_bits:
            out.append((name, got, profile.mantissa_bits))
    if rec.E.d > profile.exp_den:
        out.append(("знаменатель E", rec.E.d, profile.exp_den))
    s = support_of(rec.M)
    if s > profile.support:
        out.append(("носитель мантиссы", s, profile.support))
    return out


def fits(rec: PhiRecord, profile: Profile) -> bool:
    return not breaches(rec, profile)


def enforce(rec: PhiRecord, profile: Profile) -> PhiRecord:
    """Рабочий путь: значение внутри профиля — тождественно; иначе событие."""
    b = breaches(rec, profile)
    if b:
        raise ProfileOverflow(profile, rec, b)
    return rec


def _limit_den(q: Q, maxden: int) -> Q:
    """Наилучшее рациональное приближение при знаменателе ≤ maxden.

    Считается цепной дробью на парах пакета (`rational.limit_denominator`):
    у стандартной `Fraction` дверь для float открыта, и открывать её ради
    одного алгоритма незачем. Сверено с `Fraction.limit_denominator` на
    3 000 случайных пар — расхождений 0.
    """
    return limit_denominator(q, maxden)


def project(rec: PhiRecord, profile: Profile) -> PhiRecord:
    """ЗАМЕРОЧНОЕ: ближайшее представимое в профиле — или None, если и после
    приведения знаменателей значение в профиль не входит.

    Возвращает None, а не «что-нибудь»: сливать несливаемое — то же
    молчание, от которого профиль и защищает.
    """
    m = _limit_den(rec.M, profile.max_mantissa)
    e = _limit_den(rec.E, min(profile.exp_den, profile.max_mantissa))
    if m.n == 0:
        return None
    try:
        out = PhiRecord(rec.s, rec.alpha, m, e)
    except ValueError:
        return None
    return out if fits(out, profile) else None


def self_check(deep: bool = False) -> list:
    """Самопроверка профилей ВСТРЕЧНОЙ дорогой (ось самоприменения).

    `[вычислено здесь]` Дороги: (1) наилучшее приближение при ограниченном
    знаменателе считается ЦЕПНОЙ ДРОБЬЮ (`_limit_den`), а здесь сверяется
    ПЕРЕБОРОМ: для малого предела перебираются все знаменатели 1…maxden,
    для каждого берётся ближайший числитель, и лучший по модулю разности
    обязан совпасть с ответом цепной дроби (при равенстве расстояний
    сравниваются сами расстояния, а не записи — ничья законна); (2)
    `fits` и `breaches` — один предикат двумя выражениями: пусто ⟺ влезает;
    (3) `enforce` тождественен ровно там, где влезает, и поднимает
    `ProfileOverflow` ровно там, где нет — обе стороны, а не одна; (4)
    `project` либо честно отдаёт None, либо его ответ ВЛЕЗАЕТ в профиль
    (проверяется тем же `breaches`).

    Граница объявлена числом: перебор знаменателей идёт только при
    maxden ≤ 64 (при `deep` — ≤ 256) и на пробных парах — это цена, а не
    полнота; профили лестницы с большим пределом сверяются лишь дорогами
    (2)–(4). `project` есть ЗАМЕРОЧНЫЙ путь (он приближает), и его
    расхождение с исходной величиной находкой не считается — рабочий путь
    тождествен по построению. Возвращает список находок; исключений не
    бросает.
    """
    out = []
    probes = [PhiRecord(0, 2, Q(355, 113), Q(3)),
              PhiRecord(1, 2, Q(1, 3), Q(-2)),
              PhiRecord(0, 2, Q(2 ** 40 + 1, 3), Q(1, 7)),
              PhiRecord(0, 3, Q(7, 5), Q(1))]
    top = 256 if deep else 64
    for maxden in (2, 3, 7, 16, top):
        for rec in probes:
            for q in (rec.M, rec.E):
                got = _limit_den(q, maxden)
                if got.d > maxden:
                    out.append(f"знаменатель ответа {got.d} больше предела "
                               f"{maxden} на {q}")
                if maxden > 256:
                    continue
                best = None
                for d in range(1, maxden + 1):
                    k = (q.n * d * 2 + q.d) // (q.d * 2)      # ближайший числитель
                    for cand in (Q(k - 1, d), Q(k, d), Q(k + 1, d)):
                        dist = abs(cand - q)
                        if best is None or dist < best[0]:
                            best = (dist, cand)
                if best is not None and abs(got - q) != best[0]:
                    out.append(f"цепная дробь против перебора расходится на "
                               f"{q} при пределе {maxden}: {got} (расстояние "
                               f"{abs(got - q)}) против {best[1]} ({best[0]})")
    for prof in LADDER:
        for rec in probes:
            b = breaches(rec, prof)
            if bool(b) == fits(rec, prof):
                out.append(f"fits и breaches расходятся на {rec.format()} "
                           f"в профиле {prof.name}")
            try:
                same = enforce(rec, prof)
                if b:
                    out.append(f"enforce пропустил значение вне профиля "
                               f"{prof.name}: {rec.format()}")
                elif same is not rec:
                    out.append(f"enforce не тождествен внутри профиля "
                               f"{prof.name}: {same.format()}")
            except ProfileOverflow as e:
                if not b:
                    out.append(f"enforce поднял событие на влезающем значении "
                               f"{rec.format()} в профиле {prof.name}")
                elif list(e.breaches) != list(b):
                    out.append(f"событие назвало не те превышения: "
                               f"{e.breaches} против {b}")
            except Exception as e:                              # noqa: BLE001
                out.append(f"enforce отказал не своим событием: "
                           f"{type(e).__name__}: {e}")
            pr = project(rec, prof)
            if pr is not None and breaches(pr, prof):
                out.append(f"project отдал значение, не влезающее в профиль "
                           f"{prof.name}: {pr.format()}")
    return out


def resolution_table(samples, profiles=LADDER) -> list:
    """Таблица измеренного невежества: строка на профиль.

    Считается по парам образцов: слито — разные значения дали одно
    представление; ложно различено — равные значения дали разные (у
    канонической подстановки обязано быть нулём); выпало — образцы, для
    которых профиль поднимает событие контракта.
    """
    rows = []
    for p in profiles:
        proj = [project(s, p) for s in samples]
        dropped = sum(1 for x in proj if x is None)
        merged = false = pairs = 0
        for i in range(len(samples)):
            for j in range(i + 1, len(samples)):
                if proj[i] is None or proj[j] is None:
                    continue
                pairs += 1
                same_exact = (samples[i].s == samples[j].s
                              and samples[i].alpha == samples[j].alpha
                              and samples[i].M == samples[j].M
                              and samples[i].E == samples[j].E)
                same_repr = (proj[i].s == proj[j].s and proj[i].alpha == proj[j].alpha
                             and proj[i].M == proj[j].M and proj[i].E == proj[j].E)
                if same_repr and not same_exact:
                    merged += 1
                elif same_exact and not same_repr:
                    false += 1
        rows.append({"профиль": p.name, "пределы": p.limits(), "образцов": len(samples),
                     "выпало событием контракта": dropped, "пар сравнено": pairs,
                     "слито пар": merged, "ложно различено": false,
                     "в профиле как есть": sum(1 for s in samples if fits(s, p))})
    return rows
