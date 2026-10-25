"""Явная Ф-запись контейнера: значение ≡ (−1)^s · M · α^E.

Ограничения:
  * никаких float — только int и Q;
  * M > 0; ноль в этой записи запрещён — ноль живёт в zero_state;
  * все значения в сообщениях — явной Ф-записью (format()/record());
  * перевод (to_fraction, translate) — отдельное действие; при
    невозможности точного результата — честный отказ
    ExactRepresentationNotAvailable, не скрытое приближение.
"""

import re
from dataclasses import dataclass
from .rational import Q

from .errors import ExactRepresentationNotAvailable
# Колесо — не сосед записи, а её закон (слово автора 09.09): позиционная
# строка есть СОКРАЩЁННАЯ ЗАПИСЬ СУММЫ КОЛЁС, ноль — отсутствие колеса.
# Ввоз верхнего уровня безопасен: wheels → container → record только
# внутри функции (container.py:1611), цикла при загрузке нет.
from .wheels import Record as _WheelRecord, _DIGITS, _check_base

#: Разбор ровно формата format(); пробелы гибко, включая вокруг «/».
_PARSE_RE = re.compile(
    r"\s*Φ\(\s*s\s*=\s*(\d+)\s*,"
    r"\s*α\s*=\s*(\d+)\s*,"
    r"\s*M\s*=\s*(\d+(?:\s*/\s*\d+)?)\s*,"
    r"\s*E\s*=\s*(-?\d+(?:\s*/\s*\d+)?)\s*\)\s*"
)


@dataclass(frozen=True)
class PhiRecord:
    """Ф-запись: s ∈ {0,1}, α ≥ 2, M > 0 (Q), E — Q."""

    s: int
    alpha: int
    M: Q
    E: Q

    def __post_init__(self):
        if isinstance(self.s, bool) or not isinstance(self.s, int) \
                or self.s not in (0, 1):
            raise ValueError(f"s должен быть 0 или 1, получено {self.s!r}")
        if isinstance(self.alpha, bool) or not isinstance(self.alpha, int) \
                or self.alpha < 2:
            raise ValueError(f"α должен быть целым ≥ 2, получено "
                             f"{self.alpha!r}")
        for name in ("M", "E"):
            v = getattr(self, name)
            if isinstance(v, float):
                raise TypeError(f"{name}: float запрещён, только Q")
            if not isinstance(v, Q):
                # int допустим — нормализуем в Q (frozen: setattr).
                object.__setattr__(self, name, Q(v))
        if self.M <= 0:
            raise ValueError(
                "M должен быть > 0; ноль в Ф-записи запрещён — "
                "ноль живёт в zero_state (Zero)")

    # --- запись ---------------------------------------------------------

    def format(self) -> str:
        """Явная Ф-запись; единственная форма значения в сообщениях."""
        return f"Φ(s={self.s}, α={self.alpha}, M={self.M}, E={self.E})"

    def record(self) -> str:
        """Синоним format() — общий протокол записи пакета.

        Носители одного класса (нули, указатели, алгебраические) печатают
        Φ(C=<класс>, …); эта запись — старшая форма Φ(s, α, M, E), без
        поля C: знак вынесен отдельным полем s.
        Контейнер v6 в протоколе УЧАСТВУЕТ с 12.09 (п. 5 главы 09, дыра
        была названа прогоном главы): у PhiNumberV6 есть record(), и
        печать с записью у него один текст (`__repr__ = record`,
        `container.py`, канон соседей Φα, Notch, Zech, Ieee) — свой
        старший знак «Φv6(…)» с классом внутри числом и именем:
        `Φv6(3/2 · 2^(0/1), C=1 NATURAL_POWER)`. Прежде здесь стояло
        обратное («у PhiNumberV6 нет record(), он печатается repr»), и
        следствий было два, оба замерены: фасад ФIOS печатал контейнер
        запасной ветвью `Φ(?=…)` с приближением внутри, а
        `from_ratio(3,2) ** from_int(2)` давал голый TypeError вместо
        названного отказа.
        Две записи одной величины остаются двумя записями: `to_container`
        ниже переводит Φ(s, α, M, E) в контейнер, `from_container` —
        обратно, и ровно там названо, что сводится, а что нет.
        """
        return self.format()

    def __str__(self) -> str:
        return self.format()

    def __repr__(self) -> str:
        return self.format()

    def self_check(self, deep: bool = False) -> list:
        """Самопроверка записи ВСТРЕЧНОЙ дорогой (ось самоприменения).

        `[вычислено здесь]` Дороги здесь три, и ни одна не повторяет ту, что
        строила запись: (1) разбор своей же печати — `parse(format())`
        обязан дать ту же запись (регулярное выражение против сборки
        строки); (2) мост в СОСЕДНИЙ модуль и обратно —
        `from_container(to_container())` (контейнер v6 хранит знак в m1, а
        не полем s, и обратный перевод обязан вернуть исходную четвёрку);
        (3) величина двумя дорогами — `to_fraction()` (умножение пары на
        α^E) против `to_algebraic()` (корень степени q из M^q·α^p на
        линейке 3, `algebraic.py`), которые при целом E обязаны сойтись.

        Граница названа: при НЕЦЕЛОМ E пары целых не существует, и дорога
        (3) идёт только алгебраической стороной — расхождение там означало
        бы не ошибку записи, а отказ перевода, и оно не считается находкой;
        при α, которое контейнер не сводит к целому основанию, мост (2)
        честно отказывает — это тоже объявлено, а не находка. Сломанная
        запись возвращает находки, а не бросает исключение.
        """
        out = []
        if self.s not in (0, 1):
            out.append(f"знак вне {{0, 1}}: {self.s!r}")
        if not isinstance(self.alpha, int) or self.alpha < 2:
            out.append(f"основание не целое ≥ 2: {self.alpha!r}")
        if not isinstance(self.M, Q) or self.M <= 0:
            out.append(f"M не положительная пара: {self.M!r}")
        if not isinstance(self.E, Q):
            out.append(f"E не пара: {self.E!r}")
        if out:
            return out
        try:
            back = PhiRecord.parse(self.format())
            if back != self:
                out.append(f"разбор своей печати дал иное: {back.format()}")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"своя печать не разбирается: {type(e).__name__}: {e}")
        try:
            again = PhiRecord.from_container(self.to_container())
            if again != self:
                out.append(f"мост в контейнер и обратно дал иное: {again.format()}")
        except ExactRepresentationNotAvailable:
            pass                       # объявленная граница: α вне сводимых
        except Exception as e:                                  # noqa: BLE001
            out.append(f"мост в контейнер отказал не по границе: "
                       f"{type(e).__name__}: {e}")
        if self.E.denominator == 1:
            try:
                f, alg = self.to_fraction(), self.to_algebraic()
                if alg != f:
                    out.append(f"величина двумя дорогами расходится: "
                               f"пара {f}, линейка 3 {alg!r}")
            except Exception as e:                              # noqa: BLE001
                out.append(f"перевод величины отказал: {type(e).__name__}: {e}")
        if not deep:
            return out
        if self.E.denominator == 1:
            # третья встречная дорога: десятичная сторона собирается обратно
            # в пару целыми (усечение к нулю обязано лишь уменьшать модуль)
            try:
                digits = 6
                val = self.to_fraction()
                # Обратная дорога идёт КОЛЕСОМ и по трём сеткам сразу: запись
                # позиционна в любом основании, и свойство усечения — свойство
                # ПАРЫ (число, основание), а не десятки (слово автора 09.09).
                for base in (10, 2, 16):
                    text, exact = self.translate(digits, base)
                    back_q = _WheelRecord.parse(text, base, sep=".").value
                    if exact and back_q != val:
                        out.append(f"позиционная сторона (основание {base}) "
                                   f"объявлена точной, а расходится: {text} "
                                   f"против {val}")
                    if abs(back_q) > abs(val):
                        out.append(f"усечение к нулю дало модуль БОЛЬШЕ "
                                   f"исходного (основание {base}): {text} "
                                   f"против {val}")
                    if not exact and abs(back_q) == abs(val):
                        out.append(f"усечение объявлено неточным, а совпало "
                                   f"(основание {base}): {text}")
            except Exception as e:                              # noqa: BLE001
                out.append(f"позиционный перевод отказал: {type(e).__name__}: {e}")
        return out

    @classmethod
    def parse(cls, text: str) -> "PhiRecord":
        """Разбор ровно формата format(); пробелы гибко.

        Roundtrip обязана держаться: parse(r.format()) == r.
        """
        m = _PARSE_RE.fullmatch(text)
        if m is None:
            raise ValueError(f"не Ф-запись: {text!r}; ожидается формат "
                             f"Φ(s=..., α=..., M=..., E=...)")
        s_txt, a_txt, m_txt, e_txt = m.groups()
        return cls(int(s_txt), int(a_txt),
                   Q(re.sub(r"\s+", "", m_txt)),
                   Q(re.sub(r"\s+", "", e_txt)))

    # --- переводы (отдельные действия) ----------------------------------

    def to_fraction(self) -> Q:
        """Точный перевод в Q; только при целом E (E<0 допустимо)."""
        if self.E.denominator != 1:
            raise ExactRepresentationNotAvailable(
                f"перевод в Q невозможен: E не целое у "
                f"{self.format()}")
        value = self.M * Q(self.alpha) ** self.E.numerator
        return -value if self.s else value

    def to_algebraic(self):
        """Перевод в алгебраическое: для E=p/q значение есть корень
        q-й степени из M^q·α^p, со знаком s.

        AlgebraicNumber (Φ(C=3 POLY_ROOT, …)) выходит лишь когда радиканд
        не есть точная q-я степень в Q. При целом E (q=1) и вообще при
        рациональном корне from_root распознаёт Q и возвращает Q со
        знаком — то же, что to_fraction(): Φ(s=1, α=2, M=3/4, E=5) → −24,
        α=4, M=1, E=1/2 → 2. Это ответ, не отказ: Q ⊂ алгебраические.
        """
        # Ленивый импорт — от цикла record ↔ algebraic.
        from .algebraic import AlgebraicNumber
        p, q = self.E.numerator, self.E.denominator
        if p >= 0:
            radicand = self.M ** q * Q(self.alpha) ** p
        else:
            radicand = self.M ** q / Q(self.alpha) ** (-p)
        return AlgebraicNumber.from_root(q, radicand, sign=self.s)

    def translate(self, digits: int, base: int = 10) -> tuple[str, bool]:
        """Позиционная строка на digits знаков после запятой и флаг exact.

        Основание — ПАРАМЕТР сетки, не константа; по умолчанию 10, и при
        основании 10 ответ побайтно тот же, что давала прежняя форма
        (приёмка: proto/test_phicore_record_wheels.py, 4 000 пар).

        Слово автора 09.09 [в корпусе]: позиционное число есть сокращённая
        запись СУММЫ КОЛЁС, где 0 — ОТСУТСТВИЕ колеса. Строка потому не
        считается здесь заново, а собирается из записи колёсами
        (`wheels.Record.of`): цифра на месте e — мантисса колеса номер e, а
        напечатанный 0 — ровно то, что колеса на этом месте нет. Прежняя
        форма зашивала десятку (`rem *= 10`) и отдавала целую часть
        питоновым `str(int_part)` — то есть отрезала запись от закона дома
        и от `phicore.wheels` (замер 22.09: ввозов wheels из 54 модулей
        phicore — 0; строка 90 очереди).

        Усечённая величина лежит в ℤ[1/base] построением (знаменатель —
        степень основания), потому `Record.of` здесь не отказывает: закон
        ранга W_b(x) < ∞ ⟺ x ∈ ℤ[1/b] (wheels.py) выполнен заведомо.

        Только для целого E; перевод корней в этой версии не реализован
        и честно отказывает. Усечение к нулю; усечение ВСЕГДА объявляется
        флагом exact (остаток нулевой).
        """
        if isinstance(digits, bool) or not isinstance(digits, int) \
                or digits < 0:
            raise ValueError(f"digits должен быть целым ≥ 0, получено "
                             f"{digits!r}")
        base = _check_base(base)
        if self.E.denominator != 1:
            raise ExactRepresentationNotAvailable(
                f"десятичный перевод не реализован для нецелого E у "
                f"{self.format()}")
        f = self.to_fraction()
        scaled = abs(f) * base ** digits
        whole = scaled.numerator // scaled.denominator      # усечение к нулю
        exact = scaled.denominator == 1
        колёса = _WheelRecord.of(Q(whole, base ** digits), base,
                                 width=1).counts()
        старшее = max(max(колёса) if колёса else 0, 0)
        text = ("-" if self.s else "") + "".join(
            _DIGITS[колёса.get(e, 0)] for e in range(старшее, -1, -1))
        if digits:
            text += "." + "".join(_DIGITS[колёса.get(e, 0)]
                                  for e in range(-1, -digits - 1, -1))
        return text, exact

    @classmethod
    def from_fraction(cls, f: Q) -> "PhiRecord":
        """Q → Ф-запись: s по знаку, α=2, M=|f|, E=0.

        f == 0 запрещён — ноль живёт в zero_state (Zero).
        """
        if isinstance(f, float):
            raise TypeError("float запрещён, только Q")
        f = Q(f)
        return cls(1 if f < 0 else 0, 2, abs(f), Q(0))

    # --- мост к контейнеру v6 (миграция, этап E): знак — в m1 -----------

    def to_container(self):
        """Ф-запись → контейнер v6 без потерь: M в пару мантиссы (знак в
        m1), E в пару экспоненты; α=2 — неявная база, целое α>2 — база
        NATURAL_POWER контейнера."""
        # Ленивый импорт — record ниже контейнера в порядке сборки пакета.
        from .container import NumberClass, PhiNumberV6
        m_num = -self.M.numerator if self.s else self.M.numerator
        if self.alpha == 2:
            return PhiNumberV6.from_mantissa_exponent(
                m_num, self.M.denominator,
                self.E.numerator, self.E.denominator)
        return PhiNumberV6.exotic_base(
            m_num, self.M.denominator,
            self.E.numerator, self.E.denominator,
            NumberClass.NATURAL_POWER, {"value": self.alpha},
            NumberClass.NATURAL_POWER)

    @classmethod
    def from_container(cls, c) -> "PhiRecord":
        """Контейнер v6 → Ф-запись; сводимы только обычные числа
        (NATURAL_POWER, ALGEBRAIC_ROOT) с ненулевой мантиссой и неявной
        базой 2 либо целой базой NATURAL_POWER. Прочее честно отказывает
        ExactRepresentationNotAvailable, называя ЧТО не сведено: класс
        (ZERO, HYPER/INFINITY, IGNORANCE, иные классы), нулевую мантиссу
        обычного класса (from_int(0)) или основание. Не контейнер —
        TypeError."""
        from .container import NumberClass, PhiNumberV6
        if not isinstance(c, PhiNumberV6):
            raise TypeError("ожидается PhiNumberV6")
        if c.number_class not in (NumberClass.NATURAL_POWER,
                                  NumberClass.ALGEBRAIC_ROOT):
            raise ExactRepresentationNotAvailable(
                f"класс {c.number_class.name} контейнера не сводится к "
                f"записи Φ(s, α, M, E) с M > 0")
        if c.m1.v == 0:
            raise ExactRepresentationNotAvailable(
                f"нулевая мантисса контейнера класса {c.number_class.name} "
                f"не сводится к записи Φ(s, α, M, E) с M > 0 — ноль живёт "
                f"в zero_state (Zero)")
        if c.base_class is None:
            alpha = 2
        elif c.base_class == NumberClass.NATURAL_POWER:
            alpha = int(c.base_params.get("value", 2))
        else:
            raise ExactRepresentationNotAvailable(
                f"основание {c.base_class.name} не сводится к целому α")
        return cls(1 if c.m1.v < 0 else 0, alpha,
                   Q(abs(c.m1.v), c.m2.v),
                   Q(c.e1.v, c.e2.v))
