"""Ноль как состояние с провенансом, а не единая точка.

Четыре рода нуля (ось происхождения):
  ABSOLUTE                 — мгновение определения (аксиома A1);
  STRUCTURAL_INFINITE      — от ω−ω, поглощает;
  STRUCTURAL_INFINITESIMAL — от ε−ε, уступает;
  NOT_COMPUTED             — «вычисление не производилось», 0∉S.

Ограничения:
  * никаких float; арифметика — только возврат операндов или отказ;
  * все значения в сообщениях — явной Ф-записью (record());
  * провенанс знаком не флипается: __neg__ возвращает self.
"""

from enum import Enum

from .errors import Unresolved


class ZeroKind(Enum):
    """Род нуля — провенанс, не величина."""

    #: Мгновение определения (аксиома A1). Нейтрален и аннигилирует
    #: в умножении: 0·x = 0, включая 0·ω = 0 (правило корпуса
    #: «ε×0 = обычный ноль»).
    ABSOLUTE = "absolute"

    #: От ω−ω: неограниченная неопределённость. ПОГЛОЩАЕТ операнды
    #: (решение линии v6, подтверждено автором). Предмет решения
    #: автора №2 (окончательный статус поглощения) остаётся
    #: ОТКРЫТЫМ — здесь зафиксировано текущее правило, не итог.
    STRUCTURAL_INFINITE = "structural_infinite"

    #: От ε−ε: идеал бесконечно малых. УСТУПАЕТ в сложении
    #: (нейтрален), в умножении сохраняет свой род.
    STRUCTURAL_INFINITESIMAL = "structural_infinitesimal"

    #: «Вычисление не производилось»: 0∉S, Next(N)=1 (OP-ш Theory).
    #: Аргумент корпуса: часы — 300 лет без переполнения против
    #: счётчиков с нулём. Любая арифметика — отказ Unresolved.
    NOT_COMPUTED = "not_computed"


def _phi(x: object) -> str:
    """Ф-запись операнда для сообщений отказа."""
    rec = getattr(x, "record", None)
    return rec() if callable(rec) else str(x)


# Старшинство родов при встрече двух нулей (свежие глаза 05.09; первые
# принципы: A1 — абсолютный ноль аннигилирует всё, включая ω, потому в
# произведении он старше любого структурного; в сумме поглощающий
# структурный ноль (ω−ω) старше нейтральных). NOT_COMPUTED отвергнут раньше.
_PRODUCT_RANK = {"ABSOLUTE": 3, "STRUCTURAL_INFINITE": 2, "STRUCTURAL_INFINITESIMAL": 1}
_SUM_RANK = {"STRUCTURAL_INFINITE": 3, "ABSOLUTE": 2, "STRUCTURAL_INFINITESIMAL": 1}


def _dominant_product(a: "Zero", b: "Zero") -> "Zero":
    return a if _PRODUCT_RANK[a.kind.name] >= _PRODUCT_RANK[b.kind.name] else b


def _dominant_sum(a: "Zero", b: "Zero") -> "Zero":
    return a if _SUM_RANK[a.kind.name] >= _SUM_RANK[b.kind.name] else b


class Zero:
    """Ноль-состояние. Хранит только род (kind); величины не несёт."""

    __slots__ = ("kind",)

    def __init__(self, kind: ZeroKind):
        if not isinstance(kind, ZeroKind):
            raise TypeError("kind должен быть ZeroKind")
        self.kind = kind

    # --- конструкторы ---------------------------------------------------

    @classmethod
    def absolute(cls) -> "Zero":
        return cls(ZeroKind.ABSOLUTE)

    @classmethod
    def structural_infinite(cls) -> "Zero":
        return cls(ZeroKind.STRUCTURAL_INFINITE)

    @classmethod
    def structural_infinitesimal(cls) -> "Zero":
        return cls(ZeroKind.STRUCTURAL_INFINITESIMAL)

    @classmethod
    def not_computed(cls) -> "Zero":
        return cls(ZeroKind.NOT_COMPUTED)

    # --- свойства -------------------------------------------------------

    def absorbs(self) -> bool:
        """Поглощающий ли это ноль (ω−ω)."""
        return self.kind is ZeroKind.STRUCTURAL_INFINITE

    def _refuse(self, op: str, other: object, self_left: bool) -> None:
        left, right = (self.record(), _phi(other)) if self_left \
            else (_phi(other), self.record())
        raise Unresolved(op, left, right,
                         hint="ноль-состояние: вычисление не производилось")

    # --- сложение -------------------------------------------------------

    def _guard(self, op: str, other: object, self_left: bool) -> None:
        """Общие запреты на оба порядка операндов (свежие глаза 05.09:
        порядок операндов давал разные итоги). float — TypeError
        построением; NOT_COMPUTED с ЛЮБОЙ стороны — отказ."""
        if isinstance(other, float):
            raise TypeError("float запрещён построением: ноль не бывает «примерно»")
        if self.kind is ZeroKind.NOT_COMPUTED:
            self._refuse(op, other, self_left)
        if isinstance(other, Zero) and other.kind is ZeroKind.NOT_COMPUTED:
            other._refuse(op, self, not self_left)

    def __add__(self, other: object) -> object:
        self._guard("+", other, True)
        if isinstance(other, Zero):
            return _dominant_sum(self, other)
        if self.absorbs():
            return self
        # ABSOLUTE, STRUCTURAL_INFINITESIMAL: нейтральный элемент для
        # ВСЯКОГО x, включая указатели: ε+0 = ε (исправление замеренного
        # дефекта ε+0 → 0).
        return other

    def __radd__(self, other: object) -> object:
        self._guard("+", other, False)
        if self.absorbs():
            return self
        return other

    # --- вычитание ------------------------------------------------------

    def __sub__(self, other: object) -> object:
        """0 − x = −x через other.__neg__(); иначе — отказ."""
        self._guard("-", other, True)
        if isinstance(other, Zero):
            return _dominant_sum(self, other)
        if self.absorbs():
            return self
        neg = getattr(other, "__neg__", None)
        if callable(neg):
            return neg()
        raise Unresolved("-", self.record(), _phi(other),
                         hint="операнд не даёт __neg__")

    def __rsub__(self, other: object) -> object:
        """x − 0 = x; поглощающий ноль — self; NOT_COMPUTED — отказ."""
        self._guard("-", other, False)
        if self.absorbs():
            return self
        return other

    # --- умножение ------------------------------------------------------

    def __mul__(self, other: object) -> "Zero":
        self._guard("*", other, True)
        # ABSOLUTE: 0·x = 0, и 0·ω = 0 (правило корпуса ε×0 = обычный
        # ноль). STRUCTURAL_*: результат — self с СОХРАНЕНИЕМ kind.
        # Два нуля — по старшинству родов, не по порядку операндов
        # (свежие глаза 05.09: ABSOLUTE·(ω−ω) и (ω−ω)·ABSOLUTE расходились).
        if isinstance(other, Zero):
            return _dominant_product(self, other)
        return self

    def __rmul__(self, other: object) -> "Zero":
        self._guard("*", other, False)
        if isinstance(other, Zero):
            return _dominant_product(other, self)
        return self

    # --- знак и сравнение -----------------------------------------------

    def __neg__(self) -> "Zero":
        """Провенанс — компонент записи, знаком не флипается."""
        return self

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Zero):
            return self.kind is other.kind
        return NotImplemented

    def __hash__(self) -> int:
        return hash((Zero, self.kind))

    def __bool__(self) -> bool:
        return False

    # --- самопроверка (ось «AIAΣA спрашивает AIAΣA», 04.09) --------------

    def self_check(self, deep: bool = False) -> list:
        """Самопроверка ноль-состояния ВСТРЕЧНОЙ дорогой.

        `[вычислено здесь]` Ноль-состояние величины не несёт, поэтому
        встречная дорога у него одна и она внешняя: те же четыре рода живут
        в СОСЕДНЕМ модуле другим устройством — контейнер v6 держит их
        классами 0 и 20 с провенансом и родом (`void`, `chaos`,
        `center_zero(structural, kind)`), и правило «уступает или
        поглощает» там выведено из своего канона, не из этого. Сверяется
        именно правило: absolute и structural_infinitesimal обязаны
        УСТУПИТЬ числу (x + 0 = x), structural_infinite — ПОГЛОТИТЬ, и то же
        самое обязан дать контейнер на своём роде. Вторая дорога — запись
        (`record`) против рода: род назван в записи, и классы 0/20 в ней
        различены. Третья — контракт eq/hash и bool(ноль) = False.

        Граница: род NOT_COMPUTED арифметики не имеет вовсе (любая —
        отказ `Unresolved`), и его сверка идёт ловлей отказа, а не
        сравнением величин; у контейнера прямого двойника этому роду нет —
        расхождение здесь объявлено, а не считается находкой. Возвращает
        список находок; исключений не бросает.
        """
        out = []
        if not isinstance(self.kind, ZeroKind):
            out.append(f"род не ZeroKind: {self.kind!r}")
            return out
        if self.absorbs() != (self.kind is ZeroKind.STRUCTURAL_INFINITE):
            out.append(f"absorbs() = {self.absorbs()} не отвечает роду "
                       f"{self.kind.name}")
        rec = self.record()
        if self.kind.value not in rec and not (
                self.kind is ZeroKind.ABSOLUTE and "kind=absolute" in rec):
            out.append(f"род не назван в записи: {rec}")
        if (self.kind is ZeroKind.ABSOLUTE) != rec.startswith("Φ(C=0 ZERO"):
            out.append(f"класс в записи не отвечает роду: {rec}")
        if bool(self):
            out.append("ноль-состояние истинно в булевом чтении")
        twin = Zero(self.kind)
        if twin != self or hash(twin) != hash(self):
            out.append("контракт eq/hash нарушен: равные роды — разные ключи")
        if self.kind is ZeroKind.NOT_COMPUTED:
            try:
                self + 5
                out.append("род not_computed допустил сложение вместо отказа")
            except Unresolved:
                pass
            except Exception as e:                              # noqa: BLE001
                out.append(f"род not_computed отказал не Ф-отказом: "
                           f"{type(e).__name__}")
            return out
        try:
            yields_here = (self + 5) == 5
        except Exception as e:                                  # noqa: BLE001
            out.append(f"сложение с числом отказало: {type(e).__name__}: {e}")
            return out
        if yields_here == self.absorbs():
            out.append(f"правило сложения против рода: absorbs={self.absorbs()}, "
                       f"а 0 + 5 {'уступил' if yields_here else 'поглотил'}")
        if not deep:
            return out
        try:
            from .container import PhiNumberV6 as _P
            twin_c = {ZeroKind.ABSOLUTE: _P.center_zero(),
                      ZeroKind.STRUCTURAL_INFINITE: _P.center_zero(True, 'infinite'),
                      ZeroKind.STRUCTURAL_INFINITESIMAL: _P.center_zero(True, 'infinitesimal'),
                      }[self.kind]
            there = twin_c + _P.from_int(5)
            yields_there = not there.is_uncomputable
            if yields_there != yields_here:
                out.append(f"сосед-контейнер решает иначе: здесь "
                           f"{'уступает' if yields_here else 'поглощает'}, "
                           f"там {'уступает' if yields_there else 'поглощает'}")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"сверка с контейнером отказала: {type(e).__name__}: {e}")
        return out

    # --- запись ---------------------------------------------------------

    def record(self) -> str:
        """Явная Ф-запись: C=0 ZERO для абсолютного, C=20 UNCOMPUTABLE
        для структурных и невычисленного; kind — провенанс."""
        if self.kind is ZeroKind.ABSOLUTE:
            return "Φ(C=0 ZERO, kind=absolute)"
        return f"Φ(C=20 UNCOMPUTABLE, kind={self.kind.value})"

    def __repr__(self) -> str:
        return self.record()
