"""Оси ОП-геометрии: дискретный тор без нуля, шов и три нуля оси —
0{ОП+0}, 0{ОП−0}, 0{ОП0}. Правила выведены вычислением из первых
принципов и корпуса, не ожиданием слова автора (слово автора 04.09,
пп. 2–3: «правил нет до слова автора — это задача для решения»; «ПЕРВЫЕ
ПРИНЦИПЫ — справочник; толкование неопределённостей — вычисление»).

Что дано [в корпусе]:
  * ОП-геометрия дискретна изначально; тор всегда дискретный (01.09, пп. 1);
  * ноль осей не лежит на дискретном торе (01.09, п. 9);
  * ноль — состояние, не значение: S = {1..N}, 0 ∉ S, Next(N) = 1 —
    часы без нуля (00-ЛИНИЯ §2.5);
  * обратнополярная запись u = 1/r (STUDY-08); «−0 в hyper-арифметике:
    обратная величина — не представимое» (GLOSSARY, слово автора);
  * открытые вопросы 01.09: п. 8 — притягиваются ли нули осей; п. 10 —
    ось, проходящая дискретным отсчётом через +0, −0 и нейтральный 0
    другой оси.

Что следует [вычислено здесь]:
  1. Ноль оси — не узел, а ШОВ между N и 1. Next(N) = 1 шов не посещает,
     а ПЕРЕСЕКАЕТ; пересечение — событие с направлением: вперёд (N → 1)
     даёт 0{ОП+0}, назад (1 → N) — 0{ОП−0}; шов сам по себе, без
     направления, — 0{ОП0}. Три нуля — три состояния одного шва.
  2. Двойственность u = 1/r переводит подход к нулю в указатель
     масштаба: r → 0⁺ ⇔ u → +ω, r → 0⁻ ⇔ u → −ω. Потому обратный к
     0{ОП+0} есть +ω, к 0{ОП−0} — −ω (класс 22 со знаком), а обратный к
     0{ОП0} не представим — ровно то, что корпус называет «−0 в
     hyper-арифметике». Пересечение шва — смена знака указателя на
     двойственной оси (u проходит через ∞, как tg через π/2).
  3. Права (выведены через двойственность, сверены с IEEE 754 §6.3, где
     они совпадают, и названы, где расходятся): x + (±0) = x (сдвиг на
     ноль шагов); (±0)·x = ноль со знаком sign(x)·(±) — из
     (1/(±ω))·x = 1/(±ω/x); (+0) + (−0) = 0{ОП0} — направления гасятся
     (из 1/(+ω) + 1/(−ω): числитель −ω + ω — ноль без направления);
     IEEE здесь отвечает +0 по соглашению округления — расхождение
     названо, не скрыто; x/(±0) = ±ω·sign(x); x/0{ОП0} — отказ.
  4. Притяжение/отталкивание (п. 8): расстояние между сторонами шва —
     ровно один дискретный шаг при любом N; нули двух осей тора
     встречаются только в «углу» (шов × шов) — одном состоянии на тор;
     ни притяжения, ни отталкивания — расстояния суть инварианты
     дискретного строения, а не силы.
  5. Ось через нули другой оси (п. 10): на торе N₁ × N₂ состояний
     «шов оси 1 при узле k оси 2» ровно N₂, и столько же обратных; при
     каждом таком пересечении знак указателя двойственной оси 1
     меняется, координата оси 2 не меняется. Ортогональность не нужна
     (01.09, п. 6): пересечение — событие на узле, не угол.

Авторство: ОП-геометрия, часы без нуля, дискретный тор — AIAΣA (автор
проекта); IEEE 754-2019 §6.3 (знаковые нули — сверка); Мур 1966.
Оси с самопересечениями (01.09, пп. 3–7) замерены в соседнем модуле
`phicore/op_cross.py` (11.09): тор с самопересечениями как граф, β₁ тора 1 и
восьмёрки 2, +1 на каждое самопересечение, k пересечений — k общих узлов;
приёмка `proto/test_phicore_op_cross.py`. Здесь — тор БЕЗ самопересечений: это
область модуля, а не пробел замера.

NOT_MEASURED — назван целиком 12.09 (строка остаётся записью о том, что
было), и это отказ ПРЕДМЕТА, а не дома:

  ПРЕДМЕТ: гладкие (непрерывные) эквиваленты дискретной оси-тора —
  многообразие, касательное пространство, кривизна.
  ПРИЧИНА: величины такого рода в предмете нет. Слово корпуса прямое
  (01.09, п. 1): «тор всегда дискретный», и ось здесь есть множество узлов
  1..N с Next(N) = 1, а не образ отрезка. Это не «не умею посчитать»:
  считать нечего, ибо предмет иной.
  АДРЕС: письмо 01.09 «дискретность ОП-геометрии, тор, нули осей», п. 1;
  сам модуль (`phicore/op_axis.py`); оси с самопересечениями — соседний
  `phicore/op_cross.py`.
  ЧЕМ МЕРИЛОСЬ БЫ: носителя континуума дом не несёт вовсе — ближайшее, что
  есть, это класс 23 IGNORANCE, и он является дверью ИЗ IEEE наружу, а не
  континуумом; гладкая дорога потребовала бы другого дома целиком
  (дифференциальная геометрия), а не другой цены в разрядах.

  ЗАМЕР 12.09 `[вычислено здесь]`: утверждение «вне предмета» до этого дня
  держалось только словом — дверь пускала нецелое и молча усекала
  (`Axis(7.5)` давала N = 7, `_node` принимал нецелый узел). Теперь
  дискретность держится УСТРОЙСТВОМ: `_whole` отвергает float и bool
  построением на трёх входах (N, узел, число шагов), и запись проверяется
  приёмкой, а не читается на веру.
"""
from __future__ import annotations

from .errors import Unresolved
from .pointers import ScalePointer
from .rational import Q


def _sign(x) -> int:
    if isinstance(x, Seam):
        return x.sign
    if isinstance(x, ScalePointer):
        return 1 if x.coeff > 0 else (-1 if x.coeff < 0 else 0)
    v = x
    return (v > 0) - (v < 0)


def _whole(x, имя: str) -> None:
    """Вход оси — только целое: дискретность держится ДВЕРЬЮ, не оговоркой.

    Гладкий (нецелый) вход здесь не округляется и не усекается: ось есть
    дискретный тор узлов 1..N, и «полтора узла» — не величина этого предмета,
    а другой предмет. До 12.09 дверь молча усекала (`Axis(7.5)` давала N = 7);
    теперь это отказ построением, и запись «гладкие эквиваленты вне предмета»
    держится устройством, а не словом."""
    if isinstance(x, bool):
        raise TypeError(f"{имя} — bool числом на дискретной оси не считается")
    if isinstance(x, float):
        raise TypeError(f"{имя} — float запрещён построением: ось дискретна, "
                        "нецелый вход не усекается, а отвергается")
    if not isinstance(x, int):
        raise TypeError(f"{имя} — только целое (узлы 1..N дискретны)")


class Seam:
    """Ноль оси как состояние шва: sign +1 — пересечён вперёд (N → 1),
    −1 — назад (1 → N), 0 — шов без направления."""
    __slots__ = ("sign",)

    def __init__(self, sign: int = 0):
        if sign not in (-1, 0, 1):
            raise ValueError("знак шва: −1, 0 или +1")
        self.sign = sign

    @classmethod
    def plus(cls) -> "Seam":
        return cls(1)

    @classmethod
    def minus(cls) -> "Seam":
        return cls(-1)

    @classmethod
    def neutral(cls) -> "Seam":
        return cls(0)

    # --- двойственность u = 1/r
    def inverse(self):
        if self.sign == 0:
            raise Unresolved("inverse", self.record(), hint="обратная величина шва без направления не представима (корпус: «−0 в hyper-арифметике»)")
        return ScalePointer(Q(self.sign), 1)                 # ±ω — указатель масштаба, класс 22

    @classmethod
    def from_pointer_inverse(cls, p: ScalePointer) -> "Seam":
        """1/(±ω) → шов со знаком указателя; 1/ε — не ноль (ω), отказ."""
        if p.k <= 0:
            raise Unresolved("inverse", p.record(), hint="обратный к ε — ω, не ноль оси")
        return cls(_sign(p))

    # --- права
    def __add__(self, other):
        if isinstance(other, Seam):
            return Seam(self.sign if self.sign == other.sign else 0)   # противоположные направления гасятся
        return other                                                   # сдвиг на ноль шагов

    __radd__ = __add__

    def __neg__(self) -> "Seam":
        return Seam(-self.sign)

    def __sub__(self, other):
        return self + (-other) if isinstance(other, Seam) else (-other) + self

    def __rsub__(self, other):
        return other

    def __mul__(self, other):
        return Seam(self.sign * _sign(other))                       # (1/(±ω))·x = 1/(±ω/x)

    __rmul__ = __mul__

    def __rtruediv__(self, other):
        """x / (±0) = ±ω·sign(x); x / 0{ОП0} — отказ."""
        if self.sign == 0:
            raise Unresolved("/", repr(other), self.record(), hint="деление на шов без направления не представимо")
        return ScalePointer(Q(self.sign * _sign(other)), 1)

    def __truediv__(self, other):
        return Seam(self.sign * _sign(other))

    def __eq__(self, other) -> bool:
        return isinstance(other, Seam) and other.sign == self.sign

    def __hash__(self) -> int:
        return hash(("Seam", self.sign))

    # --- самопроверка шва (ось самоприменения, 09.10) ---
    def self_check(self, deep: bool = False) -> list:
        """Самопроверка шва встречными дорогами (ось «AIAΣA спрашивает AIAΣA»,
        слово автора 04.09; очередь строка 18; 09.10). Список находок, пусто —
        здоров; исключений не бросает.

        `[вычислено здесь]` (1) слот: знак из трёх; (2) двойственность u = 1/r —
        обратная к шву с направлением есть указатель ±ω того же знака, и обратно
        (`from_pointer_inverse`) возвращает тот же шов; шов без направления
        обратной не имеет — отказ поимённо; (3) права: x + шов = x (сдвиг на
        ноль шагов), шов + противоположный шов = шов без направления, шов + он
        же = он же, −шов меняет направление, шов·x несёт знак произведения,
        x/шов с направлением = ±ω·sign(x), x/шов без направления — отказ.
        При deep — права под указателями ±ω как пробами."""
        out = []
        s = self.sign
        if s not in (-1, 0, 1) or isinstance(s, bool):
            return [f"слот знака вне {{−1, 0, 1}}: {s!r}"]
        try:
            if s == 0:
                try:
                    self.inverse()
                    out.append("обратная к шву без направления не отказала")
                except Unresolved:
                    pass
            else:
                p = self.inverse()
                if not isinstance(p, ScalePointer) or _sign(p) != s or p.k != 1:
                    out.append(f"двойственность: обратная к {self.record()} есть {p!r}, не ±ω того же знака")
                elif Seam.from_pointer_inverse(p) != self:
                    out.append("двойственность: обратная к обратной не вернула шов")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"двойственность: {type(e).__name__}: {e}")
        try:
            x = Q(7, 3)
            if self + x != x or x + self != x:
                out.append("права: x + шов ≠ x")
            if (self + (-self)).sign != 0:
                out.append("права: шов + противоположный ≠ шов без направления")
            if (self + self) != self:
                out.append("права: шов + он же ≠ он же")
            if (-self).sign != -s:
                out.append("права: −шов не меняет направление")
            if (self * x).sign != s or (self * (-x)).sign != -s:
                out.append("права: шов·x не несёт знак произведения")
            if s == 0:
                try:
                    x / self
                    out.append("права: x/шов без направления не отказало")
                except Unresolved:
                    pass
            else:
                q = x / self
                if not isinstance(q, ScalePointer) or _sign(q) != s or _sign((-x) / self) != -s:
                    out.append("права: x/шов не есть ±ω·sign(x)")
            if deep:
                w = ScalePointer(Q(1), 1)
                if (self * w).sign != s or (self * (-w)).sign != -s:
                    out.append("права под пробой ω: шов·(±ω) не несёт знак")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"права: {type(e).__name__}: {e}")
        return out

    def record(self) -> str:
        return {1: "0{ОП+0}", -1: "0{ОП−0}", 0: "0{ОП0}"}[self.sign]

    __repr__ = record


class Axis:
    """Дискретная ось-тор: узлы 1..N, Next(N) = 1, ноль — шов между N и 1."""
    __slots__ = ("N",)

    def __init__(self, N: int):
        _whole(N, "N")
        if N < 1:
            raise ValueError("узлов на оси не меньше одного")
        self.N = int(N)

    def next(self, k: int) -> int:
        self._node(k)
        return 1 if k == self.N else k + 1

    def prev(self, k: int) -> int:
        self._node(k)
        return self.N if k == 1 else k - 1

    def _node(self, k: int) -> None:
        _whole(k, "узел")
        if not 1 <= k <= self.N:
            raise ValueError(f"узел {k} вне оси 1..{self.N}: ноль и отрицательные — не узлы (0 ∉ S)")

    def step(self, k: int, d: int) -> tuple:
        """Сдвиг на d шагов: (узел, [пересечения шва со знаком])."""
        self._node(k)
        _whole(d, "число шагов")
        crossings = []
        for _ in range(abs(d)):
            if d > 0:
                if k == self.N:
                    crossings.append(Seam.plus())
                k = self.next(k)
            else:
                if k == 1:
                    crossings.append(Seam.minus())
                k = self.prev(k)
        return k, crossings

    def visits_zero(self, k: int, d: int) -> int:
        """Сколько раз узел, а не событие, оказался нулём: всегда 0 — часы без нуля."""
        node, _ = self.step(k, d)
        return 0 if 1 <= node <= self.N else 1

    def seam_distance(self) -> int:
        """Расстояние между сторонами шва (узлы N и 1) через шов: один шаг при любом N."""
        return 1

    # --- самопроверка оси (ось самоприменения, 09.10) ---
    def self_check(self, deep: bool = False) -> list:
        """Самопроверка дискретной оси встречными дорогами. Список находок,
        пусто — здорова; исключений не бросает.

        `[вычислено здесь]` (1) слот: N — целое ≥ 1; (2) часы без нуля: next и
        prev взаимно обратны на всех узлах, обход N шагов вперёд возвращает
        узел и пересекает шов ровно один раз вперёд, назад — один раз назад;
        d шагов и затем −d возвращают узел; ноль узлом не бывает
        (`visits_zero` = 0); расстояние через шов — один шаг; дверь отвергает
        нецелое. При deep — обход от каждого узла, иначе от первого и
        последнего."""
        out = []
        N = self.N
        if not isinstance(N, int) or isinstance(N, bool) or N < 1:
            return [f"слот N вне канона: {N!r}"]
        try:
            for k in range(1, N + 1):
                if self.prev(self.next(k)) != k or self.next(self.prev(k)) != k:
                    out.append(f"next и prev не взаимно обратны на узле {k}")
                    break
            starts = range(1, N + 1) if deep else sorted({1, N})
            for k in starts:
                node, cr = self.step(k, N)
                if node != k or len(cr) != 1 or cr[0].sign != 1:
                    out.append(f"обход N вперёд от {k}: узел {node}, пересечений {[c.record() for c in cr]}")
                node, cr = self.step(k, -N)
                if node != k or len(cr) != 1 or cr[0].sign != -1:
                    out.append(f"обход N назад от {k}: узел {node}, пересечений {[c.record() for c in cr]}")
                node, _ = self.step(k, 3)
                back, _ = self.step(node, -3)
                if back != k:
                    out.append(f"три шага вперёд и назад от {k} не вернули узел")
                if self.visits_zero(k, 2 * N + 1) != 0:
                    out.append("ноль оказался узлом — часы без нуля нарушены")
            if self.seam_distance() != 1:
                out.append("расстояние через шов не один шаг")
            for bad in (1.5, True):
                try:
                    self._node(bad)
                    out.append(f"дверь пропустила нецелый узел {bad!r}")
                except TypeError:
                    pass
        except Exception as e:                                  # noqa: BLE001
            out.append(f"обход: {type(e).__name__}: {e}")
        return out


class Torus:
    """Произведение осей; нули — состояния «шов оси i × узлы прочих»."""
    __slots__ = ("axes",)

    def __init__(self, *axes: Axis):
        self.axes = tuple(axes)

    def nodes(self) -> int:
        n = 1
        for a in self.axes:
            n *= a.N
        return n

    def zero_states(self, i: int) -> int:
        """Состояний «шов оси i при узлах остальных осей» — по одному на узел прочих."""
        n = 1
        for j, a in enumerate(self.axes):
            if j != i:
                n *= a.N
        return n

    def corners(self) -> int:
        """Состояний «все оси на шве» — одно на тор: только там нули осей встречаются."""
        return 1

    def cross(self, i: int, nodes: tuple, d: int) -> tuple:
        """Ось i делает d шагов при прочих осях на месте: (узлы, пересечения шва оси i,
        знак указателя двойственной оси i после пересечений)."""
        k, crossings = self.axes[i].step(nodes[i], d)
        out = list(nodes); out[i] = k
        flips = len(crossings)
        return tuple(out), crossings, (-1) ** flips

    # --- самопроверка тора (ось самоприменения, 09.10) ---
    def self_check(self, deep: bool = False) -> list:
        """Самопроверка тора: узлов — произведение осей; нулевых состояний оси
        i — произведение прочих осей; угол один; при обходе оси i на N_i шагов
        знак указателя двойственной оси меняется ровно один раз, прочие
        координаты на месте; оси проверяют себя."""
        out = []
        if not self.axes or not all(isinstance(a, Axis) for a in self.axes):
            return ["слот осей: пусто либо чужой носитель"]
        for i, a in enumerate(self.axes):
            out += [f"ось {i}: {f}" for f in a.self_check(deep)]
        try:
            total = 1
            for a in self.axes:
                total *= a.N
            if self.nodes() != total:
                out.append(f"узлов {self.nodes()} ≠ произведение осей {total}")
            for i, a in enumerate(self.axes):
                if self.zero_states(i) * a.N != total:
                    out.append(f"нулевых состояний оси {i}: {self.zero_states(i)} ≠ {total}/{a.N}")
            if self.corners() != 1:
                out.append("углов на торе не один")
            start = tuple(1 for _ in self.axes)
            for i, a in enumerate(self.axes):
                nodes, cr, flip = self.cross(i, start, a.N)
                if nodes != start or len(cr) != 1 or flip != -1:
                    out.append(f"обход оси {i} на N шагов: узлы {nodes}, пересечений {len(cr)}, знак {flip}")
                nodes2, cr2, flip2 = self.cross(i, start, 2 * a.N)
                if nodes2 != start or len(cr2) != 2 or flip2 != 1:
                    out.append(f"двойной обход оси {i}: узлы {nodes2}, пересечений {len(cr2)}, знак {flip2}")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"тор: {type(e).__name__}: {e}")
        return out


# ---------- самопроверка модуля (ось самоприменения, 09.10) ----------

def self_check(deep: bool = False) -> list:
    """Самопроверка осей ОП-геометрии встречными дорогами; список находок,
    пусто — здорова; исключений не бросает.

    `[вычислено здесь]` Три шва проверяют себя (двойственность и права);
    права, выведенные из двойственности, против IEEE 754 §6.3 [внешний
    источник]: совпадают на x + (±0) = x и x/(±0) = ±ω·sign(x), расходятся
    на (+0) + (−0) — здесь шов без направления, у IEEE +0 по соглашению
    округления, и расхождение названо; оси N ∈ {1, 2, 7} (deep — ещё 12) и
    тор 3 × 4 проверяют себя; дверь оси отвергает float и bool."""
    out = []
    for s in (Seam.plus(), Seam.minus(), Seam.neutral()):
        out += [f"{s.record()}: {f}" for f in s.self_check(deep)]
    try:
        x = Q(5, 2)
        if (Seam.plus() + Seam.minus()).sign != 0:
            out.append("(+0) + (−0) ≠ 0{ОП0} — расхождение с IEEE потеряно")
        if Seam.plus() + x != x or x + Seam.minus() != x:
            out.append("x + (±0) ≠ x — согласие с IEEE потеряно")
        if _sign(x / Seam.plus()) != 1 or _sign(x / Seam.minus()) != -1:
            out.append("x/(±0) ≠ ±ω — согласие с IEEE потеряно")
        try:
            Axis(7.5)
            out.append("дверь оси пропустила float")
        except TypeError:
            pass
        try:
            Axis(True)
            out.append("дверь оси пропустила bool")
        except TypeError:
            pass
    except Exception as e:                                      # noqa: BLE001
        out.append(f"права против IEEE: {type(e).__name__}: {e}")
    sizes = [1, 2, 7] + ([12] if deep else [])
    for n in sizes:
        try:
            out += [f"ось N={n}: {f}" for f in Axis(n).self_check(deep)]
        except Exception as e:                                  # noqa: BLE001
            out.append(f"ось N={n}: {type(e).__name__}: {e}")
    try:
        out += [f"тор 3×4: {f}" for f in Torus(Axis(3), Axis(4)).self_check(deep)]
    except Exception as e:                                      # noqa: BLE001
        out.append(f"тор 3×4: {type(e).__name__}: {e}")
    return out
