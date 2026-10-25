"""Оси ОП-геометрии с самопересечениями и взаимными пересечениями:
дискретный тор как граф, пересечение — отождествление узлов. Закрывает
NOT_MEASURED из `op_axis.py` (01.09, пп. 3–7): здесь тор С самопересечениями.

Что дано [в корпусе] (письмо 01.09 «дискретность ОП-геометрии, тор, нули
осей», пп. 1–7, 9):
  п. 1  тор всегда дискретный;
  п. 2  незамкнутый тор — прямая/кривая;
  п. 3  тор с собственным пересечением (типа восьмёрка) топологически —
        другое многообразие;
  п. 4  ограничения на единственность пересечения нет: множественные
        пересечения равноправно допустимы;
  п. 5  в любой дискретной точке может быть пересечение с любой другой
        осью/осями;
  п. 6  оси не обязаны быть ортогональными: ось может пересекаться с
        другими осями и собой, каким бы причудливым многообразием это
        ни было;
  п. 7  следствие о единственности трёх перпендикуляров НЕ выполняется
        автоматически;
  п. 9  ноль осей не находится на дискретном торе (`op_axis`: ноль — шов
        между N и 1, событие, не узел).

Выбор модели [вычислено здесь] — то, чего письмо не говорит, названо
выбором, не словом автора:
  М1. Дискретная ось — граф: узлы 1..N, рёбра k → Next(k), Next(N) = 1
      (п. 1; часы без нуля из `op_axis`). Незамкнутая ось (п. 2) — тот же
      граф без ребра N → 1: путь.
  М2. Пересечение — ОТОЖДЕСТВЛЕНИЕ двух узлов (одной оси — самопересечение,
      разных осей — взаимное). Нового узла, угла, метрики не возникает:
      «ортогональность» в модели не параметр, требовать её не из чего
      (п. 6); отождествлять можно любые узлы (п. 5), любое число раз
      (п. 4). Восьмёрка — одно самопересечение.
  М3. Каждая ось сохраняет свой Next: шаг вдоль оси однозначен по метке
      оси, отождествление говорит лишь, что положения совпали. Потому в
      узле пересечения рёбер четыре, а ход — по-прежнему по одной оси
      (`walk` считает через `Axis.step`, швы — `Seam`).
  М4. Петли и кратные рёбра допускаются (N = 1 — петля, N = 2 — двуугольник,
      отождествление N ~ 1 — петля шва).
  М5. Мера «другого многообразия» (п. 3) — точные инварианты фактор-графа:
      узлы V, рёбра E, χ = V − E, β₀ (компоненты), β₁ = E − V + β₀ —
      целые; в приёмке β сверяются с рангом матрицы инцидентности над ℚ
      (`proto/discrete_manifold.CellComplex`). Слово «многообразие» в
      письме читается как «дискретное пространство»: в узле пересечения
      степень 4, а не 2 — замер, не спор со словом.
  М6. «Единственность трёх перпендикуляров» (п. 7) читается как
      единственность ОБЩЕЙ ТОЧКИ осей: две оси делят ровно один узел;
      три оси проходят через один общий узел. В модели число общих узлов
      k — свойство набора отождествлений, не осей: k = 1 — один случай
      из многих. Иное чтение п. 7 — не здесь.
  М7. Отождествления, задевающие концы шва N и 1, допускаются (п. 5:
      любая точка); шов остаётся ребром-событием, узлом не становится
      (п. 9); он делается петлёй ровно тогда, когда N и 1 попали в один
      класс.

Что измерено [вычислено здесь] (перебором малых окон в
`proto/test_phicore_op_cross.py`):
  (а) тор: V = E = N, χ = 0, β = (1, 1); восьмёрка (любая пара узлов):
      V = N − 1, E = N, χ = −1, β = (1, 2) — β₁ различает, «другое
      многообразие» (п. 3); незамкнутая ось: β = (1, 0), χ = 1 (п. 2);
  (б) каждое новое самопересечение — β₁ + 1 (п. 4): β₁ = 1 + число
      результативных отождествлений; все N узлов в один — букет из N петель,
      β₁ = N; повторное отождествление одного класса ничего не меняет;
  (в) две оси с k точками пересечения делят ровно k узлов, β₁ = k + 1;
      единственность общей точки — только при k = 1 (п. 7); три оси с
      попарными пересечениями: через один общий узел (β₁ = 3) — N³
      наборов из N⁶, прочие — β₁ = 4;
  (г) шов оси сохраняется при отождествлениях вне {N, 1}: классы N и 1
      различны, ход N → 1 даёт 0{ОП+0}, 1 → N — 0{ОП−0}, как без
      отождествлений; петлёй шов становится ⇔ N ~ 1 (М7).
Не измерено: п. 8 и п. 10 письма — вне этого модуля (`op_axis`, пп. 4–5
докстринга); гладкие эквиваленты не рассматриваются (корпус).

Авторство: ОП-геометрия, дискретный тор, самопересечения осей — AIAΣA
(автор проекта); χ, числа Бетти, ранг над ℚ — Эйлер 1758, Бетти 1871,
Пуанкаре 1895 [внешний источник].
"""
from __future__ import annotations

from .op_axis import Axis, Seam


class CrossedAxes:
    """Набор дискретных осей-торов с отождествлениями узлов (М1–М7).

    Узел — пара (i, k): ось i, отсчёт k ∈ 1..N_i. Класс узла — представитель
    в объединении (наименьшая пара). Ось из `open` — незамкнутая (без ребра
    N → 1)."""
    __slots__ = ("axes", "open", "pairs", "_parent")

    def __init__(self, *axes, open=(), pairs=()):
        self.axes = tuple(a if isinstance(a, Axis) else Axis(a) for a in axes)
        if not self.axes:
            raise ValueError("нужна хотя бы одна ось")
        self.open = frozenset(int(i) for i in open)
        for i in self.open:
            self._axis(i)
        self.pairs: tuple = ()
        self._parent = {(i, k): (i, k) for i, a in enumerate(self.axes) for k in range(1, a.N + 1)}
        for u, v in pairs:
            self.identify(u, v)

    # --- узлы и классы -----------------------------------------------------
    def _axis(self, i: int) -> Axis:
        if not 0 <= i < len(self.axes):
            raise ValueError(f"оси {i} нет: оси 0..{len(self.axes) - 1}")
        return self.axes[i]

    def node(self, i: int, k: int) -> tuple:
        self._axis(i)._node(k)
        return (int(i), int(k))

    def cls(self, u) -> tuple:
        """Класс узла — корень объединения (сжатие пути; арифметика только целая)."""
        u = self.node(*u)
        root = u
        while self._parent[root] != root:
            root = self._parent[root]
        while self._parent[u] != root:
            self._parent[u], u = root, self._parent[u]
        return root

    def identify(self, u, v) -> bool:
        """Отождествить узлы u ~ v (М2). True — пересечение новое (классы
        были различны), False — узлы уже в одном классе."""
        ru, rv = self.cls(u), self.cls(v)
        self.pairs += ((self.node(*u), self.node(*v)),)
        if ru == rv:
            return False
        if rv < ru:
            ru, rv = rv, ru
        self._parent[rv] = ru
        return True

    def copy(self) -> "CrossedAxes":
        c = CrossedAxes(*self.axes, open=self.open)
        c.pairs = self.pairs
        c._parent = dict(self._parent)
        return c

    def nodes(self):
        """Все исходные узлы (до отождествления)."""
        return list(self._parent)

    def vertices(self) -> list:
        """Узлы фактор-графа — классы."""
        return sorted({self.cls(u) for u in self._parent})

    def edges(self) -> list:
        """Рёбра фактор-графа: (ось i, отсчёт k, класс k, класс Next(k));
        у незамкнутой оси ребра N → 1 нет (М1). Петли и кратные рёбра
        сохраняются (М4)."""
        out = []
        for i, a in enumerate(self.axes):
            for k in range(1, a.N + 1):
                if k == a.N and i in self.open:
                    continue
                out.append((i, k, self.cls((i, k)), self.cls((i, a.next(k)))))
        return out

    # --- инварианты ---------------------------------------------------------
    def counts(self) -> tuple:
        """(V, E) фактор-графа."""
        return len(self.vertices()), len(self.edges())

    def euler(self) -> int:
        """χ = V − E."""
        V, E = self.counts()
        return V - E

    def components(self) -> int:
        """β₀ — компоненты связности фактор-графа по его рёбрам."""
        parent = {v: v for v in self.vertices()}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for _, _, a, b in self.edges():
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[rb] = ra
        return len({find(v) for v in parent})

    def betti(self) -> tuple:
        """(β₀, β₁) графа: β₁ = E − V + β₀ — точное целое (М5)."""
        V, E = self.counts()
        b0 = self.components()
        return b0, E - V + b0

    def degree(self, u) -> int:
        """Число концов рёбер в классе узла u (петля считается дважды)."""
        c = self.cls(u)
        return sum((a == c) + (b == c) for _, _, a, b in self.edges())

    def is_plain_torus(self) -> bool:
        """Тор без самопересечений: связен, все степени 2, β = (1, 1)."""
        return self.components() == 1 and all(self.degree(v) == 2 for v in self.vertices()) and self.betti() == (1, 1)

    def merges(self) -> int:
        """Число результативных отождествлений = Σ N_i − V."""
        return len(self._parent) - len(self.vertices())

    # --- самопроверка (ось «AIAΣA спрашивает AIAΣA», 04.09) ------------------
    def self_check(self, deep: bool = False) -> list:
        """Самопроверка скрещённых осей ВСТРЕЧНОЙ дорогой (ось самоприменения).

        `[вычислено здесь]` Дороги — тождества графа, каждое считается не
        тем счётом, что даёт проверяемое число: (1) Эйлер двумя дорогами —
        χ = V − E прямым счётом против χ = β₀ − β₁ через числа Бетти, где
        β₀ берётся обходом связности, а β₁ = E − V + β₀; (2) лемма о
        рукопожатиях: сумма степеней классов обязана равняться 2E, причём
        степень считается перебором концов рёбер, а E — длиной списка
        рёбер; (3) отождествления: `merges` = Σ Nᵢ − V — число слияний
        против числа классов; (4) всякое ребро ведёт из класса в класс, и
        оба конца обязаны быть среди вершин.

        Граница объявлена: у НЕЗАМКНУТОЙ оси ребра N → 1 нет (М1), поэтому
        числа рёбер у открытой и замкнутой осей разные — это условие
        построения, а не находка; петли и кратные рёбра сохраняются (М4) и
        считаются как есть (петля даёт две единицы степени). Возвращает
        список находок; исключений не бросает.
        """
        out = []
        try:
            V, E = self.counts()
            verts = self.vertices()
            edges = self.edges()
        except Exception as e:                                  # noqa: BLE001
            out.append(f"фактор-граф не собран: {type(e).__name__}: {e}")
            return out
        if (V, E) != (len(verts), len(edges)):
            out.append(f"счёт (V, E) = {(V, E)} не отвечает спискам "
                       f"{(len(verts), len(edges))}")
        b0, b1 = self.betti()
        if self.euler() != b0 - b1:
            out.append(f"Эйлер двумя дорогами расходится: V − E = "
                       f"{self.euler()}, β₀ − β₁ = {b0 - b1}")
        if b0 != self.components():
            out.append(f"β₀ = {b0} не равно числу компонент "
                       f"{self.components()}")
        deg_sum = sum(self.degree(v) for v in verts)
        if deg_sum != 2 * E:
            out.append(f"лемма о рукопожатиях нарушена: сумма степеней "
                       f"{deg_sum} ≠ 2E = {2 * E}")
        if self.merges() != len(self._parent) - V:
            out.append(f"число отождествлений {self.merges()} не равно "
                       f"Σ Nᵢ − V = {len(self._parent) - V}")
        vs = set(verts)
        for i, k, a_cls, b_cls in edges:
            if a_cls not in vs or b_cls not in vs:
                out.append(f"ребро оси {i} отсчёта {k} ведёт мимо вершин: "
                           f"{a_cls} → {b_cls}")
                break
        if not deep:
            return out
        twin = self.copy()
        if (twin.counts(), twin.betti()) != ((V, E), (b0, b1)):
            out.append("копия дала другие инварианты, чем оригинал")
        for i, axis in enumerate(self.axes):
            own = {self.cls((i, k)) for k in range(1, axis.N + 1)}
            if not own <= vs:
                out.append(f"классы оси {i} не все среди вершин")
        return out

    # --- пересечения осей ---------------------------------------------------
    def shared(self, i: int, j: int) -> int:
        """Общих узлов осей i и j — классов, содержащих узлы обеих осей."""
        self._axis(i); self._axis(j)
        ci = {self.cls((i, k)) for k in range(1, self.axes[i].N + 1)}
        cj = {self.cls((j, k)) for k in range(1, self.axes[j].N + 1)}
        return len(ci & cj)

    def common_to_all(self) -> int:
        """Классов, содержащих узлы КАЖДОЙ оси (три оси через одну точку — 1)."""
        sets = [{self.cls((i, k)) for k in range(1, a.N + 1)} for i, a in enumerate(self.axes)]
        common = sets[0]
        for s in sets[1:]:
            common = common & s
        return len(common)

    # --- шов и ход ----------------------------------------------------------
    def seam(self, i: int) -> tuple:
        """Шов оси i: (класс N, класс 1, петля ли). Шов — ребро-событие, не
        узел (п. 9); петля ⇔ N ~ 1 (М7). У незамкнутой оси шва нет."""
        a = self._axis(i)
        if i in self.open:
            raise ValueError(f"ось {i} незамкнута: шва нет (п. 2 — прямая/кривая)")
        cN, c1 = self.cls((i, a.N)), self.cls((i, 1))
        return cN, c1, cN == c1

    def walk(self, i: int, k: int, d: int) -> tuple:
        """d шагов вдоль оси i с отсчёта k (М3): (класс конечного узла,
        [пересечения шва со знаком] через `Axis.step`). Незамкнутая ось
        через край не ходит — ValueError."""
        a = self._axis(i)
        node, crossings = a.step(k, d)
        if i in self.open and crossings:
            raise ValueError(f"ось {i} незамкнута: ход через край N → 1 невозможен")
        return self.cls((i, node)), crossings

    # --- клетки для сверки аппаратом discrete_manifold ----------------------
    def cells(self) -> list:
        """Клетки фактор-графа для `CellComplex(cells, cell_dim, cell_boundary)`:
        ("v", класс) и ("e", ось, отсчёт, класс, класс)."""
        return [("v", v) for v in self.vertices()] + [("e",) + e for e in self.edges()]

    def record(self) -> str:
        V, E = self.counts()
        b0, b1 = self.betti()
        axes = "×".join(f"{'(' if i in self.open else ''}{a.N}{')' if i in self.open else ''}" for i, a in enumerate(self.axes))
        return f"оси {axes}, отождествлений {len(self.pairs)} (результативных {self.merges()}): V={V} E={E} χ={V - E} β=({b0},{b1})"

    __repr__ = record


def cell_dim(c) -> int:
    return 0 if c[0] == "v" else 1


def cell_boundary(c):
    """∂("e", i, k, a, b) = ("v", b) − ("v", a); у вершины граница пуста.
    Петля даёт нулевую границу, кратные рёбра — отдельные клетки."""
    if c[0] == "v":
        return []
    return [(1, ("v", c[4])), (-1, ("v", c[3]))]


def all_pairs(X: CrossedAxes, within: int | None = None, between: tuple | None = None) -> list:
    """Все кандидаты отождествления: `within=i` — пары узлов оси i
    (самопересечения); `between=(i, j)` — пары узел оси i × узел оси j."""
    if within is not None:
        N = X.axes[within].N
        return [((within, a), (within, b)) for a in range(1, N + 1) for b in range(a + 1, N + 1)]
    if between is not None:
        i, j = between
        return [((i, a), (j, b)) for a in range(1, X.axes[i].N + 1) for b in range(1, X.axes[j].N + 1)]
    raise ValueError("укажите within=i или between=(i, j)")
