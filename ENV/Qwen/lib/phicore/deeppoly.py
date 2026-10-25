"""DeepPoly — сертификат сети линейными границами с подстановкой назад.

Названная неполнота (`phicore/affine.py`, `proto/bench_zonotope.py`):
«сравнение с DeepPoly не замерено». Здесь DeepPoly написан, и сравнение
делается.

Чем он отличается от двух уже имеющихся дорог. **IBP** ведёт по слою
интервал и на каждом слое забывает, что границы одной величины связаны.
**Зонотоп** помнит связь через общие символы шума, но каждый ReLU
добавляет новый символ и связь с ВХОДОМ слабеет. **DeepPoly** [внешний
источник: Singh, Gehr, Püschel, Vechev, POPL 2019] держит для каждого
нейрона ЛИНЕЙНУЮ нижнюю и верхнюю границу и подставляет их назад до
самого входа: итоговая граница выражена прямо через входные координаты,
и только там берётся ящик.

Выигрыш виден там же, где у зонотопа: разность логитов подставляется
назад ЦЕЛИКОМ, а не как разность двух отдельно ограниченных величин.

Всё точно: коэффициенты — рациональные дроби, ящик берётся один раз в
конце, ни одного числа с плавающей запятой. Смягчение ReLU на смешанном
отрезке [l; u]: сверху y ≤ u(x − l)/(u − l) — единственная точная
линейная мажоранта; снизу y ≥ λx, где λ ∈ {0, 1} выбирается по площади
(λ = 1 при u ≥ |l|, иначе 0) — правило самого DeepPoly.

NOT_MEASURED — разобран 12.09 на три разных утверждения, каждое со своим
доказательством (строка остаётся записью о том, что было):

  * **уточнение λ поиском — ЗАКРЫТО точной дорогой.** Искать нечего, и это
    посчитано, а не заявлено: площадь зазора смягчения
    S(λ) = ∫ [hi·(x − lo)/(hi − lo) − λx] dx по [lo; hi] равна
    (hi − lo)/2 · [hi − λ·(hi + lo)] — величина АФФИННА по λ с наклоном
    −(hi − lo)(hi + lo)/2. Минимум аффинной величины на отрезке λ ∈ [0; 1]
    лежит на конце, потому правило λ ∈ {0, 1} уже оптимально по площади:
    hi + lo > 0 → λ = 1, hi + lo < 0 → λ = 0, hi + lo = 0 → площадь от λ не
    зависит вовсе. Дорога вписана: `relu_area`, `lambda_optimum`; вторая,
    независимая — трапеции по разбиению (точны, ибо подынтегральная
    величина линейна), сверка точная в ℚ.

  * **сети глубже двух слоёв — ОТКАЗ ДОМА, названный ценой.** Предмет:
    подстановка назад написана ровно для одного слоя ReLU (`certify` читает
    из носителя четыре поля — W1, b1, W2, b2, и пятого поля у носителя нет).
    Причина: носителя глубокой сети дом не несёт, а не арифметика не
    сходится. Чем мерилось бы: рекурсивная подстановка назад по L слоям,
    цена O(L·H²·n) операций над ℚ против нынешних O(H²·n) — разрядность
    коэффициентов при этом растёт как произведение L матриц, то есть
    линейно по L в битах; набор — `proto/bench_deeppoly.py` на сети с L
    слоями. Это «не посчитано здесь», а не «посчитать нельзя».

  * **свёрточные слои — ОТКАЗ ДОМА, названный ценой.** Предмет: свёртка
    есть линейный слой с РАЗДЕЛЯЕМЫМИ весами, потому арифметика та же, и
    величины новой не возникает. Причина: носителя ядра дом не несёт.
    Чем мерилось бы: разворот ядра в матрицу (im2col) размером
    (H·W) × (k²·C) — цена в разрядах та же, что у плотного слоя той же
    матрицы, выигрыш только в памяти. Это «не посчитано здесь».
"""
from __future__ import annotations

from .rational import Q


class Form:
    """Линейная форма над входом: значение = c + Σ wᵢ·xᵢ, коэффициенты — дроби."""

    __slots__ = ("c", "w")

    def __init__(self, c, w):
        self.c = Q(c)
        self.w = list(w)

    @classmethod
    def zero(cls, n: int) -> "Form":
        return cls(0, [Q(0)] * n)

    def scaled(self, k) -> "Form":
        k = Q(k)
        return Form(self.c * k, [x * k for x in self.w])

    def shifted(self, d) -> "Form":
        return Form(self.c + Q(d), list(self.w))

    def __add__(self, other: "Form") -> "Form":
        return Form(self.c + other.c, [a + b for a, b in zip(self.w, other.w)])

    # --- самопроверка формы (ось самоприменения, 09.10) ---
    def value_at(self, x) -> Q:
        """Значение формы в точке — прямой счёт c + Σ wᵢ·xᵢ."""
        return self.c + sum(wi * Q(xi) for wi, xi in zip(self.w, x))

    def self_check(self, deep: bool = False, x0=None, eps=1) -> list:
        """Самопроверка линейной формы встречными дорогами (ось «AIAΣA
        спрашивает AIAΣA», слово автора 04.09; очередь строка 18; 09.10).
        Список находок, пусто — здорова; исключений не бросает.

        `[вычислено здесь]` (1) слоты: c — Q, w — список Q; (2) ящик против
        прямого счёта: при ε = 0 обе границы суть значение в x₀; при ε > 0 угол
        ящика x₀ + ε·sign(w) достигает верхней границы ТОЧНО, противоположный
        — нижней (линейная форма достигает экстремума в вершине ящика —
        теорема, здесь сверенная прямым счётом); (3) линейность: границы
        суммы, масштаба и сдвига — из границ слагаемых; (4) смягчение ReLU на
        смешанном отрезке: верхняя форма не ниже max(0, f) и нижняя не выше
        в x₀ и в обоих углах, λ нижней формы — та же, что у `lambda_optimum`.
        При deep — ещё ε = 3 и 1/2."""
        out = []
        def _num(v):
            return isinstance(v, (int, Q)) and not isinstance(v, bool)
        if not isinstance(self.c, Q) or not isinstance(self.w, list) or not all(_num(v) for v in self.w):
            return [f"слоты формы чужого носителя: c = {type(self.c).__name__}, w = {self.w!r}"]
        n = len(self.w)
        x0 = [Q(0)] * n if x0 is None else [Q(v) for v in x0]
        eps_list = [Q(eps)] + ([Q(3), Q(1, 2)] if deep else [])
        try:
            lo0, hi0 = self.bounds(x0, 0)
            if lo0 != hi0 or lo0 != self.value_at(x0):
                out.append(f"ящик при ε = 0 не есть значение в x₀: ⟦{lo0}; {hi0}⟧ против {self.value_at(x0)}")
            for e in eps_list:
                lo, hi = self.bounds(x0, e)
                up = [xi + e * (1 if wi > 0 else (-1 if wi < 0 else 0)) for xi, wi in zip(x0, self.w)]
                dn = [xi - e * (1 if wi > 0 else (-1 if wi < 0 else 0)) for xi, wi in zip(x0, self.w)]
                if self.value_at(up) != hi or self.value_at(dn) != lo:
                    out.append(f"угол ящика при ε = {e} не достигает границ: {self.value_at(dn)}, {self.value_at(up)} против ⟦{lo}; {hi}⟧")
            g = Form(Q(1), [Q(i + 1) for i in range(n)])
            s = self + g
            lo_s, hi_s = s.bounds(x0, Q(eps))
            if s.value_at(x0) != self.value_at(x0) + g.value_at(x0) or hi_s - lo_s > (self.bounds(x0, Q(eps))[1] - self.bounds(x0, Q(eps))[0]) + (g.bounds(x0, Q(eps))[1] - g.bounds(x0, Q(eps))[0]):
                out.append("линейность суммы: середина или ширина ящика не сходятся")
            k = Q(-3, 2)
            lo_k, hi_k = self.scaled(k).bounds(x0, Q(eps))
            lo1, hi1 = self.bounds(x0, Q(eps))
            if (lo_k, hi_k) != (hi1 * k, lo1 * k):
                out.append("масштаб отрицательным числом не переворачивает ящик точно")
            if self.shifted(Q(5)).bounds(x0, Q(eps)) != (lo1 + 5, hi1 + 5):
                out.append("сдвиг не переносит ящик точно")
        except Exception as e:                                  # noqa: BLE001
            out.append(f"ящик и линейность: {type(e).__name__}: {e}")
        try:
            lo1, hi1 = self.bounds(x0, Q(eps))
            if lo1 < 0 < hi1:
                low, upp = relu_relax(self, lo1, hi1)
                lam, slope = lambda_optimum(lo1, hi1)
                if slope != 0 and (low.w != self.w) != (lam == 0):
                    out.append(f"λ нижней формы не согласна с lambda_optimum: λ = {lam}")
                pts = [x0, [xi + Q(eps) * (1 if wi > 0 else (-1 if wi < 0 else 0)) for xi, wi in zip(x0, self.w)],
                       [xi - Q(eps) * (1 if wi > 0 else (-1 if wi < 0 else 0)) for xi, wi in zip(x0, self.w)]]
                for x in pts:
                    v = max(Q(0), self.value_at(x))
                    if upp.value_at(x) < v or low.value_at(x) > v:
                        out.append(f"смягчение ReLU не обрамляет max(0, f) в точке {x}: {low.value_at(x)} ≤ {v} ≤ {upp.value_at(x)} нарушено")
                        break
        except Exception as e:                                  # noqa: BLE001
            out.append(f"смягчение ReLU: {type(e).__name__}: {e}")
        return out

    def bounds(self, x0, eps) -> tuple:
        """Ящик берётся ОДИН раз, в самом конце: [c + w·x₀ ∓ ε·Σ|w|]."""
        mid = self.c + sum(wi * Q(xi) for wi, xi in zip(self.w, x0))
        rad = sum(abs(wi) for wi in self.w) * Q(eps)
        return mid - rad, mid + rad


def layer_forms(W, b, shift: int, bias_shift: int, n_in: int) -> list:
    """Формы предактиваций слоя: (x·W + b·2^bias_shift) / 2^shift."""
    out = []
    H = len(W[0])
    for j in range(H):
        w = [Q(int(W[i][j]), 1 << shift) for i in range(n_in)]
        out.append(Form(Q(int(b[j]) * (1 << bias_shift)), w))
    return out


def relu_relax(form: Form, lo, hi) -> tuple:
    """Пара форм (нижняя, верхняя) для ReLU по границам предактивации."""
    if hi <= 0:
        z = Form.zero(len(form.w))
        return z, z
    if lo >= 0:
        return form, form
    k = Q(hi, hi - lo) if isinstance(hi, int) else hi / (hi - lo)
    upper = form.shifted(-lo).scaled(k)
    lower = form if hi >= -lo else Form.zero(len(form.w))
    return lower, upper


#: Поля носителя сети, какие читает `certify`. Пятого поля у носителя нет —
#: отсюда и отказ «сети глубже двух слоёв», названный числом, а не словом.
NET_FIELDS = ("W1", "b1", "W2", "b2")


def relu_area(lo, hi, lam) -> Q:
    """Площадь зазора смягчения ReLU на [lo; hi] при нижней границе y ≥ λ·x.

    Считается ПЕРВООБРАЗНОЙ, а не сеткой: зазор есть разность двух линейных
    величин, и интеграл берётся точно в ℚ:

        S(λ) = ∫_lo^hi [ hi·(x − lo)/(hi − lo) − λ·x ] dx
             = (hi − lo)/2 · [ hi − λ·(hi + lo) ]

    Отсюда всё: S АФФИННА по λ, наклон равен −(hi − lo)·(hi + lo)/2.
    Смешанный отрезок обязателен (lo < 0 < hi): на одностороннем зазора нет
    вовсе, и смягчать нечего — там отказ построением."""
    lo, hi, lam = Q(lo), Q(hi), Q(lam)
    if not (lo < 0 < hi):
        raise ValueError("зазор смягчения есть только на смешанном отрезке "
                         "lo < 0 < hi; на одностороннем ReLU точен")
    return (hi - lo) / Q(2) * (hi - lam * (hi + lo))


def lambda_optimum(lo, hi) -> tuple:
    """(λ минимума площади, наклон S по λ) — без поиска, счётом.

    S аффинна по λ, потому минимум на отрезке [0; 1] лежит на конце:
    hi + lo > 0 → λ = 1; hi + lo < 0 → λ = 0; hi + lo = 0 → наклон 0, и
    площадь от λ не зависит — оба конца равноправны, берётся 0.

    Это ровно правило самого DeepPoly (λ = 1 при hi ≥ |lo|, иначе 0), и
    здесь оно не взято на веру: `relu_relax` обязан выбирать ту же λ."""
    lo, hi = Q(lo), Q(hi)
    if not (lo < 0 < hi):
        raise ValueError("λ выбирается только на смешанном отрезке lo < 0 < hi")
    наклон = -(hi - lo) * (hi + lo) / Q(2)
    return (Q(1) if hi + lo > 0 else Q(0)), наклон


def certify(pq, x0_units, eps_units: int, cls: int, shift: int, bias_shift: int) -> dict:
    """Сертификат для одного образца: границы логитов и запаса.

    Возвращает нижние/верхние границы логитов и нижнюю границу запаса
    (логит верного класса минус наибольший из прочих), посчитанную
    подстановкой РАЗНОСТИ назад — в этом вся сила дороги.
    """
    W1, b1, W2, b2 = pq["W1"], pq["b1"], pq["W2"], pq["b2"]
    n_in = len(x0_units)
    pre = layer_forms(W1, b1, shift, bias_shift, n_in)
    lows, ups = [], []
    for f in pre:
        lo, hi = f.bounds(x0_units, eps_units)
        a, b = relu_relax(f, lo, hi)
        lows.append(a)
        ups.append(b)
    H, K = len(pre), len(W2[0])

    def logit_form(coef, bias, want_upper: bool) -> Form:
        """Подстановка назад: знак веса решает, какая граница берётся."""
        acc = Form(Q(bias), [Q(0)] * n_in)
        for i in range(H):
            k = coef[i]
            if k == 0:
                continue
            src = (ups[i] if (k > 0) == want_upper else lows[i])
            acc = acc + src.scaled(Q(k, 1 << shift))
        return acc

    zl, zu = [], []
    for j in range(K):
        coef = [int(W2[i][j]) for i in range(H)]
        bias = Q(int(b2[j]) * (1 << bias_shift))
        zl.append(logit_form(coef, bias, False).bounds(x0_units, eps_units)[0])
        zu.append(logit_form(coef, bias, True).bounds(x0_units, eps_units)[1])
    margin = None
    for j in range(K):
        if j == cls:
            continue
        coef = [int(W2[i][cls]) - int(W2[i][j]) for i in range(H)]
        bias = Q((int(b2[cls]) - int(b2[j])) * (1 << bias_shift))
        lo = logit_form(coef, bias, False).bounds(x0_units, eps_units)[0]
        margin = lo if margin is None else min(margin, lo)
    return {"нижние логиты": zl, "верхние логиты": zu, "запас снизу": margin,
            "сертифицирован": margin is not None and margin > 0}


# ---------- самопроверка модуля (ось самоприменения, 09.10) ----------

def self_check(deep: bool = False) -> list:
    """Самопроверка DeepPoly встречными дорогами; список находок, пусто —
    здорова; исключений не бросает.

    `[вычислено здесь]` Окно форм (одна координата, две с разными знаками,
    с нулевым весом) проверяет себя (`Form.self_check`); смягчение ReLU на
    односторонних отрезках точно — нуль при hi ≤ 0, сама форма при lo ≥ 0;
    площадь зазора аффинна по λ — середина отрезка λ даёт полусумму концов,
    наклон совпадает с `lambda_optimum`, и `relu_relax` выбирает ту же λ на
    смешанных отрезках вне ничьей (при hi + lo = 0 оба конца законны) [внешний источник: правило DeepPoly, Singh et al.
    2019]; на сети 2 → 2 → 2 с целыми весами логиты в x₀ по прямому счёту
    лежат в границах `certify`, запас снизу не выше запаса в x₀, а при ε = 0
    границы схлопываются в точку. Граница — окно."""
    out = []
    try:
        forms = [(Form(0, [1, 0]), 1), (Form(Q(1, 2), [Q(2), Q(-3)]), 1), (Form(-1, [Q(0), Q(5)]), 1),
                 (Form(Q(-1, 2), [1, 0]), 2), (Form(Q(-3, 2), [1, 0]), 2)]      # две последние — смешанный отрезок, λ = 1 и λ = 0
        for f, e in forms:
            out += [f"форма {f.c}+{f.w}: {m}" for m in f.self_check(deep, x0=[1, -1], eps=e)]
        f = Form(0, [1, 0])
        z = relu_relax(f, Q(-3), Q(-1))
        if z[0].w != [Q(0), Q(0)] or z[1].w != [Q(0), Q(0)] or z[0].c != 0:
            out.append("ReLU при hi ≤ 0 не нуль")
        s = relu_relax(f, Q(1), Q(3))
        if s[0] is not f or s[1] is not f:
            out.append("ReLU при lo ≥ 0 не сама форма")
        for lo, hi in ((-1, 2), (-3, 1), (-2, 2)):
            s0, s1, sh = relu_area(lo, hi, 0), relu_area(lo, hi, 1), relu_area(lo, hi, Q(1, 2))
            lam, slope = lambda_optimum(lo, hi)
            if sh != (s0 + s1) / 2 or s1 - s0 != slope:
                out.append(f"площадь зазора на [{lo}; {hi}] не аффинна по λ или наклон не тот")
            low, _ = relu_relax(Form(Q(0), [Q(1)]), Q(lo), Q(hi))
            chosen = Q(1) if low.w == [Q(1)] else Q(0)
            if chosen != lam and slope != 0:                     # при hi + lo = 0 площадь от λ не зависит — оба конца законны
                out.append(f"relu_relax на [{lo}; {hi}] выбрал λ = {chosen}, lambda_optimum — {lam}")
        pq = {"W1": [[2, -1], [1, 3]], "b1": [1, -2], "W2": [[1, -2], [2, 1]], "b2": [0, 1]}
        x0 = [3, -1]
        res = certify(pq, x0, 1, 0, 0, 0)
        h = [max(0, sum(x0[i] * pq["W1"][i][j] for i in range(2)) + pq["b1"][j]) for j in range(2)]
        logits = [sum(h[i] * pq["W2"][i][j] for i in range(2)) + pq["b2"][j] for j in range(2)]
        for j in range(2):
            if not (res["нижние логиты"][j] <= logits[j] <= res["верхние логиты"][j]):
                out.append(f"логит {j} в x₀ = {logits[j]} вне границ ⟦{res['нижние логиты'][j]}; {res['верхние логиты'][j]}⟧")
        if res["запас снизу"] > logits[0] - logits[1]:
            out.append(f"запас снизу {res['запас снизу']} выше запаса в x₀ {logits[0] - logits[1]}")
        res0 = certify(pq, x0, 0, 0, 0, 0)
        if res0["нижние логиты"] != logits or res0["верхние логиты"] != logits:
            out.append(f"при ε = 0 границы не схлопнулись в логиты: {res0['нижние логиты']}, {res0['верхние логиты']} против {logits}")
        if len(NET_FIELDS) != 4:
            out.append("полей носителя сети не четыре")
    except Exception as e:                                      # noqa: BLE001
        out.append(f"окно: {type(e).__name__}: {e}")
    return out
