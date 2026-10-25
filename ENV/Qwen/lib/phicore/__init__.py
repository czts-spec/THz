"""phicore — эталонное ядро Ф-чисел (PLAN этап 0→1, линия 00-ЛИНИЯ §3).

Первые принципы: AIAΣA, Ф-число в явной записи, ОП-геометрия (дискретная).
Три правила, вкрученные в устройство, а не в дисциплину:

1. Несведённое действие получает ОТКАЗ (`Unresolved`), не ложный класс.
2. Явная запись `Φ(s, α, M, E | C, параметры)` — единственная форма
   значения в сообщениях и представлениях; перевод — отдельное действие.
3. Ноль — состояние с провенансом; знак и происхождение — компоненты
   записи и не теряются операциями.

Фасад — ДВЕРЬ, не канон: он только выводит имя носителя наружу, поведение
и правила остаются за модулем. Правило двери: у каждого класса каталога
(1, 3–19, 20–23) хотя бы один носитель или фабрика доступны как
`phicore.<имя>`; где у соседей заведена короткая фабрика — выведена и она,
новых имён фасад не изобретает и ничего не переименовывает.

Два места, где у класса два носителя, — оба выведены, канон назван здесь:
  * классы 21/22 — канон контейнера `ScalePointer` (`OMEGA`, `EPSILON`,
    чтение указателя, 00-ЛИНИЯ §2.2); `Hyper` — чтение ЭЛЕМЕНТА с полем
    порядка (`hyper.py`), выбор чтения самого контейнера за автором.
  * класс 19 — `iteration_limit` есть направленное правило (итог живёт в
    классах 1/3), собственный носитель — структурная запись `Mixed`
    (`gridcore.py`, ключ ФIOS (STRUCT, 19)), засечка `Notch` и пара
    `ComplexPair`; `ExpMap(coeff, rate)` — отображение f = coeff·e^{rate·x}
    со своей записью C=19 (`iteration.py`), и при нём второе направленное
    правило `iteration_limit_map` — сертификат сходимости на отрезках.

Класс 17 выведен тремя дверьми: носитель `WValue`, дверь `lambert_w` и
точная дорога `exact_w` (W(k·e^k) = k, в том числе точка ветвления
W(−1/e) = −1 на обеих ветвях) вместе с тождеством `w_exp_w`, которым
двойственность 17 замыкается объектами.

Что за дверью модуля оставлено осознанно: сертификаты сетки. Публичных
имён `*_enclosure` в phicore 16 (замер прогоном), в фасаде 0 — среди них
`iteration.limit_enclosure`, `bessel.j_enclosure`,
`gamma_value.gamma_enclosure`; они вводятся по имени модуля.

`affine` и `exact_gemm` фасад НЕ ввозит: они требуют numpy, и ввоз сделал
бы внешнюю зависимость обязательной для всего пакета (замер с
заблокированным numpy: ввозятся 10 модулей из 12, эти два отказывают
ImportError). Они вводятся по имени: `from phicore import affine`.
"""
from .errors import Unresolved, ExactRepresentationNotAvailable
from .rational import Q
from .record import PhiRecord
from .pointers import ScalePointer, OMEGA, EPSILON
from .zero_state import Zero, ZeroKind
from .intpoly import IntPoly
from .algebraic import AlgebraicNumber
from .enclosure import Interval
from .pi_value import PiMonomial, PI
from .container import (PhiNumberV6, LockedInt, NumberClass,
                        CrossBaseOperationNotSupported)
from .exp_log import ExpValue, LogValue, exp, ln
from .series import PowerSeriesValue
from .dfinite import DFiniteSeries, DFiniteSeq, Recurrence
from .ieee import Ieee
from .window import Window
from .tag import AlphaTable, STATES, pack, unpack, state_of, registers_of
from .product import InfiniteProductValue
from .integral import IntegralValue, integral
from .zeta import ZetaValue, zeta, bernoulli
from .gamma_value import GammaValue, sin_pi
from .elliptic import (KValue, EValue, elliptic_k, elliptic_e,
                       complementary)
from .bessel import BesselRootValue, bessel_root
from .lambert import WValue, lambert_w, exact_w, w_exp_w
from .de_solution import DESolutionValue
from .iteration import (iteration_limit, fixed_point_poly, ExpMap,
                        iteration_limit_map)
from .unity_root import UnityRoot, I
from .arctrig import ArcTrigValue, arctrig
from .hyperbolic import HyperbolicValue, hyperbolic, pythagorean_check
from .alpha import PhiAlpha
from .profile import Profile, ProfileOverflow, LADDER, enforce, fits, project
from .factored import Factored
from . import ccore
from . import rscore
from .lattice import lll, is_reduced, gram_det
from . import lllfactor
from . import vanhoeij
from . import deeppoly
from .zech import Zech, ZechExponent
from .dna import repair, CORRESPONDENCES
from .dispatch import combine

# --- перевод в фасад (00-ЛИНИЯ §3, последний шаг каждого класса) ---------
# Носители, которые прежде требовали ввоза модуля по имени. Фасад открывает
# дверь: ни одно имя не переименовано, ни одна фабрика не заведена заново.
from .gridcore import Notch, Mixed, separate, sign_on_grid
from .complex_pair import ComplexPair
from .hyper import Hyper, laurent_inverse, laurent_div
from .ratfunc import RatFunc, ZeroInQx
from . import cyclotomic
from . import intfactor
from . import wheels
from . import zeros
from . import op_axis
from . import op_cross

#: Всё, что фасад выводит наружу: отсортировано, без дублей;
#: имена модулей-аппаратов (ccore, wheels, …) стоят рядом с
#: именами носителей — дверь у них одна.
__all__ = [
    "AlgebraicNumber", "AlphaTable", "ArcTrigValue", "BesselRootValue",
    "CORRESPONDENCES", "ComplexPair", "CrossBaseOperationNotSupported",
    "DESolutionValue", "DFiniteSeq", "DFiniteSeries", "EPSILON", "EValue",
    "ExactRepresentationNotAvailable", "ExpMap", "ExpValue", "Factored",
    "GammaValue", "Hyper", "HyperbolicValue", "I", "Ieee",
    "InfiniteProductValue", "IntPoly", "IntegralValue", "Interval",
    "KValue", "LADDER", "LockedInt", "LogValue", "Mixed", "Notch",
    "NumberClass", "OMEGA", "PI", "PhiAlpha", "PhiNumberV6", "PhiRecord",
    "PiMonomial", "PowerSeriesValue", "Profile", "ProfileOverflow", "Q",
    "RatFunc", "Recurrence", "STATES", "ScalePointer", "UnityRoot",
    "Unresolved", "WValue", "Window", "Zech", "ZechExponent", "Zero",
    "ZeroInQx", "ZeroKind", "ZetaValue", "arctrig", "bernoulli",
    "bessel_root", "ccore", "combine", "complementary", "cyclotomic",
    "deeppoly", "elliptic_e", "elliptic_k", "enforce", "exact_w", "exp",
    "fits", "fixed_point_poly", "gram_det", "hyperbolic", "integral",
    "intfactor", "is_reduced", "iteration_limit", "iteration_limit_map",
    "lambert_w", "laurent_div", "laurent_inverse", "lll", "lllfactor",
    "ln", "op_axis", "op_cross", "pack", "project", "pythagorean_check",
    "registers_of", "repair", "rscore", "separate", "sign_on_grid",
    "sin_pi", "state_of", "unpack", "vanhoeij", "w_exp_w", "wheels",
    "zeros", "zeta",
]
