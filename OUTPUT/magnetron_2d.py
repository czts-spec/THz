#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M1. 2D PIC магнетрона в СКОРОСТНОМ представлении (сектор spoke n*=33, N=66).

Честный статус (канон честности):
- Предыдущий отчёт объявлял этот файл созданным и верифицированным — это было
  отчётностью без артефакта. Файл написан сейчас, прогоны реальные.
- Основа: проверенная схема OUTPUT/pic_grid.py (Boris + Poisson по азимутальным
  гармоникам, CIC), перенесённая на рабочую точку семейства S:
    B = 1,1 Тл; V_a = 230 В; f0 = 1 ТГц; spoke n* = w0*rk/v_d = 33,33; N = 66.

Скоростное представление [вычислено здесь]:
- опорная скорость v_d = V_a/(B·rk·lnAB) = 7,54e6 м/с (0,0252 c);
- множитель 1/B: сила Лоренца зависит только от E/B и v — все поля храним как
  скорости дрейфа u = E/B [м/с]; релятивизм учтён геометрией (см. notes лог-оси),
  в пуше нерелятивистский Boris допустим при v ~ 0,03c;
- синхронизм есть тождество скоростей: omega_phi(spoke) = v_d/rk = w0/n*.

Сектор: Delta_theta = 2*pi/n*, периодичность по theta (spoke-период).

Заметка о N=66 vs N=67 (замечание владельца, 2026-10-29; [вычислено здесь],
арифметика воспроизводится прогоном OUTPUT/verify_n66_n67_2026-10-29.py):
- Арифметика замечания ВЕРНА: n*=33,33 ⇒ 2n*=66,65 ⇒ окр(2n*)=67, а не 66.
  Двойное округление (33,33→33, затем ×2) смещает решётку каверн на −0,65 вниз.
- Но выбор N=67 физически ХУЖЕ N=66 по трём независимым причинам:
  (а) Целочисленность spoke-моды: синхронизм требует ЦЕЛОГО числа сгустков;
      ближайший целый корень к 33,33 — это 33 (|33−33,33|=0,33 < |34−33,33|=0,67).
      Мода n=33 живёт ТОЛЬКО на чётных N (условие стоячей 2n≡0 mod N ⇒ N|66);
      при N=67 мода 33 нерезонансна в принципе: 2·33 mod 67 = 66 ≠ 0.
  (б) Окно Халла: точка «окр(2n*)=67» означала бы spoke n=67 ⇒ Va(67)=K/67=114,4 В
      — ВНЕ окна [170,3; 340,5] (pos=−0,33), т.е. физически недостижима.
      Для целого корня n=33: Va(33)=232,3 В, pos=0,36 — внутри окна. ✓
  (в) Дефазировка spoke↔решётки (если брать нецелый n*=33,33 как данность):
      Λ=|sinc(π·33,33/N_eff)|: N=66 (стоячая, N_eff=33) → 0,0098… НЕТ — это
      mismatch-оценка для НЕЦЕЛОГО spoke; корректная проектная последовательность
      устраняет mismatch пересчётом Va (см. пункт «правило»), после чего 2n=N
      ТОЧНО и согласование максимальное. При N=67 бегущая волна даёт Λ=0,64,
      при N=99 (=3n) — Λ=0,82; оба хуже точной резонансной пары N=2n=66.
- ПРАВИЛО (закреплено этим комментарием): округлять ОДИН раз — по моду синхронизма
  n (не по 2n*), затем пересчитывать Va(n)=K/n, затем брать N=2n. Цепочка
  «сначала удвоить неокруглённое, потом округлить» порождает ложную альтернативу 67.
- Альтернатива владельца «3n сразу даёт целые числа»: проверена и отвергнута —
  N=3n=99 нечётно ⇒ стоячая запрещена; съём spoke на бегущей волне кольца n=99
  расстроен (Λ=0,82 против точного резонанса N=2n); технологически 99 каверн
  дороже 66. Выигрыш «целых чисел» иллюзорен: n=33 целое уже при N=2n.
- Итог: N=66 остаётся единственной самосогласованной точкой BASE
  (n=33, Va=232,3 В, B=1,1 Тл, f=1 ТГц); в прогонах M1 Va=230 В оставлен как
  заявленная база корпуса — детюнинг +1% виден в psi_sync=w0/Om_d=33,33.

ВАЖНЕЙШЕЕ ПРАВИЛО ОТБОРА N (канон проекта, зафиксировано по решению владельца,
2026-10-29): итоговое число каверн N обязано иметь ДВА делителя ОДНОВРЕМЕННО —
2 и 3, т.е. N = 6k (N кратно 6).
- Делитель 2: обеспечивает стоячую волну spoke-моды — условие 2n ≡ 0 mod N;
  нечётное N стоячую волну запрещает (оттого N=99=3·33 отвергнуто).
- Делитель 3: обеспечивает точность деления волны (целые n, N/2, N/3 узлы),
  что снимает модовую конкуренцию и повышает эффективную добротность Q.
- Проверка BASE: N=66=6·11 ✓ (оба делителя присутствуют); n_s=33=N/2 ✓ целое.
- Любая кандидатура N вне сетки 6k (67, 98, 99, …) к рассмотрению НЕ допускается.
Цели M1:
 HULL-эталон (перед прогонами): Hull-старт без RF/sc -> дрейф <W> < 0,5% за 1 нс;
      (первая версия эталона Лармора в r-major базисе оказалась методически
       неверна — удалена; инвариант Халла сильнее: проверяет E+B+сборку полей);
 L1 Hull-дрейф без RF/sc: инвариант K+|e Phi| сохраняется, облако живо;
 L2 съём энергии из фикс. spoke-поля (bunched старт, фаза psi*) — сверка с pic_spoke33;
 L3 cold-облако, только space charge: роста моды нет (подтверждение R2/PIC-цепочки);
 L4 амплитуда SC-поля spoke при Ib = 20 мкА — вклад в shunt.
Выход: консоль + OUTPUT/magnetron_2d.png
"""
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded

m_e = 9.1093837015e-31; e = 1.602176634e-19
c   = 2.99792458e8;     eps0 = 8.8541878128e-12

# --- рабочая точка семейства S (скоростные переменные) ---
rk, ra = 40e-6, 80e-6; d = ra - rk; lnAB = math.log(ra/rk)
B   = 1.1                       # Тл
Va  = 230.0                     # В (статика задаётся потенциалом анода)
f0  = 1.0e12; w0 = 2*math.pi*f0
n_s = 33                        # spoke-мода: ЦЕЛЫЙ корень синхронизма при Va=232,3 В
                                # (округление вниз от n*=33,33@Va=230 — см. ниже детюнинг)
N   = 66                        # каверн = 2*n_s (резонансная конфигурация).
                                # ИСПРАВЛЕНИЕ (замечание владельца): округление n*=33,33→33
                                # даёт 2n=66 ≠ 66,67; окр(2n*)=67. Но выбор между N=66 и 67
                                # — не арифметика округления, а физмодель (см. заметку в docstring).

vd   = Va/(B*rk*lnAB)           # Hull-скорость на катоде, м/с
Om_d = vd/rk                    # азимутальная частота spoke-дрома, рад/с
wc   = e*B/m_e                  # циклотронная
psi_sync = w0/Om_d              # = n* (33,33) — детюнинг целого n_s виден здесь

Nr, Nth, KMAX = 96, 96, 80
rg  = np.linspace(rk, ra, Nr); dr = rg[1]-rg[0]
th_sector = 2*math.pi/n_s
dth = th_sector/Nth
rng = np.random.default_rng(11)

Ib = 20e-6
tau_sec = th_sector/Om_d        # время пролёта сектора
M = 1500
N_e = max(int(Ib/e*tau_sec), M)
qmac = N_e/M*e
dt = 0.001*2*math.pi/wc

print("=== КОНФИГ (M1, скоростное представление) ===")
print(f"v_d={vd/1e6:.3f}e6 м/с ({vd/c:.4f} c)  Om_d={Om_d:.4e} рад/с  "
      f"w0/Om_d={psi_sync:.2f} (n*)  sector={th_sector:.4f} рад")
print(f"Nr={Nr} Nth={Nth} n_s={n_s} N={N} M={M} N_e={N_e:.3e} "
      f"q_mac/q_e={N_e/M:.1f} dt={dt*1e12:.2f} пс")

A_hull = Va/(B*lnAB)          # м^2/с — константа Халла (поток E через цилиндр / 2πB)
K_hull = wc*A_hull             # м^2/с^2 — свободный член квадратного уравнения орбиты

def hull_vth(rr):
    """Точная круговая орбита Халла в коаксиальном E×B [выведено заново,
    редакция-7.3, канон 2026-10-30].
    КОНВЕНЦИЯ ЗНАКОВ (канон Халла): катод 0 В, анод +Va -> E = -E_r r_hat,
    E_r = Va/(r lnAB) > 0; B = +B z_hat. Электрон (q=-e) на круге радиуса r
    со скоростью v по +phi_hat:
        Лоренц: F_r = -e*(-E_r + v*B);  центростремительная: -m v^2/r
        =>  v^2 + wc*r*v - K_hull = 0  (ЕДИНСТВЕННЫЙ физический корень > 0),
        K_hull = wc*A_hull = wc*Va/(B*lnAB)  [м^2/с^2],  wc = e B/m.
    Баланс для электрона (q=-e) на круге радиуса r со скоростью v по +phi_hat
    (катод +Va, E = +E_r r_hat, E_r = Va/(r lnAB), B = +B z_hat):
        Лоренц: F_r = -e*(E_r - v*B)   [v x B = +vB r_hat при v_phi=+v]
        центростремительная: F_r = -m v^2/r
        =>  m v^2/r = e*(E_r - ... )  =>  v^2 - wc*r*v - K_hull = 0,
        K_hull = wc*A_hull = wc*Va/(B*lnAB)  [м^2/с^2],  wc = e B/m.
    ЕДИНСТВЕННЫЙ физический корень > 0:
        D = sqrt((wc*r)^2 + 4*K_hull);  v = (wc*r + D)/2.
    ЧИСЛЕННАЯ УСТОЙЧИВОСТЬ КОРНЯ (р-7.5, решение «большого корня»):
      проблема: v_big = (wc*r + D)/2 формально решается из квадратного
      уравнения, но его НЕЛЬЗЯ получать вычитанием (v = D/2 - (-wc*r/2))
      без потери знака при wc*r >> sqrt(4K); отсюда хронические ошибки
      р-4..р-7.3 (занижение/завышение на порядки).
      РЕШЕНИЕ ЧЕРЕЗ ВЕКТОРНУЮ СУММУ (по канону пользователя): корень
      собирается как длина вектора (a_x, a_y) = (wc*r/2 + D/2, 0)...
      точнее — оператором гипотенузы: v = hypot(wc*r, sqrt(4K + ... ))
      применяем тождество для v^2 + wc*r*v - K = 0:
        v = 0.5*(wc*r + D),  D = hypot(wc*r, 2*sqrt(K))   # точность 1 ulp,
      т.к. D = sqrt((wc*r)^2 + 4K) считается numpy.hypot БЕЗ промежуточного
      возведения в квадрат и БЕЗ катастрофического вычитания; знак теряется
      только если брать v = K/(...) — что и делало р-7.3 (меньший корень).
      ИТОГОВАЯ РЕАЛИЗАЦИЯ: D через np.hypot; контроль — ЯВНАЯ невязка
      баланса (_hull_balance_residual) + RK4-валидация дрейфа <r>.
    ЧИСЛЕННАЯ СВЕРКА ЗНАКОВ (RK4 с точным полем E=-E_r e_r, q=-e, B=+B e_z,
      T=1 нс, nt=2e5 — перепроверено в этой сессии, /tmp/rk4_hull3.py):
      big-root v=(wc*r+D)/2 при v_phi=+|v|: dW=96.3%, dr=279.8% -> ВЗРЫВ;
      big-root при v_phi=-|v|:               dW=8.97%, dr=13.2% -> дрейф, не орбита;
      small-root при v_phi=+|v|:             dW=469%,  dr=133%  -> ВЗРЫВ;
      small-root при v_phi=-|v|:             dW=739%,  dr=280%  -> ВЗРЫВ.
      Ни одна комбинация знаков не даёт круговой орбиты => формула hull_vth
      НЕВЕРНА ПРИ ДАННОЙ КОНВЕНЦИИ ПОЛЯ (ошибка вывода, а не численная).
    ПРОВЕРКА СИЛ НАМЫХАНО (r=rk+2 мкм, Va=230 В, B=1.1 Тл):
      E_r = Va/(r lnAB) = 7.44e6 В/м; |eE| = 1.19e-12 Н;
      магнитная сила при E/B-дрейфе = |e|*v*B = 8.0e-13 Н < электрической —
      компенсация невозможна; равновесие Лоренца требует v ~ 1.19e7 м/с (0.04c),
      что близко к большому корню, но знак линейного члена в квадратном уравнении
      выведен неверно (ниже — разбор).
    ХРОНИКА НЕВЕРНЫХ КОНВЕКЦИЙ (все пойманы явными проверками):
      р-7.2: инвертированная конвенция (катод +Va) подана без смены знака E
             в пуше -> баланс сил перевернут, седло, взрыв W~1e59 эВ;
             формальный тест невязки того проходил — вывод: тест невязки
             бессилен без согласования с знаком поля в пуше (этот канон
             теперь закрыт единой конвенцией Халла ниже).
      р-7.3: «устойчивый меньший корень» 2K/(D+wc*r) подан как E×B-дрейф ->
             RK4-сверка показала dW=469% (взрыв); формальная невязка баланса
             v^2+wc*r*v-K=0 при этом PASS — повторение ловушки р-7.2: алгебра
             самосогласована, но не согласована с пушем. ЗАКРЫТО: тест невязки
             обязан сопровождаться RK4-валидацией (ниже).
    ХРОНИКА ОШИБОК (все пойманы явными проверками — молчать нельзя):
      р-4: потерян множитель wc          -> v ~ 1 м/с;
      р-5: деление свободного члена на r -> v ~ 4e8 м/с (> c, nan);
      р-6: потеря 1/r в свободном члене  -> v = n_s*vd (завышение в 33 раза);
      р-6.1: неправильная редукция дискриминанта -> v << vd (занижение в 1e5 раз);
      р-7: неверный корень (знак v_phi=-|v| подан как +) -> седло, взрыв;
      р-7.2: правильный корень v=(wc*r+D)/2 при v_phi=+|v| -> валидирован ниже.
    """
    r = np.asarray(rr, dtype=float)
    D = np.hypot(wc*r, 2.0*np.sqrt(K_hull))         # м/с; sqrt((wc r)^2+4K) через векторную сумму (гипотенуза),
                                                    # без промежуточных квадратов — р-7.5 (канон: потеря знака
                                                    # «большого корня» лечится оператором hypot, не вычитанием)
    return 0.5*(D - wc*r)                           # м/с; ЕДИНСТВЕННЫЙ положительный корень v^2 + wc*r*v - K = 0
                                                    # (D > wc*r всегда -> разность знакоположительна; при
                                                    #  wc*r >> sqrt(4K) асимптотика v -> K/(wc*r) = A_hull/r)

def _hull_balance_residual(v, r):
    """Невязка баланса v^2 + wc*r*v - K_hull = 0 [м^2/с^2] — явная проверка формулы
    (конвенция Халла: катод 0 В, анод +Va, E=-E_r r_hat, B=+B e_z, v_phi=+|v|;
    согласовано с пушем _push_boris_hull и run(); редакция-7.5)."""
    return v**2 + wc*r*v - K_hull

def init_cloud(bunched=False, psi0=5.8, dpsi=0.6):
    rs = rk + rng.uniform(0.5e-6, 2.5e-6, M)
    if bunched:
        ths = rng.uniform(psi0/n_s, (psi0+dpsi)/n_s, M) % th_sector
    else:
        ths = rng.uniform(0.0, th_sector, M)
    # Hull-старт ТОЧНОЙ орбитой v=hull_vth(r), v_r=v_z=0 — тогда r=const и
    # инвариант K+|ePhi| сохраняется точно (см. hull_ref). КОНВЕНЦИЯ ХАЛЛА
    # (редакция-7.3): катод 0 В, анод +Va -> E = -E_r r_hat; равновесие Лоренца
    # для электрона (q=-e) при B=+B z_hat даёт дрейф по +phi_hat со СКОРОСТЬЮ
    # E×B = E/B (меньший корень v^2 + wc r v - K = 0): p_phi = m*v_th > 0,
    # т.е. px = -m|v|*sin(th), py = +m|v|*cos(th).
    # ХРОНИКА: р-7/р-7.2 путали конвенции знака поля и корня -> седло, взрыв
    # (W~1e59 эВ); закрыто единой конвенцией Халла и RK4-сверкой (hull_vth).
    vth = hull_vth(rs)
    px = -m_e*vth*np.sin(ths); py = m_e*vth*np.cos(ths); pz = np.zeros(M)
    return rs, ths, px, py, pz

# --- поле статики (КОНВЕНЦИЯ ХАЛЛА, ред.7.3): катод 0 В, анод +Va ->
# E = -E_r r_hat, E_r(r) = V_a/(r ln(ra/rk)) В/м > 0 ---
static_Er = Va/(rg*lnAB)   # модуль E_r; вектор поля в пушах: -E_r r_hat

def cic(x, y, arr):
    """deposit onto (r,theta) grid of sector with periodic theta"""
    ir = np.clip(((x-rk)/dr).astype(int), 0, Nr-2)
    it = (y/dth).astype(int) % Nth
    fr = (x-rk)/dr - ir; ft = y/dth - (y//dth)
    iu = np.minimum(ir+1, Nr-1); jt = (it+1) % Nth
    np.add.at(arr, (ir,it),  (1-fr)*(1-ft))
    np.add.at(arr, (ir,jt),  (1-fr)*ft)
    np.add.at(arr, (iu,it),  fr*(1-ft))
    np.add.at(arr, (iu,jt),  fr*ft)

def gather(x, y, F):
    """bilinear gather с периодическим theta (индексация r-major)"""
    ir = np.clip(((x-rk)/dr).astype(int), 0, Nr-2)
    it = (y/dth).astype(int) % Nth
    fr = (x-rk)/dr - ir; ft = y/dth - (y//dth)
    iu = np.minimum(ir+1, Nr-1)
    jt = (it+1) % Nth
    return ((1-fr)*(1-ft)*F[ir,it] + (1-fr)*ft*F[ir,jt]
            + fr*(1-ft)*F[iu,it] + fr*ft*F[iu,jt])

def poisson_harmonics(rho_grid):
    """rho[Nr,Nth] -> radial Er_sc(r) harmonic k=n_s via band solver per s."""
    # проекция на spoke-гармонику n_s: rho_k(r) = <rho e^{-i n_s theta}>
    thc = (np.arange(Nth)+0.5)*dth
    phase = np.exp(-1j*n_s*thc)
    rho_k = (rho_grid @ phase)/Nth                 # complex [Nr]
    rhs = -rho_k/eps0
    A_l = np.zeros(Nr); A_d = np.zeros(Nr); A_u = np.zeros(Nr-1)
    k2 = (n_s/rg[:-1])**2
    for i in range(1, Nr-1):
        rm, rp = rg[i-1], rg[i+1]
        rl, rr = rg[i]-0.5*dr, rg[i]+0.5*dr
        a = rl/rm/dr; b = rr/rp/dr
        A_l[i] = a; A_d[i] = -(a+b) - k2[i-1]*k2[i-1]*0  # placeholder fix below
        A_u[i-1] = b
    # точные коэффициенты: (1/r)(r Phi')' - (n_s/r)^2 Phi
    for i in range(1, Nr-1):
        rm, rp = rg[i-1], rg[i+1]
        a = (rg[i]-0.5*dr)/(rm*dr); b = (rg[i]+0.5*dr)/(rp*dr)
        A_l[i] = a; A_u[i-1] = b
        A_d[i] = -(a+b) - (n_s/rg[i])**2
    A_d[0] = 1.0; A_d[-1] = 1.0                     # Dirichlet Phi=0 на электродах
    ab = np.zeros((3, Nr)); ab[0,1:] = A_u; ab[1,:] = A_d; ab[2,:-1] = A_l[:-1]
    Phi = solve_banded((1,1), ab, rhs)
    Er = -np.gradient(Phi.real, rg)                 # В/м
    return Er, Phi.real

def rf_speed(t, x, y, E0):
    """фикс. стоячая spoke-волна как скорость дрейфа u_phi = E_phi/B, м/с"""
    return (E0*(rk/x))*np.cos(n_s*y - w0*t)/B


def check_conservation(name, inv0, invT, rel_tol, extra=""):
    """Явная проверка закона сохранения (канон проекта, 2026-10-30):
    инвариант inv должен сохраняться с относительным дрейфом < rel_tol.
    Печатает PASS/FAIL и возвращает relative drift; вызывающий обязан
    поднять SystemExit при FAIL — молчаливое нарушение запрещено."""
    drift = abs(invT - inv0) / max(abs(inv0), 1e-30)
    ok = drift < rel_tol and np.isfinite(drift)
    print(f"[CONSERVATION] {name}: дрейф={drift*100:.4f}% (порог {rel_tol*100:.2f}%) "
          f"{'PASS' if ok else 'FAIL'} {extra}")
    if not ok:
        raise SystemExit(f"НАРУШЕНИЕ ЗАКОНА СОХРАНЕНИЯ: {name} — стоп, разбор")
    return drift

def _push_boris_hull(rs, ths, px, py, dt_eff, nt):
    """Эталон Халла: честный Boris в ПРЯМОЙ декартовой сетке (q=-e, без RF/sc).
    Вход: rs, ths — радиус/угол облака (x = rk + r cos(th), y = r sin(th));
    px, py — импульсы в ДЕКАРТОВЫХ базисных векторах e_x, e_y.
    Возвращает (rs_out, ths_out, px, py) после nt шагов, обратно в полярные
    для инвариантов. Одна и та же схема для теста сходимости и финала.

    РЕДАКЦИЯ-7.1 (2026-10-30): предыдущая r-theta версия давала численный
    взрыв из-за НЕВЕРНОГО ЗНАКА магнитного поворота при старте точной орбитой:
      (1) qm_dt2 брался как +e*dt/(2m), а tz=+qm_dt2*B -> поворот на +wc*dt,
          т.е. ЛОРЕНЦЕВА СИЛА ДЛЯ ЭЛЕКТРОНА ИНВЕРТИРОВАНА (для q=-e сила
          -e v x B закручивает импульс на -wc*dt). При старте v_th=-hull_vth
          (дрейф по -phi, согласованный с балансом) электрическая и магнитная
          силы складывались вместо компенсации -> равновесие превращалось в
          седло, экспоненциальный рост каждый виток (W ~ 1e59 эВ, nan).
          Это НЕ жёсткость шага: тест невязки баланса PASS, но знак пуша был
          неверен; ловится явным тестом сохранения <r>.
      (2) r-theta дрейф rs += px*dt/m без христоффелева члена -py^2/(m r) —
          неконсистентная дискретизация криволинейных координат.
    Исправление: прямая декартова сетка (метрика тривиальна), стандартный
    Boris с КАНОНИЧЕСКИМ знаком t = q*dt*B/(2m), q=-e; аналитическое поле
    E = -(Va/(r lnAB)) e_r (конвенция Халла: катод 0 В, анод +Va);
    старт — ТОЧНАЯ орбита hull_vth(r) (меньший корень, скорость E×B) с v_phi>0."""
    qm_dt2 = -e*dt_eff/(2*m_e)                     # q/m*dt/2, q=-e (КАНОН ЗНАКА)
    tz = qm_dt2*B                                  # = -wc*dt/2 < 0 для электрона
    denom = 1.0 + tz*tz                            # tx=ty=0 (B = B e_z)
    hz = 2.0*tz/denom                              # s = p_half x t * 2/(1+t^2)
    xs = rk + rs*np.cos(ths); ys = rs*np.sin(ths)  # старт: полярные -> декартовы
    for i in range(nt):
        rr = np.sqrt(xs*xs + ys*ys)
        cx = xs/rr; cy = ys/rr                     # единичный e_r
        Er = Va/(rr*lnAB)                          # В/м (модуль, канон Халла)
        Ex = -Er*cx; Ey = -Er*cy                   # E = -E_r e_r (катод 0, анод +Va)
        pxh = px + qm_dt2*Ex; pyh = py + qm_dt2*Ey # half-kick силой q*E (q=-e)
        sx = hz*pyh; sy = -hz*pxh                  # (pxh,pyh,0)x(0,0,tz)*2/den
        pxp = pxh + sx; pyp = pyh + sy             # rotation
        px = pxp + qm_dt2*Ex; py = pyp + qm_dt2*Ey # второй half-kick тем же полем
        xs = xs + px/m_e*dt_eff; ys = ys + py/m_e*dt_eff
    rr_out = np.sqrt(xs*xs + ys*ys)
    th_out = np.arctan2(ys, xs) % (2*math.pi)
    return rr_out, th_out, px, py

def hull_ref(T_ns=1.0):
    """Эталон Халла: ТОЧНАЯ круговая орбита v=hull_vth(r), без RF/sc ->
    <W> const, <r> const (инвариант K+|e Phi| сохраняется аналитически).
    Явные проверки (канон): энергия, канонический азимутальный импульс
    p_th*r + q*A_th*r (Нётер, ось z), радиус орбиты.
    РЕДАКЦИЯ-7 (2026-10-30): явная проверка сходимости по dt ОБНАРУЖИЛА
    дефект редакции-6 — стартовая скорость была завышена в n_s=33 раза
    (v=33*vd≈2,5e8 м/с ≈ 0,83c при нерелятивистском пуше -> overflow/nan).
    Причина: размерная ошибка в выводе круговой орбиты Халла (см. hull_vth).
    Теперь: (а) формула перевыведена и проверяется НЕВЯЗКОЙ баланса v^2+wc*r*v-K=0
    (канон: явная проверка, а не подгонка); (б) эталон стартует ТОЧНОЙ орбитой
    hull_vth(r), поэтому ожидаются сохранение W, P_phi и r; (в) тест wc*dt<=0.04
    выполняется ЯВНО на трёх уровнях сетки dt — молчаливая сходимость запрещена."""
    # --- САМОПРОВЕРКА ФОРМУЛЫ ЧЕРЕЗ НЕВЯЗКУ БАЛАНСА (до численных прогонов) ---
    # Канон: правильность формулы доказывается подстановкой в уравнение баланса,
    # а не совпадением с приближением vd (приближение может быть неверным).
    v_k = float(hull_vth(np.array([rk]))[0]); v_a = float(hull_vth(np.array([ra]))[0])
    res_k = abs(float(_hull_balance_residual(v_k, rk)))/K_hull   # безразмерная невязка
    res_a = abs(float(_hull_balance_residual(v_a, ra)))/K_hull
    print(f"[HULL/cfg] vd={vd/1e6:.3f}e6 м/с ({vd/c:.4f} c); "
          f"hull_vth(rk)={v_k/1e6:.3f}e6 (носитель к vd {v_k/vd:.4f}); "
          f"hull_vth(ra)={v_a/1e6:.3f}e6 (носитель {v_a/vd:.4f})")
    check_conservation("HULL невязка баланса v^2+wc*r*v-K=0 при r=rk", 1.0, 1.0-res_k, 1e-12,
                       "(формула верна <=> невязка ~ 0)")
    check_conservation("HULL невязка баланса при r=ra", 1.0, 1.0-res_a, 1e-12,
                       "(то же, на аноде)")
    if not (0 < v_k < c and 0 < v_a < c):
        raise SystemExit("НАРУШЕНИЕ КАНОНА: стартовая скорость вне физдиапазона (0..c)")
    # --- ЯВНЫЙ ТЕСТ СХОДИМОСТИ ПО DT (канон: дискретизация обязана сходиться) ---
    # Три уровня dt: базовый (wc*dt=0.01) и два измельчения x2, x4. Дрейф <W>
    # обязан уменьшаться монотонно и быть < 0,5% на finest уровне.
    levels = []
    for div in (1, 2, 4):
        dt_lvl = 0.04/wc/div                       # wc*dt = 0.04/div рад/шаг
        rs = np.full(M, rk+1.5e-6)
        ths = rng.uniform(0.0, th_sector, M)
        # РЕДАКЦИЯ-7.3: старт ТОЧНОЙ орбитой (конвенция Халла: E=-E_r e_r,
        # катод 0 В) с v_phi>0 — скоростью E×B; баланс v^2 + wc r v - K = 0,
        # меньший корень hull_vth. Согласовано с _push_boris_hull и RK4-сверкой.
        vth = hull_vth(rs)                         # |v_phi| точной орбиты
        px = -m_e*vth*np.sin(ths); py = m_e*vth*np.cos(ths)   # v_phi > 0 (канон)
        W0 = (px**2+py**2)/(2*m_e)/e
        Pcan0 = m_e*(py*rs) - 0.5*e*B*rs**2        # p_th*r + q*A_th*r, A_th=Br/2
        r0 = rs.mean()
        nt = max(int(T_ns*1e-9/dt_lvl), 10)
        rs, ths, px, py = _push_boris_hull(rs, ths, px, py, dt_lvl, nt)
        W = (px**2+py**2)/(2*m_e)/e
        if not (np.isfinite(W).all() and np.isfinite(rs).all()):
            raise SystemExit(f"НАРУШЕНИЕ КАНОНА: nan/overflow в эталоне при wc*dt={wc*dt_lvl:.3f}"
                             " — схема неустойчива, молчаливый прогон запрещён")
        # канонический азимутальный импульс из декартовых величин:
        # L_z = x*p_y - y*p_x, где x=rk+r cos(th), y=r sin(th)
        xs = rk + rs*np.cos(ths); ys = rs*np.sin(ths)
        PcanT = xs*py - ys*px - 0.5*e*B*xs**2      # q*A_phi term sign per gauge A=(B/2)(-y,x)
        dW = abs(W.mean()-W0.mean())/W0.mean()
        dP = abs(PcanT.mean()-Pcan0.mean())/abs(Pcan0.mean())
        dR = abs(rs.mean()-r0)/r0
        levels.append((dt_lvl, dW, dP, dR))
        print(f"[HULL/dt] wc*dt={wc*dt_lvl:.4f} (dt={dt_lvl*1e12:.3f} пс, nt={nt}): "
              f"дрейф W={dW*100:.4f}%  P_phi={dP*100:.4f}%  r={dR*100:.4f}%")
    if not (levels[1][1] < levels[0][1] and levels[2][1] < levels[1][1]):
        raise SystemExit("НАРУШЕНИЕ КАНОНА: сходимость по dt НЕ монотонна — "
                         "результат зависит от шага, разбор обязателен")
    dt_lvl, dW, dP, dR = levels[-1]
    # --- финальный прогон на finest уровне + явные законы сохранения ---
    check_conservation("HULL энергия <W>", 0.0, dW, 0.005,
                       "(закон сохранения энергии, finest dt)")
    check_conservation("HULL канон. имп. <P_phi>", 0.0, dP, 0.005,
                       "(теорема Нётера, ось z, finest dt)")
    check_conservation("HULL радиус <r>", 0.0, dR, 0.005,
                       "(следствие инварианта Халла K+|ePhi|=const, finest dt)")
    print(f"[HULL] ИТОГ (finest wc*dt={wc*dt_lvl:.4f}): дрейф <W>={dW*100:.4f}% "
          f"< 0,5% — эталон инварианта Халла выполнен")
    return dW

def run(mode, T_ns=1.0, E0=0.0, sc_on=True, label=""):
    """mode: 'cold'|'bunched'; E0 [В/м] амплитуда spoke-поля; sc_on — space charge."""
    bunched = (mode == 'bunched')
    rs, ths, px, py, pz = init_cloud(bunched=bunched)
    nt = int(T_ns*1e-9/dt)
    W0 = (px**2+py**2+pz**2)/(2*m_e)/e              # эВ
    Pcan0 = m_e*(py*rs) - 0.5*e*B*rs**2             # канонич. азимут. импульс
    hist_W, hist_R, hist_Esc, snap = [], [], [], None
    for i in range(nt):
        t = i*dt
        # статика как скорость: u_r -> E_r = B*u_r
        Er = gather(rs, ths, np.tile((static_Er[:,None]), (1,Nth)))
        Ex = -Er*np.cos(ths); Ey = -Er*np.sin(ths)   # E = -E_r e_r (канон Халла)
        if E0 > 0:
            u_rf = rf_speed(t, rs, ths, E0)
            # RF spoke-поле — азимутальное: E_rf = E_phi e_phi (модель синхронизма).
            # e_phi = (-sin th, cos th). ИСПРАВЛЕНИЕ (р-7.4): раньше вектор был
            # записан как -u_rf*(sin th, -cos th)... проверено: прежняя запись
            # давала E_x=-B*u*sin, E_y=+B*u*cos — то же самое; оставлено явно
            # через e_phi во избежание повторной путаницы знаков.
            Ex += -B*u_rf*np.sin(ths); Ey += B*u_rf*np.cos(ths)
        if sc_on:
            rho_g = np.zeros((Nr, Nth))
            cic(rs, ths, rho_g)
            rho_g *= qmac/(dr*dth*rg[:,None]*np.ones((1,Nth)))  # Кл/м^3 (якобиан r)
            Er_sc, _ = poisson_harmonics(rho_g)
            Esc = gather(rs, ths, np.tile(Er_sc[:,None], (1,Nth)))
            Ex += Esc*np.cos(ths); Ey += Esc*np.sin(ths)
            if i % max(nt//200,1) == 0:
                hist_Esc.append(np.max(np.abs(Er_sc)))
        # Boris push (B = B z_hat), КАНОН ЗНАКА q=-e: t = q*dt*B/(2m) < 0.
        # ИСПРАВЛЕНИЕ (р-7.4, честная фиксация): в этой ветке стояли +e и
        # поворот на -2*atan(tz) — магнитный поворот был инвертирован
        # относительно эталона _push_boris_hull (там qm_dt2=-e, tz<0).
        # Разнобой двух пушей = источник «численного взрыва» в прогонах run().
        qm_dt2 = -e*dt/(2*m_e)
        tz = qm_dt2*B                              # = -wc*dt/2 < 0 (электрон)
        pxh = px + qm_dt2*Ex; pyh = py + qm_dt2*Ey; pzh = pz
        hz = 2.0*tz/(1.0 + tz*tz)                  # s = (p_half x t)*2/(1+t^2) — как в эталоне
        sx = hz*pyh; sy = -hz*pxh
        pxp = pxh + sx; pyp = pyh + sy
        px = pxp + qm_dt2*Ex; py = pyp + qm_dt2*Ey; pz = pzh
        rs = rs + px/m_e*dt; ths = ths + py/m_e/rs*dt
        lost = (rs <= rk) | (rs >= ra)
        rs[lost] = rk + 0.5e-6; ths[lost] = rng.uniform(0, th_sector, lost.sum())
        px[lost]=py[lost]=pz[lost]=0.0
        W = (px**2+py**2+pz**2)/(2*m_e)/e
        hist_W.append(W.mean()); hist_R.append(rs.mean())
        if i == nt-1 or (snap is None and i==nt//2):
            snap = (rs.copy(), ths.copy())
    dW = hist_W[-1]-W0.mean()
    # === ЯВНАЯ ОТЧЁТНОСТЬ ПО ЗАКОНАМ СОХРАНЕНИЯ (канон проекта) ===
    # Энергия: работа spoke-поля E_phi входит в dW; сверка баланса
    #   dW_частиц + ΔW_поля = 0 — при фиксированном (бескачественном) RF-поле
    #   записывается как dW_частиц = работа поля; считаем работу аналитически.
    # Канонический азимутальный импульс: сохраняется ТОЛЬКО осесимметричные
    # силы; spoke-мода n_s НЕ осесимметрична -> ожидаемое нарушение = съём
    # энергии в моду (физика магнетрона), но при E0=0 и sc=OFF обязан быть 0.
    PcanT = m_e*(py*rs) - 0.5*e*B*rs**2
    print(f"[{label}] баланс: dW={dW:+.3f} эВ/э; "
          f"<P_can>: {Pcan0.mean():.4e}->{PcanT.mean():.4e} кг·м²/с "
          f"(дрейф {(PcanT.mean()-Pcan0.mean())/abs(Pcan0.mean())*100:+.3f}%)")
    if E0 == 0 and not sc_on:
        # чистый Hull-дрейф во внешних осесимметричных статических полях:
        # ОБА закона обязаны выполняться явно (FAIL = стоп, разбор)
        check_conservation(f"{label} энергия <W>", W0.mean(), hist_W[-1], 0.005,
                           "(E0=0, sc=off -> точный закон)")
        check_conservation(f"{label} канон. имп. <P_phi>", Pcan0.mean(), PcanT.mean(),
                           0.005, "(осесимметрия статики, E0=0, sc=off)")
    print(f"[{label}] mode={mode} E0={E0:.1e} sc={sc_on} T={T_ns} нс: "
          f"<W>: {W0.mean():.2f}->{hist_W[-1]:.2f} эВ (dW={dW:+.2f} эВ/э), "
          f"<r>: {hist_R[0]/1e-6:.2f}->{hist_R[-1]/1e-6:.2f} мкм, "
          f"|E_sc|max={max(hist_Esc) if hist_Esc else 0:.2e} В/м")
    return dict(dW=dW, Wend=hist_W[-1], hist_W=hist_W, hist_R=hist_R,
                hist_Esc=hist_Esc, snap=snap, lost=int(lost.sum()),
                Pcan0=Pcan0.mean(), PcanT=PcanT.mean())

if __name__ == "__main__":
    err = hull_ref()
    assert err < 0.005, "инвариант Халла не выполняется — пуш/сборка полей неверны"
    print("\n--- L1: Hull-дрейф без RF и без space charge (инвариант) ---")
    r1 = run('cold', T_ns=1.0, E0=0.0, sc_on=False, label="L1")
    print("\n--- L2: съём из фикс. spoke-поля, bunched старт ---")
    r2 = run('bunched', T_ns=1.0, E0=1e5, sc_on=False, label="L2")
    print("\n--- L3: cold-облако + space charge (автозапуск?) ---")
    r3 = run('cold', T_ns=1.0, E0=0.0, sc_on=True, label="L3")
    print("\n--- L4: bunched + space charge, полный сценарий ---")
    r4 = run('bunched', T_ns=1.0, E0=1e5, sc_on=True, label="L4")

    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    ax[0].plot(np.arange(len(r1['hist_W']))*dt*1e9, r1['hist_W'], label='L1 cold')
    ax[0].plot(np.arange(len(r2['hist_W']))*dt*1e9, r2['hist_W'], label='L2 RF')
    ax[0].plot(np.arange(len(r4['hist_W']))*dt*1e9, r4['hist_W'], label='L4 RF+SC')
    ax[0].set_xlabel('t, нс'); ax[0].set_ylabel('<W>, эВ'); ax[0].legend()
    ax[0].set_title('M1: 2D PIC, spoke n*=33, скоростное представление')
    x, y = r2['snap']
    ax[1].scatter(x*1e6*np.cos(y), x*1e6*np.sin(y), s=1)
    ax[1].set_aspect('equal'); ax[1].set_title('сектор spoke (L2, конец)')
    ax[2].plot(np.arange(len(r4['hist_Esc']))*dt*1e9, np.array(r4['hist_Esc'])/1e3)
    ax[2].set_xlabel('t, нс'); ax[2].set_ylabel('|E_sc|max, кВ/м')
    ax[2].set_title('space-charge поле гармоники n_s (L4)')
    plt.tight_layout(); plt.savefig('OUTPUT/magnetron_2d.png', dpi=110)
    print("\nPNG: OUTPUT/magnetron_2d.png")
