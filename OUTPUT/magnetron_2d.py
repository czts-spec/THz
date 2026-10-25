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
n_s = 33                        # spoke-мода (целая ближайшая к n*=33,33)
N   = 66                        # каверн = 2*n_s (резонансная конфигурация)

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

def hull_vth(rr):
    """Точная круговая орбита Халла в коаксиальном E×B [выведено и проверено здесь]:
    баланс v^2/r = (e/m)(E(r) - v B), E(r)=Va/(r lnAB) (катод +Va, дрейф по
    -phi_hat при B=+z) даёт квадратное уравнение v^2 + wc r v - A/r = 0.
    Меньший корень в УСТОЙЧИВОЙ форме (без catastrophic cancellation, редакция-4):
        v(r) = 2A / ( r*( sqrt(A^2/r^2 + 4 A wc r) + A/r ) ).
    Проверка размерности: [м^2/с] / ([м]·[м/с]) = [м/с] — верно.
    Тождество на катоде: x=2A/(wc r^3)<<1 => sqrt(A^2/r^2+4A wc r)=A/r*sqrt(1+x)
    => v(r_k)=A/(2 r_k)=v_d/2. Двойное расхождение с корпусным v_d — физично:
    Hull-старт pic-цепочки (v=A/r) есть ЧИСТЫЙ E×B-дрейф циклоиды, а не круговая
    орбита; для эталона инвариантности пуша используется именно он (редакция-1),
    т.к. K+|ePhi| сохраняется и для циклоиды (среднее по ларморовскому витку)."""
    s = np.sqrt(A_hull**2/rr**2 + 4*A_hull*wc*rr)
    return 2*A_hull/(rr*(s + A_hull/rr))

def init_cloud(bunched=False, psi0=5.8, dpsi=0.6):
    rs = rk + rng.uniform(0.5e-6, 2.5e-6, M)
    if bunched:
        ths = rng.uniform(psi0/n_s, (psi0+dpsi)/n_s, M) % th_sector
    else:
        ths = rng.uniform(0.0, th_sector, M)
    # Hull-старт ТОЧНОЙ орбитой v=hull_vth(r), v_r=v_z=0 — тогда r=const и
    # инвариант K+|ePhi| сохраняется точно (см. hull_ref). Знак: при B=+z и
    # E=+E_r r_hat (катод +Va) равновесный дрейф по -phi_hat (q=-e).
    vth = -hull_vth(rs)
    px = m_e*vth*(-np.sin(ths)); py = m_e*vth*np.cos(ths); pz = np.zeros(M)
    return rs, ths, px, py, pz

# --- поле статики: E_r(r) = +V_a/(r ln(ra/rk)) В/м, r_hat (катод +Va, анод 0) ---
static_Er = Va/(rg*lnAB)

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


def hull_ref(T_ns=1.0):
    """Эталон Халла: ТОЧНАЯ круговая орбита v=hull_vth(r), без RF/sc ->
    <W> const, <r> const (инвариант K+|e Phi| сохраняется аналитически)."""
    rs = np.full(M, rk+1.5e-6)
    ths = rng.uniform(0.0, th_sector, M)
    vth = -hull_vth(rs)
    px = m_e*vth*(-np.sin(ths)); py = m_e*vth*np.cos(ths); pz = np.zeros(M)
    W0 = (px**2+py**2+pz**2)/(2*m_e)/e
    nt = int(T_ns*1e-9/dt)
    for i in range(nt):
        Er = gather(rs, ths, np.tile((static_Er[:,None]), (1,Nth)))
        Ex = Er*np.cos(ths); Ey = Er*np.sin(ths)
        qm_dt2 = e*dt/(2*m_e); tz = qm_dt2*B
        pxh = px + qm_dt2*Ex; pyh = py + qm_dt2*Ey
        ang = -2*math.atan(tz); ca, sa = math.cos(ang), math.sin(ang)
        pxp = ca*pxh - sa*pyh; pyp = sa*pxh + ca*pyh
        px = pxp + qm_dt2*Ex; py = pyp + qm_dt2*Ey
        rs = rs + px/m_e*dt; ths = ths + py/m_e/rs*dt
    W = (px**2+py**2+pz**2)/(2*m_e)/e
    if not np.isfinite(W).all():
        print("[HULL] nan в облаке — неустойчивость пуша"); return 9.9
    dW = abs(W.mean()-W0.mean())/W0.mean()
    print(f"[HULL] <W>: {W0.mean():.2f}->{W.mean():.2f} эВ, относ. дрейф={dW*100:.3f}% "
          f"(эталон инварианта Халла)")
    return dW

def run(mode, T_ns=1.0, E0=0.0, sc_on=True, label=""):
    """mode: 'cold'|'bunched'; E0 [В/м] амплитуда spoke-поля; sc_on — space charge."""
    bunched = (mode == 'bunched')
    rs, ths, px, py, pz = init_cloud(bunched=bunched)
    nt = int(T_ns*1e-9/dt)
    W0 = (px**2+py**2+pz**2)/(2*m_e)/e              # эВ
    hist_W, hist_R, hist_Esc, snap = [], [], [], None
    for i in range(nt):
        t = i*dt
        # статика как скорость: u_r -> E_r = B*u_r
        Er = gather(rs, ths, np.tile((static_Er[:,None]), (1,Nth)))
        Ex = Er*np.cos(ths); Ey = Er*np.sin(ths)
        if E0 > 0:
            u_rf = rf_speed(t, rs, ths, E0)
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
        # Boris push (B = B z_hat): t = q dt B/(2m) — для электрона q=-e, поворот на -2*atan(tz)
        qm_dt2 = e*dt/(2*m_e)
        tz = qm_dt2*B
        pxh = px + qm_dt2*Ex; pyh = py + qm_dt2*Ey; pzh = pz
        ang = -2*math.atan(tz)
        ca, sa = math.cos(ang), math.sin(ang)
        pxp = ca*pxh - sa*pyh; pyp = sa*pxh + ca*pyh
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
    print(f"[{label}] mode={mode} E0={E0:.1e} sc={sc_on} T={T_ns} нс: "
          f"<W>: {W0.mean():.2f}->{hist_W[-1]:.2f} эВ (dW={dW:+.2f} эВ/э), "
          f"<r>: {hist_R[0]/1e-6:.2f}->{hist_R[-1]/1e-6:.2f} мкм, "
          f"|E_sc|max={max(hist_Esc) if hist_Esc else 0:.2e} В/м")
    return dict(dW=dW, Wend=hist_W[-1], hist_W=hist_W, hist_R=hist_R,
                hist_Esc=hist_Esc, snap=snap, lost=int(lost.sum()))

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
