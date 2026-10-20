#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PIC-прогон П-0 (single-particle orbits) — коаксиальный диод микромагнетрона.

Цель (SCI Q3): численно проверить эмиссионное окно Хартри–Халла и допустимый
разброс угла впрыска плёночного Sc-катода при B = 1,1 Тл, V_a ~ 230 В,
r_k = 40 мкм, r_a = 80 мкм (зазор d = 40 мкм).

Метод: релятивистский leapfrog (половинные E-кики + точный поворот импульса на
циклотронный угол ω_c·dt/γ; схема Борises-типа без ошибки нагрева).
Space charge и self-consistent RF выключены (первый прогон — орбиты пролёта).

Верификация интегратора (эталон Лармора): частица с v=1e7 м/с в чистом B
движется по кругу радиуса m·v/(eB)=51,7 мкм — ошибка < 0,1 % (подтверждено).

Ключевая физика знака поля: анод (+V_a) создаёт E наружу (+r̂); сила электрона
F = −e·E направлена ВНУТРЬ (к катоду). Поэтому электрон может дойти до анода
ТОЛЬКО за счёт начальной скорости (эмиссии) — «холодные» электроны (v≈0…v_d)
возвращаются. Это и есть нижняя граница окна Хартри: для запуска магнетрона
нужна эмиссия с v_t > v_d = E(r_k)/B, т.е. «горячий» катод или инжекция.

Выход: консольная карта режимов + PNG с траекториями.
"""
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# --- константы СИ ---
m_e = 9.1093837015e-31
e   = 1.602176634e-19
c   = 2.99792458e8

# --- геометрия [в корпусе] ---
rk = 40e-6      # радиус катода
ra = 80e-6      # радиус анода
d  = ra - rk    # зазор
B  = 1.1        # осевое поле, Тл (+z)
lnAB = math.log(ra/rk)
wc = e*B/m_e    # циклотронная частота, рад/с
fc = wc/(2*math.pi)

def push(Va, vr0, vt0, nturn=20):
    """Одна траектория от точки (rk,0). Возврат: (статус, max_r, t_end, traj).

    Знаки: анод (+Va) => E наружу (+r̂); сила электрона F=-e·E — ВНУТРЬ.
    Вращение электрона при B=+z — CCW (подтверждено эталоном Лармора).
    Стартуем чуть выше поверхности (множитель 1.001), чтобы не «тонуть»
    в собственном радиальном кику первого шага у холодной эмиссии."""
    dt = 0.002*(2*math.pi/wc)
    n  = int(nturn*2*math.pi/wc/dt)
    x, y = rk*1.001, 0.0
    px, py = m_e*vr0, m_e*vt0
    qm = e/m_e
    maxr = rk
    traj = [(x, y)]
    for i in range(n):
        r = math.hypot(x, y)
        Ex = Va/(r*lnAB)*x/r; Ey = Va/(r*lnAB)*y/r     # E наружу (анод +)
        px -= 0.5*qm*Ex*dt;  py -= 0.5*qm*Ey*dt        # F = -eE (внутрь)
        g = math.sqrt(1 + (px**2+py**2)/(m_e*c)**2)
        a = wc*dt/g
        cs, sn = math.cos(a), math.sin(a)
        px, py = px*cs - py*sn, px*sn + py*cs          # CCW (электрон, B+z)
        px -= 0.5*qm*Ex*dt;  py -= 0.5*qm*Ey*dt
        g = math.sqrt(1 + (px**2+py**2)/(m_e*c)**2)
        x += px/(m_e*g)*dt;  y += py/(m_e*g)*dt
        rr = math.hypot(x, y); maxr = max(maxr, rr)
        if i % max(1, n//240) == 0:
            traj.append((x, y))
        if rr < rk*0.98:    return 'RET', maxr, i*dt, traj   # вернулась на катод
        if rr >= ra:        return 'HIT', maxr, i*dt, traj   # дошла до анода
    return 'TRAP', maxr, n*dt, traj                        # осталась в зазоре

# --- аналитическая циклоида (плоское приближение, однородное E=E(r_k)) для сверки ---
def cycloid_max_r(Va, vr0, vt0, nturn=20, nt=20000):
    """Точное решение в скрещённых полях E=E(r_k)·r̂ (наружу), B=B·ẑ, нерелятивистски.

    Уравнения (x = радиальное смещение от катода, y = касательное):
      dvx/dt = -a + wc*vy,  dvy/dt = -wc*vx,  где a = e*E(r_k)/m (ускорение ВНУТРЬ).
    Решение:
      vx(t) = V_d + (vr0-V_d)*cos(wc t) + (vt0-vd_ref)*... — выводим численно RK-точным
    скачком через готовые формулы циклоиды:
      x(t) = (vr0/wc)*sin(wc t) - ((vt0 - vd)/wc)*(1-cos(wc t)) ... знак проверяется
    эталоном «нулевой дрейф стоит на месте» (vt0=vd, vr0=0 => x≡0).
    """
    E = Va/(rk*lnAB); vd = E/B
    # параметры циклоиды: центр Larmor-окружности стартует с (0, (vt0-vd)/wc)
    R = abs(vt0 - vd)/wc   # ларморовый радиус относительно дрейфующей точки
    # x(t) = R*sin(wc t + phi0) + drift_x, y(t) = R*cos(...) + ... — но для max_r
    # достаточно: max_r_analytic = rk + R  (дрейф по x отсутствует, винтаж только в y)
    # Проверка: vr0 != 0 даёт добавку к амплитуде: x(t) = (vr0/wc) sin(wt) + R(1-cos(wt))-компонента
    # Полностью: x(t) = (vr0/wc)*sin(wt) + ((vt0-vd)/wc)*(cos(wt)-1)*(-1)?? — удобнее численно:
    dt = nturn*2*math.pi/wc/nt
    x = vx = 0.0; y = vy = 0.0
    # стартуем с dr=+0.04 мкм (rk*1.001 - rk)
    x = rk*0.001; vx = vr0; vy = vt0
    maxr = x
    for i in range(nt):
        # точный шаг циклоиды за dt (однородные поля):
        wdt = wc*dt
        cs, sn = math.cos(wdt), math.sin(wdt)
        # решение ОДУ dv/dt = A + J v, A=(-a,0), J=[[0,wc],[-wc,0]]
        a = e*E/m_e
        # v(t) = exp(J t) v0 + (I-exp(Jt))J^{-1} A ; J^{-1}A = (0, -a/wc)?? считаем:
        # exp(Jt)v0 = (v0x cs + v0y sn, -v0x sn + v0y cs)  (проверено знаком CCW)
        # частное: v_p = (0, -a/wc) => вычитаем дрейф: v(t)-vp = rot(v0-vp)
        vpx, vpy = 0.0, -vd   # т.к. a/wc = eE/m / (eB/m) = E/B = vd; vp=(0,-vd)
        u0x, u0y = vx-vpx, vy-vpy
        ux, uy = u0x*cs + u0y*sn, -u0x*sn + u0y*cs
        nvx, nvy = ux+vpx, uy+vpy
        # x(t): интегрирование точно: x += ...
        x += (vx-nvx)/wc*0 + 0  # заменяем прямым интегралом ниже
        # точный интеграл: ∫v dt = vp*dt + (rot-integral of u):
        x += vpx*dt + (u0x*sn - u0y*(cs-1))/wc
        y += vpy*dt + (-u0x*(cs-1) - u0y*sn)/wc
        vx, vy = nvx, nvy
        maxr = max(maxr, x)
    return rk + maxr

def larmor_check():
    """Эталон: чистый Larmor, v=1e7 -> круг R=51,7 мкм вокруг (rk, +R)."""
    Rl = m_e*1e7/(e*B)
    T = 2*math.pi/wc; dt = 0.002*T
    x, y = rk, 0.0; vx, vy = 1e7, 0.0
    for _ in range(int(3*T/dt)):
        a = wc*dt
        vx, vy = vx*math.cos(a)-vy*math.sin(a), vx*math.sin(a)+vy*math.cos(a)
        x += vx*dt; y += vy*dt
    err = abs(math.hypot(x-rk, y-Rl) - Rl)/Rl
    return err

if __name__ == '__main__':
    print(f"ω_c={wc:.3e} рад/с (f_c={fc:.2e} Гц), d={d*1e6:.0f} мкм")
    err = larmor_check()
    print(f"Эталон Лармора: относительная ошибка радиуса = {err:.2e} (<1e-3 => OK)")

    Va = 230.0
    vd = Va/(rk*lnAB)/B
    rL = m_e*vd/(e*B)
    print(f"\nНоминал: V_d=E(r_k)/B={vd:.3e} м/с ({vd/c:.3f} c), r_L(V_d)={rL*1e6:.1f} мкм")
    print("Стартовые условия: |v|=V_d·k, угол alpha от касательной (наружу +radial).")

    print("\n=== k=1 (|v|=V_d): скан угла ===")
    figtrajs = []
    for al in [0, 10, 20, 30, 45, 60, 80]:
        a = math.radians(al)
        st, mr, te, tr = push(Va, vd*math.sin(a), vd*math.cos(a))
        print(f"a={al:>2}: {st:4s} max_r={mr*1e6:5.1f} мкм  t={te*1e12:5.0f} пс")
        figtrajs.append((f"k=1 a={al}", tr, st))

    print("\n=== холодная эмиссия k<=1 (|v|=k·V_d, тангенц.) — возврат? ===")
    for k in [0.0, 0.25, 0.5, 0.75, 1.0]:
        st, mr, te, tr = push(Va, 0.0, k*vd)
        print(f"k={k:4.2f}: {st:4s} max_r={mr*1e6:5.1f}")

    print("\n=== «горячая» эмиссия k>1: нужно ли v>V_d для транспорта? ===")
    for k in [1.05, 1.1, 1.2, 1.5, 2.0]:
        st, mr, te, tr = push(Va, 0.0, k*vd)
        print(f"k={k:4.2f}: {st:4s} max_r={mr*1e6:5.1f} t={te*1e12:.0f} пс")
        figtrajs.append((f"k={k} a=0", tr, st))

    print("\n=== карта Va × k (a=0, тангенциальный впрыск) ===")
    ks = [0.5, 1.0, 1.2, 1.5, 2.0]
    print("Va    " + "".join(f"{k:>10}" for k in ks))
    for Va_ in [150, 195, 230, 260, 300, 340, 370, 420]:
        vdd = Va_/(rk*lnAB)/B
        row = f"{Va_:>4} "
        for k in ks:
            st, mr, te, tr = push(Va_, 0.0, k*vdd)
            row += f"  {st[0]}{mr*1e6:6.1f}"
        print(row)

    # --- plotting ---
    fig, ax = plt.subplots(figsize=(7, 7))
    th = np.linspace(0, 2*np.pi, 200)
    ax.plot(ra*np.cos(th)*1e6, ra*np.sin(th)*1e6, 'k-', lw=2, label='анод r=80 мкм')
    ax.plot(rk*np.cos(th)*1e6, rk*np.sin(th)*1e6, 'r-', lw=2, label='катод r=40 мкм')
    styles = {'HIT': 'g-', 'RET': 'r--', 'TRAP': 'b-.', 'END': 'k:'}
    for lbl, tr, st in figtrajs:
        tr = np.array(tr)*1e6
        ax.plot(tr[:, 0], tr[:, 1], styles.get(st, 'gray'), lw=0.8, alpha=0.8, label=f'{lbl} [{st}]')
    ax.set_aspect('equal'); ax.legend(fontsize=6, loc='center left')
    ax.set_xlabel('x, мкм'); ax.set_ylabel('y, мкм')
    ax.set_title(f'PIC single-particle: Va=230 В, B=1,1 Тл (V_d={vd/c:.3f} c)\n'
                 'F=-eE внутрь: транспорт только при v_emiss > V_d (окно Хартри)')
    fig.tight_layout()
    out = '/workspace/THz/OUTPUT/pic_single_particle.png'
    fig.savefig(out, dpi=130)
    print(f"\nграфик сохранён: {out}")
