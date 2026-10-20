#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PIC-прогон П-0: коаксиальный диод микромагнетрона в скрещённых полях E_r x B_z.
Модель: тестовая частица (space-charge выключен, self-consistent RF выключен) —
«single-particle orbit» прогон для проверки эмиссионного окна Хартри–Халла и
угла впрыска плёночного Sc-катода (SCI Q3).

Геометрия [в корпусе]: r_k = 40 мкм (катод), r_a = 80 мкм (анод), V_a = 230 В,
B = 1,1 Тл (осевая, +z). Поле E_r = V_a / (r * ln(r_a/r_k)) (коаксиал).
Эмиссия: с поверхности катода, скорость tan(alpha)*v_drift, alpha в [-60..+60]°.

Критерии окна заимствованы из справочной литературы по магнетронам
(Hartree boundary / "electric field cutoff" критерий для коаксиального диода):
  нижняя граница  V_H ~ m/e * omega_c^2 * r_k * d * f_geom,   f_geom = ln(a/b)/(a/b - 1)
  верхняя граница V_Hull ~ m * omega_c^2 * d * (r_a - d/2) / e  (циклотронный радиус = d)
Здесь проверяем численно: при Va < V_low частица не доходит до анода (возврат),
при Va > V_high Larmor-радиус > зазора -> удар в анод на первом пол-обороте
(запуск есть, но нет «десятков оборотов»).

Интегрирование: релятивистски-корректный leapfrog (Борис-подход через u=p/m),
dt выбран так, что omega_c*dt < 0.05 рад.
"""
import math
import numpy as np

# --- константы СИ ---
m_e = 9.1093837015e-31
e   = 1.602176634e-19
c   = 2.99792458e8

# --- геометрия [в корпусе] ---
rk = 40e-6      # катод
ra = 80e-6      # анод
d  = ra - rk    # зазор 40 мкм
B  = 1.1        # Тл
Ncav = 128
f    = 1.0e12   # рабочая точка 1 ТГц (таблица корпуса)
wc   = e*B/m_e                    # циклотронная частота
fc   = wc/(2*math.pi)             # ~30.8 ГГц
lnAB = math.log(ra/rk)

def Er(r, Va):
    return Va/(r*lnAB)

def run_particle(Va, alpha_deg, tmax_ps=None, dt=None):
    """Одна траектория. Возвращает dict с итогом."""
    v0t = None
    # дрейф E/B в коаксиале у катода: v_d = Er(rk)/B... точнее v_E = E/B
    vd = Er(rk, Va)/B
    v0 = vd/math.cos(math.radians(alpha_deg))   # тангенциальный впрыск: |v| cos a = v_d? 
    # физически: пучок испускается с малой скоростью; задаём долю от v_th безразмерно:
    # упрощаем: v0 = vd / cos(alpha) => v_phi = vd*tg? Нет. Стандартная постановка
    # PIC-теста Hartree: v_init tangential = vd*tan(alpha_scale). Берём простую схему:
    v_t = vd*math.tan(math.radians(alpha_deg)) if abs(alpha_deg)<89 else vd*100
    v_r = vd
    if dt is None:
        dt = 0.02/wc          # шаг по времени
    if tmax_ps is None:
        T = 2*math.pi/wc
        tmax_ps = 40*T        # 40 циклотронных периодов (~1300 волновых периодов волны?)
    nstep = int(tmax_ps/dt)
    r = rk
    th = 0.0
    vr, vt = v_r, v_t
    ux, uy = 0.0, 0.0         # нормированные импульсы p/(m c): старт нерелятив.
    # перпендикулярные компоненты в декартовых (плоскость xy), B along z
    vx, vy = vr, vt           # в момент старта: radial along +x
    gamma = 1.0/math.sqrt(max(1e-12, 1-(vx**2+vy**2)/c**2))
    px, py = m_e*gamma*vx, m_e*gamma*vy
    x, y = rk, 0.0
    hit_anode = False; returned = False; min_r = rk; max_r = rk
    turns = 0.0
    traj = []
    save_every = max(1, nstep//4000)
    for i in range(nstep):
        rr = math.hypot(x, y)
        if rr <= 0: break
        Ex = Er(rr, Va)*x/rr; Ey = Er(rr, Va)*y/rr
        qm_e = e/m_e
        # half-kick electric
        px += 0.5*qm_e*Ex*dt; py += 0.5*qm_e*Ey*dt
        # magnetic rotation (Boris): t = q dt B /(2 m gamma)
        g = math.sqrt(1+(px**2+py**2)/(m_e*c)**2)
        tvec = 0.5*qm_e*B*dt/g
        sx, sy = 2*tvec*(-py), 2*tvec*(px)   # s = 2t/(1+t^2), учтём знаменатель
        den = 1+tvec**2
        sx /= den; sy /= den
        pxp, pyp = px+sx*0 + (-sy)+px*0, py   # аккуратнее ниже
        # простая схема: p' = p + p x t ... делаем явно
        px_, py_ = px, py
        # p_minus_half уже включает E; теперь ротация: p' = pmh + (pmh + p x s)? 
        # Используем стандарт: tp = q dt B/(2m g); s = 2tp/(1+tp^2); 
        # p_plus = p_minus + cross(p_minus + cross(p_minus, tp_vec), s_vec)... 
        # Для B along z это просто поворот: 
        ang = 2*math.atan(tvec)
        cs, sn = math.cos(ang), math.sin(ang)
        px, py = px_*cs - py_*sn, px_*sn + py_*cs
        # second half E-kick
        px += 0.5*qm_e*Ex*dt; py += 0.5*qm_e*Ey*dt
        g = math.sqrt(1+(px**2+py**2)/(m_e*c)**2)
        vx, vy = px/(m_e*g), py/(m_e*g)
        x += vx*dt; y += vy*dt
        rr = math.hypot(x, y)
        if rr < rk: returned = True; break
        if rr >= ra: hit_anode = True; max_r = rr; break
        max_r = max(max_r, rr); min_r = min(min_r, rr)
        turns += abs((x*vy - y*vx))/rr**2 * dt / (2*math.pi)  # dtheta/2pi signed
        if i % save_every == 0:
            traj.append((x, y))
    return dict(hit=hit_anode, ret=returned, turns=turns, max_r=max_r,
                traj=np.array(traj) if traj else np.zeros((0,2)),
                nstep=i, t_end=(i+1)*dt)

def classify(Va, alpha):
    res = run_particle(Va, alpha)
    if res['ret']: return 'NO-TRANSPORT (вернулась на катод)'
    if res['hit'] and res['turns'] < 3: return 'HULL-OVER (удар в анод сразу)'
    if res['hit']: return 'SUSTAINED (дошла за N оборотов)'
    return 'TRAPPED-BETWEEN'

if __name__ == '__main__':
    print(f"omega_c={wc:.3e} rad/s, f_c={fc:.3e} Hz, d={d*1e6:.0f} мкм")
    # аналитические границы (справочные формулы коаксиального магнетрона)
    fgeom = lnAB/(ra/rk - 1)
    Vlow  = m_e/m_e*0  # placeholder
    Vlow  = (m_e/e)*wc**2 * rk * d * fgeom
    Vhigh = (m_e/e)*wc**2 * d * (ra - d/2)
    print(f"V_low (Hartree, коакс.)  = {Vlow:.1f} В")
    print(f"V_high (Larmor=d)        = {Vhigh:.1f} В")
    print()
    grid = []
    for Va in [100, 150, 180, 200, 230, 260, 300, 340, 380, 420]:
        row = []
        for alpha in [0, 30, 60]:
            r = run_particle(Va, alpha)
            tag = ('R' if r['ret'] else ('A-fast' if (r['hit'] and r['turns']<3) else
                   ('A-sust' if r['hit'] else 'T')))
            row.append(f"{tag}(n={r['turns']:.0f})")
        grid.append((Va,row)); print(f"Va={Va:>3} В: " + "  ".join(row))
