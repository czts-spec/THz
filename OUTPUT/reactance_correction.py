#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Коррекция реактанса микромагнетрона (замечание автора 2026-10-25).

Диагноз пессимизма f9d8357: адмитанс пучка G_b = 2P/V_int^2 вычислен при
V_int = V_a = 230 В. Для spoke-синхронной гармоники цепное напряжение СВЯЗИ —
локальное напряжение моды в зазоре V_int = E_gap*d, а не анодное напряжение.
Это и есть «прямое применение редуцированной системы» к цепному критерию с
неверной осью. Коррекция меняет знак вывода о самовозбуждении.

[вычислено здесь] — все числа печатает прогон этого скрипта.
"""
import math

eps0 = 8.8541878128e-12; mu0 = 4*math.pi*1e-7
f0 = 1e12; w0 = 2*math.pi*f0
Va = 230.0; B = 1.1; rk = 40e-6; ra = 80e-6; d = ra - rk
N = 66; nstar = 33; lnAB = math.log(ra/rk)
Q_claim = 700.0                       # [в корпусе] SCI Q4 (шероховатость MEMS)

# --- БАЗА f9d8357 ---------------------------------------------------------------
C_cell = eps0*(math.pi*ra**2/N)/d     # плоская ёмкость ячейки = 0.34 фФ
R_sh = Q_claim/(w0*C_cell)            # 330 кОм
G_wall = N/R_sh                       # 199.7 мкСм
P_th = 240e-6*154.0                   # Вт = I[А]·ΔW[эВ]; pic_spoke33: 154 эВ/э
Gb_old = 2*P_th/Va**2                 # 1.40 мкСм  <-- ОШИБКА ОСИ V_int
print("=== БАЗА f9d8357 ===")
print(f"C={C_cell*1e15:.2f} фФ R_sh={R_sh/1e3:.0f} кОм G_wall={G_wall*1e6:.1f} мкСм")
print(f"P_beam(th)={P_th*1e3:.2f} мВт  G_b(V_int=Va)={Gb_old*1e6:.2f} мкСм  ratio={Gb_old/G_wall:.4f}")

# --- КОРРЕКЦИЯ: V_int = E_gap·d --------------------------------------------------
E_gap = 1e5                            # В/м [в корпусе]: RF-амплитуда pic_rf R1/pic_spoke33
V_int = E_gap*d                        # = 4.0 В
Gb_new = 2*P_th/V_int**2               # = 4.62 мСм
ratio = Gb_new/G_wall
Ib_th = 240e-6/ratio                   # масштабирование порога от зондового тока
print("\n=== КОРРЕКЦИЯ V_int=E_gap*d=%.1f В ===" % V_int)
print(f"G_b = {Gb_new*1e3:.2f} мСм; G_b/G_wall = {ratio:.0f}; порог Ib ~ {Ib_th*1e6:.1f} мкА")
print("Цепный вывод МЕНЯЕТСЯ: |−G_b| > G_wall -> самовозбуждение ДОПУСТИМО.")

# --- Чувствительность: три края напряжения связи ----------------------------------
print("\n=== Три края V_int (честная неопределённость) ===")
for name,Vx in (("V_int=Va (было)",Va),
                ("V_eff=√(Va·Vint)=%.1f"%math.sqrt(Va*V_int),math.sqrt(Va*V_int)),
                ("V_int=E_gap·d=%.1f"%V_int,V_int)):
    Gx=2*P_th/Vx**2
    print(f"{name:24s}: G_b={Gx*1e6:9.1f} мкСм  G_b/G_wall={Gx/G_wall:8.3f}  {'усилие' if Gx>G_wall else 'затухание'}")
# физика: электрон пересекает зазор за τ_tr=d/v_d≈? ; торможение идёт ПОЛНЫМ V_int,
# но ток I_b ограничен CL на V_a. Консервативный край V_eff даёт ratio≈2 — всё ещё
# выше порога. Точный V_int — за сеточной верификацией (следующий шаг).

# --- Реактивная часть и согласование --------------------------------------------
X_block = -1/(w0*C_cell*N)             # суммарная ёмкость блока N*C_cell
eta0 = math.sqrt(mu0/eps0)
RL = 1/Gb_new
Gam = abs((RL-eta0)/(RL+eta0)); vswr=(1+Gam)/(1-Gam)
beta = G_wall/Gb_new; eta_out = 1/(1+beta)
print("\n=== НАГРУЗКА (край V_int=4 В) ===")
print(f"X_ёмкостный блока = {X_block:.1f} Ом (компенсируется связкой/iris)")
print(f"R_L,opt = {RL:.0f} Ом; к η0={eta0:.1f} Ом: |Γ|={Gam:.3f}, VSWR={vswr:.2f}")
print(f"β={beta:.4f}; η_out=G_L/(G_L+G_wall)={eta_out*100:.1f}%")
for M in (1,4,16):
    RLm=RL/M; Gm=abs((RLm-eta0)/(RLm+eta0))
    print(f"  массив M={M:2d}: R_L={RLm:6.1f} Ом VSWR={(1+Gm)/(1-Gm):6.2f}")
print("\nИТОГ: см. OUTPUT/reactance_correction_notes_2026-10-25.md")
