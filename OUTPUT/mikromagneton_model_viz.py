#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Визуализация модели микромагнетрона (масштабирование классического магнетрона).

Параметры — из корпуса Inbox/mikromagneton_concept_2026.md [в корпусе]:
  f = 300 ГГц (номинал П-2 ~1 ТГц), N = 128 каверн (π-мода, n = N/2),
  r_a = 80 мкм, r_k = 40 мкм, зазор d = 40 мкм, B = 1.0–1.5 Тл, V_a ~ 200–370 В.
Числа на панелях пересчитаны по формулам здесь [вычислено здесь].
Материал: автор проекта (czts-spec); расчёт выполнен в среде Qwen.
"""
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- константы и параметры корпуса ---
e = 1.602176634e-19; m = 9.1093837e-31; c = 2.99792458e8
f = 1.0e12; lam = c/f   # номинал корпуса: рабочая точка на 1 ТГц (таблица "Рабочая точка на 1 ТГц")
N = 128; n = N//2
ra, rk = 80e-6, 40e-6; d = ra - rk
B = 1.1
wc = e*B/m; fc = wc/(2*math.pi)
vphi = 2*math.pi*f*ra/n                    # фазовая скорость π-моды на аноде
V_Hull = m*wc**2*ra*d/(2*e)                # верхняя граница (циклотронная оценка): V = e*B^2*ra*d/(2m)
V_HB   = m*(wc*d)**2/(2*e)                 # нижняя граница: ларморовский критерий r_L(V) <= d

fig, ax = plt.subplots(2, 2, figsize=(12.5, 10))
fig.suptitle("Микромагнетрон THz — модельная визуализация (1 ТГц, N=128, B=1,1 Тл)", fontsize=13)

# --- Панель 1: сечение анодного блока ---
a = ax[0,0]; th = np.linspace(0, 2*np.pi, N+1)
a.add_patch(plt.Circle((0,0), ra*1e6, color="#b87333", alpha=.35))
a.add_patch(plt.Circle((0,0), rk*1e6, facecolor="#222", edgecolor="k"))
for t in th[:-1] + (th[1]-th[0])/2:        # каверны между лопатками
    x1,y1 = ra*1e6*math.cos(t), ra*1e6*math.sin(t)
    a.plot([x1, x1*1.6], [y1, y1*1.6], lw=3, color="#8c5a2b")
a.text(0,0,"катод\nr=40 мкм", ha="center", va="center", color="w", fontsize=9)
a.set_title(f"Сечение: 128 каверн, шаг ≈ {2*np.pi*ra*1e6/N:.1f} мкм; v_φ={vphi/c:.3f}c")
a.set_xlim(-140,140); a.set_ylim(-140,140); a.set_aspect("equal"); a.grid(alpha=.3)
a.set_xlabel("мкм"); a.set_ylabel("мкм")

# --- Панель 2: окно Хартри–Халла vs B ---
a = ax[0,1]
Bs = np.linspace(0.9, 1.5, 70); wcs = e*Bs/m
Vh = m*wcs*lam/(4*e); Vhb = m*(wcs*d)**2/(2*e)
a.fill_between(Bs, Vhb, Vh, color="gold", alpha=.5, label="расчётное окно")
# линия корпуса 200–370 В показана зелёной полосой
a.plot(Bs, Vh, "r-"); a.plot(Bs, Vhb, "b-")
a.axvline(B, ls=":", color="k"); a.axhspan(200, 370, color="green", alpha=.15)
a.text(1.12, (Vh.max()+Vhb.min())/2, f"B=1,1 Тл:\n[{Vhb[np.argmin(abs(Bs-B))]:.0f}; {Vh[np.argmin(abs(Bs-B))]:.0f}] В", fontsize=9)
a.set_xlabel("B, Тл"); a.set_ylabel("V_a, В")
a.set_title("Окно: [ларморовский порог mωc²d²/2e; циклотронная оценка mωc²·ra·d/2e]"); a.legend(fontsize=8); a.grid(alpha=.3)

# --- Панель 3: развёртка π-моды вдоль кольца ---
a = ax[1,0]
x = np.linspace(0, N, 1000)
a.plot(x, np.sin(np.pi*x), lw=2, color="crimson", label="π-мода E_z (знак чередуется)")
a.plot(x, 0.35*np.sin(np.pi*x+np.pi/2)*np.sin(np.pi*x), "--", color="gray", label="↦ пучок в захвате")
a.axhline(0, color="k", lw=.6)
a.set_xlabel("номер каверны (0…128)"); a.set_ylabel("E_z (отн.)")
a.set_title("Развёртка резонанса: период = 2 каверны; оборотов/период волны = "
            f"{fc/f:.0f}")
a.legend(fontsize=8); a.grid(alpha=.3)

# --- Панель 4: энергобаланс ячейки [в корпусе] ---
a = ax[1,1]
labels = ["P_питания\nV_a·I_a", "P_РФ выхода", "потери на аноде\n(скин 66 нм)", "КПД"]
# номинал корпуса: P_выхода 0,4–0,6 Вт при КПД 15–25 % => P_пит ~ 2,4–4 Вт; берём 0,5 Вт / 20 %
Pout, eta = 0.5, 0.20
Pin = Pout/eta; Ploss = Pin - Pout
vals = [Pin, Pout, Ploss]
bars = a.bar(range(3), vals, color=["steelblue","seagreen","indianred"])
a.set_xticks(range(4)); a.set_xticklabels(labels, fontsize=8)
a.set_ylabel("Вт (номинал 0,5 Вт)"); a.set_ylim(0,3.4)
for i,v in enumerate(vals): a.text(i, v+.06, f"{v:.1f}", ha="center")
a.text(3, 1.5, "КПД\n15–25 %\n(корпус)", ha="center", fontsize=10,
       bbox=dict(boxstyle="round", fc="gold", alpha=.6))
a.set_title("Баланс мощности ячейки (номинал П-1)"); a.grid(axis="y", alpha=.3)

fig.tight_layout(rect=[0,0,1,0.95])
out = "OUTPUT/mikromagneton_model_viz.png"
fig.savefig(out, dpi=130)
print("saved:", out)
print(f"v_phi={vphi/c:.4f}c | окно при B=1.1: [{V_HB:.0f}; {V_Hull:.0f}] В | f_c/f={fc/f:.1f}")
