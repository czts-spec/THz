#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Понижение степени системы домножением на 1/B (u = V_a/B).

Среда Qwen, THz/OUTPUT. Оси: [вычислено здесь].
Выводит решётку целых мод n=17..33 замкнутой формой (ноль итераций)
и сверяет строку n=33 с точкой корпуса (B~1,1 Тл, V_a~230 В, N=66).
"""
from math import pi, log, sqrt

m_e = 9.1093837e-31; e = 1.602176634e-19
f0 = 1e12; B0 = 1.1; Va0 = 230.0
rk = 40e-6; ra = 80e-6
lnAB = log(ra / rk)
w0 = 2 * pi * f0

C_n = w0 * rk * rk * lnAB                      # n = C_n / u, u = V_a/B
k_L = e * (ra - rk) ** 2 * B0 ** 2 / (2 * m_e) / B0 ** 2   # В/Тл^2 (V_L = k_L*B^2)
k_h = k_L * (370.0 / 195.0)                    # калибровка по окну корпуса

print(f"C_n = {C_n:.3f} В·мод/Тл | k_L = {k_L:.2f} | k_h = {k_h:.2f}")
u0 = Va0 / B0
print(f"база: u0={u0:.2f} В/Тл -> полоса B ∈ [{u0/k_h:.3f}; {u0/k_L:.3f}] Тл\n")

print(" n    u(В/Тл)   B_lo    B_mid   B_hi    V_a(В)   N=2n")
for n in range(17, 34):
    u = C_n / n
    Blo, Bhi = u / k_h, u / k_L
    Bmid = sqrt(Blo * Bhi)
    print(f"{n:2d}  {u:8.2f}  {Blo:6.3f}  {Bmid:6.3f}  {Bhi:6.3f}  {u*Bmid:7.1f}   {2*n}")

# самопроверка: n=33 воспроизводит точку корпуса
u33 = C_n / 33
B33 = sqrt((u33 / k_h) * (u33 / k_L))
V33 = u33 * B33
assert abs(V33 - 230.0) < 5.0 and abs(B33 - 1.1) < 0.05, "самопроверка не пройдена"
print(f"\n[OK] n=33: B={B33:.3f} Тл, V_a={V33:.1f} В, N=66 — совпадает с точкой корпуса.")
