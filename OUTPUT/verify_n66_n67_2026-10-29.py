import math
m_e=9.1093837015e-31; e=1.602176634e-19
rk,ra,B,f=40e-6,80e-6,1.1,1e12
lnAB=math.log(ra/rk); w0=2*math.pi*f
VL=e*B*B*(ra-rk)**2/(2*m_e); VH=e*B*B*ra*(ra-rk)/(2*m_e)
# Синхронизм: Om_d = Va/(B*rk^2*lnAB) [рад/с] (в знаменателе rk^2 — проверено),
# n* = w0/Om_d = w0*B*rk^2*lnAB/Va. Безразмерная величина; её «числитель»
# K = n*·Va = w0*B*rk^2*lnAB имеет размерность В и равен 7665,1 В в BASE:
K=w0*B*rk**2*lnAB               # В (= n*·Va); Va(n)=K/n — точка целого корня
def Va_of_n(n): return K/n
print('V_L=%.1f V_H=%.1f  K=%.2f V' % (VL, VH, K))
for n in (33, 67):
    va = Va_of_n(n)
    print('n=%d: Va=%.1f V, pos=%.3f' % (n, va, (va-VL)/(VH-VL)))
nst = K/230.0
print('n*(Va=230)=%.4f  2n*=%.4f' % (nst, 2*nst))
def lam(x_n, Neff):
    x = x_n*math.pi/Neff
    return abs(math.sin(x)/x) if abs(x)>1e-15 else 1.0
print('Lambda(spoke vs Neff=33, N=66 standing)= %.4f' % lam(nst, 33))
print('Lambda(spoke vs Neff=67, N=67 traveling)= %.4f' % lam(nst, 67))
print('Lambda(spoke vs Neff=99, N=99=3n traveling)= %.4f' % lam(nst, 99))
