#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PIC второго уровня (П-0b): пучок в π-моде ВЧ-поля микромагнетрона. SCI Q6.

Метод — lumped-cavity PIC (аналог Hutchinson, Principles of MM Wave Electronics,
ch.15 [внешний источник]): макрочастицы в полях статика+RF; амплитуда моды E0
растёт от мощности пучка и затухает с T_dec = Q/ω0:
    dE0/dt = k_c·P_beam − E0/T_dec,   k_c = 1/(ω0·C_st·γ_f·d), C_st = ε0·d
(ёмкостная модель горла каверны [вычислено здесь]).

Поля [вычислено здесь]:
  Статика: E_r = Vdc/(r·ln(ra/rk));  RF: E_φ = E0·(rk/r)·cos(n_eff·θ − ω0·t),
  n_eff = N/2 = 64 — π-мода чередующихся каверн (множитель rk/r — концентрация
  поля у катода; точная мода — за 3D-решателем, SCI Q4).

СТАРТОВАЯ КОНФИГУРАЦИЯ — строго по Халлу [вычислено здесь]:
  равновесие у катода: v_θ(r_s) = E(r_s)/(B·r_s) = Vdc/(B·r_s·lnAB), v_r=0,
  r_s = r_k·1.001. Баланс сил eE(r_s) = e·v_θ·B выполнен (проверка V1 ниже);
  возмущение v_r=0 при r_s>r_k даёт медленный дрейф Эпicycloid наружу
  (циклотронный резонанс магнетрона), а НЕ мгновенное падение на катод.
  (Ранее использованный «чисто дрейфовый» старт v_θ=v_d(r_k) при r=r_k·1.001
  не удовлетворял балансу: eE > e·v·B → часть падала на катод и набирала
  энергию из СТАТИКИ — артефакт, устранён в этой версии.)

КЛЮЧЕВЫЕ РЕЗУЛЬТАТЫ ПРОГОНА:
  R1. Одночастичный съём энергии из волны: максимум ~86 эВ/э за 2 нс при
      E0=1e5 В/м, тормозящая фаза ψ*≈5,8 рад (~330°) [вычислено здесь];
      диапазон 40–90 эВ/э по всей тормозящей половине фазы — съём УСТОЙЧИВ.
  R2. Равномерное холодное облако НЕ самовозбуждает моду (<P>≈−3 мкВт…−0,3 мВт
      при Ib до 2 мА; знак «−» = поле качает пучок, т.е. нет bunching-фазировки):
      порог запуска требует формирования spoke (как в реальных магнетронах —
      через шум/ёмкостную связь). За 3 нс рост E0 не наступил ни при Ib≤2 мА.
  R3. Заранее спокенное облако (ψ∈[ψ*, ψ*+0,6]) при фикс. E0=1e5 В/м даёт
      устойчивый тормозящий съём <W>: 160→85 эВ, η = ΔW/Vdc ≈ 0,33 и
      P(Ib=20 мкА) ≈ 1,5 мВт — верхняя оценка КПД вынужденного режима.
      Пороговая проверка потерь стенки: P_st(E0=1e5)×128 каверн ≈ 0,8 мВт —
      ТОГО ЖЕ порядка, что съём (1,5 мВт при Ib=20 мкА) ⇒ баланс достижим
      при E0≈1,4·10^5 В/м; заявка «стартовый ток 20 мкА» на границе
      выполнимости — нужен либо Q>700, либо injection locking (SCI Q6).

Верификации:
  V1 Hull-старт без RF: баланс сил, |ΔW| мало, r растёт медленно (epicycloid);
  V2 скан ψ_st: максимум съёма 86 эВ/э при ψ*≈5,8 рад;
  V3 детализация dt: Δ<5% между nt=1e5 и 2e5 в окрестности оптимума.

Выход: консоль + OUTPUT/pic_rf.png.
"""
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

m_e=9.1093837015e-31; e=1.602176634e-19; c=2.99792458e8; eps0=8.8541878128e-12
rk,ra=40e-6,80e-6; d=ra-rk; lnAB=math.log(ra/rk)
B=1.1; f0=1.0e12; w0=2*math.pi*f0; N=128; neff=N//2
Q=700.; gamma_f=0.5; Vdc0=230.
wc=e*B/m_e
Tdec=Q/w0; Cst=eps0*d; k_c=1.0/(w0*Cst*gamma_f*d)

def vd_of(Vdc): return Vdc/(rk*lnAB)/B          # E×B дрейф у катода
def vHull(r,Vdc): return Vdc/(r*lnAB)/B          # Hull: v_th=E(r)/(Br) => круг r=const

dt_step=0.002*2*math.pi/wc

def step(x,y,px,py,E0,t,Vdc):
    r=np.hypot(x,y); r=np.where(r<1e-12,1e-12,r); th=np.arctan2(y,x)
    Er=Vdc/(r*lnAB); Ep=E0*(rk/r)*np.cos(neff*th-w0*t)
    rx,ry=x/r,y/r
    Fx=-e*(Er*rx-Ep*ry); Fy=-e*(Er*ry+Ep*rx)
    px=px+Fx*dt_step/2; py=py+Fy*dt_step/2
    g=np.sqrt(1+(px*px+py*py)/(m_e*c)**2)
    a=wc*dt_step/g; cs,sn=np.cos(a),np.sin(a)
    px,py=px*cs+py*sn,-px*sn+py*cs
    vt=(-px*np.sin(th)+py*np.cos(th))/(g*m_e)
    Pu=float(np.sum(-e*Ep*vt))                   # Вт при q_mac=e
    px=px+Fx*dt_step/2; py=py+Fy*dt_step/2
    x=x+px/(g*m_e)*dt_step; y=y+py/(g*m_e)*dt_step
    return x,y,px,py,Pu

def emit_hull(M,rng,Vdc,rscale=None,bunched=False,psistar=None,dpsi=0.6):
    """Старт с Hull-орбит: r=r_s (равновесие), v_t=vHull(r_s), v_r=0."""
    if rscale is None: rs=rk*1.001
    else: rs=rscale
    if bunched:
        cell=rng.integers(0,neff,M)*2*math.pi/neff
        psi=psistar+rng.uniform(0,dpsi,M)
        th=(cell+psi)/neff
    else:
        th=rng.uniform(0,2*math.pi/neff,M)
    x=rs*np.cos(th); y=rs*np.sin(th)
    vh=Vdc/(rs*lnAB)/B
    px=np.zeros(M); py=vh*m_e*np.ones(M)
    vr=rng.normal(0,2e4,M)                        # тепловая Sc ~0,1 эВ [внешний]
    px+=vr*m_e*np.cos(th); py+=vr*m_e*np.sin(th)
    return x,y,px,py

# ---------- V1: Hull-орбита без RF ----------
print("=== V1: Hull-старт, E0=0 ===")
dt_step=0.002*2*math.pi/wc
rng=np.random.default_rng(0); x,y,px,py=emit_hull(4,rng,Vdc0)
W0=(np.sqrt(1+(px*px+py*py)/(m_e*c)**2)-1)*m_e*c*c/e
for i in range(3000): x,y,px,py,_=step(x,y,px,py,0.0,i*dt_step,Vdc0)
W1=(np.sqrt(1+(px*px+py*py)/(m_e*c)**2)-1)*m_e*c*c/e
rr=np.hypot(x,y)
print(f"max|ΔW|={float(np.max(np.abs(W1-W0))):.4f} эВ, r={np.round(rr*1e6,2)} мкм "
      f"(старт {rk*1.001*1e6:.1f}) за {3000*dt_step*1e12:.1f} пс")

# ---------- V2/V3: одночастичный съём ----------
def one_psi(psi,E0=1e5,T_ns=2.0,nt=100000,Vdc=Vdc0):
    global dt_step
    dt_step=T_ns*1e-9/nt
    rng=np.random.default_rng(7)
    x,y,px,py=emit_hull(1,rng,Vdc,rscale=rk*1.001,bunched=True,psistar=psi,dpsi=1e-9)
    Wstart=0.5*m_e*(Vdc/(rk*1.001*lnAB)/B)**2/e
    for i in range(nt): x,y,px,py,_=step(x,y,px,py,E0,i*dt_step,Vdc)
    Wend=(math.sqrt(1+(px[0]**2+py[0]**2)/(m_e*c)**2)-1)*m_e*c*c/e
    return Wstart-Wend
print("\n=== V2/V3: скан ψ_st (Hull-старт, E0=1e5 В/м; nt=1e5 | 2e5) ===")
best=None
for ps in np.arange(0,2*math.pi,math.pi/6):
    s1=one_psi(ps,nt=100000); s2=one_psi(ps,nt=200000)
    print(f"ψ={ps:5.2f}: {s1:+7.1f} | {s2:+7.1f} эВ")
    if best is None or s1>best[1]: best=(ps,s1)
print(f"=> ψ*≈{best[0]:.2f} рад, {best[1]:.1f} эВ/э")

# ---------- A: равномерное облако, самосогласование ----------
def run_beam(Ib_uA,M=256,T_ns=3.0,nt=15000,seed=1,E0_0=1e3,Vdc=Vdc0,bunched=False):
    global dt_step
    dt_step=T_ns*1e-9/nt
    rng=np.random.default_rng(seed)
    x,y,px,py=emit_hull(M,rng,Vdc,bunched=bunched,psistar=best[0])
    rs=rk*1.001; vh=Vdc/(rs*lnAB)/B
    tau_cloud=(2*math.pi/neff)/vh
    qm=Ib_uA*1e-6*tau_cloud/M
    E0=E0_0; tl=[];El=[];Wl=[];Pl=[]
    for i in range(nt):
        x,y,px,py,Pu=step(x,y,px,py,E0,i*dt_step,Vdc)
        Pr=Pu*qm/e
        E0=max(E0+(k_c*Pr-E0/Tdec)*dt_step,0.)
        rr=np.hypot(x,y); bad=(rr<rk*0.9)|(rr>ra*1.05); nb=int(bad.sum())
        if nb:
            xn,yn,pn,qn=emit_hull(nb,rng,Vdc,bunched=bunched,psistar=best[0])
            x[bad]=xn;y[bad]=yn;px[bad]=pn;py[bad]=qn
        if i%100==0:
            W=((np.sqrt(1+(px*px+py*py)/(m_e*c)**2)-1)*m_e*c*c/e).mean()
            tl.append(i*dt_step*1e12);El.append(E0);Wl.append(W);Pl.append(Pr)
    return map(np.array,(tl,El,Wl,Pl))

print("\n=== A: равномерное облако (реальный запуск), самосогласование ===")
for Ib in [20,100,500,2000]:
    t,E,W,P=run_beam(Ib,M=256,T_ns=3.0,nt=15000,seed=5)
    print(f"Ib={Ib:5d} мкА: E0max={E.max():.3e} В/м, <P>={P.mean():+.3e} Вт")

print("\n=== B: заранее spoke (ψ*..ψ*+0,6), самосогласование ===")
res={}
for Ib in [20,100,500,2000]:
    t,E,W,P=run_beam(Ib,M=256,T_ns=3.0,nt=15000,seed=5,bunched=True)
    res[Ib]=(t,E,W,P)
    Psat=w0*0.5*Cst*(gamma_f*d*E.max())**2/Q
    print(f"Ib={Ib:5d} мкА: E0max={E.max():.3e} В/м, P_end={P[-1]:.3e} Вт, "
          f"P_st@E0max={Psat:.3e} Вт/кав (×128={Psat*128:.3e} Вт)")

# ---------- C: КПД вынужденного режима (поле фикс., полный захват) ----------
print("\n=== C: оценка КПД при E0_fix=1e5 В/м, spoke-облако ===")
dt_step=3e-9/15000
rng=np.random.default_rng(11); M=512
x,y,px,py=emit_hull(M,rng,Vdc0,bunched=True,psistar=best[0])
vh=Vdc0/(rk*1.001*lnAB)/B
qm=20e-6*(3e-9/15000)/M   # вес за шаг прогона (консервативно: τ=nt*dt)
Wstart=(np.sqrt(1+(px*px+py*py)/(m_e*c)**2)-1)*m_e*c*c/e
Wlog=[]
for i in range(15000):
    x,y,px,py,_=step(x,y,px,py,1e5,i*dt_step,Vdc0)
    rr=np.hypot(x,y); bad=(rr<rk*0.9)|(rr>ra*1.05); nb=int(bad.sum())
    if nb:
        xn,yn,pn,qn=emit_hull(nb,rng,Vdc0,bunched=True,psistar=best[0])
        x[bad]=xn;y[bad]=yn;px[bad]=pn;py[bad]=qn
    if i%3000==0:
        W=((np.sqrt(1+(px*px+py*py)/(m_e*c)**2)-1)*m_e*c*c/e)
        Wlog.append((i*dt_step*1e12, float(np.mean(W))))
t0=Wlog[0][1]; tend=Wlog[-1][1]
dWe=t0-tend                       # эВ/э отдано волне
eta=dWe*e/(e*Vdc0)                # КПД = снятая энергия / e·Vdc (на электрон)
print(f"<W>: {t0:.1f} -> {tend:.1f} эВ; съём {dWe:.1f} эВ/э; "
      f"η = dW/Vdc = {eta:.3f}; P(Ib=20мкА) = {20e-6*Vdc0*eta:.3e} Вт")

# ---------- график ----------
fig,ax=plt.subplots(1,3,figsize=(15,4.5))
psi=np.arange(0,2*math.pi,math.pi/12)
dt_step=2e-9/60000
ax[0].plot(psi,[one_psi(p,nt=60000) for p in psi],'o-')
ax[0].set_xlabel("ψ_st"); ax[0].set_ylabel("съём, эВ/э")
ax[0].set_title("V2: одночастичный съём (Hull-старт)")
for Ib,col in [(20,'C0'),(2000,'C1')]:
    t,E,W,P=res.get(Ib,run_beam(Ib,M=128,T_ns=3.0,nt=8000,seed=5,bunched=True))
    ax[1].semilogy(t,np.maximum(E,1e-1),color=col,label=f"Ib={Ib} мкА")
ax[1].set_xlabel("t, пс"); ax[1].set_ylabel("E0, В/м")
ax[1].set_title("B: spoke-облако, самосогласование"); ax[1].legend()
tt=[a for a,_ in Wlog]; ww=[b for _,b in Wlog]
ax[2].plot(tt,ww,'s-'); ax[2].set_xlabel("t, пс"); ax[2].set_ylabel("<W>, эВ")
ax[2].set_title(f"C: вынужденный съём, η≈{eta:.2f}")
plt.tight_layout(); plt.savefig("OUTPUT/pic_rf.png",dpi=120)
print("\nPNG: OUTPUT/pic_rf.png")
