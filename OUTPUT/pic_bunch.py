#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PIC П-1b: формирование spoke ( bunching ) в скрещённых полях + порог самовозбуждения.
Закрывает разрыв П-1/R2: равномерное облако моду не возбуждало (нет фазировки).
Вопрос: возникает ли spoke из МАЛОГО начального сгустка и при каком токе стартует мода?

Метод тот же lumped-cavity [внешний источник: Hutchinson ch.15], поля — как в pic_rf.py
[в корпусе: ENV концепта, mikromagneton_concept_2026.md]. Новое здесь [вычислено здесь]:
  - стартовое распределение: доля frac облака в окне dpsi вокруг ψ* (остальное равномерно);
  - критерий старта: E0(t) растёт дольше 2·T_dec; порог по Ib — удвоение поиска;
  - метрика bunching: R = |<exp(i·neff·θ)>| (порядковый параметр spoke).
Запуск: python3 OUTPUT/pic_bunch.py (нужны numpy, matplotlib).
"""
import math, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt

m_e=9.1093837015e-31; e=1.602176634e-19; c=2.99792458e8; eps0=8.8541878128e-12
rk,ra=40e-6,80e-6; d=ra-rk; lnAB=math.log(ra/rk)
B=1.1; f0=1.0e12; w0=2*math.pi*f0; N=128; neff=N//2
Q=700.; gamma_f=0.5; Vdc0=230.
wc=e*B/m_e; Tdec=Q/w0; Cst=eps0*d; k_c=1.0/(w0*Cst*gamma_f*d)
PSI_STAR=5.8   # тормозящая фаза из V2 (pic_rf.py) [вычислено здесь]

def step(x,y,px,py,E0,t,Vdc):
    """Boris-шаг, идентичный верифицированному pic_rf.py [в корпусе OUTPUT]:
    половина толчка E — точный поворот в B — половина толчка E+RF; RF-мощность."""
    r=np.hypot(x,y); r=np.where(r<1e-12,1e-12,r); th=np.arctan2(y,x)
    Er=Vdc/(r*lnAB); Ep=E0*(rk/r)*np.cos(neff*th-w0*t)
    rx,ry=x/r,y/r
    Fx=-e*(Er*rx-Ep*ry); Fy=-e*(Er*ry+Ep*rx)     # сила от E (на электрон: −e·E)
    px=px+Fx*DT/2; py=py+Fy*DT/2
    g=np.sqrt(1+(px*px+py*py)/(m_e*c)**2)
    a=wc*DT/g; cs,sn=np.cos(a),np.sin(a)          # точный ротационный B-поворот
    px,py=px*cs+py*sn,-px*sn+py*cs
    vt=(-px*np.sin(th)+py*np.cos(th))/(g*m_e)
    Pu=float(np.sum(-e*Ep*vt))                    # Вт при q_mac=e (>0 = отдача волне)
    px=px+Fx*DT/2; py=py+Fy*DT/2
    x=x+px/(g*m_e)*DT; y=y+py/(g*m_e)*DT
    return x,y,px,py,Pu

def emit(M,rng,Vdc,bunched_frac=0.0,dpsi=0.6):
    """Hull-старт как в pic_rf.emit_hull: зародыш spoke занимает целые каверны
    (cell·2π/neff + окно ψ∈[ψ*,ψ*+dpsi]); фон — равномерный по 2π/neff; тепловая
    скорость Sc 0,1 эВ добавляется радиально [вычислено здесь]."""
    rs=rk*1.001
    nb=int(round(bunched_frac*M))
    cell=rng.integers(0,neff,M)*2*math.pi/neff
    psi=rng.uniform(0,2*math.pi/neff,M)            # фон: равномерный внутри одного периода π-моды
    if nb: psi[:nb]=PSI_STAR+rng.uniform(0,dpsi,nb)
    th=(cell+psi)/neff
    x=rs*np.cos(th); y=rs*np.sin(th)
    vh=Vdc/(rs*lnAB)/B
    px=np.zeros(M); py=vh*m_e*np.ones(M)
    vr=rng.normal(0,2e4,M)
    px+=vr*m_e*np.cos(th); py+=vr*m_e*np.sin(th)
    return x,y,px,py

def run(Ib_uA,frac,M=256,T_ns=6.0,nt=30000,seed=1,E0_0=1e2,Vdc=Vdc0):
    global DT
    DT=T_ns*1e-9/nt
    rng=np.random.default_rng(seed)
    x,y,px,py=emit(M,rng,Vdc,bunched_frac=frac)
    rs=rk*1.001; vh=Vdc/(rs*lnAB)/B
    tau_cloud=(2*math.pi/neff)/vh
    qm=Ib_uA*1e-6*tau_cloud/M
    E0=E0_0; tl=[];El=[];Rl=[];Wl=[]
    for i in range(nt):
        x,y,px,py,Pu=step(x,y,px,py,E0,i*DT,Vdc)
        th=np.arctan2(y,x); rr=np.hypot(x,y)
        R=float(np.abs(np.mean(np.exp(1j*neff*th))))     # bunching-параметр
        W=((np.sqrt(1+(px*px+py*py)/(m_e*c)**2)-1)*m_e*c*c/e).mean()
        Pr=Pu*qm/e                                        # Вт, вес макрочастицы
        E0=max(E0+(k_c*Pr-E0/Tdec)*DT,0.)
        bad=(rr<rk*0.9)|(rr>ra*1.05); nb=int(bad.sum())
        if nb:
            xn,yn,pn,qn=emit(nb,np.random.default_rng(seed+i),Vdc,bunched_frac=0.0)
            x[bad]=xn;y[bad]=yn;px[bad]=pn;py[bad]=qn
        if i%200==0:
            tl.append(i*DT*1e12); El.append(E0); Rl.append(R); Wl.append(W)
    return map(np.array,(tl,El,Rl,Wl))

print("=== A: рост spoke из малого зародыша (frac=0.15, Ib=20 мкА) ===")
t,E,R,W=run(20,0.15,seed=5)
print(f"E0: {E[0]:.2e} -> {E.max():.2e} В/м; R: {R[0]:.2f} -> max {R.max():.2f}; "
      f"<W>: {W[0]:.0f}->{W[-1]:.0f} эВ; T_dec={Tdec*1e12:.1f} пс")
start = E.max()>10*max(E[0],1e-9) and t[E.argmax()]>2*Tdec*1e12
print(f"старт моды за 6 нс: {'ДА' if start else 'НЕТ (затухает)'}")

print("\n=== B: порог тока (frac=0.15, поиск удвоением) ===")
thr=None
for Ib in [20,100,500,2000,10000]:
    t,E,R,W=run(Ib,0.15,seed=7)
    grow = E.max()>10*E[0] and t[E.argmax()]>2*Tdec*1e12
    print(f"Ib={Ib:6d} мкА: E0max={E.max():.3e}, Rmax={R.max():.2f}, рост={'ДА' if grow else 'нет'}")
    if grow and thr is None: thr=Ib
print(f"=> порог самовозбуждения с зародышем ~{thr} мкА (между {thr//5 if thr else '>10000'} и {thr})" if thr else "=> порога нет до 10 мА")

print("\n=== C: усиление зародыша (Ib=20 мкА, frac=0.05/0.15/0.5) ===")
fig,ax=plt.subplots(1,2,figsize=(12,4.5))
for frac,col in [(0.05,'C0'),(0.15,'C1'),(0.5,'C2')]:
    t,E,R,W=run(20,frac,seed=9)
    ax[0].semilogy(t,np.maximum(E,1e-1),color=col,label=f"frac={frac}")
    ax[1].plot(t,R,color=col,label=f"frac={frac}")
ax[0].set_xlabel("t, пс"); ax[0].set_ylabel("E0, В/м"); ax[0].set_title("Усиление зародыша: поле")
ax[1].set_xlabel("t, пс"); ax[1].set_ylabel("R=|<e^{iNθ}>|"); ax[1].set_title("Порядок spoke")
for a in ax: a.legend(); a.grid(alpha=.3)
plt.tight_layout(); plt.savefig("OUTPUT/pic_bunch.png",dpi=120)
print("PNG: OUTPUT/pic_bunch.png")
