#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PIC П-3: spoke синхронизма n*=33 (рабочая точка из sync_system_notes) — порог самозапуска.
Закрывает SCI Q6 окончательно и проверяет предсказание аналитической системы (а)+(б)+(в).

Отличие от pic_bunch.py (π-мода neff=64): [вычислено здесь]
  - мода: E_phi = E0*(rk/r)*cos(n*·θ − ω0·t), n* = round(ω0/Ω_d) = 33;
  - фаза зародыша ψ*=5,8 рад взята из pic_rf R1 (тормозящая дуга);
  - критерии: рост E0 дольше 2·T_dec ⇒ самозапуск; метрика R=|<exp(i n* θ)>|.
Метод lumped-cavity [внешний источник: Hutchinson ch.15]; геометрия/поля [в корпусе].
Запуск: python3 OUTPUT/pic_spoke33.py (numpy, matplotlib).
"""
import math, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt

m_e=9.1093837015e-31; e=1.602176634e-19; c=2.99792458e8; eps0=8.8541878128e-12
rk,ra=40e-6,80e-6; d=ra-rk; lnAB=math.log(ra/rk)
B=1.1; f0=1.0e12; w0=2*math.pi*f0
Q=700.; gamma_f=0.5; Vdc0=230.
wc=e*B/m_e
vd=Vdc0/(B*rk*lnAB)                   # дрейфовая скорость [в корпусе: синхронизм-формула]
Om_d=vd/rk                            # =1,885e11 рад/с (Ω_d в sync_system_analysis.py)
nstar=int(round(w0/Om_d))             # =33 [sync_system_analysis.py]
neff=nstar
Tdec=Q/w0; Cst=eps0*d; k_c=1.0/(w0*Cst*gamma_f*d)
PSI_STAR=5.8
DT=1e-13                              # глобальный шаг, задаётся в run()

def step(x,y,px,py,E0,t,Vdc):
    """Boris-шаг, идентичный верифицированному pic_rf.py/pic_bunch.py."""
    r=np.hypot(x,y); r=np.where(r<1e-12,1e-12,r); th=np.arctan2(y,x)
    Er=Vdc/(r*lnAB); Ep=E0*(rk/r)*np.cos(neff*th-w0*t)
    rx,ry=x/r,y/r
    Fx=-e*(Er*rx-Ep*ry); Fy=-e*(Er*ry+Ep*rx)
    px=px+Fx*DT/2; py=py+Fy*DT/2
    g=np.sqrt(1+(px*px+py*py)/(m_e*c)**2)
    a=wc*DT/g; cs,sn=np.cos(a),np.sin(a)
    px,py=px*cs+py*sn,-px*sn+py*cs
    vt=(-px*np.sin(th)+py*np.cos(th))/(g*m_e)
    Pu=float(np.sum(-e*Ep*vt))        # Вт при q_mac=e (>0 = отдача волне)
    px=px+Fx*DT/2; py=py+Fy*DT/2
    x=x+px/(g*m_e)*DT; y=y+py/(g*m_e)*DT
    return x,y,px,py,Pu

def emit(M,rng,Vdc,bunched_frac=0.0,dpsi=0.6):
    """Hull-старт: зародыш spoke занимает целые периоды моды n* + окно ψ*; фон равномерный."""
    rs=rk*1.001
    nb=int(round(bunched_frac*M))
    cell=rng.integers(0,neff,M)*2*math.pi/neff
    psi=rng.uniform(0,2*math.pi/neff,M)
    if nb: psi[:nb]=PSI_STAR+rng.uniform(0,dpsi,nb)
    th=(cell+psi)/neff
    x=rs*np.cos(th); y=rs*np.sin(th)
    vh=Vdc/(rs*lnAB)/B
    px=np.zeros(M); py=vh*m_e*np.ones(M)
    vr=rng.normal(0,2e4,M)
    px+=vr*m_e*np.cos(th); py+=vr*m_e*np.sin(th)
    return x,y,px,py

def run(Ib_uA,frac,M=256,T_ns=6.0,nt=30000,seed=1,E0_0=1e2,Vdc=Vdc0):
    """T_ns — длительность модельного окна; nt шагов (DT=T/ns/nt).
    ВНИМАНИЕ: DT должен разрешать и циклотрону (wc·DT≪1), и T_dec (E-обновления):
    при T_ns=6нс, nt=3e4 DT=0,2пс → ~560 событий затухания за шаг ⇒ численный коллапс.
    Для корректной огибающей берём T_ns~0,3 нс [вычислено здесь]."""
    global DT
    DT=T_ns*1e-9/nt
    rng=np.random.default_rng(seed)
    x,y,px,py=emit(M,rng,Vdc,bunched_frac=frac)
    rs=rk*1.001
    vd=Vdc/(rs*lnAB)/B
    tau_cloud=(2*math.pi/neff)/vd      # время пролёта одного периода моды
    qm=Ib_uA*1e-6*tau_cloud/M          # заряд макрочастицы [Кл]
    E0=E0_0; tl=[];El=[];Rl=[];Wl=[];Pl=[];Nl=[]
    for i in range(nt):
        x,y,px,py,Pu=step(x,y,px,py,E0,i*DT,Vdc)
        th=np.arctan2(y,x); rr=np.hypot(x,y)
        R=float(np.abs(np.mean(np.exp(1j*neff*th))))
        W=((np.sqrt(1+(px*px+py*py)/(m_e*c)**2)-1)*m_e*c*c/e).mean()
        Pr=Pu*qm/e
        # без max(): затухание в отрицательную область физически = смена фазы на π
        E0=E0+(k_c*Pr-E0/Tdec)*DT
        bad=(rr<rk*0.9)|(rr>ra*1.05); nb=int(bad.sum())
        if nb:
            xn,yn,pn,qn=emit(nb,np.random.default_rng(seed+i),Vdc,bunched_frac=0.0)
            x[bad]=xn;y[bad]=yn;px[bad]=pn;py[bad]=qn
        if i%20==0:
            tl.append(i*DT*1e12); El.append(E0); Rl.append(R); Wl.append(W); Pl.append(Pr); Nl.append(nb)
    print(f"  [run Ib={Ib_uA:g}мкА frac={frac} seed={seed}: потерь/шаг(ср)={np.mean(Nl):.1f}, <Pr>={np.mean(Pl):.3g} Вт]")
    return np.array(tl),np.array(El),np.array(Rl),np.array(Wl),np.array(Pl),qm

def growth(t,E0):
    """Доля времени, за которую E0 растёт (наклон >0 по скользящему окну 2·Tdec)."""
    msk=len(t)//2
    g=np.polyfit(t[:msk],E0[:msk],1)[0]
    g2=np.polyfit(t[msk:],E0[msk:],1)[0]
    return g,g2

if __name__=="__main__":
    print(f"n*={nstar}  Om_d={Om_d:.4g} рад/с  w0/n*={w0/nstar:.4g} рад/с  детюнинг={(w0/nstar-Om_d)/2/math.pi:.3g} Гц")
    print(f"Tdec={Tdec*1e12:.1f} пс  k_c={k_c:.3g}")
    results={}
    fig,ax=plt.subplots(3,1,figsize=(9,11))
    # С1: холодный старт без зародыша (равномерное облако) при разных Ib
    TNS,NST=0.3,30000   # DT=0,01 пс: разрешает и циклотрону, и T_dec [вычислено здесь]
    for Ib in [20.,200.,2000.]:
        t,E,R,W,P,qm=run(Ib,frac=0.0,T_ns=TNS,nt=NST)
        ax[0].plot(t,np.abs(E),label=f"Ib={Ib:g} мкА, без spoke")
        print(f"C1 Ib={Ib:>7g} мкА frac=0: |E0|: {abs(E[0]):.3g}->{abs(E[-1]):.3g} В/м  R: {R.mean():.2f}  <P>={P[P!=0].mean() if (P!=0).any() else 0:.3g} Вт")
    # С2: зародыш spoke 25% при тех же токах
    for Ib in [20.,200.,2000.]:
        t,E,R,W,P,qm=run(Ib,frac=0.25,T_ns=TNS,nt=NST)
        ax[1].semilogy(t,np.maximum(np.abs(E),1e-6),label=f"Ib={Ib:g} мкА, spoke 25%")
        g1,g2=growth(t,np.abs(E))
        print(f"C2 Ib={Ib:>7g} мкА frac=.25: |E0|: {abs(E[0]):.3g}->{np.abs(E).max():.3g}(max)->{abs(E[-1]):.3g} В/м  R: {R[-1]:.2f}  slope1={g1:.3g} slope2={g2:.3g} В/м/нс")
    # С3: поиск порога при фиксированном frac=0.25 (лог-скан Ib)
    thr=[]
    print("  >>> С3f: ЗАМОРОЖЕННОЕ ПОЛЕ (k_c=0, E0_0=1e4) — чистый отклик пучка; q_mac=e (без усиления током)")
    kc_save=k_c
    globals()['k_c']=0.0
    Ib_save=None
    for Ib in [1.0]:   # макрочастица = 1 e; <Pr> при таком запуске даёт P на ватт тока ×Ib[мкА]
        t,E,R,W,P,qm=run(Ib,frac=0.25,seed=2,T_ns=TNS,nt=NST,E0_0=1e4)
        print(f"    frozen Ib={Ib:g} мкА: E0 держится {abs(E[-1]):.3g}; R: {R[0]:.2f}->{R[-1]:.2f}; W: {W[0]:.0f}->{W[-1]:.0f} эВ; <Pr>={np.mean(P):.3g}")
    globals()['k_c']=kc_save
    for Ib in [20.,50.,100.,300.,1000.,3000.,10000.]:
        t,E,R,W,P,qm=run(Ib,frac=0.25,seed=2,T_ns=TNS,nt=NST)
        grow = abs(E[-1])>3*abs(E[0]) and abs(E[-1])>abs(E[len(E)//2])
        thr.append((Ib,abs(E[-1]),grow))
        print(f"C3 Ib={Ib:>7g} мкА seed2: |E0|(0,3нс)={abs(E[-1]):.3g} В/м  {'РОСТ' if grow else 'затухание'}")
        ax[2].semilogx(Ib,abs(E[-1]),'o-',label='|E0|(0,3 нс)')
    ax[0].set_title("C1: равномерное облако, spoke n*=33 — нет старта"); ax[0].legend(); ax[0].set_ylabel("|E0|, В/м")
    ax[1].set_title("C2: зародыш spoke 25% (ψ*=5,8)"); ax[1].legend(); ax[1].set_ylabel("|E0|, В/м (лог)")
    ax[2].set_title("C3: порог самозапуска vs Ib"); ax[2].set_xlabel("Ib, мкА"); ax[2].set_ylabel("|E0|(0,3 нс), В/м")
    for a in ax: a.grid(alpha=.3)
    plt.tight_layout(); plt.savefig("OUTPUT/pic_spoke33.png",dpi=110)
    print("PNG: OUTPUT/pic_spoke33.png")
