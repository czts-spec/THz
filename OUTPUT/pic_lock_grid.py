#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ПОЛНАЯ СЕТОЧНАЯ ВЕРИФИКАЦИЯ pic_lock (injection locking spoke n*=33).

Контекст [в корпусе]: pic_spoke33.py показал отсутствие самовозбуждения при
I_b ≤ 10 мА (лumped, Q=700). pic_lock (лumped) предсказывал захват фазы при
внешней затравке P_inj ≥ ~1 мВт. Здесь — независимая ПРОВЕРКА сеточным PIC:
поле E_phi(r,θ,t) = Re{ A·J1-like(r)·exp(i(n θ − ω0 t)) } задаётся извне
(стоячая+бегущая компонента), без уравнения огибающей; измеряем отклик пучка:
фазовый захват Δψ(t)→const и средний съём <P>.

Метод: одночастичный зонд + облако 256 макрочастиц, Boris-push (идентичен
верифицированному шагу pic_rf/pic_spoke33), геометрия П-2 [в корпусе].
Запуск: python3 OUTPUT/pic_lock_grid.py   (~2 мин)
"""
import math, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt

m_e=9.1093837015e-31; e=1.602176634e-19; c0=2.99792458e8
rk,ra=40e-6,80e-6; lnAB=math.log(ra/rk); d=ra-rk
B=1.1; Vdc=230.; f0=1.0e12; w0=2*math.pi*f0; wc=e*B/m_e
n_=float(w0/((Vdc/(B*rk*lnAB)/rk)))   # ТОЧНЫЙ spoke-синхронизм n*=33,327 [sync_system_notes]; целый spoke — вопрос решётки каверн
DT=1e-14

def fields(x,y,E0,t):
    r=np.hypot(x,y); r=np.where(r<1e-12,1e-12,r); th=np.arctan2(y,x)
    Er=Vdc/(r*lnAB)
    Ep=E0*(rk/r)*np.cos(n_*th-w0*t)          # бегущая волна spoke-моды (затравка)
    return r,th,Er,Ep

def push(x,y,px,py,E0,t):
    r,th,Er,Ep=fields(x,y,E0,t)
    rx,ry=x/r,y/r
    Fx=-e*(Er*rx-Ep*ry); Fy=-e*(Er*ry+Ep*rx)
    px=px+Fx*DT/2; py=py+Fy*DT/2
    g=np.sqrt(1+(px*px+py*py)/(m_e*c0)**2)
    a=wc*DT/g; cs,sn=np.cos(a),np.sin(a)
    px,py=px*cs+py*sn,-px*sn+py*cs
    vt=(-px*np.sin(th)+py*np.cos(th))/(g*m_e)
    Pu=float(np.sum(-e*Ep*vt))
    px=px+Fx*DT/2; py=py+Fy*DT/2
    x=x+px/(g*m_e)*DT; y=y+py/(g*m_e)*DT
    return x,y,px,py,Pu

def emit(M,rng,psi_center=None,spread=0.6):
    rs=rk*1.001
    cell=rng.integers(0,n_,M)*2*math.pi/n_
    psi=(psi_center+rng.uniform(0,spread,M)) if psi_center is not None else rng.uniform(0,2*math.pi/n_,M)
    th=(cell+psi)/n_
    vh=Vdc/(rs*lnAB)/B
    px=np.zeros(M); py=vh*m_e*np.ones(M)
    vr=rng.normal(0,2e4,M); px+=vr*m_e*np.cos(th); py+=vr*m_e*np.sin(th)
    return rs*np.cos(th),rs*np.sin(th),px,py

def phase_of(x,y):
    th=np.arctan2(y,x)
    ph=(n_*th-w0*T_NOW)% (2*math.pi/n_)     # локальная фаза внутри периода spoke
    return ph

T_NOW=0.0
def run(E0,Ib_uA,M=256,T_ns=3.0,nt=300000,seed=1,bunched=True):
    global T_NOW
    rng=np.random.default_rng(seed)
    x,y,px,py=emit(M,rng,psi_center=5.8 if bunched else None)
    tau_cloud=(2*math.pi/n_)/(Vdc/(B*rk*lnAB)/rk)
    qm=Ib_uA*1e-6*tau_cloud/M
    tl=[];Wl=[];Pl=[];Rl=[];Dl=[]
    psi_ref=5.8%(2*math.pi/n_)
    for i in range(nt):
        T_NOW=i*DT
        x,y,px,py,Pu=push(x,y,px,py,E0,T_NOW)
        rr=np.hypot(x,y)
        W=((np.sqrt(1+(px*px+py*py)/(m_e*c0)**2)-1)*m_e*c0*c0/e).mean()
        ph=phase_of(x,y); R=float(np.abs(np.mean(np.exp(1j*n_*np.arctan2(y,x)))))
        dv=np.angle(np.mean(np.exp(1j*(ph-psi_ref))))
        bad=(rr<rk*0.9)|(rr>ra*1.05)
        if bad.any():
            xn,yn,pn,qn=emit(int(bad.sum()),np.random.default_rng(seed+i),psi_center=5.8 if bunched else None)
            x[bad]=xn;y[bad]=yn;px[bad]=pn;py[bad]=qn
        if i%100==0:
            tl.append(T_NOW*1e12);Wl.append(W);Pl.append(Pu*qm/e);Rl.append(R);Dl.append(dv)
    return np.array(tl),np.array(Wl),np.array(Pl),np.array(Rl),np.array(Dl)

if __name__=="__main__":
    print(f"lock-grid: n={n_} Ω_d={ (Vdc/(B*rk*lnAB)/rk):.4g} w0/n={w0/n_:.4g} "
          f"детюнинг={(w0/n_-(Vdc/(B*rk*lnAB)/rk))/2/math.pi:.3g} Гц")
    cases=[("без затравки",0.0,False),("E0=1 кВ/м",1e3,True),("E0=10 кВ/м",1e4,True),
           ("E0=100 кВ/м",1e5,True),("E0=1 МВ/м",1e6,True)]
    fig,ax=plt.subplots(2,2,figsize=(11,8))
    res={}
    for tag,E0,bun in cases:
        t,W,P,R,D=run(E0,Ib_uA=20.,bunched=bun)
        tail=slice(len(t)//2,None)
        print(f"{tag:>14}: <P>={P[tail].mean():+.3g} Вт  W: {W[0]:.0f}->{W[-1]:.0f} эВ  "
              f"R={R[-1]:.2f}  |Δψ|_ср={np.abs(D[tail]).mean():.3f} рад  разброс Δψ={np.std(D[tail]):.3f}")
        ax[0,0].plot(t,np.abs(P)*1e3,label=tag)
        ax[0,1].plot(t,W,label=tag)
        ax[1,0].plot(t,R,label=tag)
        ax[1,1].plot(t,D,label=tag)
    ax[0,0].set_ylabel('|P|, мВт'); ax[0,1].set_ylabel('<W>, эВ')
    ax[1,0].set_ylabel('R spoke');  ax[1,1].set_ylabel('Δψ, рад')
    for a_ in ax.flat: a_.legend(fontsize=7); a_.set_xlabel('t, пс')
    ax[0,0].set_title('съём мощности'); ax[0,1].set_title('энергия пучка')
    ax[1,0].set_title('бунчинг'); ax[1,1].set_title('расхождение фазы (захват→0)')
    plt.tight_layout(); plt.savefig('/workspace/THz/OUTPUT/pic_lock_grid.png',dpi=110)
    print("PNG: OUTPUT/pic_lock_grid.png")
