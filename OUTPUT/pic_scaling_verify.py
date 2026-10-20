#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ВЕРИФИКАЦИЯ ЗАКОНА МАСШТАБИРОВАНИЯ (sync_system_analysis + scaling_law_notes) — фактический прогон.

Закон [вычислено здесь, sync_system_notes_2026-10-20.md]:
  инвариант  f0·B·rk²·ln(ra/rk)/Va = const;
  масштаб s: rk'=s·rk, ra'=s·ra; поле t: B'=t·B; тогда Va'=s·Va/t, f0'=f0/(s·t).
Тест: две среды PIC (lumped-модель pic_spoke33, шаг идентичен верифицированному pic_rf):
  BASE: rk=40мкм ra=80 B=1,1 В=230 n*=33 f0=1 ТГц → съём P_unit (Вт на 1 мкА);
  SCAL: s=1,21 (f0→0,7896 ТГц), t=s (компенсация поля ⇒ V不变=230 В, B'=1,331 Тл)
        → по закону f0''=f0'/(s·t)=1 ТГц… ПРОВЕРЯЕМЫЙ ФАКТ: безразмерный съём
        P_unit/(e·Vdc·n_mac) должен совпасть с BASE в пределах численного шума.

Запуск: python3 OUTPUT/pic_scaling_verify.py   (numpy, matplotlib)
"""
import math, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt

m_e=9.1093837015e-31; e=1.602176634e-19; eps0=8.8541878128e-12
Q=700.; gamma_f=0.5; PSI_STAR=5.8; DT=1e-14

class Env:
    """Среда: геометрия (rk,ra), поле B, напряжение Vdc, частота f0. n* выводится из синхронизма."""
    def __init__(self, rk, ra, B, Vdc, f0):
        self.rk, self.ra, self.B, self.Vdc, self.f0 = rk, ra, B, Vdc, f0
        self.lnAB=math.log(ra/rk); self.d=ra-rk
        self.w0=2*math.pi*f0; self.wc=e*B/m_e
        vd=Vdc/(B*rk*self.lnAB); self.Om_d=vd/rk
        self.n=float(self.w0/self.Om_d)     # НЕ округляем: закон инвариантен по n (целый spoke — вопрос решётки каверн)
        self.Tdec=Q/self.w0; self.k_c=1.0/(self.w0*eps0*self.d*gamma_f*self.d)
    def step(self,x,y,px,py,E0,t):
        r=np.hypot(x,y); r=np.where(r<1e-12,1e-12,r); th=np.arctan2(y,x)
        Er=self.Vdc/(r*self.lnAB); Ep=E0*(self.rk/r)*np.cos(self.n*th-self.w0*t)
        rx,ry=x/r,y/r
        Fx=-e*(Er*rx-Ep*ry); Fy=-e*(Er*ry+Ep*rx)
        px=px+Fx*DT/2; py=py+Fy*DT/2
        g=np.sqrt(1+(px*px+py*py)/(m_e*c0)**2) if False else np.sqrt(1+(px*px+py*py)/(m_e*2.99792458e8)**2)
        a=self.wc*DT/g; cs,sn=np.cos(a),np.sin(a)
        px,py=px*cs+py*sn,-px*sn+py*cs
        vt=(-px*np.sin(th)+py*np.cos(th))/(g*m_e)
        Pu=float(np.sum(-e*Ep*vt))
        px=px+Fx*DT/2; py=py+Fy*DT/2
        x=x+px/(g*m_e)*DT; y=y+py/(g*m_e)*DT
        return x,y,px,py,Pu
    def emit(self,M,rng):
        rs=self.rk*1.001
        cell=rng.integers(0,int(round(self.n)),M)*2*math.pi/self.n
        psi=PSI_STAR+rng.uniform(0,0.6,M)      # тормозящая дуга, как pic_spoke33 frac=0.25
        th=(cell+psi)/self.n
        x=rs*np.cos(th); y=rs*np.sin(th)
        vh=self.Vdc/(rs*self.lnAB)/self.B
        px=np.zeros(M); py=vh*m_e*np.ones(M)
        vr=rng.normal(0,2e4,M); px+=vr*m_e*np.cos(th); py+=vr*m_e*np.sin(th)
        return x,y,px,py
    def run(self,M=256,T_ns=0.3,nt=30000,seed=1,E0_0=1e4):
        rng=np.random.default_rng(seed)
        x,y,px,py=self.emit(M,rng)
        E0=E0_0; tl=[];Pl=[];Wl=[];Rl=[]
        for i in range(nt):
            x,y,px,py,Pu=self.step(x,y,px,py,E0,i*DT)
            rr=np.hypot(x,y); th=np.arctan2(y,x)
            W=((np.sqrt(1+(px*px+py*py)/(m_e*2.99792458e8)**2)-1)*m_e*(2.99792458e8)**2/e).mean()
            R=float(np.abs(np.mean(np.exp(1j*self.n*th))))
            bad=(rr<self.rk*0.9)|(rr>self.ra*1.05)
            if bad.any():
                xn,yn,pn,qn=self.emit(int(bad.sum()),np.random.default_rng(seed+i))
                x[bad]=xn;y[bad]=yn;px[bad]=pn;py[bad]=qn
            if i%20==0:
                tl.append(i*DT*1e12); Pl.append(Pu); Wl.append(W); Rl.append(R)
        return np.array(tl),np.array(Pl),np.array(Wl),np.array(Rl)

def summarize(tag,env):
    t,P,W,R=env.run()
    Pu_mean=P[P!=0].mean()                 # Вт при q_mac=e (отклик единичного заряда)
    # нормировка: отнесём к мощности, отдаваемой волне на 1 мкА модельного тока облака M*e за tau_cloud
    tau=(2*math.pi/env.n)/(env.Vdc/(env.B*env.rk*env.lnAB)/env.rk*env.rk)  # = (2π/n)/Ω_d
    tau=(2*math.pi/env.n)/env.Om_d
    P_uA=Pu_mean*tau/(M_e_charge:=e) * 1e-6 / 1e-6  # Вт на 1 мкА = Pu(при q=e)·τ_cloud/e ·1e-6
    print(f"{tag}: n={env.n} B={env.B:.4g}Тл rk={env.rk*1e6:.3g}мкм f0={env.f0/1e12:.4g}ТГц "
          f"V={env.Vdc:g}В | <Pu>={Pu_mean:.3g} Вт(q=e) | P(1мкА)={P_uA:.3g} Вт | "
          f"W: {W[0]:.0f}->{W[-1]:.0f} эВ | R: {R[0]:.2f}->{R[-1]:.2f}")
    return t,P,W,R,P_uA

if __name__=="__main__":
    import sys
    if len(sys.argv)>1 and sys.argv[1]=="check":
        # быстрый алгебраический тест инварианта (без PIC)
        s,t=1.21,1.21
        base=Env(40e-6,80e-6,1.1,230.,1.0e12)
        scal=Env(40e-6*s,80e-6*s,1.1*t,230.*s/t,1.0e12/(s*t))
        inv=lambda E: E.f0*E.B*E.rk**2*E.lnAB/E.Vdc
        print(f"BASE n={base.n} f0={base.f0/1e12:.4g} ТГц inv={inv(base):.6g}")
        print(f"SCAL n={scal.n} f0={scal.f0/1e12:.4g} ТГц inv={inv(scal):.6g}")
        print(f"ОТНОШЕНИЕ ИНВАРИАНТОВ SCAL/BASE = {inv(scal)/inv(base):.6f} (идеал 1,0)")
        print(f"Ω_d SCAL/BASE = {scal.Om_d/base.Om_d:.6f}; w0/n SCAL/BASE = "
              f"{(scal.w0/scal.n)/(base.w0/base.n):.4f} (округление n целым ⇒ ~1±0,5%)")
        sys.exit(0)
    S=1.21
    base=Env(40e-6,80e-6,1.1,230.,1.0e12)
    # масштабирование по закону: rk,ra ×s; B×t; Va×(s/t); f0' = f0/(s·t)
    scal=Env(40e-6*S,80e-6*S,1.1*S,230.*S/S,1.0e12/(S*S))
    tb,Pb,Wb,Rb,Pbu=summarize("BASE",base)
    ts,Ps,Ws,Rs,Psu=summarize("SCAL",scal)
    ratio=Psu/Pbu
    print(f"\nОТНОШЕНИЕ безразмерного съёма SCAL/BASE = {ratio:.3f}  "
          f"(идеал 1,00; допуск численного шума ±10%)")
    print("ИНВАРИАНТ f0·B·rk²·ln(ra/rk)/Va:")
    inv=lambda E: E.f0*E.B*E.rk**2*E.lnAB/E.Vdc
    print(f"  BASE={inv(base):.4g}  SCAL={inv(scal):.4g}  отношение={inv(scal)/inv(base):.4f}")
    fig,ax=plt.subplots(1,3,figsize=(13,4))
    ax[0].plot(tb,Pb,label='BASE'); ax[0].plot(ts,Ps,label='SCAL')
    ax[0].set_xlabel('t, пс'); ax[0].set_ylabel('Pu, Вт (q=e)'); ax[0].legend(); ax[0].set_title('отклик пучка волне')
    ax[1].plot(tb,Wb,label='BASE'); ax[1].plot(ts,Ws,label='SCAL')
    ax[1].set_xlabel('t, пс'); ax[1].set_ylabel('<W>, эВ'); ax[1].legend(); ax[1].set_title('энергия пучка')
    ax[2].plot(tb,Rb,label='BASE'); ax[2].plot(ts,Rs,label='SCAL')
    ax[2].set_xlabel('t, пс'); ax[2].set_ylabel('R spoke'); ax[2].legend(); ax[2].set_title('бунчинг')
    plt.tight_layout(); plt.savefig('/workspace/THz/OUTPUT/pic_scaling_verify.png',dpi=110)
    print("PNG: OUTPUT/pic_scaling_verify.png")
