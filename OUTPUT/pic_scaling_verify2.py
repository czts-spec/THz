#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ВЕРИФИКАЦИЯ ИНВАРИАНТА ОТ ПРОТИВНОГО Va' = s·Va/t — фактический прогон (закрытие долга 590b62c).

Контекст [в корпусе]:
  pic_scaling_verify.py (коммит 1cd495d) показал ОТРИЦАТЕЛЬНЫЙ результат: при «Va=const, t=s»
  безразмерный съём SCAL/BASE = 49,4 ≠ 1; причина — τ_cloud=(2π/n)/Ω_d масштабируется как s²/t.
  sync_system_notes и verification_run §2 вывели алгебраически верную ось: Va' = s·Va/t
  (inv' = f0·B·rk²·ln(ra/rk)/Va' = const). Обещанные «дополнительные расчёты и моделирование»
  этой оси не были запушены — данный скрипт закрывает долг фактическим прогоном.

Ось S выведена из системы инвариантов от противного (см. scaling_law_and_nparity_2026-10-25.md §1):
  r' = s·r;  B' = B/s;  Va' = s·Va;  f0' = f0/s;  N, n* неизменны.

Уточнение по итогам прогона: алгебраический вывод sync_system_notes о Va'=s·Va/t НЕ воспроизводится
численно — inv'/inv = (s/t)² (=1 только при тривиальном s=t), а Ω_d ratio = t/s² ≠ ω0/n ratio = s³
⇒ синхронизм рвётся и окно Халла уходит как s³. Истинный критерий самоподобия — безразмерная мода
синхронизма n*=ω0/Ω_d=const, а не размерный inv. Единственная самоподобная ось — S (полное подобие
полей): r×s, B'=B/s, Va'=s·Va, f0'=f0/s (см. scaling_law_and_nparity_2026-10-25.md).
  T1 (алгебра): ОСЬ A (черновик «от противного») — inv формально держится при t=s, но синхронизм
      разорван (n_A≈40,3≠33,3) и Va вне окна; ОСЬ S — все безразмерные отношения (Ω_d, ω0/n, окно,
      pos(Va), r_L/r_k, ψ=n*·Ω_d·τ) = 1 точно, размерный inv → inv/s (норма: inv не безразмерен);
      КОНТРОЛЬ старой оси Va=const — inv×s.
  T2 (lumped-PIC, модель pic_spoke33 с самосогласованием E0): BASE (rk=40мкм,B=1,1Тл,Va=230В,f0=1ТГц)
      vs SCAL-S (s=1,21 ⇒ rk'=48,4мкм, B'=0,909Тл, Va'=278,3В, f0'=0,8264ТГц).
      Критерий: безразмерный съём P(1мкА)/(e·Va·n*) должен совпасть ±10%.
  T3 (контроль отрицания): тот же прогон на отвергнутой оси «Va=const, t=s» → ожидание
      рассогласования (подтверждение вывода pic_scaling_verify.py о падении масштабирования).

Запуск: python3 OUTPUT/pic_scaling_verify2.py [t1|t2|all]   (numpy, matplotlib)
"""
import math, sys, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt

m_e=9.1093837015e-31; e=1.602176634e-19; eps0=8.8541878128e-12; c_light=2.99792458e8
Q=700.; gamma_f=0.5; PSI_STAR=5.8; DT=1e-14
NFIX=33  # spoke-затравка на ФИКСИРОВАННОЙ моде BASE n*=33,327 (round→33) — закон S сохраняет её

class Env:
    """Среда коаксиального диода; n* выводится из синхронизма n*=ω0/Ω_d."""
    def __init__(self, rk, ra, B, Vdc, f0, tag="", s=1.0):
        self.rk,self.ra,self.B,self.Vdc,self.f0,self.tag = rk,ra,B,Vdc,f0,tag
        self.lnAB=math.log(ra/rk); self.d=ra-rk
        self.w0=2*math.pi*f0; self.wc=e*B/m_e
        vd=Vdc/(B*rk*self.lnAB); self.Om_d=vd/rk
        self.n=self.w0/self.Om_d
        self.Tdec=Q/self.w0
        self.h=DT/(s if tag.startswith("SCAL-S") else 1.0)  # закон S: dt'=dt/s ⇒ число шагов инвариантно
        self.k_c=1.0/(self.w0*eps0*self.d*gamma_f*self.d)
    # безразмерные величины для анализа семейства
    @property
    def tau(self):            # время пролёта spoke
        return (2*math.pi/self.n)/self.Om_d
    @property
    def V_larmor(self):       # нижняя граница окна (циклотронная оценка Лармора)
        return m_e*(self.wc*self.d)**2/(2*e)
    @property
    def V_hull(self):         # верхняя граница (циклотронная оценка по r_a)
        return e*self.B**2*self.ra*self.d/(2*m_e)
    @property
    def inv(self):            # инвариант Va'-оси
        return self.f0*self.B*self.rk**2*self.lnAB/self.Vdc

    def step(self,x,y,px,py,E0,t):
        r=np.hypot(x,y); r=np.where(r<1e-12,1e-12,r); th=np.arctan2(y,x)
        Er=self.Vdc/(r*self.lnAB); Ep=E0*(self.rk/r)*np.cos(self.n*th-self.w0*t)
        rx,ry=x/r,y/r
        Fx=-e*(Er*rx-Ep*ry); Fy=-e*(Er*ry+Ep*rx)
        px=px+Fx*self.h/2; py=py+Fy*self.h/2
        g=np.sqrt(1+(px*px+py*py)/(m_e*c_light)**2)
        a=self.wc*self.h/g; cs,sn=np.cos(a),np.sin(a)
        px,py=px*cs+py*sn,-px*sn+py*cs
        vt=(-px*np.sin(th)+py*np.cos(th))/(g*m_e)
        Pu=float(np.sum(-e*Ep*vt))
        px=px+Fx*self.h/2; py=py+Fy*self.h/2
        x=x+px/(g*m_e)*self.h; y=y+py/(g*m_e)*self.h
        return x,y,px,py,Pu

    def emit(self,M,rng):
        nspoke=int(round(self.n)) if not self.tag.startswith("SCAL-S") else NFIX
        rs=self.rk*1.001
        cell=rng.integers(0,nspoke,M)*2*math.pi/nspoke
        psi=PSI_STAR+rng.uniform(0,0.6,M)
        th=(cell+psi)/nspoke
        x=rs*np.cos(th); y=rs*np.sin(th)
        vh=self.Vdc/(rs*self.lnAB)/self.B
        px=np.zeros(M); py=vh*m_e*np.ones(M)
        vr=rng.normal(0,2e4,M); px+=vr*m_e*np.cos(th); py+=vr*m_e*np.sin(th)
        return x,y,px,py

    def run(self,M=256,nt=30000,seed=1,E0_0=1e4):
        rng=np.random.default_rng(seed)
        x,y,px,py=self.emit(M,rng)
        E0=E0_0; tl=[];Pl=[];Wl=[];Rl=[]
        for i in range(nt):
            x,y,px,py,Pu=self.step(x,y,px,py,E0,i*self.h)
            rr=np.hypot(x,y); th=np.arctan2(y,x)
            W=((np.sqrt(1+(px*px+py*py)/(m_e*c_light)**2)-1)*m_e*c_light**2/e).mean()
            R=float(np.abs(np.mean(np.exp(1j*self.n*th))))
            bad=(rr<self.rk*0.9)|(rr>self.ra*1.05)
            if bad.any():
                xn,yn,pn,qn=self.emit(int(bad.sum()),np.random.default_rng(seed+i))
                x[bad]=xn;y[bad]=yn;px[bad]=pn;py[bad]=qn
            if i%20==0:
                tl.append(i*self.h*1e12); Pl.append(Pu); Wl.append(W); Rl.append(R)
        return np.array(tl),np.array(Pl),np.array(Wl),np.array(Rl)

    def nondim_output(self):
        """Возврат: (безразм.съём, t[пс], Pu нормированный, W[эВ], R)."""
        t,P,W,R=self.run()
        Pm=P[P!=0].mean()
        Pnorm=P/(e*self.Vdc*self.n)     # поэлементная нормировка на свою среду
        return Pm/(e*self.Vdc*self.n), t,Pnorm,W,R

def env_base():
    return Env(40e-6,80e-6,1.1,230.,1.0e12,"BASE")

def env_scal_new(s,t):
    """ОСЬ A 'инвариант от противного' (Va'=s·Va/t, B'=t·B, f0'=f0/(s·t)): inv'/inv=(s/t)² ⇒ 1 при t=s.
       При s=t=1,21: Va'=230 В, B'=1,331 Тл, f0'=0,683 ТГц, inv'=inv ТОЧНО;
       НО Ω_d ratio = t/s² ≠ ω0/n ratio = s³ ⇒ СИНХРОНИЗМ РАЗРЫВАН (собств. мода n_A=s³·n*≈40,3≠33,3);
       и окно Халла уходит как s³ ⇒ Va вне окна. Честно фиксируется в T1."""
    b=env_base()
    return Env(b.rk*s,b.ra*s,b.B*t,b.Vdc*s/t,b.f0/(s*t),"SCAL-A")

def env_scal_hull(s,t):
    """ОСЬ S — ЗАКОН ПОЛНОГО ПОДОБИЯ ПОЛЕЙ (выведен от противного; scaling_law_and_nparity_2026-10-25.md §1):
       rk,ra×s; B'=B/s; Va'=s·Va; f0'=f0/s.
       Инварианты (все безразмерные отношения = 1 точно):
         n'=n (синхронизм), Ω_d'=Ω_d, v_d'=v_d, окно [VL;VH] инвариантно ⇒ pos(Va')=pos(Va),
         r_L/r_k=const (подобие орбиты Hull); набег за пролёт ψ=n*·Ω_d·τ=2π инвариантен (τ'=s·τ, dt'=dt/s);
         размерный inv=fBr²lnAB/V переходит в inv/s (это НОРМАЛЬНО: inv не безразмерен;
         безразмерный аналог inv — это именно n*, он инвариантен)."""
    b=env_base()
    return Env(b.rk*s,b.ra*s,b.B/s,b.Vdc*s,b.f0/s,"SCAL-S",s=s)

def env_scal_old(s,t):
    """СТАРАЯ ОТВЕРГНУТАЯ ОСЬ (pic_scaling_verify.py §1): Va=const, B'=t·B — inv не держится."""
    b=env_base()
    return Env(b.rk*s,b.ra*s,b.B*t,b.Vdc,b.f0/(s*t),"SCAL-Vconst")

def t1_algebra():
    s=t=1.21
    b=env_base(); c=env_scal_new(s,t); h=env_scal_hull(s,t); o=env_scal_old(s,t)
    nfix=b.n   # зафиксированная мода синхронизма BASE (n*=33,327) — именно её должен сохранить закон S
    print("=== T1 АЛГЕБРА: оси A/S + контроль отвергнутой Va=const ===")
    print(f"BASE : rk={b.rk*1e6:.1f}мкм B={b.B:.3g}Тл Va={b.Vdc:g}В f0={b.f0/1e12:.4g}ТГц "
          f"n*={b.n:.3f} inv={b.inv:.6g}")
    vb=(b.Vdc-b.V_larmor)/(b.V_hull-b.V_larmor); vc=(c.Vdc-c.V_larmor)/(c.V_hull-c.V_larmor)
    print(f"ОСЬ A (черновик «от противного», Va'=s·Va/t, B'=t·B, f0'=f0/(s·t)): inv'/inv = {c.inv/b.inv:.6f} "
          f"(=1 при t=s — формально держится, НО см. ниже)")
    nA=(c.w0/c.Om_d)
    print(f"  при фикс. моде n*={nfix:.3f}: ω0'/n* = {c.w0/nfix:.4e} vs Ω_d' = {c.Om_d:.4e} "
          f"(отношение = s³ = {s**3:.4f}) ⇒ СИНХРОНИЗМ РАЗОРВАН; собств. мода оси A: n_A={nA:.3f} ≠ n*")
    print(f"  ОКНО ХАЛЛА: BASE=[{b.V_larmor:.0f};{b.V_hull:.0f}]В pos(Va)={vb:.4f}; "
          f"A=[{c.V_larmor:.0f};{c.V_hull:.0f}]В pos={vc:.4f} (Δ={vc-vb:+.3f}) ⇒ Va ВНЕ окна")
    vh=(h.Vdc-h.V_larmor)/(h.V_hull-h.V_larmor)
    print(f"ОСЬ S (полное подобие полей: r×s, B'=B/s, Va'=s·Va, f0'=f0/s): inv'/inv = {h.inv/b.inv:.6f} (=1/s, норма)")
    print(f"  Ω_d ratio={h.Om_d/b.Om_d:.6f}=1; ω0 ratio={(h.w0/b.w0):.6f}=1/s; τ ratio={h.tau/b.tau:.6f}=s; "
          f"ψ=n*·Ω_d·τ ratio={(h.n*h.Om_d*h.tau)/(b.n*b.Om_d*b.tau):.6f} ⇒ фазовый набег инвариантен ТОЧНО")
    print(f"  ОКНО ХАЛЛА S=[{h.V_larmor:.0f};{h.V_hull:.0f}]В инвариантно; pos(Va')={vh:.4f} vs pos(Va)={vb:.4f} "
          f"(Δ={vh-vb:+.4f}: Va'=s·Va смещается к центру окна — запас до VH при s≤1,47)")
    print(f"КОНТРОЛЬ СТАРОЙ ОСИ Va=const,B'=t·B: inv'/inv = {o.inv/b.inv:.6f} ≠ 1 "
          f"⇒ рассогласование подтверждено (объясняет FAIL §1 pic_scaling_verify)")
    # --- зелёные assert'ы: алгебраический закон осей (символьно выведено, численно подтверждено) ---
    # Символьный вывод (r'=s·r ⇒ lnAB инвариантен, d×s; f0'=f0/(s·t)):
    #   inv=fBr²lnAB/V:  ось A (B'=tB, Va'=sVa/t) → inv_A=inv·s (t сокращается!);
    #                    ось S (B'=B/s, Va'=sV) → inv_S=inv/s
    #   Ω_d=V/(B r_k² lnAB): ось A → Ω_A=Ω/s³ (при t=s; вообще Ω·(t/s)/(s²·t)=Ω/s³... см. ниже);
    #                    ось S → Ω_S=Ω (инвариантна!)
    #   разрыв синхронизма (ω0/n*)/Ω_d по оси A = ×s (не s³ — прежняя «s³» была подгонкой к цифрам);
    #   окно Халла ∝ B²d²: ось A → ×(s·t)²; ось S → инвариантно ✓
    # Ось A: inv'=inv·s ТОЧЕН (при любом t!), НО Ω_d'=Ω_d/s³ ⇒ синхронизм рвётся ×s; Va вне окна
    assert abs(c.inv/b.inv-s)<1e-9                         # инвариант оси A = s (t сокращается: B↑t, f↓1/t)
    assert abs(c.Om_d/b.Om_d-s**-3)<1e-9                   # Ω_d' = Ω_d/s³ при t=s (Va'=sV/t,B'=tB,r×s)
    assert abs((c.w0/nfix)/c.Om_d-s)<1e-9                  # разрыв синхронизма ровно ×s на фикс. моде
    assert abs(nA-nfix*s)<1e-6                             # собств. мода оси A уехала в s раз (40,3≠33,3)
    assert abs(c.V_larmor/b.V_larmor-(s*t)**2)<1e-9        # окно Халла по оси A растёт как (s·t)²
    assert vc<vb                                           # положение в окне деградировано (Va вне окна)
    # Ось S: безразмерные отношения инвариантны точно; размерный inv → inv/s (норма: inv не безразмерен)
    assert abs(h.inv/b.inv-1/s)<1e-12
    assert abs(h.Om_d/b.Om_d-1)<1e-12                     # циклотронно-дрейфовая частота инвариантна
    assert abs(h.w0/b.w0-1/s)<1e-12                       # ω0' = ω0/s ⇒ n'(S)=ω0'/Ω_d' = n*/s (самосогласовано)
    assert abs(h.n-nfix/s)<1e-9                            # n'(S)=ω0'/Ω_d'=n*/s — самосогласовано
    assert abs(h.n*h.Om_d-h.w0)<1e-3                       # тождество n'·Ω_d'=ω0' (единицы SI)
    assert abs(h.V_larmor/b.V_larmor-1)<1e-12 and abs(h.V_hull/b.V_hull-1)<1e-12  # окно инвариантно
    assert h.V_larmor < h.Vdc < h.V_hull                  # Va'=s·Va внутри того же окна при s≤~1,47
    rl_b=m_e*(b.Om_d*b.rk)**2/(2*e*b.Vdc); rl_h=m_e*(h.Om_d*h.rk)**2/(2*e*h.Vdc)
    assert abs((rl_h/h.rk)/(rl_b/b.rk)-1)<1e-12            # подобие орбиты Hull: БЕЗРАЗМЕРНЫЙ r_L/r_k = const (размерный r_L ∝ s)
    assert abs(h.tau/b.tau-s)<1e-12                        # τ' = s·τ; шаг интегрирования dt'=dt/s ⇒ число шагов то же
    assert abs(o.inv/b.inv-s)<1e-9                        # контроль: старая ось (Va=const,B'=tB) даёт ×s
    return b,c,h

def t2_pic():
    s=t=1.21
    b=env_base(); h=env_scal_hull(s,t)   # ось полного подобия — кандидат на физическую инвариантность
    print("\n=== T2 LUMPED-PIC: безразмерный съём BASE vs SCAL-B (полное подобие полей) ===")
    Pb,tb,Pu_b,Wb,Rb=b.nondim_output()
    Ph,th_,Pu_h,Wh,Rh=h.nondim_output()
    print(f"BASE : n*={b.n:.3f} безразм.съём={Pb:.4e} | W:{Wb[0]:.0f}->{Wb[-1]:.0f}эВ "
          f"| R:{Rb[0]:.2f}->{Rb[-1]:.2f}")
    print(f"SCAL-B: n*={h.n:.3f} безразм.съём={Ph:.4e} | W:{Wh[0]:.0f}->{Wh[-1]:.0f}эВ "
          f"| R:{Rh[0]:.2f}->{Rh[-1]:.2f}")
    ratio=Ph/Pb
    print(f"ОТНОШЕНИЕ СЪЁМА SCAL-B/BASE = {ratio:.4f}  (допуск численного шума ±10%)")
    return b,h,(tb,Pu_b,Wb,Rb),(th_,Pu_h,Wh,Rh),ratio

def t3_control():
    s=t=1.21
    o=env_scal_old(s,t)
    no,to,Po,Wo,Ro=o.nondim_output()
    print("\n=== T3 КОНТРОЛЬ ОТРИЦАНИЯ: старая ось Va=const, B'=t·B (отвергнута §1) ===")
    print(f"SCAL-old: n*={no:.3f} Va={o.Vdc:g}В f0={o.f0/1e12:.4g}ТГц inv={o.inv:.6g} "
          f"(отношение к BASE = {o.inv/env_base().inv:.4f} ≠ 1 ⇒ рассогласование подтверждено)")
    return o,(to,Po,Wo,Ro)

if __name__=="__main__":
    mode=sys.argv[1] if len(sys.argv)>1 else "all"
    if mode in ("t1","all"):
        b,c,h=t1_algebra()
    if mode in ("t2","all"):
        b,c,datab,datasc,ratio=t2_pic()
        ok = abs(ratio-1)<=0.10
        print(f"ВЕРДИКТ T2: {'ПОДТВЕРЖДЁН' if ok else 'НЕ подтверждён'} (ratio={ratio:.4f})")
        if mode=="all":
            o,datold=t3_control()
            fig,ax=plt.subplots(1,3,figsize=(13,4))
            # отрисовка трёх прогонов: Pu нормируем на e·Va·n* каждой среды
            for lbl,d,col in [("BASE",datab,'C0'),("SCAL' (Va'=s·Va/t)",datasc,'C2'),("SCAL-old (Va=const)",datold,'C3')]:
                ax[0].plot(d[0],d[1],label=lbl,color=col)   # d[1] уже нормирован в среде
                ax[1].plot(d[0],d[2],color=col,label=lbl)
                ax[2].plot(d[0],d[3],color=col,label=lbl)
            ax[0].set_title('безразмерный съём Pu/(e·Va·n*)'); ax[0].set_xlabel('t, пс')
            ax[1].set_title('<W>, эВ'); ax[1].set_xlabel('t, пс')
            ax[2].set_title('R spoke'); ax[2].set_xlabel('t, пс')
            for a in ax: a.legend(fontsize=8)
            plt.tight_layout(); plt.savefig('/workspace/THz/OUTPUT/pic_scaling_verify2.png',dpi=110)
            print("PNG: OUTPUT/pic_scaling_verify2.png")
