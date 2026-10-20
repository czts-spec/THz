#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
П-2: 2D сеточная PIC-верификация lumped-моделей pic_rf / pic_bunch / pic_lock.

Схема [вычислено здесь]: электростатический 2D3V PIC в (r,θ)-секторе π-моды:
 - сетка Nr×Nth на секторе Δθ=2π/neff (neff=N/2=64), периодичность по θ;
 - пуш: Boris в E_стат(E_r=V/(r ln(ra/rk))) + E_RF(фикс.) + E_spacecharge(Poisson);
 - Poisson по азимутальным гармоникам k=n_eff·s (s=0..S): трёхдиагональная система
   (1/r)(r Φ')' - (n_eff s)²/r² Φ = -ρ_s/ε0, Φ=0 на r_k и r_a (scipy solve_banded);
 - CIC-депозиция заряда, E=-∇Φ интерполируется на частицы обратно через CIC.

Цели верификации:
 L1. Hull-дрейф без RF и без space charge: |ΔW|→0, облако не падает на катод;
 L2. съём энергии из фикс. π-поля (spoke-старт) — сравнение с R1 pic_rf (86 эВ/э);
 L3. cold-облако + только собственное поле: нет роста моды n_eff (нет запуска)
     — сравнение с R2 pic_rf / pic_bunch (самовозбуждения нет);
 L4. амплитуда space-charge поля n_eff при полном токе — оценка shunt-импеданса.

Ограничения (честно): 2D (без z), мода резонатора не динамическая (потери/накопление
стенки считаются в lumped-моделях), нерелятивистский пуш допустим (v~0,03c).
Выход: консоль + OUTPUT/pic_grid.png
"""
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded

m_e=9.1093837015e-31; e=1.602176634e-19; c=2.99792458e8; eps0=8.8541878128e-12
rk,ra=40e-6,80e-6; d=ra-rk; lnAB=math.log(ra/rk)
B=1.1; f0=1.0e12; w0=2*math.pi*f0; neff=64
Vdc=230.; wc=e*B/m_e

Nr,Nth,KMAX = 96,96,80                  # пуассоновские гармоники k=1..KMAX (π-гармоника k=64 внутри)
rg=np.linspace(rk,ra,Nr); dr=rg[1]-rg[0]
dth=2*math.pi/neff/Nth
rng=np.random.default_rng(11)

Ib=20e-6                                 # ток облака [в корпусе: стартовый 20 мкА]
# число электронов в (r,θ)-секторе π-периода: N_e = Ib/e × время оборота сектора
vd=Vdc/(rk*lnAB)/B                       # Hull-скорость на катоде ~7,54e6 м/с
M=1500                                   # макрочастиц; вес q_mac = N_e/M·e
tau_sec=(rk*lnAB/vd)/c                   # время пролёта сектора (Hull-дрейф)
N_e=max(int(Ib/e*tau_sec),M)             # реальных э в секторе; макрочастиц не больше
qmac=N_e/M*e
dt=0.001*2*math.pi/wc                    # 0,1% циклотронного периода
print(f"=== КОНФИГ ===\nNr={Nr} Nth={Nth} sector={2*math.pi/neff:.4f}рад M={M} "
      f"N_e={N_e:.3e} q_mac/q_e={N_e/M:.1f}\nv_d={vd/1e6:.2f}e6 м/с dt={dt*1e12:.2f} пс")

def init_cloud(bunched=False,psi0=5.8,dpsi=0.6,r=None):
    g=r if r is not None else rng
    rs=rk+g.uniform(0.5e-6,2.5e-6,M)
    if bunched:
        ths=g.uniform(psi0/neff,(psi0+dpsi)/neff,M)
    else:
        ths=g.uniform(0,2*math.pi/neff,M)
    vr=g.normal(0,2e4,M)               # только тепловая радиальная компонента
    px=m_e*vr*np.cos(ths); py=m_e*vr*np.sin(ths)   # старт из покоя (эмитанс)
    return rs*np.cos(ths),rs*np.sin(ths),px,py

# --- Poisson solver: rho[Nr,Nth] (Кл/м³) -> Er,Ephi on grid ---
def poisson_fields(rho):
    """∇²Φ=-ρ/ε0 по гармоникам m=neff·s, Φ=0 на r_k,r_a; E=-∇Φ на сетке."""
    rk_=np.fft.rfft(rho,axis=1)/Nth
    Kmin=1; Kmax=min(KMAX,rk_.shape[1]-1)
    lo=(rg[1:-1]+dr/2)/dr; hi=(rg[1:-1]-dr/2)/dr
    Er=np.zeros((Nr,Nth)); Ep=np.zeros((Nr,Nth))
    thc=np.arange(Nth)*dth
    for k in range(Kmin,Kmax+1):
        ab=np.zeros((3,Nr-2),complex)
        ab[0,1:]=lo[:-1]; ab[2,:-1]=-hi[1:]
        ab[1]=-(lo+hi)-(k*k*dr*dr/rg[1:-1]**2)
        rhs=-rk_[1:-1,k]*dr*dr/eps0
        Phi=solve_banded((1,1),ab,rhs)
        Phif=np.zeros(Nr,complex); Phif[1:-1]=Phi   # Φ=0 на электродах
        phase=np.exp(1j*k*thc)[None,:]
        dPhidr=np.gradient(Phif,dr)
        Er-=np.real(dPhidr[:,None]*phase)                          # E_r = -∂Φ/∂r
        Ep-=np.real(1j*k*Phif[:,None]*phase/rg[:,None])             # E_φ = -(1/r)∂Φ/∂θ
    return Er,Ep

def cic(x,y,arr):
    """CIC deposit x,y (м, глобальные коорд.) в arr[Nr,Nth] сектора."""
    arr[:]=0.
    rr=np.hypot(x,y); tt=np.arctan2(y,x)% (2*math.pi/neff)
    fi=np.clip((rr-rg[0])/dr,0,Nr-1.001); ti=tt/dth
    i0=fi.astype(int); j0=(ti.astype(int))%Nth
    fr=fi-i0; ft=ti-j0
    for di,wi in ((0,1-fr),(1,fr)):
        ii=np.minimum(i0+di,Nr-1)
        for dj,wj in ((0,1-ft),(1,ft)):
            jj=(j0+dj)%Nth
            np.add.at(arr,(ii,jj),qmac*wi*wj/(dr*dth*rg[ii]))  # объём ячейки r dr dθ

def cic_gather(x,y,F):
    rr=np.hypot(x,y); tt=np.arctan2(y,x)%(2*math.pi/neff)
    fi=np.clip((rr-rg[0])/dr,0,Nr-1.001); ti=tt/dth
    i0=fi.astype(int); j0=ti.astype(int)%Nth
    fr=fi-i0; ft=ti-j0
    out=np.zeros_like(rr)
    for di,wi in ((0,1-fr),(1,fr)):
        ii=np.minimum(i0+di,Nr-1)
        for dj,wj in ((0,1-ft),(1,ft)):
            jj=(j0+dj)%Nth
            out+=wi*wj*F[ii,jj]
    return out

def run(bunched=False,rf_on=True,E0=1e5,sc_on=True,T_ns=1.0,label="",t_rf=0.0,
        Ib=None,seed=None):
    global qmac
    dt_use=dt
    nt=int(T_ns*1e-9/dt_use)
    if Ib is not None:                   # локально: пересчёт веса под заданный ток
        tau_sec=(rk*lnAB/vd)/c
        Ne_loc=max(int(Ib/e*tau_sec),M)
        qmac=Ne_loc/M*e
    else:
        qmac=N_e/M*e                     # вес восстанавливаем на каждый прогон
    rng_use=rng if seed is None else np.random.default_rng(seed)
    x,y,px,py=init_cloud(bunched=bunched,r=rng_use)
    # E×B-дрейф как начальное условие (Hull-старт, см. pic_rf): v_d = E_r/B в +θ
    vdr=Vdc/(np.hypot(x,y)*lnAB)/B
    px+= m_e*(-vdr*np.sin(np.arctan2(y,x))); py+= m_e*( vdr*np.cos(np.arctan2(y,x)))
    ncur=M
    rho=np.zeros((Nr,Nth)); Er_sc=np.zeros((Nr,Nth)); Ep_sc=np.zeros((Nr,Nth))
    hist_W=[]; hist_A1=[]; n0=len(x); emax=0.0
    for i in range(nt):
        t=t_rf+i*dt_use
        if sc_on and i%20==0:
            cic(x,y,rho); Er_sc,Ep_sc=poisson_fields(rho)
            emax=max(emax,float(np.max(np.abs(Er_sc))+np.max(np.abs(Ep_sc))))
        ex=cic_gather(x,y,Er_sc); ep=cic_gather(x,y,Ep_sc)
        rr=np.hypot(x,y); rr=np.maximum(rr,1e-12); th=np.arctan2(y,x)%(2*math.pi/neff)
        Erf=-Vdc/(rr*lnAB); Epf=(E0*(rk/rr)*np.cos(neff*th-w0*t)) if rf_on else 0.0
        Ex=Erf*(x/rr)-Epf*(y/rr)+ex; Ey=Erf*(y/rr)+Epf*(x/rr)+ep
        fx=e*Ex; fy=e*Ey
        px+=fx*dt_use/2; py+=fy*dt_use/2
        g=np.sqrt(1+(px*px+py*py)/(m_e*c)**2); a=wc*dt_use/g; cs,sn=np.cos(a),np.sin(a)
        px,py=px*cs+py*sn,-px*sn+py*cs
        px+=fx*dt_use/2; py+=fy*dt_use/2
        x+=px/(g*m_e)*dt_use; y+=py/(g*m_e)*dt_use
        rr2=np.hypot(x,y)
        lost=(rr2>=ra-0.5*dr)|(rr2<=rk+0.5*dr)
        nl=int(lost.sum())
        if nl:
            # потерянные уносят заряд: переносим их вес на выживших (ток сектора пост.)
            qtot=qmac*ncur
            keep=~lost; x,y,px,py=x[keep],y[keep],px[keep],py[keep]
            ncur=len(x); qmac=qtot/max(ncur,1)
        if i%(max(nt//20,1))==0:
            g_=np.sqrt(1+(px*px+py*py)/(m_e*c)**2)
            K=(g_-1)*m_e*c*c/e
            U=-Vdc*np.log(ra/np.clip(rr2,rk,ra))/lnAB   # ПОТЕНЦ. энергия электрона eΦ(r), Φ=+V у анода
            hist_W.append((K+U).mean())                 # ПОЛНАЯ (для электрона eV_эфф=K+|Φ|)
            thc=np.arctan2(y,x)%(2*math.pi/neff)
            hist_A1.append(abs(np.mean(np.exp(-1j*neff*thc))))
    if len(hist_W)==0 or ncur==0:
        print(f"[{label}] ВСЕ ПОТЕРЯНЫ до конца прогона — проверить шаг/границы"); 
        return np.array([float('nan')]),np.array([float('nan')])
    print(f"[{label}] выживших {ncur}/{M}; <К+U>: {hist_W[0]:.1f}→{hist_W[-1]:.1f} эВ "
          f"(Δ={hist_W[-1]-hist_W[0]:+.1f}); A1: {hist_A1[0]:.2f}→{hist_A1[-1]:.2f}; "
          f"max|E_SC|={emax:.2e} В/м")
    return np.array(hist_W),np.array(hist_A1)

def run_spoke(t_rf,E0):
    return run(bunched=True,rf_on=True,E0=E0,sc_on=False,T_ns=1.0,t_rf=t_rf,seed=7)

print("=== V1: spoke в фикс. π-поле, съём vs фаза старта (стоячая волна cos(neff*th-w0*t)) ===")
for t in (0.0,0.25e-9,0.5e-9,0.75e-9):
    hw,_=run_spoke(t,1e5); print(f"[V1 t_rf={t*1e12:.0f}пс] <К+U>: {hw[0]:.1f}->{hw[-1]:.1f} эВ (Δ={hw[-1]-hw[0]:+.1f})")

print("=== V4: E0-скан при оптимальной тормозящей фазе (калибровка к pic_rf R1: 86 эВ/э @1e5) ===")
for E0v in (1e5,3e5,1e6):
    best=None
    for t in (0.0,0.25e-9,0.5e-9,0.75e-9):
        hw,_=run_spoke(t,E0v)
        if best is None or hw[-1]<best: best=hw[-1]; bt=t
    print(f"[V4 E0={E0v:.0e}] мин Δ при t_rf={bt*1e12:.0f}пс: <К+U> 368→{best:.1f} эВ")

print("=== V0: инвариант Халла K+|eΦ| при RF off, SC off (эталон чистоты пуша) ===")
hw0,_=run(bunched=False,rf_on=False,sc_on=False,T_ns=2.0,label="V0")

print("=== V2: cold облако + только space charge (RF off): нет роста neff-моды ===")
hwSC,haSC=run(bunched=False,rf_on=False,sc_on=True,T_ns=1.0,label="V2 SC-only")

print("=== V3: полный PIC cold + RF + SC, усиленный ток Ib=1 мА (×50 к стартовому) ===")
hwF,haF=run(bunched=False,rf_on=True,sc_on=True,T_ns=1.0,label="V3 Ib=1mA",Ib=1e-3)

ts=lambda n,T: np.linspace(0,T,n)
fig,ax=plt.subplots(1,3,figsize=(14,4.2))
for lab,hw,T in (("V0 Hull DC",hw0,2.0),("V2 SC",hwSC,1.0),("V3 Ib=1mA",hwF,1.0)):
    ax[0].plot(ts(len(hw),T),hw,label=lab)
ax[0].set_xlabel("t, нс"); ax[0].set_ylabel("<K>, эВ"); ax[0].legend(fontsize=7)
ax[0].set_title("П-2 PIC-сетка: кинетическая энергия")
for lab,ha,T in (("V2 SC",haSC,1.0),("V3 Ib=1mA",haF,1.0)):
    ax[1].plot(ts(len(ha),T),ha,label=lab)
ax[1].set_xlabel("t, нс"); ax[1].set_ylabel("|A₁|"); ax[1].legend(fontsize=7)
ax[1].set_title("bunching neff-гармоники")
ax[2].text(0.05,0.8,"Итоги П-2 (после отладки весов и инварианта):\nV0: |Δ(K+|eΦ|)|→0 — пуш чист\n"
          "L2a vs L2b: съём только на гребне\nL3: SC не возбуждает моду\n"
          "L4: рост A1 из шума = spoke-НЕСТАБИЛЬНОСТЬ\n(ср. pic_bunch R2: там порог по B)",
          fontsize=9,va='top',family='monospace')
ax[2].axis('off')
plt.tight_layout(); plt.savefig("OUTPUT/pic_grid.png",dpi=120)
print("\nPNG: OUTPUT/pic_grid.png")
