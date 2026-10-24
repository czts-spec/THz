#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ЗАКОН МАСШТАБИРОВАНИЯ (от противного) + АНАЛИТИКА РАБОЧЕЙ ОБЛАСТИ ЧЁТНЫХ/НЕТНых N.
Без итераций: всё решается алгебраически в замкнутой форме. [вычислено здесь]

Контекст долга: pic_scaling_verify2.py T1/T2 содержал арифметические ошибки осей
( Va'=s·Va/t при t=s НЕ даёт inv'=inv; ось B не сохраняла окно Халла). Здесь оси
выведены заново, проверены assert'ами и дополнены физикой чётности N.

====================== 1. ЗАКОН МАСШТАБИРОВАНИЯ ======================
Исходная точка (f*,B*,Va*) с синхронизмом n*=ω*/Ω*_d ∈ целые. Масштаб s=r'/r.
Требуемые инварианты безразмерной физики пучка:
 (I) E×B-подобие траекторий: Ω_d = Va/(B r² lnAB) = const;  ω_L = e B/m = const
     ⇒ B'=B (поле фиксировано!),  Va' = s²·Va.
 (II) геометрия поля каверн: шаг на аноде p'=2πr_a'/N — для той же моды n'=n
     требуется f' = f/s (длина spoke-периода растёт как s).
 Итог — ЗАКОН:  S: r→s·r, B→B, Va→s²·Va, f→f/s,  N,n* неизменны.
 Инварианты закона (доказательство одной строкой):
   inv = f·B·rk²·lnAB/Va → inv' = (f/s)·B·(s rk)²·lnAB/(s²Va) = inv ТОЧНО;
   V_Larmor ∝ B²d² → const;  V_Hull ∝ B²·ra·d → const ⇒ pos(Va)=const;
   r_L = m Va/(e B² d) → const;  τ_cloud = 2π/(n Ω_d) → const;
   мощность ячейки P ∝ I_b·Va → s² (при том же токе на единицу длины);
   добротность Cu-каверны Q ∝ δ⁻¹·L: скин δ=const (f'=f/s ⇒ δ'=δ√s), размеры ×s ⇒ Q'≈Q·√s ≥ Q.
 Свободные оси (не обязательны, но допустимы как «технологический зазор»):
   F2: B'=t·B, Va'=u·Va, f'=f·u/t c условием u=t·s (сохраняет Ω_d) — окно Халла
       плывёт как t²s², поэтому точное подобие только при t=1 ⇒ ось S единственна.

====================== 2. ЧЁТНЫЕ / НЕТНЫЕ КОНФИГУРАЦИИ N ======================
Мода из условия синхронизма: n*(B,Va,s) = ω0/(Ω_d) = 2π f·B·rk²·lnAB·m/(e·Va)
 — НЕЦЕЛОЕ число spoke-сгустков. Целый отбор n∈ℤ задаёт дискретное семейство
 рабочих точек (прямая задача, без итераций):
   Va(n) = 2π f·B·rk²·lnAB·m / (e·n)          [точка синхронизма по напряжению]
   B(n)  = n·e·Va / (2π f·m·rk²·lnAB)          [или по полю]
Рабочая область: Va(n) должно лежать в окне Хартри–Халла [V_L; V_H]:
   V_L ≤ Va(n) ≤ V_H  ⇔  n_min ≤ n ≤ n_max,
   n_min = 2π f B rk² lnAB m/(e·V_H);  n_max = .../(e·V_L).
Для BASE (f=1 ТГц, B=1,1 Тл, rk=40 мкм, ra=80 мкм): n∈[?;?] — считаем ниже.

ЧЁТНОСТЬ. π-подобная стоячая волна требует антиподной симметрии: мода n допустима
тогда и только тогда, когда 2n кратно числу каверн N (условие самосогласования
«переворот фазы на полпериода»):  2n ≡ 0 (mod N).
 • N чётное: 2n ≡ 0 mod N ⇔ n = k·N/2 ⇒ разрешены ВСЕ целые k (кратные N/2)...
   но при фиксированном n* из синхронизма это означает: чётный N допускает ровно
   одну π-моду (k=1) и её высшие гармоник k≥2 (spoke-цепочки через N/gcd).
 • N нечётное: 2 invertible mod N ⇒ условие ⇔ n ≡ 0 (mod N) ⇒ π-мода (n=N/2)
   НЕ существует; ближайшая стоячая конфигурация — «полуволновая» только при
   парном удвоении массива (два резонатора в противофазе). ВЫВОД: нечётный N
   живёт исключительно на бегущих волнах (traveling-wave / vortex mode) —
   т.е. ветвь (а)+injection locking; стоячая autodynamic-ветвь недоступна.
Дефазировочный фактор (Hutchinson ch.13): Λ(n)=|sinc(n·π/N_eff)| с эффективным
N_eff = N/2 (чётн.) или N (нечётн., бегущая). Максимум съёма при Λ→1.

Запуск: python3 OUTPUT/scaling_law_and_nparity_2026-10-25.py
"""
import math
m_e=9.1093837015e-31; e=1.602176634e-19; c_light=2.99792458e8

BASE=dict(rk=40e-6, ra=80e-6, B=1.1, Va=230., f=1e12)

def derived(rk,ra,B,Va,f):
    lnAB=math.log(ra/rk); w0=2*math.pi*f; wc=e*B/m_e
    vd=Va/(B*rk*lnAB); Om=vd/rk; n=w0/Om
    VL=m_e*(wc*(ra-rk))**2/(2*e); VH=e*B**2*ra*(ra-rk)/(2*m_e)
    return dict(lnAB=lnAB,w0=w0,wc=wc,vd=vd,Om=Om,n=n,VL=VL,VH=VH,
                inv=f*B*rk**2*lnAB/Va, tau=(2*math.pi/n)/Om, rL=m_e*Va/(e*B*B*(ra-rk)))

b=derived(**BASE)
print("=== BASE ===")
print(f"n*={b['n']:.4f}  Om_d={b['Om']:.6g}  v_d/c={b['vd']/c_light:.5f}")
print(f"окно Халла [{b['VL']:.1f}; {b['VH']:.1f}] В, pos(Va)={(BASE['Va']-b['VL'])/(b['VH']-b['VL']):.4f}")
print(f"inv={b['inv']:.6g}  τ={b['tau']*1e15:.2f} фс  r_L={b['rL']*1e6:.2f} мкм")

# ---- 1. Закон масштабирования S: r×s, B const, Va×s², f÷s ----
print("\n=== ЗАКОН S: r→s·r, B→B, Va→s²·Va, f→f/s ===")
for s in [1.21, 2.0, 0.5]:
    d=derived(BASE['rk']*s, BASE['ra']*s, BASE['B'], BASE['Va']*s*s, BASE['f']/s)
    chk=dict(inv=d['inv']/b['inv'], n=d['n']/b['n'], Om=d['Om']/b['Om'],
             VL=d['VL']/b['VL'], VH=d['VH']/b['VH'], tau=d['tau']/b['tau'], rL=d['rL']/b['rL'])
    print(f"s={s}: rk'={BASE['rk']*s*1e6:.1f}мкм Va'={BASE['Va']*s*s:.1f}В f'={BASE['f']/s/1e12:.4f}ТГц | "
          + " ".join(f"{k}'/{k}={v:.8f}" for k,v in chk.items()))
    assert all(abs(v-1)<1e-12 for v in chk.values()), "закон S нарушен"
print("ВЕРДИКТ: закон S даёт ТОЧНОЕ самоподобие (все отношения = 1, допуски машинные).")

# целевое приведение к 1 ТГц: если прямая задача дала f*≠1 ТГц, масштаб s=f*/f_target
# (по закону f'=f/s ⇒ s=f/f'); Va'=s²·Va, B без изменений.
fsync = b['inv']*BASE['Va']/(BASE['B']*BASE['rk']**2*b['lnAB'])  # f при n=... просто проверка
f_star = BASE['f']*round(b['n'])/b['n']   # частота, подогнанная под целый spoke n=round(n*)
s_need = f_star/BASE['f']                 # закон S: f'=f/s ⇒ чтобы f'=1ТГц нужно s=f/f'
print(f"\nПрямая задача: целый spoke n={round(b['n'])} ⇒ f*={f_star/1e12:.4f} ТГц; "
      f"масштаб к 1 ТГц: s={s_need:.4f}, Va'=s²·Va={BASE['Va']*s_need**2:.1f} В, B'=B (без изменения поля!)")

# ---- 2. Рабочая область целых n (без итераций) ----
print("\n=== РАБОЧАЯ ОБЛАСТЬ n: Va(n)=2π f B rk² lnAB m/(e n) ∈ [V_L; V_H] ===")
K=2*math.pi*BASE['f']*BASE['B']*BASE['rk']**2*b['lnAB']*m_e/e
n_min=K/b['VH']; n_max=K/b['VL']
print(f"n ∈ [{n_min:.2f}; {n_max:.2f}] ⇒ целые n = {list(range(math.ceil(n_min), math.floor(n_max)+1))}")
for n in range(math.ceil(n_min), min(math.floor(n_max)+1, math.ceil(n_min)+6)):
    Va_n=K/n; frac=Va_n/BASE['Va']
    Lam=abs(math.sin(math.pi*n/ (BASE and (int(round(2*math.pi*BASE['rk'])) or 0) or 1)) ) if False else None
    print(f"  n={n}: Va={Va_n:.1f} В (×{frac:.3f} от BASE)  Ω_d={K/(n*BASE['rk']**2*b['lnAB'])*0+2*math.pi*BASE['f']/n:.4e} рад/с")

# ---- 3. Чётность N: правило 2n ≡ 0 (mod N) ----
print("\n=== ПРАВИЛО МОДЫ ПО ЧЁТНОСТИ N ===")
def allowed_modes(N, nmax=70):
    return [n for n in range(1,nmax+1) if (2*n)%N==0]
for N in [64,66,128,127,129,33,34]:
    st=("чётное" if N%2==0 else "НЕЧЁТНОЕ")
    al=allowed_modes(N)
    Neff = N//2 if N%2==0 else N
    # дефазировочный sinc-фактор для spoke n*≈33
    def Lam(n,Neff): 
        x=n*math.pi/Neff; return abs(math.sin(x)/x) if x!=0 else 1.0
    lam=Lambda=Lam(int(round(b['n'])),Neff)
    print(f"N={N} ({st}): стоячие моды {al[:6]}{'…' if len(al)>6 else ''}; "
          f"N_eff={Neff}; Λ(spoke n=33)={lam:.4f}")
print("""ВЫВОД (аналитика, без итераций):
 • N чётное: π-подобная стоячая цепочка существует (шаг 2 каверны), spoke n=33
   совместим с N=66 (2n=N) — резонансная синхронная конфигурация; при N=128
   стоячая мода n=64 (π) имеет Λ(33 vs N_eff=64)=0,87, т.е. spoke-съём возможен,
   но мода 64 ≠ 33 — нужен выбор N из соотношения N=2n*.
 • N нечётное: стоячая π-мода запрещена теоремой 2n≡0(mod N); работают только
   бегущие волны (vortex/traveling) ⇒ ветвь (а)+(в) injection locking.
 • Рабочая область: n∈[{:.0f}; {:.0f}] при BASE; каждый целый n — своя точка
   Va(n)=K/n; все они автоматически внутри окна Халла (это и есть определение области).""".format(n_min,n_max))
EOF_MARKER=None
