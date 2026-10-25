//! Ядро значения на Rust — этап 4 плана.
//!
//! Тот же алгоритм, что в `ccore.c`: значение ± ∏ p^(num/den), члены
//! отсортированы по простому, общий знаменатель показателей один на всё
//! значение, ёмкость фиксирована.
//!
//! Смысл этапа — не скорость (она та же), а ОТСУТСТВИЕ НЕОПРЕДЕЛЁННОГО
//! ПОВЕДЕНИЯ в путях с широкими целыми и переполнением профиля,
//! проверяемое компилятором. Потому арифметика вынесена в модуль
//! `math`, объявленный `#![forbid(unsafe_code)]`: там небезопасного кода
//! нет не по обещанию, а по запрету компилятора. `unsafe` остаётся
//! только в тонком слое связи, где разыменовываются указатели ctypes, и
//! он ничего не считает.
//!
//! Переполнение: все промежуточные величины в i64, запись в i32 — через
//! явную проверку диапазона. В отладке Rust паникует на переполнении, в
//! выпуске заворачивает; здесь не полагаемся ни на то, ни на другое —
//! проверка стоит до записи и в обоих режимах одна.

#![allow(non_camel_case_types)]

pub const PHI_CAP: usize = 16;

pub const PHI_OK: i32 = 0;
pub const PHI_E_SUPPORT: i32 = 1;
pub const PHI_E_DEN: i32 = 2;
pub const PHI_E_EXP: i32 = 3;
pub const PHI_E_EVEN_ROOT: i32 = 4;
pub const PHI_E_INPUT: i32 = 5;

#[repr(C)]
#[derive(Clone, Copy)]
pub struct PhiT {
    pub sign: i32,
    pub den: i32,
    pub n: i32,
    pub prime: [i32; PHI_CAP],
    pub num: [i32; PHI_CAP],
}

impl PhiT {
    pub fn zeroed() -> Self {
        PhiT { sign: 1, den: 1, n: 0, prime: [0; PHI_CAP], num: [0; PHI_CAP] }
    }
}

/// Вся арифметика. Небезопасного кода здесь нет по запрету компилятора.
mod math {
    #![forbid(unsafe_code)]

    use super::*;

    pub fn gcd64(mut a: i64, mut b: i64) -> i64 {
        if a < 0 { a = -a; }
        if b < 0 { b = -b; }
        while b != 0 {
            let t = a % b;
            a = b;
            b = t;
        }
        a
    }

    pub fn fits32(v: i64) -> bool {
        v >= i32::MIN as i64 && v <= i32::MAX as i64
    }

    /// Канон на месте: выброс нулей, сокращение знаменателя с НОД.
    pub fn canon(v: &mut PhiT) {
        let mut k = 0usize;
        for i in 0..(v.n as usize) {
            if v.num[i] != 0 {
                v.prime[k] = v.prime[i];
                v.num[k] = v.num[i];
                k += 1;
            }
        }
        v.n = k as i32;
        let mut g = v.den as i64;
        for i in 0..k {
            g = gcd64(g, v.num[i] as i64);
        }
        if g > 1 {
            v.den = (v.den as i64 / g) as i32;
            for i in 0..k {
                v.num[i] = (v.num[i] as i64 / g) as i32;
            }
        }
    }

    pub fn set(out: &mut PhiT, sign: i32, den: i32, primes: &[i32], nums: &[i32]) -> i32 {
        if sign != 1 && sign != -1 || den < 1 {
            return PHI_E_INPUT;
        }
        if primes.len() != nums.len() {
            return PHI_E_INPUT;
        }
        if primes.len() > PHI_CAP {
            return PHI_E_SUPPORT;
        }
        *out = PhiT::zeroed();
        out.sign = sign;
        out.den = den;
        let mut k = 0usize;
        for i in 0..primes.len() {
            if primes[i] < 2 {
                return PHI_E_INPUT;
            }
            if i > 0 && primes[i] < primes[i - 1] {
                return PHI_E_INPUT;
            }
            if k > 0 && out.prime[k - 1] == primes[i] {
                let s = out.num[k - 1] as i64 + nums[i] as i64;
                if !fits32(s) {
                    return PHI_E_EXP;
                }
                out.num[k - 1] = s as i32;
            } else {
                out.prime[k] = primes[i];
                out.num[k] = nums[i];
                k += 1;
            }
        }
        out.n = k as i32;
        canon(out);
        PHI_OK
    }

    /// Слияние двух отсортированных списков; sgn = +1 умножение, −1 деление.
    pub fn merge(a: &PhiT, b: &PhiT, sgn: i64, out: &mut PhiT) -> i32 {
        let g = gcd64(a.den as i64, b.den as i64);
        let d = (a.den as i64 / g) * b.den as i64;
        if !fits32(d) {
            return PHI_E_DEN;
        }
        let (ka, kb) = (d / a.den as i64, d / b.den as i64);
        let mut r = PhiT::zeroed();
        r.sign = a.sign * b.sign;
        r.den = d as i32;
        let (mut i, mut j, mut k) = (0usize, 0usize, 0usize);
        let (na, nb) = (a.n as usize, b.n as usize);
        while i < na || j < nb {
            if k >= PHI_CAP {
                return PHI_E_SUPPORT;
            }
            let (p, e): (i32, i64) = if j >= nb || (i < na && a.prime[i] < b.prime[j]) {
                let v = (a.prime[i], a.num[i] as i64 * ka);
                i += 1;
                v
            } else if i >= na || b.prime[j] < a.prime[i] {
                let v = (b.prime[j], sgn * (b.num[j] as i64 * kb));
                j += 1;
                v
            } else {
                let v = (a.prime[i], a.num[i] as i64 * ka + sgn * (b.num[j] as i64 * kb));
                i += 1;
                j += 1;
                v
            };
            if !fits32(e) {
                return PHI_E_EXP;
            }
            if e == 0 {
                continue;
            }
            r.prime[k] = p;
            r.num[k] = e as i32;
            k += 1;
        }
        r.n = k as i32;
        canon(&mut r);
        *out = r;
        PHI_OK
    }

    pub fn root(a: &PhiT, k: i32, out: &mut PhiT) -> i32 {
        if k < 1 {
            return PHI_E_INPUT;
        }
        if a.sign < 0 && k % 2 == 0 {
            return PHI_E_EVEN_ROOT;
        }
        let d = a.den as i64 * k as i64;
        if !fits32(d) {
            return PHI_E_DEN;
        }
        let mut r = *a;
        r.den = d as i32;
        canon(&mut r);
        *out = r;
        PHI_OK
    }

    pub fn power(a: &PhiT, k: i32, out: &mut PhiT) -> i32 {
        let mut r = PhiT::zeroed();
        r.sign = if a.sign > 0 || k % 2 == 0 { 1 } else { -1 };
        r.den = a.den;
        r.n = a.n;
        for i in 0..(a.n as usize) {
            let e = a.num[i] as i64 * k as i64;
            if !fits32(e) {
                return PHI_E_EXP;
            }
            r.prime[i] = a.prime[i];
            r.num[i] = e as i32;
        }
        canon(&mut r);
        *out = r;
        PHI_OK
    }

    pub fn inverse(a: &PhiT, out: &mut PhiT) -> i32 {
        let mut r = *a;
        for i in 0..(r.n as usize) {
            r.num[i] = -r.num[i];
        }
        *out = r;
        PHI_OK
    }

    pub fn eq(a: &PhiT, b: &PhiT) -> i32 {
        if a.sign != b.sign || a.den != b.den || a.n != b.n {
            return 0;
        }
        for i in 0..(a.n as usize) {
            if a.prime[i] != b.prime[i] || a.num[i] != b.num[i] {
                return 0;
            }
        }
        1
    }

    /// Партия слияний: n независимых пар за один переход границы.
    /// Коды отдаются поэлементно — партия не обрывается на первом событии.
    pub fn merge_pairs(a: &[PhiT], b: &[PhiT], out: &mut [PhiT], codes: &mut [i32],
                       sgn: i64) -> i32 {
        if a.len() != b.len() || a.len() != out.len() || a.len() != codes.len() {
            return PHI_E_INPUT;
        }
        if sgn != 1 && sgn != -1 {
            return PHI_E_INPUT;
        }
        for i in 0..a.len() {
            codes[i] = merge(&a[i], &b[i], sgn, &mut out[i]);
            if codes[i] != PHI_OK {
                out[i] = PhiT::zeroed();
                out[i].n = 0;
            }
        }
        PHI_OK
    }

    /// Цепочка зависимых действий: acc = (acc·b)/b, reps раз.
    pub fn bench_chain(a: &PhiT, b: &PhiT, out: &mut PhiT, reps: i32) -> i32 {
        let mut acc = *a;
        let mut tmp = PhiT::zeroed();
        for _ in 0..reps {
            let code = merge(&acc, b, 1, &mut tmp);
            if code != PHI_OK {
                return code;
            }
            let code = merge(&tmp, b, -1, &mut acc);
            if code != PHI_OK {
                return code;
            }
        }
        *out = acc;
        PHI_OK
    }
}

// ── тонкий слой связи: только здесь есть unsafe, и он ничего не считает ──

#[no_mangle]
pub extern "C" fn phir_capacity() -> i32 {
    PHI_CAP as i32
}

#[no_mangle]
pub extern "C" fn phir_sizeof() -> i32 {
    core::mem::size_of::<PhiT>() as i32
}

#[no_mangle]
pub unsafe extern "C" fn phir_set(out: *mut PhiT, sign: i32, den: i32,
                                  primes: *const i32, nums: *const i32, n: i32) -> i32 {
    if out.is_null() || n < 0 {
        return PHI_E_INPUT;
    }
    let n = n as usize;
    if n > PHI_CAP {
        return PHI_E_SUPPORT;
    }
    let (p, q) = if n == 0 {
        (&[][..], &[][..])
    } else {
        if primes.is_null() || nums.is_null() {
            return PHI_E_INPUT;
        }
        (core::slice::from_raw_parts(primes, n), core::slice::from_raw_parts(nums, n))
    };
    math::set(&mut *out, sign, den, p, q)
}

macro_rules! binary {
    ($name:ident, $sgn:expr) => {
        #[no_mangle]
        pub unsafe extern "C" fn $name(a: *const PhiT, b: *const PhiT, out: *mut PhiT) -> i32 {
            if a.is_null() || b.is_null() || out.is_null() {
                return PHI_E_INPUT;
            }
            math::merge(&*a, &*b, $sgn, &mut *out)
        }
    };
}

binary!(phir_mul, 1);
binary!(phir_div, -1);

#[no_mangle]
pub unsafe extern "C" fn phir_root(a: *const PhiT, k: i32, out: *mut PhiT) -> i32 {
    if a.is_null() || out.is_null() {
        return PHI_E_INPUT;
    }
    math::root(&*a, k, &mut *out)
}

#[no_mangle]
pub unsafe extern "C" fn phir_pow(a: *const PhiT, k: i32, out: *mut PhiT) -> i32 {
    if a.is_null() || out.is_null() {
        return PHI_E_INPUT;
    }
    math::power(&*a, k, &mut *out)
}

#[no_mangle]
pub unsafe extern "C" fn phir_inverse(a: *const PhiT, out: *mut PhiT) -> i32 {
    if a.is_null() || out.is_null() {
        return PHI_E_INPUT;
    }
    math::inverse(&*a, &mut *out)
}

#[no_mangle]
pub unsafe extern "C" fn phir_eq(a: *const PhiT, b: *const PhiT) -> i32 {
    if a.is_null() || b.is_null() {
        return 0;
    }
    math::eq(&*a, &*b)
}

#[no_mangle]
pub unsafe extern "C" fn phir_merge_pairs(a: *const PhiT, b: *const PhiT,
                                          out: *mut PhiT, codes: *mut i32,
                                          n: i32, sgn: i32) -> i32 {
    if a.is_null() || b.is_null() || out.is_null() || codes.is_null() || n < 0 {
        return PHI_E_INPUT;
    }
    let n = n as usize;
    math::merge_pairs(core::slice::from_raw_parts(a, n),
                      core::slice::from_raw_parts(b, n),
                      core::slice::from_raw_parts_mut(out, n),
                      core::slice::from_raw_parts_mut(codes, n),
                      sgn as i64)
}

#[no_mangle]
pub unsafe extern "C" fn phir_bench_chain(a: *const PhiT, b: *const PhiT,
                                          out: *mut PhiT, reps: i32) -> i32 {
    if a.is_null() || b.is_null() || out.is_null() {
        return PHI_E_INPUT;
    }
    math::bench_chain(&*a, &*b, &mut *out, reps)
}

// ── Колёса и наследственная форма: ОБРАТНОЕ преобразование на ядре Rust ──────
// Слово автора 24.09: обратное преобразование (машинное целое → сумма колёс,
// показатели — рекурсивно) на ядре. Счёт — в модуле без unsafe; тонкий слой
// связи ниже только переводит указатели в срезы. Коды как у C: ≥ 0 — число
// колёс/членов; −1 негодный вход; −2 буфер мал; −3 переполнение i64.
mod wheels_math {
    #![forbid(unsafe_code)]

    pub fn wheels(mut n: i64, base: i32, m: &mut [i32], e: &mut [i32]) -> i32 {
        if base < 2 || n < 0 {
            return -1;
        }
        let b = base as i64;
        let mut k = 0usize;
        let mut pos: i32 = 0;
        while n != 0 {
            let d = n % b;
            n /= b;
            if d != 0 {
                if k >= m.len() || k >= e.len() {
                    return -2;
                }
                m[k] = d as i32;
                e[k] = pos;
                k += 1;
            }
            pos += 1;
        }
        k as i32
    }

    pub fn wheels_batch(ns: &[i64], base: i32, m: &mut [i32], e: &mut [i32], cnt: &mut [i32]) -> i64 {
        if base < 2 {
            return -1;
        }
        let mut off = 0usize;
        for (i, &n) in ns.iter().enumerate() {
            let k = wheels(n, base, &mut m[off..], &mut e[off..]);
            if k < 0 {
                return k as i64;
            }
            cnt[i] = k;
            off += k as usize;
        }
        off as i64
    }

    fn pow_checked(b: i64, e: i64) -> Option<i64> {
        let mut r: i64 = 1;
        let mut k = e;
        while k > 0 {
            r = r.checked_mul(b)?;
            k -= 1;
        }
        Some(r)
    }

    /// (значение дерева в newbase, глубина, колёс всех уровней) — как у
    /// phicore.wheels.Hereditary; None при переполнении i64.
    pub fn hereditary(mut n: i64, base: i32, newbase: i32) -> Option<(i64, i32, i64)> {
        let b = base as i64;
        let nb = newbase as i64;
        let mut total: i64 = 0;
        let mut wc: i64 = 0;
        let mut dmax: i32 = 0;
        let mut any = false;
        let mut pos: i32 = 0;
        while n != 0 {
            let d = n % b;
            n /= b;
            if d != 0 {
                any = true;
                wc += 1;
                let mut e: i64 = 0;
                let mut sd: i32 = 0;
                if pos >= 1 {
                    let (ev, sdd, sw) = hereditary(pos as i64, base, newbase)?;
                    e = ev;
                    sd = sdd;
                    wc += sw;
                }
                let p = pow_checked(nb, e)?;
                let term = d.checked_mul(p)?;
                total = total.checked_add(term)?;
                if sd > dmax {
                    dmax = sd;
                }
            }
            pos += 1;
        }
        Some((total, if any { dmax + 1 } else { 0 }, wc))
    }

    pub fn goodstein(n: i64, out: &mut [i64]) -> i32 {
        if n < 0 || out.is_empty() {
            return -1;
        }
        out[0] = n;
        let mut k = 1usize;
        let mut base: i32 = 2;
        while k < out.len() && out[k - 1] != 0 {
            match hereditary(out[k - 1], base, base + 1) {
                Some((v, _, _)) => out[k] = v - 1,
                None => return -3,
            }
            k += 1;
            base += 1;
        }
        k as i32
    }
}

#[no_mangle]
pub unsafe extern "C" fn phir_wheels(n: i64, base: i32, m: *mut i32, e: *mut i32, cap: i32) -> i32 {
    if m.is_null() || e.is_null() || cap < 0 {
        return -1;
    }
    let (ms, es) = (core::slice::from_raw_parts_mut(m, cap as usize),
                    core::slice::from_raw_parts_mut(e, cap as usize));
    wheels_math::wheels(n, base, ms, es)
}

#[no_mangle]
pub unsafe extern "C" fn phir_wheels_batch(n: *const i64, count: i32, base: i32, m: *mut i32,
                                           e: *mut i32, cnt: *mut i32, cap: i64) -> i64 {
    if n.is_null() || m.is_null() || e.is_null() || cnt.is_null() || count < 0 || cap < 0 {
        return -1;
    }
    let ns = core::slice::from_raw_parts(n, count as usize);
    let (ms, es) = (core::slice::from_raw_parts_mut(m, cap as usize),
                    core::slice::from_raw_parts_mut(e, cap as usize));
    let cs = core::slice::from_raw_parts_mut(cnt, count as usize);
    wheels_math::wheels_batch(ns, base, ms, es, cs)
}

#[no_mangle]
pub unsafe extern "C" fn phir_hereditary(n: i64, base: i32, newbase: i32, value: *mut i64,
                                         depth: *mut i32, wheels: *mut i64) -> i32 {
    if value.is_null() || depth.is_null() || wheels.is_null() || base < 2 || newbase < 2 || n < 0 {
        return -1;
    }
    match wheels_math::hereditary(n, base, newbase) {
        Some((v, d, w)) => {
            *value = v;
            *depth = d;
            *wheels = w;
            0
        }
        None => -3,
    }
}

#[no_mangle]
pub unsafe extern "C" fn phir_goodstein(n: i64, steps: i32, out: *mut i64) -> i32 {
    if out.is_null() || steps < 0 {
        return -1;
    }
    let os = core::slice::from_raw_parts_mut(out, steps as usize + 1);
    wheels_math::goodstein(n, os)
}
