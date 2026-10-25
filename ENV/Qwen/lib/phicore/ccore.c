/* Ядро значения на C — этап 3 плана.
 *
 * Значение: ± prod p^(num_p / den), члены отсортированы по простому,
 * общий знаменатель показателей один на всё значение. Ёмкость
 * фиксирована (PHI_CAP): профиль ограниченной ширины — условие
 * синтезируемости (этап 5), а не оптимизация.
 *
 * Правила, которые здесь важнее скорости:
 *  - в горячем пути НЕТ выделения памяти: все структуры на стеке
 *    вызывающего;
 *  - переполнение НЕ заворачивается молча: промежуточные величины
 *    считаются в int64, и выход за int32 есть код возврата, а не
 *    испорченное значение;
 *  - канон один и тот же, что у эталона на Python: сортировка по
 *    простому, нулевые показатели выброшены, знаменатель сокращён с
 *    НОД всех числителей. Отсюда равенство есть побайтовое совпадение
 *    полей, а не вызов упрощателя.
 *
 * Собирается без заголовков Python: связь — ctypes. Так ядро остаётся
 * чистым алгоритмом и годится этапам 4 (Rust) и 5 (FPGA).
 */
#include <stdint.h>
#include <string.h>

#define PHI_CAP 16

#define PHI_OK            0
#define PHI_E_SUPPORT     1  /* носитель шире ёмкости */
#define PHI_E_DEN         2  /* общий знаменатель за int32 */
#define PHI_E_EXP         3  /* показатель за int32 */
#define PHI_E_EVEN_ROOT   4  /* чётный корень из отрицательного */
#define PHI_E_INPUT       5  /* негодный вход */

typedef struct {
    int32_t sign;              /* +1 или -1 */
    int32_t den;               /* общий знаменатель показателей, >= 1 */
    int32_t n;                 /* число членов */
    int32_t prime[PHI_CAP];
    int32_t num[PHI_CAP];
} phi_t;

static int64_t gcd64(int64_t a, int64_t b)
{
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) { int64_t t = a % b; a = b; b = t; }
    return a;
}

static int fits32(int64_t v)
{
    return v >= INT32_MIN && v <= INT32_MAX;
}

/* Канон на месте: выброс нулей, сокращение знаменателя. Порядок членов
 * входом уже задан (все построители держат его отсортированным). */
static int canon(phi_t *v)
{
    int32_t k = 0;
    for (int32_t i = 0; i < v->n; i++) {
        if (v->num[i] != 0) {
            v->prime[k] = v->prime[i];
            v->num[k] = v->num[i];
            k++;
        }
    }
    v->n = k;
    int64_t g = v->den;
    for (int32_t i = 0; i < v->n; i++) g = gcd64(g, v->num[i]);
    if (g > 1) {
        v->den = (int32_t)(v->den / g);
        for (int32_t i = 0; i < v->n; i++) v->num[i] = (int32_t)(v->num[i] / g);
    }
    return PHI_OK;
}

int phi_capacity(void) { return PHI_CAP; }
int phi_sizeof(void) { return (int)sizeof(phi_t); }

/* Построение из отсортированного списка пар. Повторы простых
 * складываются; порядок обязан быть неубывающим — иначе PHI_E_INPUT. */
int phi_set(phi_t *out, int32_t sign, int32_t den,
            const int32_t *primes, const int32_t *nums, int32_t n)
{
    if (!out || (sign != 1 && sign != -1) || den < 1) return PHI_E_INPUT;
    if (n < 0 || n > PHI_CAP) return PHI_E_SUPPORT;
    memset(out, 0, sizeof(*out));
    out->sign = sign;
    out->den = den;
    int32_t k = 0;
    for (int32_t i = 0; i < n; i++) {
        if (primes[i] < 2) return PHI_E_INPUT;
        if (i > 0 && primes[i] < primes[i - 1]) return PHI_E_INPUT;
        if (k > 0 && out->prime[k - 1] == primes[i]) {
            int64_t s = (int64_t)out->num[k - 1] + nums[i];
            if (!fits32(s)) return PHI_E_EXP;
            out->num[k - 1] = (int32_t)s;
        } else {
            out->prime[k] = primes[i];
            out->num[k] = nums[i];
            k++;
        }
    }
    out->n = k;
    return canon(out);
}

/* Слияние двух отсортированных списков; sgn = +1 умножение, -1 деление. */
static int merge(const phi_t *a, const phi_t *b, int sgn, phi_t *out)
{
    if (!a || !b || !out) return PHI_E_INPUT;
    int64_t g = gcd64(a->den, b->den);
    int64_t d = (int64_t)a->den / g * b->den;
    if (!fits32(d)) return PHI_E_DEN;
    int64_t ka = d / a->den, kb = d / b->den;
    phi_t r;
    memset(&r, 0, sizeof(r));
    r.sign = a->sign * b->sign;
    r.den = (int32_t)d;
    int32_t i = 0, j = 0, k = 0;
    while (i < a->n || j < b->n) {
        if (k >= PHI_CAP) return PHI_E_SUPPORT;
        int32_t p;
        int64_t e;
        if (j >= b->n || (i < a->n && a->prime[i] < b->prime[j])) {
            p = a->prime[i];
            e = (int64_t)a->num[i] * ka;
            i++;
        } else if (i >= a->n || b->prime[j] < a->prime[i]) {
            p = b->prime[j];
            e = sgn * ((int64_t)b->num[j] * kb);
            j++;
        } else {
            p = a->prime[i];
            e = (int64_t)a->num[i] * ka + sgn * ((int64_t)b->num[j] * kb);
            i++; j++;
        }
        if (!fits32(e)) return PHI_E_EXP;
        if (e == 0) continue;              /* нулевой показатель не хранится */
        r.prime[k] = p;
        r.num[k] = (int32_t)e;
        k++;
    }
    r.n = k;
    canon(&r);
    *out = r;
    return PHI_OK;
}

int phi_mul(const phi_t *a, const phi_t *b, phi_t *out) { return merge(a, b, 1, out); }
int phi_div(const phi_t *a, const phi_t *b, phi_t *out) { return merge(a, b, -1, out); }

/* Корень степени k — деление показателей: знаменатель умножается на k. */
int phi_root(const phi_t *a, int32_t k, phi_t *out)
{
    if (!a || !out || k < 1) return PHI_E_INPUT;
    if (a->sign < 0 && (k % 2) == 0) return PHI_E_EVEN_ROOT;
    int64_t d = (int64_t)a->den * k;
    if (!fits32(d)) return PHI_E_DEN;
    phi_t r = *a;
    r.den = (int32_t)d;
    canon(&r);
    *out = r;
    return PHI_OK;
}

/* Целая степень — умножение показателей. */
int phi_pow(const phi_t *a, int32_t k, phi_t *out)
{
    if (!a || !out) return PHI_E_INPUT;
    phi_t r;
    memset(&r, 0, sizeof(r));
    r.sign = (a->sign > 0 || (k % 2) == 0) ? 1 : -1;
    r.den = a->den;
    r.n = a->n;
    for (int32_t i = 0; i < a->n; i++) {
        int64_t e = (int64_t)a->num[i] * k;
        if (!fits32(e)) return PHI_E_EXP;
        r.prime[i] = a->prime[i];
        r.num[i] = (int32_t)e;
    }
    canon(&r);
    *out = r;
    return PHI_OK;
}

int phi_inverse(const phi_t *a, phi_t *out)
{
    if (!a || !out) return PHI_E_INPUT;
    phi_t r = *a;
    for (int32_t i = 0; i < r.n; i++) r.num[i] = -r.num[i];
    *out = r;
    return PHI_OK;
}

/* Партия действий за один переход границы: n независимых слияний одним
 * вызовом. Замер §445 показал, что граница ctypes стоит 15 мкс, а сама
 * операция — 0,04 мкс: поштучный вызов меряет границу, а не алгоритм.
 * Здесь на партию приходится один переход, и коды возврата отдаются
 * поэлементно — партия НЕ обрывается на первом событии ширины, иначе
 * пропали бы итоги всех последующих пар. */
int phi_merge_pairs(const phi_t *a, const phi_t *b, phi_t *out, int32_t *codes,
                    int32_t n, int32_t sgn)
{
    if (!a || !b || !out || !codes || n < 0) return PHI_E_INPUT;
    if (sgn != 1 && sgn != -1) return PHI_E_INPUT;
    for (int32_t i = 0; i < n; i++) {
        codes[i] = merge(&a[i], &b[i], sgn, &out[i]);
        if (codes[i] != PHI_OK) memset(&out[i], 0, sizeof(phi_t));
    }
    return PHI_OK;
}

/* Цепочка действий внутри C — чтобы цена САМОГО ядра отделялась от цены
 * связи. Просто «повторить одно и то же reps раз» здесь не годится:
 * итерации независимы, и компилятор вправе оставить одну — замер тогда
 * лжёт. Потому каждая итерация ЗАВИСИТ от предыдущей: acc = (acc·b)/b.
 * Значение при этом возвращается к исходному, и это проверяется —
 * замер обязан быть ещё и верным. Действий 2·reps. */
int phi_bench_chain(const phi_t *a, const phi_t *b, phi_t *out, int32_t reps)
{
    phi_t acc = *a, tmp;
    for (int32_t i = 0; i < reps; i++) {
        int code = merge(&acc, b, 1, &tmp);
        if (code != PHI_OK) return code;
        code = merge(&tmp, b, -1, &acc);
        if (code != PHI_OK) return code;
    }
    *out = acc;
    return PHI_OK;
}

/* Равенство — совпадение канона, поле в поле. */
int phi_eq(const phi_t *a, const phi_t *b)
{
    if (!a || !b) return 0;
    if (a->sign != b->sign || a->den != b->den || a->n != b->n) return 0;
    for (int32_t i = 0; i < a->n; i++)
        if (a->prime[i] != b->prime[i] || a->num[i] != b->num[i]) return 0;
    return 1;
}

/* ── Колёса и наследственная форма: ОБРАТНОЕ преобразование на ядре ────────
 *
 * Слово автора 24.09: «обратное преобразование на акселераторах не выполнено,
 * и что самое странное не предполагалось». Здесь оно есть: машинное целое
 * (позиционная сокращённая запись) раскладывается на колёса m·b^e внутри
 * ядра, а не сводится к десятичной форме снаружи. Наследственная форма
 * (Goodstein 1944) — то же разложение, применённое к собственным показателям.
 * Коды: ≥ 0 — число колёс; −1 негодный вход; −2 буфер мал; −3 переполнение
 * int64 при чтении дерева в новой базе. Эталон — phicore/wheels.py, сверка
 * побитовая (proto/test_phicore_hereditary.py).
 */
int phi_wheels(int64_t n, int32_t base, int32_t *m, int32_t *e, int32_t cap)
{
    if (base < 2 || n < 0) return -1;
    int32_t k = 0;
    for (int32_t pos = 0; n; ++pos) {
        int64_t d = n % base;
        n /= base;
        if (d) {
            if (k >= cap) return -2;
            m[k] = (int32_t)d;
            e[k] = pos;
            ++k;
        }
    }
    return k;
}

int64_t phi_wheels_batch(const int64_t *n, int32_t count, int32_t base,
                         int32_t *m, int32_t *e, int32_t *cnt, int64_t cap)
{
    if (base < 2 || count < 0) return -1;
    int64_t off = 0;
    for (int32_t i = 0; i < count; ++i) {
        int64_t room = cap - off;
        int32_t k = phi_wheels(n[i], base, m + off, e + off,
                               room > INT32_MAX ? INT32_MAX : (int32_t)room);
        if (k < 0) return k;
        cnt[i] = k;
        off += k;
    }
    return off;
}

static int pow_checked(int64_t b, int64_t e, int64_t *out)
{
    int64_t r = 1;
    while (e-- > 0) {
        if (r > INT64_MAX / b) return -3;
        r *= b;
    }
    *out = r;
    return 0;
}

/* значение наследственной формы n (записанной в base), прочитанной в newbase;
 * depth и wheels — как у phicore.wheels.Hereditary: depth 0 для n = 0,
 * иначе 1 + max глубины поддеревьев показателей; wheels — колёса всех уровней */
static int hered(int64_t n, int32_t base, int32_t newbase,
                 int64_t *value, int32_t *depth, int64_t *wheels)
{
    int64_t total = 0, wc = 0;
    int32_t dmax = 0, any = 0;
    for (int32_t pos = 0; n; ++pos) {
        int64_t d = n % base;
        n /= base;
        if (!d) continue;
        any = 1;
        wc += 1;
        int64_t e = 0;
        int32_t sd = 0;
        if (pos >= 1) {
            int64_t ev, sw;
            int rc = hered(pos, base, newbase, &ev, &sd, &sw);
            if (rc) return rc;
            e = ev;
            wc += sw;
        }
        int64_t p;
        if (pow_checked(newbase, e, &p)) return -3;
        if (d > INT64_MAX / p) return -3;
        int64_t term = d * p;
        if (total > INT64_MAX - term) return -3;
        total += term;
        if (sd > dmax) dmax = sd;
    }
    *value = total;
    *depth = any ? dmax + 1 : 0;
    *wheels = wc;
    return 0;
}

int phi_hereditary(int64_t n, int32_t base, int32_t newbase,
                   int64_t *value, int32_t *depth, int64_t *wheels)
{
    if (base < 2 || newbase < 2 || n < 0) return -1;
    return hered(n, base, newbase, value, depth, wheels);
}

/* Последовательность Гудстейна: out[0] = n, далее шаг (base k+2 → k+3) − 1,
 * останов на нуле; возвращает число записанных членов или −3 при переполнении */
int32_t phi_goodstein(int64_t n, int32_t steps, int64_t *out)
{
    if (n < 0 || steps < 0) return -1;
    out[0] = n;
    int32_t k = 1;
    int32_t base = 2;
    while (k <= steps && out[k - 1] != 0) {
        int64_t v, w;
        int32_t d;
        int rc = hered(out[k - 1], base, base + 1, &v, &d, &w);
        if (rc) return rc;
        out[k++] = v - 1;
        base += 1;
    }
    return k;
}
