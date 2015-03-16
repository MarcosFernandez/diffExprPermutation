#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <pthread.h>
#include "perm.h"
#include <x86intrin.h>

#define lfact(x) lfact_store[(int)(x)]

#define HAVE_SSE 1

/**
 *\brief Compares two float values
 *\param s1 First number to compare
 *\param s2 Second number to compare
 *\return -1 if s2 is the biggest, 1 if s1 is the biggest, 0 if both are equal
 */
static int cmp_float(const void* s1, const void* s2)
{
    const float* f1, *f2;
    f1 = s1;
    f2 = s2;
    if (*f1 < *f2) return -1;
    if (*f1 > *f2) return 1;
    return 0;
}

#ifdef HAVE_SSE

/* Streaming SIMD Extensions (SSE) is a SIMD instruction set extension to the x86 architecture, designed by Intel and introduced in 1999 in their Pentium III series processors as a reply to AMD's 3DNow. SSE contains 70 new instructions, most of which work on single precision floating point data. SIMD instructions can greatly increase performance when exactly the same operations are to be performed on multiple data objects. */

#define SET_MASKS(z)                            \
    const __m128 sign_mask = _mm_set1_ps(-0.f); \
    __m128d pmask[4];                           \
    pmask[0] = _mm_set_pd(0.0, 0.0);            \
    pmask[1] = _mm_set_pd(0.0, z);              \
    pmask[2] = _mm_set_pd(z, 0.0);              \
    pmask[3] = _mm_set_pd(z, z);

#define UPDATE_PVAL(x, i, ms, pv)    \
    int msk = _mm_movemask_ps(x);    \
    __m128d p = _mm_load_pd(pv + i); \
    p = _mm_add_pd(p, ms[msk & 3]);  \
    _mm_store_pd(pv + i, p);         \
    p = _mm_load_pd(pv + 2 + i);     \
    p = _mm_add_pd(p, ms[msk >> 2]); \
    _mm_store_pd(pv + 2 + i, p);

/**
 *\brief Percentile selection for two ranks
 *\extends low level operation using x86 intrisic library. More info at: http://www-db.in.tum.de/~finis/x86-intrin-cheatsheet-v2.1.pdf
 */
static inline void check_pc_2(double z, float* xp, float* yp, int nv,
    float* diff, double* pval)
{
    SET_MASKS(z);
    for (int i = 0; i < nv; i += 4) 
    {
        __m128 tx = _mm_load_ps(xp + i);
        tx = _mm_sub_ps(tx, _mm_load_ps(yp + i));
        tx = _mm_andnot_ps(sign_mask, tx);
        tx = _mm_cmpnle_ps(tx, _mm_load_ps(diff + i));
        UPDATE_PVAL(tx, i, pmask, pval);
    }
}

/**
 *\brief Percentile selection for three ranks
 *\extends low level operation using x86 intrisic library. More info at: http://www-db.in.tum.de/~finis/x86-intrin-cheatsheet-v2.1.pdf
 */
static inline void check_pc_3(double z, float* xp, float* xp1, float* yp,
    float xpc, int nv, float* diff, double* pval)
{
    SET_MASKS(z);
    __m128 frac = _mm_set1_ps(xpc);
    for (int i = 0; i < nv; i += 4) 
    {
        __m128 tx = _mm_load_ps(xp + i);
        __m128 tx1 = _mm_load_ps(xp1 + i);
        tx1 = _mm_sub_ps(tx1, tx);
        tx1 = _mm_mul_ps(tx1, frac);
        tx = _mm_add_ps(tx, tx1);
        tx = _mm_sub_ps(tx, _mm_load_ps(yp + i));
        tx = _mm_andnot_ps(sign_mask, tx);
        tx = _mm_cmpnle_ps(tx, _mm_load_ps(diff + i));
        UPDATE_PVAL(tx, i, pmask, pval);
    }
}

/**
 *\brief Percentile selection for four ranks
 *\extends low level operation using x86 intrisic library. More info at: http://www-db.in.tum.de/~finis/x86-intrin-cheatsheet-v2.1.pdf
 */
static inline void check_pc_4(double z, float* xp, float* xp1, float* yp,
    float* yp1, float xpc, float ypc, int nv,
    float* diff, double* pval)
{
    SET_MASKS(z);
    __m128 xfrac = _mm_set1_ps(xpc);
    __m128 yfrac = _mm_set1_ps(ypc);
    for (int i = 0; i < nv; i += 4) {
        __m128 tx = _mm_load_ps(xp + i);
        __m128 tx1 = _mm_load_ps(xp1 + i);
        tx1 = _mm_sub_ps(tx1, tx);
        tx1 = _mm_mul_ps(tx1, xfrac);
        tx = _mm_add_ps(tx, tx1);
        __m128 tx2 = _mm_load_ps(yp + i);
        tx1 = _mm_load_ps(yp1 + i);
        tx1 = _mm_sub_ps(tx1, tx2);
        tx1 = _mm_mul_ps(tx1, yfrac);
        tx1 = _mm_add_ps(tx1, tx2);
        tx = _mm_sub_ps(tx, tx1);
        tx = _mm_andnot_ps(sign_mask, tx);
        tx = _mm_cmpnle_ps(tx, _mm_load_ps(diff + i));
        UPDATE_PVAL(tx, i, pmask, pval);
    }
}

#else

/*NO SSE instructions */

/**
 *\brief Percentile selection for two ranks
 */
static inline void check_pc_2(double z, float* xp, float* yp, int nv,
    float* diff, double* pval)
{
    for (int i = 0; i < nv; i++) {
        float df = fabs(*xp++ - *yp++);
        if (df >= diff[i])
            pval[i] += z;
    }
}

/**
 *\brief Percentile selection for three ranks
 */
static inline void check_pc_3(double z, float* xp, float* xp1, float* yp,
    float xpc, int nv, float* diff, double* pval)
{
    for (int i = 0; i < nv; i++) {
        float tx = *xp++;
        float df = fabs(tx + xpc * (*xp1++ - tx) - *yp++);
        if (df >= diff[i])
            pval[i] += z;
    }
}

/**
 *\brief Percentile selection for four ranks
 */
static inline void check_pc_4(double z, float* xp, float* xp1, float* yp,
    float* yp1, float xpc, float ypc, int nv,
    float* diff, double* pval)
{
    for (int i = 0; i < nv; i++) {
        float tx = *xp++;
        float ty = *yp++;
        float df = fabs(tx + xpc * (*xp1++ - tx) - (ty + ypc * (*yp1++ - ty)));
        if (df >= diff[i])
            pval[i] += z;
    }
}

#endif

/**
 *\brief Cycle through all posible values for the percintile in bouth group of samples
 *\param int sx Number of values before th required percentile for first group
 *\param int sx1
 *\param float xpct Precentile value first sample
 *\param int sy Number of values before th required percentile for second group
 *\param int sy1
 *\param int nv Number of genes
 *\param int n Number of Samples
 *\param int konst
 *\param float** expr Array of expression values per each gene
 *\param float* diff Array of differences per gene
 *\param double* pval Array of pvalues per gene
 *\param double * lfact_store
 *\param int step,
 *\param int start
 */ 
static double cycle_perm(int sx, int sx1, float xpct, int sy, int sy1,
    float ypct, int nv, int n, double konst, float** expr,
    float* diff, double* pval, double* lfact_store,
    int step, int start)
{
    double tot = 0.0;
    if (xpct > 1.0e-6) {
        if (ypct > 1.0e-6) {
            // Both groups require averaging two values
            for (int ix = sx + start; ix < n - sx1 - 1; ix += step) {
                if (ix - sx <= sy) {
                    // Median of second group comes after that of first group
                    double z1 = lfact(ix) - lfact(sx) - lfact(ix - sx) - konst;
                    for (int iy = sx + sy + 1; iy < n - sy1 - 1; iy++) {
                        int d = iy - ix - 1;
                        int cx = d - sy + ix - sx;
                        if (cx) {
                            double z2 = z1 - lfact(cx - 1) - lfact(sy1);
                            int j = sy + sx - ix;
                            for (int ix1 = ix + 1; ix1 <= iy - cx; ix1++, j--) {
                                int k = n - iy - 2 - sy1;
                                double z = exp(z2 + lfact(cx - 1 + j) - lfact(j) + lfact(k + sy1) - lfact(k));
                                int iy1;
                                for (iy1 = iy + 1; iy1 < n - sy1 - 1 && z > 0.0; iy1++, k--) {
                                    tot += z;
                                    check_pc_4(z, expr[ix], expr[ix1], expr[iy], expr[iy1], xpct,
                                        ypct, nv, diff, pval);
                                    z *= (double)k / (double)(k + sy1);
                                }
                                check_pc_4(z, expr[ix], expr[ix1], expr[iy], expr[iy1], xpct,
                                    ypct, nv, diff, pval);
                                tot += z;
                            }
                        }
                        else {
                            int d = iy - ix - 1;
                            int k1 = sy + sx - ix;
                            int k = d - k1;
                            int ix1 = iy + 1;
                            double z2 = z1 + lfact(k + k1) - lfact(k) - lfact(k1);
                            double z3 = z2 - lfact(sy1);
                            k = n - iy - 3 - sy1;
                            double z = exp(z3 + lfact(k + sy1) - lfact(k));
                            int iy1;
                            for (iy1 = iy + 2; iy1 < n - sy1 - 1 && z > 0.0; iy1++, k--) {
                                tot += z;
                                check_pc_4(z, expr[ix], expr[ix1], expr[iy], expr[iy1], xpct,
                                    ypct, nv, diff, pval);
                                z *= (double)k / (double)(sy1 + k);
                            }
                            check_pc_4(z, expr[ix], expr[ix1], expr[iy], expr[iy1], xpct,
                                ypct, nv, diff, pval);
                            tot += z;
                            z3 = z2 - lfact(sx1);
                            iy1 = iy + 1;
                            k = n - iy - 3 - sx1;
                            z = exp(z3 + lfact(sx1 + k) - lfact(k));
                            for (ix1 = iy + 2; ix1 < n - sx1 - 1 && z > 0.0; ix1++, k--) {
                                tot += z;
                                check_pc_4(z, expr[ix], expr[ix1], expr[iy], expr[iy1], xpct,
                                    ypct, nv, diff, pval);
                                z *= (double)k / (double)(sx1 + k);
                            }
                            check_pc_4(z, expr[ix], expr[ix1], expr[iy], expr[iy1], xpct,
                                ypct, nv, diff, pval);
                            tot += z;
                        }
                    }
                }
                else {
                    // Median of second group comes before that of first group
                    for (int iy = sy; iy <= sx + sy; iy++) {
                        int d = ix - iy - 1;
                        int cy = d - sx + iy - sy;
                        if (cy) {
                            int j = sx + sy - iy;
                            double z1 = lfact(iy) - lfact(sy) - lfact(iy - sy) - lfact(cy - 1) - lfact(sx1) - konst;
                            for (int iy1 = iy + 1; iy1 <= ix - cy; iy1++, j--) {
                                int k = n - ix - 2 - sx1;
                                double z = exp(z1 + lfact(j + cy - 1) + lfact(sx1 + k) - lfact(k) - lfact(j));
                                int ix1;
                                for (ix1 = ix + 1; ix1 < n - sx1 - 1 && z > 0.0; ix1++, k--) {
                                    tot += z;
                                    check_pc_4(z, expr[ix], expr[ix1], expr[iy], expr[iy1], xpct,
                                        ypct, nv, diff, pval);
                                    z *= (double)k / (double)(sx1 + k);
                                }
                                check_pc_4(z, expr[ix], expr[ix1], expr[iy], expr[iy1], xpct,
                                    ypct, nv, diff, pval);
                                tot += z;
                            }
                        }
                        else {
                            int d = ix - iy - 1;
                            int k = sx + sy - iy;
                            int k1 = d - k;
                            int iy1 = ix + 1;
                            double z1 = lfact(iy) - lfact(sy) - lfact(iy - sy) + lfact(k + k1) - lfact(k) - lfact(k1) - konst;
                            k = n - ix - 3 - sx1;
                            double z = exp(z1 - lfact(sx1) + lfact(sx1 + k) - lfact(k));
                            int ix1;
                            for (ix1 = ix + 2; ix1 < n - sx1 - 1 && z > 0.0; ix1++, k--) {
                                tot += z;
                                check_pc_4(z, expr[ix], expr[ix1], expr[iy], expr[iy1], xpct,
                                    ypct, nv, diff, pval);
                                z *= (double)k / (double)(sx1 + k);
                            }
                            check_pc_4(z, expr[ix], expr[ix1], expr[iy], expr[iy1], xpct,
                                ypct, nv, diff, pval);
                            tot += z;
                            k = n - ix - 3 - sy1;
                            z = exp(z1 - lfact(sy1) + lfact(sy1 + k) - lfact(k));
                            ix1 = ix + 1;
                            for (iy1 = ix + 2; iy1 < n - sy1 - 1 && z > 0.0; iy1++, k--) {
                                tot += z;
                                check_pc_4(z, expr[ix], expr[ix1], expr[iy], expr[iy1], xpct,
                                    ypct, nv, diff, pval);
                                z *= (double)k / (double)(sy1 + k);
                            }
                            check_pc_4(z, expr[ix], expr[ix1], expr[iy], expr[iy1], xpct,
                                ypct, nv, diff, pval);
                            tot += z;
                        }
                    }
                }
            }
        }
        else {
            // Group 1 only requires averaging two values
            for (int ix = sx + start; ix < n - sx1 - 1; ix += step) {
                if (ix - sx <= sy) {
                    // Median of second group comes after that of first group
                    double z1 = lfact(ix) - lfact(sx) - lfact(ix - sx) - konst;
                    for (int iy = sx + sy + 1; iy < n - sy1; iy++) {
                        int d = iy - ix - 1;
                        int cx = d - sy - sx + ix;
                        if (cx) {
                            int j = n - iy - 1 - sy1;
                            int k = sx + sy - ix;
                            double z = exp(z1 + lfact(n - iy - 1) - lfact(j) - lfact(sy1) + lfact(cx - 1 + k) - lfact(cx - 1) - lfact(k));
                            int ix1;
                            for (ix1 = ix + 1; ix1 <= iy - cx - 1 && z > 0.0; ix1++, k--) {
                                tot += z;
                                check_pc_3(z, expr[ix], expr[ix1], expr[iy], xpct, nv, diff,
                                    pval);
                                z *= (double)k / (double)(cx - 1 + k);
                            }
                            check_pc_3(z, expr[ix], expr[ix1], expr[iy], xpct, nv, diff,
                                pval);
                            tot += z;
                        }
                        else {
                            int k = n - iy - 2 - sx1;
                            double z = exp(z1 + lfact(sx1 + k) - lfact(sx1) - lfact(k));
                            int ix1;
                            for (ix1 = iy + 1; ix1 < n - sx1 - 1 && z > 0.0; ix1++, k--) {
                                tot += z;
                                check_pc_3(z, expr[ix], expr[ix1], expr[iy], xpct, nv, diff,
                                    pval);
                                z *= (double)k / (double)(sx1 + k);
                            }
                            tot += z;
                            check_pc_3(z, expr[ix], expr[ix1], expr[iy], xpct, nv, diff,
                                pval);
                        }
                    }
                }
                else {
                    // Median of second group comes before that of first group
                    for (int iy = sy; iy <= sx + sy; iy++) {
                        int d = ix - iy - 1;
                        int j = sx + sy - iy;
                        int k = n - ix - 2 - sx1;
                        double z = exp(lfact(iy) - lfact(sy) - lfact(iy - sy) + lfact(d) - lfact(j) - lfact(d - j) + lfact(n - ix - 2) - lfact(sx1) - lfact(k) - konst);
                        int ix1;
                        for (ix1 = ix + 1; ix1 < n - sx1 - 1 && z > 0.0; ix1++, k--) {
                            tot += z;
                            check_pc_3(z, expr[ix], expr[ix1], expr[iy], xpct, nv, diff,
                                pval);
                            z *= (double)k / (double)(sx1 + k);
                        }
                        check_pc_3(z, expr[ix], expr[ix1], expr[iy], xpct, nv, diff, pval);
                        tot += z;
                    }
                }
            }
        }
    }
    else if (ypct > 1.0e-6) {
        fprintf(stderr, "Internal error\n");
        exit(-1);
    }
    else {
        // Both groups require single values
        int ix, iy;
        int k = sy - start;
        // Median of second group comes after that of first group
        for (int j = start; j <= sy; j += step, k -= step) {
            ix = sx + j;
            int k1 = n - sx - 2 - sy - sy1;
            int k2 = 0;
            double z = exp(lfact(ix) - lfact(sx) - lfact(ix - sx) + lfact(k1 + sy1) - lfact(k1) - lfact(sy1) - konst);
            for (iy = sx + sy + 1; iy < n - sy1 - 1 && z > 0.0; iy++, k1--, k2++) {
                tot += z;
                check_pc_2(z, expr[ix], expr[iy], nv, diff, pval);
                z *= (double)(k1 * (k + k2 + 1)) / (double)((sy1 + k1) * (k2 + 1));
            }
            check_pc_2(z, expr[ix], expr[iy], nv, diff, pval);
            tot += z;
        }
        // Median of second group comes before that of first group
        for (ix = sx + sy + 1 + start; ix < n - sx1; ix += step) {
            int k = ix - sx - sy - 1;
            int k2 = sx;
            double z = exp(lfact(n - ix - 1) - lfact(sx1) - lfact(n - ix - sx1 - 1) + lfact(sy + k2) - lfact(sy) - lfact(k2) - konst);
            iy = sx + sy;
            for (int k1 = 0; k1 < sx && z > 0.0; k2--, k1++, iy--) {
                tot += z;
                check_pc_2(z, expr[ix], expr[iy], nv, diff, pval);
                z *= (double)(k2 * (k + k1 + 1)) / (double)((sy + k2) * (k1 + 1));
            }
            check_pc_2(z, expr[ix], expr[iy], nv, diff, pval);
            tot += z;
        }
    }
    return tot;
}

typedef struct {
    int sx, sx1, sy, sy1;
    float xpct, ypct;
    int nv;
    int n;
    double konst;
    float** expr;
    float* diff;
    double* pval;
    double* lfact;
    int step;
    int start;
    double tot;
} cycle_perm_data;

static void* call_cycle_perm(void* s)
{
    cycle_perm_data* cpd = s;
    cpd->tot = cycle_perm(cpd->sx, cpd->sx1, cpd->xpct, cpd->sy, cpd->sy1, cpd->ypct,
        cpd->nv, cpd->n, cpd->konst, cpd->expr, cpd->diff, cpd->pval,
        cpd->lfact, cpd->step, cpd->start);
    return 0;
}

/**
 *\brief Calculate exact permutation p-values
 *\param double pct - percentile (0.5 median, 0.25 percentile 25%,...)
 *\param cp_mode mode - Percentile selection [Nearest/Linear]
 *\param n_threas - number of threads
 */
int check_percentile(double pct, cp_mode mode, int n_threads,
    struct perm_data* pd)
{
    int a = pd->n_cases;
    int b = pd->n_controls;

    /* sx and sy are the number of values *before* the
     required percentile for group 1 and 2.  If
     xpct/ypct >0 then we need to combine two adjacent
     values so z = (1-xpct)*X[sx]+xpct*X[sx+1]; */

    int sx, sy;
    float xpct, ypct;
    double z = (double)(pct * a) / 100 - 0.5;
    if (z < 0.0)
        z = 0.0;
    sx = (int)z;
    if (sx >= a - 1) {
        sx = a - 1;
        xpct = 0.0;
    }
    else {
        xpct = z - (double)sx;
    }
    z = (double)(pct * b) / 100 - 0.5;
    if (z < 0.0)
        z = 0.0;
    sy = (int)z;
    if (sy >= b - 1) {
        sy = b - 1;
        ypct = 0.0;
    }
    else {
        ypct = z - (double)sy;
    }
    int n = a + b;
    if (mode == AUTO)
        mode = n > 100 ? NEAREST_RANK : LINEAR_INTERPOLATION;
    if (mode == NEAREST_RANK) {
        sx += (int)(0.5 + xpct);
        sy += (int)(0.5 + ypct);
        xpct = ypct = 0.0;
    }
    int nv = pd->n_var;
    int nv1 = (((nv + 3) >> 2) << 2);

    // Now we calculate the original percentila differences for each gene
    // and after, we create a new array of expression values which is
    // sorted independently for each gene and transposed from the original order
    // to make it faster to access all the values corresponding to a given rank
    float** expr = malloc(sizeof(void*) * n);
    for (int i = 0; i < n; i++) {
        expr[i] = malloc(sizeof(float) * nv1);
    }
    float* diff = malloc(sizeof(float) * nv1);
    for (int i = nv; i < nv1; i++) {
        diff[i] = FLT_MAX;
    }
    float* tmp_cases = malloc(sizeof(float) * n);
    float* tmp_controls = tmp_cases + a;
    float* tp = pd->data;
    for (int i = 0; i < nv; i++) {
        int k1 = 0, k2 = 0;
        for (int j = 0; j < n; j++) {
            if (pd->grp[j] == CASE) {
                if (k1 == a)
                    break;
                tmp_cases[k1++] = *tp++;
            }
            else {
                if (k2 == b)
                    break;
                tmp_controls[k2++] = *tp++;
            }
        }
        if (k1 != a || k2 != b) {
            fprintf(
                stderr,
                "Bad numbers of cases & controls! (found %d,%d - expected %d,%d)\n",
                k1, k2, a, b);
            exit(-1);
        }
        qsort(tmp_cases, a, sizeof(float), cmp_float);
        pd->pc_cases[i] = xpct > 1.0e-6
            ? tmp_cases[sx] + xpct * (tmp_cases[sx + 1] - tmp_cases[sx])
            : tmp_cases[sx];
        qsort(tmp_controls, b, sizeof(float), cmp_float);
        pd->pc_controls[i] = ypct > 1.0e-6
            ? tmp_controls[sy] + ypct * (tmp_controls[sy + 1] - tmp_controls[sy])
            : tmp_controls[sy];
        pd->pc_diff[i] = pd->pc_cases[i] - pd->pc_controls[i];
        diff[i] = fabs(pd->pc_diff[i]);
        pd->fold_change[i] = pd->pc_cases[i] / pd->pc_controls[i];
        int j = 0, k = 0;
        for (int k1 = 0; k1 < n; k1++) {
            float z;
            if (j < a && (k >= b || tmp_cases[j] <= tmp_controls[k]))
                z = tmp_cases[j++];
            else
                z = tmp_controls[k++];
            expr[k1][i] = z;
        }
    }
    free(tmp_cases);
    if (ypct > 1.0e-6 && xpct <= 1.0e-6) {
        int t;
        t = a;
        a = b;
        b = t;
        t = sx;
        sx = sy;
        sy = t;
        xpct = ypct;
        ypct = 0.0;
    }
    int sx1 = (xpct > 1.0e-6) ? a - sx - 2 : a - sx - 1;
    int sy1 = (ypct > 1.0e-6) ? b - sy - 2 : b - sy - 1;
    //    printf("%d %d %g  %d %d %g  nv = %d\n", sx, sx1, xpct, sy, sy1, ypct, nv);
    static double* lfact_store;

    lfact_store = malloc(sizeof(double) * (n + 1));
    for (int i = 0; i <= n; i++)
        lfact_store[i] = lgamma((double)(i + 1));

    double konst = lfact(n) - lfact(a) - lfact(b);
    double tot;
    double* pval;
    if (n_threads == 1) {
        pval = malloc(sizeof(double) * nv1);
        for (int i = 0; i < nv1; i++)
            pval[i] = 0.0;
        tot = cycle_perm(sx, sx1, xpct, sy, sy1, ypct, nv, n, konst, expr, diff,
            pval, lfact_store, 1, 0);
    }
    else {
        pval = malloc(sizeof(double) * nv1 * n_threads);
        cycle_perm_data* cpd = malloc(sizeof(cycle_perm_data) * n_threads);

        cpd[0].sx = sx;
        cpd[0].sx1 = sx1;
        cpd[0].sy = sy;
        cpd[0].sy1 = sy1;
        cpd[0].xpct = xpct;
        cpd[0].ypct = ypct;
        cpd[0].nv = nv;
        cpd[0].n = n;
        cpd[0].konst = konst;
        cpd[0].expr = expr;
        cpd[0].diff = diff;
        cpd[0].lfact = lfact_store;
        cpd[0].step = n_threads;

        pthread_t* pth = malloc(sizeof(pthread_t) * n_threads);
        for (int th = 0; th < n_threads; th++) {
            for (int i = 0; i < nv1; i++)
                pval[nv1 * th + i] = 0.0;

            if (th)
                memcpy(cpd + th, cpd, sizeof(cycle_perm_data));
            cpd[th].pval = pval + nv1 * th;
            cpd[th].start = th;
            int i = pthread_create(pth + th, NULL, call_cycle_perm, cpd + th);
            if (i) {
                fprintf(stderr, "Thread creation failed: aborting\n");
                exit(-1);
            }
        }
        for (int th = 0; th < n_threads; th++)
            pthread_join(pth[th], NULL);
        tot = cpd[0].tot;
        for (int th = 1; th < n_threads; th++) {
            tot += cpd[th].tot;
            double* tp = cpd[th].pval;
            for (int i = 0; i < nv; i++)
                pval[i] += tp[i];
        }
        free(cpd);
        free(pth);
    }
    //    printf("tot = %g, %g\n", log(tot), tot);
    for (int i = 0; i < nv; i++) {
        //        printf("%d %f %g\n", i, pd->pc_diff[i], pval[i]);
        pd->p_values[i] = pval[i];
    }
    free(pval);
    free(lfact_store);
    free(diff);
    for (int i = 0; i < n; i++)
        free(expr[i]);
    free(expr);
    return 0;
}
