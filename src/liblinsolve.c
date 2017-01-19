#include <vsu/linsolve.h>
#include <stdlib.h>
#include <err.h>

#define _PROTO_SWAP(TYP,FUNCNAME) \
    static inline void \
    FUNCNAME(TYP *arr, int i, int j) \
{ \
    TYP t = arr[i]; \
    arr[i] = arr[j]; \
    arr[j] = t; \
}

_PROTO_SWAP(void*, swap_ptr);
_PROTO_SWAP(RR, swap_val);
_PROTO_SWAP(int, swap_int);

static void
swap_rows(struct lin_eqn *eq, int i, int j)
{
    swap_ptr((void**)eq->A, i, j);
    swap_val(eq->b, i, j);
}
/* layout is as following:
 * +--------------------+----+-----+
 * | vector of row ptrs | RR* | n   |
 * +--------------------+----+-----+
 * | vector of values   | RR  | n   |
 * +--------------------+----+-----+
 * | matrix data        | RR  | n^2 |
 * +--------------------+----+-----+
 */

/* allocate and initialize system */
lin_eqn_ptr
lin_eqn_alloc(size_t n)
{
    size_t buflen = sizeof(struct lin_eqn) +
        sizeof(RR)*(n*(n+1));
    void *buf = malloc(buflen);
    if (buf == NULL)
        err(EXIT_FAILURE, "lin_eqn_init(): cannot malloc()");
    lin_eqn_ptr eq = buf;
    eq->buf = buf;
    eq->n = n;
    eq->fwd_n = eq->bck_n = 0;
    eq->A = eq->buf + sizeof(struct lin_eqn);
    eq->b = (RR*)eq->A + n;
    RR *col = eq->b + n;
    for (int i = 0; i < n; ++i) {
        eq->A[i] = col;
        eq->b[i] = 0;
        col += n;
    }
    return eq;
}

/* sets up pivot for forward sweep */
/* result: A[k][k] has maximal absoute value in its row */
static void
setup_pivot(struct lin_eqn *eq)
{
    int max = eq->fwd_n;
    RR maxabs = tabs(eq->fwd_n);
    for (int i = eq->fwd_n; i < eq->n; ++i)
    {
        register RR curabs = tabs(eq->A[i][eq->fwd_n]);
        if (curabs > maxabs) {
            max = i;
            maxabs = curabs;
        }
    }
    if (max != eq->fwd_n)
        swap_rows(eq, eq->fwd_n, max);
}

static void
lin_mul_row(struct lin_eqn *eq, int row, RR mul)
{
    for (int i = 0; i < eq->n; ++i)
        eq->A[row][i] *= mul;
    eq->b[row] *= mul;
}

static void
lin_add_row(lin_eqn_ptr eq, int src, int dst, RR mul)
{
    for (int i = 0; i < eq->n; ++i)
        eq->A[dst][i] += mul * eq->A[src][i];
    eq->b[dst] += mul * eq->b[src];
}

int
lin_fwd_sweep(lin_eqn_ptr eq)
{
    int k = eq->fwd_n;
    while (k < eq->n)
    {
        setup_pivot(eq);
        if (eq->A[k][k] == 0) {
            return LINSOLVE_DEGENERATE;
        }
        for (int i = k+1; i < eq->n; ++i) {
            if (eq->A[i][k] == 0) continue;
            lin_add_row(eq, k, i, - eq->A[i][k]/eq->A[k][k]);
            eq->A[i][k] = 0;
        }
        /* to avoid values growing too much in process */
        lin_mul_row(eq, k, 1/eq->A[k][k]);
        k++;
        eq->fwd_n = k;
    }
    return LINSOLVE_SUCCESS;
}

int
lin_bck_sweep(lin_eqn_ptr eq)
{
    int k = eq->n - eq->bck_n - 1;
    while (k > 0)
    {
        for (int i = k - 1; i >= 0; --i) {
            if (eq->A[i][k] == 0) continue;
            lin_add_row(eq, k, i, -1*eq->A[i][k]/eq->A[k][k]);
            eq->A[i][k] = 0;
        }
        lin_mul_row(eq, k, 1/eq->A[k][k]);
        ++eq->bck_n;
        --k;
    }
    return LINSOLVE_SUCCESS;
}

int linsolve(lin_eqn_ptr eq)
{
    int ret = lin_fwd_sweep(eq);
    if (ret != LINSOLVE_SUCCESS)
        return ret;
    ret = lin_bck_sweep(eq);
    /* For columns it'd be necessary to:
     *
     * for (int i = 0; i < eq->n; ++i) {
     *     register RR* base = (RR*) (eq->buf + _OFFSET_MATRDATA(eq->n));
     *     int orig_row = (eq->A[i] - base)/eq->n;
     *     if (i != orig_row)
     *         swap_rows(eq, i, orig_row);
     * }
     */
    eq->x = eq->b;
    return ret;
}
