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
_PROTO_SWAP(T, swap_val);
_PROTO_SWAP(int, swap_int);

static void
swap_rows(struct lin_eqn *eq, int i, int j)
{
    swap_ptr((void**)eq->A, i, j);
    swap_val(eq->b, i, j);
}
/* layout is as following:
 * +--------------------+----+-----+
 * | vector of row ptrs | T* | n   |
 * +--------------------+----+-----+
 * | vector of values   | T  | n   |
 * +--------------------+----+-----+
 * | matrix data        | T  | n^2 |
 * +--------------------+----+-----+
 */

/* offset of a vector of pointers to rows */
#define _OFFSET_ROWS(n) 0
/* offset of a vector of values */
#define _OFFSET_VALS(n) ( (n) * sizeof(T*) )
/* offset to matrix contents */
#define _OFFSET_MATRDATA(n) ( _OFFSET_VALS(n) + (n)*sizeof(T) )
/* size of buffer to allocate in bytes */
#define _LIN_EQN_BUFSIZE(n) ( _OFFSET_MATRDATA(n) + (n)*(n)*sizeof(T) )

/* allocate and initialize system */
void
lin_eqn_init(struct lin_eqn *eq, size_t n)
{
    eq->n = n;
    eq->fwd_n = eq->bck_n = 0;
    /* 
     * 
     */
    size_t buflen = _LIN_EQN_BUFSIZE(n);
    eq->buf = malloc(buflen);
    if (eq->buf == NULL)
        err(EXIT_FAILURE, "lin_eqn_init(): cannot malloc()");
    eq->A = (T**) eq->buf;
    eq->b = (T*) eq->buf + _OFFSET_VALS(eq->n);
    for (int i = 0; i < n; ++i) {
        eq->A[i] = (T*)(eq->buf + _OFFSET_MATRDATA(eq->n)) + i*eq->n;
        eq->b[i] = 0;
    }
}

/* sets up pivot for forward sweep */
/* result: A[k][k] has maximal absoute value in its row */
static void
setup_pivot(struct lin_eqn *eq)
{
    int max = eq->fwd_n;
    T maxabs = tabs(eq->fwd_n);
    for (int i = eq->fwd_n; i < eq->n; ++i)
    {
        register T curabs = tabs(eq->A[i][eq->fwd_n]);
        if (curabs > maxabs) {
            max = i;
            maxabs = curabs;
        }
    }
    if (max != eq->fwd_n)
        swap_rows(eq, eq->fwd_n, max);
}

static void
lin_mul_row(struct lin_eqn *eq, int row, T mul)
{
    for (int i = 0; i < eq->n; ++i)
        eq->A[row][i] *= mul;
    eq->b[row] *= mul;
}

static void
lin_add_row(struct lin_eqn *eq, int src, int dst, T mul)
{
    for (int i = 0; i < eq->n; ++i)
        eq->A[dst][i] += mul * eq->A[src][i];
    eq->b[dst] += mul * eq->b[src];
}

int
lin_fwd_sweep(struct lin_eqn *eq)
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
lin_bck_sweep(struct lin_eqn *eq)
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

int linsolve(struct lin_eqn *eq)
{
    int ret = lin_fwd_sweep(eq);
    if (ret != LINSOLVE_SUCCESS)
        return ret;
    ret = lin_bck_sweep(eq);
    /* For columns it'd be necessary to:
     *
     * for (int i = 0; i < eq->n; ++i) {
     *     register T* base = (T*) (eq->buf + _OFFSET_MATRDATA(eq->n));
     *     int orig_row = (eq->A[i] - base)/eq->n;
     *     if (i != orig_row)
     *         swap_rows(eq, i, orig_row);
     * }
     */
    eq->x = eq->b;
    return ret;
}
