#include <inttypes.h>
#include <stddef.h>
#include <vsu/num.h>

#ifndef _VSU_LINSOLVE_H_
#define _VSU_LINSOLVE_H_

#define LINSOLVE_SUCCESS 0
#define LINSOLVE_DEGENERATE 1

typedef double T;
#define TFORMAT "%f"

/* GAUSSIAN ELIMINATION */

struct lin_eqn
{
    /* internal buffer; lin_eqn MUST BE disposed via free(buf) */
    void *buf;

    /* matrix is represented as
     * an array of rows (allows fastest swapping)
     */
    T **A;

    T *b; /* values vector */
    T *x; /* after system's solved, `x` points to an answer */

    /* permutation; initialized with [1, ..., n] */
    /* unnecessary? */
    /* int *p; */

    /* forward-sweep counter;
     * degree of largest (yet) lower-triangulated principal minor
     */
    size_t fwd_n;

    size_t bck_n; /* backward-sweep counter; */
    size_t n; /* dimension of matrix */
};

static inline T
tabs(T t)
{
    if (t < 0)
        return -t;
    return t;
}

void
lin_eqn_init(struct lin_eqn *eq, size_t n);

int
lin_bck_sweep(struct lin_eqn *eq);

int
lin_fwd_sweep(struct lin_eqn *eq);

int
linsolve(struct lin_eqn *eq);

/* TRIDIAGONAL */
/* a_i x_{i-1} + b_i x_i + c_i x_{i+1} = f_i */
struct tridiag_eqn
{
    RR a, b, c, f;
    RR x;
};

int
tridiag_eqn_init(struct tridiag_eqn *eq, size_t n);

int
tridiag_eqn_solve(struct tridiag_eqn *eq);

#endif
