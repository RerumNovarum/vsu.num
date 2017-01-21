#include <inttypes.h>
#include <stdbool.h>
#include <stddef.h>
#include <vsu/num.h>

#ifndef _VSU_LINSOLVE_H_
#define _VSU_LINSOLVE_H_

#define LINSOLVE_SUCCESS 0
#define LINSOLVE_DEGENERATE 1

/* GAUSSIAN ELIMINATION */

struct lin_eqn
{
    /* internal buffer; lin_eqn MUST BE disposed via free(buf) */
    void *buf;

    /* matrix is represented as
     * an array of rows (allows fastest swapping)
     */
    RR **A;

    RR *b; /* values vector */
    RR *x; /* once system's solved, `x` points to an answer */
    

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

typedef struct lin_eqn * lin_eqn_ptr;

static inline RR
tabs(RR t)
{
    if (t < 0)
        return -t;
    return t;
}

lin_eqn_ptr
lin_eqn_alloc(size_t n);

int
lin_bck_sweep(lin_eqn_ptr eq);

int
lin_fwd_sweep(lin_eqn_ptr eq);

int
linsolve(lin_eqn_ptr eq);

/* TRIDIAGONAL */
/* a_i x_{i-1} + b_i x_i + c_i x_{i+1} = f_i */
struct tridiag_eqn
{
    RR *a, *b, *c, *f;
    RR *x;
    size_t n;
};

typedef struct tridiag_eqn * tridiag_eqn_ptr;

tridiag_eqn_ptr
tridiag_eqn_alloc(size_t n);

int
tridiag_eqn_solve(tridiag_eqn_ptr eq);

void
tridiag_eqn_init(tridiag_eqn_ptr eq);

/* disturbed tridiagonal system
 * derived by adding $k$ and $(k+1)$th columns
 */
struct tridiag_twocol
{
    RR *a, *b, *c, *p, *q, *f;
    RR *x;
    uint32_t k1, k2;
    size_t n;
    bool free_self;
};

typedef struct tridiag_twocol * tridiag_twocol_ptr;

tridiag_twocol_ptr
tridiag_twocol_alloc(size_t n);

void
tridiag_twocol_free(tridiag_twocol_ptr eq);

void
tridiag_twocol_solve(tridiag_twocol_ptr eq);

void
tridiag_twocol_compute_f(tridiag_twocol_ptr eq, RR *x, RR *f);

tridiag_twocol_ptr
tridiag_twocol_copyof(tridiag_twocol_ptr eq);

#endif
