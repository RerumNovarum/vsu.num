#include <vsu/linsolve.h>
#include <err.h>

tridiag_twocol_ptr
tridiag_twocol_alloc(size_t n)
{
    size_t a_no = n,
           b_no = n,
           c_no = n, /* TODO: -1 */
           p_no = n,
           q_no = n, 
           f_no = n;
    void *buf = malloc(
            sizeof(struct tridiag_twocol) +
            sizeof(RR)*(a_no + b_no + c_no + p_no + q_no + f_no));
    if (buf == NULL) return NULL;
    tridiag_twocol_ptr eq = buf;
    eq->a = buf + sizeof(struct tridiag_twocol);
    eq->b = eq->a + a_no;
    eq->c = eq->b + b_no;
    eq->p = eq->c + c_no;
    eq->q = eq->p + p_no;
    eq->f = eq->q + q_no;
    eq->x = NULL;
    eq->n = n;
    eq->free_self = true;
    return eq;
}

void
tridiag_twocol_free(tridiag_twocol_ptr eq)
{
    if (eq->free_self) free(eq);
}

static inline int
_tridiag_twocol_is_consistent(tridiag_twocol_ptr eq)
{
    uint32_t k1, k2, n;
    k1 = eq->k1;
    k2 = eq->k2;
    n  = eq->n;
    return (eq->a[0] == 0) &&
        (eq->c[eq->n-1] == 0)  &&
        (eq->p[k1] == eq->b[k1]) &&
        ((k1 > 0) || eq->c[k1-1] == eq->p[k1-1]) &&
        ((k1 < n) || eq->p[k1] == eq->b[k1]) &&
        (((k1+1) < n) || eq->a[k1+1] == eq->p[k1+1]) &&
        ((k2 > 0) || eq->q[k2-1] == eq->c[k2-1]) &&
        ((k2 < n) || eq->q[k2] == eq->b[k2]) &&
        ((k2+1 < n) || eq->q[k2+1] == eq->a[k2+1]);
}

void
tridiag_twocol_solve(tridiag_twocol_ptr eq)
{
    if (!_tridiag_twocol_is_consistent(eq)) {
        err(1, "tridiag_twocol_solve: system is inconsistent");
    }

    uint32_t k1, k2, n;
    k1 = eq->k1;
    k2 = eq->k2;
    n  = eq->n;

    for (int i = 0; (i+1) < n; ++i) {
        if (eq->b[i] == 0) {
            err(1, "tridiag_twocol_solve: left-to-mid sweep: b[%d]=0", i);
        }
        /* a[i] = 0 */
        eq->c[i] /= eq->b[i];
        eq->p[i] /= eq->b[i];
        eq->q[i] /= eq->b[i];
        eq->f[i] /= eq->b[i];
        eq->b[i] = 1;
        eq->b[i+1] -= eq->a[i+1]*eq->c[i];
        eq->f[i+1] -= eq->a[i+1]*eq->f[i];
        eq->p[i+1] -= eq->a[i+1]*eq->p[i];
        eq->q[i+1] -= eq->a[i+1]*eq->q[i];
        if ((i+2) == k1)
            eq->c[i+1] = eq->p[i+1]; /* -= eq->a[i+1]*eq->p[i]; */
        if ((i+2) == k2)
            eq->c[i+1] = eq->q[i+1]; /* -= eq->a[i+1]*eq->q[i]; */
        eq->a[i+1] = 0;
    }
    for (int i = n-1; i > 0; --i) {
        if (eq->b[i] == 0) {
            err(1, "tridiag_twocol_solve: right-to-mid sweep: b[%d]=0", i);
        }
        /* c[i] = 0 */
        eq->a[i] /= eq->b[i];
        eq->p[i] /= eq->b[i];
        eq->q[i] /= eq->b[i];
        eq->f[i] /= eq->b[i];
        eq->b[i] = 1;
        eq->f[i-1] -= eq->c[i-1]*eq->f[i];
        eq->p[i-1] -= eq->c[i-1]*eq->p[i];
        eq->q[i-1] -= eq->c[i-1]*eq->q[i];
        if (i == (k2+1))
            eq->b[i-1] = eq->q[k2];
        if (i == (k1+1))
            eq->b[i-1] = eq->p[k1];
        if (i == (k2+2))
            eq->a[i-1] = eq->q[i-1];
        if (i == (k1+2))
            eq->a[i-1] = eq->p[i-1];
        eq->c[i-1] = 0;
    }
    eq->f[0] /= eq->b[0];
    eq->p[0] /= eq->b[0];
    eq->q[0] /= eq->b[0];
    eq->b[0] = 1.0;
    if (k1 < eq->n) {
        for (int i = 0; i < eq->n; ++i) {
            if (i == k1) continue;
            eq->f[i] -= eq->p[i]*eq->f[k1];
            eq->q[i] -= eq->p[i]*eq->q[k1];
            /* if (i == (k1+1)) eq->a[i] = 0; */
            eq->p[i] = 0;
        }
    }
    if (k2 < eq->n) {
        for (int i = 0; i < eq->n; ++i) {
            if (i == (k2)) continue;
            eq->f[i] -= eq->q[i]*eq->f[k2];
            eq->q[i] = 0;
        }
    }
    eq->x = eq->f;
    eq->f = NULL;
}

void
tridiag_twocol_compute_f(tridiag_twocol_ptr eq, RR *x, RR *f)
{
    uint32_t k1 = eq->k1;
    uint32_t k2 = eq->k2;
    uint32_t n  = eq->n;
    for (int i = 0; i < n; ++i) {
        f[i] = eq->b[i]*x[i];
        if (i > 0) f[i] += eq->a[i]*x[i-1];
        if ((i+1) < n) f[i] +=  eq->c[i]*x[i+1];
    }
    if (k1 < n) {
        for (int i = 0; i < n; ++i) {
            if ((i+1) < k1 || i > (k1+1)) {
                f[i] += eq->p[i]*x[k1];
            }
        }
    }
    if (k2 < n) {
        for (int i = 0; i < n; ++i) {
            if ((i+1) < k2 || i > k2+1) {
                f[i] += eq->q[i]*x[k2];
            }
        }
    }
}

tridiag_twocol_ptr
tridiag_twocol_copyof(tridiag_twocol_ptr eq)
{
    size_t n = eq->n;
    tridiag_twocol_ptr cpy = tridiag_twocol_alloc(n);
    RR *f = cpy->f, *eqf = eq->f;
    if (eqf == NULL) {
        cpy->f = NULL;
        cpy->x = f;
        f = eq->x;
    }
    for (int i = 0; i < n; ++i) {
        cpy->a[i] = eq->a[i];
        cpy->b[i] = eq->b[i];
        cpy->c[i] = eq->c[i];
        cpy->p[i] = eq->p[i];
        cpy->q[i] = eq->q[i];
        f[i] = eqf[i];
    }
    cpy->k1 = eq->k1;
    cpy->k2 = eq->k2;
    return cpy;
}
