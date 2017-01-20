#include <vsu/linsolve.h>

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
    return eq;
}

int
tridiag_twocol_solve(tridiag_twocol_ptr eq)
{
    uint32_t k = eq->k;
    uint32_t n = eq->n;
    if (eq->a[0] != 0 ||
            eq->c[eq->n-1] != 0 ||
            eq->p[k] != eq->b[k] ||
            ((k > 0) && eq->p[k-1] != eq->p[k-1]) ||
            ((k+1 < n) && eq->q[k] != eq->c[k]) || /* if k+1>=n, there is no q */
            ((k+1 < n) && eq->q[k+1] != eq->b[k+1]) ||
            ((k+1 < n) && eq->p[k+1] != eq->a[k+1]) ||
            ((k+2 < n) && eq->q[k+2] != eq->a[k+2])) {
        return LINSOLVE_DEGENERATE;
    }
    for (int i = 0; i < eq->n-1; ++i) {
        if (eq->b[i] == 0) return LINSOLVE_DEGENERATE;
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
        eq->a[i+1] = 0;
    }
    for (int i = eq->n-1; i > 0; --i) {
        if (eq->b[i] == 0) return LINSOLVE_DEGENERATE;
        eq->c[i] /= eq->b[i];
        eq->p[i] /= eq->b[i];
        eq->q[i] /= eq->b[i];
        eq->f[i] /= eq->b[i];
        eq->b[i] = 1;
        eq->f[i-1] -= eq->c[i-1]*eq->f[i];
        eq->p[i-1] -= eq->c[i-1]*eq->p[i];
        eq->q[i-1] -= eq->c[i-1]*eq->q[i];
        eq->c[i-1] = 0;
    }
    eq->f[0] /= eq->b[0];
    eq->p[0] /= eq->b[0];
    eq->q[0] /= eq->b[0];
    eq->b[0] = 1.0;
    for (int i = 0; i < eq->n; ++i) {
        if (i == k) continue;
        eq->f[i] -= eq->p[i];
        eq->p[i] = 0;
    }
    k += 1;
    if (k < eq->n) {
        for (int i = 0; i < eq->n; ++i) {
            if (i == k) continue;
            eq->f[i] -= eq->q[i];
            eq->q[i] = 0;
        }
    }
    eq->x = eq->f;
    eq->f = NULL;
    return LINSOLVE_SUCCESS;
}
