#include <vsu/linsolve.h>
#include <stdlib.h>
#include <err.h>

void
tridiag_eqn_init(tridiag_eqn_ptr eq)
{
    eq-> x = NULL;
}

tridiag_eqn_ptr
tridiag_eqn_alloc(size_t n)
{
    size_t a_no = n /* TODO: -1 */,
           b_no = n,
           c_no = n /* TODO: -1 */,
           f_no = n;
    void *buf = malloc(
            sizeof(struct tridiag_eqn) +
            sizeof(RR)*(a_no + b_no + c_no + f_no));
    if (buf == NULL) return -1;
    tridiag_eqn_ptr eq = buf;
    eq->a = buf + sizeof(struct tridiag_eqn);
    eq->b = eq->a + a_no;
    eq->c = eq->b + b_no;
    eq->f = eq->c + c_no;
    eq->x = NULL;
    eq->n = n;
    return eq;
}

int
tridiag_eqn_solve(struct tridiag_eqn *eq)
{
    for (int i = 0; i < eq->n-1; ++i) {
        if (eq->b[i] == 0) return LINSOLVE_DEGENERATE;
        eq->c[i] /= eq->b[i];
        eq->f[i] /= eq->b[i];
        eq->b[i] = 1;
        eq->b[i+1] -= eq->a[i+1]*eq->c[i];
        eq->f[i+1] -= eq->a[i+1]*eq->f[i];
        eq->a[i+1] = 0;
    }
    for (int i = eq->n-1; i > 0; --i) {
        if (eq->b[i] == 0) return LINSOLVE_DEGENERATE;
        eq->f[i] /= eq->b[i];
        eq->c[i] /= eq->b[i];
        eq->b[i] = 1;
        eq->f[i-1] -= eq->c[i-1] * eq->f[i];
        eq->c[i-1] = 0;
    }
    eq->x = eq->f;
    eq->f = NULL;
    return LINSOLVE_SUCCESS;
}
