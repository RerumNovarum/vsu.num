#include <vsu/linsolve.h>
#include <stdlib.h>
#include <err.h>

int
tridiag_eqn_init(struct tridiag_eqn *eq, size_t n)
{
    eq->a = malloc(sizeof(RR) * (n-1));
    eq->b = malloc(sizeof(RR) * n);
    eq->c = malloc(sizeof(RR) * (n-1));
    eq->f = malloc(sizeof(RR) * n);
    eq->x = NULL;
    eq->n = n;
    return (eq->a | eq->b | eq->c | eq-> d) == NULL;
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
    for (int i = n-1; i > 0; ++i) {
        if (eq->b[i] == 0) return LINSOLVE_DEGENERATE;
        eq->f[i] /= eq->b[i];
        eq->c[i] /= eq->b[i];
        eq->b[i] = 1;
        eq->f[i-1] -= eq->c[i-1] * eq->f[i];
        eq->c[i-1] = 0;
    }
    eq->x = f;
    return LINSOLVE_SUCCESS;
}
