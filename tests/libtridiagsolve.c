#include <vsu/num.h>
#include <vsu/linsolve.h>
#include <stdio.h>

int ret = 0;

int main()
{
    tridiag_eqn_ptr eq = tridiag_eqn_alloc(3);
    RR sol1[] = { 1, 2, 3 };
    eq->b[0] = 1; eq->c[0] = -1;
    eq->a[1] = 1; eq->b[1] = 2; eq->c[1] = 1;
    eq->a[2] = 1; eq->b[2] = -1;
    eq->f[0] = -1;
    eq->f[1] = 8;
    eq->f[2] = -1;
    if (tridiag_eqn_solve(eq) != LINSOLVE_SUCCESS) {
        puts("cannot solve eq1");
        ret |= 1;
    }
    for (int i = 0; i < sizeof(sol1)/sizeof(*sol1); ++i) {
        if (eq->x[i] != sol1[i]) {
            printf("eq1, got x[%d]=%Lf, expected %Lf\n",
                    i, eq->x[i], sol1[i]);
            ret |= 1;
        }
    }
    fflush(stdout);
    return ret;
}
