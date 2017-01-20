#include <vsu/num.h>
#include <vsu/linsolve.h>

int ret = 0;

int tc1(size_t n, uint32_t k)
{
    tridiag_twocol_ptr eq = tridiag_twocol_alloc(n);
    int j = 0;
    for (int i = 0; i < n; ++i)
        eq->p[i] = j++;
    for (int i = 0; i < n; ++i)
        eq->q[i] = j++;
    for (int i = 0; i < n; ++i)
        eq->a[i] = j++;
    for (int i = 0; i < n; ++i)
        eq->c[i] = j++;
    for (int i = 0; i < n; ++i)
        eq->b[i] = j++;
    eq->a[0] = eq->c[n-1] = 0;

    eq->k = k = k-1;
    eq->p[k] = eq->b[k];
    if (k < n) {
        eq->q[k] = eq->c[k];
        eq->q[k+1]   = eq->b[k+1];
        if ((k+2) < n)
            eq->q[k+2] = eq->a[k+2];
        eq->p[k+1] = eq->a[k+1];
    }

    for (int i = 0; i < n; ++i) {
        eq->f[i] = eq->a[i] + eq->b[i] + eq->c[i] +
            + eq->p[i] + eq->q[i];
    }

    if (tridiag_twocol_solve(eq) != 0) {
        ret |= 1;
        printf("failed to solve tc1(%zu, %u)\n", n, k+1);
    }
    RR dev = 0;
    for (int i = 0; i < n; ++i) {
        RR newdev = fabsl(eq->x[i] - 1.0);
        if (dev < newdev) dev = newdev;
    }

    printf("tc1(%zu,%u): error=%Lf\n", n, k+1, dev);
}

int main()
{
    tc1(6, 2);
    tc1(8, 2);
    tc1(16, 2);
    return ret;
}
