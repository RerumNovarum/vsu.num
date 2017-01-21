#include <vsu/num.h>
#include <vsu/linsolve.h>

int ret = 0;

int
test(tridiag_twocol_ptr eq)
{
    uint32_t n = eq->n;
    uint32_t k1 = eq->k1, k2 = eq->k2;
    tridiag_twocol_ptr cpy  = tridiag_twocol_copyof(eq);
    tridiag_twocol_solve(eq);
    if (eq->x == NULL) {
        ret |= 1;
        free(cpy);
        return 1;
    }
    tridiag_twocol_compute_f(cpy, eq->x, eq->a); /* using eq->a as a temporary storage */

    printf("test(%zu, %u, %u):\n", n, k1, k2);
    printf("  answer: (");
    for (int i = 0; i < n; ++i) {
        printf("%Lf, ", eq->a[i]);
    }
    printf(")\n  nullity vector: (");
    RR dev = 0;
    for (int i = 0; i < n; ++i) {
        RR newdev = fabsl(eq->a[i] - cpy->f[i]);
        printf("%Lf, ", newdev);
        if (dev < newdev) dev = newdev;
    }
    printf(")\n  error=%Lf\n", dev);
    free(cpy);
    return 0;
}

int tc0(size_t n, uint32_t k) {
    tridiag_twocol_ptr eq = tridiag_twocol_alloc(n);
    for (int i = 0; i < n; ++i)
        eq->p[i] = 1;
    for (int i = 0; i < n; ++i)
        eq->q[i] = 1;
    for (int i = 0; i < n; ++i)
        eq->a[i] = 1;
    for (int i = 0; i < n; ++i)
        eq->c[i] = 1;
    for (int i = 0; i < n; ++i)
        eq->b[i] = 2;
    int k1 = k-1;
    int k2 = k;
    eq->a[0] = eq->c[n-1] = 0;
    eq->k1 = k1;
    eq->k2 = k2;
    if (k1 > 0) eq->p[k1-1] = eq->c[k1-1];
    eq->p[k1] = eq->b[k1];
    if ((k1+1) < n) eq->p[k1+1] = eq->a[k1+1];
    if (k2 > 0) eq->q[k2-1] = eq->c[k2-1];
    eq->q[k2] = eq->b[k2];
    if ((k2+1) < n) eq->p[k2+1] = eq->a[k2+1];

    RR u[n];
    for (int i = 0; i < n; ++i) u[i] = 1.0;
    tridiag_twocol_compute_f(eq, u, eq->f);

    if (test(eq) != 0)
        printf("tc0(%zu, %u) failed\n", n, k);
}

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
    int k1 = k-1;
    int k2 = k;
    eq->a[0] = eq->c[n-1] = 0;
    eq->k1 = k1;
    eq->k2 = k2;
    if (k1 > 0) eq->p[k1-1] = eq->c[k1-1];
    eq->p[k1] = eq->b[k1];
    if ((k1+1) < n) eq->p[k1+1] = eq->a[k1+1];
    if (k2 > 0) eq->q[k2-1] = eq->c[k2-1];
    eq->q[k2] = eq->b[k2];
    if ((k2+1) < n) eq->p[k2+1] = eq->a[k2+1];

    RR u[n];
    for (int i = 0; i < n; ++i) u[i] = 1.0;
    tridiag_twocol_compute_f(eq, u, eq->f);

    if (test(eq) != 0)
        printf("tc0(%zu, %u) failed\n", n, k);
}

int main()
{
    tc0(6,2);
    tc0(8,2);
    tc0(16,2);
    tc1(6, 2);
    tc1(8, 2);
    tc1(16, 2);
    return ret;
}
