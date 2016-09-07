#include <math.h>

#include "approx_newton.h"

NUMBER
newton_approx_equidist_kth_dd(
        int k,
        NUMBER_R h,
        NUMBER const *const vals,
        size_t pt_no)
{
    NUMBER dd = 0;
    for (int i = 0; i <= k; ++i)
    {
        NUMBER term = vals[i];
        for (int j = 2; j <= i; ++j)
            term /= j;
        for (int j = 2; j <= (k-i); ++j)
            term /= j;
        term /= pow(h, k);
        if ((k-i)%2 == 0)
            dd += term;
        else
            dd -= term;
    }
    return dd;
}

void
newton_approx_equidist(
        NUMBER_R a,
        NUMBER_R b,
        NUMBER const *const vals,
        size_t pt_no,
        NUMBER **out_dds)
{
    *out_dds = malloc(sizeof(NUMBER) * pt_no);
    NUMBER_R h = (b - a);
    if (pt_no > 1) h /= (pt_no-1);
    else h = 0;
    for (int k = 0; k < pt_no; ++k)
    {
        (*out_dds)[k] =
            newton_approx_equidist_kth_dd(k, h, vals, pt_no);
    }
}
