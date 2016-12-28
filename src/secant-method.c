#include <vsu/num.h>

RR root_secant_method(RR_TO_RR f, RR x0, RR x1)
{
    RR last_er = 1e9;
    RR er;
    while ((er = fabsl(x0-x1)) < last_er)
    {
        last_er = er;
        RR t = x1;
        x1 = x0 - (*f)(x0) * (x1-x0) / ((*f)(x1) - (*f)(x0));
        x0 = t;
    }
    return x1;
}
