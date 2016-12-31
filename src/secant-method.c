#include <vsu/num.h>
#include <stdio.h>

RR root_secant_method(RR_TO_RR f, RR x0, RR x1, RR eps)
{
    /* current iteration's error */
    RR er;
    /* next approximation */
    RR x2;
    /* flag indicating that initial approximation (for which er < eps) has been found */
    int ini_found = 0;
    /* stopping criterion -- flag indicating that divergence started */
    int done = 0;
    /* flag indicating that er<eps on current iteration */
    int improved;
    while (!done)
    {
        /* next approximation is the intersection of secant line with the X axis */
        x2 = x0 -  (x1-x0) * (*f)(x0)/((*f)(x1) - (*f)(x0));
        /* an usual stopping criterion is |x2-x1| < eps */
        /* though we consider it to be the initial approximation criterion */
        er = fabsl(x2 - x1);
        improved = er < eps;
        /* once we had er<eps at least once, initial approximation is found */
        ini_found |= improved;

        /* if we have found initial approximation */
        /* but the error grew the last iteration */
        /* we conclude that divergence has begun, */
        /* stop iterations and return previous approximations */
        if (ini_found && !improved) {
            done = 1;
        } else {
            /* otherwise we accept the new approximation */
            x0 = x1;
            x1 = x2;
        }
        /* and if we already have the initial approximation */
        /* we raise (well, lower) the bar for acceptable errors */
        if (improved) {
            eps = er;
        }
    }
    return x1;
}
