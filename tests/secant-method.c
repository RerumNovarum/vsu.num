#include <vsu/num.h>

int ret = 0;

RR f1(RR x)
{
    return 2*x*x + x;
}

int main()
{
    RR eps = 1e-8;

    RR x_0;
    RR er;
    x_0 = root_secant_method(f1, 0.5, 4, eps);
    er = f1(x_0) - 0;
    printf("2x^2 + x + 3, x_0=%Lf, x_1=%Lf: x_root=%Lf error=%Lf\n", 0.5l, 4.0l, x_0, er);

    x_0 = root_secant_method(log10l, 0.2, 2, eps);
    er = log10l(x_0) - 0;
    printf("log10, x_0=%Lf, x_1=%Lf: x_root=%Lf error=%Lf\n", 0.2l, 2.0l, x_0, er);
    return ret;
}
