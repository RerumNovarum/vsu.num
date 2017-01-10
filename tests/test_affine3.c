#include <string.h>
#include <vsu/num.h>

int ret = 0;

static int
test_composition2(char *tcase_name, struct affine3 A, struct affine3 B, RR X, RR Y, RR Z)
{
    struct affine3 T = affine3mul(A, B);

    struct affine3 Tv = affine3mul_n(2, A, B);
    if (T.a11 != Tv.a11)
    {
        ret = 1;
        printf("variadic composition (%s) failed; a11=%Lf, expected %Lf\n", tcase_name, Tv.a11, T.a11);
    }
    if (T.a12 != Tv.a12)
    {
        ret = 1;
        printf("variadic composition (%s) failed; a12=%Lf, expected %Lf\n", tcase_name, Tv.a12, T.a12);
    }
    if (T.a13 != Tv.a13)
    {
        ret = 1;
        printf("variadic composition (%s) failed; a13=%Lf, expected %Lf\n", tcase_name, Tv.a13, T.a13);
    }
    if (T.a21 != Tv.a21)
    {
        ret = 1;
        printf("variadic composition (%s) failed; a21=%Lf, expected %Lf\n", tcase_name, Tv.a21, T.a21);
    }
    if (T.a22 != Tv.a22)
    {
        ret = 1;
        printf("variadic composition (%s) failed; a22=%Lf, expected %Lf\n", tcase_name, Tv.a22, T.a22);
    }
    if (T.a23 != Tv.a23)
    {
        ret = 1;
        printf("variadic composition (%s) failed; a23=%Lf, expected %Lf\n", tcase_name, Tv.a23, T.a23);
    }
    if (T.a31 != Tv.a31)
    {
        ret = 1;
        printf("variadic composition (%s) failed; a31=%Lf, expected %Lf\n", tcase_name, Tv.a31, T.a31);
    }
    if (T.a32 != Tv.a32)
    {
        ret = 1;
        printf("variadic composition (%s) failed; a32=%Lf, expected %Lf\n", tcase_name, Tv.a32, T.a32);
    }
    if (T.a33 != Tv.a33)
    {
        ret = 1;
        printf("variadic composition (%s) failed; a33=%Lf, expected %Lf\n", tcase_name, Tv.a33, T.a33);
    }
    if (T.b1 != Tv.b1)
    {
        ret = 1;
        printf("variadic composition (%s) failed; b1=%Lf, expected %Lf\n", tcase_name, Tv.b1, T.b1);
    }
    if (T.b2 != Tv.b2)
    {
        ret = 1;
        printf("variadic composition (%s) failed; b2=%Lf, expected %Lf\n", tcase_name, Tv.b2, T.b2);
    }
    if (T.b3 != Tv.b3)
    {
        ret = 1;
        printf("variadic composition (%s) failed; b3=%Lf, expected %Lf\n", tcase_name, Tv.b3, T.b3);
    }
    if (T.alpha != Tv.alpha)
    {
        ret = 1;
        printf("variadic composition (%s) failed; alpha=%Lf, expected %Lf\n", tcase_name, Tv.alpha, T.alpha);
    }


    RR x0 = X, y0 = Y, z0 = Z;
    affine3apply_rr(B, &x0, &y0, &z0);
    affine3apply_rr(A, &x0, &y0, &z0);
    RR x1 = X, y1 = Y, z1 = Z;
    affine3apply_rr(T, &x1, &y1, &z1);
    if (fabsl(x0 - x1) > 1e-9)
    {
        ret = 1;
        printf("composition (%s) failed on vector (%Lf, %Lf, %Lf); x diverges\n", tcase_name, X, Y, Z);
    }
    if (fabsl(y0 - y1) > 1e-9)
    {
        ret = 1;
        printf("composition (%s) failed on vector (%Lf, %Lf, %Lf); y diverges\n", tcase_name, X, Y, Z);
    }
    return 0;
}

int main()
{
    struct affine3 rotz_pi_2 = affine3rotz(M_PI_2);
    struct affine3 tr_1_2 = affine3tr(1, 2, 3);
    struct affine3 inv_tr_1_2 = affine3tr(-1, -2, -3);
    struct affine3 scale_2_2 = affine3scale(2, 2, 2);
    struct affine3 mirror_x = affine3scale(-1, 1, -1);

    RR x = 1, y = 1, z = 1;
    affine3apply_rr(rotz_pi_2, &x, &y, &z);
    if (fabsl(x - (-1)) > 1e-9)
    {
        printf("rot pi/2 failed: x=%Lf\n", x);
        ret = 1;
    }
    if (fabsl(y - (1)) > 1e-9)
    {
        printf("rot pi/2 failed: y=%Lf\n", y);
        ret = 1;
    }

    if (ret == 0)
        puts("rot pi/2 seems to be ok\n");

    test_composition2("rotz[pi/2] tr[1, 2, 3]", rotz_pi_2, tr_1_2, 1, 2, 3);
    test_composition2("rotz[pi/2] tr[1, -2, 3]", rotz_pi_2, tr_1_2, 1, -2, 3);
    test_composition2("rotz[pi/2] tr[-1, 2, -3]", rotz_pi_2, tr_1_2, -1, 2, -3);
    test_composition2("rotz[pi/2] tr[-1, -2, -3]", rotz_pi_2, tr_1_2, -1, -2, -3);

    /* TODO: make some real tests */
    if (ret == 0)
        puts("variadic affine3mul_n seems to be ok\n");

    puts("TODO: make some real tests you lazy moron\n");
    
    return ret;
}
