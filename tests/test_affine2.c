#include <string.h>
#include <vsu/num.h>

int ret = 0;

static int
test_composition2(char *tcase_name, struct affine2 A, struct affine2 B, RR X, RR Y)
{
    struct affine2 T = affine2mul(A, B);

    struct affine2 Tv = affine2mul_n(2, A, B);
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
    if (T.alpha != Tv.alpha)
    {
        ret = 1;
        printf("variadic composition (%s) failed; alpha=%Lf, expected %Lf\n", tcase_name, Tv.alpha, T.alpha);
    }

    RR x0 = X, y0 = Y;
    affine2apply_rr(B, &x0, &y0);
    affine2apply_rr(A, &x0, &y0);
    RR x1 = X, y1 = Y;
    affine2apply_rr(T, &x1, &y1);
    if (fabsl(x0 - x1) > 1e-9)
    {
        ret = 1;
        printf("composition (%s) failed on vector (%Lf, %Lf); x diverges\n", tcase_name, X, Y);
    }
    if (fabsl(y0 - y1) > 1e-9)
    {
        ret = 1;
        printf("composition (%s) failed on vector (%Lf, %Lf); y diverges\n", tcase_name, X, Y);
    }
    return 0;
}

int main()
{
    struct affine2 rot_pi_2 = affine2rot(M_PI_2);
    struct affine2 tr_1_2 = affine2tr(1, 2);
    struct affine2 inv_tr_1_2 = affine2tr(-1, -2);
    struct affine2 scale_2_2 = affine2scale(2, 2);
    struct affine2 mirror_x = affine2scale(-1, 1);

    RR x = 1, y = 1;
    affine2apply_rr(rot_pi_2, &x, &y);
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

    test_composition2("rot[pi/2] tr[1, 2]", rot_pi_2, tr_1_2, 1, 2);
    test_composition2("rot[pi/2] tr[1, -2]", rot_pi_2, tr_1_2, 1, -2);
    test_composition2("rot[pi/2] tr[-1, 2]", rot_pi_2, tr_1_2, -1, 2);
    test_composition2("rot[pi/2] tr[-1, -2]", rot_pi_2, tr_1_2, -1, -2);

    /* TODO: make some real tests */
    if (ret == 0)
        puts("variadic affine2mul_n seems to be ok\n");

    puts("TODO: make some real tests you lazy moron\n");
    
    return ret;
}
