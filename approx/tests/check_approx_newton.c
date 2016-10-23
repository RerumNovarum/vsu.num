#include <stdlib.h>
#include <check.h>
#include <complex.h>
#include <math.h>

#include "vsu/num.h"
#include "vsu/approx_newton.h"

#define EPS 1e-8

CC _x3_3x2_1(NUMBER x)
{
    return (x + 3)*x*x + 1;
}
CC _log_re_x(NUMBER x)
{
    return logl(REAL(x));
}

START_TEST(test_newton_eval)
{
    CC dds[] = 
    {
        (CC) 1,
        (CC) 2,
        (CC) 3,
        (CC) 4*I
    };
    size_t n = 4;
    const size_t tc_no = 1;
    CC x[] =
    {
        1,
        0,
        4
    };
    CC X[][4] =
    {
        { 1, 2, 3, 4 }
    };
    CC y[] =
    {
        1,
        -19,
        4*I
    };

    for (int t = 0; t < tc_no; ++t)
    {
        CC f_x = newton_eval(x[t], X[t], dds, n);
        if (!(y[t] == f_x))
        {
            char *ys, *fxs;
            size_t len;
            num_ntos(y[t], &ys, &len);
            num_ntos(f_x, &fxs, &len);
            ck_abort_msg("newton_eval failed; expected %s, got %s\n",
                    ys, fxs);
            free(ys);
            free(fxs);
        }
    }
}
END_TEST

START_TEST(test_newton_kth_dd)
{
    CC X[] = {0, 1.0/3, 2.0/3, 1};
    CC Y[] = {0, 1.0/3, 2.0/3, 1};
    CC DD[] = { 0, 1, 0, 0 };
    RR h = 1.0/3;
    size_t pt_no = 4;
    size_t n = 4;

    for (int k = 0; k < 4; ++k)
    {
        CC dd = newton_approx_equidist_kth_dd(k, h, Y, n);
        if (!num_eq(DD[k], dd, EPS))
        {
            char *dds, *exps;
            size_t len;
            num_ntos(dd, &dds, &len);
            num_ntos(DD[k], &exps, &len);
            fprintf(stderr, "newton_kth_dd: got %s, expected %s\n", dds, exps);
            free(dds);
            free(exps);
            ck_abort();
        }
    }
}
END_TEST

START_TEST(test_approx_newton_equidist)
{
    TABLE table;

    CC (*ff[])(NUMBER)=
    {
        NUM_ZERO,
        NUM_IDENTITY,
        NUM_SQR,
        _x3_3x2_1
    };
    int abn[][3] =
    {
        {0, 1, 4},
        {0, 1, 4},
        {0, 2, 2},
        {0, 2, 5}
    };

    size_t tc_no = sizeof(ff)/sizeof(*ff);

    for (int t = 0; t < tc_no; ++t)
    {
        int a = abn[t][0];
        int b = abn[t][1];
        size_t n = abn[t][2];
        num_table_equidistt(ff[t], abn[t][0], abn[t][1], abn[t][2], &table);
        CC *dds;
        newton_approx_equidist(a, b, table.y, n, &dds);
        for (int j = 0; j < n; ++j)
        {
            CC y = newton_eval(table.x[j], table.x, dds, table.pt_no);
            CC d =  y - table.y[j];
            if (! (cabsl(d) < EPS))
            {
                char *sgot, *sexp;
                size_t len;
                num_ntos(y, &sgot, &len);
                num_ntos(table.y[j], &sexp, &len);
                ck_abort_msg("deviation must not exceed %e; y=%s expected %s\n",
                        EPS, sgot, sexp);
                free(sgot);
                free(sexp);
            }
        }
        free(table.x);
        free(table.y);
        free(dds);
    }
}
END_TEST

Suite * num_io_suite(void)
{
    Suite *s;
    TCase *tc_core;

    s = suite_create("num_approx_newton");
    tc_core = tcase_create("Core");

    tcase_add_test(tc_core, test_newton_eval);
    tcase_add_test(tc_core, test_newton_kth_dd);
    tcase_add_test(tc_core, test_approx_newton_equidist);
    suite_add_tcase(s, tc_core);

    return s;
}

int main()
{
    int num_failed;
    Suite *s;
    SRunner *sr;

    s = num_io_suite();
    sr = srunner_create(s);

    srunner_run_all(sr, CK_NORMAL);
    num_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return (num_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
