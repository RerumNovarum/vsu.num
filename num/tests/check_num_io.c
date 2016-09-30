#include <stdlib.h>
#include <check.h>

#include "vsu/num.h"

static char * sample_encoded_nums[] =
{
    "1",            /* */
    "2.0i",         /* */
    "1.0+2.0i",     /* */
    "0.25",         /* */
    "-0.25",        /* */
    "1-i",          /* */
    "-i",           /* */
    "-2i",          /* */
    "\n2",          /* nextnum should skip preceding whitespaces */
    "\n2 2i"        /* whitespace after nonempty number terminates parsing */
};
static NUMBER sample_expected_nums[] =
{
    (NUMBER) 1.0 + 0.0*I,
    (NUMBER) 2.0*I,
    (NUMBER) 1.0 + 2.0*I,
    (NUMBER) .25,
    (NUMBER) -.25,
    (NUMBER) 1.0 - I,
    (NUMBER) -I,
    (NUMBER) -2*I,
    (NUMBER) 2,
    (NUMBER) 2
};

START_TEST(test_num_fnextn_valid_input)
{
    size_t n = sizeof(sample_encoded_nums) / sizeof(*sample_encoded_nums);
    for (int i = 0; i < n; ++i)
    {
        NUMBER n;
        ck_assert_int_gt(num_snextn(sample_encoded_nums[i], strlen(sample_encoded_nums[i]), &n), -1);
        if (n != sample_expected_nums[i])
        {
            ck_abort_msg(
                    "got %Lf+%Lfi; expected %Lf+%Lfi",
                    REAL(n),
                    IMAG(n),
                    REAL(sample_expected_nums[i]),
                    IMAG(sample_expected_nums[i]));
        }
    }
}
END_TEST

START_TEST(test_num_fnextn)
{
    // TODO: use different samples
    size_t n = sizeof(sample_encoded_nums) / sizeof(*sample_encoded_nums);
    for (int i = 0; i < n; ++i)
    {
        char *enc;
        size_t len;
        num_ntos(sample_expected_nums[i], &enc, &len);
        ck_assert_str_eq(enc, sample_encoded_nums[i]);
        free(enc);
    }
}
END_TEST

START_TEST(test_num_readtable)
{
    char *encoded[] = 
    {
        "0",
        "1\n3.25\ni",
        "2\n1 2\n3 4"
    };
    NUMBER x1[] = { 3.25 };
    NUMBER y1[] = { I };
    NUMBER x2[] = {1, 2};
    NUMBER y2[] = {3, 4};
    TABLE tables[] =
    {
        {0, 0, 0},
        {1, x1, y1},
        {2, x2, y2}
    };
    size_t tc_no = sizeof(encoded) / sizeof(*encoded);

    for (int i = 0; i < tc_no; ++i) 
    {
        TABLE t;
        FILE *in = fmemopen(encoded[i], strlen(encoded[i]), "r");
        num_table_read(in, &t.x, &t.y, &t.pt_no);

        ck_assert_int_eq(tables[i].pt_no, t.pt_no);
        for (int j = 0; j < t.pt_no; ++j)
        {
            ck_assert_msg(tables[i].x[j] == t.x[j],
                    "%d'th table, %d'th x element: expected %Lf+%Lfi, got %Lf+%Lfi",
                    i, j,
                    REAL(tables[i].x[j]), IMAG(tables[i].x[j]),
                    REAL(t.x[j]), IMAG(t.x[j]));
            ck_assert_msg(tables[i].y[j] == t.y[j],
                    "%d'th table, %d'th y element: expected %Lf+%Lfi, got %Lf+%Lfi",
                    i, j,
                    REAL(tables[i].y[j]), IMAG(tables[i].y[j]),
                    REAL(t.y[j]), IMAG(t.y[j]));
        }
        fclose(in);
    }
}
END_TEST

START_TEST(test_num_table_eq)
{
    NUMBER x[] = { 0, 1, 2, 3 };
    NUMBER y[] = { 1, 1, 1, 1 };
    TABLE table = { 4, x, y };
    TABLE zero = { 0, 0, 0 };
    ck_assert(num_table_eqt(&table, &table));
    ck_assert(!num_table_eqt(&table, &zero));
}
END_TEST

NUMBER _RETURN_COMPLEX_ONE(NUMBER x)
{
    return (NUMBER) 1 + 0*I;
}

NUMBER _IDENTITY(NUMBER x)
{
    return x;
}

START_TEST(test_num_table_equidist)
{
    NUMBER x[] = { 0, 1, 2, 3 };
    NUMBER y[] = { 1, 1, 1, 1 };
    TABLE reference = { 4, x, y };
    TABLE t;
    num_table_equidistt(_RETURN_COMPLEX_ONE, x[0], x[3], 4, &t);
    bool eq = num_table_eqt(&reference, &t);
    ck_assert(eq);
}
END_TEST

START_TEST(test_num_table_print)
{
    NUMBER x[] = { 0, 1, 2, 3 };
    NUMBER y[] = { 3+I, 2-I, 1, 0 };
    size_t pt_no = 4;

    char *buff;
    size_t len;
    FILE *mem = open_memstream(&buff, &len);
    num_table_print(x, y, pt_no, false, mem);
    fclose(mem);

    mem = fmemopen(buff, len, "r");
    TABLE decoded;
    num_table_readt(mem, &decoded);
    fclose(mem);
    ck_assert(num_table_eq(x, y, pt_no, decoded.x, decoded.y, decoded.pt_no));
}
END_TEST

START_TEST(test_echo_sample_numtable)
{
    NUMBER x[] = { 0, 1, 2, 3 };
    NUMBER y[] = { 3+3*I, 2+2*I, 1+I, 0 };
    TABLE table = { 4, x, y };

    fputs("SAMPLE NUMGRID:\n", stderr);
    num_table_printt(table, false, stderr);
    fputs("\nPRETTYPRINTED:\n", stderr);
    num_table_printt(table, true, stderr);
}
END_TEST

Suite * num_io_suite(void)
{
    Suite *s;
    TCase *tc_core;

    s = suite_create("num_io");
    tc_core = tcase_create("Core");

    tcase_add_test(tc_core, test_num_table_eq);
    tcase_add_test(tc_core, test_num_fnextn_valid_input);
//    tcase_add_test(tc_core, test_num_fnextn);
    tcase_add_test(tc_core, test_num_table_equidist);
    tcase_add_test(tc_core, test_num_readtable);
    tcase_add_test(tc_core, test_num_table_print);
    tcase_add_test(tc_core, test_echo_sample_numtable);
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
