#include <vsu/num.h>

void inline static
_num_fput_table_cc(
        const CC *x,
        const CC *y,
        size_t n,
        bool pretty,
        FILE *out)
{
    flockfile(out);
    if (pretty)
    {
        fprintf(out, "n=%d\n", n);
        for (int i = 0; i < n; ++i)
        {
            num_fput_cc(x[i], out);
            fputs_unlocked("\t|\t", out);
            num_fput_cc(y[i], out);
            fputc_unlocked('\n', out);
        }
    } else
    {
        fprintf(out, "%zu\n", n);
        for (int i = 0; i < n; ++i)
        {
            num_fput_cc(x[i], out);
            if (i != n)
                fputs_unlocked(" ", out);
        }
        fputc_unlocked('\n', out);
        for (int i = 0; i < n; ++i)
        {
            num_fput_cc(y[i], out);
            if (i != n)
                fputs_unlocked(" ", out);
        }
    }
    funlockfile(out);
}

void
num_fput_table_cc(
        const CC *x,
        const CC *y,
        size_t n,
        FILE *out)
{
    _num_fput_table_cc(x, y, n, false, out);
}
void
num_fputpretty_table_cc(
        const CC *x,
        const CC *y,
        size_t n,
        FILE *out)
{
    _num_fput_table_cc(x, y, n, true, out);
}

#define _NUM_TABLE_EQ(TYP, TBLTYP, FEQNAME, FEQTNAME) \
bool FEQNAME( \
        const TYP *x1, \
        const TYP *y1, \
        size_t n1, \
        const TYP *x2, \
        const TYP *y2, \
        size_t n2) \
{ \
    if (n1 != n2) return false; \
    for (int i = 0; i < n1; ++i) \
    { \
        if (x1[i] != x2[i]) return false; \
        if (y1[i] != y2[i]) return false; \
    } \
    return true; \
} \
bool \
FEQTNAME( \
        TBLTYP *t1, \
        TBLTYP *t2) \
{ \
    return FEQNAME( \
            t1->x, t1->y, t1->pt_no, \
            t2->x, t2->y, t2->pt_no); \
}

_NUM_TABLE_EQ(CC, TABLE_CC, num_table_cc_eq, num_table_cc_eqt);
#undef _NUM_TABLE_EQ

RR num_max_deviation_cc(CC *X, CC *Y, size_t pt_no)
{
    RR max = -INFINITY;
    for (int i = 0; i < pt_no; ++i)
    {
        RR d = CABS(X[i] - Y[i]);
        max = fmaxl(max, d);
    }
    return max;
}
CC
NUM_ZERO_CC(CC x) { return 0; }
CC
NUM_IDENTITY_CC(CC x) { return x; }
CC
NUM_SQR_CC(CC x) { return x*x; }
