#include "vsu/num.h"

static char *FORMAT_REAL = "%12.4Lf";
static char *FORMAT_IMAG = "%12.4Lfi";
static char *FORMAT_COMPLEX_SUM = "%5.2Lf+%05.2Lfi";
static char *FORMAT_COMPLEX_DIF = "%6.2Lf%05.2Lfi";

inline CC num_fgetn(FILE *in)
{
    CC x;
    num_fnextn(in, &x);
    return x;
}

inline CC num_sgetn(char *str, size_t sz)
{
    CC x;
    num_snextn(str, sz, &x);
    return x;
}

int num_snextn(char *str, size_t sz, CC *result)
{
    FILE *in = fmemopen(str, sz, "r");
    int i = num_fnextn(in, result);
    fclose(in);
    return i;
}

enum num_fnextn_state
{
    STATE_PREINIT, /* initial */
    STATE_INIT, /* begin reading another number */
    STATE_NINT, /* reading integral part */
    STATE_NFRA, /* reading fractional part */
    STATE_SAVR, /* save real part */
    STATE_SAVI, /* save imaginary part */
    STATE_DONE /* done */
};

static inline bool isdelimiter(char c)
{
    return c == ',' || isspace(c);
}

int num_fnextn(FILE *in, CC *result)
{
    /* parse complex numbers
     * encoded as a string of a strict form;
     * (originally)
     * three forms of a CC are allowed
     * (denoted as regex, as forward slashes suggest):
     * /-?[0-9]*\.[0-9]+/
     * /-?[0-9]*\.[0-9]+i/
     * /-?[0-9]*\.[0-9]+[+-][0-9]*\.[0-9]+i/
     *
     * e.g.:
     * 1
     * 2i
     * 1+2i
     *
     * at the moment these constrains are not properly checked
     * and long sums are allowed;
     *
     * TODO: constrains
     * TODO: do you ever need number of bytes read? just return a CC;
     */

    const RR BASE  = 10;

    enum num_fnextn_state state = STATE_PREINIT;
    RR order = 1.0;
    RR sign = 1.0;
    RR x = 0.0;
    int i = 0;

    flockfile(in);
    while (1)
    {
        int c = fgetc_unlocked(in);
        if (isdelimiter(c))
        {
            ++i;
        } else
        {
            ungetc(c, in);
            break;
        }
            
    }

    *result = 0;

    bool global_empty = true;
    bool cur_empty = true;
    while(state != STATE_DONE)
    {
        int c;
#ifdef DEBUG
        c = fgetc_unlocked(in);
        fprintf(stderr,
                "i=%d s=%d c='%c'=%d x=%Lf*(%Lf+%Lfi)\n",
                i, state, c, c,
                sign, REAL(x), IMAG(x));
        ungetc(c, in);
#endif
        switch (state)
        {
            case STATE_PREINIT:
                x = 0.0;
                order = 1.0;
                sign = 1.0;
                state = STATE_INIT;
                cur_empty = true;
                break;
            case STATE_INIT:
                c = fgetc_unlocked(in);
                i += 1;

                if (c == EOF || isdelimiter(c))
                {
                    if (global_empty)
                        return -1;
                    state = STATE_DONE;
                } else if (c == '+')
                {
                } else if (c == '-')
                {
                    sign = -sign;
                } else if (c == '.')
                {
                    state = STATE_NFRA;
                } else if (isdigit(c) || c == 'i')
                {
                    state = STATE_NINT;
                    ungetc(c, in);
                    i -= 1;
                } else
                {
                    return -1;
                }
                break;
            case STATE_NINT:
                c = fgetc_unlocked(in);
                i += 1;
                if (c == '+' || c == '-' || isdelimiter(c) || c == EOF)
                {
                    state = STATE_SAVR;
                    ungetc(c, in);
                    i -= 1;
                } else if (c == 'i')
                {
                    if (cur_empty) x = 1.0;
                    state = STATE_SAVI;
                } else if (c == '.')
                {
                    state = STATE_NFRA;
                } else if (isdigit(c))
                {
                    cur_empty = false;
                    RR digit = (CC_R) (c - '0');
                    x *= BASE;
                    x += digit;
                } else
                {
                    return -1;
                }
                break;
            case STATE_NFRA:
                c = fgetc_unlocked(in);
                i += 1;
                if (c == '+' || c == '-' || isdelimiter(c) || c == EOF)
                {
                    state = STATE_SAVR;
                    ungetc(c, in);
                    i -= 1;
                } else if (c == 'i')
                {
                    if (cur_empty)
                        return -1;
                    state = STATE_SAVI;
                } else if (isdigit(c))
                {
                    cur_empty = false;
                    RR digit = (CC_R) (c - '0');
                    order /= BASE;
                    x += order * digit;
                } else
                {
                    return -1;
                }
                break;
            case STATE_SAVR:
                (*result) += sign*x;
                state = STATE_PREINIT;
                global_empty = false;
                break;
            case STATE_SAVI:
                (*result) += sign*x*I;
                state = STATE_PREINIT;
                global_empty = false;
                break;
        }
    }
    funlockfile(in);

    return i-1; /* -1 excludes closing EOF or space */
}

int num_fputn(CC num, FILE *out)
{
    RR real = REAL(num);
    RR imag = IMAG(num);
    if (imag == 0)
    {
        return fprintf(out, FORMAT_REAL, real);
    }
    else if (real == 0)
    {
        return fprintf(out, FORMAT_IMAG, imag);
    }
    else if (imag < 0)
    {
        return fprintf(out, FORMAT_COMPLEX_DIF, real, imag);
    }
    else
    {
        return fprintf(out, FORMAT_COMPLEX_SUM, real, imag);
    }
}

int num_ntos(CC num, char **deststr, size_t *len)
{
    FILE *out = open_memstream(deststr, len);
    int bytes_no = num_fputn(num, out);
    fclose(out);
    return bytes_no;
}

bool num_eq(CC a, NUMBER b, RR eps)
{
    CC d = b-a;
    /* return cabsl(d) < eps; */
    if (REAL(d) >= eps) return false;
    if (IMAG(d) >= eps) return false;
    return true;
}

void 
num_fnext_table_cc(FILE *in, CC **x, NUMBER **y, size_t *n)
{
    fscanf(in, "%zu", n);
    *x = malloc((*n)*sizeof(CC));
    *y = malloc((*n)*sizeof(CC));

    for (int i = 0; i < *n; ++i) {
        num_fnextn(in, (*x)+i);
    }
    for (int i = 0; i < *n; ++i) {
        num_fnextn(in, (*y)+i);
    }
}

void 
DELETEME(FILE *in, TABLE *t)
{
    num_fnext_table_cc(in, &(t->x), &(t->y), &(t->pt_no));
}

void
num_fill_vals(
        CC_TO_CC f,
        CC *dom,
        size_t n,
        CC **out)
{
    *out = malloc(n * sizeof(CC));
    for (int i = 0; i < (n); ++i)
    {
        (*out)[i] = (*f)(dom[i]);
    }
}

void
num_fill_vals_r(
        RR_TO_RR f,
        RR *dom,
        size_t n,
        RR **out)
{
    *out = malloc(n * sizeof(RR));
    for (int i = 0; i < (n); ++i)
    {
        (*out)[i] = (*f)(dom[i]);
    }
}

void
num_fill_vals_frproj(
        CC_TO_CC f,
        RR *dom,
        size_t n,
        RR **out)
{
    *out = malloc(n * sizeof(RR));
    for (int i = 0; i < (n); ++i)
    {
        (*out)[i] = REAL((*f)(dom[i]));
    }
}
#define NUM_GRID_EQUIDIST(a, b, n, row, type) \
{ \
    if ((n) > 0)  \
    {  \
        size_t bytes = (n) * sizeof(type);  \
        (row) = malloc(bytes);  \
        for (int k = 1; k < (n); ++k)  \
        {  \
            (row)[k] = (type)(a + (b - a)*k/(n-1));  \
        }  \
        (row)[0] = (type) a;  \
        (row)[n-1] = (type) b;  \
    }  \
    else  \
    {  \
        (row) = 0;  \
    } \
}
    

void
num_grid_gen_eqdst_cc(
        CC a,
        CC b,
        size_t n,
        CC **row)
{
    NUM_GRID_EQUIDIST(a, b, n, *row, CC);
}

void
num_grid_gen_eqdst_rr(
        RR a,
        RR b,
        size_t n,
        RR **row)
{
    NUM_GRID_EQUIDIST(a, b, n, *row, RR);
}

void
 num_table_equidist(
        CC (*f)(NUMBER x),
        CC a,
        CC b,
        size_t n,
        CC **x,
        CC **y)
{
    num_grid_gen_eqdst_cc(a, b, n, x);
    *y = malloc(n * sizeof(CC));
    num_fill_vals(f, *x, n, y);
}

void
num_table_equidistt(
        CC (*f)(NUMBER x),
        CC a,
        CC b,
        size_t n,
        TABLE *table)
{
    table->pt_no = n;
    num_table_equidist(f, a, b, n, &(table->x), &(table->y));
}

void
num_table_print(
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
            num_fputn(x[i], out);
            fputs_unlocked("\t|\t", out);
            num_fputn(y[i], out);
            fputc_unlocked('\n', out);
        }
    } else
    {
        fprintf(out, "%zu\n", n);
        for (int i = 0; i < n; ++i)
        {
            num_fputn(x[i], out);
            if (i != n)
                fputs_unlocked(" ", out);
        }
        fputc_unlocked('\n', out);
        for (int i = 0; i < n; ++i)
        {
            num_fputn(y[i], out);
            if (i != n)
                fputs_unlocked(" ", out);
        }
    }
    funlockfile(out);
}

void
num_table_printt(
        TABLE table,
        bool pretty,
        FILE *out)
{
    num_table_print(
            table.x, table.y, table.pt_no,
            pretty, out);
}

bool num_table_eq(
        const CC *x1,
        const CC *y1,
        size_t n1,
        const CC *x2,
        const CC *y2,
        size_t n2)
{
    if (n1 != n2) return false;
    for (int i = 0; i < n1; ++i)
    {
        if (x1[i] != x2[i]) return false;
        if (y1[i] != y2[i]) return false;
    }
    return true;
}

bool
num_table_eqt(
        TABLE *t1,
        TABLE *t2)
{
    return num_table_eq(
            t1->x, t1->y, t1->pt_no,
            t2->x, t2->y, t2->pt_no);
}

CC NUM_ZERO(NUMBER x)
{
    return (CC) 0;
}
CC NUM_IDENTITY(NUMBER x)
{
    return x;
}
CC NUM_SQR(NUMBER x)
{
    return x*x;
}

RR num_max_deviation(CC *X, NUMBER *Y, size_t pt_no)
{
    RR max = -INFINITY;
    for (int i = 0; i < pt_no; ++i)
    {
        RR d = ABS(X[i] - Y[i]);
        max = fmaxl(max, d);
    }
    return max;
}
