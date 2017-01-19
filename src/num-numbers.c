#include <stdio.h>
#include <err.h>
#include <vsu/num.h>

void inline static fskipspaces(FILE *f)
{
    int c;
    do
    {
        c = fgetc_unlocked(f);
    } while(isspace(c));
    ungetc(c, f);
}

int num_fnextn_unlocked(FILE *f, NN *n)
{
    *n = 0;
    int d;
    fskipspaces(f);
    d = fgetc_unlocked(f);
    if (!isdigit(d))
        return NUM_INVALID_INPUT;
    do
    {
        *n = (*n * 10) + (d - '0');
        d = fgetc_unlocked(f);
    }
    while(isdigit(d));
    ungetc(d, f);
    return 0;
}

int num_fnextz_unlocked(FILE *f, ZZ *z)
{
    *z = 1;
    int d;
    do d = fgetc_unlocked(f); while (isspace(d));
    if (d == '-') *z = -1;
    else ungetc(d, f);
    NN n;
    int ret = num_fnextn_unlocked(f, &n);
    if (ret == 0)
        *z *= n;
    return ret;
}

enum num_fnextr_state
{
    INIT, INT, FRAC, DONE
};

static inline bool isdelimiter(char c)
{
    return c == ',' || isspace(c);
}

int num_fnextr_unlocked(FILE *f, RR *res)
{
    RR r = 0;
    RR sign = 1.0;
    int d;
    enum num_fnextr_state state = INIT;

    while (state != DONE)
    {
        d = fgetc_unlocked(f);
        switch (state)
        {
            case INIT:
                if (isdelimiter(d)) {
                }
                else if (d == '+')
                {
                    state = INT;
                }
                else if (d == '-')
                {
                    sign = -1.0;
                    state = INT;
                }
                else if (d == '.')
                    state = FRAC;
                else if (isdigit(d))
                {
                    state = INT;
                    ungetc(d, f);
                }
                else return NUM_INVALID_INPUT;
                continue;
            case INT:
                if (isdigit(d))
                    r = 10*r + sign*(d-'0');
                else if (d == '.')
                    state = FRAC;
                else state = DONE;
                continue;
            case FRAC:
                if (isdigit(d))
                {
                    sign /= 10;
                    r = r + sign*(d - '0');
                } else state = DONE;
                continue;
        }
    }
    ungetc(d, f);
    *res = r;
    return 0;
}

int num_fnextc_unlocked(FILE *f, CC *c)
{
    *c = 0;
    bool done = false;
    int sign = 0;
    while (!done)
    {
        RR r;
        int ret = num_fnextr_unlocked(f, &r);
        if (ret != 0) return ret;
        int d = fgetc_unlocked(f);
        if (tolower(d) == 'i') {
            *c += sign*r*I;
            r = 0;
            d = fgetc_unlocked(f);
        }
        if (d == '+') {
            *c += r;
            sign = 1;
        } else if (d == '-') {
            *c += r+0*I;
            sign = -1;
        } else {
            done = true;
            ungetc(d, f);
        }
    }
    return 0;
}

#define _WRAP(TYP, FNEXT, FGET, SGET, SNEXT, F_UNLOCKED) \
    int FNEXT(FILE *f, TYP *x) \
{ \
    flockfile(f); \
    int ret = F_UNLOCKED(f, x); \
    funlockfile(f); \
    return ret; \
} \
    TYP FGET(FILE *f) \
{ \
    TYP x; \
    int ret = FNEXT(f, &x); \
    if (ret != 0) \
        err(ret, "error reading number"); \
    return x; \
} \
    int SNEXT(char *str, size_t slen, TYP *result) \
{ \
    FILE *f = fmemopen(str, slen, "r"); \
    int ret = F_UNLOCKED(f, result); \
    fclose(f); \
    return ret; \
} \
    TYP SGET(char *str, size_t slen) \
{ \
    TYP x; \
    int ret = SNEXT(str, slen, &x); \
    if (ret != 0) err(ret, "error reading number"); \
    return x; \
}

_WRAP(NN, num_fnextn, num_fgetn, num_sgetn, num_snextn, num_fnextn_unlocked);
_WRAP(ZZ, num_fnextz, num_fgetz, num_sgetz, num_snextz, num_fnextz_unlocked);
_WRAP(RR, num_fnextr, num_fgetr, num_sgetr, num_snextr, num_fnextr_unlocked);
_WRAP(CC, num_fnextc, num_fgetc, num_sgetc, num_snextc, num_fnextc_unlocked);
#undef _WRAP

bool inline static num_eqr0(RR a, RR eps)
{
    return (a > 0) ? (a < eps) : (-a < eps);
}
bool num_eqr(RR a, RR b, RR eps)
{
    a -= b;
    return num_eqr0(a, eps);
}
bool num_eqc(CC a, CC b, RR eps)
{
    a -= b;
    return num_eqr0(REAL(a), eps)
        && num_eqr0(IMAG(a), eps);
}

void
num_fill_vals_cc(
        CC_TO_CC f,
        CC *dom,
        size_t n,
        CC *out)
{
    for (int i = 0; i < (n); ++i)
        out[i] = (*f)(dom[i]);
}

void
num_fill_vals_rr(
        RR_TO_RR f,
        RR *dom,
        size_t n,
        RR *out)
{
    for (int i = 0; i < (n); ++i)
        out[i] = (*f)(dom[i]);
}

void
num_fill_vals_frproj(
        CC_TO_CC f,
        RR *dom,
        size_t n,
        RR *out)
{
    for (int i = 0; i < (n); ++i)
        out[i] = REAL((*f)(dom[i]));
}

#define NUM_GRID_EQUIDIST(FNAME, TYP) \
    void \
FNAME( \
        TYP a, \
        TYP b, \
        size_t n, \
        TYP *row) \
{ \
    if ((n) > 0)  \
    { \
        for (int k = 1; k < n; ++k) \
        { \
            row[k] = (TYP)(a + (b - a)*k/(n-1)); \
        } \
        row[0] = (TYP) a;  \
        row[n-1] = (TYP) b;  \
    } \
    else \
    { \
        row = 0;  \
    } \
}
    

NUM_GRID_EQUIDIST(num_grid_eqdst_cc, CC);
NUM_GRID_EQUIDIST(num_grid_eqdst_rr, RR);
#undef NUM_GRID_EQUIDIST

static char *FORMAT_REAL = "%12.4Lf";
static char *FORMAT_IMAG = "%12.4Lfi";
static char *FORMAT_COMPLEX_SUM = "%5.2Lf+%05.2Lfi";
static char *FORMAT_COMPLEX_DIF = "%6.2Lf%05.2Lfi";

int num_fput_cc(CC num, FILE *out)
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
