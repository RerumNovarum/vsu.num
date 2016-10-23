#ifndef _VSU_NUM_H_
#define _VSU_NUM_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <ctype.h>
#include <complex.h>
#include <math.h>


#ifdef DEBUG
#define DPRINTF(...) \
    do { fprintf(stderr, __VA_ARGS__); } while(0)
#else
#define DPRINTF(...)
#endif

typedef uint64_t NN;
typedef int64_t ZZ;
typedef long double RR;
typedef long double complex CC;


typedef NN (*NN_TO_NN)(NN);
typedef ZZ (*ZZ_TO_ZZ)(ZZ);
typedef RR (*RR_TO_RR)(RR);
typedef CC (*CC_TO_CC)(CC);

#define REAL(x) creall((x))
#define IMAG(x) cimagl((x))
#define CABS(x) cabsl((x))

typedef struct numpoint_r
{
    RR x, y;
} POINT_R;

typedef struct numpoint
{
    CC x, y;
} POINT;

NN num_fgetn(FILE *f);
ZZ num_fgetz(FILE *f);
RR num_fgetr(FILE *f);
CC num_fgetc(FILE *f);

NN num_sgetn(char *str, size_t sz);
ZZ num_sgetz(char *str, size_t sz);
RR num_sgetr(char *str, size_t sz);
CC num_sgetc(char *str, size_t sz);

int num_fnextn(FILE *f, NN *result);
int num_fnextz(FILE *f, ZZ *result);
int num_fnextr(FILE *f, RR *result);
int num_fnextc(FILE *f, CC *result);

int num_snextn(char *s, size_t slen, NN *result);
int num_snextz(char *s, size_t slen, ZZ *result);
int num_snextr(char *s, size_t slen, RR *result);
int num_snextc(char *s, size_t slen, CC *result);

int num_ntos(CC num, char **deststr, size_t *len);
int num_fputn(CC num, FILE *out);
bool num_eq(CC a, CC b, RR eps);

typedef struct num_vec
{
    CC *x;
    size_t n;
} VECTOR;

typedef struct num_table
{
    size_t pt_no;
    CC *x;
    CC *y;
} TABLE;

void 
num_fnext_table_cc(FILE *in, CC **x, CC **y, size_t *n);
void 
num_fnext_table_rr(FILE *in, RR **x, RR **y, size_t *n);

void
num_grid_gen_eqdst_cc(
        CC a,
        CC b,
        size_t n,
        CC **outgrid);
void
num_grid_gen_eqdst_rr(
        RR a,
        RR b,
        size_t n,
        RR **outgrid);
void
num_fill_vals(
        CC_TO_CC f,
        CC *domain,
        size_t n,
        CC **out);
void
num_fill_vals_r(
        RR_TO_RR f,
        RR *domain,
        size_t n,
        RR **out);
void
num_fill_vals_frproj(
        CC_TO_CC f,
        RR *domain,
        size_t n,
        RR **out);
void
num_table_eqdst_cc(
        CC_TO_CC f,
        CC a,
        CC b,
        size_t n,
        CC **x,
        CC **y);
void
num_fput_table(
        const CC *x,
        const CC *y,
        size_t n,
        bool pretty,
        FILE *out);
bool
num_table_eq(
        const CC *x1,
        const CC *y1,
        size_t n1,
        const CC *x2,
        const CC *y2,
        size_t n2);
bool
num_table_eqt(
        TABLE *t1,
        TABLE *t2);
RR
num_max_deviation(CC *X, CC *Y, size_t pt_no);
CC
NUM_ZERO(CC x);
CC
NUM_IDENTITY(CC x);
CC
NUM_SQR(CC x);
#endif
