#ifndef _VSU_NUM_H_
#define _VSU_NUM_H_

#include <stdio.h> /* FILE &c */
#include <stdlib.h>
#include <inttypes.h>
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

#define NUM_PERIOD_C '.'

/* error codes */
#define NUM_INVALID_INPUT 1

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

NN num_fgetn(FILE *f);
ZZ num_fgetz(FILE *f);
RR num_fgetr(FILE *f);
CC num_fgetc(FILE *f);

NN num_sgetn(char *str, size_t sz);
ZZ num_sgetz(char *str, size_t sz);
RR num_sgetr(char *str, size_t sz);
CC num_sgetc(char *str, size_t sz);

int num_fnextn_unlocked(FILE *f, NN *result);
int num_fnextz_unlocked(FILE *f, ZZ *result);
int num_fnextr_unlocked(FILE *f, RR *result);
int num_fnextc_unlocked(FILE *f, CC *result);
int num_fnextn(FILE *f, NN *result);
int num_fnextz(FILE *f, ZZ *result);
int num_fnextr(FILE *f, RR *result);
int num_fnextc(FILE *f, CC *result);

int num_snextn(char *s, size_t slen, NN *result);
int num_snextz(char *s, size_t slen, ZZ *result);
int num_snextr(char *s, size_t slen, RR *result);
int num_snextc(char *s, size_t slen, CC *result);

int num_ntos(CC num, char **deststr, size_t *len);
int num_fput_cc(CC num, FILE *out);
bool num_eq(CC a, CC b, RR eps);

typedef struct num_vec
{
    CC *x;
    size_t n;
} VECTOR_CC;

typedef struct num_table
{
    size_t pt_no;
    CC *x;
    CC *y;
} TABLE_CC;

void 
num_fnext_table_cc(FILE *in, CC **x, CC **y, size_t *n);
void 
num_fnext_table_rr(FILE *in, RR **x, RR **y, size_t *n);

void
num_grid_eqdst_cc(
        CC a,
        CC b,
        size_t n,
        CC *out);
void
num_grid_eqdst_rr(
        RR a,
        RR b,
        size_t n,
        RR *out);
void
num_fill_vals_cc(
        CC_TO_CC f,
        CC *domain,
        size_t n,
        CC *out);
void
num_fill_vals_rr(
        RR_TO_RR f,
        RR *domain,
        size_t n,
        RR *out);
void
num_fill_vals_rproj(
        CC_TO_CC f,
        RR *domain,
        size_t n,
        RR *out);
void
num_table_eqdst_cc(
        CC_TO_CC f,
        CC a,
        CC b,
        size_t n,
        CC **x,
        CC **y);
void
num_fput_table_cc(
        const CC *x,
        const CC *y,
        size_t n,
        FILE *out);
void
num_fputpretty_table_cc(
        const CC *x,
        const CC *y,
        size_t n,
        FILE *out);
bool
num_table_cc_eq(
        const CC *x1,
        const CC *y1,
        size_t n1,
        const CC *x2,
        const CC *y2,
        size_t n2);
bool
num_table_cc_eqt(
        TABLE_CC *t1,
        TABLE_CC *t2);
RR
num_max_deviation_cc(CC *X, CC *Y, size_t pt_no);
CC
NUM_ZERO_CC(CC x);
CC
NUM_IDENTITY_CC(CC x);
CC
NUM_SQR_CC(CC x);

/* affine transformations */

struct affine2
{
    RR a11;
    RR a12;
    RR a21;
    RR a22;
    RR b1;
    RR b2;
    RR alpha;
};

struct vec2hom
{
    RR x, y, alpha;
};
struct vec2rr
{
    RR x, y;
};

#define AFFINE2_ID \
    ((struct affine2){ 1, 0, 0, 1, 0, 0, 1 })

void
affine2apply_hom(struct affine2, RR *x, RR *y, RR *alpha);

struct vec2rr
affine2apply_rr_immut(struct affine2, struct vec2rr);

void
affine2apply_rr(struct affine2, RR *x, RR *y);

struct affine2
affine2mul(struct affine2 A, struct affine2 B);

struct affine2
affine2mul_n(size_t n, struct affine2 sentinel, ...);

struct affine2
affine2rot(RR phi);

struct affine2
affine2tr(RR x, RR y);

struct affine2
affine2scale(RR x, RR y);

/* affine transformations in RR^3 */

struct affine3
{
    RR a11, a12, a13,
       a21, a22, a23,
       a31, a32, a33;
    RR b1, b2, b3;
    RR alpha;
};

struct affine3
affine3mul(struct affine3 T1, struct affine3 T2);

struct affine3
affine3mul_n(size_t n, struct affine3 T1, ...);

void
affine3apply_rr(struct affine3 T, RR *X, RR *Y, RR *Z);

#define AFFINE3_ID \
    ((struct affine3){ .a11=1, .a12=0, .a13=0, .b1=0, \
                       .a21=0, .a22=1, .a23=0, .b2=0, \
                       .a31=0, .a32=0, .a33=1, .b3=0, \
                                               .alpha=1})

struct affine3
affine3rotx(RR phi);

struct affine3
affine3roty(RR phi);

struct affine3
affine3rotz(RR phi);

struct affine3
affine3scale(RR x, RR y, RR z);

struct affine3
affine3tr(RR x, RR y, RR z);

/* Spline interpolation */


/* Common equations solution methods */

RR root_secant_method(RR_TO_RR f, RR x0, RR x1, RR eps);

#endif /* ifdef _VSU_NUM_H_ */
