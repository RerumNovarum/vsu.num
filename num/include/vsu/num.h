#ifndef _VSU_NUM_H_
#define _VSU_NUM_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <ctype.h>
#include <complex.h>


#ifdef DEBUG
#define DPRINTF(...) \
    do { fprintf(stderr, __VA_ARGS__); } while(0)
#else
#define DPRINTF(...)
#endif

typedef long double NUMBER_R;
typedef long double complex NUMBER;

typedef NUMBER (*NUM_TO_NUM)(NUMBER);
typedef NUMBER (*NUMR_TO_NUMR)(NUMBER_R);

#define REAL(x) creall((x))
#define IMAG(x) cimagl((x))

typedef struct numpoint_r
{
    NUMBER_R x, y;
} POINT_R;

typedef struct numpoint
{
    NUMBER x, y;
} POINT;

NUMBER
num_fgetn(FILE *f);
NUMBER
num_sgetn(char *str, size_t sz);
int num_fnextn(FILE *f, NUMBER *result);
int num_snextn(char *str, size_t sz, NUMBER *result);
int num_ntos(NUMBER num, char **deststr, size_t *len);
int num_fputn(NUMBER num, FILE *out);
bool num_eq(NUMBER a, NUMBER b, NUMBER_R eps);

typedef struct numvector
{
    NUMBER *x;
    size_t n;
} VECTOR;

typedef struct numtable
{
    size_t pt_no;
    NUMBER *x;
    NUMBER *y;
} TABLE;

void 
num_table_read(FILE *in, NUMBER **x, NUMBER **y, size_t *n);

void 
num_table_readt(FILE *in, TABLE *t);

void
num_grid_equidist(
        NUMBER a,
        NUMBER b,
        size_t n,
        NUMBER **grid);

void
num_grid_equidist_r(
        NUMBER_R a,
        NUMBER_R b,
        size_t n,
        NUMBER_R **grid);
void
num_fill_vals(
        NUM_TO_NUM f,
        NUMBER *domain,
        size_t n,
        NUMBER **out);
void
num_fill_vals_r(
        NUMR_TO_NUMR f,
        NUMBER_R *domain,
        size_t n,
        NUMBER_R **out);

void
num_fill_vals_frproj(
        NUM_TO_NUM f,
        NUMBER_R *domain,
        size_t n,
        NUMBER_R **out);
void
num_table_equidist(
        NUM_TO_NUM f,
        NUMBER a,
        NUMBER b,
        size_t n,
        NUMBER **x,
        NUMBER **y);
void
num_table_equidistt(
        NUM_TO_NUM f,
        NUMBER a,
        NUMBER b,
        size_t n,
        TABLE *table);

void
num_table_print(
        const NUMBER *x,
        const NUMBER *y,
        size_t n,
        bool pretty,
        FILE *out);

void
num_table_printt(
        TABLE table,
        bool pretty,
        FILE *out);

bool
num_table_eq(
        const NUMBER *x1,
        const NUMBER *y1,
        size_t n1,
        const NUMBER *x2,
        const NUMBER *y2,
        size_t n2);
bool
num_table_eqt(
        TABLE *t1,
        TABLE *t2);

NUMBER
NUM_ZERO(NUMBER x);

NUMBER
NUM_IDENTITY(NUMBER x);

NUMBER
NUM_SQR(NUMBER x);
#endif
