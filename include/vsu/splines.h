#ifndef _SPLINES_H_
#define _SPLINES_H_

#include <vsu/num.h>
#include <vsu/linsolve.h>

enum vsunum_cubspline_disposable_field
{
    FREE_X, FREE_A, FREE_B, FREE_C, FREE_D
};

/* represents a cubic spline
 * built on the grid (x_0, ..., x_n)
 * in the form
 *   s(x) = s_i(x) = 
 *   = a_i + b_i (x-x_i) + 1/2 c_i (x-x_i)^2 + 1/6 d_i (x-x_i)^3,
 * where x\in [x_i, x_{i+1}].
 */
struct vsunum_cubspline
{
    RR *x,
       *a, *b, *c, *d; /* output coefficients */
    size_t n; /* = no_points - 1 */
    enum vsunum_cubspline_disposable_field must_free;
};

typedef struct vsunum_cubspline * vsunum_cubspline_ptr;

RR
vsunum_cubspline_eval(vsunum_cubspline_ptr spline, RR x);

struct vsunum_cubspline_builder
{
    RR *x, *f,
       *a, *b, *c, *d;
    size_t n;
    struct tridiag_eqn ce; /* equations for c */
};

typedef struct vsunum_cubspline_builder * vsunum_cubspline_builder_ptr;

void
vsunum_cubspline_eqns_init(vsunum_cubspline_builder_ptr s,
        vsunum_cubspline_ptr spline,
        RR *x, RR *f);

void
vsunum_cubspline_eqns_leftderiv(vsunum_cubspline_builder_ptr s, RR df0);

void
vsunum_cubspline_eqns_rightderiv(vsunum_cubspline_builder_ptr s, RR dfn);

int
vsunum_cubspline_eqns_solve(vsunum_cubspline_builder_ptr s);

int
vsunum_cubspline_from_table_and_derivs(vsunum_cubspline_ptr spline,
        RR *x, RR *f,
        RR df0, RR dfn, size_t n);

vsunum_cubspline_ptr
vsunum_cubspline_alloc(size_t no_points);

void
vsunum_cubspline_free(vsunum_cubspline_ptr spline);

/* Load a spline from a ascii file of the form:
 * N
 * x_0, ..., x_N
 * a_0 b_0 c_0 d_0
 * ...
 * a_{N-1} b_{N-1} c_{N-1} d_{N-1}
 */
vsunum_cubspline_ptr
vsunum_cubspline_load_coeffs(char *filename);

void
vsunum_cubspline_print(FILE *out, vsunum_cubspline_ptr spline);

/* Load a spline from knot-value table
 * and boundary constraints.
 * The input file is ascii-encoded,
 * begins with a number N of knots,
 * then N lines of knot-value records
 * of the form "x_i f_i",
 * and finally two lines with boundary conditions
 * encoded as "key=val"
 * where val is a real number
 * and key like "df0" or "dfn"
 * for derivatives at left and right boundaries
 * respectively. Other constraints
 * are to be implemented.
 */
vsunum_cubspline_ptr
vsunum_cubspline_load_table(char *filename);

#endif /* ifndef  _SPLINES_H_ */
