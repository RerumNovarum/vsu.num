#include "vsu/plot.h"
#define _VSU_NUM_PLOT_TICK_WIDTH_COEF_ 1e-2

struct _numplot_space
{
    RR l, r, b, t;
};

cairo_user_data_key_t
static _numplot_space_key;

void
num_cairo_begin_plot(
        cairo_t *cr,
        RR X_0, CC_R Y_0,
        RR W, CC_R H,
        RR l, CC_R r,
        RR b, CC_R t)
{
    cairo_save(cr);
    cairo_translate(cr, X_0, Y_0 + H);
    cairo_scale(cr, W/(r-l), -H/(t-b));
    cairo_translate(cr, -l, -b);

    struct _numplot_space *space = malloc(sizeof(struct _numplot_space));
    space->l = l;
    space->r = r;
    space->b = b;
    space->t = t;
    cairo_status_t _st =
        cairo_set_user_data(cr, &_numplot_space_key, space, free);
    if (_st != CAIRO_STATUS_SUCCESS)
        err(NUMPLOT_NO_MEMORY, "cairo can not set userdata: no memory");
}

void
num_cairo_end_plot(cairo_t *cr)
{
    struct _numplot_space *s = cairo_get_user_data(cr, &_numplot_space_key);
    double l = s->l, r = s->r, b = s->b, t = s->t;

    cairo_rectangle(cr, l, b, r-l, t-b);
    RR tick_width;
    tick_width = (t - b)*_VSU_NUM_PLOT_TICK_WIDTH_COEF_;
    for (RR x = l; x < r; x += 1.0)
    {
        cairo_move_to(cr, x, b);
        cairo_line_to(cr, x, b + tick_width);
        cairo_move_to(cr, x, t);
        cairo_line_to(cr, x, t - tick_width);
    }

    tick_width = (r - l)*_VSU_NUM_PLOT_TICK_WIDTH_COEF_;
    for (RR y = b; y < t; y += 1.0)
    {
        cairo_move_to(cr, l, y);
        cairo_line_to(cr, l + tick_width, y);
        cairo_move_to(cr, r, y);
        cairo_line_to(cr, r - tick_width, y);
    }
    cairo_restore(cr);
}

void
num_cairo_plot(
        cairo_t *cr,
        const RR const *const x,
        const RR const *const y,
        size_t n)
{
    if (n <= 0) return;

    cairo_move_to(cr, x[0], y[0]);
    for (int i = 0; i < n; ++i)
    {
        cairo_line_to(cr, x[i], y[i]);
    }
}

void
num_cairo_plot_func_real(
        cairo_t *cr,
        RR a, CC_R b, 
        const CC_TO_CC f,
        size_t pt_no)
{
    if (pt_no <= 0) return;

    CC *X, *Y;

    num_table_equidist(f, a, b, pt_no, &X, &Y);
    cairo_move_to(cr, REAL(X[0]), REAL(Y[0]));
    for (int i = 0; i < pt_no; ++i)
    {
        RR x, y;
        x = REAL(X[i]);
        y = REAL(Y[i]);
        DPRINTF("line_to(%Lf, %Lf)\n", x, y);
        cairo_line_to(cr, REAL(X[i]), REAL(Y[i]));
    }
}

void
num_cairo_plot_newton_real(
        cairo_t *cr,
        const CC const *const x,
        size_t plt_pt_no,
        const CC const *const x_0,
        const CC const *const dds,
        size_t pt_no)
{
    if (pt_no <= 0) return;

    cairo_move_to(cr, REAL(x[0]), REAL(newton_eval(x[0], x_0, dds, pt_no)));
    for (int i = 0; i < plt_pt_no; ++i)
    {
        cairo_line_to(
                cr,
                REAL(x[i]),
                REAL(newton_eval(x[i], x_0, dds, pt_no)));
    }
}

void
num_cairo_plot_newton_approx_equidist_vals(
        cairo_t *cr,
        RR a, CC_R b, 
        RR const *const y,
        size_t pt_no,
        size_t plt_pt_no)
{
    if (pt_no <= 0) return;

    /* yeah, thats waste of space */
    CC *vals = malloc(sizeof(NUMBER) * pt_no);
    for (int i = 0; i < pt_no; ++i)
        vals[i] = (CC) y[i];

    CC *dds;
    CC *grid;
    CC *plt_grid;
    newton_approx_equidist(a, b, vals, pt_no, &dds);
    num_grid_gen_eqdst_cc(a, b, pt_no, &grid);
    num_grid_gen_eqdst_cc(a, b, plt_pt_no, &plt_grid);
    num_cairo_plot_newton_real(cr, plt_grid, plt_pt_no, grid, dds, pt_no);
    free(grid);
    free(plt_grid);
    free(dds);
    free(vals);
    
}

void
num_cairo_plot_newton_approx_equidist_real(
        cairo_t *cr,
        RR a, CC_R b, 
        CC_TO_CC f,
        size_t pt_no,
        size_t plt_pt_no)
{
    CC *plt_grid;
    CC *grid;
    CC *vals;
    CC *dds;

    num_table_equidist(f, a, b, pt_no, &grid, &vals);
    num_grid_gen_eqdst_cc(a, b, plt_pt_no, &plt_grid);
    newton_approx_equidist(a, b, vals, pt_no, &dds);
    num_cairo_plot_newton_real(cr, plt_grid, plt_pt_no, grid, dds, pt_no);

    free(plt_grid);
    free(grid);
    free(vals);
    free(dds);
}
