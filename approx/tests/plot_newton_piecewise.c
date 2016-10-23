#include <cairo.h>

#include "vsu/num.h"
#include "vsu/plot.h"

CC
_num_01x_else_1_sub_xSqr(CC x)
{
    if (ABS(x) < 1) return x;
    return 1 - x*x;
}

int
main()
{
    cairo_surface_t *srf = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, 1024, 1024);
    cairo_t *cr = cairo_create(srf);
    cairo_set_source_rgb(cr, 0,0,0);
    cairo_paint(cr);
    num_cairo_begin_plot(cr, 0, 0, 1024, 1024, -3.2, 3.2, -3.2, 3.2);
    cairo_set_source_rgb(cr, 0, 255, 0);
    cairo_set_line_width(cr, .04);
    num_cairo_plot_func_real(cr, -3, 3, _num_01x_else_1_sub_xSqr, 1024);
    cairo_stroke(cr);
    cairo_set_source_rgb(cr, 255, 0, 0);
    cairo_set_line_width(cr, .02);
    num_cairo_plot_newton_approx_equidist_real(
            cr,
            -3,
            3,
            _num_01x_else_1_sub_xSqr,
            4,
            1024);
    cairo_stroke(cr);
    num_cairo_end_plot(cr);
    cairo_stroke(cr);
    cairo_surface_write_to_png(srf, "piecewise.png");
    cairo_destroy(cr);
    cairo_surface_destroy(srf);
    return 0;
}
