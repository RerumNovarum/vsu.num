#include <cairo.h>

#include "vsu/num.h"
#include "vsu/plot.h"

int
main()
{
    cairo_surface_t *srf = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, 400, 400);
    cairo_t *cr = cairo_create(srf);
    num_cairo_begin_plot(cr, 0, 0, 400, 400, -1, 1, -1, 1);
    num_cairo_plot_func_real(cr, -1, 1, NUM_SQR, 128);
    num_cairo_end_plot(cr);
    cairo_stroke(cr);
    cairo_surface_write_to_png(srf, "sqr.png");
    cairo_destroy(cr);
    cairo_surface_destroy(srf);
    return 0;
}

