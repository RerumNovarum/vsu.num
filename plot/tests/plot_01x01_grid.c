#include <cairo.h>

#include "vsu/num.h"
#include "vsu/plot.h"

int
main()
{
    cairo_surface_t *srf = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, 320, 240);
    cairo_t *cr = cairo_create(srf);
    num_cairo_begin_plot(cr, 0, 0, 320, 240, 0, 1, 0, 1);
    num_cairo_end_plot(cr);
    cairo_stroke(cr);
    cairo_surface_write_to_png(srf, "01x01_grid.png");
    cairo_destroy(cr);
    cairo_surface_destroy(srf);
    return 0;
}
