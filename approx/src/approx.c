#include <stdio.h>
#include <argp.h>
#include <err.h>
#include <cairo.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

#include "vsu/num.h"
#include "vsu/plot.h"
#include "vsu/approx_newton.h"

#define SUCCESS 0;
#define INVALID_INPUT 1
#define NO_MEMORY 2

#define O_INPUT_FILE 'i'
#define O_OUTPUT_FILE 'o'
#define O_PRINT 'P'
#define O_WIDTH 'w'
#define O_HEIGHT 'h'
#define O_Y_MIN 'y'
#define O_Y_MAX 'Y'
#define O_X_MIN 'x'
#define O_X_MAX 'X'
#define O_FUNC 'f'
#define O_PT_NO 'n'
#define O_PLT_PT_NO 1
#define O_LINEWIDTH 'L'
#define O_PLOT 'p'
#define O_PLOT_NEWTON_EQUIDIST 'N'

const char *argp_program_version = "approx 0.01";
const char *argp_program_bug_address = "<rerumnovarum@openmailbox.org>";

static char doc[] =
"TODO: Description";

static char args_doc[] =
"";

static struct argp_option options[] = {
    { "input", O_INPUT_FILE, "INPUT_FILE", 0, "Input file; Defaults to '-' (stdin)" },
    { "output", O_OUTPUT_FILE, "OUTPUT_FILE", 0, "Output file; Defaults to '-' (stdout)" },
    { "print", O_PRINT, 0, 0, "Print table" },
    { "width", O_WIDTH, "INTEGER", 0, "Output image width" },
    { "height", O_HEIGHT, "INTEGER", 0, "Output image height" },
    { "y_min", O_Y_MIN, "NUMBER", 0,
        "Defines lower bound of range to be drawn; "
            "Defaults to 0"},
    { "y_max", O_Y_MAX, "NUMBER", 0,
        "Defines upper bound of range to be drawn; "
        "Defaults to 1" },
    { "x_min", O_X_MIN, "NUMBER", 0, "Left bound of domain [a,b]; Defaults to 0" },
    { "x_max", O_X_MAX, "NUMBER", 0, "Right bound of domain [a,b]; Defaults to 1" },
    { "func", O_FUNC, "FUNCTION_NAME", 0, "Function to be examined" },
    { "pt_no", O_PT_NO, "UNSIGNED", 0, "Number of points" },
    { "plt_pt_no", O_PLT_PT_NO, "UNSIGNED", 0,
        "Number of pivot points for a plot; "
        "[x,X] will be split into `plt_pt_no` parts "
        "at which the graph will be drawn as a line segments between pivot points `x[i],f(x[i])`" },
    { "linewidth", O_LINEWIDTH, "REAL", 0, "Linewidth" },
    { "plot", O_PLOT, 0, 0, "Plot a function graph" },
    { "plot_newton_equidist", PLOT_NEWTON_EQUIDIST, 0,
        "Plot a graph of interpolation polynomial built in Newton form on equidistant grid" }
    { 0 }
};

NUMBER  _log_re(NUMBER x)
{
    return logl(REAL(x));
}

static struct
{
    char * name;
    size_t name_len;
    NUM_TO_NUM func;
} functions[] =
{
    /* lexicographically ordered */
    /* shorter string is always */
    /* lexicographically smaller than longer */
    { "id", 2, NUM_IDENTITY },
    { "sqr", 3, NUM_SQR },
    { "log", 6, _log_re }
};

static NUM_TO_NUM function_by_name(char *name, size_t len)
{
    /* TODO: binary search it */
    size_t f_no = sizeof(functions)/sizeof(*functions);

    for (int i = 0; i < f_no; ++i)
    {
        if (strcmp(name, functions[i].name) == 0)
        {
            return functions[i].func;
        }
    }
    return 0;
}

struct arguments
{
    char *input_fname;
    char *output_fname;
    NUM_TO_NUM func;
    size_t width, height;
    NUMBER_R *pltgrid;
    bool pltgrid_pending;  /* need to recompute the pltgrid */
    NUMBER_R a, b;  /* function domain [a,b] */
    bool grid_equidist_pending;
    NUMBER *grid;
    NUMBER *vals;
    NUMBER_R y_min, y_max;
    size_t pt_no;
    size_t plt_pt_no;
    double linewidth;
};

void
static cmd_deviation_equidist(struct arguments *args)
{
    cairo_surface_t *surface;
    cairo_t *cr;

    surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, args->width, args->height);
    cr = cairo_create(surface);

    /* TODO: background */

    cairo_set_source_rgb(cr, 0, 0, 0);
    num_cairo_begin_plot(
            cr,
            0, 0,
            args->width, args->height,
            args->a, args->b,
            args->y_min, args->y_max);
    cairo_set_line_width(cr, args->linewidth);
    cairo_stroke(cr);

    cairo_set_source_rgb(cr, 0, 0xff, 0);
    num_cairo_plot_func_real(
            cr,
            args->a, args->b,
            args->func,
            args->plt_pt_no);
    cairo_stroke(cr);

    cairo_set_source_rgb(cr, 0xff, 0, 0);
    cairo_set_line_width(cr, .5 * args->linewidth);
    num_cairo_plot_newton_approx_equidist_real(
            cr,
            args->a,
            args->b,
            args->func,
            args->pt_no,
            args->plt_pt_no);
    cairo_stroke(cr);

    cairo_set_line_width(cr, args->linewidth);
    cairo_set_source_rgb(cr, 0, 0, 0);
    num_cairo_end_plot(cr);
    cairo_stroke(cr);
    cairo_set_source_rgba(cr, 0, 0, 0, 0);

    cairo_surface_write_to_png(surface, args->output_fname);
    cairo_destroy(cr);
    cairo_surface_destroy(surface);
}

void
static recompute_grid(arguments *args)
{
    if (args->grid_equidist_pending)
    {
        if (b <= a)
        {
            err(INVALID_INPUT, "$a$ should be strictly less than $b$");
        }
        num_grid_equidist(args->a, args->b, args->grid, args->pt_no);
        args->grid_pending = false;
    }
}

void
static recompute_pltgrid(arguments *args)
{
    if (args->pltgrid_pending)
    {
        if (b <= a)
        {
            err(INVALID_INPUT, "$a$ should be strictly less than $b$");
        }
        num_grid_equidist_r(args->a, args->b, args->pltgrid, args->plt_pt_no);
        args->pltgrid_pending = false;
    }
}



static error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
    struct arguments *args = state->input;
    switch (key)
    {
        case O_INPUT_FILE:
            args->input_fname = arg;
            break;
        case O_OUTPUT_FILE:
            args->output_fname = arg;
            break;
            /*
             * TODO: replace $REAL\circ num_sgetn$ with something real-specific
             *       or, at least, check if imaginary part is zero
             */
        case O_WIDTH:
            args->width = (size_t) REAL(num_sgetn(arg, strlen(arg)));
            args->pltgrid_pending = true;
            break;
        case O_HEIGHT:
            args->height = (size_t) REAL(num_sgetn(arg, strlen(arg)));
            break;
        case O_Y_MIN:
            args->y_min = (double) REAL(num_sgetn(arg, strlen(arg)));
            break;
        case O_Y_MAX:
            args->y_max = (double) REAL(num_sgetn(arg, strlen(arg)));
            break;
        case O_X_MIN:
            args->a = REAL(num_sgetn(arg, strlen(arg)));
            args->grid_pending = true;
            break;
        case O_X_MAX:
            args->b = REAL(num_sgetn(arg, strlen(arg)));
            args->grid_pending = true;
            break;
        case O_FUNC:
            args->func = function_by_name(arg, strlen(arg));
            if (args->func == 0)
            {
                err(INVALID_INPUT, "unknown function: `%s`", arg);
            }
            break;
        case O_PT_NO:
            /* TODO: replace with posix functions? */
            args->pt_no = (size_t) REAL(num_sgetn(arg, strlen(arg)));
            break;
        case O_LINEWIDTH:
            args->linewidth = REAL(num_sgetn(arg, strlen(arg)));
            break;
        case O_PLT_PT_NO:
            args->plt_pt_no = (size_t) REAL(num_sgetn(arg, strlen(arg)));
            break;
        case O_PLOT:
            recompute_grid(args);
            recompute_pltgrid(args);
            break;
        case O_PLOT:
            recompute_grid(args);
            recompute_pltgrid(args);
            break;
        case ARGP_KEY_ARG:
            break;
        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc };

int
main(int argc, char **argv)
{
    struct arguments args;
    args.input_fname = "-";
    args.output_fname = "-";
    args.func = 0;
    args.a = 0;
    args.b = 1;
    args.grid_equidist_pending = true;
    args.table = 0;
    args.width = 1536;
    args.height = 1536;
    args.pltgrid_pending = true;
    args.y_min = 0;
    args.y_max = 1;
    args.pt_no = 3;
    args.plt_pt_no = 128;
    args.linewidth = .1;

    // error_t argp_parse (argp, argc, argv, flags, arg_index, input)
    argp_parse(&argp, argc, argv, 0, 0, &args);
}
