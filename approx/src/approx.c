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

const char *argp_program_version = "approx 0.01";
const char *argp_program_bug_address = "<rerumnovarum@openmailbox.org>";

static char doc[] =
"TODO: Description";

static char args_doc[] =
"";

static struct argp_option options[] = {
    { "input", 'i', "INPUT_FILE", 0, "Input file; Defaults to '-' (stdin)" },
    { "output", 'o', "OUTPUT_FILE", 0, "Output file; Defaults to '-' (stdout)" },
    { "print", 'p', 0, 0, "Print table" },
    { "width", 'w', "INTEGER", 0, "Output image width" },
    { "height", 'h', "INTEGER", 0, "Output image height" },
    { "y_min", 'y', "NUMBER", 0,
        "Defines lower bound of range to be drawn; "
            "Defaults to 0"},
    { "y_max", 'Y', "NUMBER", 0,
        "Defines upper bound of range to be drawn; "
        "Defaults to 1" },
    { "x_min", 'x', "NUMBER", 0, "Left bound of domain [a,b]; Defaults to 0" },
    { "x_max", 'X', "NUMBER", 0, "Right bound of domain [a,b]; Defaults to 1" },
    { "func", 'f', "FUNCTION_NAME", 0, "Function to be examined" },
    { "pt_no", 'n', "UNSIGNED", 0, "Number of points" },
    { "plt_pt_no", 'P', "UNSIGNED", 0,
        "Number of pivot points for a plot; "
        "[x,X] will be split into `plt_pt_no` parts "
        "at which the graph will be drawn as a line segments between pivot points `x[i],f(x[i])`" },
    { "linewidth", 'L', "REAL", 0, "Linewidth" },
    { "dev-equidist", 'd', 0, 0,
        "Plot function and its global interpolation polynomial "
        "built on equidistant grid on [a,b]; Requires --func; "
        "Uses interpolation polynomial of degree `pt_no`; "
        "Set number of plot pivot points with `plt_pt_no`" },
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
    { "log_re", 6, _log_re }
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
    NUMBER_R a, b;
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

static error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
    struct arguments *args = state->input;
    switch (key)
    {
        case 'i':
            args->input_fname = arg;
            break;
        case 'o':
            args->output_fname = arg;
            break;
            /*
             * TODO: replace $REAL\circ num_sgetn$ with something real-specific
             *       or, at least, check if imaginary part is zero
             */
        case 'w':
            args->width = (size_t) REAL(num_sgetn(arg, strlen(arg)));
            break;
        case 'h':
            args->height = (size_t) REAL(num_sgetn(arg, strlen(arg)));
            break;
        case 'y':
            args->y_min = (double) REAL(num_sgetn(arg, strlen(arg)));
            break;
        case 'Y':
            args->y_max = (double) REAL(num_sgetn(arg, strlen(arg)));
            break;
        case 'x':
            args->a = REAL(num_sgetn(arg, strlen(arg)));
            break;
        case 'X':
            args->b = REAL(num_sgetn(arg, strlen(arg)));
            break;
        case 'f':
            args->func = function_by_name(arg, strlen(arg));
            if (args->func == 0)
            {
                err(INVALID_INPUT, "unknown function: `%s`", arg);
            }
            break;
        case 'n':
            /* TODO: replace with posix functions? */
            args->pt_no = (size_t) REAL(num_sgetn(arg, strlen(arg)));
            break;
        case 'L':
            args->linewidth = REAL(num_sgetn(arg, strlen(arg)));
            break;
        case 'P':
            args->plt_pt_no = (size_t) REAL(num_sgetn(arg, strlen(arg)));
            break;
        case 'd':
            if (args->func == 0)
            {
                err(INVALID_INPUT, "function is not defined");
            }
            cmd_deviation_equidist(args);
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
    args.width = 1536;
    args.height = 1536;
    args.y_min = 0;
    args.y_max = 1;
    args.pt_no = 4;
    args.plt_pt_no = 128;
    args.linewidth = .1;

    // error_t argp_parse (argp, argc, argv, flags, arg_index, input)
    argp_parse(&argp, argc, argv, 0, 0, &args);
}
