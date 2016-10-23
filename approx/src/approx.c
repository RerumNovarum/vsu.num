#include <stdio.h>
#include <argp.h>
#include <err.h>
#include <cairo.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

#include "vsu/errors.h"
#include "vsu/num.h"
#include "vsu/plot.h"
#include "vsu/approx_newton.h"

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
#define O_WRITE_PLOT 2
#define O_NEWTON_EQUIDIST_PRINT_ERR 3
#define O_WHITE 4
#define O_BLACK 5
#define O_AB 6
#define O_YMIN_YMAX 7

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
    { "y_min", O_Y_MIN, "CC", 0,
        "Defines lower bound of range to be drawn; "
            "Defaults to 0"},
    { "y_max", O_Y_MAX, "CC", 0,
        "Defines upper bound of range to be drawn; "
        "Defaults to 1" },
    { "x_min", O_X_MIN, "CC", 0, "Left bound of domain [a,b]; Defaults to 0" },
    { "x_max", O_X_MAX, "CC", 0, "Right bound of domain [a,b]; Defaults to 1" },
    { "func", O_FUNC, "FUNCTION_NAME", 0, "Function to be examined" },
    { "pt_no", O_PT_NO, "UNSIGNED", 0, "Number of points" },
    { "plt_pt_no", O_PLT_PT_NO, "UNSIGNED", 0,
        "Number of pivot points for a plot; "
        "[x,X] will be split into `plt_pt_no` parts "
        "at which the graph will be drawn"
        " as a line segments between pivot points `x[i],f(x[i])`" },
    { "linewidth", O_LINEWIDTH, "REAL", 0, "Linewidth" },
    { "plot", O_PLOT, 0, 0, "Plot a function graph" },
    { "plot_newtoneq", O_PLOT_NEWTON_EQUIDIST, 0, 0,
        "Plot a graph of interpolation polynomial built in Newton form on equidistant grid" },
    { "write", O_WRITE_PLOT, 0, 0, "Flush current image to current OUTPUT_FILE" },
    { "newtoneq_err", O_NEWTON_EQUIDIST_PRINT_ERR, "POINTS", 0,
        "Print maximum absolute deviation computed at points "
            "of equidistant grid covering [x,X] with POINTS points" },
    { "white", O_WHITE, 0, 0, "Clear background white" },
    { "black", O_BLACK, 0, 0, "Paint background black" },
    { "ab", O_AB, "A,B", 0, "Interval [a,b] in the form of two comma-delimited numbers" },
    { "yY", O_YMIN_YMAX, "y,Y", 0, "Comma-delimited y_min,y_max" },
    { 0 }
};

CC  _num_log_re(NUMBER x)
{
    return logl(REAL(x));
}

CC  _num_random(NUMBER x)
{
	return rand()*1.0/RAND_MAX;
}

CC _num_poly_4(NUMBER x)
{
	CC x2 = x*x;
	return 1 + x2*(-10 + x2);
}

static struct
{
    char * name;
    size_t name_len;
    CC_TO_CC func;
} functions[] =
{
    /* lexicographically ordered */
    /* shorter string is always */
    /* lexicographically smaller than longer */
    { "id", 2, NUM_IDENTITY },
    { "sqr", 3, NUM_SQR },
    { "log", 3, _num_log_re },
    { "poly_4", 3, _num_poly_4 },
	{ "random", 6, _num_random }
};

static CC_TO_CC function_by_name(char *name, size_t len)
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
    CC_TO_CC func;
    size_t width, height;
    RR a, b;  /* function domain [a,b] */
    CC *grid;
    CC *vals;
    RR *pltgrid;
    RR y_min, y_max;
    CC *dds;
    bool dds_pending;
    size_t pt_no;
    size_t plt_pt_no;
    double linewidth;
    cairo_t *cr;
    cairo_surface_t *surface;
    bool vals_pending;  /* need to recompute the vals */
    bool pltgrid_pending; /* &c */
    bool grid_pending;
    bool cairo_pending; /* need to create new cairo context */
    bool write_pending;
};

void
static compute_pendings(struct arguments *args)
{
    if (args->grid_pending)
    {
        if (args->b <= args->a)
        {
			fprintf(stderr, "a=%Ld, b=%Ld", args->a, args->b);
            err(INVALID_INPUT, "$a$ should be strictly less than $b$");
        }
        if (args->grid != NULL)
        {
            free(args->grid);
            args->grid = NULL;
        }
        num_grid_gen_eqdst_cc(args->a, args->b, args->pt_no, &(args->grid));
        if (args->vals_pending)
        {
            if (args->vals != NULL)
            {
                free(args->vals); /* not realloc since it'd require error handling */
                args->vals = NULL;
            }
            num_fill_vals(args->func, args->grid, args->pt_no, &(args->vals));
        }
    }
    args->grid_pending = false;
    args->vals_pending = false;
}

void
static compute_pendings_plt(struct arguments *args)
{
    if (args->cairo_pending)
    {
        if (args->cr != 0) cairo_destroy(args->cr);
        if (args->surface != 0) cairo_surface_destroy(args->surface);
        cairo_surface_t *surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, args->width, args->height);
        cairo_t *cr = cairo_create(surface);
        args->surface = surface;
        args->cr = cr;
        /* num_cairo_end_plot is in O_WRITE */
        num_cairo_begin_plot(
                args->cr,
                0, 0,
                args->width, args->height,
                args->a, args->b,
                args->y_min, args->y_max);
    }
    if (args->pltgrid_pending)
    {
        if (args->b <= args->a)
        {
            err(INVALID_INPUT, "$a$ should be strictly less than $b$");
        }
        if (args->pltgrid != NULL)
        {
            free(args->pltgrid);
			args->pltgrid = NULL;
        }
        num_grid_gen_eqdst_rr(args->a, args->b, args->plt_pt_no, &(args->pltgrid));
    }
    args->pltgrid_pending = false;
    args->cairo_pending = false;
}

void
static compute_pendings_dds(struct arguments *args)
{
    if (!args->dds_pending) return;
    compute_pendings(args);
    args->dds_pending = false;
    if (args->dds != NULL) free(args->dds);
    newton_approx_equidist(
            args->a, args->b,
            args->vals,
            args->pt_no,
            &(args->dds));
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
			args->cairo_pending = true;
            break;
            /*
             * TODO: replace $REAL\circ num_sgetn$ with something real-specific
             *       or, at least, check if imaginary part is zero
             */
        case O_WIDTH:
            args->width = (size_t) REAL(num_sgetn(arg, strlen(arg)));
			args->pltgrid_pending = true;
			args->cairo_pending = true;
            break;
        case O_HEIGHT:
			args->height = (size_t) REAL(num_sgetn(arg, strlen(arg)));
			args->cairo_pending = true;
			break;
		case O_Y_MIN:
			args->y_min = (double) REAL(num_sgetn(arg, strlen(arg)));
			args->cairo_pending = true;
			break;
		case O_Y_MAX:
			args->y_max = (double) REAL(num_sgetn(arg, strlen(arg)));
            break;
        case O_X_MIN:
            args->a = REAL(num_sgetn(arg, strlen(arg)));
            args->grid_pending = true;
			args->pltgrid_pending = true;
			args->cairo_pending = true;
            break;
        case O_X_MAX:
            args->b = REAL(num_sgetn(arg, strlen(arg)));
            args->grid_pending = true;
			args->pltgrid_pending = true;
			args->cairo_pending = true;
            break;
        case O_AB:
			{
				char *a = arg;
				size_t len = strlen(arg);
				CC z;
				int read = num_snextn(a, len, &z);
				a += read;
				len -= read;
				args->a = REAL(z);
				num_snextn(a, len, &z);
				args->b = REAL(z);
				args->grid_pending = true;
				args->pltgrid_pending = true;
				args->cairo_pending = true;
				break;
			}
        case O_YMIN_YMAX:
			{
				char *a = arg;
				size_t len = strlen(arg);
				CC z;
				int read = num_snextn(a, len, &z);
				a += read;
				len -= read;
				args->y_min = REAL(z);
				num_snextn(a, len, &z);
				args->y_max = REAL(z);
				args->cairo_pending = true;
				break;
			}
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
			args->grid_pending = true;
			args->pltgrid_pending = true;
			args->cairo_pending = true;
			break;
        case O_LINEWIDTH:
            args->linewidth = REAL(num_sgetn(arg, strlen(arg)));
            break;
        case O_BLACK:
			compute_pendings(args);
            compute_pendings_plt(args);
            cairo_set_source_rgb(args->cr, 0, 0, 0);
			cairo_paint(args->cr);
            break;
        case O_WHITE:
			compute_pendings(args);
            compute_pendings_plt(args);
            cairo_set_source_rgb(args->cr, 0xff, 0xff, 0xff);
			cairo_paint(args->cr);
            break;
        case O_PLT_PT_NO:
            args->plt_pt_no = (size_t) REAL(num_sgetn(arg, strlen(arg)));
            break;
        case O_PLOT:
			{
				compute_pendings(args);
				compute_pendings_plt(args);
				args->write_pending = true;
				cairo_set_source_rgb(args->cr, 0, 0xff, 0);
				cairo_set_line_width(args->cr, args->linewidth);
				RR *pltvals;
				num_fill_vals_frproj(args->func, args->pltgrid, args->plt_pt_no, &(pltvals));
				num_cairo_plot(
						args->cr,
						args->pltgrid,
						pltvals,
						args->plt_pt_no);
				free(pltvals);
				cairo_stroke(args->cr);
				break;
			}
		case O_PLOT_NEWTON_EQUIDIST:
			{
				compute_pendings(args);
				compute_pendings_dds(args);
				compute_pendings_plt(args);
				args->write_pending = true;
				cairo_set_source_rgb(args->cr, 0xff, 0, 0);
				cairo_set_line_width(args->cr, args->linewidth);
				RR *pltvals = malloc(sizeof(CC_R) * args->plt_pt_no);;
				for (int i = 0; i < args->plt_pt_no; ++i)
				{
					pltvals[i] = REAL(newton_eval(args->pltgrid[i], args->grid, args->dds, args->pt_no));
				}
				num_cairo_plot(args->cr, args->pltgrid, pltvals, args->plt_pt_no);
				free(pltvals);
				cairo_stroke(args->cr);
				break;
			}
		case O_WRITE_PLOT:
            compute_pendings_plt(args);
			cairo_set_line_width(args->cr, args->linewidth);
			cairo_set_source_rgb(args->cr, 0, 0, 0);
			/* num_cairo_begin_plot is in compute_pendings() */
			num_cairo_end_plot(args->cr);
			cairo_stroke(args->cr);
			cairo_set_source_rgba(args->cr, 0, 0, 0, 0);
			cairo_surface_write_to_png(args->surface, args->output_fname);
            break;
        case O_NEWTON_EQUIDIST_PRINT_ERR:
			compute_pendings(args);
            compute_pendings_dds(args);
            size_t err_pt_no = (size_t) REAL(num_sgetn(arg, strlen(arg)));
            CC *grid;
            CC *expected;
            CC *result;
            num_grid_gen_eqdst_cc(args->a, args->b, err_pt_no, &grid);
            num_fill_vals(args->func, grid, err_pt_no, &expected);
            newton_fill(
                    args->grid, args->dds, args->pt_no,
                    grid, &result, err_pt_no);
            RR dev = num_max_deviation(expected, result, err_pt_no);
            free(grid);
            free(expected);
            free(result);
            FILE *f;
            bool use_std= (strcmp(args->output_fname, "-") == 0);
            if (use_std) f = stdout;
            else f = fopen(args->output_fname, "w");
            fprintf(f, "maxdev=%Lf\n", dev);
            if (!use_std) fclose(f);
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
    args.func = NULL;
    args.a = 0;
    args.b = 1;
    args.width = 2048;
    args.height = 2048;
    args.y_min = 0;
    args.y_max = 1;
    args.pt_no = 4;
    args.plt_pt_no = 1536;
    args.linewidth = .01;
    args.grid = NULL;
    args.vals = NULL;
    args.pltgrid = NULL;
    args.cr = NULL;
    args.surface = NULL;
    args.cairo_pending = true;
    args.pltgrid_pending = true;
    args.vals_pending = true;
    args.grid_pending = true;
	args.dds = NULL;
	args.dds_pending = true;

    // error_t argp_parse (argp, argc, argv, flags, arg_index, input)
    argp_parse(&argp, argc, argv, 0, 0, &args);
    if (args.write_pending)
    {
        struct argp_state state;
        state.input = &args;
        parse_opt(O_WRITE_PLOT, "write", &state);
    }
}
