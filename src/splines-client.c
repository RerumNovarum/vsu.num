#include <vsu/splines.h>
#include <string.h>
#include <argp.h>

struct argp_option options[] = {
    { "table", 't', "FILE", 0,
        "Input file with a table of knot-value entries (of the form \"x_i f_i\") "
        "one at a line." },
//    { "coeffs", 'c', "FILE", OPTION_ARG_OPTIONAl,
//        "Compute and print coefficients of the spline" },
    { "x", 'x', "NUMBER", 0, "Evaluate the value of the spline at x=NUMBER" },
    { "print", 'p', 0, 0, "Print the grid and coefficients" },
    { 0 }
};

struct arguments
{
    vsunum_cubspline_ptr spline;
};

static int
parse_opt(int key, char *arg,
          struct argp_state *state)
{
    struct arguments *args = state->input;
    switch(key)
    {
        case 't':
            if (args->spline != NULL)
                vsunum_cubspline_free(args->spline);
            args->spline = vsunum_cubspline_load_table(arg);
            if (args->spline == NULL) {
                argp_usage(state);
                return EINVAL;
            }
            return 0;
//      'c':
//          args->spline = vsunum_spline_load_coeffs(arg);
//          break;
        case 'x':
            if (args->spline == NULL) {
                argp_usage(state);
                return EINVAL;
            }
            RR v = vsunum_cubspline_eval(args->spline, num_sgetr(arg, strlen(arg)));
            printf("%Lf\n", v);
            return 0;
        case 'p':
            vsunum_cubspline_print(stdout, args->spline);
            return 0;
        default:
            return ARGP_ERR_UNKNOWN;
    }
}

int main(int argc, char **argv)
{
    struct argp argp = { options, parse_opt, 0, 0 };
    struct arguments arguments = { 0 };
    argp_parse(&argp, argc, argv, 0, 0, &arguments);
}
