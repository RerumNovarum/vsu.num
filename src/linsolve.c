#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vsu/linsolve.h>

#define PROGRAM_VERSION \
    "linsolve-001"
#define PROGRAM_AUTHOR \
    "rerumnovarum@openmailbox.org"

static int exit_status = EXIT_SUCCESS;

static inline size_t
count_numbers(char *line)
{
    if (line == NULL)
        return 0;
    size_t n = 0;
    char *tail = NULL;
    while (1)
    {
        strtod(line, &tail);
        if (line == tail) return n;
        line = tail;
        ++n;
    }
}

static inline size_t
read_vector(char *line, size_t n, RR *out)
{
    if (line == NULL)
        return 0;
    size_t n_read = 0;
    char *tail = NULL;
    while (n_read < n)
    {
        out[n_read] = strtod(line, &tail);
        if (line == tail) return n;
        line = tail;
        ++n_read;
    }
    return tail == line ? (n_read - 1) : n_read;
}

static void
solve_file(char *file)
{
    FILE *f = fopen(file, "r");
    if (f == NULL) {
        fprintf(stderr, "cannot open file \"%s\"", file);
        exit_status |= 1;
        return;
    }
    char *line = NULL;
    size_t len = 0;
    getline(&line, &len, f);

    size_t n = count_numbers(line);
    lin_eqn_ptr eq = lin_eqn_alloc(n);
    for (int i = 0; i < n; ++i) {
        read_vector(line, n, eq->A[i]);
        free(line);
        line = NULL;
        len = 0;
        getline(&line, &len, f);
    }
    read_vector(line, n, eq->b);
    free(line);
    fclose(f);

    linsolve(eq);
    size_t origlen = strlen(file);
    char *solfile = alloca(origlen + 5);
    strcpy(solfile, file);
    strcpy(solfile + origlen, ".sol");

    f = fopen(solfile, "w");
    for (int i = 0; i < (n-1); ++i) {
            fprintf(f, "%Lf ", eq->x[i]);
    }
    fprintf(f, "%Lf" "\n", eq->x[n-1]);
    fclose(f);
}

static void
print_usage(void)
{
    puts(PROGRAM_VERSION "\n\
            linsolve FILE [FILE ...]\n\
\n\
            Solves a system of linear equations Ax=b\n\
            Layout of FILE is as following:\n\
            a_11 ... a_1n\n\
            .............\n\
            a_n1 ... a_nn\n\
            b_1  ... b_n\n\
email bugs to <" PROGRAM_AUTHOR ">");
}

int main(int argc, char **argv)
{
    if (argc == 1)
        print_usage();
    for (int i = 1; i < argc; ++i) {
        solve_file(argv[i]);
    }
    return exit_status;
}
