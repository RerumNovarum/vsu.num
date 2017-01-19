#include <vsu/splines.h>
#include <vsu/linsolve.h>
#include <vsu/num.h>
#include <string.h>
#include <stdio.h>
#include <err.h>

vsunum_cubspline_ptr
vsunum_cubspline_alloc(size_t n)
{
    size_t pt_n = n+1;
    size_t bytes_needed = 5*sizeof(RR)*pt_n + sizeof(struct vsunum_cubspline);
    /* TODO: consider memset()-ing zeroes */
    void *buf = malloc(bytes_needed);
    if (buf == NULL) return NULL;
    vsunum_cubspline_ptr spline = buf;
    spline->x = buf + sizeof(struct vsunum_cubspline);
    spline->a = spline->x + pt_n;
    spline->b = spline->a + (pt_n-1);
    spline->c = spline->b + (pt_n-1);
    spline->d = spline->c + (pt_n-1);
    spline->n = n;
    spline->must_free = FREE_X;
    return spline;
}

void
vsunum_cubspline_free(vsunum_cubspline_ptr spline)
{
    if (spline->must_free & FREE_X)
        free(spline->x);
    if (spline->must_free & FREE_A)
        free(spline->a);
    if (spline->must_free & FREE_B)
        free(spline->b);
    if (spline->must_free & FREE_C)
        free(spline->c);
    if (spline->must_free & FREE_D)
        free(spline->d);
}

void
vsunum_cubspline_eqns_init(vsunum_cubspline_builder_ptr s,
        vsunum_cubspline_ptr spline,
        RR *x, RR *f)
{
    size_t n = spline->n;
    /* TODO: consider memset()-ing zeroes into coeffs vectors */
    *s = (struct vsunum_cubspline_builder) {
        .a = spline->a,
        .b = spline->b,
        .c = spline->c,
        .d = spline->d,
        .x = x,
        .f = f,
        .n = spline->n,
        .ce = (struct tridiag_eqn) {
            .a = spline->a,
            .b = spline->b,
            .c = spline->d,
            .f = spline->c,
            .n = spline->n
        }
    };
    RR *cea = s->ce.a;
    RR *ceb = s->ce.b;
    RR *cec = s->ce.c;
    RR *cef = s->ce.f;
    RR h, hnext;
    h = x[1] - x[0];
    for (uint32_t i = 0; i <= n-2; ++i) {
        hnext = x[i+2] - x[i+1];
        cea[i+1] = h/6;
        ceb[i+1] = (h + hnext)/3;
        cec[i+1] = hnext/6;
        cef[i+1] = (f[i+2] - f[i+1])/hnext - (f[i+1] - f[i])/h;
        h = hnext;
    }
}

void
vsunum_cubspline_eqns_leftderiv(vsunum_cubspline_builder_ptr s, RR df0)
{
    RR h0 = s->x[1] - s->x[0];
    RR f1 = s->f[1];
    RR f0 = s->f[0];
    s->ce.a[0] = 0;
    s->ce.b[0] = 1.0/3;
    s->ce.c[0] = 1.0/6;
    s->ce.f[0] = ((f1 - f0)/h0 - df0)/h0;
}

void
vsunum_cubspline_eqns_rightderiv(vsunum_cubspline_builder_ptr s, RR dfn)
{

    RR hnm2, hnm1;
    size_t n;
    n = s->n;
    hnm2 = s->x[n-1] - s->x[n-2];
    hnm1 = s->x[n] - s->x[n-1];
    s->ce.a[n-1] = hnm2/3;
    s->ce.b[n-1] = hnm2*2/3 + hnm1/2;
    s->ce.c[n-1] = 0; 
    s->ce.f[n-1] = 3*(s->f[n] - s->f[n-1])/hnm1
                    - 2*(s->f[n-1] - s->f[n-2])/hnm2
                    - dfn;
}

int
vsunum_cubspline_eqns_solve(vsunum_cubspline_builder_ptr s)
{
    RR *a, *b, *c, *d, *x, *f;
    a = s->a; b = s->b; c = s->c; d = s->d;
    x = s->x; f = s->f;
    size_t n = s->n;

    /* s->ce.f coincides with s->c
     * so the solution of tridiagonal system s->ce
     * is stored in the vector s->c.
     */
    int ret = tridiag_eqn_solve(&s->ce);
    if (ret != 0)
        return ret;

    for (uint32_t i = 0; i <= n-2; ++i) {
        RR h_i = x[i+1] - x[i];
        a[i] = f[i];
        b[i] = (f[i+1] - f[i])/h_i - h_i*(c[i] + c[i+1]/2)/3;
        d[i] = (c[i+1] - c[i])/h_i;
    }

    RR hnm2, hnm1;
    hnm2 = x[n-1] - x[n-2];
    hnm1 = x[n] - x[n-1];
    a[n-1] = f[n-1];
    b[n-1] = b[n-2] + hnm2*c[n-2] +hnm2/2*d[n-2];
    d[n-1] = (f[n] - f[n-1] - hnm1*b[n-1] - hnm1*c[n-1]/2)/hnm1*2/hnm1*3/hnm1;
}

int
vsunum_cubspline_from_table_and_derivs(vsunum_cubspline_ptr spline,
        RR *x, RR *f,
        RR df0, RR dfn, size_t n)
{
    struct vsunum_cubspline_builder s = { 0 };
    vsunum_cubspline_eqns_init(&s, spline, x, f);
    vsunum_cubspline_eqns_leftderiv(&s, df0);
    vsunum_cubspline_eqns_rightderiv(&s, dfn);
    int ret = vsunum_cubspline_eqns_solve(&s);
    return ret;
}

static inline uint32_t
_spline_find_interval(vsunum_cubspline_ptr spline, RR x0)
{
    uint32_t l, r, m;
    RR *x = spline->x;
    if (x0 < x[0] || x0 > x[spline->n])
        err(1, "_spline_find_interval: out of bounds");
    l = 0;
    r = spline->n-1;
    while ((r-l) > 1)
    {
        m = (l+r)>>1;
        if (x[m] <= x0) l = m;
        if (x[m] >= x0) r = m;
    }
    if (x0 == x[r])
        return r;
    return l;
}

RR
vsunum_cubspline_eval(vsunum_cubspline_ptr spline, RR x)
{
    uint32_t i;
    RR v;

    i = _spline_find_interval(spline, x);
    x = x - spline->x[i];
    v = spline->a[i] +
            x*(spline->b[i] +
                    x*(spline->c[i]/2.0 +
                        x*spline->d[i]/6.0));
    return v;
}

vsunum_cubspline_ptr
vsunum_cubspline_load_table(char *filename)
{
    int ret;

    FILE *f = fopen(filename, "r");
    if (f == NULL) return NULL;

    size_t n;
    fscanf(f, "%zu", &n);

    vsunum_cubspline_ptr spline =
        vsunum_cubspline_alloc(n);
    RR *vals = malloc(sizeof(RR)*n);
    if (vals == NULL || spline == NULL) {
        vsunum_cubspline_free(spline);
        return NULL;
    }
    for (uint32_t i = 0; i <= n; ++i) {
        num_fnextr_unlocked(f, spline->x + i);
        num_fnextr_unlocked(f, vals + i);
    }

    struct vsunum_cubspline_builder builder = { 0 };
    vsunum_cubspline_eqns_init(&builder, spline, spline->x, vals);

    struct {
        char *key;
        void (*handler)(vsunum_cubspline_builder_ptr, RR);
    } boundaries[] = {
        { "df0", vsunum_cubspline_eqns_leftderiv },
        { "dfn", vsunum_cubspline_eqns_rightderiv }
    };
    size_t boundaries_no = sizeof(boundaries)/sizeof(*boundaries);
    for (int i = 0; i < 2; ++i) {
        char key[4];
        RR val;
        ret = fscanf(f, "%*[ \n\r\t]%3[^=]=%Lf", key, &val);
        if (ret == EOF)
            err(1, "vsunum_cubspline_load_table: bad input format");
        for (int b = 0; b < boundaries_no; ++b) {
            if (strcmp(boundaries[b].key, key) == 0) {
                (boundaries[b]).handler(&builder, val);
            }
        }
    }
    vsunum_cubspline_eqns_solve(&builder);
    spline->x = builder.x;
    spline->a = builder.a;
    spline->b = builder.b;
    spline->c = builder.c;
    spline->d = builder.d;
    free(vals);
    return spline;
}

void
vsunum_cubspline_print(FILE *out, vsunum_cubspline_ptr spline)
{
    flockfile(out);
    fprintf(out, "%zu\n", spline->n);
    for (int i = 0; i < spline->n; ++i)
        fprintf(out, "%Lf ", spline->x[i]);
    fprintf(out, "%Lf\n", spline->x[spline->n]);
    for (int i = 0; i < spline->n; ++i) {
        fprintf(out, "%Lf %Lf %Lf %Lf\n",
                spline->a[i], spline->b[i],
                spline->c[i], spline->d[i]);
    }
    funlockfile(out);
}
