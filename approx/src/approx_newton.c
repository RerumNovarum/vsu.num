#include "vsu/approx_newton.h"

NUMBER
newton_eval(
        NUMBER x,
        NUMBER const *const pts,
        NUMBER const *const dds,
        size_t pt_no)
{
    NUMBER y = dds[pt_no-1];
    for (int i = pt_no-1; i > 0; --i)
    {
        y =
            y * (x - pts[i - 1])
            + dds[i - 1];
    }
    return y;
}

NUMBER
newton_evalg(
        NUMBER x,
        const NEWTON_DDS dds,
        const TABLE table)
{
    VECTOR pts = { table.x, table.pt_no };
    return newton_eval(x, pts.x, dds.x, dds.n);
}

void
newton_fill(
        NUMBER const *const pts,
        NUMBER const *const dds,
        size_t dds_pt_no,
        NUMBER const *const dom,
        NUMBER **out,
        size_t pt_no)
{
    (*out) = malloc(pt_no * sizeof(NUMBER));
    for (int i = 0; i < pt_no; ++i)
    {
        (*out)[i] = newton_eval(dom[i], pts, dds, dds_pt_no);
    }
}
