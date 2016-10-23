#include "vsu/approx_newton.h"

CC
newton_eval(
        CC x,
        CC const *const pts,
        CC const *const dds,
        size_t pt_no)
{
    CC y = dds[pt_no-1];
    for (int i = pt_no-1; i > 0; --i)
    {
        y =
            y * (x - pts[i - 1])
            + dds[i - 1];
    }
    return y;
}

CC
newton_evalg(
        CC x,
        const NEWTON_DDS dds,
        const TABLE table)
{
    VECTOR pts = { table.x, table.pt_no };
    return newton_eval(x, pts.x, dds.x, dds.n);
}

void
newton_fill(
        CC const *const pts,
        CC const *const dds,
        size_t dds_pt_no,
        CC const *const dom,
        CC **out,
        size_t pt_no)
{
    (*out) = malloc(pt_no * sizeof(CC));
    for (int i = 0; i < pt_no; ++i)
    {
        (*out)[i] = newton_eval(dom[i], pts, dds, dds_pt_no);
    }
}
