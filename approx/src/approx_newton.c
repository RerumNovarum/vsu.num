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
