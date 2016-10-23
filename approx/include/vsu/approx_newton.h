#include "vsu/num.h"

#ifndef vsu_num_approx_newton_h
#define vsu_num_approx_newton_h

typedef struct numvector NEWTON_DDS;

CC
newton_eval(
        CC x,
        CC const *const grid,
        CC const *const dds,
        size_t pt_no);
CC
newton_evalg(
        CC X,
        const NEWTON_DDS dds,
        const TABLE grid);
CC
newton_approx_equidist_kth_dd(
        int k,
        RR h,
        CC const *const vals,
        size_t pt_no);
void
newton_approx_equidist(
        RR a,
        RR b,
        CC const *const vals,
        size_t pt_no,
        CC **out_dds);

void
newton_fill(
        CC const *const pts,
        CC const *const dds,
        size_t dds_pt_no,
        CC const *const dom,
        CC **out,
        size_t pt_no);
#endif
