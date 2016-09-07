#include "vsunum.h"

#ifndef vsu_num_approx_newton_h
#define vsu_num_approx_newton_h

typedef struct numvector NEWTON_DDS;

NUMBER
newton_eval(
        NUMBER x,
        NUMBER const *const grid,
        NUMBER const *const dds,
        size_t pt_no);
NUMBER
newton_evalg(
        NUMBER X,
        const NEWTON_DDS dds,
        const TABLE grid);
NUMBER
newton_approx_equidist_kth_dd(
        int k,
        NUMBER_R h,
        NUMBER const *const vals,
        size_t pt_no);
void
newton_approx_equidist(
        NUMBER_R a,
        NUMBER_R b,
        NUMBER const *const vals,
        size_t pt_no,
        NUMBER **out_dds);
#endif
