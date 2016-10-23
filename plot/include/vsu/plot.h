#ifndef NUMPLOT_H
#define NUMPLOT_H

#include <cairo.h>
#include <err.h>

#include "vsu/num.h"
#include "vsu/approx_newton.h"
#include "vsu/errors.h"

#define NUMPLOT_NO_MEMORY NO_MEMORY

void
num_cairo_begin_plot(
        cairo_t *cr,
        RR L, CC_R R,
        RR B, CC_R T,
        RR a, CC_R b,
        RR y_min, CC_R y_max);
void
num_cairo_end_plot(cairo_t *cr);

void
num_cairo_plot(
        cairo_t *cr,
        RR const *const x,
        RR const *const y,
        size_t pt_no);

void
num_cairo_plot_func_real(
        cairo_t *cr,
        RR a, CC_R b,
        CC_TO_CC f,
        size_t pt_no);
/*
void
num_cairo_plot_fun(
        cairo_t *cr,
        double L, double R, double B, double T,
        RR a, CC_R b, 
        const RR_TO_RR,
        size_t pt_no);
        */

void
num_cairo_plot_newton_real(
        cairo_t *cr,
        const CC const *const x,
        size_t plt_pt_no,
        const CC const *const x_0,
        const CC const *const dds,
        size_t pt_no);

void
num_cairo_plot_newton_approx_equidist_vals(
        cairo_t *cr,
        RR a, CC_R b, 
        const RR const *const y,
        size_t pt_no,
        size_t plt_pt_no);

void
num_cairo_plot_newton_approx_equidist_real(
        cairo_t *cr,
        RR a, CC_R b, 
        CC_TO_CC f,
        size_t pt_no,
        size_t plt_pt_no);

#endif
