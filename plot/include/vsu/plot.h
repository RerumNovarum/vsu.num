#ifndef NUMPLOT_H
#define NUMPLOT_H

#include <cairo.h>
#include <err.h>

#include "vsu/num.h"
#include "vsu/approx_newton.h"

#define NUMPLOT_NO_MEMORY 2

void
num_cairo_begin_plot(
        cairo_t *cr,
        NUMBER_R L, NUMBER_R R,
        NUMBER_R B, NUMBER_R T,
        NUMBER_R a, NUMBER_R b,
        NUMBER_R y_min, NUMBER_R y_max);
void
num_cairo_end_plot(cairo_t *cr);

void
num_cairo_polychain(
        cairo_t *cr,
        NUMBER_R const *const x,
        NUMBER_R const *const y,
        size_t pt_no);

void
num_cairo_plot_func_real(
        cairo_t *cr,
        NUMBER_R a, NUMBER_R b,
        NUM_TO_NUM f,
        size_t pt_no);
/*
void
num_cairo_plot_fun(
        cairo_t *cr,
        double L, double R, double B, double T,
        NUMBER_R a, NUMBER_R b, 
        const NUMR_TO_NUMR,
        size_t pt_no);
        */

void
num_cairo_plot_newton_real(
        cairo_t *cr,
        const NUMBER const *const x,
        size_t plt_pt_no,
        const NUMBER const *const x_0,
        const NUMBER const *const dds,
        size_t pt_no);

void
num_cairo_plot_newton_approx_equidist_vals(
        cairo_t *cr,
        NUMBER_R a, NUMBER_R b, 
        const NUMBER_R const *const y,
        size_t pt_no,
        size_t plt_pt_no);

void
num_cairo_plot_newton_approx_equidist_real(
        cairo_t *cr,
        NUMBER_R a, NUMBER_R b, 
        NUM_TO_NUM f,
        size_t pt_no,
        size_t plt_pt_no);

#endif
