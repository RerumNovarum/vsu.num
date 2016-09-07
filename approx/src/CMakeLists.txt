find_package(PkgConfig REQUIRED)

find_library(CAIRO cairo)
find_library(M_LIB m)
find_path(CAIRO_INCLUDE_DIR NAMES cairo)
find_path(CAIRO_INCLUDES NAMES cairo.h PATHS /usr/include/cairo)
include_directories(${CAIRO_INCLUDE_DIR} ${CAIRO_INCLUDES})

add_library(vsu_num STATIC vsunum.c)
add_library(vsu_num_approx_newton STATIC approx_newton.c approx_newton_equidist.c)
add_library(vsu_num_plot STATIC numplot.c)
set(LIBS ${LIBS} ${CAIRO} ${M_LIB} vsu_num vsu_num_approx_newton vsu_num_plot)

add_executable(approx approx.c)
target_link_libraries(approx ${LIBS})
target_link_libraries(vsu_num_plot vsu_num_approx_newton)
