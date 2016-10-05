project(vsu_plot)

find_package(PkgConfig)

pkg_check_modules(CAIRO REQUIRED cairo)

add_library(vsu_plot STATIC src/numplot.c)
target_include_directories(vsu_plot PUBLIC include ${CAIRO_INCLUDE_DIRS})

target_link_libraries(vsu_plot vsu_num vsu_approx ${CAIRO_LIBRARIES})

enable_testing()
set(LIBS vsu_plot vsu_num vsu_approx)

add_executable(plot_01x01_grid tests/plot_01x01_grid.c)
target_link_libraries(plot_01x01_grid ${LIBS})
add_test(plot_01x01_grid plot_01x01_grid)

add_executable(plot_sqr tests/plot_sqr.c)
target_link_libraries(plot_sqr ${LIBS})
add_test(plot_sqr plot_sqr)

add_executable(plot_log tests/plot_log.c)
target_link_libraries(plot_log ${LIBS})
add_test(plot_log plot_log)
