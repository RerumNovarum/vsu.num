project(vsu_plot)

find_package(PkgConfig)

pkg_check_modules(CAIRO REQUIRED cairo)

include_directories(${CAIRO_INCLUDE_DIRS})
add_library(vsu_plot STATIC src/numplot.c)
target_include_directories(vsu_plot PUBLIC ./include)

target_link_libraries(vsu_plot vsu_num vsu_approx ${CAIRO_LIBRARIES})
