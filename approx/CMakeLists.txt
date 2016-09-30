cmake_minimum_required(VERSION 3.5)
project(vsu_num_approx)

find_package(PkgConfig REQUIRED)
pkg_check_modules(CAIRO REQUIRED cairo)
find_library(M_LIBRARY m)

include_directories(PRIVATE ${CAIRO_INCLUDE_DIRS})

add_library(vsu_approx STATIC src/approx_newton.c src/approx_newton_equidist.c)
target_include_directories(vsu_approx PUBLIC ./include)
target_link_libraries(vsu_approx vsu_num ${CAIRO_LIBRARIES})

add_executable(approx src/approx.c)
target_link_libraries(approx ${M_LIBRARY} vsu_num vsu_plot ${CAIRO_LIBRARIES})

enable_testing()

pkg_check_modules(CHECK REQUIRED check)

add_executable(check_approx_newton tests/check_approx_newton.c)
target_link_libraries(check_approx_newton vsu_approx vsu_num ${M_LIBRARY} ${CHECK_LIBRARIES})

add_test(NAME approx_newton COMMAND check_approx_newton)
