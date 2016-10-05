cmake_minimum_required(VERSION 3.5)
project(vsu_num_approx)

find_package(PkgConfig REQUIRED)
find_library(M_LIBRARY m)


add_library(vsu_approx STATIC src/approx_newton.c src/approx_newton_equidist.c)
target_include_directories(vsu_approx PUBLIC include)
target_link_libraries(vsu_approx vsu_plot vsu_num ${M_LIBRARY})

add_executable(approx src/approx.c)
target_link_libraries(approx vsu_approx vsu_plot)

enable_testing()

pkg_check_modules(CHECK REQUIRED check)

add_executable(check_approx_newton tests/check_approx_newton.c)
target_link_libraries(check_approx_newton vsu_approx ${CHECK_LIBRARIES})

add_test(approx_newton check_approx_newton)
