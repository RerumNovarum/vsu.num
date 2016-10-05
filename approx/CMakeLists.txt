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

add_executable(plot_newton_piecewise tests/plot_newton_piecewise.c)
target_link_libraries(plot_newton_piecewise vsu_approx vsu_plot)
add_test(plot_newton_piecewise plot_newton_piecewise)

add_test(approx_newton check_approx_newton)
add_test(NAME approx_cli_sqr_plot COMMAND approx -f sqr -o sqr.png --white --plot)
add_test(NAME approx_cli_sqr_plot_newton COMMAND approx -f sqr -o sqr-newton.png --white --plot_newtoneq)
add_test(NAME approx_cli_sqr_plot_orig_newton COMMAND approx -f sqr -o sqr-orig-newton.png --black -L0.05 --plot -L0.02 --plot_newtoneq)
add_test(NAME approx_cli_print_sqr_maxdev COMMAND approx -f sqr --newtoneq_err=100)
add_test(NAME approx_cli_random COMMAND approx -f random --ab 0,1 --yY 0,1 --newtoneq_err=100 -o random.png --pt_no=4 --white -L0.001 --plot -L0.005 --plt_pt_no=1536 --plot_newtoneq)
add_test(NAME approx_cli_poly_4 COMMAND approx -f poly_4 --ab 0,1 --yY 0,1 --newtoneq_err=100 -o poly_4.png --pt_no=4 --black -L0.001 --plot -L0.005 --plt_pt_no=1536 --plot_newtoneq)
add_test(NAME approx_cli_log_plot_orig_newton COMMAND approx --ab 0.0000000000000000000001,2 --yY -2,2 --pt_no 16 -f log --newtoneq_err=100 -o log-orig-newton.png --white -L0.005 --plot -L0.01 --plot_newtoneq)
