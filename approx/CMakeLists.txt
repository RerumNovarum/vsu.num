project(vsu_num_approx)

add_subdirectory (src)
add_subdirectory (tests)

enable_testing()
add_test(NAME check_num_io COMMAND check_num_io)
add_test(NAME check_newton COMMAND check_approx_newton)
