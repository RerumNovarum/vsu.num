find_package(PkgConfig)

find_library(CHECK_LIBRARY check)
find_library(M_LIBRARY m)

set(LIBS ${LIBS} ${CHECK_LIBRARY} ${M_LIBRARY} vsu_num vsu_num_approx_newton)

set(ENV{CTEST_OUTPUT_ON_FAILURE} "1")
add_executable(check_num_io check_num_io.c)
add_executable(check_approx_newton check_approx_newton.c)
target_link_libraries(check_num_io ${LIBS})
target_link_libraries(check_approx_newton ${LIBS})
