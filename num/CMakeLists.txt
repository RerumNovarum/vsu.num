project(vsu_num)

find_package(PkgConfig)

find_library(M_LIBRARY m)
add_library(vsu_num STATIC src/vsunum.c)
target_include_directories(vsu_num PUBLIC ./include)
target_link_libraries(vsu_num ${M_LIBRARY})

enable_testing()

pkg_check_modules(CHECK REQUIRED check)

add_executable(check_num_io tests/check_num_io.c)
target_link_libraries(check_num_io vsu_num ${CHECK_LIBRARIES}  ${M_LIBRARY})

add_test(NAME num_io COMMAND check_num_io)
