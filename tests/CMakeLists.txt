project(vsu_tests)

find_package(PkgConfig)
pkg_check_modules(CHECK REQUIRED check)

set(LIBS ${CHECK_LIBRARIES} vsu_plot vsu_num vsu_approx)

enable_testing()
