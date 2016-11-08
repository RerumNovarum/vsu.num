CMAKE_MINIMUM_REQUIRED(VERSION 3.2)
project(vsunum)

find_package(PkgConfig)

find_library(M_LIBRARY m)
set(vsunum_SRC include/vsu/num.h src/num.c src/num-numbers.c)
add_library(vsunum STATIC ${vsunum_SRC})
target_include_directories(vsunum PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(vsunum ${M_LIBRARY})

install(TARGETS vsunum EXPORT vsunum-config DESTINATION lib)
install(EXPORT vsunum-config DESTINATION share/libvsu/)
install(DIRECTORY include DESTINATION include)

enable_testing()
add_executable(test_numbers tests/test_numbers.c)
target_link_libraries(test_numbers vsunum)
add_test(test_numbers test_numbers)

pkg_check_modules(CHECK REQUIRED check)
