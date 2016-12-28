CMAKE_MINIMUM_REQUIRED(VERSION 3.2)
project(vsunum)

find_package(PkgConfig)

find_library(M_LIBRARY m)
set(vsunum_SRC include/vsu/num.h src/num.c src/num-numbers.c src/affine2.c src/secant-method.c)
add_library(vsunum STATIC ${vsunum_SRC})
target_include_directories(vsunum PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(vsunum ${M_LIBRARY})

install(TARGETS vsunum EXPORT vsunum-config DESTINATION lib)

add_library(liblinsolve STATIC src/liblinsolve.c)
target_include_directories(liblinsolve PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

add_executable(linsolve src/linsolve.c)
target_link_libraries(linsolve liblinsolve)

install(TARGETS liblinsolve EXPORT vsunum-config DESTINATION lib)
install(TARGETS linsolve DESTINATION bin)

install(EXPORT vsunum-config DESTINATION share/libvsu/)
install(DIRECTORY include DESTINATION ${CMAKE_INSTALL_PREFIX})


enable_testing()

add_executable(test_numbers tests/test_numbers.c)
target_link_libraries(test_numbers vsunum)
add_test(test_numbers test_numbers)

add_executable(test_affine2 tests/test_affine2.c)
target_link_libraries(test_affine2 vsunum)
add_test(test_affine2 test_affine2)

add_executable(test_secant tests/secant-method.c)
target_link_libraries(test_secant vsunum)
add_test(test_secant test_secant)
# pkg_check_modules(CHECK REQUIRED check)
