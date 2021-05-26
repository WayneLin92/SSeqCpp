find_package(Catch2 REQUIRED)

add_executable(mytest test.cpp)
target_compile_features(mytest PRIVATE cxx_std_17)
target_include_directories(mytest PRIVATE ../include)
target_link_libraries(mytest PRIVATE Catch2::Catch2 algebras)

add_test(NAME algebraslibtest COMMAND mytest) # Command can be a target
