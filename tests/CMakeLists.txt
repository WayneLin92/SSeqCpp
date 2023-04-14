find_package(Catch2 CONFIG REQUIRED)
if(MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /Zi")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_LINK_FLAGS_RELEASE} /DEBUG")
endif()

# test_groebner
add_executable(test_groebner test_groebner.cpp)
target_compile_features(test_groebner PRIVATE cxx_std_17)
target_include_directories(test_groebner PRIVATE ../include)
target_link_libraries(test_groebner PRIVATE Catch2::Catch2 algebras)

add_test(NAME test_groebner COMMAND test_groebner)

# test_database
add_executable(test_database test_database.cpp)
target_compile_features(test_database PRIVATE cxx_std_17)
target_include_directories(test_database PRIVATE ../include)
target_link_libraries(test_database PRIVATE Catch2::Catch2 algebras)

add_test(NAME test_database COMMAND test_database)

# A temporary target
add_executable(tmp tmp.cpp)
target_compile_features(tmp PRIVATE cxx_std_17)
target_include_directories(tmp PRIVATE ../include ../thirdparty/fmt)
target_link_libraries(tmp PRIVATE algebras fmt)

# target for benchmark
add_executable(benchmark benchmark.cpp)
target_compile_features(benchmark PRIVATE cxx_std_17)
target_include_directories(benchmark PRIVATE ../include)
target_link_libraries(benchmark PRIVATE algebras)