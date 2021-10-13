find_package(Catch2 CONFIG REQUIRED)

# A temporary target
add_executable(tmp tmp.cpp)
target_compile_features(tmp PRIVATE cxx_std_17)
target_include_directories(tmp PRIVATE ../include)
target_link_libraries(tmp PRIVATE algebras)

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