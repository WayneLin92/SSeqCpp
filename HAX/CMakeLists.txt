add_executable(MayE2 MayE2.cpp utilities2.cpp MayE2.h)
target_compile_features(MayE2 PRIVATE cxx_std_17)
target_include_directories(MayE2 PRIVATE ../include)
target_link_libraries(MayE2 PRIVATE algebras sqlite3)

add_executable(B9 B9.cpp B9.h)
target_compile_features(B9 PRIVATE cxx_std_17)
target_include_directories(B9 PRIVATE ../include)
target_link_libraries(B9 PRIVATE algebras sqlite3)
