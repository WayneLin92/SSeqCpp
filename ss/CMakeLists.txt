add_executable(main main.cpp)
target_compile_features(main PRIVATE cxx_std_17)

target_include_directories(main PRIVATE ../include)
target_link_libraries(main PRIVATE algebras)
