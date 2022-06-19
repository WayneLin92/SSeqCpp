add_executable(ss_basis_prod basis_prod.cpp)
target_compile_features(ss_basis_prod PRIVATE cxx_std_17)
target_include_directories(ss_basis_prod PRIVATE ../include)
target_link_libraries(ss_basis_prod PRIVATE algebras)

IF (WIN32)
    # VS 二进制运行时的current working directory
    set_target_properties(ss_basis_prod PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Release)
ENDIF()