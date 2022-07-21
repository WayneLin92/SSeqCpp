add_executable(ss main.h table_prod.cpp table_ss.cpp add_diff.cpp deduce.cpp SS.cpp main.cpp)
target_compile_features(ss PRIVATE cxx_std_17)
target_include_directories(ss PRIVATE ../include)
target_link_libraries(ss PRIVATE algebras)

IF (WIN32)
    # VS 二进制运行时的current working directory
    set_target_properties(ss PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Release)
ENDIF()