add_executable(ss main.h pigroebner.h mylog.h basis.cpp database_ss.cpp add_diff.cpp pigroebner.cpp add_ext.cpp deduce.cpp ss.cpp diagram.cpp main.cpp homotopy.cpp htpy_def.cpp config.cpp migrate.cpp mylog.cpp staircase.cpp plot.cpp cofseq.cpp)
target_compile_features(ss PRIVATE cxx_std_17)
target_include_directories(ss PRIVATE ../include ../thirdparty/nlohmann_json ../thirdparty/fmt)
target_link_libraries(ss PRIVATE algebras fmt)

IF (WIN32)
    # VS 二进制运行时的current working directory
    set_target_properties(ss PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Release)
ENDIF()