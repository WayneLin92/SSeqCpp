if(MSVC)
    add_compile_options(/permissive-)
endif()

set(CMAKE_CXX_STANDARD 20)
add_executable(ss main.h pigroebner.h mylog.h database_ss.cpp add_diff.cpp pigroebner.cpp deduce.cpp ss.cpp category.cpp main.cpp homotopy.cpp loadsave.cpp migrate.cpp mylog.cpp staircase.cpp plot.cpp cofseq.cpp utility.cpp naming.cpp)
target_compile_features(ss PRIVATE cxx_std_20)
target_include_directories(ss PRIVATE ../include ../thirdparty/nlohmann_json ../thirdparty/fmt)
target_link_libraries(ss PRIVATE algebras fmt)

IF (WIN32)
    # VS 二进制运行时的current working directory
    set_target_properties(ss PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Release)
ENDIF()

if(MSVC)
    target_compile_options(ss PRIVATE /W4 /WX)
    # set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /Zi")
    # message(STATUS ${CMAKE_CXX_FLAGS_RELEASE})
endif()