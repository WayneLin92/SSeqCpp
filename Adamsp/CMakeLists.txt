add_executable(Adamsp main.h main.cpp steenrodp.h steenrodp.cpp)
target_compile_features(Adamsp PRIVATE cxx_std_17)
target_include_directories(Adamsp PRIVATE ../include ../thirdparty/fmt ../thirdparty/nlohmann_json)
target_link_libraries(Adamsp PRIVATE algebras fmt)

IF (WIN32)
    # VS 二进制运行时的current working directory
    set_target_properties(Adamsp PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Release)
ENDIF()