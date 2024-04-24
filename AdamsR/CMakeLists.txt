add_executable(AdamsR main.h main.cpp steenrodR.h steenrodR.cpp)
target_compile_features(AdamsR PRIVATE cxx_std_17)
target_include_directories(AdamsR PRIVATE ../include ../thirdparty/fmt ../thirdparty/nlohmann_json)
target_link_libraries(AdamsR PRIVATE algebras fmt)

IF (WIN32)
    # VS 二进制运行时的current working directory
    set_target_properties(AdamsR PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Release)
ENDIF()
