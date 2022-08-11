MESSAGE(STATUS "Building target AdamsE2")

add_executable(Adams AdamsRes.cpp groebner_steenrod.h groebner_steenrod.cpp groebner_steenrod_const.h AdamsResProd.cpp AdamsResExport.cpp main.h main.cpp)
target_compile_features(Adams PRIVATE cxx_std_17)
target_include_directories(Adams PRIVATE ../include)
target_link_libraries(Adams PRIVATE algebras)

IF (WIN32)
    # VS 二进制运行时的current working directory
    set_target_properties(Adams PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Release)
ENDIF()
