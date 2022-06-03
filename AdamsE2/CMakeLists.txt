MESSAGE(STATUS "Building target AdamsE2")

add_executable(AdamsE2 AdamsE2.cpp groebner_steenrod.h groebner_steenrod.cpp)
target_compile_features(AdamsE2 PRIVATE cxx_std_17)
target_include_directories(AdamsE2 PRIVATE ../include)
target_link_libraries(AdamsE2 PRIVATE algebras)

add_executable(AdamsE2Prod AdamsE2Prod.cpp groebner_steenrod_const.h groebner_steenrod_const.cpp)
target_compile_features(AdamsE2Prod PRIVATE cxx_std_17)
target_include_directories(AdamsE2Prod PRIVATE ../include)
target_link_libraries(AdamsE2Prod PRIVATE algebras)

add_executable(AdamsE2Export AdamsE2Export.cpp)
target_compile_features(AdamsE2Export PRIVATE cxx_std_17)
target_include_directories(AdamsE2Export PRIVATE ../include)
target_link_libraries(AdamsE2Export PRIVATE algebras)

IF (WIN32)
    # VS 二进制运行时的current working directory
    set_target_properties(AdamsE2 PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Release)
    set_target_properties(AdamsE2Prod PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Release)
    set_target_properties(AdamsE2Export PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Release)
ENDIF()
