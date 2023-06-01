add_executable(Adams AdamsRes.cpp groebner_res.h groebner_res.cpp groebner_res_const.h groebner_res_const.cpp AdamsResProd.cpp AdamsResExport.cpp AdamsCell.cpp cells.cpp main.h main.cpp AdamsUtilities.cpp)
target_compile_features(Adams PRIVATE cxx_std_17)
target_include_directories(Adams PRIVATE ../include ../thirdparty/fmt)
target_link_libraries(Adams PRIVATE algebras fmt)

IF (WIN32)
    # VS 二进制运行时的current working directory
    set_target_properties(Adams PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Release)
ENDIF()
