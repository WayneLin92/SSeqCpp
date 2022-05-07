if (WIN32)
    find_package(unofficial-sqlite3 CONFIG REQUIRED)
else()
    find_package(sqlite3 REQUIRED)
endif()
if(MSVC)
    add_compile_options(/permissive-)
endif()

set(HEADER_LIST "${AlgTop_SOURCE_DIR}/include/algebras/algebras.h"
    "${AlgTop_SOURCE_DIR}/include/algebras/groebner.h"
    "${AlgTop_SOURCE_DIR}/include/algebras/myio.h"
    "${AlgTop_SOURCE_DIR}/include/algebras/myexception.h"
    "${AlgTop_SOURCE_DIR}/include/algebras/database.h"
    "${AlgTop_SOURCE_DIR}/include/algebras/dbalg.h"
    "${AlgTop_SOURCE_DIR}/include/algebras/linalg.h"
    "${AlgTop_SOURCE_DIR}/include/algebras/benchmark.h"
    "${AlgTop_SOURCE_DIR}/include/algebras/utility.h"
    "${AlgTop_SOURCE_DIR}/include/algebras/steenrod.h"
)

add_library(algebras STATIC algebras.cpp database.cpp dbalg.cpp groebner.cpp linalg.cpp myio.cpp steenrod.cpp utility.cpp ${HEADER_LIST})
if(MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /Zi")
    # message(STATUS ${CMAKE_CXX_FLAGS_RELEASE})
endif()

target_include_directories(algebras PRIVATE ../include/algebras)

if (WIN32)
    target_link_libraries(algebras PRIVATE unofficial::sqlite3::sqlite3)
else()
    target_link_libraries(algebras PRIVATE sqlite3)  # ?
endif()

target_compile_features(algebras PRIVATE cxx_std_17)

# IDEs should put the headers in a nice place
source_group(
    TREE "${PROJECT_SOURCE_DIR}/include"
    PREFIX "Header Files"
    FILES ${HEADER_LIST}
)