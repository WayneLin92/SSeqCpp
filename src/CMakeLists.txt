if (WIN32)
    find_package(unofficial-sqlite3 REQUIRED)
else()
    find_package(sqlite3 REQUIRED)
endif()

set(HEADER_LIST "${AlgTop_SOURCE_DIR}/include/algebras/algebras.h"
    "${AlgTop_SOURCE_DIR}/include/algebras/groebner.h"
    "${AlgTop_SOURCE_DIR}/include/algebras/myio.h"
    "${AlgTop_SOURCE_DIR}/include/algebras/myparser.h"
    "${AlgTop_SOURCE_DIR}/include/algebras/myexception.h"
    "${AlgTop_SOURCE_DIR}/include/algebras/database.h"
    "${AlgTop_SOURCE_DIR}/include/algebras/dbalg.h"
    "${AlgTop_SOURCE_DIR}/include/algebras/linalg.h"
    "${AlgTop_SOURCE_DIR}/include/algebras/benchmark.h"
    "${AlgTop_SOURCE_DIR}/include/algebras/utility.h"
)

add_definitions(-D_CRT_SECURE_NO_WARNINGS)

add_library(algebras STATIC algebras.cpp database.cpp dbalg.cpp groebner.cpp linalg.cpp myio.cpp myparser.cpp ${HEADER_LIST})

target_include_directories(algebras PRIVATE ../include/algebras)

target_link_libraries(algebras PRIVATE unofficial::sqlite3::sqlite3)

target_compile_features(algebras PRIVATE cxx_std_17)

# IDEs should put the headers in a nice place
source_group(
    TREE "${PROJECT_SOURCE_DIR}/include"
    PREFIX "Header Files"
    FILES ${HEADER_LIST}
)