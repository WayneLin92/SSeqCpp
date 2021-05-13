set(HEADER_LIST "${AlgTop_SOURCE_DIR}/include/algebras/algebras.h"
    "${AlgTop_SOURCE_DIR}/include/algebras/myio.h"
    "${AlgTop_SOURCE_DIR}/include/algebras/myparser.h"
)

add_library(algebras STATIC algebras.cpp database.cpp groebner.cpp linalg.cpp myio.cpp myparser.cpp ${HEADER_LIST})

target_include_directories(algebras PUBLIC ../include/algebras ../thirdparty/sqlite3/include)

target_link_libraries(algebras PUBLIC sqlite3)

target_compile_features(algebras PUBLIC cxx_std_17)

# IDEs should put the headers in a nice place
source_group(
    TREE "${PROJECT_SOURCE_DIR}/include"
    PREFIX "Header Files"
    FILES ${HEADER_LIST}
)