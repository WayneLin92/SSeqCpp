if(MSVC)
    add_compile_options(/permissive-)
endif()

set(HEADER_LIST 
    "${SSeqCpp_SOURCE_DIR}/include/algebras/algebras.h"
    "${SSeqCpp_SOURCE_DIR}/include/algebras/algebrasZ.h"
    "${SSeqCpp_SOURCE_DIR}/include/algebras/benchmark.h"
    "${SSeqCpp_SOURCE_DIR}/include/algebras/database.h"
    "${SSeqCpp_SOURCE_DIR}/include/algebras/dbAdamsSS.h"
    "${SSeqCpp_SOURCE_DIR}/include/algebras/groebner.h"
    "${SSeqCpp_SOURCE_DIR}/include/algebras/linalg.h"
    "${SSeqCpp_SOURCE_DIR}/include/algebras/myexception.h"
    "${SSeqCpp_SOURCE_DIR}/include/algebras/myhash.h"
    "${SSeqCpp_SOURCE_DIR}/include/algebras/myio.h"
    "${SSeqCpp_SOURCE_DIR}/include/algebras/steenrod.h"
    "${SSeqCpp_SOURCE_DIR}/include/algebras/groebner_steenrod.h"
    "${SSeqCpp_SOURCE_DIR}/include/algebras/utility.h"
)

set(CMAKE_CXX_STANDARD 20)
add_library(algebras STATIC algebras.cpp algebrasZ.cpp groebner.cpp database.cpp dbAdamsSS.cpp linalg.cpp myio.cpp steenrod.cpp groebner_steenrod.cpp utility.cpp benchmark.cpp ${HEADER_LIST})
if(MSVC)
    target_compile_options(algebras PRIVATE /W4 /WX)
    # message(STATUS ${CMAKE_CXX_FLAGS_RELEASE})
endif()

target_include_directories(algebras PRIVATE ../include/algebras ../thirdparty/fmt ../thirdparty/nlohmann_json ../thirdparty/sqlite3)
target_link_libraries(algebras PRIVATE sqlite3 fmt)