cmake_minimum_required(VERSION 3.11)

project(
    SSeqCpp
    VERSION 1.0.1
    DESCRIPTION "Spectral sequences in C++"
    LANGUAGES C CXX
)

set(BUILD_DOXYGEN OFF)


if(MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /Zi")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_LINK_FLAGS_RELEASE} /DEBUG")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
endif()

# The compiled library code is here
add_subdirectory(src)

# Only do these if this is the main project, and not if it is included through add_subdirectory
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    # Optionally set things like CMAKE_CXX_STANDARD, CMAKE_POSITION_INDEPENDENT_CODE here
    set (CMAKE_CXX_STANDARD 20)

    # Let's ensure -std=c++xx instead of -std=g++xx
    set(CMAKE_CXX_EXTENSIONS OFF)

    # Support folders in IDEs
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)

    # Set the output directory
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

    # Testing
    include(CTest)

    add_subdirectory(thirdparty/fmt)
    add_subdirectory(thirdparty/sqlite3)

    # The executable code is here
    add_subdirectory(ss)
    # add_subdirectory(MayE2)
    # add_subdirectory(MayE4)
    add_subdirectory(Adams)
    # add_subdirectory(AdamsR)
    # add_subdirectory(Adamsp)
    # add_subdirectory(tests)

    # Docs only available if this is the main app
    find_package(Doxygen)
    if(Doxygen_FOUND AND BUILD_DOXYGEN)
        add_subdirectory(docs)
    else()
        message(STATUS "Not building docs")
    endif()
endif()