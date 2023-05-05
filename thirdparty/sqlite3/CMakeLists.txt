add_definitions(-DSQLITE_THREADSAFE=0)
add_definitions(-DSQLITE_DEFAULT_MEMSTATUS=0)
add_definitions(-DSQLITE_OMIT_DEPRECATED)

add_library(sqlite3 STATIC sqlite3.h sqlite3.c)
if(UNIX)
target_link_libraries(sqlite3 PRIVATE pthread dl)
endif()