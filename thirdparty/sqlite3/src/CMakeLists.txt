add_library(sqlite3 STATIC shell.c sqlite3.c ../include/sqlite3/sqlite3.h ../include/sqlite3/sqlite3ext.h)

target_include_directories(sqlite3 PUBLIC ../include/sqlite3)
target_compile_features(sqlite3 PUBLIC cxx_std_17)