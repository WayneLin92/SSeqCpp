set(DOXYGEN_EXTRACT_ALL YES)
set(DOXYGEN_BUILTIN_STL_SUPPORT YES)
set(DOXYGEN_GENERATE_HTML YES)
set(DOXYGEN_GENERATE_TREEVIEW YES)
set(DOXYGEN_HTML_EXTRA_STYLESHEET "${PROJECT_SOURCE_DIR}/thirdparty/doxygen-awesome-css/doxygen-awesome.css" "${PROJECT_SOURCE_DIR}/thirdparty/doxygen-awesome-css/doxygen-awesome-sidebar-only.css" "${PROJECT_SOURCE_DIR}/thirdparty/doxygen-awesome-css/custom.css")

set(DOXYGEN_EXCLUDE_SYMBOLS "_*")

doxygen_add_docs(
    docs
    ALL
    ${PROJECT_SOURCE_DIR}/include
    COMMENT "Generate html pages"
)
