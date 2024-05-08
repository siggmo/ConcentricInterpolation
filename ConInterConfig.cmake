set(CMAKE_MODULE_PATH_save "${CMAKE_MODULE_PATH}")
list(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_CURRENT_LIST_DIR}/customFindScripts")

include(CMakeFindDependencyMacro)
find_dependency(GSL)
find_dependency(LAPACK)
find_dependency(LAPACKE)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH_save}")
unset(CMAKE_MODULE_PATH_save)

include(${CMAKE_CURRENT_LIST_DIR}/ConInter.cmake)
