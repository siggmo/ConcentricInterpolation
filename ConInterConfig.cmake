include(CMakeFindDependencyMacro)
find_dependency(GSL)
find_dependency(LAPACKE)

include(${CMAKE_CURRENT_LIST_DIR}/ConInter.cmake)
