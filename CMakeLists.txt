cmake_minimum_required(VERSION 3.0...3.28)

# ##############################################################################
# GENERAL SETTINGS
# ##############################################################################

project(
  ConInter
  VERSION 0.1.0
  LANGUAGES CXX
  DESCRIPTION "Concentric Interpolation")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
# set(CMAKE_COLOR_DIAGNOSTICS ON) # not supported below 3.24

add_compile_options(-fdiagnostics-color=always -fmax-errors=1)

# handle rpath weirdness
if(NOT APPLE)
  set(CMAKE_INSTALL_RPATH $ORIGIN)
endif()

option(CMAKE_INTERPROCEDURAL_OPTIMIZATION
       "Enable interprocedural optimization for all targets." ON)

# ##############################################################################
# REQUIRED PACKAGES
# ##############################################################################

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

option(CONINTER_DEPS_FALLBACK "Skip automatic dependency detection and use hardcoded linker paths." OFF)

if (NOT CONINTER_DEPS_FALLBACK)
  set(BLA_STATIC ON)
  set(LAPACKE_STATIC TRUE)
  find_package(LAPACK REQUIRED)
  find_package(LAPACKE REQUIRED)
  find_package(GSL REQUIRED)
endif()

# ##############################################################################
# SOURCES
# ##############################################################################

option(CONINTER_BUILD_SHARED "Build shared library" OFF)
if(NOT CONINTER_BUILD_SHARED)
  add_library(ConInter STATIC)
else()
  add_library(ConInter SHARED)
endif()
add_library(ConInter::ConInter ALIAS ConInter)
# since API and internal headers are equivalent for this project, just set the
# public header paths to src/
target_include_directories(
  ConInter PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>)
add_subdirectory(src)

option(CONINTER_BUILD_EXAMPLES "Build examples" ${PROJECT_IS_TOP_LEVEL})
if(CONINTER_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

option(CONINTER_MKL "Use MKL headers for lapacke and cblas" OFF)
if(CONINTER_MKL)
  target_compile_definitions(ConInter PUBLIC CONINTER_MKL)
endif()

# ##############################################################################
# LINKING
# ##############################################################################

if (NOT CONINTER_DEPS_FALLBACK)
  target_link_libraries(ConInter PRIVATE GSL::gsl ${LAPACKE_LIBRARIES})
else()
  target_link_directories(ConInter PRIVATE "/usr/lib/x86_64-linux-gnu/")
  target_link_libraries(ConInter PRIVATE gsl lapacke)
endif()

# ##############################################################################
# INSTALLING
# ##############################################################################
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
install(
  TARGETS ConInter
  EXPORT ConInter
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT ConInter_Development
  INCLUDES
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ConInter
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ConInter
  PRIVATE_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ConInter)

# install custom cmake find scripts
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindLAPACKE.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ConInter/customFindScripts)

install(
  EXPORT ConInter
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ConInter
  NAMESPACE ConInter::
  FILE ConInter.cmake
  COMPONENT ConInter_Development)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/ConInterConfigVersion.cmake
  COMPATIBILITY SameMajorVersion)
install(FILES ConInterConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/ConInterConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ConInter)
