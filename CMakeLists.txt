cmake_minimum_required(VERSION 3.16)

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

# add flags from presets
add_compile_options(${PRESET_CXX_FLAGS})

# handle rpath weirdness
if(NOT APPLE)
  set(CMAKE_INSTALL_RPATH $ORIGIN)
endif()

# ##############################################################################
# REQUIRED PACKAGES
# ##############################################################################

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(BLA_STATIC ON)
set(LAPACKE_STATIC TRUE)
find_package(LAPACK REQUIRED)
find_package(LAPACKE REQUIRED)
find_package(GSL REQUIRED)

# ##############################################################################
# HEADER SEARCH PATHS
# ##############################################################################

# include_directories(src)

# ##############################################################################
# SOURCES
# ##############################################################################

add_library(ConInter STATIC)
add_library(ConInter::ConInter ALIAS ConInter)
# since API and internal headers are equivalent for this project, just set the
# public header paths to src/
target_include_directories(
  ConInter PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>)
add_subdirectory(src)

add_subdirectory(examples)

# ##############################################################################
# LINKING
# ##############################################################################

target_link_libraries(ConInter PRIVATE GSL::gsl lapacke)

# ##############################################################################
# NAMING
# ##############################################################################

# ##############################################################################
# INSTALLING
# ##############################################################################
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
install(
  TARGETS ConInter
  EXPORT ConInter
#   RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
#   LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT ConInter_Development
  INCLUDES
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ConInter
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ConInter
  PRIVATE_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ConInter)

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
