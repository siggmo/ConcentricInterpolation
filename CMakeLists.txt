cmake_minimum_required(VERSION 3.22)

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

include_directories(src)

# ##############################################################################
# SOURCES
# ##############################################################################

add_library(ConInter STATIC)
add_library(ConInter::ConInter ALIAS ConInter)
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
