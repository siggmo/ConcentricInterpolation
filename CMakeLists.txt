cmake_minimum_required(VERSION 3.0...3.28)

# ##############################################################################
# GENERAL SETTINGS
# ##############################################################################

project(
  ConInter
  VERSION 0.1.0
  LANGUAGES C CXX
  DESCRIPTION "Concentric Interpolation")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ##############################################################################
# Options
# ##############################################################################

option(CONINTER_BUILD_SHARED
       "Build ConcentricInterpolation as a shared library." OFF)

option(CONINTER_MKL "Use Intel MKL library for CBLAS and LAPACKE." OFF)

option(CONINTER_DEPS_FALLBACK
       "Skip automatic dependency detection and use hardcoded dependency paths."
       OFF)

option(CONINTER_DEPS_LINK_STATIC
       "Link statically against dependencies LAPACK, LAPACKE and BLAS." ON)

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set(CONINTER_BUILD_EXAMPLES_DEFAULT ON)
endif()
option(CONINTER_BUILD_EXAMPLES "Build examples"
       ${CONINTER_BUILD_EXAMPLES_DEFAULT})

# IPO
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.9")
  cmake_policy(SET CMP0069 NEW)
  include(CheckIPOSupported)
  get_property(languages GLOBAL PROPERTY ENABLED_LANGUAGES)
  check_ipo_supported(RESULT result LANGUAGES ${languages})
  if(result AND NOT DEFINED CMAKE_INTERPROCEDURAL_OPTIMIZATION)
    option(CMAKE_INTERPROCEDURAL_OPTIMIZATION
           "Enable interprocedural optimization for all targets." ON)
    message(STATUS "IPO activated.")
  elseif(NOT result)
    message(STATUS "IPO not supported.")
  endif()
endif()

# ##############################################################################
# RPATH settings
# ##############################################################################

set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(RPATH_BASE "$ORIGIN")
file(RELATIVE_PATH REL_PATH_LIB "/${CMAKE_INSTALL_BINDIR}/"
     "/${CMAKE_INSTALL_LIBDIR}/")
set(CMAKE_INSTALL_RPATH
    "${RPATH_BASE};${RPATH_BASE}/${REL_PATH_LIB}"
    CACHE STRING "Install RPATH")

# ##############################################################################
# REQUIRED PACKAGES
# ##############################################################################

if(CONINTER_MKL)
  if(NOT "$ENV{SETVARS_COMPLETED}" STREQUAL "1")
    message(STATUS "Intel OneAPI setvars script has not been run. Might not be \
able to find MKL.")
  endif()
  set(BLA_VENDOR "Intel10_64lp")
else()
  if("$ENV{SETVARS_COMPLETED}" STREQUAL "1")
    message(
      STATUS
        "Intel OneAPI environment is active, but CONINTER_MKL is not set. This \
might impact which LAPACK and BLAS libraries are found.")
  endif()
  set(BLA_VENDOR "Generic")
endif()

if(NOT CONINTER_DEPS_FALLBACK)
  # The CMake find scripts automatically detect the MKL variants of CBLAS and
  # LAPACKE if setvars has been run. TODO: Sometimes if both
  # intel-oneapi-mkl-devel (from Intel's repo, in /opt) and intel-mkl (from
  # Debian's repo, in /usr/lib/) are installed, find_package uses the ones in
  # /usr/lib/, even though the setvars script has been run.
  set(BLA_STATIC ${CONINTER_DEPS_LINK_STATIC})
  find_package(LAPACK REQUIRED)
  set(LAPACKE_STATIC ${CONINTER_DEPS_LINK_STATIC})
  find_package(LAPACKE REQUIRED)

  # some older FindLAPACK.cmake scripts don't create a target
  if(NOT TARGET LAPACK::LAPACK)
    add_library(LAPACK::LAPACK INTERFACE IMPORTED)
    set_target_properties(
      LAPACK::LAPACK
      PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${LAPACK_INCLUDE_DIRS}"
                 INTERFACE_LINK_DIRECTORIES "${LAPACK_LIBRARY_DIRS}"
                 INTERFACE_LINK_LIBRARIES "${LAPACK_LIBRARIES}")
  endif()

  # create target for LAPACKE since our find script doesn't do that.
  add_library(LAPACKE::LAPACKE INTERFACE IMPORTED)
  set_target_properties(
    LAPACKE::LAPACKE
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${LAPACKE_INCLUDE_DIRS}"
               INTERFACE_LINK_DIRECTORIES "${LAPACKE_LIBRARY_DIRS}"
               INTERFACE_LINK_LIBRARIES "${LAPACKE_LIBRARIES}")

  target_link_libraries(LAPACKE::LAPACKE INTERFACE LAPACK::LAPACK)
else()
  # if the cmake version doesn't ship with (recent) find scripts for LAPACK and
  # BLAS, use legacy hardcoded paths:
  if(NOT CONINTER_MKL)
    set(CONINTER_LAPACKE_INC "/usr/include")
    set(CONINTER_LAPACKE_LIB "/usr/lib/x86_64-linux-gnu/")
    if(CONINTER_DEPS_LINK_STATIC)
      set(CONINTER_LAPACKE_LIBS gfortran m liblapacke.a liblapack.a libblas.a)
    else()
      set(CONINTER_LAPACKE_LIBS liblapacke.so liblapack.so libblas.so)
    endif()
  else()
    set(INTEL "$ENV{ONEAPI_ROOT}")
    set(CONINTER_LAPACKE_INC "${INTEL}/mkl/latest/include")
    set(CONINTER_LAPACKE_LIB
        "${INTEL}/mkl/latest/lib;${INTEL}/compiler/latest/lib")

    if(CONINTER_DEPS_LINK_STATIC)
      set(CONINTER_LAPACKE_LIBS m dl libmkl_intel_lp64.a libmkl_intel_thread.a
                                libmkl_core.a libiomp5.a)
    else()
      set(CONINTER_LAPACKE_LIBS
          m dl libmkl_intel_lp64.so libmkl_intel_thread.so libmkl_core.so
          libiomp5.so)
    endif()
  endif()

  add_library(LAPACKE::LAPACKE INTERFACE IMPORTED)
  set_target_properties(
    LAPACKE::LAPACKE
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${CONINTER_LAPACKE_INC}"
               INTERFACE_LINK_DIRECTORIES "${CONINTER_LAPACKE_LIB}"
               INTERFACE_LINK_LIBRARIES "${CONINTER_LAPACKE_LIBS}")

  # on some older systems one still has to explicitly link against pthreads
  find_package(Threads REQUIRED)
  target_link_libraries(LAPACKE::LAPACKE INTERFACE Threads::Threads)
endif()

# ##############################################################################
# SOURCES
# ##############################################################################

if(NOT CONINTER_BUILD_SHARED)
  add_library(ConInter STATIC)
else()
  add_library(ConInter SHARED)
endif()
add_library(ConInter::ConInter ALIAS ConInter)

if(CONINTER_MKL)
  target_compile_definitions(ConInter PUBLIC CONINTER_MKL)
endif()

# since API and internal headers are equivalent for this project, just set the
# public header paths to src/
target_include_directories(
  ConInter PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>)
add_subdirectory(src)

if(CONINTER_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

# ##############################################################################
# LINKING
# ##############################################################################

target_link_libraries(ConInter PUBLIC LAPACKE::LAPACKE)

# ##############################################################################
# INSTALLING
# ##############################################################################

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
install(
  TARGETS ConInter
  EXPORT ConInter
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
          COMPONENT ConInter_Runtime
          NAMELINK_COMPONENT ConInter_Development
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT ConInter_Development
  INCLUDES
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ConInter
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ConInter
                COMPONENT ConInter_Development
  PRIVATE_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ConInter
                 COMPONENT ConInter_Development)

install(
  EXPORT ConInter
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ConInter
  NAMESPACE ConInter::
  FILE ConInterTargets.cmake
  COMPONENT ConInter_Development)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/ConInterConfigVersion.cmake
  COMPATIBILITY SameMajorVersion)
configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/ConInterConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/ConInterConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/ConInter"
  NO_SET_AND_CHECK_MACRO NO_CHECK_REQUIRED_COMPONENTS_MACRO)
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/ConInterConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/ConInterConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ConInter
  COMPONENT ConInter_Development)

install(
  FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/FindLAPACKE.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ConInter/modules
  COMPONENT ConInter_Development)
