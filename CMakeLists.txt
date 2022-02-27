cmake_minimum_required(VERSION 3.14)
project(levmar VERSION 2.6 LANGUAGES C CXX)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

################################### Options ####################################
include(utils)
option_with_default(CMAKE_BUILD_TYPE "Build type" Release)
option_with_print(BUILD_SHARED_LIBS "Build final library as shared, not static" OFF)
option_with_default(BUILD_FPIC "Libraries will be compiled with position independent code" ON)
if(${BUILD_SHARED_LIBS} AND NOT ${BUILD_FPIC})
    message(FATAL_ERROR "BUILD_SHARED_LIBS ON and BUILD_FPIC OFF are incompatible, as shared library requires position independent code")
endif()

option(ENABLE_LAPACK "Use vendor Lapack implementation. If disabled, then use internal implementation (self-contained)." OFF)
if (ENABLE_LAPACK)
  set(HAVE_LAPACK 1)
else()
  set(HAVE_LAPACK 0)
endif()

option(BUILD_DEMO "Build demo program?" TRUE)

option(ENABLE_SINGLE_PREC "enable single precision computation." ON)
option(ENABLE_DOUBLE_PREC "enable double precision computation." ON)

######################### Process & Validate Options ###########################
include(autocmake_safeguards)
#include(autocmake_static_library)

################################# Main Project #################################

# gnu compatibility,
# see https://cmake.org/cmake/help/latest/module/GNUInstallDirs.html
include(GNUInstallDirs)

# build for double precision
if(ENABLE_DOUBLE_PREC)
  set(LM_DBL_PREC 1)
endif()

# build for single precision
if(ENABLE_SINGLE_PREC)
  set(LM_SNGL_PREC 1)
endif()

# For add_library STATIC/SHARED is governed by BUILD_SHARED_LIBS
add_library(levmar)
target_sources(levmar PRIVATE
  src/lm.c
  src/Axb.c
  src/misc.c
  src/lmlec.c
  src/lmbc.c
  src/lmblec.c
  src/lmbleic.c
  src/misc.h
  src/compiler.h
  )

configure_file(
  ${CMAKE_SOURCE_DIR}/src/levmar.h.in
  ${CMAKE_BINARY_DIR}/src/levmar.h
  )

target_include_directories(levmar PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/src>
    $<INSTALL_INTERFACE:include>
  )

set_target_properties(levmar PROPERTIES POSITION_INDEPENDENT_CODE ${BUILD_FPIC})

set_target_properties(levmar PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(levmar PUBLIC m)

if (${BUILD_SHARED_LIBS})
    target_link_libraries(levmar PRIVATE ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${LIBC_INTERJECT})
endif()

add_library(levmar::levmar ALIAS levmar)

# libtool hook (so that the build can used inside ssbu with autotools)
include(libtool_utils)

if(${BUILD_SHARED_LIBS})
	create_libtool_file(levmar ${CMAKE_BINARY_DIR}/liblevmar.so ${CMAKE_BINARY_DIR} )
else()
	create_libtool_file(levmar ${CMAKE_BINARY_DIR}/liblevmar.a ${CMAKE_BINARY_DIR} )
endif()


set(INSTALL_TARGETS levmar)

if(BUILD_DEMO)

  add_executable(lmdemo)
  target_sources(lmdemo PRIVATE src/lmdemo.c)
  target_link_libraries(lmdemo PRIVATE levmar::levmar)

  add_executable(lm_demo_fit)
  target_sources(lm_demo_fit PRIVATE src/lm_demo_fit.c)
  target_link_libraries(lm_demo_fit PRIVATE levmar::levmar)

  add_executable(lm_demo_fit2)
  target_sources(lm_demo_fit2 PRIVATE src/lm_demo_fit2.c)
  target_link_libraries(lm_demo_fit2 PRIVATE levmar::levmar)

endif(BUILD_DEMO)


################################# EXPORT INTERFACE #################################

#target_compile_definitions(levmar INTERFACE USING_${PROJECT_NAME})
#target_include_directories(levmar INTERFACE
#                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

################################# EXPORT CONFIG #################################
include(CMakePackageConfigHelpers)

# setup some variables
set(version_config ${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake)
set(project_config_src ${PROJECT_SOURCE_DIR}/${PROJECT_NAME}-config.cmake.in)
set(project_config_dst ${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config.cmake)
set(targets_export_name ${PROJECT_NAME}-targets)


# important variables
set(INSTALL_LIBDIR ${CMAKE_INSTALL_LIBDIR} CACHE STRING
  "Installation directory for libraries, relative to ${CMAKE_INSTALL_PREFIX}.")

set(INSTALL_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR} CACHE STRING
  "Installation directory for include files, relative to ${CMAKE_INSTALL_PREFIX}.")

set(INSTALL_PKGCONFIG_DIR ${CMAKE_INSTALL_LIBDIR}/pkgconfig CACHE PATH
  "Installation directory for pkgconfig (.pc) files, relative to ${CMAKE_INSTALL_PREFIX}.")

set(INSTALL_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/cmake CACHE STRING
  "Installation directory for cmake files, relative to ${CMAKE_INSTALL_PREFIX}.")

# Generate the version, config and target files into the build directory.
write_basic_package_version_file(
  ${version_config}
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion)

# Generate cmake my_package-config.cmake file
configure_package_config_file(
  ${project_config_src}
  ${project_config_dst}
  INSTALL_DESTINATION ${INSTALL_CMAKE_DIR})
  # Use a namespace because CMake provides better diagnostics
  # for namespaced imported targets.
export(
  TARGETS ${INSTALL_TARGETS} NAMESPACE levmar::
  FILE ${PROJECT_BINARY_DIR}/${targets_export_name}.cmake)

# macro helper to generate pkg-config file
include(generate_pkg_config)

# generate the pkg-config file
generate_pkg_config(src/levmar)

################################# INSTALL LIBRARY #################################

# install header
install(
  FILES ${CMAKE_BINARY_DIR}/src/levmar.h
  DESTINATION ${INSTALL_INCLUDEDIR}
  )

# install binaries
install(
  TARGETS ${INSTALL_TARGETS}
  EXPORT ${targets_export_name}
  ARCHIVE DESTINATION ${INSTALL_LIBDIR}
  LIBRARY DESTINATION ${INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

# install cmake config and targets
install(
  FILES ${project_config_dst} ${version_config}
  DESTINATION ${INSTALL_CMAKE_DIR})

install(
  EXPORT ${targets_export_name}
  DESTINATION ${INSTALL_CMAKE_DIR}
  NAMESPACE levmar::)

# install pkgconfig
install(
  FILES ${CMAKE_BINARY_DIR}/src/${PROJECT_NAME}.pc
  DESTINATION "${INSTALL_PKGCONFIG_DIR}")

add_custom_target(check
	COMMAND echo "make check : Nothing really necessary to do"
        COMMAND echo "Just make sure target all is run"
        COMMAND make all
	COMMENT "Just to make upstream check target happy.")

### Uninstall target
# https://gitlab.kitware.com/cmake/community/-/wikis/FAQ#can-i-do-make-uninstall-with-cmake
if(NOT TARGET uninstall)
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

  add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake)
endif()
