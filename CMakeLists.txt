cmake_minimum_required (VERSION 3.20)

#define the project
project (szp
  VERSION 1.0.0
  DESCRIPTION "SZp Error Bounded Lossy Compressor"
  LANGUAGES C CXX
  )

#correct was to set a default build type
# https://blog.kitware.com/cmake-and-the-default-build-type/
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "No build type was set. Setting build type to ${default_build_type}.")
  set(CMAKE_BUILD_TYPE ${default_build_type} CACHE
    STRING "Choose the type to build" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()


set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
option(BUILD_SHARED_LIBS "build shared libraries over static libraries" ON)

include(CTest)
include(GNUInstallDirs)

find_package(OpenMP REQUIRED COMPONENTS C CXX)
find_package(PkgConfig)
pkg_search_module(ZSTD  IMPORTED_TARGET libzstd)


add_library(
    SZp
        SZp/src/szp.cc
        SZp/src/szp_CompressionToolkit.cc
        SZp/src/szp_TypeManager.cc
        SZp/src/szp_double.cc
        SZp/src/szp_float.cc
        SZp/src/szp_helper.cc
        SZp/src/szp_rw.cc
        SZp/src/szpd_double.cc
        SZp/src/szpd_float.cc
)
target_compile_features(SZp PUBLIC c_std_99 cxx_std_14)
target_link_libraries(SZp PUBLIC PkgConfig::ZSTD OpenMP::OpenMP_CXX)
target_include_directories(SZp
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/SZp/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDE_DIR}/SZp>
)


install(TARGETS SZp EXPORT SZpTargets)
install(EXPORT SZpTargets NAMESPACE SZp:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SZp/)
export(EXPORT SZpTargets FILE ${CMAKE_CURRENT_BINARY_DIR}/cmake/SZpTargets.cmake NAMESPACE SZp::)
include(CMakePackageConfigHelpers)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/SZpConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/SZpConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SZp
)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/SZpConfigVersion.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY AnyNewerVersion
)
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/SZpConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/SZpConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SZp
)
install(DIRECTORY SZp/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SZp)
