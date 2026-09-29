cmake_minimum_required(VERSION 3.20)
project(pramana VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(MSVC)
  add_compile_options(/W3 /permissive- /utf-8 /bigobj /EHsc /Zc:__cplusplus)
  add_compile_definitions(_CRT_SECURE_NO_WARNINGS NOMINMAX WIN32_LEAN_AND_MEAN)
  # /fp:precise is the default: no FMA contraction, IEEE round-to-nearest.
  # The certifier's outward rounding relies on this.
  set(CMAKE_CXX_FLAGS_RELEASE "/O2 /DNDEBUG")
else()
  add_compile_options(-Wall -Wextra -ffp-contract=off)
  set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")
endif()

# ---------------------------------------------------------------------------
# Solver core: every algorithm lives here. The ONLY external dependencies are
# the C++ standard library and the OS loader (the CUDA driver/NVRTC are loaded
# dynamically at runtime, never linked). See docs/FROM_SCRATCH.md.
# ---------------------------------------------------------------------------
file(GLOB_RECURSE PRAMANA_CORE_SOURCES CONFIGURE_DEPENDS
  ${CMAKE_SOURCE_DIR}/src/util/*.cpp
  ${CMAKE_SOURCE_DIR}/src/core/*.cpp
  ${CMAKE_SOURCE_DIR}/src/io/*.cpp
  ${CMAKE_SOURCE_DIR}/src/presolve/*.cpp
  ${CMAKE_SOURCE_DIR}/src/lp/*.cpp
  ${CMAKE_SOURCE_DIR}/src/ipm/*.cpp
  ${CMAKE_SOURCE_DIR}/src/pdhg/*.cpp
  ${CMAKE_SOURCE_DIR}/src/gpu/*.cpp
  ${CMAKE_SOURCE_DIR}/src/mip/*.cpp
  ${CMAKE_SOURCE_DIR}/src/cert/*.cpp
  ${CMAKE_SOURCE_DIR}/src/router/*.cpp
  ${CMAKE_SOURCE_DIR}/src/parametric/*.cpp
  ${CMAKE_SOURCE_DIR}/src/api/solver.cpp
)
add_library(pramana_core STATIC ${PRAMANA_CORE_SOURCES})
target_include_directories(pramana_core PUBLIC ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/include)
set_target_properties(pramana_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
if(NOT WIN32)
  find_package(Threads REQUIRED)
  target_link_libraries(pramana_core PUBLIC Threads::Threads ${CMAKE_DL_LIBS})
endif()

# C API shared library (used by the Python bindings through ctypes).
add_library(pramana SHARED ${CMAKE_SOURCE_DIR}/src/api/c_api.cpp)
target_link_libraries(pramana PRIVATE pramana_core)
target_compile_definitions(pramana PRIVATE PRAMANA_BUILDING_DLL)
set_target_properties(pramana PROPERTIES OUTPUT_NAME pramana)

# Command line interface.
add_executable(pramana_cli ${CMAKE_SOURCE_DIR}/src/cli/main.cpp)
target_link_libraries(pramana_cli PRIVATE pramana_core)
set_target_properties(pramana_cli PROPERTIES OUTPUT_NAME pramana)
if(WIN32)
  # Executable and DLL would both be pramana.*; keep the import lib names distinct.
  set_target_properties(pramana_cli PROPERTIES OUTPUT_NAME pramana)
  set_target_properties(pramana PROPERTIES OUTPUT_NAME libpramana)
endif()

# Unit tests (own minimal framework, no third-party test library).
file(GLOB PRAMANA_TEST_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/tests/cpp/*.cpp)
add_executable(pramana_tests ${PRAMANA_TEST_SOURCES})
target_link_libraries(pramana_tests PRIVATE pramana_core)
target_compile_definitions(pramana_tests PRIVATE PRAMANA_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

enable_testing()
add_test(NAME unit COMMAND pramana_tests)
