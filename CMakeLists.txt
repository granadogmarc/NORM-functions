cmake_minimum_required(VERSION 3.15)

project(norm-functions LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --------------------------------------------------
# OpenMP
# --------------------------------------------------
find_package(OpenMP REQUIRED)

# --------------------------------------------------
# ROOT
# --------------------------------------------------
find_package(ROOT REQUIRED COMPONENTS
  Core
  RIO
  Tree
  Hist
  Graf
  Gpad
  Physics
)

# --------------------------------------------------
# Executable
# --------------------------------------------------
add_executable(own-normFactors
  src/own-normFactors.cxx
  src/normFunctions.cxx
  src/scannerGeometry.cxx
)

target_include_directories(own-normFactors
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${ROOT_INCLUDE_DIRS}
)

target_link_libraries(own-normFactors
  PRIVATE
    ROOT::Core
    ROOT::RIO
    ROOT::Tree
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad
    ROOT::Physics
    OpenMP::OpenMP_CXX
)

# --------------------------------------------------
# Tests (no ROOT dependency)
# --------------------------------------------------
enable_testing()

# CASToR .geom parser and crystal positions (reads the files in scanners/)
add_executable(test_scanner_geometry
  tests/test_scanner_geometry.cxx
  src/scannerGeometry.cxx
)
target_include_directories(test_scanner_geometry PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(test_scanner_geometry PRIVATE
  SCANNER_DIR="${CMAKE_CURRENT_SOURCE_DIR}/scanners")
target_compile_options(test_scanner_geometry PRIVATE -Wall -Wextra)
add_test(NAME scanner_geometry COMMAND test_scanner_geometry)

# Only built when the source file is present in the checkout
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_indexing_32x16x2.cxx)
  add_executable(test_indexing_32x16x2
    tests/test_indexing_32x16x2.cxx
  )
  target_compile_options(test_indexing_32x16x2 PRIVATE -Wall -Wextra)
  add_test(NAME indexing_32x16x2 COMMAND test_indexing_32x16x2)
endif()
