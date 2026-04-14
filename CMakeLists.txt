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
  XMLIO
)

# --------------------------------------------------
# Executable
# --------------------------------------------------
add_executable(own-normFactors
  src/own-normFactors.cxx
  src/normFunctions.cxx
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
    ROOT::XMLIO
    OpenMP::OpenMP_CXX
)

# --------------------------------------------------
# Tests (no ROOT dependency)
# --------------------------------------------------
enable_testing()

add_executable(test_indexing_32x16x2
  tests/test_indexing_32x16x2.cxx
)
target_compile_options(test_indexing_32x16x2 PRIVATE -Wall -Wextra)
add_test(NAME indexing_32x16x2 COMMAND test_indexing_32x16x2)
