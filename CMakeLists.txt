cmake_minimum_required(VERSION 3.15)

project(norm-functions LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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
)
