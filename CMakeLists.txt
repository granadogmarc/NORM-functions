cmake_minimum_required(VERSION 3.0)
project(own-normFactors)


set(ROOT_DIR /home/granado/Applications/root/root-6.32.02/install)


# Locate ROOT and define a variable for the include directory
find_package(ROOT REQUIRED COMPONENTS MathCore Hist Tree)
message(STATUS "ROOT libraries: ${ROOT_LIBRARIES}")
include_directories(${ROOT_INCLUDE_DIRS})


link_directories(/home/granado/Applications/root/root-6.32.02/install/lib)
#add_executable(first_test first_test.cxx normFunctions.cxx ${CASToR_DATAFILES_SRC} ${CASToR_SRC})
add_executable( own-normFactors own-normFactors.cxx normFunctions.cxx ${CASToR_DATAFILES_SRC} ${CASToR_SRC})

#target_link_libraries(first_test ${ROOT_LIBRARIES})
target_link_libraries(own-normFactors ${ROOT_LIBRARIES})
#add_executable(createLUT createLUT.cxx ${CASToR_DATAFILES_SRC} ${CASToR_SRC})
# Link against ROOT libraries
#target_link_libraries(createLUT ${ROOT_LIBRARIES})
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})



