cmake_minimum_required(VERSION 2.6)

if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

set(projName "skullStripping")
project(${projName})

configure_file (
  "${PROJECT_SOURCE_DIR}/cmakeConfig.h.in"
  "${PROJECT_BINARY_DIR}/cmakeConfig.h"
  )
 
# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}")

option(USE_OMP "Use OMP in multi-atlas?" ON)
if (USE_OMP)
    find_package(OpenMP)
    if(OPENMP_FOUND)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}") 
    endif()
endif (USE_OMP)

################################################################################
################################
# include ITK, to avoid setting it everytime, do:
# set the environment variable ITK_DIR to the dir containing ITKConfig.cmake in .bashrc
#set(ITK_DIR $ENV{ITK_DIR})
find_package (ITK)
if (ITK_FOUND)
    include(${ITK_USE_FILE})
endif(ITK_FOUND) 

set(itkLibraries itkvnl ${ITK_LIBRARIES})

set(namicLibraries ${itkLibraries})

set(headerAndLibDir
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/lib)

include_directories(${headerAndLibDir})

include_directories(
    ${ANTS_SRC}/Examples/include
    ${ANTS_SRC}/ImageRegistration
    ${ANTS_SRC}/Utilities
    )
link_directories(${ANTS_BUILD}/lib/)

# regular multi-atlas segmentation based on ANTS registration
add_subdirectory(multiAtlasANTS)
