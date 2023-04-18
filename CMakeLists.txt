cmake_minimum_required(VERSION 3.0.1)

project(GeneRaxCoreProject
        VERSION 0.0)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer")

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG(-std=c++14 CXX_STANDARD_14_SUPPORT -g)
if(CXX_STANDARD_14_SUPPORT)
    set(CMAKE_CXX_STANDARD 14)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
else()
    message(FATAL_ERROR "Compiler ${CMAKE_CXX_COMPILER} has no C++14 support.")
endif()

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE
      "Debug"
            CACHE STRING "Build type" FORCE
            )
          message(STATUS "CMAKE_BUILD_TYPE not set, defaulting to Debug")
endif ()


add_subdirectory(ext)

if (NOT DISABLE_MPI)
  include(FindMPI)
  if(MPI_FOUND)
    set(WITH_MPI TRUE)
    add_definitions(-DWITH_MPI)
    message(STATUS "WITH MPI ${OpenMP_EXE_LINKER_FLAGS}")
  endif()
endif()


set(GENERAXCORE_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src/
   ${CMAKE_CURRENT_SOURCE_DIR}/ext
   ${MPI_INCLUDE_PATH}
   ${CORAX_INCLUDE_DIRS}
   ${MPISCHEDULER_INCLUDE_DIRS}
   CACHE INTERNAL "${PROJECT_NAME}: Include dirs"
  )
  

add_subdirectory(src)

