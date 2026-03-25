cmake_minimum_required(VERSION 3.10)
project(DistroMatch)

if(NOT CMAKE_BUILD_TYPE)
        message("Defaulting CMAKE_BUILD_TYPE to Release")
        set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# src-directory of the project
include_directories(
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/src/algorithm
    ${PROJECT_SOURCE_DIR}/src/datastructre
    ${PROJECT_SOURCE_DIR}/src/io
    ${PROJECT_SOURCE_DIR}/src/programs
    ${PROJECT_SOURCE_DIR}/src/util
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lpthread  -flto")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -D_MPI")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -D_MPI -g -g3 -march=native")

include(FetchContent)

FetchContent_Declare(
    argtable3
    GIT_REPOSITORY https://github.com/argtable/argtable3.git
    GIT_TAG        4368929685a6634efff59678dda5ba79e7c50c9a
)

FetchContent_Declare(
    absl
    GIT_REPOSITORY https://github.com/abseil/abseil-cpp
    GIT_TAG 20240116.2
)
set(ABSL_PROPAGATE_CXX_STD ON)
set(BUILD_TESTING OFF)

FetchContent_MakeAvailable(argtable3 absl)

set(LIBS
    argtable3
    absl::flat_hash_map)

set(DistroMatch_SRC
    ${PROJECT_SOURCE_DIR}/src/distroMatch.cpp
    )

add_executable(DistroMatch ${DistroMatch_SRC})
target_link_libraries(DistroMatch PRIVATE ${LIBS})
target_include_directories(DistroMatch PRIVATE ${CMAKE_BINARY_DIR}/_deps/argtable3-src/src/)

find_package(MPI)
if(MPI_CXX_FOUND)
    target_link_libraries(DistroMatch PUBLIC MPI::MPI_CXX)
endif()

message(STATUS "CXX Flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "Linker flags: ${CMAKE_EXE_LINKER_FLAGS}")

if (CMAKE_BUILD_TYPE STREQUAL Release)
    message (STATUS "Release build!")
    message (STATUS "Release CXX Flags: ${CMAKE_CXX_FLAGS_RELEASE}")
endif()

if (CMAKE_BUILD_TYPE STREQUAL Debug)
    message (STATUS "Debug build!")
    message (STATUS "Debug CXX Flags: ${CMAKE_CXX_FLAGS_DEBUG}")
endif()
