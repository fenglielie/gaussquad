cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project(gaussquad VERSION 1.0)
include(cmake/zero.cmake)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

zero_setup()
zero_check()

add_library(gaussquad INTERFACE)
target_include_directories(gaussquad INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_subdirectory(example)

include(CTest)
if(BUILD_TESTING)
    include(FetchContent)

    FetchContent_Declare(
        googletest
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        URL https://github.com/google/googletest/releases/download/v1.17.0/googletest-1.17.0.tar.gz
    )

    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    set(BUILD_GMOCK OFF CACHE BOOL "" FORCE) # Disable GMock

    FetchContent_MakeAvailable(googletest)

    add_subdirectory(test)
endif()
