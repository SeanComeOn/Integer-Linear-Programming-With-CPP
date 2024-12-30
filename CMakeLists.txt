cmake_minimum_required(VERSION 3.10)
project(simplex-project)

set(CMAKE_CXX_STANDARD 14)

include_directories(include)

set(SOURCES
    src/main.cpp
    src/Simplex.cpp
)

add_executable(simplex ${SOURCES})

# # 设置编译选项
# set(CMAKE_CXX_FLAGS_DEBUG "-g")
# set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# # 设置默认构建类型为 Debug
# if(NOT CMAKE_BUILD_TYPE)
#     set(CMAKE_BUILD_TYPE Debug)
# endif()
