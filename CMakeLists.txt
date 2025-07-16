cmake_minimum_required(VERSION 3.10)

project(suffix)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")   # 对C代码生效
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2") # 对C++代码生效 

add_compile_options(-g)

# 添加头文件包含目录
include_directories(
  "${PROJECT_SOURCE_DIR}/src"
  "${PROJECT_SOURCE_DIR}/src/edlib/include"  # 添加 edlib 的头文件目录
)

# 添加edlib库的源文件
set(EDLIB_SRC "${PROJECT_SOURCE_DIR}/src/edlib/src/edlib.cpp")

file(GLOB SRC_FILES
  "${PROJECT_SOURCE_DIR}/src/*.h"
  "${PROJECT_SOURCE_DIR}/src/*.cc"
  "${PROJECT_SOURCE_DIR}/src/*.cpp"
  "${PROJECT_SOURCE_DIR}/src/*.hpp"
  "${PROJECT_SOURCE_DIR}/src/ksw2/*.c"
  "${PROJECT_SOURCE_DIR}/src/ksw2/*.h"
)

# 构建edlib库
add_library(edlib STATIC ${EDLIB_SRC})

add_subdirectory("${PROJECT_SOURCE_DIR}/src/spoa")

add_executable(${CMAKE_PROJECT_NAME} ${SRC_FILES})

# 链接edlib库
target_link_libraries(${CMAKE_PROJECT_NAME} "spoa" edlib)