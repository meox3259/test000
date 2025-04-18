cmake_minimum_required(VERSION 3.10)

project(suffix)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")   # 对C代码生效
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2") # 对C++代码生效 [[1]][[6]]

file(GLOB SRC_FILES
  "${PROJECT_SOURCE_DIR}/src/*.h"
  "${PROJECT_SOURCE_DIR}/src/*.cc"
  "${PROJECT_SOURCE_DIR}/src/*.cpp"
  "${PROJECT_SOURCE_DIR}/src/*.hpp"
  "${PROJECT_SOURCE_DIR}/src/ksw2/*.c"
  "${PROJECT_SOURCE_DIR}/src/ksw2/*.h"
  "${PROJECT_SOURCE_DIR}/src/edlib/include/*.h"
  "${PROJECT_SOURCE_DIR}/src/edlib/src/*.c"
)

add_subdirectory("${PROJECT_SOURCE_DIR}/src/spoa")

add_executable(${CMAKE_PROJECT_NAME} ${SRC_FILES})

target_link_libraries(${CMAKE_PROJECT_NAME} "spoa")