# ~/c/BBBrem/CMakeLists.txt

cmake_minimum_required (VERSION 2.8)

project(BBBrem)
SET(CMAKE_C_COMPILER /opt/local/bin/gcc)
SET(CMAKE_CXX_COMPILER /opt/local/bin/g++)
set(CMAKE_CXX_STANDARD 11)
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11"   CHECK_CXX_FLAG)
message("c++11    CHECK_CXX_FLAG    = ${CHECK_CXX_FLAG}")
#SET(CMAKE_OSX_SYSROOT /opt/local)
#SET(CMAKE_PREFIX_PATH /opt/local)

message("PROJECT_NAME = ${PROJECT_NAME} CMAKE_SOURCE_DIR=${CMAKE_SOURCE_DIR} CMAKE_BINARY_DIR=${CMAKE_BINARY_DIR} CMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")

cmake_policy(SET CMP0008 NEW) # Libraries linked by full-path must have a valid library file name.

#source files
file (GLOB my_src_files   RELATIVE ${CMAKE_SOURCE_DIR}     "*.C" )    # all .C
file (GLOB my_headers     RELATIVE ${CMAKE_SOURCE_DIR}     "*.h" )    # all .h
file (GLOB my_txt         RELATIVE ${CMAKE_SOURCE_DIR}     "*.txt" )  # all .txt, only shown, not compiled

message("my_src_files = ${my_src_files}")
message("my_headers   = ${my_headers}")
message("my_txt       = ${my_txt}")

add_executable (BBBrem ${my_src_files} ${my_headers} ${my_txt} )

target_link_libraries (BBBrem )

install(TARGETS BBBrem RUNTIME DESTINATION ${PROJECT_NAME}/${CMAKE_BUILD_TYPE} )


