cmake_minimum_required(VERSION 3.10)

project(MyProject)

include_directories(handle_file.h sam.h sequence_handler.h utils.h)

add_executable(Myproject main.cc sam.cc handle_file.cc sam.cc sequence_handler.cc utils.cc)