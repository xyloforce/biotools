cmake_minimum_required(VERSION 3.0)

project(bioTools)

add_library(bio_entry bio_entry.cpp)
add_library(bio_file bio_file.cpp bio_entry.cpp)
add_library(bed_tools bed_tools.cpp bio_file.cpp bio_entry.cpp)
add_library(wiggle_tools wiggle_tools.cpp bio_file.cpp bio_entry.cpp)
add_library(tools tools.cpp)

add_executable(countFeatures.bin countFeatures.cpp)
add_executable(intersectKeepingNames.bin intersectKeepingNames.cpp)

target_link_libraries(countFeatures.bin PUBLIC bed_tools tools)
target_link_libraries(intersectKeepingNames.bin PUBLIC bed_tools tools)
