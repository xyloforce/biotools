cmake_minimum_required(VERSION 3.5)
set(CMAKE_CXX_STANDARD 17)

project(bioTools)

add_library(bio_entry bio_entry.cpp)
add_library(bio_file bio_file.cpp bio_entry.cpp)
add_library(bed_tools bed_tools.cpp bio_file.cpp bio_entry.cpp)
add_library(fasta_tools fasta_tools.cpp bio_file.cpp bio_entry.cpp)
add_library(vcf_tools vcf_tools.cpp bio_file.cpp bio_entry.cpp)
add_library(wiggle_tools wiggle_tools.cpp bio_file.cpp bio_entry.cpp)
add_library(tools tools.cpp)

add_executable(countFeatures.bin countFeatures.cpp)
add_executable(intersectKeepingNames.bin intersectKeepingNames.cpp)
add_executable(validateMuts.bin validateMuts.cpp)

target_link_libraries(countFeatures.bin PUBLIC bed_tools tools)
target_link_libraries(intersectKeepingNames.bin PUBLIC vcf_tools bed_tools tools)
target_link_libraries(validateMuts.bin PUBLIC vcf_tools fasta_tools tools)
