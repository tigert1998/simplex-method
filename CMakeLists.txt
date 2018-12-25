cmake_minimum_required(VERSION 3.1)
project(simplex-method)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(src_dir "src")
set(bin_dir "bin")
set(include_dir "include")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${bin_dir})

include_directories(${include_dir})
aux_source_directory(${src_dir} srcs)

add_executable("main" src/main.cc src/linear_program_solver.cc include/timer.h)
add_executable("demo" src/demo.cc src/linear_program_solver.cc include/timer.h)
