# Copyright (c) 2015 Thomas Heller
# Copyright (c) 2015 Dominic Marcello
# Copyright (c) 2018-2019 Parsa Amini
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

cmake_minimum_required(VERSION 3.12.4)
cmake_policy(VERSION 3.12.4)

# Libraries linked via full path no longer produce linker search paths.
cmake_policy(SET CMP0003 NEW)
cmake_policy(SET CMP0074 NEW)

project(meshless CXX)

# Search path for CMake modules to be loaded by include() and find_package()
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

find_package(HPX REQUIRED NO_CMAKE_PACKAGE_REGISTRY)
find_package(Silo REQUIRED)
find_package(Boost REQUIRED)

set(source_files
	src/delaunay.cpp
	src/main.cpp
    src/math.cpp
	src/particle.cpp
	src/tree.cpp
	src/silo.cpp
	src/state.cpp
)

set(header_files
	src/dim.hpp  
	src/delaunay.hpp
	src/math.hpp  
	src/particle.hpp  
	src/real.hpp  
	src/silo.hpp
	src/tree.hpp  
	src/vect.hpp
)

add_hpx_executable(
  meshless
  DEPENDENCIES
    Silo::silo Boost::boost
  SOURCES
    ${source_files}
  HEADERS
    ${header_files}
)



