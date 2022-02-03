cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

project(json LANGUAGES CXX)

# ISO standard compliance (C++14)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Optimization
set(CMAKE_CXX_FLAGS "-O2")

# Compiler warnings
if (MSVC)
    # warning level 4
    add_compile_options(/W4)
else()
    # lots of warnings
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# external libraries
find_library(jansson jansson)
find_library(curl curl)

add_executable(json json.cc)

target_link_libraries(json
			PUBLIC jansson
		    PUBLIC curl
		  )

# INSTALL target (make install)
install(TARGETS json)

