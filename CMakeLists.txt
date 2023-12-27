cmake_minimum_required(VERSION 3.15.0)
project(GD)
set(CMAKE_CXX_STANDARD 17)
add_subdirectory(pybind11)
pybind11_add_module(GD 
	main.cpp 
	src/feature.cpp
	src/raster.cpp
	src/sampler.cpp
	src/image_fusion.cpp
	src/image_segment.cpp
	src/rand_Engine.cpp
	include/feature.h
	include/raster.h 
	include/feature.h 
	include/sampler.h 
	include/image_fusion.h 
	include/image_segment.h 
	include/templates.hpp
	include/pythonAPI.hpp

)
