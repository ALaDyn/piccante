cmake_minimum_required (VERSION 3.1)
project(piccante)
enable_language(CXX)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

## use C++11 if you can, otherwise, add the NO_CXX11 definition to correctly build piccante (in this case, boost 1.47+ is required for RNG)
set(CMAKE_CXX_STANDARD 11)
#add_definitions(-DNO_CXX11)

# make sure that the default is a RELEASE
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE)
endif()

if(MSVC)
  add_definitions(-D_SCL_SECURE_NO_WARNINGS)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
  add_definitions(-D_USE_MATH_DEFINES)
endif()

include(ExternalProject)

## use Boost if you can (mandatory if C++11 is not available); if found, a USE_BOOST flag is declared - remove definition manually if you enable C++11!
find_package(Boost 1.47.0 COMPONENTS system filesystem)
if(Boost_FOUND)
  add_definitions(-DUSE_BOOST)
  include_directories(${Boost_INCLUDE_DIR})
endif()
if(NOT Boost_FOUND AND NO_CXX11)
  message(FATAL_ERROR "Without C++11 available, the system must provide Boost Libraries 1.47+")
endif()

find_package(MKL)
if(NOT MKL_FOUND)
  find_package(FFTW)
  if(FFTW_FOUND)
    include_directories(${FFTW_INCLUDES})
  endif()
else()
  include_directories(${MKL_INCLUDE_DIRS})
endif()
if(FFTW_MPI_FOUND OR MKL_FOUND)
  add_definitions(-D_USE_FFTW_FILTER)
endif()

find_package(HDF5)
if(HDF5_FOUND)
#  DO NOT ENABLE, since it is broken
#  add_definitions(-DUSE_HDF5)
#  include_directories(${HDF5_INCLUDE_DIRS})
endif()

find_package(JsonCpp)
if(NOT JsonCpp_FOUND)
	message(STATUS "  --> JsonCpp will be built as an external project")
	ExternalProject_Add(json-cpp
		URL https://github.com/open-source-parsers/jsoncpp/archive/1.8.1.zip
		CMAKE_ARGS
			-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
			-DJSONCPP_WITH_TESTS:BOOL=OFF
		# Disable the install step.
		INSTALL_COMMAND ""
		# Wrap the download, configure and build steps in a script to log the output.
		LOG_DOWNLOAD ON
		LOG_CONFIGURE ON
		LOG_BUILD ON
	)
	ExternalProject_Get_Property(json-cpp source_dir)
	ExternalProject_Get_Property(json-cpp binary_dir)
	set(JsonCpp_INCLUDE_DIR "${source_dir}/include")
	set(JsonCpp_LIBRARY_DIR "${binary_dir}/src/lib_json")
	set(JsonCpp_LIBRARY "jsoncpp")
endif()

find_package(MPI REQUIRED)

file(DOWNLOAD "https://raw.githubusercontent.com/cenit/jburkardt/master/sobol/sobol.cpp"
  ${CMAKE_SOURCE_DIR}/src/sobol.cpp
  TIMEOUT 3
  STATUS download_status
)
file(DOWNLOAD "https://raw.githubusercontent.com/cenit/jburkardt/master/sobol/sobol.hpp"
  ${CMAKE_SOURCE_DIR}/src/sobol.hpp
  TIMEOUT 3
  STATUS download_status
)
execute_process(
  COMMAND git --work-tree=. --git-dir=.git apply ${CMAKE_SOURCE_DIR}/cmake/sobol.patch --ignore-whitespace --whitespace=nowarn --verbose
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  RESULT_VARIABLE error_code
)

set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS}")

include_directories(${JsonCpp_INCLUDE_DIR} ${jsoncpp_DIR} ${MPI_INCLUDE_PATH} ${CMAKE_SOURCE_DIR}/src/)
link_directories(${JsonCpp_LIBRARY_DIR})

set (COMMON_PICCANTE_FILES 	
	${CMAKE_SOURCE_DIR}/src/current.cpp
	${CMAKE_SOURCE_DIR}/src/em_field.cpp
	${CMAKE_SOURCE_DIR}/src/grid.cpp
	${CMAKE_SOURCE_DIR}/src/jsonparser.cpp
	${CMAKE_SOURCE_DIR}/src/output_manager.cpp
	${CMAKE_SOURCE_DIR}/src/particle_species.cpp
	${CMAKE_SOURCE_DIR}/src/sobol.cpp
	${CMAKE_SOURCE_DIR}/src/structures.cpp
	${CMAKE_SOURCE_DIR}/src/utilities.cpp
)

add_executable (piccante                                 ${CMAKE_SOURCE_DIR}/src/main-piccante.cpp ${COMMON_PICCANTE_FILES})
add_executable (piccante_dev                             ${CMAKE_SOURCE_DIR}/src/main-devel.cpp ${COMMON_PICCANTE_FILES})

add_executable (piccante_clusters                        ${CMAKE_SOURCE_DIR}/examples/MAINS/clusters.cpp ${COMMON_PICCANTE_FILES})
add_executable (piccante_landau_damping                  ${CMAKE_SOURCE_DIR}/examples/MAINS/landau_damping.cpp ${COMMON_PICCANTE_FILES})
add_executable (piccante_1d_TNSA                         ${CMAKE_SOURCE_DIR}/examples/MAINS/1d_TNSA.cpp ${COMMON_PICCANTE_FILES})
add_executable (piccante_1d_solitons                     ${CMAKE_SOURCE_DIR}/examples/MAINS/1d_solitons.cpp ${COMMON_PICCANTE_FILES})
add_executable (piccante_1d_wakefield                    ${CMAKE_SOURCE_DIR}/examples/MAINS/1d_wakefield.cpp ${COMMON_PICCANTE_FILES})
add_executable (piccante_2d_grating                      ${CMAKE_SOURCE_DIR}/examples/MAINS/2d_grating.cpp ${COMMON_PICCANTE_FILES})
add_executable (piccante_2d_twostream                    ${CMAKE_SOURCE_DIR}/examples/MAINS/2d_twostream.cpp ${COMMON_PICCANTE_FILES})
add_executable (piccante_2d_twostream_no_output          ${CMAKE_SOURCE_DIR}/examples/MAINS/2d_twostream_no_output.cpp ${COMMON_PICCANTE_FILES})
add_executable (piccante_2d_twostream_scaling            ${CMAKE_SOURCE_DIR}/examples/MAINS/2d_twostream_scaling.cpp ${COMMON_PICCANTE_FILES})
add_executable (piccante_2d_twostream_scaling_no_noutput ${CMAKE_SOURCE_DIR}/examples/MAINS/2d_twostream_scaling_no_noutput.cpp ${COMMON_PICCANTE_FILES})
add_executable (piccante_2d_wakefield                    ${CMAKE_SOURCE_DIR}/examples/MAINS/2d_wakefield.cpp ${COMMON_PICCANTE_FILES})
add_executable (piccante_2d_wakefield_no_output          ${CMAKE_SOURCE_DIR}/examples/MAINS/2d_wakefield_no_output.cpp ${COMMON_PICCANTE_FILES})
add_executable (piccante_3d_twostream                    ${CMAKE_SOURCE_DIR}/examples/MAINS/3d_twostream.cpp ${COMMON_PICCANTE_FILES})
add_executable (piccante_3d_twostream_no_output          ${CMAKE_SOURCE_DIR}/examples/MAINS/3d_twostream_no_output.cpp ${COMMON_PICCANTE_FILES})
add_executable (piccante_3d_twostream_outputParamScan    ${CMAKE_SOURCE_DIR}/examples/MAINS/3d_twostream_outputParamScan.cpp ${COMMON_PICCANTE_FILES})
add_executable (piccante_3d_twostream_scaling            ${CMAKE_SOURCE_DIR}/examples/MAINS/3d_twostream_scaling.cpp ${COMMON_PICCANTE_FILES})
add_executable (piccante_3d_twostream_scaling_no_output  ${CMAKE_SOURCE_DIR}/examples/MAINS/3d_twostream_scaling_no_output.cpp ${COMMON_PICCANTE_FILES})
add_executable (piccante_3d_wakefield                    ${CMAKE_SOURCE_DIR}/examples/MAINS/3d_wakefield.cpp ${COMMON_PICCANTE_FILES})
add_executable (piccante_3d_wakefield_no_output          ${CMAKE_SOURCE_DIR}/examples/MAINS/3d_wakefield_no_output.cpp ${COMMON_PICCANTE_FILES})

if(NOT JsonCpp_FOUND)
  add_dependencies(piccante                                  json-cpp)
  add_dependencies(piccante_dev                              json-cpp)
  add_dependencies(piccante_clusters                         json-cpp)
  add_dependencies(piccante_landau_damping                   json-cpp)
  add_dependencies(piccante_1d_TNSA                          json-cpp)
  add_dependencies(piccante_1d_solitons                      json-cpp)
  add_dependencies(piccante_1d_wakefield                     json-cpp)
  add_dependencies(piccante_2d_grating                       json-cpp)
  add_dependencies(piccante_2d_twostream                     json-cpp)
  add_dependencies(piccante_2d_twostream_no_output           json-cpp)
  add_dependencies(piccante_2d_twostream_scaling             json-cpp)
  add_dependencies(piccante_2d_twostream_scaling_no_noutput  json-cpp)
  add_dependencies(piccante_2d_wakefield                     json-cpp)
  add_dependencies(piccante_2d_wakefield_no_output           json-cpp)
  add_dependencies(piccante_3d_twostream                     json-cpp)
  add_dependencies(piccante_3d_twostream_no_output           json-cpp)
  add_dependencies(piccante_3d_twostream_outputParamScan     json-cpp)
  add_dependencies(piccante_3d_twostream_scaling             json-cpp)
  add_dependencies(piccante_3d_twostream_scaling_no_output   json-cpp)
  add_dependencies(piccante_3d_wakefield                     json-cpp)
  add_dependencies(piccante_3d_wakefield_no_output           json-cpp)
endif()

target_link_libraries(piccante                                 ${MPI_LIBRARIES} ${JsonCpp_LIBRARY}) 
target_link_libraries(piccante_dev                             ${MPI_LIBRARIES} ${JsonCpp_LIBRARY}) 
target_link_libraries(piccante_clusters                        ${MPI_LIBRARIES} ${JsonCpp_LIBRARY}) 
target_link_libraries(piccante_landau_damping                  ${MPI_LIBRARIES} ${JsonCpp_LIBRARY}) 
target_link_libraries(piccante_1d_TNSA                         ${MPI_LIBRARIES} ${JsonCpp_LIBRARY}) 
target_link_libraries(piccante_1d_solitons                     ${MPI_LIBRARIES} ${JsonCpp_LIBRARY}) 
target_link_libraries(piccante_1d_wakefield                    ${MPI_LIBRARIES} ${JsonCpp_LIBRARY}) 
target_link_libraries(piccante_2d_grating                      ${MPI_LIBRARIES} ${JsonCpp_LIBRARY}) 
target_link_libraries(piccante_2d_twostream                    ${MPI_LIBRARIES} ${JsonCpp_LIBRARY}) 
target_link_libraries(piccante_2d_twostream_no_output          ${MPI_LIBRARIES} ${JsonCpp_LIBRARY}) 
target_link_libraries(piccante_2d_twostream_scaling            ${MPI_LIBRARIES} ${JsonCpp_LIBRARY}) 
target_link_libraries(piccante_2d_twostream_scaling_no_noutput ${MPI_LIBRARIES} ${JsonCpp_LIBRARY}) 
target_link_libraries(piccante_2d_wakefield                    ${MPI_LIBRARIES} ${JsonCpp_LIBRARY}) 
target_link_libraries(piccante_2d_wakefield_no_output          ${MPI_LIBRARIES} ${JsonCpp_LIBRARY}) 
target_link_libraries(piccante_3d_twostream                    ${MPI_LIBRARIES} ${JsonCpp_LIBRARY}) 
target_link_libraries(piccante_3d_twostream_no_output          ${MPI_LIBRARIES} ${JsonCpp_LIBRARY}) 
target_link_libraries(piccante_3d_twostream_outputParamScan    ${MPI_LIBRARIES} ${JsonCpp_LIBRARY}) 
target_link_libraries(piccante_3d_twostream_scaling            ${MPI_LIBRARIES} ${JsonCpp_LIBRARY}) 
target_link_libraries(piccante_3d_twostream_scaling_no_output  ${MPI_LIBRARIES} ${JsonCpp_LIBRARY}) 
target_link_libraries(piccante_3d_wakefield                    ${MPI_LIBRARIES} ${JsonCpp_LIBRARY}) 
target_link_libraries(piccante_3d_wakefield_no_output          ${MPI_LIBRARIES} ${JsonCpp_LIBRARY}) 

if(Boost_FOUND)
  target_link_libraries(piccante                                 ${Boost_LIBRARIES}) 
  target_link_libraries(piccante_dev                             ${Boost_LIBRARIES}) 
  target_link_libraries(piccante_clusters                        ${Boost_LIBRARIES}) 
  target_link_libraries(piccante_landau_damping                  ${Boost_LIBRARIES}) 
  target_link_libraries(piccante_1d_TNSA                         ${Boost_LIBRARIES}) 
  target_link_libraries(piccante_1d_solitons                     ${Boost_LIBRARIES}) 
  target_link_libraries(piccante_1d_wakefield                    ${Boost_LIBRARIES}) 
  target_link_libraries(piccante_2d_grating                      ${Boost_LIBRARIES}) 
  target_link_libraries(piccante_2d_twostream                    ${Boost_LIBRARIES}) 
  target_link_libraries(piccante_2d_twostream_no_output          ${Boost_LIBRARIES}) 
  target_link_libraries(piccante_2d_twostream_scaling            ${Boost_LIBRARIES}) 
  target_link_libraries(piccante_2d_twostream_scaling_no_noutput ${Boost_LIBRARIES}) 
  target_link_libraries(piccante_2d_wakefield                    ${Boost_LIBRARIES}) 
  target_link_libraries(piccante_2d_wakefield_no_output          ${Boost_LIBRARIES}) 
  target_link_libraries(piccante_3d_twostream                    ${Boost_LIBRARIES}) 
  target_link_libraries(piccante_3d_twostream_no_output          ${Boost_LIBRARIES}) 
  target_link_libraries(piccante_3d_twostream_outputParamScan    ${Boost_LIBRARIES}) 
  target_link_libraries(piccante_3d_twostream_scaling            ${Boost_LIBRARIES}) 
  target_link_libraries(piccante_3d_twostream_scaling_no_output  ${Boost_LIBRARIES}) 
  target_link_libraries(piccante_3d_wakefield                    ${Boost_LIBRARIES}) 
  target_link_libraries(piccante_3d_wakefield_no_output          ${Boost_LIBRARIES}) 
endif()

if(HDF5_FOUND)
  target_link_libraries(piccante                                 ${HDF5_CXX_LIBRARIES}) 
  target_link_libraries(piccante_dev                             ${HDF5_CXX_LIBRARIES}) 
  target_link_libraries(piccante_clusters                        ${HDF5_CXX_LIBRARIES}) 
  target_link_libraries(piccante_landau_damping                  ${HDF5_CXX_LIBRARIES}) 
  target_link_libraries(piccante_1d_TNSA                         ${HDF5_CXX_LIBRARIES}) 
  target_link_libraries(piccante_1d_solitons                     ${HDF5_CXX_LIBRARIES}) 
  target_link_libraries(piccante_1d_wakefield                    ${HDF5_CXX_LIBRARIES}) 
  target_link_libraries(piccante_2d_grating                      ${HDF5_CXX_LIBRARIES}) 
  target_link_libraries(piccante_2d_twostream                    ${HDF5_CXX_LIBRARIES}) 
  target_link_libraries(piccante_2d_twostream_no_output          ${HDF5_CXX_LIBRARIES}) 
  target_link_libraries(piccante_2d_twostream_scaling            ${HDF5_CXX_LIBRARIES}) 
  target_link_libraries(piccante_2d_twostream_scaling_no_noutput ${HDF5_CXX_LIBRARIES}) 
  target_link_libraries(piccante_2d_wakefield                    ${HDF5_CXX_LIBRARIES}) 
  target_link_libraries(piccante_2d_wakefield_no_output          ${HDF5_CXX_LIBRARIES}) 
  target_link_libraries(piccante_3d_twostream                    ${HDF5_CXX_LIBRARIES}) 
  target_link_libraries(piccante_3d_twostream_no_output          ${HDF5_CXX_LIBRARIES}) 
  target_link_libraries(piccante_3d_twostream_outputParamScan    ${HDF5_CXX_LIBRARIES}) 
  target_link_libraries(piccante_3d_twostream_scaling            ${HDF5_CXX_LIBRARIES}) 
  target_link_libraries(piccante_3d_twostream_scaling_no_output  ${HDF5_CXX_LIBRARIES}) 
  target_link_libraries(piccante_3d_wakefield                    ${HDF5_CXX_LIBRARIES}) 
  target_link_libraries(piccante_3d_wakefield_no_output          ${HDF5_CXX_LIBRARIES}) 
endif()

if (FFTW_FOUND AND FFTW_MPI_FOUND AND NOT MKL_FOUND)
  target_link_libraries(piccante                                 ${FFTW_LIBRARIES} ${FFTW_MPI_LIBRARIES}) 
  target_link_libraries(piccante_dev                             ${FFTW_LIBRARIES} ${FFTW_MPI_LIBRARIES}) 
  target_link_libraries(piccante_clusters                        ${FFTW_LIBRARIES} ${FFTW_MPI_LIBRARIES}) 
  target_link_libraries(piccante_landau_damping                  ${FFTW_LIBRARIES} ${FFTW_MPI_LIBRARIES}) 
  target_link_libraries(piccante_1d_TNSA                         ${FFTW_LIBRARIES} ${FFTW_MPI_LIBRARIES}) 
  target_link_libraries(piccante_1d_solitons                     ${FFTW_LIBRARIES} ${FFTW_MPI_LIBRARIES}) 
  target_link_libraries(piccante_1d_wakefield                    ${FFTW_LIBRARIES} ${FFTW_MPI_LIBRARIES}) 
  target_link_libraries(piccante_2d_grating                      ${FFTW_LIBRARIES} ${FFTW_MPI_LIBRARIES}) 
  target_link_libraries(piccante_2d_twostream                    ${FFTW_LIBRARIES} ${FFTW_MPI_LIBRARIES}) 
  target_link_libraries(piccante_2d_twostream_no_output          ${FFTW_LIBRARIES} ${FFTW_MPI_LIBRARIES}) 
  target_link_libraries(piccante_2d_twostream_scaling            ${FFTW_LIBRARIES} ${FFTW_MPI_LIBRARIES}) 
  target_link_libraries(piccante_2d_twostream_scaling_no_noutput ${FFTW_LIBRARIES} ${FFTW_MPI_LIBRARIES}) 
  target_link_libraries(piccante_2d_wakefield                    ${FFTW_LIBRARIES} ${FFTW_MPI_LIBRARIES}) 
  target_link_libraries(piccante_2d_wakefield_no_output          ${FFTW_LIBRARIES} ${FFTW_MPI_LIBRARIES}) 
  target_link_libraries(piccante_3d_twostream                    ${FFTW_LIBRARIES} ${FFTW_MPI_LIBRARIES}) 
  target_link_libraries(piccante_3d_twostream_no_output          ${FFTW_LIBRARIES} ${FFTW_MPI_LIBRARIES}) 
  target_link_libraries(piccante_3d_twostream_outputParamScan    ${FFTW_LIBRARIES} ${FFTW_MPI_LIBRARIES}) 
  target_link_libraries(piccante_3d_twostream_scaling            ${FFTW_LIBRARIES} ${FFTW_MPI_LIBRARIES}) 
  target_link_libraries(piccante_3d_twostream_scaling_no_output  ${FFTW_LIBRARIES} ${FFTW_MPI_LIBRARIES}) 
  target_link_libraries(piccante_3d_wakefield                    ${FFTW_LIBRARIES} ${FFTW_MPI_LIBRARIES}) 
  target_link_libraries(piccante_3d_wakefield_no_output          ${FFTW_LIBRARIES} ${FFTW_MPI_LIBRARIES}) 
endif()

if (MKL_FOUND AND NOT FFTW_FOUND AND NOT FFTW_MPI_FOUND)
  target_link_libraries(piccante                                 ${MKL_LIBRARIES}) 
  target_link_libraries(piccante_dev                             ${MKL_LIBRARIES}) 
  target_link_libraries(piccante_clusters                        ${MKL_LIBRARIES}) 
  target_link_libraries(piccante_landau_damping                  ${MKL_LIBRARIES}) 
  target_link_libraries(piccante_1d_TNSA                         ${MKL_LIBRARIES}) 
  target_link_libraries(piccante_1d_solitons                     ${MKL_LIBRARIES}) 
  target_link_libraries(piccante_1d_wakefield                    ${MKL_LIBRARIES}) 
  target_link_libraries(piccante_2d_grating                      ${MKL_LIBRARIES}) 
  target_link_libraries(piccante_2d_twostream                    ${MKL_LIBRARIES}) 
  target_link_libraries(piccante_2d_twostream_no_output          ${MKL_LIBRARIES}) 
  target_link_libraries(piccante_2d_twostream_scaling            ${MKL_LIBRARIES}) 
  target_link_libraries(piccante_2d_twostream_scaling_no_noutput ${MKL_LIBRARIES}) 
  target_link_libraries(piccante_2d_wakefield                    ${MKL_LIBRARIES}) 
  target_link_libraries(piccante_2d_wakefield_no_output          ${MKL_LIBRARIES}) 
  target_link_libraries(piccante_3d_twostream                    ${MKL_LIBRARIES}) 
  target_link_libraries(piccante_3d_twostream_no_output          ${MKL_LIBRARIES}) 
  target_link_libraries(piccante_3d_twostream_outputParamScan    ${MKL_LIBRARIES}) 
  target_link_libraries(piccante_3d_twostream_scaling            ${MKL_LIBRARIES}) 
  target_link_libraries(piccante_3d_twostream_scaling_no_output  ${MKL_LIBRARIES}) 
  target_link_libraries(piccante_3d_wakefield                    ${MKL_LIBRARIES}) 
  target_link_libraries(piccante_3d_wakefield_no_output          ${MKL_LIBRARIES}) 
endif()

install (TARGETS 
	piccante 
	piccante_dev 
	DESTINATION ${CMAKE_SOURCE_DIR}/bin/
)

install (TARGETS 
	piccante_clusters 
	piccante_landau_damping 
	piccante_1d_TNSA 
	piccante_1d_solitons 
	piccante_1d_wakefield 
	piccante_2d_grating 
	piccante_2d_twostream 
	piccante_2d_twostream_no_output 
	piccante_2d_twostream_scaling 
	piccante_2d_twostream_scaling_no_noutput 
	piccante_2d_wakefield 
	piccante_2d_wakefield_no_output 
	piccante_3d_twostream 
	piccante_3d_twostream_no_output 
	piccante_3d_twostream_outputParamScan 
	piccante_3d_twostream_scaling 
	piccante_3d_twostream_scaling_no_output 
	piccante_3d_wakefield 
	piccante_3d_wakefield_no_output 
	DESTINATION ${CMAKE_SOURCE_DIR}/bin/examples/
)
