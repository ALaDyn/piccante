cmake_minimum_required (VERSION 3.1)
project(piccante)
enable_language(CXX)

## use C++11 if you can, otherwise, add the NO_CXX11 definition to correctly build piccante (in this case, boost 1.47+ is required for RNG)
#set(CMAKE_CXX_STANDARD 11)
add_definitions(-DNO_CXX11)

if(MSVC)
  add_definitions(-D_SCL_SECURE_NO_WARNINGS)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
  add_definitions(-D_USE_MATH_DEFINES)
endif()

## use Boost if you can (mandatory if C++11 is not available); if found, a USE_BOOST flag is declared - remove definition manually if you enable C++11!
find_package(Boost 1.47.0 COMPONENTS system filesystem)
if(Boost_FOUND)
  add_definitions(-DUSE_BOOST)
endif()

find_package(MPI REQUIRED)

# make sure that the default is a RELEASE
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE)
endif()

set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS}")

include_directories(${MPI_INCLUDE_PATH} ${Boost_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/src/)

set (COMMON_PICCANTE_FILES 	
	${CMAKE_SOURCE_DIR}/src/current.cpp
	${CMAKE_SOURCE_DIR}/src/em_field.cpp
	${CMAKE_SOURCE_DIR}/src/grid.cpp
	${CMAKE_SOURCE_DIR}/src/jsoncpp.cpp
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

target_link_libraries(piccante                                 ${MPI_LIBRARIES} ${Boost_LIBRARIES}) 
target_link_libraries(piccante_dev                             ${MPI_LIBRARIES} ${Boost_LIBRARIES}) 

target_link_libraries(piccante_clusters                        ${MPI_LIBRARIES} ${Boost_LIBRARIES}) 
target_link_libraries(piccante_landau_damping                  ${MPI_LIBRARIES} ${Boost_LIBRARIES}) 
target_link_libraries(piccante_1d_TNSA                         ${MPI_LIBRARIES} ${Boost_LIBRARIES}) 
target_link_libraries(piccante_1d_solitons                     ${MPI_LIBRARIES} ${Boost_LIBRARIES}) 
target_link_libraries(piccante_1d_wakefield                    ${MPI_LIBRARIES} ${Boost_LIBRARIES}) 
target_link_libraries(piccante_2d_grating                      ${MPI_LIBRARIES} ${Boost_LIBRARIES}) 
target_link_libraries(piccante_2d_twostream                    ${MPI_LIBRARIES} ${Boost_LIBRARIES}) 
target_link_libraries(piccante_2d_twostream_no_output          ${MPI_LIBRARIES} ${Boost_LIBRARIES}) 
target_link_libraries(piccante_2d_twostream_scaling            ${MPI_LIBRARIES} ${Boost_LIBRARIES}) 
target_link_libraries(piccante_2d_twostream_scaling_no_noutput ${MPI_LIBRARIES} ${Boost_LIBRARIES}) 
target_link_libraries(piccante_2d_wakefield                    ${MPI_LIBRARIES} ${Boost_LIBRARIES}) 
target_link_libraries(piccante_2d_wakefield_no_output          ${MPI_LIBRARIES} ${Boost_LIBRARIES}) 
target_link_libraries(piccante_3d_twostream                    ${MPI_LIBRARIES} ${Boost_LIBRARIES}) 
target_link_libraries(piccante_3d_twostream_no_output          ${MPI_LIBRARIES} ${Boost_LIBRARIES}) 
target_link_libraries(piccante_3d_twostream_outputParamScan    ${MPI_LIBRARIES} ${Boost_LIBRARIES}) 
target_link_libraries(piccante_3d_twostream_scaling            ${MPI_LIBRARIES} ${Boost_LIBRARIES}) 
target_link_libraries(piccante_3d_twostream_scaling_no_output  ${MPI_LIBRARIES} ${Boost_LIBRARIES}) 
target_link_libraries(piccante_3d_wakefield                    ${MPI_LIBRARIES} ${Boost_LIBRARIES}) 
target_link_libraries(piccante_3d_wakefield_no_output          ${MPI_LIBRARIES} ${Boost_LIBRARIES}) 

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
