cmake_minimum_required (VERSION 3.1)
project(piccante)
enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)

find_package(Boost COMPONENTS system filesystem REQUIRED)
find_package(MPI REQUIRED)

# make sure that the default is a RELEASE
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE)
endif()

set (CMAKE_C_FLAGS_RELEASE "-O3")
set (CMAKE_CXX_FLAGS_RELEASE "-O3")
set (CMAKE_C_FLAGS_DEBUG "-O0 -g")
set (CMAKE_CXX_FLAGS_DEBUG "-O0 -g")

set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS}")

include_directories(${MPI_INCLUDE_PATH} ${Boost_INCLUDE_DIR})

add_executable (piccante 
	src/current.cpp
	src/em_field.cpp
	src/grid.cpp
	src/jsoncpp.cpp
	src/jsonparser.cpp
	src/main-piccante.cpp
	src/output_manager.cpp
	src/particle_species.cpp
	src/sobol.cpp
	src/structures.cpp
	src/utilities.cpp
)

add_executable (piccante_dev
	src/current.cpp
	src/em_field.cpp
	src/grid.cpp
	src/jsoncpp.cpp
	src/jsonparser.cpp
	src/main-devel.cpp
	src/output_manager.cpp
	src/particle_species.cpp
	src/sobol.cpp
	src/structures.cpp
	src/utilities.cpp
)


target_link_libraries(piccante 
	${MPI_LIBRARIES} 
	${Boost_LIBRARIES}
) 

target_link_libraries(piccante_dev
	${MPI_LIBRARIES} 
	${Boost_LIBRARIES}
) 

install (TARGETS piccante piccante_dev DESTINATION ${CMAKE_SOURCE_DIR}/bin/)
