cmake_minimum_required(VERSION 2.6)
project(SCOREMPI)

#default compiler is intel
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -openmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -openmp") 

# further optional arguments
# -Wno-deprecated-declarations -Wno-write-strings")  #"-pg" option for debugging

#change compiler flags if gcc/g++ compiler is utilized
if (${CMAKE_COMPILER_IS_GNUCXX})
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -fopenmp")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -fopenmp")
endif()

if (${CMAKE_COMPILER_IS_GNUCC})
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -fopenmp")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -fopenmp")
endif()

#MPI is required so find the corresponding libraries
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

#which source files require compilation
add_executable(scorempi
			../src/lodepng.cpp
			../src/SCORE_Definitions.cpp
			../src/SCORE_Io.cpp
			../src/SCORE_Main.cpp
			../src/SCORE_Math.cpp
			../src/SCORE_Random.cpp
			../src/SCORE_Kernel.cpp
)

#linking process with additional mounting of MPI libraries
target_link_libraries(scorempi ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(scorempi PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(scorempi PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()