cmake_minimum_required(VERSION 3.0)

# project name
project (OTToolbox)

# C++11 is required to build this project
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_STANDARD 11)

# add src root as include directory
include_directories(.)

# Position Independent Code
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

# Misc compiler flags
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wall -Wextra -pedantic-errors")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra -pedantic-errors")


# Verbose mode
option(SET_VERBOSE "Verbose output" ON)
if( SET_VERBOSE )
	set (CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -DVERBOSE")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVERBOSE")
endif ( SET_VERBOSE )



# configure CPLEX
option(USE_CPLEX "Use the CPLEX Network Flow Solver as backend" OFF)

if( USE_CPLEX )
	set(CPLEX_LIBRARY /media/daten/lib/cplex12.6.1/cplex/lib/x86-64_linux/static_pic/libcplex.a pthread)
	set(CPLEX_INCLUDE_DIRECTORY "/media/daten/lib/cplex12.6.1/cplex/include")
	set(CPLEX_COMPILE_OPTIONS -DNDEBUG -DIL_STD -DILOSTRICTPOD)
	# define USE_CPLEX compiler flag
	set(CPLEX_COMPILE_OPTIONS ${CPLEX_COMPILE_OPTIONS} -DUSE_CPLEX)

	add_subdirectory (LP_CPLEX)

endif ( USE_CPLEX )



# configure Lemon
option(USE_LEMON "Use the Lemon Network Flow Solver as backend" OFF)

if( USE_LEMON )
	set(LEMON_LIBRARY "/media/daten/lib/lemon/libemon.a")
	set(LEMON_INCLUDE_DIRECTORY "/media/daten/lib/")
	set(LEMON_COMPILE_OPTIONS "")
	set(LEMON_COMPILE_OPTIONS ${LEMON_COMPILE_OPTIONS} -DUSE_LEMON)


	add_subdirectory (LP_Lemon)

endif ( USE_LEMON )


# configure lp_solve
option(USE_LPSOLVE "Use lp_solve as backend" OFF)

if( USE_LPSOLVE )
	set(LPSOLVE_LIBRARY "/media/daten/lib/lp_solve_5.5.2.5/liblpsolve55.a" "dl")
	set(LPSOLVE_INCLUDE_DIRECTORY "/media/daten/lib/lp_solve_5.5.2.5")
	set(LPSOLVE_COMPILE_OPTIONS "")
	set(LPSOLVE_COMPILE_OPTIONS ${LPSOLVE_COMPILE_OPTIONS} -DUSE_LPSOLVE)


	add_subdirectory (LP_lp_solve)

endif ( USE_LPSOLVE )


# configure Sinkhorn
option(USE_SINKHORN "Compile Sinkhorn solver" ON)
if( USE_SINKHORN )
	# configure eigen library
	set(EIGEN_LIBRARY "/usr/include/eigen3")
	set(SINKHORN_COMPILE_OPTIONS "")
	set(SINKHORN_INCLUDE_DIRECTORY ${EIGEN_LIBRARY})
	add_subdirectory (Sinkhorn)
endif ( USE_SINKHORN )


# find all project files
add_subdirectory (Common)
add_subdirectory (ShortCutSolver)
add_subdirectory (Examples)


install (FILES Common.h LP_CPLEX.h LP_Lemon.h ShortCutSolver.h Sinkhorn.h DESTINATION ${CMAKE_BINARY_DIR}/../include)
