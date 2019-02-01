if( USE_SINKHORN )

	add_library(Sinkhorn STATIC TSinkhornSolver.cpp TSinkhornKernel.cpp TSinkhornSolverBarycenter.cpp)

	target_include_directories(Sinkhorn PUBLIC ${SINKHORN_INCLUDE_DIRECTORY})
	target_compile_options(Sinkhorn PUBLIC ${SINKHORN_COMPILE_OPTIONS})

	install (TARGETS Sinkhorn DESTINATION ${CMAKE_BINARY_DIR}/../bin)

	install (FILES
		TSinkhornKernel.h
		TSinkhornSolver.h
		TSinkhornSolverBarycenter.h
		DESTINATION ${CMAKE_BINARY_DIR}/../include/Sinkhorn)


endif( USE_SINKHORN )
