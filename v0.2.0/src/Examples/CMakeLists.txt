####################################################################################################################
## ShortCut

## CPLEX
if( USE_CPLEX )

	add_executable(ShortCut_CPLEX
		ShortCut.cpp)

	target_compile_options(ShortCut_CPLEX PUBLIC ${CPLEX_COMPILE_OPTIONS})
	target_include_directories(ShortCut_CPLEX PUBLIC ${CPLEX_INCLUDE_DIRECTORY})
	target_link_libraries(ShortCut_CPLEX
		LP_CPLEX
		ShortCutSolver
		${CPLEX_LIBRARY}
		Common
		)

	install (TARGETS ShortCut_CPLEX DESTINATION ${CMAKE_BINARY_DIR}/../examples/)

	# dense example
	add_executable(Dense_LP_CPLEX
		Dense_LP.cpp)

	target_compile_options(Dense_LP_CPLEX PUBLIC ${CPLEX_COMPILE_OPTIONS})
	target_include_directories(Dense_LP_CPLEX PUBLIC ${CPLEX_INCLUDE_DIRECTORY})
	target_link_libraries(Dense_LP_CPLEX
		LP_CPLEX
		${CPLEX_LIBRARY}
		Common
		)

	install (TARGETS Dense_LP_CPLEX DESTINATION ${CMAKE_BINARY_DIR}/../examples/)


endif( USE_CPLEX )


## Lemon
if( USE_LEMON )
	add_executable(ShortCut_Lemon
		ShortCut.cpp)

	target_compile_options(ShortCut_Lemon PUBLIC ${LEMON_COMPILE_OPTIONS})
	target_include_directories(ShortCut_Lemon PUBLIC ${LEMON_INCLUDE_DIRECTORY})
	target_link_libraries(ShortCut_Lemon
		${LEMON_LIBRARY}
		ShortCutSolver
		LP_Lemon
		Common
		)

	install (TARGETS ShortCut_Lemon DESTINATION ${CMAKE_BINARY_DIR}/../examples/)


	# dense example
	add_executable(Dense_LP_Lemon
		Dense_LP.cpp)

	target_compile_options(Dense_LP_Lemon PUBLIC ${LEMON_COMPILE_OPTIONS})
	target_include_directories(Dense_LP_Lemon PUBLIC ${LEMON_INCLUDE_DIRECTORY})
	target_link_libraries(Dense_LP_Lemon
		${LEMON_LIBRARY}
		LP_Lemon
		Common
		)

	install (TARGETS Dense_LP_Lemon DESTINATION ${CMAKE_BINARY_DIR}/../examples/)


endif( USE_LEMON )


if( USE_LPSOLVE )

	add_executable(ShortCut_lp_solve
		ShortCut.cpp)

	target_compile_options(ShortCut_lp_solve PUBLIC ${LPSOLVE_COMPILE_OPTIONS})
	target_include_directories(ShortCut_lp_solve PUBLIC ${LPSOLVE_INCLUDE_DIRECTORY})
	target_link_libraries(ShortCut_lp_solve
		ShortCutSolver
		LP_lp_solve
		${LPSOLVE_LIBRARY}
		Common
		)

	install (TARGETS ShortCut_lp_solve DESTINATION ${CMAKE_BINARY_DIR}/../examples/)


	# dense example	
	add_executable(Dense_LP_lp_solve
		Dense_LP.cpp)

	target_compile_options(Dense_LP_lp_solve PUBLIC ${LPSOLVE_COMPILE_OPTIONS})
	target_include_directories(Dense_LP_lp_solve PUBLIC ${LPSOLVE_INCLUDE_DIRECTORY})
	target_link_libraries(Dense_LP_lp_solve
		${LPSOLVE_LIBRARY}
		LP_lp_solve
		Common
		)

	install (TARGETS Dense_LP_lp_solve DESTINATION ${CMAKE_BINARY_DIR}/../examples/)

endif( USE_LPSOLVE )




####################################################################################################################
## Sinkhorn
if ( USE_SINKHORN )
	add_executable(Sinkhorn_Standard
		Sinkhorn_Standard.cpp)

	target_include_directories(Sinkhorn_Standard PUBLIC ${SINKHORN_INCLUDE_DIRECTORY})
	target_compile_options(Sinkhorn_Standard PUBLIC ${SINKHORN_COMPILE_OPTIONS})
	target_link_libraries(Sinkhorn_Standard
		Sinkhorn
		Common
		)

	install (TARGETS Sinkhorn_Standard DESTINATION ${CMAKE_BINARY_DIR}/../examples/)

	## Sinkhorn-Barycenter
	add_executable(Sinkhorn_Barycenter
		Sinkhorn_Barycenter.cpp)

	target_include_directories(Sinkhorn_Barycenter PUBLIC ${SINKHORN_INCLUDE_DIRECTORY})
	target_compile_options(Sinkhorn_Barycenter PUBLIC ${SINKHORN_COMPILE_OPTIONS})
	target_link_libraries(Sinkhorn_Barycenter
		Sinkhorn
		Common
		)

	install (TARGETS Sinkhorn_Barycenter DESTINATION ${CMAKE_BINARY_DIR}/../examples/)

endif ( USE_SINKHORN )


####################################################################################################################


####################################################################################################################
## Spherical Reflector


if( USE_CPLEX )

	add_executable(ReflectorSpherical_CPLEX
		ReflectorSpherical_CPLEX.cpp)

	target_compile_options(ReflectorSpherical_CPLEX PUBLIC ${CPLEX_COMPILE_OPTIONS})
	target_include_directories(ReflectorSpherical_CPLEX PUBLIC ${CPLEX_INCLUDE_DIRECTORY})
	target_link_libraries(ReflectorSpherical_CPLEX
		LP_CPLEX
		ShortCutSolver
		${CPLEX_LIBRARY}
		Common
		)

	install (TARGETS ReflectorSpherical_CPLEX DESTINATION ${CMAKE_BINARY_DIR}/../examples/)

endif( USE_CPLEX )


if ( USE_SINKHORN )
	add_executable(ReflectorSpherical_Sinkhorn
		ReflectorSpherical_Sinkhorn.cpp)

	target_include_directories(ReflectorSpherical_Sinkhorn PUBLIC ${SINKHORN_INCLUDE_DIRECTORY})
	target_compile_options(ReflectorSpherical_Sinkhorn PUBLIC ${SINKHORN_COMPILE_OPTIONS})
	target_link_libraries(ReflectorSpherical_Sinkhorn
		Sinkhorn
		Common
		)

	install (TARGETS ReflectorSpherical_Sinkhorn DESTINATION ${CMAKE_BINARY_DIR}/../examples/)


endif( USE_SINKHORN )



