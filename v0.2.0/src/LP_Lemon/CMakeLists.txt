if( USE_LEMON )
	add_library(LP_Lemon STATIC TLemonSolver.cpp Interfaces-Lemon.cpp)
	target_include_directories(LP_Lemon PUBLIC ${LEMON_INCLUDE_DIRECTORY})

	install (TARGETS LP_Lemon DESTINATION ${CMAKE_BINARY_DIR}/../bin)
	
	
	install (FILES
		Interfaces-Lemon.h
		TLemonSolver.h
		DESTINATION ${CMAKE_BINARY_DIR}/../include/LP_Lemon)

endif( USE_LEMON )



