add_library(ShortCutSolver STATIC
	Interfaces.cpp
	MultiScaleSolver.cpp
	TShortCutSolver.cpp
	TShieldGenerator.cpp
	TShieldGenerator_Models.cpp
	)

install (TARGETS ShortCutSolver DESTINATION ${CMAKE_BINARY_DIR}/../bin)

install (FILES
	Interfaces.h
	MultiScaleSolver.h
	TShieldGenerator.h
	TShieldGenerator_Models.h
	TShortCutSolver.h
	DESTINATION ${CMAKE_BINARY_DIR}/../include/ShortCutSolver)

