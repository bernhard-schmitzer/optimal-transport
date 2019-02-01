add_subdirectory (Models)

add_library(Common STATIC
	GridTools.cpp Tools.cpp
	TVarListHandler.cpp TCouplingHandler.cpp
	THierarchicalPartition.cpp THierarchyBuilder.cpp
	MultiScaleTools.cpp
	Verbose.cpp
	TEpsScaling.cpp THierarchicalCostFunctionProvider.cpp
	$<TARGET_OBJECTS:Models>
	)


install (TARGETS Common DESTINATION ${CMAKE_BINARY_DIR}/../bin)

install (FILES
	ErrorCodes.h
	GridTools.h
	MultiScaleTools.h
	PythonTypes.h
	TCostFunctionProvider.h
	TCouplingHandler.h
	TEpsScaling.h
	THierarchicalCostFunctionProvider.h
	THierarchicalPartition.h
	THierarchyBuilder.h
	Tools.h
	TVarListHandler.h
	Verbose.h
	DESTINATION ${CMAKE_BINARY_DIR}/../include/Common)
