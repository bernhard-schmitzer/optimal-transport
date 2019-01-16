#!/bin/bash


. compile_constants.sh


echo "###########################################################"
echo "# Compiling CPLEX Basic"
echo "###########################################################"

echo $FLAGS_CPLEX $FLAGS_CPLEX_INCLUDE $FLAGS_INCLUDE $FLAGS_MISC -c -o OT_CPLEX/Release/OT_CPLEX.o \
	OT_CPLEX/OT_CPLEX.cpp
$CC $FLAGS_CPLEX $FLAGS_CPLEX_INCLUDE $FLAGS_INCLUDE $FLAGS_MISC -c -o OT_CPLEX/Release/OT_CPLEX.o \
	OT_CPLEX/OT_CPLEX.cpp


echo $FLAGS_CPLEX $FLAGS_CPLEX_INCLUDE $FLAGS_INCLUDE $FLAGS_MISC -c -o OT_CPLEX/Release/TCPLEXNetSolver.o \
	OT_CPLEX/TCPLEXNetSolver.cpp
$CC $FLAGS_CPLEX $FLAGS_CPLEX_INCLUDE $FLAGS_INCLUDE $FLAGS_MISC -c -o OT_CPLEX/Release/TCPLEXNetSolver.o \
	OT_CPLEX/TCPLEXNetSolver.cpp


echo -shared -o OT_CPLEX/Release/libOT_CPLEX.so \
	OT_CPLEX/Release/OT_CPLEX.o \
	OT_CPLEX/Release/TCPLEXNetSolver.o \
	Common/TCouplingHandler.o \
	Common/TVarListHandler.o \
	Common/tools.o \
	$FLAGS_CPLEX_LINKER $FLAGS_LINKER
$CC -shared -o OT_CPLEX/Release/libOT_CPLEX.so \
	OT_CPLEX/Release/OT_CPLEX.o \
	OT_CPLEX/Release/TCPLEXNetSolver.o \
	Common/TCouplingHandler.o \
	Common/TVarListHandler.o \
	Common/tools.o \
	$FLAGS_CPLEX_LINKER $FLAGS_LINKER


echo "###########################################################"
echo "# Compiling CPLEX ShortCuts"
echo "###########################################################"


echo $FLAGS_CPLEX $FLAGS_CPLEX_INCLUDE $FLAGS_SHORTCUT_INCLUDE $FLAGS_INCLUDE $FLAGS_MISC -c \
	-o ShortCutSolver_CPLEX/Release/Interfaces-CPLEX.o \
	ShortCutSolver_CPLEX/Interfaces-CPLEX.cpp
$CC $FLAGS_CPLEX $FLAGS_CPLEX_INCLUDE $FLAGS_SHORTCUT_INCLUDE $FLAGS_INCLUDE $FLAGS_MISC -c \
	-o ShortCutSolver_CPLEX/Release/Interfaces-CPLEX.o \
	ShortCutSolver_CPLEX/Interfaces-CPLEX.cpp


echo $FLAGS_CPLEX $FLAGS_CPLEX_INCLUDE $FLAGS_SHORTCUT_INCLUDE $FLAGS_INCLUDE $FLAGS_MISC -c \
	-o ShortCutSolver_CPLEX/Release/ShortCutSolver_CPLEX.o \
	ShortCutSolver_CPLEX/ShortCutSolver_CPLEX.cpp
$CC $FLAGS_CPLEX $FLAGS_CPLEX_INCLUDE $FLAGS_SHORTCUT_INCLUDE $FLAGS_INCLUDE $FLAGS_MISC -c \
	-o ShortCutSolver_CPLEX/Release/ShortCutSolver_CPLEX.o \
	ShortCutSolver_CPLEX/ShortCutSolver_CPLEX.cpp

echo -shared -o ShortCutSolver_CPLEX/Release/libShortCutSolver_CPLEX.so \
	ShortCutSolver_CPLEX/Release/Interfaces-CPLEX.o \
	ShortCutSolver_CPLEX/Release/ShortCutSolver_CPLEX.o \
	ShortCutSolver/Release/Interfaces.o \
	ShortCutSolver/Release/ShortCutSolver-Tools.o \
	OT_CPLEX/Release/TCPLEXNetSolver.o \
	Common/TCouplingHandler.o \
	Common/TVarListHandler.o \
	Common/tools.o \
	$FLAGS_CPLEX_LINKER $FLAGS_LINKER
$CC -shared -o ShortCutSolver_CPLEX/Release/libShortCutSolver_CPLEX.so \
	ShortCutSolver_CPLEX/Release/Interfaces-CPLEX.o \
	ShortCutSolver_CPLEX/Release/ShortCutSolver_CPLEX.o \
	ShortCutSolver/Release/Interfaces.o \
	ShortCutSolver/Release/ShortCutSolver-Tools.o \
	OT_CPLEX/Release/TCPLEXNetSolver.o \
	Common/TCouplingHandler.o \
	Common/TVarListHandler.o \
	Common/tools.o \
	$FLAGS_CPLEX_LINKER $FLAGS_LINKER
