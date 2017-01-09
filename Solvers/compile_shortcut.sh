#!/bin/bash

. compile_constants.sh

echo "###########################################################"
echo "# Compiling ShortCutSolver"
echo "###########################################################"

echo $FLAGS_INCLUDE $FLAGS_MISC -c -o ShortCutSolver/Release/ShieldGenerator/TShieldGenerator.o \
	ShortCutSolver/ShieldGenerator/TShieldGenerator.cpp
$CC $FLAGS_INCLUDE $FLAGS_MISC -c -o ShortCutSolver/Release/ShieldGenerator/TShieldGenerator.o \
	ShortCutSolver/ShieldGenerator/TShieldGenerator.cpp

echo $FLAGS_INCLUDE $FLAGS_MISC -c -o ShortCutSolver/Release/Interfaces.o ShortCutSolver/Interfaces.cpp
$CC $FLAGS_INCLUDE $FLAGS_MISC -c -o ShortCutSolver/Release/Interfaces.o ShortCutSolver/Interfaces.cpp

echo $FLAGS_INCLUDE $FLAGS_MISC -c -o ShortCutSolver/Release/ShortCutSolver-Interface.o ShortCutSolver/ShortCutSolver-Interface.cpp
$CC $FLAGS_INCLUDE $FLAGS_MISC -c -o ShortCutSolver/Release/ShortCutSolver-Interface.o ShortCutSolver/ShortCutSolver-Interface.cpp

echo $FLAGS_INCLUDE $FLAGS_MISC -c -o ShortCutSolver/Release/ShortCutSolver-Tools.o ShortCutSolver/ShortCutSolver-Tools.cpp
$CC $FLAGS_INCLUDE $FLAGS_MISC -c -o ShortCutSolver/Release/ShortCutSolver-Tools.o ShortCutSolver/ShortCutSolver-Tools.cpp

echo $FLAGS_INCLUDE $FLAGS_MISC -c -o ShortCutSolver/Release/TShortCutSolver.o ShortCutSolver/TShortCutSolver.cpp
$CC $FLAGS_INCLUDE $FLAGS_MISC -c -o ShortCutSolver/Release/TShortCutSolver.o ShortCutSolver/TShortCutSolver.cpp



echo -shared -o ShortCutSolver/Release/libShortCutSolver.so \
	ShortCutSolver/Release/ShieldGenerator/TShieldGenerator.o \
	ShortCutSolver/Release/Interfaces.o \
	ShortCutSolver/Release/ShortCutSolver-Interface.o \
	ShortCutSolver/Release/ShortCutSolver-Tools.o \
	ShortCutSolver/Release/TShortCutSolver.o \
	Common/THierarchicalPartition.o \
	Common/TCouplingHandler.o \
	Common/TVarListHandler.o \
	Common/tools.o \
	$FLAGS_LINKER


$CC -shared -o ShortCutSolver/Release/libShortCutSolver.so \
	ShortCutSolver/Release/ShieldGenerator/TShieldGenerator.o \
	ShortCutSolver/Release/Interfaces.o \
	ShortCutSolver/Release/ShortCutSolver-Interface.o \
	ShortCutSolver/Release/ShortCutSolver-Tools.o \
	ShortCutSolver/Release/TShortCutSolver.o \
	Common/THierarchicalPartition.o \
	Common/TCouplingHandler.o \
	Common/TVarListHandler.o \
	Common/tools.o \
	$FLAGS_LINKER


