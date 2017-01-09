#!/bin/bash


. compile_constants.sh


echo "###########################################################"
echo "# Compiling Lemon Basic"
echo "###########################################################"

echo $FLAGS_LEMON_INCLUDE $FLAGS_INCLUDE $FLAGS_MISC -c -o OT_Lemon/Release/OT_Lemon.o \
	OT_Lemon/OT_Lemon.cpp
$CC $FLAGS_LEMON_INCLUDE $FLAGS_INCLUDE $FLAGS_MISC -c -o OT_Lemon/Release/OT_Lemon.o \
	OT_Lemon/OT_Lemon.cpp

echo $FLAGS_LEMON_INCLUDE $FLAGS_INCLUDE $FLAGS_MISC -c -o OT_Lemon/Release/TLemonSolver.o \
	OT_Lemon/TLemonSolver.cpp
$CC $FLAGS_LEMON_INCLUDE $FLAGS_INCLUDE $FLAGS_MISC -c -o OT_Lemon/Release/TLemonSolver.o \
	OT_Lemon/TLemonSolver.cpp

echo -shared -o OT_Lemon/Release/libOT_Lemon.so \
	OT_Lemon/Release/OT_Lemon.o \
	OT_Lemon/Release/TLemonSolver.o \
	Common/TCouplingHandler.o \
	Common/TVarListHandler.o \
	$FLAGS_LEMON_LINKER $FLAGS_LINKER
$CC -shared -o OT_Lemon/Release/libOT_Lemon.so \
	OT_Lemon/Release/OT_Lemon.o \
	OT_Lemon/Release/TLemonSolver.o \
	Common/TCouplingHandler.o \
	Common/TVarListHandler.o \
	$FLAGS_LEMON_LINKER $FLAGS_LINKER

echo "###########################################################"
echo "# Compiling Lemon ShortCuts"
echo "###########################################################"


echo $FLAGS_LEMON_INCLUDE $FLAGS_INCLUDE $FLAGS_SHORTCUT_INCLUDE $FLAGS_MISC -c \
	-o ShortCutSolver_Lemon/Release/Interfaces-Lemon.o \
	ShortCutSolver_Lemon/Interfaces-Lemon.cpp

$CC $FLAGS_LEMON_INCLUDE $FLAGS_INCLUDE $FLAGS_SHORTCUT_INCLUDE $FLAGS_MISC -c \
	-o ShortCutSolver_Lemon/Release/Interfaces-Lemon.o \
	ShortCutSolver_Lemon/Interfaces-Lemon.cpp

echo $FLAGS_LEMON_INCLUDE $FLAGS_INCLUDE $FLAGS_SHORTCUT_INCLUDE $FLAGS_MISC -c \
	-o ShortCutSolver_Lemon/Release/ShortCutSolver_Lemon.o \
	ShortCutSolver_Lemon/ShortCutSolver_Lemon.cpp

$CC $FLAGS_LEMON_INCLUDE $FLAGS_INCLUDE $FLAGS_SHORTCUT_INCLUDE $FLAGS_MISC -c \
	-o ShortCutSolver_Lemon/Release/ShortCutSolver_Lemon.o \
	ShortCutSolver_Lemon/ShortCutSolver_Lemon.cpp

echo -shared -o ShortCutSolver_Lemon/Release/libShortCutSolver_Lemon.so \
	ShortCutSolver_Lemon/Release/Interfaces-Lemon.o \
	ShortCutSolver_Lemon/Release/ShortCutSolver_Lemon.o \
	OT_Lemon/Release/TLemonSolver.o \
	ShortCutSolver/Release/Interfaces.o \
	ShortCutSolver/Release/ShortCutSolver-Tools.o \
	Common/TCouplingHandler.o \
	Common/TVarListHandler.o \
	Common/tools.o \
	$FLAGS_LEMON_LINKER $FLAGS_LINKER

$CC -shared -o ShortCutSolver_Lemon/Release/libShortCutSolver_Lemon.so \
	ShortCutSolver_Lemon/Release/Interfaces-Lemon.o \
	ShortCutSolver_Lemon/Release/ShortCutSolver_Lemon.o \
	OT_Lemon/Release/TLemonSolver.o \
	ShortCutSolver/Release/Interfaces.o \
	ShortCutSolver/Release/ShortCutSolver-Tools.o \
	Common/TCouplingHandler.o \
	Common/TVarListHandler.o \
	Common/tools.o \
	$FLAGS_LEMON_LINKER $FLAGS_LINKER
