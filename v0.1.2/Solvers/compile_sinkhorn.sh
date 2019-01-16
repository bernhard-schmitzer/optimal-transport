#!/bin/bash

CC="gcc"

FLAGS_INCLUDE="-ICommon/"
FLAGS_MISC="-fPIC -O3 -Wall -Wextra"
FLAGS_LINKER="-lm -lstdc++"


echo "###########################################################"
echo "# Compiling Sinkhorn Solver"
echo "###########################################################"


echo $FLAGS_INCLUDE $FLAGS_MISC -c -o Sinkhorn/Release/TSinkhorn.o Sinkhorn/TSinkhorn.cpp
$CC $FLAGS_INCLUDE $FLAGS_MISC -c -o Sinkhorn/Release/TSinkhorn.o Sinkhorn/TSinkhorn.cpp

echo $FLAGS_INCLUDE $FLAGS_MISC -c -o Sinkhorn/Release/Sinkhorn.o Sinkhorn/Sinkhorn.cpp
$CC $FLAGS_INCLUDE $FLAGS_MISC -c -o Sinkhorn/Release/Sinkhorn.o Sinkhorn/Sinkhorn.cpp



echo -shared -o Sinkhorn/Release/libSinkhorn.so \
	Sinkhorn/Release/Sinkhorn.o \
	Sinkhorn/Release/TSinkhorn.o \
	Common/THierarchicalPartition.o \
	Common/TMultiVarListHandler.o \
	Common/TVarListHandler.o \
	$FLAGS_LINKER


$CC -shared -o Sinkhorn/Release/libSinkhorn.so \
	Sinkhorn/Release/Sinkhorn.o \
	Sinkhorn/Release/TSinkhorn.o \
	Common/THierarchicalPartition.o \
	Common/TMultiVarListHandler.o \
	Common/TVarListHandler.o \
	$FLAGS_LINKER


