#!/bin/bash

. compile_constants.sh

echo "###########################################################"
echo "# Compiling Common Modules"
echo "###########################################################"


echo "# Common"

echo $FLAGS_INCLUDE $FLAGS_MISC -c -o Common/TCouplingHandler.o Common/TCouplingHandler.cpp 
$CC $FLAGS_INCLUDE $FLAGS_MISC -c -o Common/TCouplingHandler.o Common/TCouplingHandler.cpp

echo $FLAGS_INCLUDE $FLAGS_MISC -c -o Common/TVarListHandler.o Common/TVarListHandler.cpp
$CC $FLAGS_INCLUDE $FLAGS_MISC -c -o Common/TVarListHandler.o Common/TVarListHandler.cpp

echo $FLAGS_INCLUDE $FLAGS_MISC -c -o Common/TMultiVarListHandler.o Common/TMultiVarListHandler.cpp
$CC $FLAGS_INCLUDE $FLAGS_MISC -c -o Common/TMultiVarListHandler.o Common/TMultiVarListHandler.cpp

echo $FLAGS_INCLUDE $FLAGS_MISC -c -o Common/THierarchicalPartition.o Common/THierarchicalPartition.cpp
$CC $FLAGS_INCLUDE $FLAGS_MISC -c -o Common/THierarchicalPartition.o Common/THierarchicalPartition.cpp

echo $FLAGS_INCLUDE $FLAGS_MISC -c -o Common/THierarchicalPartitionInterface.o Common/THierarchicalPartitionInterface.cpp
$CC $FLAGS_INCLUDE $FLAGS_MISC -c -o Common/THierarchicalPartitionInterface.o Common/THierarchicalPartitionInterface.cpp

echo $FLAGS_INCLUDE $FLAGS_MISC -c -o Common/tools.o Common/tools.cpp
$CC $FLAGS_INCLUDE $FLAGS_MISC -c -o Common/tools.o Common/tools.cpp


echo "# CostFunctionComputation"


echo $FLAGS_INCLUDE $FLAGS_MISC -c -o CostFunctionComputation/Release/CostFunctionComputation.o \
		CostFunctionComputation/CostFunctionComputation.cpp
$CC $FLAGS_INCLUDE $FLAGS_MISC -c -o CostFunctionComputation/Release/CostFunctionComputation.o \
		CostFunctionComputation/CostFunctionComputation.cpp
	
echo -shared -o CostFunctionComputation/Release/libCostFunctionComputation.so \
	CostFunctionComputation/Release/CostFunctionComputation.o \
	Common/THierarchicalPartition.o \
	Common/THierarchicalPartitionInterface.o \
	Common/TVarListHandler.o \
	$FLAGS_LINKER

$CC -shared -o CostFunctionComputation/Release/libCostFunctionComputation.so \
	CostFunctionComputation/Release/CostFunctionComputation.o \
	Common/THierarchicalPartition.o \
	Common/THierarchicalPartitionInterface.o \
	Common/TVarListHandler.o \
	$FLAGS_LINKER

