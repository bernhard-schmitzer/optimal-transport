#!/bin/bash

CC="gcc"

FLAGS_INCLUDE="-ICommon/"
FLAGS_MISC="-fPIC -O3 -Wall -Wextra"
FLAGS_LINKER="-lm -lstdc++"
FLAGS_SHORTCUT_INCLUDE="-IShortCutSolver/"



#############################################
# CPLEX
#############################################

FLAGS_CPLEX="-DNDEBUG -DIL_STD -DILOSTRICTPOD"

# Please replace XXXX by the actual path to the libraries on your system.
# Be careful that also the version number and the architecture type (x86-64_linux or x86-64_osx) is correct.
# In the following there are some examples. Make sure to remove the comments in front of the setting that you want to use!


# Linux

# v12.5
#FLAGS_CPLEX_INCLUDE="-I/XXXX/cplex12.5/cplex/include -IOT_CPLEX/"
#FLAGS_CPLEX_LINKER="-L/XXXX/cplex12.5/cplex/lib/x86-64_sles10_4.1/static_pic/ -lcplex"

# v12.61
#FLAGS_CPLEX_INCLUDE="-I/XXXX/cplex12.6.1/cplex/include -IOT_CPLEX/"
#FLAGS_CPLEX_LINKER="-L/XXXX/cplex12.6.1/cplex/lib/x86-64_linux/static_pic/ -lcplex"

# OS X

# v12.61
#FLAGS_CPLEX_INCLUDE="-I/XXXX/cplex12.6.1/cplex/include -IOT_CPLEX/"
#FLAGS_CPLEX_LINKER="-L/XXXX/cplex12.6.1/cplex/lib/x86-64_osx/static_pic/ -lcplex"

#############################################
# LEMON
#############################################

# Please replace XXXX by the path to the LEMON library build on your system.
# Relative to XXXX/ there should be (among others) the files lemon/list_graph.h and lemon/libemon.a.
# Make sure to remove the comments after adjusting the location!

#FLAGS_LEMON_INCLUDE="-I/XXXX/ -IOT_Lemon/"
#FLAGS_LEMON_LINKER="/XXXX/lemon/libemon.a"

