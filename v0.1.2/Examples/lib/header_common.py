import sys
import os
import numpy as np

import scipy
import scipy.io as sciio

import datetime

from .header_rootdir import *

import ScriptTools
from ScriptTools import ParameterType as ptype

import OTTools
import OTTools.HierarchicalPartition as HierarchicalPartition
import Solvers.CostFunctionComputation as SolverCFC

