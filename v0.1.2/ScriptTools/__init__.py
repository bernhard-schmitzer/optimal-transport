#import glob
import numpy as np
#import scipy.io as sciio
from enum import Enum
import sys
import csv
csv.field_size_limit(1000000)
#import os
#import re


##########################################################################################################################################
# Script Parameter Handling

class ParameterType(Enum):
	"""Enum type for parameter type specification."""
	integer=1
	real=2
	string=3
	integerList=4
	realList=5
	stringList=6


def parseParameter(v,t):
	"""The parameter value v (given as str) is parsed according to type specification t (ParameterType)."""
	if t==ParameterType.string:
		return v
	elif t==ParameterType.integer:
		return np.fromstring(v,dtype=np.int,count=1,sep=" ")[0]
	elif t==ParameterType.real:
		return np.fromstring(v,dtype=np.double,count=1,sep=" ")[0]
	elif t==ParameterType.integerList:
		return np.fromstring(v,dtype=np.int,sep=",")
	elif t==ParameterType.realList:
		return np.fromstring(v,dtype=np.double,sep=",")
	elif t==ParameterType.stringList:
		return v.split(",")

def readParameters(filename,paramsType):
	"""File <filename> is opened and read as tab separated CSV. The dictionary paramsType specifies which parameters are read: each entry
		specifies one parameter. Name given by the key, type by the value of type ParameterType. If a key is found in the file, its
		value is read, parsed and added to the result dictionary under the given key."""
	stream=open(filename,"r")
	params_reader=csv.reader(stream,delimiter="\t",strict=True)
	params={}
	for dat in params_reader:
		if len(dat)==2:
			key,val=dat
			if len(key)>0:
				if key[0]!="#":
					if key in paramsType.keys():
						params[key]=parseParameter(val,paramsType[key])
	return params


def getCommandLineParameters(params,newVariables,offset=1,src=sys.argv):
	"""List of command line arguments is parsed according to newVariables: each variable gives name and type. Variables are added to params."""
	if len(src)-offset<len(newVariables):
		raise ValueError("Not enough command line parameters given.")
	for n,t in newVariables:
		v=src[offset]
		params[n]=parseParameter(v,t)
		offset+=1


def saveParameters(filename,params):
	"""Save dictionary to CSV file."""
	w=csv.writer(open(filename,"w"),delimiter="\t")
	for key in sorted(params):
		w.writerow([key, params[key]])

def execFromFile(filename,glob_vars=None,loc_vars=None):
	"""Load code from a file and run it in given scope."""
	f=open(filename)
	content=f.read()
	f.close()
	code = compile(content, filename, 'exec')
	exec(code, glob_vars, loc_vars)
