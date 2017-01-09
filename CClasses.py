import ctypes as ct
#import numpy as np

############################################################################################################################################
# Auxiliary Typecasting Functions
############################################################################################################################################

c_double_p = ct.POINTER(ct.c_double)
c_int_p = ct.POINTER(ct.c_int)
c_long_p = ct.POINTER(ct.c_long)

class TDoubleMatrix(ct.Structure):
	_fields_ = [("data",c_double_p),("depth",ct.c_int),("dimensions",c_int_p)]

class TInteger64Matrix(ct.Structure):
	_fields_ = [("data",c_long_p),("depth",ct.c_int),("dimensions",c_int_p)]


class TInteger32Matrix(ct.Structure):
	_fields_ = [("data",c_int_p),("depth",ct.c_int),("dimensions",c_int_p)]

    
def _toCInt(a):
	return ct.c_int(a)

def _toCDouble(a):
	return ct.c_double(a)

###

def _toCDoubleMatrix(a):
	return ct.byref(_toCDoubleMatrix_pre(a) )

def _toCInteger64Matrix(a):
	return ct.byref(_toCInteger64Matrix_pre(a) )

def _toCInteger32Matrix(a):
	return ct.byref(_toCInteger32Matrix_pre(a) )

###

def _toCDoubleMatrix_pre(a):
	return TDoubleMatrix(a.ctypes.data_as(c_double_p),_toCInt(a.ndim),a.ctypes.shape_as(ct.c_int))

def _toCInteger64Matrix_pre(a):
	return TInteger64Matrix(a.ctypes.data_as(c_long_p),_toCInt(a.ndim),a.ctypes.shape_as(ct.c_int))

def _toCInteger32Matrix_pre(a):
	return TInteger32Matrix(a.ctypes.data_as(c_int_p),_toCInt(a.ndim),a.ctypes.shape_as(ct.c_int))

