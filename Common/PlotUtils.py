import numpy as np

def cm2inch(*tupl):
	inch = 2.54
	if isinstance(tupl[0], tuple):
		return tuple(i/inch for i in tupl[0])
	else:
		return tuple(i/inch for i in tupl)


def VecToHSV(field,maxval="max"):
	hsvField=np.zeros((field.shape[0],3))
	angle=np.arctan2(field[:,1],field[:,0])
	absval=np.linalg.norm(field,axis=1)
	hsvField[:,0]=np.maximum(0,np.minimum(1,(angle/(2*np.pi)+0.5)))
	hsvField[:,1]=1.
	if maxval=="max":
		maxval=np.max(absval)
	hsvField[:,2]=np.minimum(1,absval/maxval)
	
	return hsvField


