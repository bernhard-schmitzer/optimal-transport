def addParams(base,extraNames,extraVals):
	for n,v in zip(extraNames,extraVals):
		base[n]=v

def loopAtLevel(base,extraNames,extraVals,level,posList,callback=None,mask=None):
	# check via mask, if current level is to be looped

	doLoop=True
	if mask is not None:
		if mask[level]==0:
			doLoop=False

	if doLoop:
		# if looping is required, go over full range in level
		loopRange=range(len(extraVals[level]))
	else:
		# otherwise just go over value 0
		loopRange=[0]

	# loop through all entries at current level
	for i in loopRange:
		#print("level: "+str(level) + "\tvalue:" + str(i))
		posList[level]=i
		addParams(base,extraNames[level],extraVals[level][i])
		# while not at finest level
		if level<len(extraNames)-1:
			# set entries at finer level to -1
			for j in range(level+1,len(extraNames)):
				posList[j]=-1
			# call routine at next finer level
			loopAtLevel(base,extraNames,extraVals,level+1,posList,callback,mask)
		# at finest level
		else:
			if callback is None:
				print(posList)
				print(base)
			else:
				callback(posList,base)

def loopParams(base,extraNames,extraVals,callback=None,mask=None):
	tmpBase=base.copy()
	posList=[0 for i in range(len(extraNames))]
	loopAtLevel(base,extraNames,extraVals,0,posList,callback=callback,mask=mask)

def getTagFileName(base_name,posList,extraFilenames):
	result=base_name
	for i,p in enumerate(posList):
		result=result+"_"+extraFilenames[i][p]
	return result
