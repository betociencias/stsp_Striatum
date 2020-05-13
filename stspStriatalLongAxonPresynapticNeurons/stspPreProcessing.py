import matplotlib.pylab as gr
import scipy as sc
from scipy import fftpack
#sc.test("all")
import os
import csv
import glob

#Exploring the synaptic plasticity between Cortico-Striatal and Thalamo-striatal synapses
#srcDir="$HOME/data/synapticPlasticityMario/STP_Cortico-Str_Talamo-Str/"

def getFileList(dataDir, prefix, suffix, includeDataDir=1):
	"""
	Example:
	files=getFileList(dataDir, prefix, suffix)
	"""
	files=list()
	for f in os.listdir(dataDir):
		a= sc.int32(os.path.isfile(os.path.join(dataDir,f)))
		b= sc.int32(str.find(f,prefix)>-1)
		c= sc.int32(str.find(f,suffix)>0)
		#print(a,b,c)
		if (a*b*c):
			if includeDataDir:
				files.append(dataDir+f)
			else:
				files.append(f)
	nFiles = len(files)
	print("Found %d files with the indicated string"%nFiles)
	print(files)
	return files

def getDirNames(myPath):
	dirNames = [d for d in os.listdir(myPath) if os.path.isdir(os.path.join(myPath, d))]
	return dirNames

def getFileNames(myPath):
	fNames = [f for f in os.listdir(myPath) if os.path.isfile(os.path.join(myPath, f))]
	return fNames

def walkDir(srcDir):
	f = []
	for (dirpath, dirnames, filenames) in os.walk(srcDir):
		f.extend(filenames)
	return f

def readCSV(csvF,nHeaderRows=5, delimiter=" "):
	data=list()
	header=list()
	n=0
	f=open(csvF, 'r')
	reader = csv.reader(f, delimiter=delimiter)
	n=0
	for row in reader:
		if n<nHeaderRows:
			#print(row)
			header.append(row)
			n=n+1
		else:			
			break

	for row in reader:
		#print(row)
		data.append(sc.float32(row))
	f.close()
	return sc.array(data),header


def cleanLHDerivative(epsc,trial, upThresh=1e5,dnThresh=-1e5,sampInterval=1/6000.0, offsetPoints=10,showJob=0, showDeriv=0):
	"""
	Example:
	cc, stimStarts=cleanLHDerivative(epsc,upThresh,dnThresh,sampInterval,showJob,showDeriv)
	Arguments:
	upThresh, dnThresh ~ parameters are upper and lower bounds for the derivative (to keep)
	sampInterval ~ the time between sampling points
	offsetPoints ~ # number of points to take into account after the end of the high-derivative intervals.
	"""
	dcdt=sc.zeros(len(epsc))
	dcdt[1:]=(epsc[1:]-epsc[:-1])/sampInterval
	badInds= sc.where( (dcdt<dnThresh) | (dcdt>upThresh))[0]
	#print(badInds)
	# Make sure the pieces are contiguous
	aa=sc.where(badInds[1:]-badInds[:-1]>10)[0]
	ab=sc.hstack((badInds[0],badInds[aa+1]))
	bb=sc.hstack((badInds[aa],badInds[-1]))+offsetPoints
	# print(aa,"\n",ab,"\n",bb)
	# 
	cleanEPSC=sc.copy(epsc)
	if ((len(bb)==len(ab)) & (len(bb)>0)):
		for nn in sc.arange(len(bb)):
			cleanEPSC[ab[nn]:bb[nn]]=cleanEPSC[ab[nn]-1]
	if showJob:
		f0=gr.figure(figsize=(15,5))
		r=1; c=1; gr.ioff()
		ax0=f0.add_subplot(r,c,1)
		ax0.plot(epsc,"b",lw=1, alpha=0.9, label=r"$I_{%d}$"%(trial))
		ax0.plot(cleanEPSC,"k",lw=3, alpha=0.4, label=r"$I_{%d}$"%(trial))
		ax0.set_xlabel("indices")
	if showDeriv:
		ax0.plot(badInds,dcdt[badInds]/1000.0,"wo",ms=5,label=r"$\partial_t I_{%d}(bad)$"%(trial))
		ax0.plot(dcdt/1000.0,"k.",ms=3,label=r"$\partial_t I_{%d}$"%(trial))
	if showJob+showDeriv:
		ax0.legend()
		gr.ion(); gr.draw()
	return cleanEPSC, ab, bb

def cleanListLHDerivative(trials,upThresh=1e5,dnThresh=-1e5,sampInterval=1/6000.0,offsetPoints=10,showJob=0,showDeriv=0):
	nTrials=len(trials)
	cleanRecs= list(); artStarts=list(); artEnds=list()
	for n in sc.arange(nTrials):
		epsc=trials[n]
		cc,ab,bb=cleanLHDerivative(epsc,upThresh,dnThresh,sampInterval,offsetPoints,showJob,showDeriv)
		cleanRecs.append(cc)
		artStarts.append(ab)
		artEnds.append(bb)
	return cleanRecs,artStarts,artEnds

# FFT low-high cut
def lowHighCutFreqSmoothing(x,lowCut=400,hiCut=600):
	rft= fftpack.rfft(x); 
	rft[:lowCut]=0; 
	rft[hiCut:]=0; 
	#y_smooth = sc.ifft(rft, N)
	y=fftpack.irfft(rft); 
	return y #,y_smooth

def getSamplesFromList(fNames,srcDir="./",sampRate=6000.0, 
	upThresh=1e5,dnThresh=-1e5, lowCut=0, hiCut=1000, offsetPoints=10,
	cleanArtifacts=1, smooth=1):
	nFiles = len(fNames)
	data=list()
	for n in sc.arange(nFiles): 
		data.append(dict())
		fName=srcDir+fNames[n]
		#print(fName); 
		dd,header=readCSV(fName,nHeaderRows=5,delimiter=",")
		data[n]["file"]=fNames[n][:-4].split()[-1]
		data[n]["srcDir"]=srcDir
		data[n]["trains"]=dd.transpose()
		data[n]["header"]=header
		data[n]["nTrains"]= len(data[n]["trains"])
		#sampRate=1000*sc.float32(header[0][0].split()[-1])
		data[n]["Hz"]=sampRate
		data[n]["sampInterval"]=1/sampRate
		data[n]["sampTimes"]=sc.arange(0,len(dd.transpose()[0])/sampRate,1/sampRate)
		if cleanArtifacts:
			cc, artStarts,artEnds=cleanListLHDerivative(data[n]["trains"],showJob=0,showDeriv=0)
			data[n]["cleanTrains"]= cc
			data[n]["stimStartInds"]=artStarts
			data[n]["stimEndInds"]=artEnds		
		if smooth:
			sRecs=list()
			for m in sc.arange(data[n]["nTrains"]):
				sRecs.append(lowHighCutFreqSmoothing(cc[m],lowCut,hiCut))
			data[n]["smoothTrains"]=sc.array(sRecs)
	return data

def smoothTrainList(trainList, lowCut=0, hiCut=1000):
	nTrains=len(trainList)
	sRecs=list()
	for m in sc.arange(nTrains):
		sRecs.append(lowHighCutFreqSmoothing(trainList[m],lowCut,hiCut))
	return sRecs

def extractRecording(prefix, suffix,srcDir="./.",sampRate=6000.0, 
	upThresh=1e5,dnThresh=-1e5,offsetPoints=10, lowCut=0, hiCut=1000,
	cleanArtifacts=1, smooth=1):
    fNames= getFileList(srcDir, prefix, suffix,includeDataDir=0)
    data= getSamplesFromList(fNames,srcDir,sampRate, upThresh,dnThresh,
    	lowCut,hiCut,offsetPoints,cleanArtifacts,smooth)
    return data

def splitTrain(train,spliceInds, nPtsPerPiece=10):
	"""
	Syntax:
	pieces=splitTrain(train,spliceInds, nPtsPerPiece=1)
	"""
	nPieces = len(spliceInds)
	pieces=list()
	ii= sc.minimum(spliceInds[0]+1,nPieces)
	baseLine=train[ii]
	if nPtsPerPiece>1:
		for n in sc.arange(nPieces):
			a = spliceInds[n]+1
			b = a+nPtsPerPiece
			print(a,b)
			pieces.append(train[a:b]-baseLine)
	else:
		print("The number of points per piece is too small\n")
	return sc.array(pieces)

def splitTimeSeries(x,y,spliceTimes, nPtsPerPiece=1):
	"""
	Syntax:
	xPieces,yPieces=splitTimeSeries(x,y,spliceTimes, nPtsPerPiece=1)
	"""
	nPieces = len(spliceTimes)
	xPieces= list(); yPieces=list()
	for n in sc.arange(nPieces-1):
		a = sc.int32(sc.where(x>=spliceTimes[n])[0].min())
		b = sc.int32(a + nPtsPerPiece)
		xPieces.append(x[a:b])
		yPieces.append(y[a:b])
	return sc.array(xPieces), sc.array(yPieces)

def calcMultipleLocalPeaks(timeSeriesList):
	"""
	calcMultipleLocalExtrema(timeSeriesList)
	"""
	M=list(); 
	nItems = len(timeSeriesList)
	for n in sc.arange(nItems):
		M.append(timeSeriesList[n].max())
	return M


def smoothSplitTrainList(trainList, sampTimes, startTimes, nPts = 300.0, low=0, hi=1000, dt =1.0/6000.0,calcMin=1):
	"""
	smoothSplitTrainList(trainList, sampTimes, startTimes, nPts = 300.0, low=0, hi=1000, dt =1.0/recV["Hz"],)
	takes list of responses to stimulation trains, 
	smooths the responses, and splits them with respect to the stimulus times. 
	The function returns a list with as many elements as stimuli, 
	each element is an array containing the responses to each stimulus. 
	"""
	x_smooth=lowFreqSmoothing(trainList,lowCut=low,hiCut=hi)
	nTrains = len(trainList)
	print("Splitting %d trains"%nTrains)
	splitTrains=list()
	smoothSplitTrains=list()
	for n in sc.arange(nTrains):
		print("Splitting train %d"%n)
		x = trainList[n]
		x_smooth=lowFreqSmoothing(x,lowCut=low,hiCut=hi)
		pieces,m=splitTrace(sampTimes, x,startTimes, nPts)
		sm_pieces,sm=splitTrace(sampTimes, x_smooth, startTimes, nPts,calcMin)
		splitTrains.append(pieces)
		smoothSplitTrains.append(sm_pieces)
	return splitTrains,m, smoothSplitTrains
