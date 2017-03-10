import os
import csv
import multiprocessing

import settings
import segmentNucAndGFP
import cellProfilerGetRelation
import measureGFPSize
import plotFeatures
import fusionCSV
import clusterDroplets
import computeZprime




	


############################# MAIN #####################################3 


settings.init() 

CPPath=settings.pathList[0]
inputDataPath=settings.pathList[1]
resultPath=settings.pathList[2]
outputDetPath=settings.pathList[3]
inputCellProfilerPath=settings.pathList[4]
outputCellProfilerPath=settings.pathList[5]



listPlates= [x for x in os.listdir(inputDataPath) if os.path.isdir(inputDataPath+x) and x.startswith('plate')]


csv.register_dialect('unixpwd', delimiter=',',quotechar = '"', doublequote = True, skipinitialspace = False,lineterminator = '\n',  quoting = csv.QUOTE_NONE)


## Auto
nProc = multiprocessing.cpu_count()
# listPlates=listPlates[21:]

## Manual
# nProc=7
# listPlates=['plate1', 'plate2', 'plate3', 'plate4', 'plate5', 'plate6', 'plate7', 'plate8', 'plate9', 'plate10', 'plate11']



print listPlates


segFlag=False
CPFlag=False
measFlag=False

clusterFlag=False
fusionCSVFlag=False

DBFlag=False
featureFlag= True
zprimeFlag=True


if segFlag:
	segmentNucAndGFP.segmentFatDroplet(listPlates)		

if CPFlag:
	cellProfilerGetRelation.runCellProfilerRelationship()

if measFlag:
	measureGFPSize.measureGFP()

if clusterFlag:
	clusterDroplets.getClusterOfDroplets(nbClass=2)

if fusionCSVFlag:
	fusionCSV.getPerImageMeasurements()

if featureFlag:
	plotFeatures.plotFeat()

if zprimeFlag:
	computeZprime.getZprime()



print "Done"


############ END MAIN #################




