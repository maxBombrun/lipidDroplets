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



settings.init() 

CPPath=settings.pathList[0]
inputDataPath=settings.pathList[1]
resultPath=settings.pathList[2]
outputDetPath=settings.pathList[3]
inputCellProfilerPath=settings.pathList[4]
outputCellProfilerPath=settings.pathList[5]



nProc = multiprocessing.cpu_count()
listPlates= [x for x in os.listdir(inputDataPath) if os.path.isdir(inputDataPath+x) and x.startswith('plate')]

csv.register_dialect('unixpwd', delimiter=',',quotechar = '"', doublequote = True, skipinitialspace = False,lineterminator = '\n',  quoting = csv.QUOTE_NONE)

print listPlates

## Segmentation of the nuclei and the lipid droplets
segmentNucAndGFP.segmentFatDroplet(listPlates)

## Cells approximation based on the previous segmentation
## and features extraction through CellProfiler
cellProfilerGetRelation.runCellProfilerRelationship()

## Individual lipid droplet measurements
## Creation of size distribution vectors
measureGFPSize.measureGFP()

## Classification of the cells based on the vectors
clusterDroplets.getClusterOfDroplets(nbClass=2)

## Creation of CSV output, summarizing the measurements per-cell and per-well
fusionCSV.getPerImageMeasurements()

## Plotting of the features
plotFeatures.plotFeat()

## Validation of the features through the computation of the Zprime factor
computeZprime.getZprime()


