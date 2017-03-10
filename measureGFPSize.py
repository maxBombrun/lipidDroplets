# measureGFPSize.py
#######################################################
#
#	Creation of a new feature:
#	- Load images of the detected LDs and cells
#	- For each image, for each cell, 
#		-> measure each LDs
#		-> create a vector of LDs size distribution based on fixed binning
#	- Export the collection per image
#
#######################################################

import os

import cv2
import numpy as np

import matplotlib.pyplot as plt

from skimage.measure import label
from skimage.measure import regionprops

import multiprocessing
import Queue
import threading

import csv

import settings


## Thread
## Load images of cells and LDs, and  
## Measure the diameter of each LDs, cell-wise
## Create a vector of size distribution
def workerMeasure():
	while not q.empty():
		[idxPlate,idxWell, idxImg] = q.get()

		print ('Measuring P%s_%s'% (idxPlate, idxWell))

		imgGFP=cv2.imread(inputCellProfilerPath+'P'+idxPlate+'_'+idxWell+"_GFP.tif",0)
		imgMaskGFP=cv2.imread(outputDetPath +'P'+idxPlate+'_'+idxWell+"_GFP_MASK.tif",0)
		imgMaskCell=plt.imread(outputCellProfilerPath + 'P'+idxPlate+'_'+idxWell+"_GFP_Cell.tif", 0)  

		HistGFP=measureGFPSizePerCell(imgGFP,imgMaskGFP,imgMaskCell)

		currName= 'P'+idxPlate+'_'+idxWell
		# CellList 0 is the background !
		for cellList in range(1, len(HistGFP)):
			content= [idxPlate] + [idxImg] + [idxWell] + [cellList] + list(HistGFP[cellList])
			list2write.append(content)


## Function to detect lipid droplets in BODIPY channel ##
## Input: imgGFP: Original image of LDs
## imgMaskGFP: Mask of LDs
## labelCell: Mask of cells 
## Output: Collection of size distribution vectors for the current image
def measureGFPSizePerCell(imgGFP, imgMaskGFP, labelCell):
	## Create a label map of the masks
	labelGFP, numGFP=label(imgMaskGFP,return_num=True)
	numCell= np.amax(labelCell)

	listSizeGFPperCELL = [[]]*(numCell+1)

	## We place the size every GFP in the correct list based on the corresponding nucleus label
	for gfpObj in regionprops(labelGFP):
		xcurr,ycurr= gfpObj.centroid
		# if ycurr>debX and ycurr< debX+width and xcurr>debY and xcurr< debY+height:
		labCellOfCurrGFP= labelCell[int(xcurr),int(ycurr)] 
		listSizeGFPperCELL[labCellOfCurrGFP]=listSizeGFPperCELL[labCellOfCurrGFP]+[gfpObj.equivalent_diameter]        
	histoBins=[0, 5, 10, 15, 20, 25, 30, 35, 40]

	Hcell=[histPerso(x, histoBins) for x in listSizeGFPperCELL]

	return Hcell


## Function to create the vector of size distribution ##
## Input: inputList: size of LDs in the current cell
## listSize: Binning
## Output: res: number of LDs per bin
def histPerso(inputList,listSize):
	res= [0]*(len(listSize)-1)

	for i in inputList:
		if i>listSize[len(listSize)-1]:
			print ('%s is so big' % str(i))
			res[len(res)-1]=res[len(res)-1]+1
		else:
			j=0
			while i>= listSize[j]:
				j=j+1
			j=j-1
			res[j]=res[j]+1

	return res




#################### Main ####################

def measureGFP():
## Measurements !

	## Initialise variables
	global q
	global list2write
	global inputCellProfilerPath
	global outputDetPath
	global outputCellProfilerPath
	
	## Initialise paths		
	outputDetPath=settings.pathList[3]
	inputCellProfilerPath=settings.pathList[4]
	outputCellProfilerPath=settings.pathList[5]

	## Initialise threads
	nProc = multiprocessing.cpu_count()
	threadsList = [threading.Thread(target=workerMeasure) for i in range(0, nProc)]
	q = Queue.Queue()

	list2write=[]
	idxImg=0
 	listFiles= [x for x in os.listdir(outputCellProfilerPath) if x.endswith('Cell.tif')]

	## Fill up the queue for threading
	for currImg in listFiles:
		idxPlate=currImg.rsplit('_')
		idxWell=idxPlate[1]
		idxPlate=idxPlate[0][1:]

		idxImg=idxImg+1	
		q.put([idxPlate,idxWell,idxImg])

	# start all threads        
	for thread in threadsList:
		thread.start()
			
	# join all threads
	for thread in threadsList:
		thread.join()


	## Export size of the LDs cell-wise, non-fat included
	print('Writing csv...')

	list2write=sorted(list2write, key=lambda x: (x[0],x[1]))
	with open(outputDetPath+ "/GFPsize.csv",'wb') as csvfile:
		writer= csv.writer(csvfile, dialect='unixpwd')	
		header2write= ['Plate','ImageNumber', 'Well', 'ObjectNumber', 'bin_0_5', 'bin_5_10', 'bin_10_15', 'bin_15_20', 'bin_20_25', 'bin_25_30', 'bin_30_35', 'bin_35_40']
		writer.writerow(header2write)
		for content in list2write:
			writer.writerow(content)

	## Export size of the LDs cell-wise, only fat cells	 
	with open(outputDetPath+ "/GFPSize_Clean.csv",'wb') as mainFile:
		writer= csv.writer(mainFile, dialect='unixpwd')	
		writer.writerow(header2write)
		for currRow in list2write:
			if sum(map(int,currRow[4:]))!=0:
				writer.writerow(currRow)
