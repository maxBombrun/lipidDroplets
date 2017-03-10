# fusionCSV.py
#######################################################
#
#	Fusion the CSV from the parallelisation:
#	- Add the different measurement and the LDs size distribution together
#	- Export in a per-cell csv and per-well csv files
#
#######################################################

import os

import math
import csv

import settings

from scipy.stats import gaussian_kde
import numpy as np
import matplotlib.pyplot as plt

import copy

## Proof for the manual intensity threshold
def plotDensityDistribution(list2plot, name2give):
	colorTab=['r','g','b']
	plt.figure(figsize=(10,10))
	for idxDim in range(len(list2plot)):
		currList=list2plot[idxDim]
		currList= [item for sublist in currList for item in sublist]
		density = gaussian_kde(currList)
		if name2give=='Integrated':
			xs = np.linspace(-0.5 ,30,200)
		else:
			xs = np.linspace(-0.05 ,0.05,200)
		density.covariance_factor = lambda : .25
		density._compute_covariance()

		plt.plot(xs,density(xs), colorTab[idxDim])
	plt.savefig(outputDetPath+ "/gfpIntDis/"+str(name2give)+".png", bbox_inches='tight')
	plt.close()	



#################### Main ####################

def getPerImageMeasurements():
## Fusion !

	## Initialise global variables
	global q
	global outputDetPath

	## Initialise paths		
	outputDetPath=settings.pathList[3]
	outputCellProfilerPath=settings.pathList[5]

	## Initialise local variables
	idxCell=0
	idxImg=0
	idxImgMax=0
	idxColumnBef=[0, 1, 4, 8, 10]
	idxColumnAfter=[0, 1, 6, 6]	
	listOfNum=[]
	listOfIntBef=[]
	listOfIntAfter=[]
	listOfCell=[]
	listForCell=[]


	print('Fusion in progress...')

	## Load the LDs size distribution
	with open(outputDetPath+ "/GFPSize.csv",'r') as csvfile:
		reader= csv.reader(csvfile)	
		content2Store = list(reader)

	with open(outputDetPath+ "/GFPSize_UNCLASSED.csv",'r') as csvfile:
		reader= csv.reader(csvfile)	
		listOfClass = list(reader)

	## Load the measurements, output of the CellProfiler pipeline
	listAnalysis= [x for x in os.listdir(outputCellProfilerPath+ "temp/") if os.path.isdir(outputCellProfilerPath+ "temp/"+x) and x.startswith('Analysis_')]
	
	for dbIdx in range(len(listAnalysis)):
		currSubPath=outputCellProfilerPath+ "temp/Analysis_"+str(dbIdx)+"/"
		if os.path.isdir(currSubPath):	
			with open(currSubPath+ "/MyExpt_Cells.csv",'r') as csvfile:
				reader= csv.reader(csvfile)	
				listOfCell2read = list(reader)	

		for idxLine in range(len(listOfCell2read[1:])):
			line2write=[float(listOfCell2read[idxLine+1][item]) for item in [2, 3]]  + [np.median([float(listOfCell2read[idxLine+1][itemMed]) for itemMed in [4, 5, 6, 7]]) ] + [np.median([float(listOfCell2read[idxLine+1][itemMed]) for itemMed in [8, 9, 10, 11]]) ]
			content2Store[idxCell+idxLine+1] += line2write 
		idxCell+=len(listOfCell2read[1:])

	## Export the cell measurements
	with open(outputDetPath+ "/Cell_Measures.csv",'wb') as mainFile:
		writer= csv.writer(mainFile, dialect='unixpwd')
		content2Store[0]+= [listOfCell2read[0][item] for item in [2, 3]]
		content2Store[0]+=['Average_Texture_SumEntropy', 'Average_Texture_SumVariance', 'Class']
		writer.writerow(content2Store[0])
		listOfCell=content2Store[1:]
		listOfCell=sorted(listOfCell, key=lambda x: x[0])
		idxLine=1
		for currRow in listOfCell:
			if sum(map(int,currRow[4:12]))!=0:
				currRow+= [listOfClass[idxLine][0]]
				idxLine+=1
				writer.writerow(currRow)

	## Load the cell measurements
	## Can be avoid by storing in previous step !! TODO
	with open(outputDetPath+ "/Cell_Measures.csv",'r') as csvfile:
		reader= csv.reader(csvfile)	
		content2Store = list(reader)

	## Load the fat cells LDs size distribution
	with open(outputDetPath+ "/GFPSize_Clean.csv",'r') as csvfile:
		reader= csv.reader(csvfile)	
		listOfCleanedNum2read = list(reader)

	listForCell.append([item for item in content2Store[0][0:3]]+ ['Average_Intensity_IntegratedIntensity_GFP', 'Average_Intensity_IntegratedIntensity_NUC'] + [item for item in content2Store[0][14:16]]+ ['Percent_Class0','Percent_Class1', 'BeforeFiltering', 'AfterFiltering', 'AfterCleaning', 'NucW/Drop', 'Percent With',	'Percent Without'])

	## Prepare the header for writing
	header2Write= content2Store[0]
	header2Write=[header2Write[item] for item in [1,3,12, 12]]
	listOfIntAfter.append(header2Write)
	for content in content2Store[1:]:
		content = [float(content[item]) for item in [1,3,12, 12]]
		content[3]=math.log(content[2],2)
		listOfIntAfter.append(content) 

	## Prepare the content for writing
	nbClass=int(max([item[0] for item in listOfClass[1:]]))
	for dbIdx in range(len(listAnalysis)):
		currSubPath=outputCellProfilerPath+ "temp/Analysis_"+str(dbIdx)+"/"
		if os.path.isdir(currSubPath):	
			with open(currSubPath+ "/MyExpt_FilteredNuc.csv",'r') as csvfile:
				reader= csv.reader(csvfile)	
				list2read = list(reader)
			with open(currSubPath+ "/MyExpt_ExpNuclei.csv",'r') as csvfile:
				reader= csv.reader(csvfile)	
				listOfNum2read = list(reader)


		if dbIdx==0:
			header2Write= listOfNum2read[0]
			header2Write=[header2Write[item] for item in idxColumnBef]
			listOfIntBef.append(header2Write)

		idxImgMax=int(list2read[len(list2read)-1][0])

		for content in listOfNum2read[1:]:
			content = [float(content[item]) for item in idxColumnBef]
			content[0]= content[0]+ idxImg
			listOfIntBef.append(content) 


		for currIdxImg in range(1, idxImgMax+1):
			if [item[1] for item in content2Store[1:] if item[1]==str(currIdxImg+ idxImg)]:
				allInfo= [item[0:3] for item in content2Store[1:] if item[1]==str(currIdxImg+ idxImg)][0]
				for i in range(12,16):
					currInfo= [float(item[i]) for item in content2Store[1:] if item[1]==str(currIdxImg+ idxImg)]
					allInfo+=[float(sum(currInfo))/float(len(currInfo))]

				currInfo= [item[0] for item in listOfClass[1:] if item[2]==str(currIdxImg+ idxImg)]
				for i in range(nbClass+1):
					allInfo+=[float(len([int(item) for item in currInfo if item==str(i)])*100)/float(len(currInfo))]	

				SumNucBef= max([int(item[1]) for item in listOfNum2read[1:] if item[0]==str(currIdxImg)])
				SumNucWith= max([int(item[1]) for item in list2read[1:] if item[0]==str(currIdxImg)])
				SumNucClean= len([int(item[1]) for item in listOfCleanedNum2read[1:] if item[1]==str(currIdxImg + idxImg)])
				PercWith= 100*float(SumNucClean)/SumNucBef
			 	PercWithout= 100*float(SumNucBef-SumNucClean)/SumNucBef
				listOfNum.append([currIdxImg + idxImg]+ [SumNucBef] + [SumNucWith] + [SumNucClean] + [SumNucBef-SumNucClean] + [PercWith] + [PercWithout])
				

				allInfo+=[SumNucBef] + [SumNucWith] + [SumNucClean] + [SumNucBef-SumNucClean] + [PercWith] + [PercWithout]
				listForCell.append(allInfo)
			else:
				print ('%s saturated'%(currIdxImg+ idxImg))

		idxImg=idxImg + idxImgMax

	## Write the per-well measurements
	with open(outputDetPath+ "/Well_Measures.csv",'wb') as mainFile:
		writer= csv.writer(mainFile, dialect='unixpwd')
		for i in listForCell:
			writer.writerow(i)

	## Write the measurements relative to the intensity of the fat cells
	with open(outputDetPath+ "/NUC_Intensity.csv",'wb') as mainFile:
		writer= csv.writer(mainFile, dialect='unixpwd')
		writer.writerow(listOfIntAfter[0])
		listOfIntAfter=listOfIntAfter[1:]
		listOfIntAfter=sorted(listOfIntAfter, key=lambda x: x[0])
		for i in listOfIntAfter:
			writer.writerow(i)

	## Write the measurements relative to the intensity of all the cells
	with open(outputDetPath+ "/NUC_IntensityBefore.csv",'wb') as mainFile:
		writer= csv.writer(mainFile, dialect='unixpwd')
		writer.writerow(listOfIntBef[0])
		listOfIntBef=listOfIntBef[1:]
		listOfIntBef=sorted(listOfIntBef, key=lambda x: x[0])
		for i in listOfIntBef:
			writer.writerow(i)

	## Write the measurements relative to the number of nuclei
	listOfNum=sorted(listOfNum, key=lambda x: x[0])
	with open(outputDetPath+ "/NUC_Numeration.csv",'wb') as mainFile:
		writer= csv.writer(mainFile, dialect='unixpwd')	
		content= ['ImageNumber', 'BeforeFiltering', 'AfterFiltering', 'AfterCleaning', 'NucW/Drop', 'Percent With',	'Percent Without']
		writer.writerow(content)
		for i in listOfNum:
			writer.writerow(i)

		

	print('Fusion done')






	# listOfPlate=[range(1, 136), range(136, 280), range(280, 396)]

	# firstList2plot=[[] for i in range(len(listOfPlate))]
	# secondList2plot=[[] for i in range(len(listOfPlate))]
	# thirdList2plot=[[] for i in range(len(listOfPlate))]

	# for nbPat in range(len(listOfPlate)):
	# 	idxPat=listOfPlate[nbPat]
	# 	for currIdx in idxPat:
	# 		firstList2plot[nbPat].append([float(item[2]) for item in listOfIntBef if float(item[0])==currIdx])
	# 		secondList2plot[nbPat].append([float(item[3]) for item in listOfIntBef if float(item[0])==currIdx])
	# 		thirdList2plot[nbPat].append([float(item[4]) for item in listOfIntBef if float(item[0])==currIdx])

	# plotDensityDistribution(firstList2plot, "Integrated")
	# plotDensityDistribution(secondList2plot, "Max")
	# plotDensityDistribution(thirdList2plot, "Mean")