# plotFeatures.py
#######################################################
#
#	Plot the three main features:
#	- Percentage of MS/FL, Cell classification
#	- Percentage of Fat/non-fat cell
#	- Intensity of nuclei in the Hoescht channel
#
#######################################################

import os

import multiprocessing
import Queue
import threading

import csv
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


from scipy.stats import gaussian_kde
import numpy as np
import math

import settings

from PIL import Image


import sys


def workerOnList():
	while not q.empty():
		# Get the index of the plate from the queue
		idxPlate = q.get()
		# Get the total number of images concerning this plate from the csv file

		listofImg=[int(item[1]) for item in listOfWellFeatures if str(item[0]).zfill(2)==str(idxPlate).zfill(2)]
		nbImg= max(listofImg)
		debImg= min(listofImg)

		for idxImg in range(debImg, nbImg+1):
			print ('P%s, Img%s'%(idxPlate,idxImg))
			if [item[1] for item in listOfWellFeatures if item[1]==str(idxImg)]:
				idxWell= [str(item[2]) for item in listOfWellFeatures if item[1]==str(idxImg)][0]	

				slicesOfSize= [item[7:9] for item in listOfWellFeatures if item[1]==str(idxImg)][0]
				slicesOfSize= [float(item) for item in slicesOfSize]

				slicesOfCount= [item[13:15] for item in listOfWellFeatures if item[1]==str(idxImg)][0]
				slicesOfCount= [float(item) for item in slicesOfCount]

				# listOfIntensity= [math.log(float(item[13]),2) for item in listOfCellFeatures if item[1]==str(idxImg)]
				listOfIntensity= [float(item[13]) for item in listOfCellFeatures if item[1]==str(idxImg)]
				if len(listOfIntensity)>1:
					density = gaussian_kde(listOfIntensity)
					xs = np.linspace(min(listOfIntensity)-0.5 ,max(listOfIntensity)+0.5,200)
					density.covariance_factor = lambda : .25
					density._compute_covariance()
					finalListToPlot[idxPlate-1].append([idxWell, slicesOfSize, slicesOfCount, [xs, density(xs)]])
				else:
					finalListToPlot[idxPlate-1].append([idxWell, slicesOfSize, slicesOfCount, [0,0]])


				idxPatient=  next(index for index, value in enumerate(listOfPatient) if idxPlate in value)
				if idxWell in listOfPos:
					list2plotIntPos[idxPatient].append(listOfIntensity)
				else:	
					list2plotIntNeg[idxPatient].append(listOfIntensity)
			else:
				print ('P%s, Img%s saturated'%(idxPlate,idxImg))



			


			

def plotPieOfPlate(idxPlate, finalListIndex, colorsToUse, labelsToUse):
		f, axarr = plt.subplots(6, 12, sharex=True, sharey=True, figsize=(80, 40), squeeze=True)
		for idxImg in range(len(finalListToPlot[idxPlate-1])):
			slizesOfPie= finalListToPlot[idxPlate-1][idxImg][finalListIndex]
			idxWell = finalListToPlot[idxPlate-1][idxImg][0]
			if idxWell in ['C02', 'E02','G02']:
				colorsToUse = ['gold', 'darkgray']
			elif idxWell in ['B10', 'E05','G08']:
				colorsToUse = ['mediumblue', 'darkgray']
			elif idxWell in ['B02', 'C11', 'D02', 'D06', 'F02', 'F09', 'C02', 'E02','G02']:
				colorsToUse = ['darkred', 'darkgray']
			axarr[int(listOfLetters.index(idxWell[0])), int(idxWell[1:])-1].pie(slizesOfPie, labels=labelsToUse[:len(slizesOfPie)+1], colors=colorsToUse, autopct='%1.1f%%', shadow=True, startangle=90)
		# Set aspect ratio to be equal so that pie is drawn as a circle.
		plt.axis('equal')
		plt.savefig(outputDetPath+ "/piePlots/pie"+ labelsToUse[0] + "_" + str(idxPlate).zfill(2)+".png", bbox_inches='tight')
		plt.close()


def plotDoubleBarStats(list2plotPos, list2plotNeg, color2plot, label2plot, legendFlag ):

	font = {'size'   : 80}

	matplotlib.rc('font', **font)

	list2plotPos=[item for item in list2plotPos if item ]
	list2plotNeg=[item for item in list2plotNeg if item ]
	
	meanListPos= [float(sum(item))/len(item) for item in list2plotPos]
	meanListNeg= [float(sum(item))/len(item) for item in list2plotNeg]

	stdPos=[]*len(list2plotPos)
	stdNeg=[]*len(list2plotPos)

	for item in range(len(list2plotPos)):
		variance = map(lambda x: (x - meanListPos[item])**2, list2plotPos[item])
		stdPos.append(math.sqrt(sum(variance)/float(len(variance))))
		variance = map(lambda x: (x - meanListNeg[item])**2, list2plotNeg[item])
		stdNeg.append(math.sqrt(sum(variance)/float(len(variance))))


	ind = np.arange(len(list2plotPos))
	width = 0.35      
	fig, ax = plt.subplots(figsize=(30, 30),dpi=600)
	ax.tick_params(axis='both', which='major',  size=80)

	rects1 = ax.bar(ind, meanListNeg,                  # data
	                width,                          # bar width
	                color= color2plot,        # bar colour
	                yerr=stdNeg,                  # data for error bars
	                error_kw={'ecolor':'darkgray',    # error-bars colour
	                          'linewidth':6})       # error-bar width

	rects2 = ax.bar(ind + width, meanListPos, 
	                width, 
	                color='darkgray', 
	                yerr=stdPos, 
	                error_kw={'ecolor':color2plot,
	                          'linewidth':6})



	axes = plt.gca()

	botLimit=0
	topLimit= max([max([item for item in meanListPos])]+ [max([item for item in meanListNeg])]) + abs(2*max([max([item for item in stdPos])]+[max([item for item in stdNeg])]))

	ax.set_ylabel( label2plot, fontsize=100)

	ax.set_xticks(ind + width)
	ax.set_xticklabels(('Patient 1', 'Patient 2', 'Patient 3'), fontsize=100)
	
	if legendFlag:
		ax.legend((rects1[0], rects2[0]), ('Neg Control', 'Pos Control'), fontsize=80)


	data_x=[]
	for rect in rects1:
		data_x.append(rect.get_x() + rect.get_width()/2)

	markPos= [max(item) for item in list2plotNeg]
	data_y = [markPos[0], markPos[1],markPos[2]]
	print ('maxNeg:')
	print data_y
	if topLimit< max(data_y):
		topLimit= max(data_y)
	ax.plot(data_x, data_y, linestyle=' ', color='darkgray', marker='o', markersize=40)

	markPos= [min(item) for item in list2plotNeg]
	data_y = [markPos[0], markPos[1],markPos[2]]
	print ('minNeg:')
	print data_y
	if botLimit> min(data_y):
		botLimit= min(data_y)
	ax.plot(data_x, data_y, linestyle=' ', color='darkgray', marker='s', markersize=40)

	data_x=[]
	for rect in rects2:
		data_x.append(rect.get_x() + rect.get_width()/2)
	markPos= [max(item) for item in list2plotPos]
	data_y = [markPos[0], markPos[1],markPos[2]]
	print ('maxPos:')
	print data_y
	if topLimit< max(data_y):
		topLimit= max(data_y)
	ax.plot(data_x, data_y, linestyle=' ', color='lightskyblue', marker='o', markersize=40)

	markPos= [min(item) for item in list2plotPos]
	data_y = [markPos[0], markPos[1],markPos[2]]
	print ('minPos:')
	print data_y
	if botLimit> min(data_y):
		botLimit= min(data_y)
	ax.plot(data_x, data_y, linestyle=' ', color='lightskyblue', marker='s', markersize=40)

	axes.set_ylim([math.floor(botLimit), math.ceil(topLimit)]) 



def autolabel(ax, rects):
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
                '%.2f' % float(height),
                ha='center',            # vertical alignment
                va='bottom',             # horizontal alignment
                fontsize=80
                )



#################### Main ####################
def plotFeat():
## Plotting !

	## Initialise global variables
	global q

	global listOfWellFeatures
	global listOfCellFeatures

	global listOfPatient
	global listOfLetters
	global listOfWells	


	global outputDetPath
	global finalListToPlot
	

	global listOfWells
	global colorsOfPie
	global listOfPos
	global gsSizePie
	global gsCountPie
	global gsHisto

	global list2plotIntPos
	global list2plotIntNeg

	global nbClass

	## Initialise paths
	outputDetPath=settings.pathList[3]

	## Initialise variables
	listOfPos= ['C02', 'E02','G02']

	listOfLetters=['B', 'C', 'D', 'E', 'F', 'G']
	listOfWells=[]
	for letters in listOfLetters:
		for val in range(2,12):
			listOfWells.append(letters+ str(val).zfill(2))

	listOfWells=['B02', 'C11', 'D02', 'D06', 'F02', 'F09', 'C02', 'E02','G02']

	colorsOfPie = ['yellowgreen', 'chocolate', 'lightskyblue', 'mistyrose' ]
	
	listOfPatient=[range(1, 16), range(16, 32), range(32, 45)]


	list2plotNumPos=[[] for i in range(len(listOfPatient))]
	list2plotNumNeg=[[] for i in range(len(listOfPatient))]

	list2plotSizePos=[[] for i in range(len(listOfPatient))]
	list2plotSizeNeg=[[] for i in range(len(listOfPatient))]

	list2plotIntPos=[[] for i in range(len(listOfPatient))]
	list2plotIntNeg=[[] for i in range(len(listOfPatient))]

	list2plotVarPos=[[] for i in range(len(listOfPatient))]
	list2plotVarNeg=[[] for i in range(len(listOfPatient))]

	## Initialise threads
	nProc = multiprocessing.cpu_count()
	q = Queue.Queue()
	threadsList = [threading.Thread(target=workerOnList) for i in range(0, nProc)]


	if not os.path.isdir(outputDetPath+ "/featuresPlot/"):
		os.mkdir(outputDetPath+ "/featuresPlot/")
	if not os.path.isdir(outputDetPath+ "/piePlots/"):
		os.mkdir(outputDetPath+ "/piePlots/")


### OpenFile ###
	# Plate	ImageNumber	Well	Average_Intensity_IntegratedIntensity_GFP	Average_Intensity_IntegratedIntensity_NUC	Average_Texture_SumEntropy	Average_Texture_SumVariance	Percent_Class0	Percent_Class1	BeforeFiltering	AfterFiltering	AfterCleaning	NucW/Drop	Percent With	Percent Without
	with open(outputDetPath+ "/Well_Measures.csv",'r') as csvfile:
		reader= csv.reader(csvfile)	
		list2read = list(reader)
	headerOfSizeDistribution= list2read[0]
	listOfWellFeatures= list2read[1:]
	# Plate	ImageNumber	Well	ObjectNumber	bin_0_5	bin_5_10	bin_10_15	bin_15_20	bin_20_25	bin_25_30	bin_30_35	bin_35_40	Intensity_IntegratedIntensity_GFP	Intensity_IntegratedIntensity_NUC	Average_Texture_SumEntropy	Average_Texture_SumVariance	Class
	with open(outputDetPath+ "/Cell_Measures.csv",'r') as csvfile:
		reader= csv.reader(csvfile)	
		list2read = list(reader)
	headerOfSizeDistribution= list2read[0]
	listOfCellFeatures= list2read[1:]

	nbPlate= max([int(item[0]) for item in listOfWellFeatures])
	nbClass= max([int(item[16]) for item in listOfCellFeatures]) +1

	finalListToPlot=[[] for i in range(nbPlate)]


### Measure/Analyse ###
	for idxPlate in range(1, nbPlate+1):	
		q.put(idxPlate)

	# start all threads        
	for thread in threadsList:
		thread.start()
			
	# join all threads
	for thread in threadsList:
		thread.join()


	## FinalListToPlot = [Patient1:[[Img, bla,], [Img, bla,]], [Img, bla,]], [Patient2] , [Patient3]]
	for idxPatient in range(len(listOfPatient)):
		plateList= [item for item in finalListToPlot if finalListToPlot.index(item)+1 in listOfPatient[idxPatient]]
		list2plotSizePos[idxPatient]= [item[1][0] for featureList in plateList for item in featureList if item[0] in ['C02', 'E02','G02']]
		list2plotSizeNeg[idxPatient]= [item[1][0] for featureList in plateList for item in featureList if item[0] not in ['C02', 'E02','G02']]
		list2plotNumPos[idxPatient]= [item[2][0] for featureList in plateList for item in featureList if item[0] in ['C02', 'E02','G02']]
		list2plotNumNeg[idxPatient]= [item[2][0] for featureList in plateList for item in featureList if item[0] not in ['C02', 'E02','G02']]


	plotFlag=False
### Plot ###
	if plotFlag:
		for idxPlate in range(1, nbPlate+1):
			f, axarr = plt.subplots(1, len(listOfWells)*3, figsize=(340, 20))
			for idxImg in range(len(finalListToPlot[idxPlate-1])):
				idxWell=finalListToPlot[idxPlate-1][idxImg][0]
				slicesOfSize=finalListToPlot[idxPlate-1][idxImg][1]
				slicesOfCount=finalListToPlot[idxPlate-1][idxImg][2]
				[X, D]=finalListToPlot[idxPlate-1][idxImg][3]

				axarr[ int(listOfWells.index(idxWell))].pie(slicesOfSize, colors=colorsOfPie[:2], autopct='%1.1f%%', shadow=True, startangle=90)
				axarr[ int(len(listOfWells)+ listOfWells.index(idxWell))].pie(slicesOfCount, labels=['With', 'Without'], colors=colorsOfPie[2:], autopct='%1.1f%%', shadow=True, startangle=90)

				axarr[ int(2*len(listOfWells)+ listOfWells.index(idxWell))].plot(X, D)
			plt.savefig(outputDetPath+ "/featuresPlot/Features_Pl"+str(idxPlate).zfill(2)+".png", bbox_inches='tight')
			plt.close()	
		
			plotPieOfPlate(idxPlate, 1, [ 'lightskyblue', 'mistyrose'], ['MS', 'FL'])
			plotPieOfPlate(idxPlate, 2, ['yellowgreen', 'chocolate'], ['With', 'Without'])
		

		## Fusion the plates Features
		listFiles= [outputDetPath+ "/featuresPlot/"+x for x in os.listdir(outputDetPath+ "/featuresPlot/") if x.startswith('Features_')]
		images = map(Image.open, listFiles)
		widths, heights = zip(*(i.size for i in images))
		total_width = max(widths)
		max_height = sum(heights)
		new_im = Image.new('RGB', (total_width, max_height))
		y_offset = 0
		for im in images:
		  new_im.paste(im, (0,y_offset))
		  y_offset += im.size[1]

		new_im.save(outputDetPath+ "/featuresPlot/Features.png")

	for idxPat in range(len(listOfPatient)):
		with open(outputDetPath+ "/featuresPlot/Size_Pat" + str(idxPat) + ".csv",'wb') as csvfile:
			writer= csv.writer(csvfile, dialect='unixpwd')	
			writer.writerow(list2plotSizePos[idxPat]) 	
			writer.writerow(list2plotSizeNeg[idxPat]) 	
		with open(outputDetPath+ "/featuresPlot/Num_Pat" + str(idxPat) + ".csv",'wb') as csvfile:
			writer= csv.writer(csvfile, dialect='unixpwd')	
			writer.writerow(list2plotNumPos[idxPat]) 	
			writer.writerow(list2plotNumNeg[idxPat]) 	


	for i in range(len(listOfPatient)):
		list2plotIntPos[i]= [item for sublist in list2plotIntPos[i] for item in sublist]
		list2plotIntNeg[i]= [item for sublist in list2plotIntNeg[i] for item in sublist]

	plotDoubleBarStats(list2plotNumPos, list2plotNumNeg, 'lightskyblue', 'Percentage fat cell',0)
	plt.savefig(outputDetPath+ "/featuresPlot/withFeatures.png")
	plotDoubleBarStats(list2plotSizePos, list2plotSizeNeg, 'lightskyblue', 'Percentage FL class',0)
	plt.savefig(outputDetPath+ "/featuresPlot/sizeFeatures.png")
	plotDoubleBarStats(list2plotIntPos, list2plotIntNeg, 'lightskyblue', 'Average integrated intensity',1)
	plt.savefig(outputDetPath+ "/featuresPlot/intFeatures.png")

