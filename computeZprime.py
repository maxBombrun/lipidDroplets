# computeZprime.py
#######################################################
#
#	For all features:
#		- Compute the ZPrime factor
#		- Compute the Receiver Operator Characteristic (ROC) and Area Under the Curve (AUC) 
#		- Plot the ROC with AUC and ZPrime of the feature
#
#######################################################

import os
import csv
import matplotlib
import matplotlib.pyplot as plt
import settings
import math
import numpy as np

from sklearn.metrics import roc_curve, auc
from scipy.stats import gaussian_kde


## Compute Zprime factor ##
## Input: listOfPos: list of data from the positive control
##		  listOfNeg: list of data from the negative control
## Output: Z: Zprime factor
def computeZPrime(listOfPos, listOfNeg):
	listOfPos=[item for item in listOfPos if item is not None]
	listOfNeg=[item for item in listOfNeg if item is not None]

	meanPos= sum(listOfPos)/float(len(listOfPos))
	variance = map(lambda x: (x - meanPos)**2, listOfPos)
	stdPos = math.sqrt(sum(variance)/float(len(variance)))

	meanNeg= sum(listOfNeg)/float(len(listOfNeg))
	variance = map(lambda x: (x - meanNeg)**2, listOfNeg)
	stdNeg = math.sqrt(sum(variance)/float(len(variance)))

	Z= 1- (float(3*(stdPos+stdNeg))/abs(meanPos-meanNeg))

	return Z



#################### Main ####################

def getZprime():

### Init ###	

	outputDetPath=settings.pathList[3]

	with open(outputDetPath+ "/Well_Measures.csv",'r') as csvfile:
		reader= csv.reader(csvfile)	
		listOfWellFeatures = list(reader)
	listOfHeader= listOfWellFeatures[0]
	listOfWellFeatures = listOfWellFeatures[1:]
	
	listOfHeader= listOfHeader[3:9] + listOfHeader[13:15]

	listIdxPos = [item[3:9] + item[13:15] for item in listOfWellFeatures if item[2] in ['C02', 'E02','G02']]
	listIdxNeg = [item[3:9] + item[13:15] for item in listOfWellFeatures if not item[2] in ['C02', 'E02','G02']]
	
	font = {'size'   : 60}
	matplotlib.rc('font', **font)

	for i in range(len(listIdxNeg[0])):
		fig, ax = plt.subplots(figsize=(40, 40))

		listNegFeature = [float(item[i]) for item in listIdxNeg]
		listPosFeature = [float(item[i]) for item in listIdxPos]

		maxX= max(max(listNegFeature, listPosFeature))
		minX= min(min(listNegFeature, listPosFeature))

		densityNeg = gaussian_kde(listNegFeature)
		xsNeg = np.linspace(minX-0.5 ,maxX+0.5,200)
		densityNeg.covariance_factor = lambda : .25
		densityNeg._compute_covariance()


		densityPos = gaussian_kde(listPosFeature)
		xsPos = np.linspace(minX-0.5 ,maxX+0.5,200)
		densityPos.covariance_factor = lambda : .25
		densityPos._compute_covariance()
		
		maxY= max(max(densityNeg(xsNeg).tolist(), densityPos(xsPos).tolist()))
		minY= min(min(densityNeg(xsNeg).tolist(), densityPos(xsPos).tolist()))

		Z= computeZPrime(listPosFeature, listNegFeature)
		print str(listOfHeader[i])
		print('Z= %s'%Z)

		red_line, blue_line =plt.plot(xsNeg, densityNeg(xsNeg),'r', xsPos, densityPos(xsPos), 'b')
		plt.text(3*maxX/4, 3*maxY/4, r'$Z= %s$'%Z, fontsize=40)
		plt.legend([red_line, blue_line],["Neg", "Pos"], fontsize=40)
		if not os.path.isdir(outputDetPath+ "/features/"):
			os.mkdir(outputDetPath+ "/features/")	
		plt.savefig(outputDetPath+ "/features/feat_" + str(listOfHeader[i]) +".png", bbox_inches='tight')
		plt.close()	

		FPR=[]
		TPR=[]
		for seuil in range(0, 100, 2):
			FPR.append(float(len([item for item in listNegFeature if float(item) > seuil]))/ float(len(listNegFeature)))
			TPR.append(float(len([item for item in listPosFeature if float(item) > seuil]))/ float(len(listPosFeature)))

		FPR=np.asarray(FPR)
		TPR=np.asarray(TPR)
		AUC = auc(FPR,TPR)
		
		if i==1:
			nameFeat='Av. integ. intensity'
		elif i==4:
			nameFeat='FL class'
		elif i==7:
			nameFeat='Fat cell'
		else:
			nameFeat='ROC curve'

		font = {'size'   : 80}

		matplotlib.rc('font', **font)
		plt.figure(figsize=(30, 30), dpi=600)
		plt.xlabel("FPR", fontsize=100)
		plt.ylabel("TPR", fontsize=100)
		plt.title("%s, AUC=%0.2f" % (str(nameFeat), (AUC)), fontsize=100)
# 
		plt.plot(FPR, TPR, color='green', linewidth=8)

		x = [0.0, 1.0]
		plt.plot(x, x, linestyle='dashed', color='red', linewidth=8, label='random')
		plt.xlim(0.0, 1.0)
		plt.ylim(0.0, 1.0)
		plt.tight_layout()
		plt.savefig(outputDetPath+ "/features/roc_" + str(listOfHeader[i]) +".png", bbox_inches='tight')

		plt.close()	