# clusterDroplets.py
#######################################################
#
#	Classification based on the size distribution:
#	- Load collection of vectors
#	- For each patient, take all images, and cluster every cells in two classes
#		-> Many small (MS), lot of LDs in small bins [0-15]
#		-> Few large (FL), heterogeneous size distribution with presence of big LDs
#	- Export the collection per image with new column for the class
#
#######################################################

import settings
import numpy as np
import sys
import os
import csv
import random 
import Queue
import threading
import math

from scipy import stats
from matplotlib.mlab import PCA
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import cv2
from pyemd import emd


## Thread
## Add the class value in the first column of the correct row 
## based on the classification
def workerCluster():
	while not q.empty():
		idxClass = q.get()
	 	currClass=[item.tolist() for item in classifK[idxClass]]

		idxRowInList= [item for item, sublist in enumerate(listOfSize) if sublist in currClass]
		for itemIterate in idxRowInList:
			content=  [str(idxClass)]+ list2Process[itemIterate]
			list2write.append(content)


## Function to find centroids of the current cluster ##
## Input: X: Data to cluster
## 		K: Number of class
## 		oldmu: Position of the centroids 
## Output: mu: New evaluated position of the centroids
## 		clusters: New evaluated cluster
def find_centers(X, K, oldmu):
    # Initialize to K random centers
	if not oldmu:
		oldmu = random.sample(X, K)
	mu = random.sample(X, K)

	while not has_converged(mu, oldmu):
		oldmu = mu
		# Assign all points in X to clusters
		clusters = cluster_points(X, mu)
		# Reevaluate centers
		mu = reevaluate_centers(oldmu, clusters)
	return(mu, clusters)


## Measure intra-cluster distance based on Earth's Moving##
## Input: X: Old cluster
## 		 mu: Position of the centroids 
## Output: clusters: New cluster
def cluster_points(X, mu):
	clusters  = {}
	for x in X:
		bestmukey= min([(i[0], emd(np.asarray(x, np.float), np.asarray(mu[i[0]], np.float) , distance_matrix)) for i in enumerate(mu)], key=lambda t:t[1])[0]
		try:
			clusters[bestmukey].append(x)
		except KeyError:
			clusters[bestmukey] = [x]
	return clusters


## Measure new positions of the centroids##
## Input: mu: Position of the centroids 
## 		 clusters: Data of the cluster 
## Output: newmu: New evaluated position of the centroids
def reevaluate_centers(mu, clusters):
    newmu = []
    keys = sorted(clusters.keys())
    for k in keys:
        newmu.append(np.mean(clusters[k], axis = 0))

    return newmu


## Check if the position of the centroids changed, i.e., if the classification converge##
## Input: mu: New position of the centroids 
## 		 oldmu: Old position of the centroids  
## Output: boolean, either the position changed (True) or not (False)
def has_converged(mu, oldmu):
    return (set([tuple(a) for a in mu]) == set([tuple(a) for a in oldmu]))



## Create a distance matrix for the Earth Moving metric##
## Input: NBO: Size of the vectors for the distance, i.e., number of bins = SIFT's Orientation-dimension
## Output: D: matrix of distance
def defineDistanceMatrix(NBO):
	XNBP= 1 # SIFT's X-dimension
	YNBP= 1 # SIFT's Y-dimension

	N= XNBP*YNBP*NBO
	thresh= 2

	D= np.zeros((N, N), dtype=np.float)
	i= 0 
	for y1 in range(YNBP):
		for x1 in range(XNBP):
			for nbo1 in range(NBO): 
				j= 0 
				for y2 in range(YNBP):
					for x2 in range(XNBP):
						for nbo2 in range(NBO):
							D[i,j] = math.sqrt(math.pow((y1-y2),2) + math.pow((x1-x2),2))  +  min(abs(nbo1-nbo2), NBO-abs(nbo1-nbo2))
							j=j+1 
				i= i+1 

	D[D>thresh]= thresh

	return D



#################### Main ####################

def getClusterOfDroplets(nbClass):
## Clustering !

	## Initialise variables
	global q
	global list2write
	global classifK
	global listOfSize
	global list2Process
	global distance_matrix

	## Initialise paths	
	outputDetPath=settings.pathList[3]

	listOfPatient=[range(1, 16), range(16, 32), range(32, 45)]
	listOfPatient=[map(str, map("{:02d}".format, item)) for item in listOfPatient]

	## Read size distribution
	with open(outputDetPath+ "/GFPSize_Clean.csv",'r') as csvfile:
		reader= csv.reader(csvfile)	
		list2read = list(reader)
	header2Write= list2read[0]
	list2read= list2read[1:]

	list2write=[]
	centrK=[]
	## Initialise a "Many small" class instance
	prevClass0=[10, 5, 0, 0, 0, 0, 0, 0]

	## Initialise distance matrix for Earth's Moving metric
	distance_matrix = defineDistanceMatrix(len(prevClass0))


	savecentrK=[[] for i in range(len(listOfPatient))]
	# Cluster Per patient
	for idxPatient in range(len(listOfPatient)):

		list2Process= [item for item in list2read if item[0] in ([str(idxItem) for idxItem in listOfPatient[idxPatient]])]
		listOfSizeOri=[map(int,item[4:]) for item in list2Process]

		listOfSize= np.asarray(listOfSizeOri)
		listOfSize= listOfSize.tolist()

		if listOfSize:
			X = np.asarray(listOfSize)

			(centrK, classifK)=find_centers(X,nbClass, centrK)

			nProc=len(classifK)
			## Assure that class0 is the "MS" class
			if prevClass0:
				if  sum([abs(a_i - b_i) for a_i, b_i in zip(centrK[1].tolist(), prevClass0)]) < sum([abs(a_i - b_i) for a_i, b_i in zip(centrK[0].tolist(), prevClass0)]):
					temp= centrK[0]
					centrK[0]= centrK[1]
					centrK[1]=temp
					temp= classifK[0]
					classifK[0]= classifK[1]
					classifK[1]=temp

			print centrK
			savecentrK[idxPatient].append(centrK)
			prevClass0= centrK[0].tolist()
	
			## Initialise and fill up the queue for threading
			threadsList = [threading.Thread(target=workerCluster) for i in range(0, nProc)]
			q = Queue.Queue()

			[q.put(query) for query in range(len(classifK))]

			# start all threads        
			for thread in threadsList:
				thread.start()
					
			# join all threads
			for thread in threadsList:
				thread.join()

	## If the sizes of the input and output lists are different
	if len(list2write) != len(list2read):
		print('Clustering error')
		sys.exit(1)



	print('Writing csv...')


	## Export classification of the cell, sorted by class
	list2write=sorted(list2write, key=lambda x: (x[0],x[1]))
	with open(outputDetPath+ "/GFPSize_CLASSED.csv",'wb') as mainFile:
		writer= csv.writer(mainFile, dialect='unixpwd')	
		writer.writerow(['Class'] +header2Write)
		for content in list2write:
			writer.writerow(content) 		

	## Export classification of the cell, sorted by plate/ image/ cell
	list2write=sorted(list2write, key=lambda x: (x[1], x[2], str(x[4]).zfill(6)))
	with open(outputDetPath+ "/GFPSize_UNCLASSED.csv",'wb') as mainFile:
		writer= csv.writer(mainFile, dialect='unixpwd')	
		writer.writerow(['Class'] +header2Write)
		for content in list2write:
			writer.writerow(content) 						

	## Export positions of the centroid for the classes
	# with open(outputDetPath+ "/cluster_Center.csv",'wb') as mainFile:
	# 	writer= csv.writer(mainFile, dialect='unixpwd')	
	# 	for content in savecentrK:
	# 		for subContent in content:
	# 			writer.writerow(subContent) 			