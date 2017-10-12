# segmentNucAndGFP.py
#######################################################
#
#	Segmentation pipeline for :
#	- Load original images
#	- Detect lipid droplets (LDs) and nuclei
#	- Save the "cleaned" images
#
#######################################################

import os

import cv2
import numpy as np
from scipy import ndimage

import multiprocessing
import Queue
import threading

import settings

import csv


## Thread
## Get the paths of the images and  
## Process them to extract clean masks of the LDs and the nuclei
## Save the images as tif for later use in CellProfiler
def workerSeg():
	while not q.empty():
		currSubPath= q.get()
		print ("%s\n" % currSubPath)
		listFiles= [x for x in os.listdir(currSubPath) if x.endswith('.tif')]
		gfpList=[x for x in listFiles if x.startswith("GFP")]
		nucList=[x for x in listFiles if x.startswith("Hoechst")]
		for gfpName in gfpList:
			gfpPath=currSubPath+gfpName

			idxPlate=gfpPath.rsplit('plate')
			idxPlate= idxPlate[1].rsplit('/')
			idxPlate= idxPlate[0].zfill(2)
			idxWell=gfpPath.rsplit('Well ')
			idxWell= idxWell[1].rsplit('/')
			idxWell= idxWell[0]

			print("Detecting P%s_%s\n"% (idxPlate, idxWell))
			imgGFP = cv2.imread(gfpPath, -1)
			cv2.imwrite(inputCellProfilerPath+'P'+idxPlate+'_'+idxWell+"_GFP.tif",imgGFP, [cv2.CV_16U])	

			imgMaskGFP=detectGFP(imgGFP, 20)
			imgDet=imgMaskGFP*imgGFP
			imgDet=imgDet.astype(np.uint16)

			cv2.imwrite(outputDetPath +'P'+idxPlate+'_'+idxWell+"_GFP_CL.tif",imgDet, [cv2.CV_16U])		
			cv2.imwrite(outputDetPath +'P'+idxPlate+'_'+idxWell+"_GFP_MASK.tif",imgMaskGFP*255)


		for nucName in nucList:
			nucPath=currSubPath+nucName
			print("%s in progress"% nucPath)
			imgNUC=cv2.imread(nucPath,  -1 )

			imgMaskNUC=detectNucs(imgNUC, 40)  
			imgDet=imgMaskNUC*imgNUC   
			imgDet=imgDet.astype(np.uint16)

			cv2.imwrite(outputDetPath +'P'+idxPlate+'_'+idxWell+"_NUC_MASK.tif",imgMaskNUC*255) 
			cv2.imwrite(outputDetPath +'P'+idxPlate+'_'+idxWell+"_NUC_CL.tif",imgDet, [cv2.CV_16U])
			cv2.imwrite(inputCellProfilerPath+'P'+idxPlate+'_'+idxWell+"_NUC.tif",imgDet, [cv2.CV_16U])


## Function to detect lipid droplets in BODIPY channel ##
## Input: imgToProces: Image of LD
## typicDiam: Average diameter of lipid droplets /!\ Depend on the resolution
## Output: Image with cleaned LD
def detectGFP(imgToProces, typicDiam): 
	halfDiam= round(typicDiam/2)
	
	# Clean the background
	imgToProcNorm= norm255(imgToProces)
	imgToProcNorm=np.uint8(imgToProcNorm)   
	ret,imgToProcOtsu = cv2.threshold(imgToProcNorm,0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)   
	imgToProcOtsu[imgToProcOtsu>0]=1
	cleanImage=imgToProcOtsu*imgToProces
 
	# Enhance the edges
	sigma = 0.6
	img = cleanImage.max()-cleanImage    
	eigenimage = enhanceContrast(img, sigma)
	eigenimage = eigenimage.max()-eigenimage
	eigenimage[eigenimage<100]=0
	eigenimage[eigenimage>0]=1  
	eigenimage= norm255(cleanImage)-norm255(eigenimage)
	eigenimage[eigenimage<=0]=0  
	
	maskEigen=eigenimage.copy()
	maskEigen[maskEigen>0]=1  
	

	# Enhance the droplets         
	TH1=TopHatBasic(cleanImage, typicDiam- halfDiam)
	TH2=TopHatBasic(cleanImage, typicDiam )   
	TH3=TopHatBasic(cleanImage, typicDiam+ halfDiam) 
	THimage= TH1 +TH2 +TH3
	THimage[THimage>0]=1 
	
	maskTH=THimage.copy()
	maskTH[maskTH>0]=1 


	## Compute mask of the 2 enhancements; 1 values are the uncertainty
	sumMask=maskTH+maskEigen



	## We cleaned the uncertainty based on the intensity
	T= sumMask.copy()
	T[T==2]=0
	T= norm255(T)
	T=np.uint8(T)  
	Topen=ndimage.binary_opening(T, structure=np.ones((2,2))).astype(np.int) 
	Topen=Topen*imgToProces
	Topen= norm255(Topen)
	Topen=np.uint8(Topen) 
	H=(Topen[Topen>0])
	ret,imgToProcOtsu = cv2.threshold(H,0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU) 
	TOtsu=Topen.copy()
	TOtsu[TOtsu<ret]=0
	TOtsu[TOtsu!=0]=1
	
	finalMask=sumMask.copy()
	finalMask[finalMask<2]=0
	finalMask[finalMask!=0]=1
	finalMask=finalMask+TOtsu

	return finalMask

## Function to detect nuclei in Hoescht channel ##
## Input: imgToProces: Image of nuclei
## typicDiam: Average diameter of nuclei /!\ Depend on the resolution
## Output: Image with cleaned nuclei
def detectNucs(imgToProces, typicDiam):
	halfDiam= round(typicDiam/2)
	imgProc=TopHatBasic(imgToProces, typicDiam- halfDiam) + TopHatBasic(imgToProces, typicDiam) + TopHatBasic(imgToProces, typicDiam+ halfDiam)
	imgProcNorm= norm255(imgProc)
	imgProcNorm=np.uint8(imgProcNorm)
	imgProcNorm[imgProcNorm>0]=1
	imgProcOpen=ndimage.binary_opening(imgProcNorm, structure=np.ones((4,4))).astype(np.int)
	
	return imgProcOpen

## Function of normalization [0;255] ##
## Input: Image to normalize
## Output: Normalized image
def norm255(ImgIn):
	ImgIn= ImgIn.astype(float)
	imgOut=((ImgIn-ImgIn.min())/(ImgIn.max()-ImgIn.min()))*255

	return imgOut

## Modified New White Top Hat ##
## Enhanced dim bright spot and reduce local background
## Input: f : Image to process
##		  diamVal: Typical diameter for the filtering
## Output: Processed image
def TopHatBasic(f, diamVal):       
	nL=int(diamVal)
	nW=int(diamVal)
	nM=int(round(nW/7))

	Bb = np.ones((nL,nL))
	dB= np.ones((nW,nW))
	dB[nM:nW-nM,nM:nW-nM]=0;
  
	# White
	dilatedI = ndimage.grey_dilation(f, footprint=dB)
	blackSquare= ndimage.grey_erosion (dilatedI, footprint=Bb)
	NCoi= np.minimum (blackSquare, f)
	MNWTH= f- NCoi

	return MNWTH


## Enhanced contrast curvature ##
## Enhanced edges and gap between LD
## Input: forimg : Image to process
##		  sigma: Sigma for the gaussian kernel
## Output: Processed image
def enhanceContrast(orimg,sigma):
	Img = np.copy(orimg)
	gau = gauss_kern2(Img, sigma)
	Imgfft = np.fft.rfft2(Img)
	gfft = np.fft.rfft2(gau)
	fftimage = np.multiply(Imgfft, gfft)
	Img_smooth =np.real(np.fft.ifftshift( np.fft.irfft2(fftimage)))

	Iy, Ix = np.gradient(Img_smooth)
	Ixy, Ixx = np.gradient(Ix)
	Iyy, Iyx = np.gradient(Iy)
	
	eigvallam = np.zeros((2,Ixx.shape[0],Ixx.shape[1] ))
	trhessian = Ixx+Iyy
	dethessian = Ixx*Iyy-Ixy*Ixy
	
	eigvallam[0,:,:] = 0.5*(trhessian + np.sqrt(trhessian*trhessian - (4*dethessian) ))
	eigvallam[1,:,:] = 0.5*(trhessian - np.sqrt(trhessian*trhessian - (4*dethessian) ))
	eigh1 = eigvallam.min(0)
	eigh2 = eigvallam.max(0)
	
	etemp = np.copy(eigh1)
	etemp = etemp-etemp.min()
	etemp = etemp/etemp.max()
	etemp = np.uint8(etemp*255)

	return(etemp)   


def gauss_kern2(Img, sigma):
	""" Returns a normalized 2D gauss kernel array for convolutions """
	h2,h1 = Img.shape    
	x, y = np.mgrid[0:h2, 0:h1]
	x = x-h2/2
	y = y-h1/2

	g = np.exp( -( x**2 + y**2 ) / (2*sigma**2) );
	return g / g.sum()   




#################### Main ####################

def segmentFatDroplet(listPlates):
## Detection !

	## Initialise variables
	global q
	global inputCellProfilerPath
	global outputDetPath
	global list4CP

	## Initialise paths	
	inputDataPath=settings.pathList[1]
	outputDetPath=settings.pathList[3]
	inputCellProfilerPath=settings.pathList[4]

	## Initialise threads
	nProc = multiprocessing.cpu_count()
	threadsList = [threading.Thread(target=workerSeg) for i in range(0, nProc)]
	q = Queue.Queue()
	
	## Fill up the queue for threading
	for currFolder in listPlates:
		currPath=inputDataPath+str(currFolder)+"/"
		idxSpace=currFolder.rsplit('plate')
		currFolder=idxSpace[1]
		  
		listSubFolder= [x for x in os.listdir(currPath) if os.path.isdir(currPath+x)]
		for currSubFolder in listSubFolder:
			currSubPath=currPath+ str(currSubFolder)+"/"
			q.put(currSubPath)

	# start all threads        
	for thread in threadsList:
		thread.start()
	# join all threads
	for thread in threadsList:
		thread.join()

	## Get the paths of the processed images
	list4CP=[]
	listImages= [x for x in os.listdir(inputCellProfilerPath) if x.endswith('_GFP.tif')]
	listSubImages=[[x.rsplit('_')[0]+'_'+x.rsplit('_')[1]+'_GFP.tif'] + [x.rsplit('_')[0]+'_'+x.rsplit('_')[1]+'_NUC.tif'] for x in listImages]
	listSubImages= [item for sublist in listSubImages for item in sublist]

	## Export paths for input in CellProfiler
	list4CP = [[os.path.abspath(inputCellProfilerPath+x)] for x in listSubImages] 
	list4CP=sorted(list4CP, key=lambda x: x[0])
	with open(inputCellProfilerPath + "CellProfilerInput.csv",'wb') as csvfile:
		writer= csv.writer(csvfile, dialect='unixpwd')	
		for i in range(0, len(list4CP)):
			content=list4CP[i]
			writer.writerow(content) 	