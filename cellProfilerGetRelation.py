# cellProfilerGetRelation.py
#######################################################
#
#	Relationship pipeline for :
#	- Load path of images
#	- Separate image in thread and run CellProfiler
#	- The current CellProfiler creates a relationship between the location of
#	the nuclei and the lipid droplets i.e. cells approximation, extract measurements and create a database
#	for CellProfilerAnalyst
#
#######################################################

import os
import sys

import settings

import multiprocessing
import Queue
import threading
import subprocess
import shlex
import shutil
import csv

import math


## Thread
## Get the paths of the images and the CellProfiler pipeline 
## Run CellProfiler
## Check the output and the completion of the process
def workerCP():
	while not q.empty():
		idxFile = q.get()

		projPath=(os.path.abspath("./FatPipelineV3.cppipe")).replace('\\', '/')

		doneTempPath=os.path.abspath(inputCellProfilerPath + "temp/done_"+ str(idxFile) + ".txt")
		doneTempPath=doneTempPath.replace('\\', '/')
		csvTempPath=os.path.abspath(inputCellProfilerPath + "/temp/CellProfilerInput_"+ str(idxFile) +".csv")
		csvTempPath=csvTempPath.replace('\\', '/')
		outTempPath= os.path.abspath(outputCellProfilerPath+ "/temp/Analysis_"+ str(idxFile) +"/")
		outTempPath=outTempPath.replace('\\', '/')

		command=  CPPath +" --jvm-heap-size=4g -c -r -o " + outTempPath + " -p "+ projPath + " --file-list=" + csvTempPath +" -d " + doneTempPath


		print ('command:%s\n'%command)
		command=shlex.split(command)		

		p=subprocess.Popen(command, stdout=Fout, stderr=Ferr)

		while not os.path.isfile(doneTempPath):
			pass
		while (os.stat(doneTempPath).st_size==0):
			pass
		with open(doneTempPath, 'r') as fOutput:
			first_line = fOutput.readline()
		if first_line.rstrip() != 'Complete':
			print('Fail')
			os.remove(inputCellProfilerPath + "temp/done_"+ str(idxFile) + ".txt")
			q.put(idxFile)
		else:
			if p:
				p.kill()




#################### Main ####################

def runCellProfilerRelationship():
## CellProfiler !

	## Initialise variables
	global q 
	global inputCellProfilerPath
	global outputCellProfilerPath
	global CPPath

	global Ferr
	global Fout

	if not os.path.exists('./../logFiles/'):
		os.makedirs('./../logFiles/')
	Ferr=open('./../logFiles/Ferr.txt','w')
	Fout=open('./../logFiles/Fout.txt','w')

	## Initialise paths
	CPPath=settings.pathList[0]
	inputCellProfilerPath=settings.pathList[4]
	outputCellProfilerPath=settings.pathList[5]

	outPath= os.path.abspath(outputCellProfilerPath)
	outPath=outPath.replace('\\', '/')
	csvPath=  os.path.abspath(inputCellProfilerPath + "CellProfilerInput.csv")
	csvPath=csvPath.replace('\\', '/')

	## Initialise threads
	nProc = multiprocessing.cpu_count()
	q = Queue.Queue()

	## Load mask image paths
	with open(inputCellProfilerPath + "CellProfilerInput.csv", 'r') as mainCsv:
		reader = csv.reader(mainCsv)
		readerData = list(reader)
	
	## Define the number of images to process in parallel
	num_seqImg=len(readerData)
	if (num_seqImg%2 != 0 or num_seqImg==0):
		print 'List broken, run again'
		sys.exit(1)
	if nProc>num_seqImg:
		nProc=num_seqImg

	threadsList = [threading.Thread(target=workerCP) for i in range(0, nProc)]
	num_seqImg=len(readerData)/2

	stepImg=50
	nProc=int(round(float(num_seqImg*2)/stepImg))

	## Create temporary folder for thread outputs
	for i in range(0, nProc):
		processFlag=True
		if not os.path.exists(outputCellProfilerPath + "temp/"):
			os.makedirs(outputCellProfilerPath + "temp/")
		if not os.path.exists(inputCellProfilerPath + "temp/"):
			os.makedirs(inputCellProfilerPath + "temp/")
		else:
			doneTemp=os.path.abspath(inputCellProfilerPath + "temp/done_"+ str(i) + ".txt")	
			if os.path.exists(doneTemp):
				with open(doneTemp, 'r') as fOutput:
					first_line = fOutput.readline()
				if first_line.rstrip() != 'Complete':
					os.remove(inputCellProfilerPath + "temp/done_"+ str(i) + ".txt")
				else:
					processFlag=False
					print first_line.rstrip()	

		## Fill up the queue with paths
		if processFlag:	
			with open(inputCellProfilerPath + "temp/CellProfilerInput_"+ str(i) +".csv",'wb') as csvfile:
				writer= csv.writer(csvfile)	
				content=readerData[i*stepImg: min(((i+1)*stepImg), (num_seqImg*2))]
				writer.writerows(content)
				if (content):
					q.put(i)


	# start all threads        
	for thread in threadsList:
		thread.start()
			
	# join all threads
	for thread in threadsList:
		thread.join()

	## Move all tif out of /temp
	for idxProc in range(0, nProc):
		currSubPath=outputCellProfilerPath+ "temp/Analysis_"+ str(idxProc) +"/"
		if os.path.isdir(currSubPath):
			listFiles= [x for x in os.listdir(currSubPath) if x.endswith('.tif')]
			[shutil.move(currSubPath+x,outputCellProfilerPath) for x in listFiles]