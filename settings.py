# settings.py
#######################################################
#
#	Definition of the different paths:
#	- CellProfiler (Software, input, output)
#	- Input
#	- Output
#
#######################################################
import os

def init():
	global pathList

	CPPath= "D:/Logiciel/CellProfiler2.2/CellProfiler.exe"
	inputDataPath= "C:/Data/Granulometry/Data/"

	resultPath= "./../Results/"
	colorDisplayPath =resultPath +"colorDisplay/"
	outputDetPath = resultPath + "outputResults/"
	inputCellProfilerPath =resultPath +"inputCP/"
	outputCellProfilerPath =resultPath +"CPResults/"
	if not os.path.isdir(resultPath):
		os.mkdir(resultPath)
	if not os.path.isdir(colorDisplayPath):
		os.mkdir(colorDisplayPath)
	if not os.path.isdir(outputDetPath):
		os.mkdir(outputDetPath)
	if not os.path.isdir(inputCellProfilerPath):
		os.mkdir(inputCellProfilerPath)
	if not os.path.isdir(outputCellProfilerPath):
		os.mkdir(outputCellProfilerPath)	

	pathList = [CPPath]+ [inputDataPath] + [resultPath] + [outputDetPath] + [inputCellProfilerPath] + [outputCellProfilerPath]
