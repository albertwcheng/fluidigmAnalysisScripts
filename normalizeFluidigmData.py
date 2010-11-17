#!/usr/bin/env python

'''

Copyright 2010 Wu Albert Cheng <albertwcheng@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

Usage: For normalizing fluidigm output to endogenous controls

'''

from optparse import OptionParser
from albertcommon import *
from sys import *

''' 
tab file usual format

[:::::			R 1			:::::]
Index			Excel			Field
-----			-----			-----
1			A			ChamberID
2			B			SampleName
3			C			SampleType
4			D			SamplerConc
5			E			FAMMGBName
6			F			FAMMGBType
7			G			CtValue
8			H			CtQuality
9			I			CtCall
10			J			CtThreshold

'''

def emptyAllTheseArrays(listOfArrays):
	for L in listOfArrays:
		del L[:]


def readFluidigmFile(filename,SampleNameCol,geneNameCol,CtValueCol,CtQualityCol,CtCallCol,headerRow,startRow,fs,geneNameMap):
	header,prestarts=getHeader(filename,headerRow,startRow,fs)
	SampleNameCol=getCol0ListFromCol1ListStringAdv(header,SampleNameCol)[0]
	GeneNameCol=getCol0ListFromCol1ListStringAdv(header,GeneNameCol)[0]
	CtValueCol=getCol0ListFromCol1ListStringAdv(header,CtValueCol)[0]
	CtQualityCol=getCol0ListFromCol1ListStringAdv(header,CtQualityCol)[0]
	CtCallCol=getCol0ListFromCol1ListStringAdv(header,CtCallCol)[0]
	
	#return (sampleNames,sampleStruct)
	# sampleNames = [names...]
	# sampleStruct = [CtValues,CtQuals,CtCalls]
	#  CtValues=[...]
	#  CtQuals=[...]
	#  CtCalls=[...]
	
	sampleNames=[] #such that the order of sample appearance is retained

	
	sampleNameToStructs=dict()
	
	#now load the values
	#open file first
	fil=open(filename)
	for lin in fil:
		fields=lin.rstrip("\r\n").split(fs)
		SampleName=fields[SampleNameCol]
		GeneName=fields[GeneNameCol]
		CtValue=float(fields[CtValueCol])
		CtQuality=float(fields[CtQualityCol])
		CtCall=fields[CtCallCol]
		
		
		if geneNameMap and len(geneNameMap)>0: #if geneName Map is present
			try:
				GeneName=geneNameMap[GeneName] #map it to proper name
			except:
				pass
		
		#now find the sample struct for that sample name
		try:
			sampleStruct=sampleNameToStructs[SampleName]
		except KeyError:
			sampleStruct=dict()
			sampleNameToStructs[SampleName]=sampleStruct
			sampleNames.append(SampleName) #new sample add it
		
		#now find the gene
		try:
			CtValuesForThisGene,CtQualsForThisGene,CtCallsForThisGene=sampleStruct[GeneName]
		except:
			CtValuesForThisGene=[]
			CtQualsForThisGene=[]
			CtCallsForThisGene=[]
			dataForThisGene=[CtValuesForThisGene,CtQualsForThisGene,CtCallsForThisGene]
			sampleStruct[GeneName]=dataForThisGene
			
		CtValuesForThisGene.append(CtValue)
		CtQualsForThisGene.append(CtQuality)
		CtCallsForThisGene.append(CtCall)
		
	fil.close()
	
	sampleStructs=[]
	
	for sampleName in sampleNames:  #such that the order of sample appearance is retained
		sampleStructs.append(sampleNameToStructs[sampleName])
	
	return sampleNames,sampleStructs


def deleteAnElementFromAllArrays(listOfArrays,idx):
	for arr in listOfArrays:
		del arr[idx]

def filterFluidigmData(sampleNames,sampleStructs,options): #inplace, return a list of invalid samples [to be ignored for normalization and NA'ed for output
	'''
	parser.add_option("--keep-Ct-call-failed",dest="discardCtCallFailed",default=True,action="store_false",help="keep Ct Call failed entries  (applied to data and controls) [False]")
	parser.add_option("--Ct-quality-threshold",dest="CtQualityThreshold",default=-1.0,type=float,help="specify the CtQualityThreshold under which the value will be discarded (applied to data and controls) [-1.0:No treshold]")
	parser.add_option("--Ct-value-threshold",dest="CtValueThreshold",default=50.0,type=float,help="specify the max Ct value to be valid (applied to data and controls) [50.0]")
	parser.add_option("--max-Ct-deviation-between-replicates",dest="MaxCtRepDev",default=5.0,type=float,help="set the maximal Ct deviation between replicates allowed  (applied to data and controls) [5.0]")
	parser.add_option("--min-number-of-valid-data-point-per-control",dest="minValidReplicatesControl",default=1,type=int,help="set a minimum number of data points for that control [1]")
	parser.add_option("--Ct-value-threshold-for-per-control-average",dest="CtValueThresholdPerControl",default=50.0,type=float,help="set a Ct value threshold per control (average among replicates of that control) higher than which the whole row of data will be discarded [50.0]")
	parser.add_option("--min-number-of-valid-data-point-per-gene",dest="minValidReplicates",default=2,type=int,help="set a minimum number of data points for that gene (not including controls) in that sample for calling normalized value [2]")	
	parser.add_option("--Ct-value-threshold-for-data-average",dest="CtValueThresholdPerData",default=50.0,type=float,help="set a Ct value threshold of the data (average among replicates of that data before normalization to controls) higher than which the data will be discarded [50.0]")


	options.controlNames

	'''
	
	sampleToRemove=[]
	
	for idx in range(len(sampleNames),-1,-1):
		sampleName=sampleNames[idx]
		sampleStruct=sampleStructs[idx]

		invalidSample=False

		#first apply Ct quality and Ct Call threshold in all data including controls
		for geneName,dataForThisGene in sampleStruct.items():
			CtValuesForThisGene,CtQualsForThisGene,CtCallForThisGene=dataForThisGene
			numValuesForThisGene=len(CtValuesForThisGene)
			
			isControl=(geneName in options.controlNames)
			
			for j in range(numValuesForThisGene,-1,-1): #from back for convenient deletion
				thisCtValue=CtValuesForThisGene[j]
				thisCtQuality=CtQualsForThisGene[j]
				thisCtCall=CtCallForThisGene[j]
				if options.discardCtCallFailed && thisCtCall.strip().upper()=="FAIL":
					#discard this
					deleteAnElementFromAllArrays(dataForThisGene,j)
					continue
					
				if thisCtQuality<options.CtQualityThreshold):
					deleteAnElementFromAllArrays(dataForThisGene,j)
					continue
				
				if thisCtValue>=options.CtValueThreshold:
					deleteAnElementFromAllArrays(dataForThisGene,j)
					continue
					
			#now filter for deviation
			dev=max(CtValuesForThisGene)-min(CtValuesForThisGene)
			if dev>options.MaxCtRepDev:
				#remove all elements!!
				emptyAllTheseArrays(dataForThisGene)
					
				
		#now for each control genes:
		for controlName in options.controlNames:
			try:
				CtValuesForThisGene,CtQualsForThisGene,CtCallForThisGene=sampleStruct[controlName]
			except KeyError:
				#control not found, remove the whole row
				sampleToRemove.append(idx)
				invalidSample=True
				break #don't care about other stuff
				
			if len(CtValuesForThisGene)==0 or len(CtValuesForThisGene)<options.minValidReplicatesControl:
				#this whole row is to be removed
				sampleToRemove.append(idx)
				invalidSample=True
				break #don't care about other stuff
			
			
			controlCtMean=mean(CtValuesForThisGene)
			if controlCtMean>options.CtValueThresholdPerControl:
				sampleToRemove.append(idx)
				invalidSample=True
				break #don't care about other stuff
				
		#before checking all data points, ignore if sample is already to be removed
		if not invalidSample:
			#now check each points
			for geneName,dataForThisGene in sampleStruct.items():
				CtValuesForThisGene,CtQualsForThisGene,CtCallForThisGene=dataForThisGene
				numValuesForThisGene=len(CtValuesForThisGene)
				
				if numValuesForThisGene<options.minValidReplicates:
					#remove this particular gene
					emptyAllTheseArrays(dataForThisGene)				
					
				isControl=(geneName in options.controlNames)
				
				if isControl:
					continue #ignore controls
				
				dataMean=mean(CtValuesForThisGene)
				if dataMean>options.CtValueThresholdPerData:
					emptyAllTheseArrays(dataForThisGene)	
	
	
	return sampleToRemove			
		
def normalizeFluidigmData(sampleNames,sampleStructs,options,invalidSampleRows):		
	for idx in range(0,len(sampleNames)):
		if idx in invalidSampleRows:
			continue #no need to process these
		
		sampleName=sampleNames[idx]
		sampleStruct=sampleStructs[idx]
		
		#now get control values:
		controlValues=[]
		for controlName in options.controlNames:
			controlValues.append(mean(sampleStruct[controlName][0])) #column 0 is CtValues
		
		controlMean=mean(controlValues)
	
		#now go to data and normalize by data-controlMean
		for geneName,dataForThisGene in sampleStructs.items():
			CtValuesForThisGene,CtQualsForThisGene,CtCallForThisGene=dataForThisGene
			#not collapsed yet
			for j in range(0,len(CtValuesForThisGene)):
				CtValuesForThisGene[j]-=controlMean

def getMeanCtValueForGene(sampleStruct,geneName):
	return mean(sampleStruct[geneName][0]) #0 stores CtValues
	
def printFluidigmData(sampleNames,sampleStructs,options,invalidSampleRows,geneNamesInOrder,printControl):
	#now print header
	fieldsToPrint=["SampleName"]
	if printControl:
		fieldsToPrint.extend(options.controlNames)
	
	geneNamesInOrderNoControls=[]
	for geneName in geneNamesInOrder:
		if geneName not in options.controlNames:
			geneNamesInOrderNoControls.append(geneName)
	
	fieldsToPrint.extend(geneNamesInOrderNoControls)
	
	numColsToPrint=len(fieldsToPrint)
	
	print >> stdout,"\t".join(fieldsToPrint)
	
	
	
	#now print content
	for idx in range(0,len(sampleNames)):
		sampleName=sampleNames[idx]
		sampleStruct=sampleStructs[idx]
		
		fieldsToPrint=[sampleName]
		if idx in invalidSampleRows: #print all NAs
			fieldsToPrint.extend(["NA"]*(numColsToPrint-1))
		else:
			#valid!
			#print controls?
			if options.printControl:
				for controlName in options.controlNames:
					fieldsToPrint.append(getMeanCtValueForGene(controlName))
			
			#now print data
			for geneName in geneNamesInOrderNoControls:
				fieldsToPrint.append(getMeanCtValueForGene(geneName))
			
	print >> stdout,"\t".join(fieldsToPrint)
	
		
if __name__=='__main__':
	parser=OptionParser("usage: %prog [options] infile > outfile")
	parser.add_option("--control-names",dest="controlNames",default="gapdh,hprt",help="specify control gene names separate by comma [gapdh,hprt]")

	#parser.add_option("--not-collapse-replicate",dest="collapseReplicates",default=True,action="store_false",help="do not collapse replicates")
	parser.add_option("--gene-name-map",dest="geneNameMap",default=None,help="specify mapping of gene ID to gene Name [None] formatted: [geneID]tab[geneName]")
	
	#filters:
	parser.add_option("--keep-Ct-call-failed",dest="discardCtCallFailed",default=True,action="store_false",help="keep Ct Call failed entries  (applied to data and controls) [False]")
	parser.add_option("--Ct-quality-threshold",dest="CtQualityThreshold",default=-1.0,type=float,help="specify the CtQualityThreshold under which the value will be discarded (applied to data and controls) [-1.0:No treshold]")
	parser.add_option("--Ct-value-threshold",dest="CtValueThreshold",default=50.0,type=float,help="specify the max Ct value to be valid (applied to data and controls) [50.0]")
	parser.add_option("--max-Ct-deviation-between-replicates",dest="MaxCtRepDev",default=5.0,type=float,help="set the maximal Ct deviation between replicates allowed  (applied to data and controls) [5.0]")
	parser.add_option("--min-number-of-valid-data-point-per-control",dest="minValidReplicatesControl",default=1,type=int,help="set a minimum number of data points for that control [1]")
	parser.add_option("--Ct-value-threshold-for-per-control-average",dest="CtValueThresholdPerControl",default=50.0,type=float,help="set a Ct value threshold per control (average among replicates of that control) higher than which the whole row of data will be discarded [50.0]")
	parser.add_option("--min-number-of-valid-data-point-per-gene",dest="minValidReplicates",default=2,type=int,help="set a minimum number of data points for that gene (not including controls) in that sample for calling normalized value [2]")	
	parser.add_option("--Ct-value-threshold-for-data-average",dest="CtValueThresholdPerData",default=50.0,type=float,help="set a Ct value threshold of the data (average among replicates of that data before normalization to controls) higher than which the data will be discarded [50.0]")
	
			
	#formatting modifyers
	parser.add_option("--fs",dest="fs",default="\t",help="specify the field separator of the infile [tab]")
	parser.add_option("--start-row",dest="startRow",default=2,type=int,help="specify the start row of data (excluding header) [2]")
	parser.add_option("--header-row",dest="headerRow",default=1,type=int,help="specify the header row of data [1]")
	parser.add_option("--sample-name-col",dest="sampleNameCol",default="2",help="specify the column of the infile that contains the sample name [2]")
	parser.add_option("--gene-name-col",dest="geneNameCol",default="5",help="specify the column of the infile that contains the gene name [5]")	
	parser.add_option("--Ct-value-col",dest="CtValueCol",default="7",help="specify the column of the infile that contains the Ct value [7]")
	parser.add_option("--Ct-quality-col",dest="CtQualityCol",default="8",help="specify the column of the infile that contains the Ct Call  [8]")
	parser.add_option("--Ct-call-col",dest="CtCallCol",default="9",help="specify the column of the infile that contains the Ct Call  [9]")
	
	parser.add_option("--print-controls",dest="printControl",default=False,action="store_true",help="print also raw control Ct values [NO]")
		
	(options, args) = parser.parse_args(argv)
	
	try:
		infile,=args
	except:
		parser.print_help()
		exit()
	
	geneNamesInOrder=[]
	
	geneNameMap=dict()

	if parser.geneNameMap:
		#if a gene name map is specified, usu is, then load it
		print >> stderr,"load gene map from",parser.geneNameMap
		fil=open(parser.geneNameMap)
		for lin in fil:
			fields=lin.strip().split("\t")
			geneNameMap[fields[0]]=fields[1]
			geneNamesInOrder.append(fields[1])
			
		fil.close()	
	
	#now process arguments
	parser.controlNames=parser.controlNames.split(",")
	
	#now load the infile
	sampleNames,sampleStructs=readFluidigmFile(infile,parser.SampleNameCol,parser.geneNameCol,parser.CtValueCol,parser.CtQualityCol,parser.CtCallCol,parser.headerRow,parser.startRow,parser.fs,geneNameMap)
	
	#filter data
	invalidSampleRows=filterFluidigmData(sampleNames,sampleStructs,options)
	
	#now normalize data
	normalizeFluidigmData(sampleNames,sampleStructs,options,invalidSampleRows)
	
	#now output
	printFluidigmData(sampleNames,sampleStructs,options,invalidSampleRows,geneNamesInOrder,printControl)
	