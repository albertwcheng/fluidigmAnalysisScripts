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
def mean(numberList):
    floatNums = [float(x) for x in numberList]
    return sum(floatNums) / len(numberList)
    
def emptyAllTheseArrays(listOfArrays):
	for L in listOfArrays:
		del L[:]


def readFluidigmFile(filename,SampleNameCol,geneNameCol,CtValueCol,CtQualityCol,CtCallCol,headerRow,startRow,fs,geneNameMap):
	header,prestarts=getHeader(filename,headerRow,startRow,fs)
	SampleNameCol=getCol0ListFromCol1ListStringAdv(header,SampleNameCol)[0]
	GeneNameCol=getCol0ListFromCol1ListStringAdv(header,geneNameCol)[0]
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
	lino=0
	for lin in fil:
		lino+=1
		if lino<startRow:
			continue
		fields=lin.rstrip("\r\n").split(fs)
		SampleName=fields[SampleNameCol]
		GeneName=fields[GeneNameCol]
		CtValue=float(fields[CtValueCol])
		CtQuality=float(fields[CtQualityCol])
		CtCall=fields[CtCallCol]
		
		
		if geneNameMap and len(geneNameMap)>0: #if geneName Map is present
			try:
				#print >> stderr,GeneName,"to",
				GeneName=geneNameMap[GeneName] #map it to proper name
				#print >> stderr,GeneName
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

def filterFluidigmData(sampleNames,sampleStructs,options,geneMapID): #inplace, return a list of invalid samples [to be ignored for normalization and NA'ed for output
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
	#filterMessages=[] #construct similar matrix structure as sampleStructsToRememberFilterStatus
	
	for idx in range(len(sampleNames)-1,-1,-1):
		sampleName=sampleNames[idx]
		sampleStruct=sampleStructs[idx]

		invalidSample=False

		#first apply Ct quality and Ct Call threshold in all data including controls
		for geneName,dataForThisGene in sampleStruct.items():
			CtValuesForThisGene,CtQualsForThisGene,CtCallForThisGene=dataForThisGene
			numValuesForThisGene=len(CtValuesForThisGene)
			
			isControl=(geneName in options.controlNames)
			
			for j in range(numValuesForThisGene-1,-1,-1): #from back for convenient deletion
				thisCtValue=CtValuesForThisGene[j]
				thisCtQuality=CtQualsForThisGene[j]
				thisCtCall=CtCallForThisGene[j]
				if options.discardCtCallFailed and thisCtCall.strip().upper()=="FAIL":
					#discard this
					print >> stderr,"Filter1: Discard data point due to FAIL Ct Call","for sample",sampleName,"and gene",geneName,"ID:",geneMapID[geneName]
					deleteAnElementFromAllArrays(dataForThisGene,j)
					continue
					
				if thisCtQuality<options.CtQualityThreshold:
					print >> stderr,"Filter2: Discard data point due to CtQuality lower than Threshold. CtQuality=",thisCtQuality,"<",options.CtQualityThreshold,"=threshold","for sample",sampleName,"and gene",geneName,"ID:",geneMapID[geneName]
					deleteAnElementFromAllArrays(dataForThisGene,j)
					continue
				
				if thisCtValue>=options.CtValueThreshold:
					print >> stderr,"Filter3: Discard data point due to CtValue >= Threshold. CtValue=",thisCtValue,">=",options.CtValueThreshold,"=threshold","for sample",sampleName,"and gene",geneName,"ID:",geneMapID[geneName]
					deleteAnElementFromAllArrays(dataForThisGene,j)
					continue
			
			if len(CtValuesForThisGene)>0:		
				#now filter for deviation
				dev=max(CtValuesForThisGene)-min(CtValuesForThisGene)
				if dev>options.MaxCtRepDev:
					#remove all elements!!
					print >> stderr,"Filter4: Discard data point due to deviation > Maximal. dev=",dev,">=",options.MaxCtRepDev,"=threshold","for sample",sampleName,"and gene",geneName,"ID:",geneMapID[geneName]
					#emptyAllTheseArrays(dataForThisGene)
					for k in range(0,len(CtCallForThisGene)):
						CtCallForThisGene[k]="INC" #inconsistent, mark that so that we know later not to 0 this gene
					
					
				
		#now for each control genes:
		for controlName in options.controlNames:
			try:
				CtValuesForThisGene,CtQualsForThisGene,CtCallForThisGene=sampleStruct[controlName]
			except KeyError:
				#control not found, remove the whole row
				print >> stderr,"Filter5: Discard data row (the whole sample) because control",controlName,"not found for sample",sampleName
				sampleToRemove.append(idx)
				invalidSample=True
				break #don't care about other stuff
				
			#added Dec 4th.
			if CtCallForThisGene[0]=="INC":
				#inconsistent replicates of control
				print >> stderr,"Filter10: Discard data row (the whole sample) because max-min of control",controlName,"was >maximum=",options.MaxCtRepDev,"for sample",sampleName
				sampleToRemove.append(idx)
				invalidSample=True
				break
				
			if len(CtValuesForThisGene)==0 or len(CtValuesForThisGene)<options.minValidReplicatesControl:
				#this whole row is to be removed
				print >> stderr,"Filter6: Discard data row (the whole sample) because number of replicates for",controlName,"was",len(CtValuesForThisGene),"<minimum=",options.minValidReplicatesControl,"for sample",sampleName
				sampleToRemove.append(idx)
				invalidSample=True
				break #don't care about other stuff
			
			
			controlCtMean=mean(CtValuesForThisGene)
			if controlCtMean>options.CtValueThresholdPerControl:
				print >> stderr,"Filter7: Discard data row (the whole sample) because mean control Ct for",controlName,"was",controlCtMean,">maximum=",options.CtValueThresholdPerControl,"for sample",sampleName
				sampleToRemove.append(idx)
				invalidSample=True
				break #don't care about other stuff
		
			
		#before checking all data points, ignore if sample is already to be removed
		if not invalidSample:
			#now check each points
			for geneName,dataForThisGene in sampleStruct.items():
			
				isControl=(geneName in options.controlNames)
				
				if isControl:
					continue #ignore controls
								
				CtValuesForThisGene,CtQualsForThisGene,CtCallForThisGene=dataForThisGene
				numValuesForThisGene=len(CtValuesForThisGene)
				
				if CtCallForThisGene[0]=="INC":
					#inconsistent replicates of gene : don't do anything
					continue
				
				if numValuesForThisGene<options.minValidReplicates:
					#remove this particular gene
					print >> stderr,"Filter8: Discard gene because number of replicates for",geneName,"was",numValuesForThisGene,"<minimum=",options.minValidReplicates,"for sample",sampleName,"and gene",geneName,"ID:",geneMapID[geneName]
					emptyAllTheseArrays(dataForThisGene)				
					

				
				if len(CtValuesForThisGene)>0:
					dataMean=mean(CtValuesForThisGene)
					if dataMean>options.CtValueThresholdPerData:
						print >> stderr,"Filter9: Discard gene because mean Ct for",geneName,"was",dataMean,">maximum=",options.CtValueThresholdPerData,"for sample",sampleName
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
			#print >> stderr,idx,sampleName,sampleStruct[controlName]
			controlValues.append(mean(sampleStruct[controlName][0])) #column 0 is CtValues
		
		controlMean=mean(controlValues)
	
		#now go to data and normalize by data-controlMean
		for geneName,dataForThisGene in sampleStruct.items():
		
			if geneName in options.controlNames:
				continue #not normalize the controls
		
			CtValuesForThisGene,CtQualsForThisGene,CtCallForThisGene=dataForThisGene
			#not collapsed yet
			for j in range(0,len(CtValuesForThisGene)):
				CtValuesForThisGene[j]-=controlMean 

def getMeanCtValueForGene(sampleStruct,geneName,NAString):
	try:
		return mean(sampleStruct[geneName][0]) #0 stores CtValues
	except:
		return "NA"

def getStrArray(L):
	S=[]
	for x in L:
		S.append(str(x))
	return S
		
def printFluidigmData(sampleNames,sampleStructs,options,invalidSampleRows,geneNamesInOrder,printControl,geneMapID):
	#now print header
	
	outputMatrix=[]

	if options.useConventionalDeltaCt:
		blankValue=options.missingDataPolicy[2]
	else:
		blankValue=options.missingDataPolicy[1]	
	
	fieldsToPrint=["SampleName"]
	if printControl:
		#fieldsToPrint.extend(options.controlNames)
		for controlName in options.controlNames:
			fieldsToPrint.append(options.controlPrefix+controlName)
	
	geneNamesInOrderNoControls=[]
	for geneName in geneNamesInOrder:
		if geneName not in options.controlNames and geneName not in geneNamesInOrderNoControls: #no duplicate names and no control names
			geneNamesInOrderNoControls.append(geneName)
	
	fieldsToPrint.extend(geneNamesInOrderNoControls)
	
	numColsToPrint=len(fieldsToPrint)
	
	#print >> stdout,"\t".join(fieldsToPrint)
	outputMatrix.append(fieldsToPrint)
	

	#now print content
	for idx in range(0,len(sampleNames)):
		sampleName=sampleNames[idx]
		sampleStruct=sampleStructs[idx]
		
		fieldsToPrint=[sampleName]
		if idx in invalidSampleRows: #print all NAs
			fieldsToPrint.extend([options.NAString]*(numColsToPrint-1))
		else:
			#valid!
			#print controls?
			if options.printControl:
				for controlName in options.controlNames:
					fieldsToPrint.append(str(getMeanCtValueForGene(sampleStruct,controlName,options.NAString)))
			
			#now print data
			for geneName in geneNamesInOrderNoControls:
				try:
					CtValuesForThisGene,CtQualsForThisGene,CtCallForThisGene=sampleStruct[geneName]
				except:
					#gene not found
					print >> stderr,"warning: gene data not found for sample",sampleName,"and gene",geneName,"ID:",geneMapID[geneName]
					fieldsToPrint.append(options.NAString)
					
				meanCt=getMeanCtValueForGene(sampleStruct,geneName,options.NAString)
				if not options.useConventionalDeltaCt and meanCt!=options.NAString:
					meanCt=options.outputOffset-meanCt
				
				if meanCt!=options.NAString:
					if meanCt<0:
						print >> stderr,"warning: output < 0 output=",meanCt,"for sample",sampleName,"and gene",geneName,"ID:",geneMapID[geneName]
				
				#is this gene consistent?
				consistent=True
				for CtCall in CtCallForThisGene:
					if CtCall=="INC":
						consistent=False
			
				if consistent:
					if meanCt==options.NAString:
						if options.missingDataPolicy[0] in ["fill","fillifonevalid"]:
							meanCt=blankValue
				
					fieldsToPrint.append(str(meanCt))	
				else:
					fieldsToPrint.append(options.NAString)	
		#print >> stdout,"\t".join(fieldsToPrint)
		outputMatrix.append(fieldsToPrint)
		
	#TODO:whether to NA out whole col blank values?
	if options.missingDataPolicy[0]=='fillifonevalid': #because all are filled with blank values a priori

		#find all columns with all blank values and NA them.
		for c in range(1,len(outputMatrix)):
			allblanks=True
			for r in range(1,len(outputMatrix)):
				if outputMatrix[r][c]!=blankValue and outputMatrix[r][c]!=options.NAString:
					allblanks=False
					break
			
			if allblanks:
				#now set all row except header to NAString
				for r in range(1,len(outputMatrix)):
					outputMatrix[r][c]=options.NAString
	
	
	for outputRow in outputMatrix:
		print >> stdout,"\t".join(outputRow)
	
	
		
if __name__=='__main__':
	parser=OptionParser("usage: %prog [options] infile > outfile")
	parser.add_option("--control-names",dest="controlNames",default="gapdh,hprt",help="specify control gene names separate by comma [gapdh,hprt]")

	#parser.add_option("--not-collapse-replicate",dest="collapseReplicates",default=True,action="store_false",help="do not collapse replicates")
	parser.add_option("--gene-name-map",dest="geneNameMap",default=None,help="specify mapping of gene ID to gene Name [None] formatted: [geneID]tab[geneName]")
	
	#filters:
	parser.add_option("--keep-Ct-call-failed",dest="discardCtCallFailed",default=True,action="store_false",help="keep Ct Call failed entries  (applied to data and controls) [False]")
	parser.add_option("--Ct-quality-threshold",dest="CtQualityThreshold",default=-1.0,type=float,help="specify the CtQualityThreshold under which the value will be discarded (applied to data and controls) [-1.0:No treshold]")
	parser.add_option("--Ct-value-threshold",dest="CtValueThreshold",default=50.0,type=float,help="specify the max Ct value to be valid (applied to data and controls) [50.0]")
	parser.add_option("--max-Ct-deviation-between-replicates",dest="MaxCtRepDev",default=2.0,type=float,help="set the maximal Ct deviation between replicates allowed  (applied to data and controls) [2.0]")
	parser.add_option("--min-number-of-valid-data-point-per-control",dest="minValidReplicatesControl",default=1,type=int,help="set a minimum number of data points for that control [1]")
	parser.add_option("--Ct-value-threshold-for-per-control-average",dest="CtValueThresholdPerControl",default=30.0,type=float,help="set a Ct value threshold per control (average among replicates of that control) higher than which the whole row of data will be discarded [30.0]")
	parser.add_option("--min-number-of-valid-data-point-per-gene",dest="minValidReplicates",default=2,type=int,help="set a minimum number of data points for that gene (not including controls) in that sample for calling normalized value [2]")	
	parser.add_option("--Ct-value-threshold-for-data-average",dest="CtValueThresholdPerData",default=50.0,type=float,help="set a Ct value threshold of the data (average among replicates of that data before normalization to controls) higher than which the data will be discarded [50.0]")
	
	#normalization method
	parser.add_option("--output-ACx",dest="useConventionalDeltaCt",action="store_false",default=False,help="[default] output ACx = mean(Ct(Controls))-Ct(gene)+x where x is specified by --offset-output. This gives a higher value for higher expression and make the average control expression to as if it is x. i.e., the Ct value of a gene if the Ct value of the control is scaled to x. For getting Ratio Sample2/Sample1, do 2^(ACx(S2)-ACx(S1)) or operate as if in log2 space") 
	parser.add_option("--offset-output",dest="outputOffset",default=20.0,type=float,help="specify the output offset. i.e., the x in ACx values. [20.0]")	
	parser.add_option("--output-conventional-delta-ct",dest="useConventionalDeltaCt",default=False,action="store_true",help="output conventional delta Ct = Ct(gene)- mean(Ct(Controls)) such that higher delta Ct means lower expression relative to control. When getting Ratio Sample2/Sample1, then by delta delta Ct method, from this type of value it would be: 2^(DeltaCt(sample1)-DeltaCt(sample2))")
	

	
	
			
	#formatting modifiers
	parser.add_option("--fs",dest="fs",default="\t",help="specify the field separator of the infile [tab]")
	parser.add_option("--start-row",dest="startRow",default=2,type=int,help="specify the start row of data (excluding header) [2]")
	parser.add_option("--header-row",dest="headerRow",default=1,type=int,help="specify the header row of data [1]")
	parser.add_option("--sample-name-col",dest="sampleNameCol",default="2",help="specify the column of the infile that contains the sample name [2]")
	parser.add_option("--gene-name-col",dest="geneNameCol",default="5",help="specify the column of the infile that contains the gene name [5]")	
	parser.add_option("--Ct-value-col",dest="CtValueCol",default="7",help="specify the column of the infile that contains the Ct value [7]")
	parser.add_option("--Ct-quality-col",dest="CtQualityCol",default="8",help="specify the column of the infile that contains the Ct Call  [8]")
	parser.add_option("--Ct-call-col",dest="CtCallCol",default="9",help="specify the column of the infile that contains the Ct Call  [9]")
	
	#output modifiers:
	parser.add_option("--print-controls",dest="printControl",default=False,action="store_true",help="print also raw control Ct values [NO]")
	parser.add_option("--NA-string",dest="NAString",default="NA",help="set the output string for invalid data [NA]")
	parser.add_option("--control-prefix",dest="controlPrefix",default="control.",help="add a prefix to control genes [control.]")	
	parser.add_option("--missing-data-policy",dest="missingDataPolicy",default=["fill","0.0","50.0"],nargs=3,help="set how missing data are outputed for the final matrix (3 args): <nofill|*fill|fillifonevalid> <valueACxToFill> <valueDeltaCtToFill> default:fill 0.0 50.0")
	

	
	
	(options, args) = parser.parse_args(argv)
	
	try:
		infile,=args[1:]
		#print >> stderr,infile
	except:
		parser.print_help()
		exit()
	
	geneNamesInOrder=[]
	
	geneNameMap=dict()
	geneMapID=dict()

	if options.geneNameMap:
		#if a gene name map is specified, usu is, then load it
		print >> stderr,"load gene map from",options.geneNameMap
		fil=open(options.geneNameMap)
		for lin in fil:
			fields=lin.strip().split("\t")
			geneNameMap[fields[0]]=fields[1]
			if fields[1] in geneMapID:
				geneMapID[fields[1]]+=","+fields[0]
			else:
				geneMapID[fields[1]]=fields[0]
			geneNamesInOrder.append(fields[1])
			
		fil.close()	
	
	#now process arguments
	options.controlNames=options.controlNames.split(",")
	
	#now load the infile
	sampleNames,sampleStructs=readFluidigmFile(infile,options.sampleNameCol,options.geneNameCol,options.CtValueCol,options.CtQualityCol,options.CtCallCol,options.headerRow,options.startRow,options.fs,geneNameMap)
	
	#filter data
	invalidSampleRows=filterFluidigmData(sampleNames,sampleStructs,options,geneMapID)
	
	#print >> stderr,"invalid sample rows",sorted(invalidSampleRows)
	
	#now normalize data
	normalizeFluidigmData(sampleNames,sampleStructs,options,invalidSampleRows)
	
	#now output
	printFluidigmData(sampleNames,sampleStructs,options,invalidSampleRows,geneNamesInOrder,options.printControl,geneMapID)
	