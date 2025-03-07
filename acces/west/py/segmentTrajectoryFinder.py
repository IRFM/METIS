import os
import sys
import xml.etree.ElementTree as Et
import numpy as np
import matplotlib.pyplot as plt



def segmentTrajectoryFinder(infile):
	"""
	SCENARIO = segmentTrajectoryFinder(FILE)
	
	Determines the nominal scenario of a pulse from a DCS-format FILE.
	
	FILE is is the path of a DCS-formatted supervision (typically Sup.xml) file.
	
	SCENARIO is a n x 3 array describing the nominal scenario. 
	
	Example : 
	>> SCENARIO = segmentTrajectoryFinder(FILE)
	
	>> SCENARIO = 
	[['Init' '34.0' '34.0']
	 ['segment1' '34.1' '0.1']]
	
	The first column is the name of the segment as defined in Xedit by the session leader
	The second column is the absolute time at which the segment ends (the first segment starts 
	at 0.0s. The third column is the duration of the segment. Note that the duration information 
	is theorticaly redundant but is kept for practical purposes to avoid making substractions
	every time the duration information is needed.
	
	Note that the nominal scenario is, by convention, the sequence of segments which contain
	the most segments with number less than a 100. Segments with number > 100 are (also by 
	convention) either backup segments or stop segments.
	
	
	Bug reports and suggestions are welcome.
	
	---
	C. Reux - January 2017. Based on initial work by H. Joshi
	"""
	

	## getting root of xml file 
	tree= Et.parse(infile)
	root=tree.getroot()


	# arrMainSeg to store main segment's name. Only temporary structure.
	arrMainSeg = np.array([])
	
	# arrTargetSeg stores all the segments and all the ways to transit from them to other segments.
	# Contains the segment name, the target segment, the name of the jump condition, the (maximum) time
	# at which the jump can take place.
	arrTargetSeg = np.array([])
	
	arrSearchSeg={}
	# arrsearchSeg draws a list containing each segment (once) and all the segments towards a jump
	# is possible from this particular segment.
	
	# arrsearch stores all the possible scenarios (=sequence of segments). 
	# It contains a list of string arrays with segment names interleaved with a '->'
	# For the moment, the search does not take into account the possibility to transit at times others
	# than the maximum duration.
	arrSearch = np.array([])
	
	# arrSearch2 contains the segments names of the nominal scenario.
	arrSearch2 = np.array([])
	
	#arrSegment stores the final result (nominal scenario) as described in the function help text.
	#arrSegment - final result stored in this array 
	arrSegment = np.array([]) 
	
	
	###  functions are called as below 
	
	# Determination of all segments and transition ways
	arrTargetSeg=set_target_segment(arrMainSeg,arrTargetSeg,root)

	#print_target_segment(arrTargetSeg)

	### reshape arrTargetSeg from a single line vector to an array.
	arrColumns=5
	arrTargetSegElements=len(arrTargetSeg)
	arrRows=int(arrTargetSegElements/arrColumns)
	arrTargetSeg = np.reshape(arrTargetSeg,(arrRows,arrColumns))


	### create segment data into searchable form with a unique instance of each originating segment
	# in the array.
	arrSearchSeg={}
	arrSearchSeg=create_seg_traj(arrTargetSeg)
	
	
	# Search all possible scenarios
	arrSearch = search_all(arrSearchSeg,'Init','TheEnd')


	
	# Search the nominal scenario among all possible scenarios
	# The nominal scenario is defined as the one containing the most segments
	# with ID numbers < 100
	
	longest=0
	# Loop over all possible scenarios
	for i in range(0,len(arrSearch)):
		arrSearch_temp = prepare_arrSearch(arrSearch,i)
		count_nominal = 1
		# Loop over the segments in the chosen scenario
		for ii in range(1,len(arrSearch_temp)):
			# Count the number of segments with ID < 100
			segment_id = float(arrSearch_temp[ii][7:])
			if segment_id<100:
				count_nominal = count_nominal+1
		shortest = count_nominal
		if( shortest > longest):
			longest = shortest
			choice = i


	# Select the nominal scenario
	arrSearch2 = prepare_arrSearch(arrSearch,choice)
	

	arrSegment=prepare_arrSegment(arrSearch2,arrTargetSeg)

	arrSegment=format_arrSegment(arrSegment)	

	return(arrSegment)


def set_target_segment(arrMainSeg,arrTargetSeg,root):
	segSr=1
	#for SegmentList in root.findall('{http://www.ipp.mpg.de/2004/Schema/SUP_Desc}SegmentList'):
	# Loop on the segment list
	for SegmentList in root.findall('SegmentList'):
		for SEGMENT in SegmentList:
			# Get the segment name and append it in the list
			name = SEGMENT.get('name')
			arrMainSeg=np.append(arrMainSeg,name)
			targetList = np.array([])
			# Search for the watchdog jump condition (i.e. the longest possible duration
			# of the segment)
			for Control in SEGMENT:
				if (Control.get('program')=='Watchdog'):
					# Careful: startTime is actually the time at which the watchdog
					# jump condition kicks in. This is actually the end of the segment, 
					# not its beginning 
					startTime = Control.get('startTime')
			
			# Loop over other jump conditions
			for Control in SEGMENT:
				targetSegment = Control.get('targetSegment')
				program = Control.get('program')
				
				# Avoid redundant scenarios if several jump conditions lead to the same 
				# segments. Note that in reality, their startTime may be different. This
				# is not taken into account.
				if not targetSegment in targetList:
	#adding target segment detail in arrTargetSeg array
				#if (program == "Watchdog" or program == "StartUp" ):
					targetList = np.append(targetList,targetSegment)
					arrTargetSeg=np.append(arrTargetSeg,name)
					arrTargetSeg=np.append(arrTargetSeg,targetSegment)
					arrTargetSeg=np.append(arrTargetSeg,program)
					arrTargetSeg=np.append(arrTargetSeg,startTime)
					arrTargetSeg=np.append(arrTargetSeg,segSr)
					segSr+=1
	return arrTargetSeg



def search_all(arrSearchSeg, start, end, path=[]):
        # Recursive search of all scenarios
        ar=['->']
        path = path + [start]+ar
        if start == end:
            return [path]

        if start not in arrSearchSeg:
            return []

        paths = []
        for node in arrSearchSeg[start]:
            if node not in path:
                newpaths = search_all(arrSearchSeg, node, end, path)
                for newpath in newpaths:
                    paths.append(newpath)
        return paths


def create_seg_traj(arrTargetSeg):
	# Creates an array which segment names as indexes to avoid only one occurrence of each segment.
	# arrSearSeg gives the segments which it is possible to jump to starting from the index segment.
	arrSearchSeg = {}
	for i in range(0,len(arrTargetSeg)):
		key = arrTargetSeg[i][0]
		#print(key)
		#arrSearchSeg.setdefault(key,[]).append("%s via %s"%(arrTargetSeg[i][1],arrTargetSeg[i][2]))
		arrSearchSeg.setdefault(key,[]).append(arrTargetSeg[i][1])

	return arrSearchSeg

def prepare_arrSearch(arrSearch,choice):
	"""
	ARRSEARCH2 = prepare_arrSearch(ARRSEARCH,CHOICE)
	Selects a given scenario among the list of possible scenarios 
	ARRSEARCH is an array of scenarios formatted by the search_all function 
	CHOICE is the number of the scenario in arrSearch
	"""
	arrSearch2 = np.asarray(arrSearch[choice])

	## delete -> and "TheEnd" from arrSearch
	index = np.where(arrSearch2 == "->")
	arrSearch2 = np.delete(arrSearch2,index)
	index = np.where(arrSearch2 == "TheEnd")
	arrSearch2 = np.delete(arrSearch2,index)
	return arrSearch2

def prepare_arrSegment(arrSearch,arrTargetSeg):
	arrSegment = np.array([])	
	#### After selection of trajectory  and fill arrSegment with segments and time information
	#### re declare arrSegment array to freshly append selected choice of trajectory
	startTime=0.0
	for i in range(0,len(arrSearch)):
		segment_done = 0
		for j in range(0,len(arrTargetSeg)):	
			# Note: only one transition is considered at the moment
			# This assumes that all transitions happen at Watchdog time.
			if (arrSearch[i] == arrTargetSeg[j][0]) and not segment_done:
				startTime = startTime+float(arrTargetSeg[j][3])
				arrSegment=np.append(arrSegment,arrSearch[i])
				arrSegment=np.append(arrSegment,startTime)
				arrSegment=np.append(arrSegment,float(arrTargetSeg[j][3]))
				segment_done = 1

	return arrSegment
	
def format_arrSegment(arrSegment):
	##### to check number of elements in arrSegment array
	countarrSegmentElements = len(arrSegment)
	### arrColumns = 3 because arrSegment array has two columns : 1. name of segment and 2. start time
	arrColumns =3
	arrRows = (countarrSegmentElements/arrColumns)
	arrSegment = np.reshape(arrSegment,(arrRows,arrColumns))

	return arrSegment


### Functions declaration over





