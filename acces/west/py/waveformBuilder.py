import os
import xml.etree.ElementTree as Et
import numpy as np
import matplotlib.pyplot as plt
import sys
import math
import waveform_pack


def waveformBuilder(arrSegment,signal_list,infile):
	"""
	OUTPUT_WAVEFORM_PACK = waveformBuilder(ARRSEGMENT,SIGNAL_LIST,INFILE)
	
	Builds waveforms of a given signals following a chosen scenario.
	
	ARRSEGMENT is the chosen scenario as provided by segmentTrajectoryFinder.py
	It is a numpy string array and should be formatted this way:
	[[segmentname1 AbsoluteEndtime1 Duration1]
	[segmentname2 AbolsuteEndtime2 Duration2]
	[segmentname3 AbsoluteEndtime3 Duration3]
	[... ... ...]]
	
	Example:
	[['Init' '34.0' '34.0']
	['segment1' '35.5' '1.5']
	['segment2' '37.0' '2.0']]
	
	All times in seconds.
	See the documentation of the segmentTrajectoryFinder.py for more details.
	
	SIGNAL_LIST is a numpy string array of the signals for which the waveforms are to be built.
	Signal names are DCS names.
	 
	INFILE is the path of the DCS xml containing the waveform information
	Usual name is DP.xml
	
	OUTPUT_WAVFORM_PACK is an array of waveform_packs objects. It contains the reconstructed 
	waveform data as well as their position in the segment sequence.
	The class is defined in waveform_pack.py. 
	
	OUTPUT_WAVEFORM_PACK.times is the absolute time vector. Numpy float array
	OUTPUT_WAVEFORM_PACK.reltimes is the relative time vector (i.e. the time of the point inside 
	the segment it belongs to. Numpy float array.
	OUTPUT_WAVEFORM_PACK.values is the data vector. Units are defined by DCS and are 
	usually SI. Numpy float array.
	OUTPUT_WAVEFORM_PACK.segments is a vector indicating which segment each data point 
	belongs to. Numpy string array.
	
	Features currently implemented:
	- Linear or Step-executed waveforms
	- Absolute or relative waveforms (note: in DCS, relative means relative to the start of the segment)
	
	Features not yet implemented:
	- Global forward-transition-times for waveforms (JET-like feature). Not yet implemented in Xedit anyway.
	- Proportional control for Inh/b coils (proportional means proportional to Ip). 
	  The reconstructed waveform in this case will only be the proportional factor F in kA/MA, not Ip*F in kA
	- Envelopes proportional to the reference theyr are related are not functional.
	
	Bug reports and suggestions are welcome.
	
	--
	C. Reux - Januray 2017. Based on initial work by H. Joshi
	
	"""
	
	# xml file parsing.
	tree = Et.parse(infile)
	xmlroot = tree.getroot()
	
	# Create the array of waveform_pack objects for final result storage.
	output_waveform_pack=np.array([],dtype='object')
		
	### Declare array to store signal data before casting it in a waveform_pack object. 
	arrsignalTrajectory = np.array([])
	
	### Declare array to store final result
	#arrPlot=np.array([])
	
	### call functions
	
	### Loop over all signals requested by the user.
	for i in range(0,len(signal_list)):

		# Initialization of a single waveform_pack object, to be appended with others in case several
		# signals are requested.
		single_waveform_pack = waveform_pack.waveform_pack()
		
		### Declare array to store signal default data from Signal definition in DP file
		arrSignal= np.array([])
		
		# Look for signal default and general values in the DCS xml file
		arrSignal=search_signal(arrSignal,signal_list[i],xmlroot)
		

		# Find all data points for the requested signal in all segments.
		arrsignalTrajectory=search_signal_trajectory(signal_list[i],xmlroot)


		
		# Reshape the long 1D array arrsignalTrajectory into a multi-line array for easier handling.
		arrsignalTrajectory=format_signal_trajectory(arrsignalTrajectory)


		# Waveform assembly from all info gathered previously.
		[single_waveform_pack.segments,single_waveform_pack.times,single_waveform_pack.reltimes,single_waveform_pack.values]=\
		prepare_final(xmlroot,arrSegment,arrsignalTrajectory,arrSignal)
		
		# Append the current waveform for the current signal to those already calculated and
		# stored in the waveform_pack array.
		output_waveform_pack =  np.append(output_waveform_pack,single_waveform_pack)
			
	return(output_waveform_pack)
	



### Functions Declaration


def search_signal(arrSignal,searchSignal,xmlroot):
	"""
	Retrieves the signal declaration in the xml file: default values, type, dimension, etc.)
	"""
	 
	for SIGNAL in xmlroot.findall('DECLARATIONS/SIGNALS/SIGNAL'):
		name = SIGNAL.attrib['name']
		if (name == searchSignal):
			signal_type    = SIGNAL.attrib['signal_type']
			type           = SIGNAL.attrib['type']
			dimension      = SIGNAL.attrib['dimension']
			for DEFAULT_VALUE in SIGNAL:
				default_value  = (DEFAULT_VALUE.text).strip('\n\ \r\t')
				
				arrSignal = np.append(arrSignal, name)
				arrSignal = np.append(arrSignal, signal_type)
				arrSignal = np.append(arrSignal, type)
				arrSignal = np.append(arrSignal, dimension)
				arrSignal = np.append(arrSignal, default_value)
				
	return (arrSignal)




def search_signal_trajectory(searchSignal,xmlroot):
	"""
	Finds all occurrences of a signal in all segments.
	"""

	# Initialization of the output array for all occurences of a signal in segments.
	arrsignalTrajectory = np.array([])
	
	# Loop over the segments
	for SEGMENT in xmlroot.findall("SEGMENTS/SEGMENT"):
		# Enable the following to check segment detail
		# print(SEGMENT.attrib,SEGMENT.tag)
		segment_name = SEGMENT.attrib['Id']
		
		# Initial confusion between segment Id and segment name. 
		# Aims at keeping the identification set by SLs in Xedit
		# Rename the first segment using the identifier
		if segment_name=='0':
			segment_name='Init'
		else:
			segment_name='segment%s' %(segment_name)
		
		# Loop over all the waveforms available for a given the segment to find the one 
		# the user is interested in.		
		for SIGNAL_TRAJECTORIES in SEGMENT:
			for SIGNAL_TRAJECTORY in SIGNAL_TRAJECTORIES:
				# Assign signal name to signalTrajectory variable
				signalTrajectory = SIGNAL_TRAJECTORY.attrib['Name']
				# Execution rule is step (mostly for flags) or Linear (for normal references)
				executionRule = SIGNAL_TRAJECTORY.find('EXECUTION_RULE').attrib['is']
				# Look for the signal we are interested in.
				if (searchSignal == signalTrajectory):
					# For the desired signal, gather all time points in the segment and their related information
					for REFERENCE in SIGNAL_TRAJECTORY:
						for POINT in REFERENCE:
							for VALUE in POINT:
								# Extract all information about a data point
								value = VALUE.text
								arrsignalTrajectory = np.append(arrsignalTrajectory,segment_name)
								arrsignalTrajectory = np.append(arrsignalTrajectory,POINT.get("time"))
								arrsignalTrajectory = np.append(arrsignalTrajectory,VALUE.get("type"))
								arrsignalTrajectory = np.append(arrsignalTrajectory,VALUE.get("dimension"))
								arrsignalTrajectory = np.append(arrsignalTrajectory,(VALUE.text).strip('\n\r\ \t'))
								arrsignalTrajectory = np.append(arrsignalTrajectory,signalTrajectory)
								arrsignalTrajectory = np.append(arrsignalTrajectory,executionRule)

	return arrsignalTrajectory


### Printing signal detail stored in arrsignalTrajectory array

def format_signal_trajectory(arrsignalTrajectory):
	"""
	Reshapes arrsignalTrajectory from a 1D array into a mult-line array for easier handling.
	"""
	countTrajectoryElements = len(arrsignalTrajectory)
	### Value 7 is taken because there are six (6) different values stored in arrsignalTrajectory
	arrColumns = 7
	arrRows = int(countTrajectoryElements/arrColumns)

	#reshape signal trajectory with rows and columns
	arrsignalTrajectory  = np.reshape(arrsignalTrajectory,(arrRows,arrColumns))

	return arrsignalTrajectory


### preparing waveform_pack with final result combining arrSegment and arrsignalTrajectory
def prepare_final(xmlroot,arrSegment,arrsignalTrajectory,arrSignal):
	"""
	Assembly of the final waveform. Features: Takes into account step or linear executed waveforms. 
	"""


	#PCS time cycle in seconds. Used for jumps at the beginning of a segment	
	time_cycle = 0.002
	
	# Number of explicit data points in all segments.
	rowsSignal = int(len(arrsignalTrajectory))
	
	# Number of segments in the scenario.
	rowsSegment = int(len(arrSegment))

	#  Output arrays containing the waveform data points 
	# (absolute time, relative time, data and associated segment).
	times_output = np.array([])
	relative_times_output = np.array([])
	values_output = np.array([])
	segments_output = np.array([])
	segmentTime = 0.0

	

	
	# Iterate segments and build the waveform segment by segment.
	for j in range(0,rowsSegment):
			
		segmentName = arrSegment[j,0]
		
		# Guard test to make sure there is at least one explicit data point.
		# Otherwise there is no point determining the waveform type (relative/absolute)
		# or try to build the waveform (default value only)
		if (len(arrsignalTrajectory)>0):
			waveformType = waveformTypefinder(xmlroot,segmentName,arrsignalTrajectory[0][5])
		else:
			break
		
		# First point of the first segment. Whatever signalTrajectory contains, the start of a segment is considered
		# as a real point by the PCS. 
		# On any other segment than the first, the signal value at the end of the previous segment is used.
		# On the first segment, the signal default value is used.
		# Note: LastPointIn is used in case points have been defined outside the segment duration. See below.			
		if j==0:
			# Start time of the segment
			segmentTime = 0.0
			
			# Default value from the signal declaration.
			waveformValue = float(arrSignal[4])
			
			# Segment the point belongs to (i.e. the current segment)
			segments_output = np.append(segments_output,segmentName)
			
			# Absolute time
			times_output = np.append(times_output,segmentTime)
			
			#Relative time
			relative_times_output = np.append(relative_times_output,segmentTime)
			
			# Append the value to the vector
			values_output = np.append(values_output,waveformValue)
			
			# LastValue is used for relative waveforms (can be on the segment immediately after)
			lastValue = waveformValue
			
			# Starting value of the segment. Used for relative waveforms
			startValueSegment = waveformValue
			
			# Same remark as lastValue
			lastTime = segmentTime
			
			# Flag to check if the lastPoint in the segment is inside its time boundaries.
			# Used to ignore any subsequent point after the single one which could be defined outside.
			# See below.
			lastPointIn = True
			#print('First point of first segment')
		
		# First point of any segment except the first one.	
		else:

			
			segmentTime = float(arrSegment[j-1,1])
			segments_output = np.append(segments_output,segmentName)
			times_output = np.append(times_output,segmentTime)
			# First tpoint of the segment. Therefore the relative time is 0.0s.
			relative_times_output = np.append(relative_times_output,0.0)
			values_output = np.append(values_output,waveformValue)
			lastValue = waveformValue
			startValueSegment = waveformValue
			lastTime = segmentTime
			lastPointIn = True
			#print('First point of any segment except first')
		
		# Other points in the segment (i.e. not the first point)
		# Loop over the explicit points defined in the signal trajectory to search if
		# one of them belongs to the current segment.		
		for i in range(0,rowsSignal):
			
			#print("arrsignaltraject = %s"%arrsignalTrajectory[i,0])
			#print("arrSegment = %s"%arrSegment[j,0])	
			
			# Case where an explicit point exists in the current segment.		
			if (arrsignalTrajectory[i,0] == arrSegment[j,0]):
				#print('Found a point in segment %s. Point value %s. Point time raw rel %s '%(arrSegment[j,0],arrsignalTrajectory[i,4],arrsignalTrajectory[i,1]))
					
				# Case for the first segment (absolute time is zero)
				if (j==0):
					segment_duration = float(arrSegment[j,1])
				# Any other segment
				else:
					#segment_duration = float(arrSegment[j,1])-float(arrSegment[j-1,1])
					segment_duration = float(arrSegment[j,2]) #Mod CR 20170517. Old calculation had a rounding error.
					
				
					
				# Case of a normal point defined within the segment time boundaries.	
				if (float(arrsignalTrajectory[i,1])<=segment_duration):
				
					# Special case if the first defined point in the segment is also the beginning of the segment
					# (or one less than one time step away from it)
					# In this case, a jump is done from the value in the last segment to the current value within
					# the duration of one time cycle.				
					if float(arrsignalTrajectory[i,1])<time_cycle:
						waveformTime_rel = float(arrsignalTrajectory[i,1])+time_cycle
						#print('First point in segment is the beginning of the segment')
						#print('FPITBOS: %d'%(waveformTime_rel))
						
					else:
						waveformTime_rel = float(arrsignalTrajectory[i,1])	
						#print('First point in segment is NOT the beginning of the segment')

					segmentName = arrSegment[j,0]
					
					# Case for the first segment: the end time of the previous segment (=start time of the 
					# current one) is not yet defined.
					# In this case, the starting time is imply 0.0s
					if j==0:
						segmentTime = 0.0
					# Other segments. The starting segment time is used extensively to build the absolute time vectors.
					else:
						segmentTime = float(arrSegment[j-1,1])
					#print('Normalpoint: segmentTime = %d' %(segmentTime))

					# Treatment of execution mode: case for Linear interpolation in between points					
					if (arrsignalTrajectory[i,6]=='Linear') or (arrsignalTrajectory[i,6]=='linear'):
						
						# Treatment of the wavefor type (absolute/relative)
						# Case for absolute waveform
						if waveformType:
							waveformValue = float(arrsignalTrajectory[i,4])
						# Case for relative wavefor(i.e. relative to the value at the beginning of the segment
						#In the PCS, the measured value is used rather than the reference value (obviously not available
						#in the pre-analysis tool).
						else:
							waveformValue = float(arrsignalTrajectory[i,4])*startValueSegment
							
						
						waveformTime = segmentTime + waveformTime_rel
						segments_output = np.append(segments_output,segmentName)
						times_output = np.append(times_output,waveformTime)
						relative_times_output = np.append(relative_times_output,waveformTime_rel)
						values_output = np.append(values_output,waveformValue)
						lastValue = waveformValue
						lastTime = waveformTime
						lastPointIn = True
						
						# Case 2: Step interpolation in between points. The jump is done in one time cycle.
						# An escape condition is added in case the point time coincides with the end of the segment.
						# Keeping such points would finish the step in the next segment, which builds a non-monotonous
						# time vector. The point is therefore rejected. TO BE COMPARED WITH WHAT DCS ACTUALLY DOES.
						# Note that two points are added to the waveform: the one juste before the step, and 
						# another one at the top of the step.
					elif ((arrsignalTrajectory[i,6]=='Step') or (arrsignalTrajectory[i,6]=='step')) and (waveformTime_rel+segmentTime)<(float(arrSegment[j,1])-time_cycle):
						waveformValue = lastValue
						waveformTime = segmentTime + waveformTime_rel
						segments_output = np.append(segments_output,segmentName)
						times_output = np.append(times_output,waveformTime)
						relative_times_output = np.append(relative_times_output,waveformTime_rel)
						values_output = np.append(values_output,waveformValue)
						if waveformType:
							#Case Absolute waveform
							waveformValue2 = float(arrsignalTrajectory[i,4])
						else:
							#Case Relative waveform
							waveformValue2 = float(arrsignalTrajectory[i,4])*startValueSegment
							
							
						waveformTime2 = waveformTime+time_cycle
						segments_output = np.append(segments_output,segmentName)
						times_output = np.append(times_output,waveformTime2)
						relative_times_output = np.append(relative_times_output,waveformTime_rel+time_cycle)
						values_output = np.append(values_output,waveformValue2)
						lastValue = waveformValue2
						lastTime = waveformTime2
						lastPointIn = True
						
				# Case where a linear-executed point is defined after the end of the segment. The linear ramp 
				# between the last two points is started but interrupted before the end of the segment.
				# A point is created at the very end of the segment. 
				# Any subsequent point is ignored.
				# Step-executed points defined after the end of the segment are obviously ignored.				

				
				elif ((float(arrsignalTrajectory[i,1])>segment_duration) and lastPointIn and ((arrsignalTrajectory[i,6]=='Linear') or (arrsignalTrajectory[i,6]=='linear'))):
					segmentName = arrSegment[j,0]
					segmentTime = float(arrSegment[j-1,1])
					
					waveformTime_rel = float(arrsignalTrajectory[i,1])
					waveformOutTime = segmentTime + waveformTime_rel
					waveformOutValue = float(arrsignalTrajectory[i,4])
					
					waveformTime = float(arrSegment[j,1])
					
					# Interpolation at the segment end time between the last available point and the point defined 
					# outside the current segment.
					if waveformType:
						#Case Absolute value
						waveformValue = np.interp(waveformTime, [lastTime, waveformOutTime],[lastValue,waveformOutValue])
					else:
						#Case relative value
						waveformValue = np.interp(waveformTime, [lastTime, waveformOutTime],[lastValue,waveformOutValue])*startValueSegment
						print('relative Linear detected')
					
					segments_output = np.append(segments_output,segmentName)
					times_output = np.append(times_output,waveformTime)
					relative_times_output = np.append(relative_times_output,waveformTime-segmentTime)
					values_output = np.append(values_output,waveformValue)
					lastValue = waveformValue
					lastTime = waveformTime
					#LastPointIn is flagged false to ignore all subsequent points.
					lastPointIn = False
					
		# Last point in the last segment. In case no point is defined at the end, a extra point is added
		# only for waveform completeness.
		if (j==(rowsSegment-1) and not (lastTime==float(arrSegment[j,1]))):
			segmentName = arrSegment[j,0]
			waveformTime = float(arrSegment[j,1])
			waveformValue = lastValue
			segments_output = np.append(segments_output,segmentName)
			times_output = np.append(times_output,waveformTime)
			if j==0:
				#Case where the last segment is also the first one (only one segment)
				relative_times_output = np.append(relative_times_output,waveformTime)
			else:
				relative_times_output = np.append(relative_times_output,waveformTime-float(arrSegment[j-1,1]))
			values_output = np.append(values_output,waveformValue)

		

	return(segments_output,times_output,relative_times_output,values_output)
	



def waveformTypefinder(xmlroot,current_segment,sigTrajectoryName):
	"""
	Determines the waveform type (absolute or relative to the start value of the segment).
	Returns a flag for the waveform type: absolute (1), relative (0) or proportional (2).	
	The value is stored in a an int value for each segment and each component
	The int is a bitmask storing the type of a number of waveforms in the component:
	absolute (1), relative (0). Relative means actually proportional for PF coils.
	Which bit corresponds to which waveform is given in the comments of the DP.xml file.
	NOTE: the bit 0 in the DP.xml comments is the bit of least weight.
	NOTE: the proportional type (2) is correctly detected in the present function
	but not yet implemented in the waveform builder.
	"""
	
	TypeBit_SignalNameList = ['rts:WEST_PCS/Actuators/Poloidal/ref_type.ref',\
	'rts:WEST_PCS/Actuators/Gas/ref_type.ref',\
	'rts:WEST_PCS/Actuators/Heating/ref_type.ref',\
	'rts:WEST_PCS/Actuators/Poloidal/env_type.ref',\
	'rts:WEST_PCS/Actuators/Gas/env_type.ref',\
	'rts:WEST_PCS/Actuators/Heating/env_type.ref']
	# The size of the bitmask cannot be guessed by the algorithm. 
	# It is needed in order not assign values to extra zeroes 
	# if the integer was capped to fixed number of bits.
	# Each of the bitmask sizes refers to the 6 signals referenced above.
	TypeBit_SizeBit = [27,24,10,27,28,10]

	
	TypeBit_SignalList = np.array([])
	# Iterate on the number of waveform type bitmasks.
	k=0
	for TypeBit_SignalName in TypeBit_SignalNameList:

		# Loop over all segments in the scenario to find the one requested by the user.
		for SEGMENT in xmlroot.findall("SEGMENTS/SEGMENT"):
			# Location of signal in segment
			# Enable following to check segment detail, which helps in error solving
			#print(SEGMENT.attrib,SEGMENT.tag)
			segment_name = SEGMENT.attrib['Id']
			# As usual, renaming of the first segment.
			if segment_name=='0':
				segment_name='Init'
			else:
				segment_name='segment%s' %(segment_name)
			
			# If the segment is the one requested by the arguments, loop over the signal trajectories
			# to find the signal requested by the arguments (in Signal Trajectory). No point going 
			# all the time-consuming mapping of only a few waveforms are present in each segment.
			if (current_segment==segment_name):	
				for SIGNAL_TRAJECTORIES in SEGMENT:
					for SIGNAL_TRAJECTORY in SIGNAL_TRAJECTORIES:
						# Assign the signal name explored by the loop to signalTrajectory
						signalTrajectory = SIGNAL_TRAJECTORY.attrib['Name']
						
						# matching searchSignal's waveform name (nnnn/env_tye.ref) with signalTrajectory's value
						# When it is found, extract he bit mask
						if (TypeBit_SignalName == signalTrajectory):
							for REFERENCE in SIGNAL_TRAJECTORY:
								# Theoretically, there is only one point for the whole segment, but the
								# whole structure (point, value, etc.) is still used.
								for POINT in REFERENCE:
									for VALUE in POINT:
										
										# Extract the decimal value of the bitmask
										value_decimal = VALUE.text	
					
										# Format the string taking into account the size of the bitmask.
										# Example: 244 should be actually 00000244 if 8 bits are treated.
										format_str = '#0%ib'%(TypeBit_SizeBit[k]+2)
										
										# Convert to binary array
										value_binary = format(int(value_decimal),format_str)[2:]

										TypeBit_SignalList = np.append(TypeBit_SignalList,value_binary)
										
										# Enable the following if command lines to display the bitmask name (debug	
										#if (TypeBit_SignalName == 'rts:WEST_PCS/Actuators/Heating/ref_type.ref'):
										#	print(TypeBit_SignalName)
										#	print(VALUE.text)
										#	print(format_str)
										#	print(value_binary)
										#	print(TypeBit_SizeBit[k]+2)
										#	print(k)
										# End debug lines
										
										
		k=k+1
	
	
	# Now the bitmask is known, search the relevant bit for the signal given as input.
	
	if (sigTrajectoryName=='rts:WEST_PCS/Plasma/Ip/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[0][11])
	elif (sigTrajectoryName=='rts:WEST_PCS/Plasma/dXLow/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[0][12])
	elif (sigTrajectoryName=='rts:WEST_PCS/Plasma/Z/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[0][14])	
	elif (sigTrajectoryName=='rts:WEST_PCS/Plasma/rext/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[0][13])		
	elif (sigTrajectoryName=='rts:WEST_PCS/Plasma/Zgeo/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[0][14])
	elif (sigTrajectoryName=='rts:WEST_PCS/Plasma/Rgeo/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[0][15])
	# Special "relative" Flag for poloidal coil currents: 2 because they are proportional to Ip
	# Arithmetic Trick so that the Absolute flag (1) stays the same but 0 becomes 2 
	# NOTE: although correctly detected here, Proportional waveforms are not yet implemented in
	# the waveform builder.
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Poloidal/IXh/waveform.ref'):
		AbsRelBitValue = abs(int(TypeBit_SignalList[0][17])-2)		
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Poloidal/IXb/waveform.ref'):
		AbsRelBitValue = abs(int(TypeBit_SignalList[0][16])-2)	
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Poloidal/IBh/waveform.ref'):
		AbsRelBitValue = abs(int(TypeBit_SignalList[0][25])-2)	
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Poloidal/IDh/waveform.ref'):
		AbsRelBitValue = abs(int(TypeBit_SignalList[0][24])-2)	
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Poloidal/IEh/waveform.ref'):
		AbsRelBitValue = abs(int(TypeBit_SignalList[0][23])-2)	
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Poloidal/IFh/waveform.ref'):
		AbsRelBitValue = abs(int(TypeBit_SignalList[0][22])-2)	
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Poloidal/IFb/waveform.ref'):
		AbsRelBitValue = abs(int(TypeBit_SignalList[0][21])-2)	
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Poloidal/IEb/waveform.ref'):
		AbsRelBitValue = abs(int(TypeBit_SignalList[0][20])-2)	
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Poloidal/IDb/waveform.ref'):
		AbsRelBitValue = abs(int(TypeBit_SignalList[0][19])-2)
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Poloidal/IBb/waveform.ref'):
		AbsRelBitValue = abs(int(TypeBit_SignalList[0][18])-2)
	# Not for A coil: its relative mode is a real relative mode.
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Poloidal/IA/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[0][26])
		
		
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Gas/REF3/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[1][0])
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Gas/REF2/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[1][1])		
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Gas/REF1/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[1][2])
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Gas/valve21/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[1][3])		
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Gas/valve20/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[1][4])						
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Gas/valve19/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[1][5])			
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Gas/valve18/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[1][6])			
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Gas/valve17/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[1][7])		
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Gas/valve16/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[1][8])		
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Gas/valve15/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[1][9])		
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Gas/valve14/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[1][10])			
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Gas/valve13/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[1][11])		
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Gas/valve12/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[1][12])			
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Gas/valve11/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[1][13])		
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Gas/valve10/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[1][14])		
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Gas/valve09/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[1][15])		
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Gas/valve08/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[1][16])		
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Gas/valve07/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[1][17])	
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Gas/valve06/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[1][18])	
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Gas/valve05/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[1][19])	
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Gas/valve04/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[1][20])	
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Gas/valve03/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[1][21])	
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Gas/valve02/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[1][22])	
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Gas/valve01/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[1][23])	
		
	
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Heating/ICRH/phase/3/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[2][0])
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Heating/ICRH/power/3/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[2][1])		
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Heating/ICRH/phase/2/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[2][2])
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Heating/ICRH/power/2/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[2][3])	
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Heating/ICRH/phase/1/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[2][4])
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Heating/ICRH/power/1/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[2][5])		
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Heating/LHCD/phase/2/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[2][6])
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Heating/LHCD/power/2/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[2][7])			
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Heating/LHCD/phase/1/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[2][8])
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Heating/LHCD/power/1/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[2][9])

	
	elif (sigTrajectoryName=='rts:WEST_PCS/Plasma/Ip/min_env_waveform.ref')\
	or (sigTrajectoryName=='rts:WEST_PCS/Plasma/Ip/max_env_waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[3][11])				
	elif (sigTrajectoryName=='rts:WEST_PCS/Plasma/dXLow/min_env_waveform.ref')\
	or (sigTrajectoryName=='rts:WEST_PCS/Plasma/dXLow/max_env_waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[3][12])
	elif (sigTrajectoryName=='rts:WEST_PCS/Plasma/Z/min_envwaveform.ref')\
	or (sigTrajectoryName=='rts:WEST_PCS/Plasma/Z/max_env_waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[3][14])	
	elif (sigTrajectoryName=='rts:WEST_PCS/Plasma/rext/min_env_waveform.ref')\
	or (sigTrajectoryName=='rts:WEST_PCS/Plasma/rext/max_env_waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[3][13])		
	elif (sigTrajectoryName=='rts:WEST_PCS/Plasma/zgeo/min_env_waveform.ref')\
	or (sigTrajectoryName=='rts:WEST_PCS/Plasma/zgeo/min_env_waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[3][14])
	elif (sigTrajectoryName=='rts:WEST_PCS/Plasma/rgeo/min_env_waveform.ref')\
	or (sigTrajectoryName=='rts:WEST_PCS/Plasma/rgeo/min_env_waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[3][15])
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Poloidal/IXh/min_env_waveform.ref')\
	or (sigTrajectoryName=='rts:WEST_PCS/Actuators/Poloidal/IXh/min_env_waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[3][17])		
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Poloidal/IXb/min_env_waveform.ref')\
	or (sigTrajectoryName=='rts:WEST_PCS/Actuators/Poloidal/IXb/min_env_waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[3][16])	
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Poloidal/IBh/min_env_waveform.ref')\
	or (sigTrajectoryName=='rts:WEST_PCS/Actuators/Poloidal/IBh/min_env_waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[3][25])	
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Poloidal/IDh/min_env_waveform.ref')\
	or (sigTrajectoryName=='rts:WEST_PCS/Actuators/Poloidal/IDh/min_env_waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[3][24])	
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Poloidal/IEh/min_env_waveform.ref')\
	or (sigTrajectoryName=='rts:WEST_PCS/Actuators/Poloidal/IEh/min_env_waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[3][23])	
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Poloidal/IFh/min_env_waveform.ref')\
	or (sigTrajectoryName=='rts:WEST_PCS/Actuators/Poloidal/IFh/min_env_waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[3][22])	
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Poloidal/IFb/min_env_waveform.ref')\
	or (sigTrajectoryName=='rts:WEST_PCS/Actuators/Poloidal/IFb/min_env_waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[3][21])	
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Poloidal/IEb/min_env_waveform.ref')\
	or (sigTrajectoryName=='rts:WEST_PCS/Actuators/Poloidal/IEb/min_env_waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[3][20])	
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Poloidal/IDb/min_env_waveform.ref')\
	or (sigTrajectoryName=='rts:WEST_PCS/Actuators/Poloidal/IDb/min_env_waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[3][19])
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Poloidal/IBb/min_env_waveform.ref')\
	or (sigTrajectoryName=='rts:WEST_PCS/Actuators/Poloidal/IBb/min_env_waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[3][18])
	elif (sigTrajectoryName=='rts:WEST_PCS/Actuators/Poloidal/IA/min_env_waveform.ref')\
	or (sigTrajectoryName=='rts:WEST_PCS/Actuators/Poloidal/IA/min_env_waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[3][26])	
	
	
	elif (sigTrajectoryName=='rts\WEST_PCS/Actuators/Gas/envelope4/min_env_waveform.ref')\
	or (sigTrajectoryName=='rts\WEST_PCS/Actuators/Gas/envelope4/max_env_waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[4][0])
	elif (sigTrajectoryName=='rts\WEST_PCS/Actuators/Gas/envelope3/min_env_waveform.ref')\
	or (sigTrajectoryName=='rts\WEST_PCS/Actuators/Gas/envelope3/max_env_waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[4][1])
	elif (sigTrajectoryName=='rts\WEST_PCS/Actuators/Gas/envelope2/min_env_waveform.ref')\
	or (sigTrajectoryName=='rts\WEST_PCS/Actuators/Gas/envelope2/max_env_waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[4][2])
	elif (sigTrajectoryName=='rts\WEST_PCS/Actuators/Gas/envelope1/min_env_waveform.ref')\
	or (sigTrajectoryName=='rts\WEST_PCS/Actuators/Gas/envelope1/max_env_waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[4][3])	
	elif (sigTrajectoryName=='rts\WEST_PCS/Actuators/Gas/REF3/waveform.ref')\
	or (sigTrajectoryName=='rts\WEST_PCS/Actuators/Gas/REF3/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[4][4])
	elif (sigTrajectoryName=='rts\WEST_PCS/Actuators/Gas/REF2/min_env_waveform.ref')\
	or (sigTrajectoryName=='rts\WEST_PCS/Actuators/Gas/REF2/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[4][5])
	elif (sigTrajectoryName=='rts\WEST_PCS/Actuators/Gas/REF1/min_env_waveform.ref')\
	or (sigTrajectoryName=='rts\WEST_PCS/Actuators/Gas/REF1/waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[4][6])
		
	
	
	elif (sigTrajectoryName=='rts\WEST_PCS/Actuators/Heating/ICRH/phase/3/min_env_waveform.ref')\
	or (sigTrajectoryName=='rts\WEST_PCS/Actuators/Heating/ICRH/phase/3/max_env_waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[5][0])
	elif (sigTrajectoryName=='rts\WEST_PCS/Actuators/Heating/ICRH/power/3/min_env_waveform.ref')\
	or (sigTrajectoryName=='rts\WEST_PCS/Actuators/Heating/ICRH/phase/3/max_env_waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[5][1])		
	elif (sigTrajectoryName=='rts\WEST_PCS/Actuators/Heating/ICRH/phase/2/min_env_waveform.ref')\
	or (sigTrajectoryName=='rts\WEST_PCS/Actuators/Heating/ICRH/phase/3/max_env_waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[5][2])
	elif (sigTrajectoryName=='rts\WEST_PCS/Actuators/Heating/ICRH/power/2/min_env_waveform.ref')\
	or (sigTrajectoryName=='rts\WEST_PCS/Actuators/Heating/ICRH/phase/3/max_env_waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[5][3])	
	elif (sigTrajectoryName=='rts\WEST_PCS/Actuators/Heating/ICRH/phase/1/min_env_waveform.ref')\
	or (sigTrajectoryName=='rts\WEST_PCS/Actuators/Heating/ICRH/phase/3/max_env_waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[5][4])
	elif (sigTrajectoryName=='rts\WEST_PCS/Actuators/Heating/ICRH/power/1/min_env_waveform.ref')\
	or (sigTrajectoryName=='rts\WEST_PCS/Actuators/Heating/ICRH/phase/3/max_env_waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[5][5])		
	elif (sigTrajectoryName=='rts\WEST_PCS/Actuators/Heating/LHCD/phase/2/min_env_waveform.ref')\
	or (sigTrajectoryName=='rts\WEST_PCS/Actuators/Heating/ICRH/phase/3/max_env_waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[5][6])
	elif (sigTrajectoryName=='rts\WEST_PCS/Actuators/Heating/LHCD/power/2/min_env_waveform.ref')\
	or (sigTrajectoryName=='rts\WEST_PCS/Actuators/Heating/ICRH/phase/3/max_env_waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[5][7])			
	elif (sigTrajectoryName=='rts\WEST_PCS/Actuators/Heating/LHCD/phase/1/min_env_waveform.ref')\
	or (sigTrajectoryName=='rts\WEST_PCS/Actuators/Heating/ICRH/phase/3/max_env_waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[5][8])
	elif (sigTrajectoryName=='rts\WEST_PCS/Actuators/Heating/LHCD/power/1/min_env_waveform.ref')\
	or (sigTrajectoryName=='rts\WEST_PCS/Actuators/Heating/ICRH/phase/3/max_env_waveform.ref'):
		AbsRelBitValue = int(TypeBit_SignalList[5][9])	
	
	else:
		AbsRelBitValue = 0				
				
	return(AbsRelBitValue)
										
	

