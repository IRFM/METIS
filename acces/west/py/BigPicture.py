#! /Applications/Python-3.3.5/bin/python3.3
# /usr/bin/python 

import os,sys
from xml.dom import minidom
import numpy as np
import re
import check_result
import waveform_pack
from waveformBuilder import *
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from matplotlib import gridspec
from itertools import cycle

def BigPicture_disp(segmentTrajectory,dpFile):
	
	
	#Figure window size parameters. Note that NX adjusts the vertical size if too tall.
	fig_Xposition = 0
	fig_Yposition = 0
	fig_width = 1200
	fig_height = 1000
	

	signal_array = np.array([])		
	
	#print(dpFile)
	
	#Build the signal list together with the plot number they are supposed to be displayed in.
	#First stored in a linear list for memory management efficiency
	
	signal_list = []
	signal_list.append('rts:WEST_PCS/Plasma/Ip/waveform.ref')
	signal_list.append('1')
	signal_list.append('Ip')
	
	signal_list.append('rts:WEST_PCS/Actuators/Gas/REF1/waveform.ref')
	signal_list.append('5')
	signal_list.append('ne')
	
	signal_list.append('rts:WEST_PCS/Actuators/Gas/valve01/waveform.ref')
	signal_list.append('6')
	signal_list.append('valve 1')	
	signal_list.append('rts:WEST_PCS/Actuators/Gas/valve02/waveform.ref')
	signal_list.append('6')
	signal_list.append('valve 2')		
	signal_list.append('rts:WEST_PCS/Actuators/Gas/valve03/waveform.ref')
	signal_list.append('6')
	signal_list.append('valve 3')		
	signal_list.append('rts:WEST_PCS/Actuators/Gas/valve04/waveform.ref')
	signal_list.append('6')
	signal_list.append('valve 4')		
	signal_list.append('rts:WEST_PCS/Actuators/Gas/valve05/waveform.ref')
	signal_list.append('6')
	signal_list.append('valve 5')		
	signal_list.append('rts:WEST_PCS/Actuators/Gas/valve06/waveform.ref')
	signal_list.append('6')
	signal_list.append('valve 6')		
	signal_list.append('rts:WEST_PCS/Actuators/Gas/valve07/waveform.ref')
	signal_list.append('6')
	signal_list.append('valve 7')		
	signal_list.append('rts:WEST_PCS/Actuators/Gas/valve08/waveform.ref')
	signal_list.append('6')
	signal_list.append('valve 8')		
	signal_list.append('rts:WEST_PCS/Actuators/Gas/valve09/waveform.ref')
	signal_list.append('6')
	signal_list.append('valve 9')		
	signal_list.append('rts:WEST_PCS/Actuators/Gas/valve10/waveform.ref')
	signal_list.append('6')
	signal_list.append('valve 10')		
	signal_list.append('rts:WEST_PCS/Actuators/Gas/valve11/waveform.ref')
	signal_list.append('6')
	signal_list.append('valve 11')		
	signal_list.append('rts:WEST_PCS/Actuators/Gas/valve12/waveform.ref')
	signal_list.append('6')
	signal_list.append('valve 12')		
	signal_list.append('rts:WEST_PCS/Actuators/Gas/valve13/waveform.ref')	
	signal_list.append('6')
	signal_list.append('valve 13')		
	signal_list.append('rts:WEST_PCS/Actuators/Gas/valve14/waveform.ref')	
	signal_list.append('6')
	signal_list.append('valve 14')		
	signal_list.append('rts:WEST_PCS/Actuators/Gas/valve15/waveform.ref')	
	signal_list.append('6')
	signal_list.append('valve 15')		
	signal_list.append('rts:WEST_PCS/Actuators/Gas/valve16/waveform.ref')	
	signal_list.append('6')
	signal_list.append('valve 16')		
	signal_list.append('rts:WEST_PCS/Actuators/Gas/valve17/waveform.ref')	
	signal_list.append('6')
	signal_list.append('valve 17')		
	signal_list.append('rts:WEST_PCS/Actuators/Gas/valve18/waveform.ref')	
	signal_list.append('6')
	signal_list.append('valve 18')		
	signal_list.append('rts:WEST_PCS/Actuators/Gas/valve19/waveform.ref')	
	signal_list.append('6')
	signal_list.append('valve 19')		
	signal_list.append('rts:WEST_PCS/Actuators/Gas/valve20/waveform.ref')	
	signal_list.append('6')
	signal_list.append('valve 20')		
	signal_list.append('rts:WEST_PCS/Actuators/Gas/valve21/waveform.ref')
	signal_list.append('6')
	signal_list.append('valve 21')	
			
	signal_list.append('rts:WEST_PCS/Actuators/Poloidal/IA/waveform.ref')
	signal_list.append('3')
	signal_list.append('IA')	
	signal_list.append('rts:WEST_PCS/Actuators/Poloidal/IBh/waveform.ref')
	signal_list.append('3')		
	signal_list.append('IBh')	
	signal_list.append('rts:WEST_PCS/Actuators/Poloidal/IDh/waveform.ref')
	signal_list.append('3')		
	signal_list.append('IDh')	
	signal_list.append('rts:WEST_PCS/Actuators/Poloidal/IEh/waveform.ref')
	signal_list.append('3')		
	signal_list.append('IEh')	
	signal_list.append('rts:WEST_PCS/Actuators/Poloidal/IFh/waveform.ref')
	signal_list.append('3')		
	signal_list.append('IFh')	
	signal_list.append('rts:WEST_PCS/Actuators/Poloidal/IBb/waveform.ref')
	signal_list.append('3')		
	signal_list.append('IBb')	
	signal_list.append('rts:WEST_PCS/Actuators/Poloidal/IDb/waveform.ref')
	signal_list.append('3')		
	signal_list.append('IDb')	
	signal_list.append('rts:WEST_PCS/Actuators/Poloidal/IEb/waveform.ref')
	signal_list.append('3')		
	signal_list.append('IEb')	
	signal_list.append('rts:WEST_PCS/Actuators/Poloidal/IFb/waveform.ref')
	signal_list.append('3')		
	signal_list.append('IFb')	
	
	signal_list.append('rts:WEST_PCS/Actuators/Poloidal/IXh/waveform.ref')
	signal_list.append('4')		
	signal_list.append('IXh')	
	signal_list.append('rts:WEST_PCS/Actuators/Poloidal/IXb/waveform.ref')
	signal_list.append('4')		
	signal_list.append('IXb')	
	
	signal_list.append('rts:WEST_PCS/Plasma/Rgeo/waveform.ref')
	signal_list.append('2')
	signal_list.append('Rgeo')		
	signal_list.append('rts:WEST_PCS/Plasma/Zgeo/waveform.ref')
	signal_list.append('2')	
	signal_list.append('Zgeo')		
	signal_list.append('rts:WEST_PCS/Plasma/Rext/waveform.ref')
	signal_list.append('2')	
	signal_list.append('Rext')		
	signal_list.append('rts:WEST_PCS/Plasma/dXLow/waveform.ref')
	signal_list.append('2')		
	signal_list.append('dXLow')
			
	signal_list.append('rts:WEST_PCS/Actuators/Heating/LHCD/power/1/waveform.ref')
	signal_list.append('7')
	signal_list.append('P_LH1')			
	signal_list.append('rts:WEST_PCS/Actuators/Heating/LHCD/phase/1/waveform.ref')
	signal_list.append('8')	
	signal_list.append('Phi_LH1')		
	signal_list.append('rts:WEST_PCS/Actuators/Heating/LHCD/power/2/waveform.ref')
	signal_list.append('7')	
	signal_list.append('P_LH2')		
	signal_list.append('rts:WEST_PCS/Actuators/Heating/LHCD/phase/2/waveform.ref')			
	signal_list.append('8')	
	signal_list.append('Phi_LH2')		
	signal_list.append('rts:WEST_PCS/Actuators/Heating/ICRH/power/1/waveform.ref')
	signal_list.append('7')	
	signal_list.append('P_ICRH1')		
	signal_list.append('rts:WEST_PCS/Actuators/Heating/ICRH/phase/1/waveform.ref')
	signal_list.append('8')	
	signal_list.append('Phi_ICRH1')		
	signal_list.append('rts:WEST_PCS/Actuators/Heating/ICRH/power/2/waveform.ref')
	signal_list.append('7')	
	signal_list.append('P_ICRH2')		
	signal_list.append('rts:WEST_PCS/Actuators/Heating/ICRH/phase/2/waveform.ref')
	signal_list.append('8')	
	signal_list.append('Phi_ICRH2')		
	signal_list.append('rts:WEST_PCS/Actuators/Heating/ICRH/power/3/waveform.ref')
	signal_list.append('7')	
	signal_list.append('P_ICRH3')		
	signal_list.append('rts:WEST_PCS/Actuators/Heating/ICRH/phase/3/waveform.ref')	
	signal_list.append('8')	
	signal_list.append('Phi_ICRH3')	

	signal_number = len(signal_list)/3
	
	#Reshaping of the signal list into an array.
	#signal_array(:,1) = DCS signal name
	#signal_array(:,2) = subplot to be plotted in. 
	#signal_array(:,3) = name of the waveform, to be displayed in legend
	
	signal_array = np.array([])
	signal_array = np.array(np.array(signal_list))
	signal_array = np.reshape(signal_array,(signal_number,3))
	

	#Build the waveforms from the segment trajectory, the signal names and the name of the DP.xml file
			
	wform = waveformBuilder(segmentTrajectory,signal_array[:,0],dpFile)
	

	#Array containing axes labels. 
	#Should be the same size as the number of subplots times 2 (2 axes per subplot)
	
	label_array = np.array([])
	label_array = np.append(label_array,'') #No label on patch subplot, just to keep the correct numbers for the following ones.
	label_array = np.append(label_array,'[A]')
	label_array = np.append(label_array,'[m]')
	label_array = np.append(label_array,'[A]')
	label_array = np.append(label_array,'[A]')
	label_array = np.append(label_array,'[m^-2]')
	label_array = np.append(label_array,'[Pa.m^3.s^-1]')	
	label_array = np.append(label_array,'[W]')
	label_array = np.append(label_array,'[deg]')				


	#Figure initialization
	
	fig = plt.figure()
	fig.canvas.set_window_title("The Big Picture")
	
	#Axes creation. 2 y axes per suplot using twinx.
	axarr = np.array([])
	gs = gridspec.GridSpec(5,1,height_ratios=[1,6,6,6,6])
	axarr = np.append(axarr,plt.subplot(gs[0]))
	axarr = np.append(axarr,plt.subplot(gs[1]))	
	axarr = np.append(axarr,axarr[-1].twinx())
	axarr = np.append(axarr,plt.subplot(gs[2]))	
	axarr = np.append(axarr,axarr[-1].twinx())
	axarr = np.append(axarr,plt.subplot(gs[3]))
	axarr = np.append(axarr,axarr[-1].twinx())
	axarr = np.append(axarr,plt.subplot(gs[4]))	
	axarr = np.append(axarr,axarr[-1].twinx())	
	
	
	#axarr = np.append(axarr,plt.subplot(5,1,1))
	#axarr = np.append(axarr,plt.subplot(5,1,2))	
	#axarr = np.append(axarr,axarr[-1].twinx())
	#axarr = np.append(axarr,plt.subplot(5,1,3))	
	#axarr = np.append(axarr,axarr[-1].twinx())
	#axarr = np.append(axarr,plt.subplot(5,1,4))
	#axarr = np.append(axarr,axarr[-1].twinx())
	#axarr = np.append(axarr,plt.subplot(5,1,5))	
	#axarr = np.append(axarr,axarr[-1].twinx())
	
	
	
	
	#Link the x axes to be able to zoom all suplots simultaneously
	axarr[0].get_shared_x_axes().join(axarr[0], axarr[1], axarr[3], axarr[5], axarr[7])
	
	#Reduce spaces between subplots to compactify the display
	fig.subplots_adjust(left=0.12,bottom=0.02,right=0.90,top=0.98,wspace=0.15,hspace=0.10)
	
	#Black background for better readbility in the control room			
	for i in np.arange(9):
		axarr[i].set_axis_bgcolor('black')
		
	#Color maps for curves.
	#Default pyplot colormap has a black iteration which is invisible in a black background
	#n_colors is the number of different colors in the colormap
	
		n_colors=7
		plot_color_map_base = plt.cm.prism
		#plot_color_map = plot_color_map_base(np.linspace(0.02,1.0,n_colors))
		plot_color_map = ['#FF0000','#00FF00','#FF00FF','#0000FF','#FFFF00','#00FFFF','#669900','#cccccc','#996600',]
		axarr[i].set_color_cycle(plot_color_map)

	
	#Prepare each subplot. 
	#Could have been done with a loop, but some subplots display all-zero waveforms (heating).
	#Others do not, for compacity purposes (21 gas valves)
	#Even axis numbers: left-hand y axes
	#Odd axis numbers: right-hand y axes
	
	#Patch plot (segments)
	
	axarr[0].plot([0],[0],'x-',label=signal_array[i,2],linewidth=3,markersize=10)
	
	segments_patches = []
	for ii in (np.arange(len(segmentTrajectory))):
		
		segment_duration = float(segmentTrajectory[ii,2])
		segment_start_time = float(segmentTrajectory[ii,1])-segment_duration
		rect = mpatches.Rectangle([segment_start_time,0.0],segment_duration,1.0,ec='none')
		segments_patches.append(rect)
	#p_collec = PatchCollection(segments_patches, cmap=mpl.cm.jet, alpha=0.4)
	colors = np.linspace(0, 1, len(segments_patches))
	p_collec = PatchCollection(segments_patches,cmap=plt.cm.Dark2)
	p_collec.set_array(np.array(colors))
	
	
	axarr[0].add_collection(p_collec)
	axarr[0].set_ylim([0,1])
	
	for ii in (np.arange(len(segmentTrajectory))):
		segment_duration = float(segmentTrajectory[ii,2])
		segment_start_time = float(segmentTrajectory[ii,1])-segment_duration
		if (segmentTrajectory[ii,0]=='Init'):
			text_segment = 'Init'
		else:
			text_segment = segmentTrajectory[ii,0][7:]
		
		axarr[0].text(segment_start_time+segment_duration/2.0,0.2,text_segment,color='white',ha='center')
	axarr[0].tick_params(axis='x',bottom='off',labelbottom='off')
	axarr[0].tick_params(axis='y',left='off',labelleft='off')
				
	
	
	
	# Plasma current (subplot 1, left axis)
	for i in np.where(signal_array[:,1]=='1')[0]:
		
		#Plot signal
		axarr[1].plot(wform[i].times,wform[i].values,'x-',label=signal_array[i,2],linewidth=3,markersize=10)
		
	#Get y axis limits to find the adequate size of vertical bars for segment indicators
	ylim_inf = axarr[1].get_ylim()[0]
	ylim_sup = axarr[1].get_ylim()[1]
		
	#Set yaxis labels just beyond automatic limits to avoid having horizontal lines stuck to the top
	# or the bottom of the subplot.	
	axarr[1].set_ylim([ylim_inf-0.1*abs(ylim_sup-ylim_inf),ylim_sup+0.1*abs(ylim_sup-ylim_inf)])	
		
	#Y axis label taken from table defined above.
	axarr[1].set_ylabel(label_array[1])
		
	#Display legend near the relevant axis (left or right). To be improved.
	axarr[1].legend(loc=2,fontsize=11)
		
	#Display white vertical lines to mark segment transitions.		
	for ii in (np.arange(len(segmentTrajectory))):
		segment_time = segmentTrajectory[ii,1]
		axis_range = axarr[1].get_ylim()
		axarr[1].plot([segment_time, segment_time],axis_range,'w',linewidth=1.0)

	#Repeat for right-hand axis on first subplot.
	# R,Z (subplot 1, right axis)
	for i in np.where(signal_array[:,1]=='2')[0]:

		#Do not display non-existent waveforms (unused shape control modes).
		if len(wform[i].values!=0)>0:
		
			axarr[2].plot(wform[i].times,wform[i].values,'x--',label=signal_array[i,2],linewidth=3)
	
	axarr[2].set_ylabel(label_array[2])		
	ylim_inf = axarr[2].get_ylim()[0]
	ylim_sup = axarr[2].get_ylim()[1]
	axarr[2].set_ylim([ylim_inf-0.1*abs(ylim_sup-ylim_inf),ylim_sup+0.1*abs(ylim_sup-ylim_inf)])	
	axarr[2].legend(loc=1,fontsize=11)	
	
	for ii in (np.arange(len(segmentTrajectory))):
		segment_time = segmentTrajectory[ii,1]
		axis_range = axarr[2].get_ylim()
		axarr[2].plot([segment_time, segment_time],axis_range,'w',linewidth=0.5)
			
	
	
			
	#Repeat for next subplot.
	# Coil currents (subplot 2, left axis)			
	for i in np.where(signal_array[:,1]=='3')[0]:

		axarr[3].plot(wform[i].times,wform[i].values,'x-',label=signal_array[i,2],linewidth=3)
	
	axarr[3].set_ylabel(label_array[3])			
	ylim_inf = axarr[3].get_ylim()[0]
	ylim_sup = axarr[3].get_ylim()[1]
	axarr[3].set_ylim([ylim_inf-0.1*abs(ylim_sup-ylim_inf),ylim_sup+0.1*abs(ylim_sup-ylim_inf)])	
	axarr[3].legend(loc=2,fontsize=10)
	for ii in (np.arange(len(segmentTrajectory))):
		segment_time = segmentTrajectory[ii,1]
		axis_range = axarr[3].get_ylim()
		axarr[3].plot([segment_time, segment_time],axis_range,'w',linewidth=1.0)

	
	
	# Empty plot (subplot 2, right axis)
	for i in np.where(signal_array[:,1]=='4')[0]:

		axarr[4].plot(wform[i].times,wform[i].values,'x--',label=signal_array[i,2],linewidth=3)
	
	axarr[4].set_ylabel(label_array[4])			
	ylim_inf = axarr[4].get_ylim()[0]
	ylim_sup = axarr[4].get_ylim()[1]
	axarr[4].set_ylim([ylim_inf-0.1*abs(ylim_sup-ylim_inf),ylim_sup+0.1*abs(ylim_sup-ylim_inf)])			
	axarr[4].legend(loc=1,fontsize=11)		
		
	for ii in (np.arange(len(segmentTrajectory))):
		segment_time = segmentTrajectory[ii,1]
		axis_range = axarr[4].get_ylim()
		axarr[4].plot([segment_time, segment_time],axis_range,'w',linewidth=0.5)	
		
	#axarr[4].tick_params(axis='y',right='off',labelright='off')
	
	
	
	# Density request (subplot 3, left axis)
	for i in np.where(signal_array[:,1]=='5')[0]:
		if len(np.where(wform[i].values!=0)[0])>0:
			axarr[5].plot(wform[i].times,wform[i].values,'x-',label=signal_array[i,2],linewidth=3)
	
	axarr[5].set_ylabel(label_array[5])	
	ylim_inf = axarr[5].get_ylim()[0]
	ylim_sup = axarr[5].get_ylim()[1]
	axarr[5].set_ylim([ylim_inf-0.1*abs(ylim_sup-ylim_inf),ylim_sup+0.1*abs(ylim_sup-ylim_inf)])	
	axarr[5].legend(loc=2,fontsize=11)				
	
	for ii in (np.arange(len(segmentTrajectory))):
		segment_time = segmentTrajectory[ii,1]
		axis_range = axarr[5].get_ylim()
		axarr[5].plot([segment_time, segment_time],axis_range,'w',linewidth=0.5)
	
	
	
	# Gas valve opening requests (subplot 3, right axis)
	for i in np.where(signal_array[:,1]=='6')[0]:
		if len(np.where(wform[i].values!=0)[0])>0:
			axarr[6].plot(wform[i].times,wform[i].values,'x--',label=signal_array[i,2],linewidth=3)
	
	axarr[6].set_ylabel(label_array[6])	
	ylim_inf = axarr[6].get_ylim()[0]
	ylim_sup = axarr[6].get_ylim()[1]			
	axarr[6].set_ylim([ylim_inf-0.1*abs(ylim_sup-ylim_inf),ylim_sup+0.1*abs(ylim_sup-ylim_inf)])	
	axarr[6].legend(loc=1,fontsize=10)
		
	for ii in (np.arange(len(segmentTrajectory))):
		segment_time = segmentTrajectory[ii,1]
		axis_range = axarr[6].get_ylim()
		axarr[6].plot([segment_time, segment_time],axis_range,'w',linewidth=0.5)			
	
	
			
	# Heating powers (subplot 4, left axis)
	for i in np.where(signal_array[:,1]=='7')[0]:
		axarr[7].plot(wform[i].times,wform[i].values,'x-',label=signal_array[i,2],linewidth=3)
	
	ylim_inf = axarr[7].get_ylim()[0]
	ylim_sup = axarr[7].get_ylim()[1]
	axarr[7].set_ylabel(label_array[7])
	axarr[7].set_ylim([ylim_inf-0.1*abs(ylim_sup-ylim_inf),ylim_sup+0.1*abs(ylim_sup-ylim_inf)])	
	axarr[7].legend(loc=2,fontsize=10)			

	
	for ii in (np.arange(len(segmentTrajectory))):
		segment_time = segmentTrajectory[ii,1]
		axis_range = axarr[7].get_ylim()
		axarr[7].plot([segment_time, segment_time],axis_range,'w',linewidth=0.5)
	
	
	
	# Heating phases (subplot 4, right axis)
	for i in np.where(signal_array[:,1]=='8')[0]:
		
		axarr[8].plot(wform[i].times,wform[i].values,'x--',label=signal_array[i,2],linewidth=3)
		
	ylim_inf = axarr[8].get_ylim()[0]
	ylim_sup = axarr[8].get_ylim()[1]
	axarr[8].set_ylabel(label_array[8])
	axarr[8].set_ylim([ylim_inf-0.1*abs(ylim_sup-ylim_inf),ylim_sup+0.1*abs(ylim_sup-ylim_inf)])			
	axarr[8].legend(loc=1,fontsize=12)	
	for ii in (np.arange(len(segmentTrajectory))):
		segment_time = segmentTrajectory[ii,1]
		axis_range = axarr[8].get_ylim()
		axarr[8].plot([segment_time, segment_time],axis_range,'w',linewidth=0.5)		
						


	#Set figure geometry
	mngr = plt.get_current_fig_manager()
	mngr.window.setGeometry(fig_Xposition,fig_Yposition,fig_width,fig_height)	
	
	plt.show()
