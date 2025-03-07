#! /Applications/Python-3.3.5/bin/python3.3
# /usr/bin/python 

import os,sys
from xml.dom import minidom
import numpy as np
import re

class waveform_pack():
	"""
	Class defining a waveform data package
	Used by waveformBuilder.py to send its data.
	"""
	def __init__(self):
		self.times = np.array([])	# absolute time vector
		self.reltimes = np.array([])	# relative (=inside the segment) time vector
		self.values = np.array([])	# Value vector
		self.segments = np.array([])	# Segment number of each point 
		# self.ylabel = ''		# ylabel (units)
	
