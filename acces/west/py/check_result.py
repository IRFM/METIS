#! /Applications/Python-3.3.5/bin/python3.3
# /usr/bin/python 

import os,sys
from xml.dom import minidom
import numpy as np
import re

#class defnition for check result data. Given as output by check_functions.py.

class check_result():
	def __init__(self):
		self.check_name	= ''			# Name of the check (for WOIs containing more than one check)	
		self.check_result_code = 0 		# 0 = error; 1 = warning; 2 = data unavailable (on/offline); 3 = OK
		self.check_result_text = '' 		# Free explanation text
		self.check_fail_abs_times = []		# Absolute times in the scenario at which the check is failed
		self.check_fail_values = []		# Value at which the check is failed
		self.check_fail_values_unit = '' 	# Physical unit of the parameter
		self.check_fail_rel_times = []		# Relative times in the segment at which the check is failed
		self.check_fail_segments = []		# Segments corresponding to the check failures
		self.check_fail_limit = 0		# Limit tested
	

