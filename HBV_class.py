#!/usr/bin/env python
# -*- coding: utf-8 -*-
# HBV-96 model class


class HBV96(object):
	"""docstring for HBV96"""

	"""Public Static Variables"""
	# This part remains to be implenmented as a module or codeblock

	# Parameter index, 18 modifiable and 2 not
	_ind = ['ltt',
			'utt',
			'ttm',
			'cfmax',
			'fc',
			'ecorr',
			'etf',
    		'lp',
    		'k',
    		'k1',
    		'alpha',
    		'beta',
    		'cwh',
    		'cfr',
    		'c_flux',
    		'perc',
    		'rfcf',
    		'sfcf'
    		'tfac', #Non-modifiable, nm of hrs in time step
    		'area'	#Non-modifiable, catchment area
    		]

	# Lower boundary parameters
	P_LB = [-1.5, 		#ltt
	        0.001, 		#utt
	        0.001, 		#ttm
	        0.04, 		#cfmax [mm c^-1 h^-1]
	        50.0, 		#fc
	        0.6, 		#ecorr
	        0.001, 		#etf
	        0.2, 		#lp
	        0.00042, 	#k [h^-1] upper zone
	        0.0000042, 	#k1 lower zone
	        0.001, 		#alpha
	        1.0, 		#beta
	        0.001, 		#cwh
	        0.01, 		#cfr
	        0.0, 		#c_flux
	        0.001, 		#perc mm/h
	        0.6, 		#rfcf
	        0.4] 		#sfcf

	# Upper boundary parameters
	P_UB = [2.5, 		#ltt
	        3.0, 		#utt
	        2.0, 		#ttm
	        0.4, 		#cfmax [mm c^-1 h^-1]
	        500.0, 		#fc
	        1.4, 		#ecorr
	        5.0, 		#etf
	        0.5, 		#lp
	        0.0167, 	#k upper zone
	        0.00062, 	#k1 lower zone
	        1.0, 		#alpha
	        6.0, 		#beta
	        0.1, 		#cwh
	        1.0, 		#cfr
	        0.08, 		#c_flux - 2mm/day
	        0.125, 		#perc mm/hr
	        1.4, 		#rfcf
	        1.4]		#sfcf

	# Initial status
	DEF_ST = {	'sp': 0.0,	#sp: Snow pack
				'sm': 30.0,	#sm: Soil moisture
				'uz': 30.0,	#uz: Upper zone direct runoff
				'lz': 30.0,	#lz: Lower zone direct runoff
				'wc': 0.0	#wc: Water content
				}		
					
	# Boundary conditions DataFrame
	BOUND = {"low": dict(zip(_ind[:18], P_LB)),
				"up": dict(zip(_ind[:18], P_UB))}

	# Initial flow rate		
	DEF_q0 = 10.0

	# HBV96 model initializer, only Mparameter enough?
	def __init__(self):
		self._par = dict(zip(self._ind, len(self._ind)*[0.0]))

	@property
	def par(self):
		return self._par

	# Render parameter function
	def render_par(self, obj_name, par_name):
		return	eval('self.'+obj_name+'[par_name]')

	# Set parameter function, external accessibility for Django
	def set_par(self, obj_name, par_name, value):
		exec('self.'+obj_name+'[par_name] = value')

# Test how private statement works in nested functions
x = {'a':1, 'b':2}
def wraper(x):
	x['c'] = 3
	_bc = x['b']*x['c']
	def _body():
		x['d'] = 4
		x['bcd'] = _bc*x['d']
		def _guts():
			x['e'] = 5
			x['bce'] = _bc*x['e']
			return x
		return _guts
	return _body

'''
Test result is all a big BRAVOS !
In [8]: wr = wraper(x)

In [9]: wr()()
Out[9]: {'a': 1, 'b': 2, 'bcd': 24, 'bce': 30, 'c': 3, 'd': 4, 'e': 5}

'''