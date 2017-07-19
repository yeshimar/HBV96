#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This file consists of all tests dedicated to the implentetion of the HBV96 model
import unittest

class AboutTests(unittest.TestCase):
	"""docstring for AboutTests"""

	def test_nested_functions_and_scopes(self):
	# Test how private statement works in nested functions
		x = {'a':1, 'b':2}
		def wraper():
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
		_result = {'a': 1, 'b': 2, 'bcd': 24, 'bce': 30, 'c': 3, 'd': 4, 'e': 5}
		'''
		Test result is a big BRAVOS !

		In [8]: wr = wraper(x)
		In [9]: wr()()
		Out[9]: {'a': 1, 'b': 2, 'bcd': 24, 'bce': 30, 'c': 3, 'd': 4, 'e': 5}

		'''
		self.assertEqual(_result, wraper()()())