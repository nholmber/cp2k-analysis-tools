#!/usr/bin/env python3
# -*- coding: <utf8> -*-

import numpy as np
import re
from itertools import combinations
from scipy.interpolate import RectBivariateSpline, InterpolatedUnivariateSpline

class GridPotential:
	'''
	Manipulate a potential defined on a one or two dimensional rectangular grid.
	'''

	def __init__(self, fname, **kwargs):
		self._y = None
		self._dimension = 1
		self.extrema = None
		self.read(fname, **kwargs)

	def read(self, fname, **kwargs):
		if 'dimension' in kwargs:
			try:
				if kwargs['dimension'] == 1:
					pass
				elif kwargs['dimension'] == 2:
					self._dimension = 2
					y = []
				else:
					raise ValueError
			except ValueError:
				print("Potential dimension must be +1 or +2.")

		# Parse header lines
		with open(fname, 'r') as fp:
			line = fp.readline()
			while re.search('^\s*#', line) is not None:
				# Try looking for vector containing y grid points in the header lines
				if self._dimension == 2:
					rx = re.findall("[-+]?\d+\.\d+", line)
					if  rx != []: y.append(rx)

				line = fp.readline()

			if 'usecols' in kwargs:
				self._npot = len(kwargs['usecols'])
				usecols = kwargs['usecols']
			else:
				self._npot = len(line.split())
				usecols = tuple(i for i in range(self._npot))

			if self._dimension == 2:
				#TODO: Parse 2d potential in gnuplot format
				y = [i for j in y for i in j]
				if y == []:
					raise ValueError("Tried to interpret file as a two dimensional potential but couldn't find vector for "
				 					 "second dimension grid points. Pass the vector in the header lines that start with '#'.")
				if len(y) != self._npot-1:
					raise ValueError("The number of columns and the size of the grid point vector for the second dimension don't match. "
									 "Expected {0:d} columns, got {1:d}.".format(len(y)+1, self._npot))
				self._y = np.array(y, dtype=np.float64)

		data = np.loadtxt(fname, usecols=usecols, unpack=True)

		self._x = data[0]
		if self._dimension == 1:
			self._values = data[1:].reshape(data.shape[1])
		else:
			self._values = data[1:]

	def absolute_to_relative(self, unit=None, reference=None):
		if reference is None:
			self._values -= np.amin(self._values)
			if hasattr(self, '_valuesinterp') and self._valuesinterp is not None: self._valuesinterp -= np.amin(self._valuesinterp)
		else:
			self._values -= reference
			if hasattr(self, '_valuesinterp') and self._valuesinterp is not None: self._valuesinterp -= reference

		if unit is not None:
			try:
				assert isinstance(unit, float)
				self._values *= unit
				if hasattr(self, '_valuesinterp') and self._valuesinterp is not None: self._valuesinterp *= unit
			except AssertionError:
				print('A floating point number must be provided for unit conversion.')
				raise

	@property
	def potential(self):
		return (self._x, self._y, self._values)

	@property
	def interpolated_potential(self):
		return (self._xinterp, self._yinterp, self._valuesinterp)

	@property
	def extrema(self):
		return (self._minima, self._maxima, self._saddle)

	@extrema.setter
	def extrema(self, value):
		try:
			if value is None:
				val1 = val2 = val3 = None
			else:
				# Unpack
				val1, val2, val3 = value
		except ValueError:
			raise ValueError("Pass an iterable with three items")
		else:
			self._minima = val1
			self._maxima = val2
			self._saddle = val3

	@property
	def global_minimum(self):
		return (np.amin(self._values))

	def interpolate(self, order=3, grid=None, forced_grid=False):
		'''
		Interpolate potential using splines of the given order (default: cubic)
		and tighter grid (default: add midpoint between each original grid point).
		'''
		self._create_interpolate_grid(grid, forced_grid)

		if self._dimension == 1:
			self._spline = InterpolatedUnivariateSpline(self._x, self._values, k=order)
			self._valuesinterp = self._spline(self._xinterp)

		elif self._dimension == 2:
			self._spline = RectBivariateSpline(self._y, self._x, self._values, kx=order, ky=order)
			self._valuesinterp = self._spline(self._yinterp, self._xinterp)

	def _create_interpolate_grid(self, grid, forced_grid):
		self._xinterp = None
		self._yinterp = None
		self._valuesinterp = None

		if grid is not None :
			if forced_grid:
				if self._dimension == 1:
					self._xinterp =  grid
				else:
					self._xinterp = grid[0]
					self._yinterp = grid[1]
				return

			grid = np.array(grid, dtype=np.float64)
			try:
				assert grid.ndim == self._dimension
			except AssertionError:
				print('The number of dimensions for argument grid is wrong.')
				raise
			if self._dimension == 1:
				self._xinterp = np.concatenate((self._x, grid))
			else:
				self._xinterp = np.concatenate((self._x, grid[0]))
				self._yinterp = np.concatenate((self._y, grid[1]))

		else:
			stride = np.diff(self._x)/2
			self._xinterp = np.concatenate((self._x,  self._x+np.concatenate((stride, np.zeros((1,))))))
			if self._dimension == 2:
				stride = np.diff(self._y)/2
				self._yinterp = np.concatenate((self._y,  self._y+np.concatenate((stride, np.zeros((1,))))))

		self._xinterp.sort()
		self._xinterp = np.unique(self._xinterp)
		if self._dimension == 2:
				self._yinterp.sort()
				self._yinterp = np.unique(self._yinterp)

	def find_extrema(self, use_interpolated=False):
		'''
		Find local minima/maxima and saddle points of grid potential.
		'''
		if use_interpolated:
			data = self._valuesinterp
			grid = self._xinterp
		else:
			data = self._values

		if self._dimension == 2:
			self._find_extrema_2d(data)
		else:
			try:
				self._find_extrema_1d(data, grid=grid)
			except UnboundLocalError:
				self._find_extrema_1d(data)

		return self.extrema

	def _find_extrema_2d(self, data):
		# Find min/max along data rows
		row_maxima = set()
		row_minima = set()
		for i, row in enumerate(data):
			minima = []
			maxima = []
			for j, col in enumerate(row):
				if j == 0: continue
				try:
					if col < row[j-1] and col < row[j+1]:
						minima.append(j)
					if col > row[j-1] and col > row[j+1]:
						maxima.append(j)
				except IndexError:
					# Minima/Maxima cannot lie on data boundaries
					pass
			[row_maxima.add(tupl) for tupl in [(j, i) for j in maxima]]
			[row_minima.add(tupl) for tupl in [(j, i) for j in minima]]

		col_maxima = set()
		col_minima = set()
		# Find min/max along data columns
		for i, row in enumerate(data.T):
			minima = []
			maxima = []
			for j, col in enumerate(row):
				if j == 0: continue
				try:
					if col < row[j-1] and col < row[j+1]:
						minima.append(j)
					if col > row[j-1] and col > row[j+1]:
						maxima.append(j)
				except IndexError:
					# Minima/Maxima cannot lie on data boundaries
					pass

			[col_maxima.add(tupl) for tupl in [(i, j) for j in maxima]]
			[col_minima.add(tupl) for tupl in [(i, j) for j in minima]]

		# Construct local extrema as unions of row/col local extrema
		self._minima = row_minima & col_minima
		self._maxima = row_maxima & col_maxima
		self._saddle = row_minima & col_maxima | row_maxima & col_minima

	def _find_extrema_1d(self, data, **kwargs):
		self._minima = []
		self._maxima = []
		for j, val in enumerate(data):
			if j == 0: continue
			try:
				if val < data[j-1] and val < data[j+1]:
					self._minima.append(j)
				if val > data[j-1] and val > data[j+1]:
					self._maxima.append(j)
			except IndexError:
				# Minima/Maxima cannot lie on data boundaries
				pass

		self._saddle = None
		try:
			grid = kwargs['grid']
		except:
			grid = None

		if len(self._minima) > 1:
			self._saddle = []
			for pair in combinations(self._minima, 2):
				index = np.argmax(data[pair[0]:pair[1]])+pair[0]
				if grid is None:
					self._saddle.append(index)
				else:
					self._saddle.append(grid[index])