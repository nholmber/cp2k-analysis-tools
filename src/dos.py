#!/usr/bin/env python3
# -*- coding: <utf8> -*-

import numpy as np
import functools
from math import sqrt, pi

# Constants
HARTREE_TO_EV=27.2113838565563
EPS = 1.0e-1
THRESHOLD = 0.45
HOMO_THRESHOLD = 0.1
VACUUM_THRESHOLD = 1e-3
ENERGY_WINDOW = 1.0
ABSOLUTE_POTENTIAL = -4.44 # Volts
EPS_SMALL = 1e-6

def check_array_type(len_array=2, my_shape=-1, num_args=1):
	"""Decorator that checks that all the arguments of a function are valid numeric numpy arrays."""
	def check_array_type_decorator(func):
		@functools.wraps(func)
		def wrapper(*args, **kwargs):
			if len(args) != num_args+1:
				raise TypeError("Function should be called with exactly {num_args} argument(s).")
			for array in args[1:]:
				if not isinstance(array, np.ndarray) or len(array.shape) != len_array:
					raise TypeError(f"Input array to function should be a {len_array}-dimensional numpy array.")
				if my_shape > 0 and array.shape[-1] != my_shape:
					raise TypeError(f"Input array to function should be a {len_array}-dimensional, {my_shape}-column numpy array.")
				if array.dtype.kind not in {'u', 'i', 'f'}:
					raise TypeError(f"Input array to function should a numeric numpy array.")
			func(*args, *kwargs)
		return wrapper
	return check_array_type_decorator


def is_number(num_args=1):
	"""Decorator that checks that all the arguments of a function are valid numbers."""
	def is_number_decorator(func):
		@functools.wraps(func)
		def wrapper(*args, **kwargs):
			if len(args) != num_args+1:
				raise TypeError(f"Function should be called with exactly {num_args} argument(s).")
			for value in args[1:]:
				if not isinstance(value, (float, int)):
					raise TypeError("Input should be a float or an int.")
			func(*args, *kwargs)
		return wrapper
	return is_number_decorator


class DOSFrame:
	"""Class to store density of states (DOS) and electrode potential
	information from a single frame of a CP2K molecular dynamics simulation."""
	def __init__(self, dos, fermi):
		self.dos = dos 						# Projected density of states
		self.fermi = fermi 					# Fermi energy (i.e. HOMO energy with OT)
		self._potential = None				# Electrostatic potential averaged along surface normal vector (V)
		self._vacuum = None					# Vacuum potential (V)
		self._voltage = None				# Electrode potential [absolute, vs SHE] (V)
		self._homo = None					# Energy of the HOMO
		self._lumo = None					# Energy of the LUMO
		self._smear = None					# Gaussian smeared DOS
		self._tpdos = None					# Total projected DOS
		self._eigenvalues = None			# List of eigenvalues (the energy range of a smeared DOS)

	@property
	def dos(self):
	    return self._dos

	@dos.setter
	@check_array_type()
	def dos(self, array):
		self._dos = array

	@property
	def eigenvalues(self):
		return self._eigenvalues

	@eigenvalues.setter
	@check_array_type(len_array=1)
	def eigenvalues(self, value):
		self._eigenvalues = value

	@property
	def fermi(self):
		return self._fermi

	@fermi.setter
	@is_number()
	def fermi(self, fermi):
		if hasattr(self, "_fermi"):
			offset = self.fermi - fermi
			self.dos[:,0] -= offset
		else:
			self.dos[:,0] -= fermi

		self._fermi = fermi

		if hasattr(self, "_voltage"):
		 	self.voltage = np.array([self.vacuum[1]-self.fermi, self.vacuum[1]-self.fermi+ABSOLUTE_POTENTIAL])

	@property
	def potential(self):
	    return self._potential

	@potential.setter
	@check_array_type(my_shape=2)
	def potential(self, array):
		self._potential = array

	@property
	def homo(self):
		return self._homo

	@homo.setter
	@is_number()
	def homo(self, value):
		self._homo = value

	@property
	def lumo(self):
	    return self._lumo

	@lumo.setter
	@is_number()
	def lumo(self, value):
		self._lumo = value

	def calculate_homo_lumo(self):
		"""Calculates the HOMO and LUMO levels of a given DOS."""
		# Determine LUMO
		candidates, indices = [], []
		if self.tpdos is None:
			raise ValueError("TPDOS must be calculated prior to evaluating HOMO/LUMO.")

		for i,x in enumerate(self.tpdos):
			if candidates != [] and candidates[0] + ENERGY_WINDOW < self.dos[i,0]: break
			if self.dos[i,1] == 0.0 and x > THRESHOLD:
				candidates.append(x)
				indices.append(i)

		if candidates == []: raise ValueError('Unable to determine LUMO')
		maxval = indices[candidates.index(max(candidates))]
		self.lumo = self.dos[maxval,0]
		self._lumo_occupation = self.tpdos[maxval]

		# Determine HOMO
		candidates, indices = [], []
		i = len(self.tpdos)-1
		for x in self.tpdos[::-1]:
			if candidates != [] and candidates[0] - ENERGY_WINDOW > self.dos[i,0]: break
			if self.dos[i,1] > HOMO_THRESHOLD and x > THRESHOLD:
				candidates.append(x)
				indices.append(i)
			i -= 1

		if candidates == []: raise ValueError('Unable to determine HOMO')
		maxval = indices[candidates.index(max(candidates))]
		self.homo = self.dos[maxval,0]
		self._homo_occupation = self.tpdos[maxval]

	@property
	def smear(self):
	    return self._smear

	@smear.setter
	@check_array_type(my_shape=2)
	def smear(self, array):
		self._smear = array

	def smear_dos(self, width):
		"""Smears the total projected DOS with a Gaussian function.

		Parameters
		----------
		width:
			width of Gaussian to use for smearing.
		"""
		if not hasattr(self, "_width") or self._width != width:
			emin, emax = np.amin(self.dos[:, 0]), np.amax(self.dos[:, 0])
			self._dos_smeared = np.zeros(len(self.dos[:,0]))
			self.eigenvalues = np.linspace(emin, emax, len(self.dos[:,0]))
			if not self.tpdos is not None:
				raise ValueError("TPDOS must be calculated prior to smearing")
			for e, pd in self.tpdos:
				self._dos_smeared += pd*delta(self.eigenvalues, e, width)

			self.smear = np.hstack((self.eigenvalues.reshape((len(self._dos[:,0]),1)), self._dos_smeared.reshape((len(self._dos[:,0]),1))))
			self._width = width

	@property
	def voltage(self):
	    return self._voltage

	@voltage.setter
	@check_array_type(len_array=1)
	def voltage(self, value):
		self._voltage = value

	def get_voltage(self, relative=False):
		return self.voltage[1] if relative else self.voltage[0]

	def calculate_voltage(self):
		if self.potential is None:
			raise ValueError("The electrostatic potential must be defined before calculating the voltage.")

		if self.vacuum is None:
			raise ValueError("The vacuum potential must be defined before calculating the voltage.")

		self.voltage =  np.array([self.vacuum[1]-self.fermi, self.vacuum[1]-self.fermi+ABSOLUTE_POTENTIAL])

	@property
	def vacuum(self):
	    # Returns a two element list [position, potential]
	    return self._vacuum

	@vacuum.setter
	@check_array_type(len_array=1)
	def vacuum(self, value):
		self._vacuum = value

	@is_number()
	def calculate_vacuum(self, x):
		"""Compute value of vacuum potential which acts a reference point to calculate the electrode potential.

		Parameters
		----------
		x:
			Point where to set vacuum reference. The reference point is automatically detected
			if x is a negative number.
		"""
		# Negative X means automatic detection of vacuum potential
		if self.potential is None:
			raise ValueError("Electrostatic potential must be defined before setting the vacuum potential")
		try:
			x = float(x)
		except TypeError:
			raise ValueError(f'Input coordinate "{x}" could not be recognized as a valid float.')

		if x < 0:
			loc = self._find_vacuum_potential()
			self.vacuum = self.potential[loc, :]
		else:
			loc = np.where( abs(self.potential[:,0]-x) < EPS )
			try:
				if not loc[0][0]: raise ValueError
			except (IndexError, ValueError):
				raise ValueError(f'Given coordinate {x} outside the range (within EPS {EPS}) of allowed coordinates. Unable to determine vacuum potential')

			self.vacuum = self.potential[loc[0][0]]

		self.calculate_voltage()

	def _find_vacuum_potential(self):
		'''
		Determine vacuum potential by:
			1) Numerically differentiating electrostatic potential
			2) Find the continuous vacuum region where the electrostatic potential
			   remains roughly constant (i.e. derivative small)
			3) Find the discontinuity due to dipole correction
			4) Set vacuum potential to
				* [vacuum_start+x_max]/2 if dipole correction plane beneath surface
				* [vacuum_start+dipole_plane]/2 if dipole correction plane above surface
		'''
		# Compute finite difference derivative of electrostatic potential
		diff = np.diff(self.potential[:,1])/(self.potential[1,0]-self.potential[0,0])

		# Find where vacuum region starts/ends
		xright, xleft = 5, 5
		# List of grid points where derivative is below VACUUM_THRESHOLD
		vacuum = np.where(abs(diff)<VACUUM_THRESHOLD)[0]
		start = None
		end = None
		# Iterate over vacuum grid points to find vacuum boundaries (likely to extend over the cell due to PBC)
		for i in vacuum:
			# Skip to the next grid point if this grid point is not actually within a vacuum region, i.e.,
			# the neighboring (xleft-1)/(xright-1) grid points to the left/right of the grid point
			# are not included in the vacuum list
			if all([j not in vacuum for j in range(i+1, i+xright)]) \
			   and all([j not in vacuum for j in range(i-1, i-(xleft+1), -1)]): continue

			# Left side boundary found when the (xright-1) grid points to right of i are not within the vacuum list
			if all([j not in vacuum for j in range(i+1, i+xright)]) and i+xright < len(self.potential[:,0]):
				if end is None:
					end = i
				else:
					raise ValueError("Could not locate vacuum end boundary.")

			# Right side boundary found when the (xleft-1) grid points to left of i are not within the vacuum list
			if all([j not in vacuum for j in range(i-1, i-(xleft+1), -1)]) and i-(xleft+1) > 0:
				if start is None:
					start = i
				else:
					raise ValueError("Could not locate vacuum start boundary.")

		# Find grid point in vacuum where the derivative of the potential is largest
		idx = np.argmax(abs(np.hstack((diff[start:], diff[:end+1]))))
		# Use this point as the dipole plane
		dipole_plane = np.where(abs(diff-np.hstack((diff[start:], diff[:end+1]))[idx])<EPS_SMALL)[0][0]
		if dipole_plane > start:
			return int((dipole_plane+start)/2)
		else:
			return int((start+len(diff)+1)/2)

		if start is None or end is None: raise ValueError("Could not locate vacuum region.")

	@property
	def tpdos(self):
	    return self._tpdos

	@tpdos.setter
	@check_array_type(my_shape=2)
	def tpdos(self, value):
		self._tpdos = value

	def calculate_tpdos(self):
		"""Calculate total projected DOS as a sum over the DOS."""
		self.tpdos = np.vstack((self.dos[:,0], self.dos[:,2:].sum(axis=1))).T


class DOS(object):
	"""Container for multiple DOSFrame instances i.e. the results from a molecular dynamics simulation"""
	def __init__(self):
		self.dos = None

	@property
	def dos(self):
		return self._dos

	@dos.setter
	def dos(self, value):
		self._dos = value

	def read_dos(self, fname, t):
		"""Reads a CP2K-formatted DOS output file and returns a dictionary where each
		element is DOSFrame type indexed by the MD time step.

		Parameters
		----------
		fname:
			Name of the DOS output file. Multiple frames should be concatenated to a single file.

		t:
			A two element list that determines the initial timestep and the stride of the MD simulation.

		"""
		self.dos = {frame: DOSFrame(dos, fermi) for (frame, dos, fermi) in yield_dos(fname, t)}

	def read_potential(self, fname):
		"""Reads a file that contains the averaged electrostatic potential along
		the surface normal direction. The potential is added to the elements of
		the dictionary that contains the DOSFrame objects.

		Parameters
		----------
		fname:
			Name of the averaged electrostatic potential file. Multiple frames should be concatenated to a single file.
			The electrostatic potential can be averaged with the CP2K utility cubecruncher, which takes as
			input a cube file containing the potential (Keyword V_HARTREE_CUBE in the manual).

		"""
		frames = list(self.dos)
		frames.reverse()
		nframes = 0
		try:
			for potential in yield_potential(fname):
				nframes += 1
				self.dos[frames.pop()].potential = potential

		except IndexError:
			raise IndexError('Incompatible trajectories. Indexation different.')

		if nframes != len(self.dos.keys()):
			raise IndexError("Incompatible trajectory lengths.")

	def evaluate(self, method, frame, *args, **kwargs):
		"""Function wrapper that calls a method iteratively on all frames"""
		if not hasattr(self, '_method_list'):
			# Store a list of public functions
			self._method_list =  [func for func in dir(DOSFrame) if callable(getattr(DOSFrame, func)) and not func.startswith("_")]

		# Raise error if requested method not in list of public functions
		if method not in self._method_list:
			raise ValueError(f"Function {method} is not a callable function of type DOSFrame. "
							 f"Accepted functions: {self._method_list}.")

		if frame == 'all':
			# Iteratively call function on all frames
			return_val = []
			for frame, obj in self.dos.items():
				func = getattr(obj, method)
				my_return = func(*args, *kwargs)
				# Collect function return values in case they differ from None
				if my_return is not None: return_val.append(my_return)

			# Return non-None return vals
			if return_val != []: return return_val

		elif isinstance(frame, (float, int)):
			# Single frame
			try:
				frame = float(frame)
			except ValueError:
				pass
			obj = self.dos[frame]
			func = getattr(obj, method)
			my_return = func(*args, *kwargs)

			if my_return is not None: return my_return

		elif isinstance(frame, list):
			# Custom list of frames
			return_val = []
			for fr in frame():
				try:
					fr = float(fr)
				except ValueError:
					raise ValueError(f"Requested frame {fr} was not recognized as a valid float.")

				obj = self.dos[fr]
				func = getattr(obj, method)
				my_return = func(*args, *kwargs)
				if my_return is not None: return_val.append(my_return)

			if return_val != []: return return_val

		else:
			raise ValueError(f"Requested frame {frame} could not be parsed. Expected string 'all', a float, or a list of floats.")


	def get_property(self, property, frame):
		"""Get the value of a property of a DOSFrame object at timestep frame.

		Parameters
		----------
		property:
			The property that should be returned.

		frame:
			The timestep where to get the property. Accepts a single number,
			a list of number, or 'all' which will return the value of the
			property over the entire trajectory.
		"""
		if frame == 'all':
			# Iteratively call get_property_frame on all frames
			result = []
			for frame, obj in self.dos.items():
				result.append(self._get_property_frame(property, frame))

			return result

		elif isinstance(frame, (float, int)):
			# Single frame
			frame = float(frame)
			return self._get_property_frame(property, frame)

		elif isinstance(frame, list):
			# A custom list of frames
			result = []
			for fr in frame():
				try:
					fr = float(fr)
				except ValueError:
					raise ValueError(f"Requested frame {fr} was not recognized as a valid float.")

				result.append(self._get_property_frame(property, fr))

			return result

		else:
			raise ValueError(f"Requested frame {frame} could not be parsed. Expected string 'all', a float, or a list of floats.")


	def _get_property_frame(self, property, frame):
		"""Get attribute 'property' from DOSFrame at given timestep 'frame'."""
		if frame not in self.dos.keys():
			raise ValueError(f"Tried to get property {property} from invalid frame {frame}.")

		if not hasattr(self, '._variable_list'):
			 self._variable_list =  [var for var in dir(DOSFrame) if not callable(getattr(DOSFrame, var)) and not var.startswith("__")]

		if property not in self._variable_list:
			raise ValueError(f"Tried to get invalid property {property}. Available properties {self._variable_list}.")

		return getattr(self.dos[frame], property)


	def set_property(self, property, value, frame):
		"""A setter that modifies the value of a DOSFrame property at timestep frame.
		Parameters
		----------
		property:
			The property that should be modified.

		value:
			The new value of the DOSFrame property.

		frame:
			The timestep where to get the property. Accepts a single number,
			a list of number, or 'all' which will return the value of the
			property over the entire trajectory.
		"""
		if frame == 'all':
			# Iteratively call set_property_frame on all frames
			for frame, obj in self.dos.items():
				self._set_property_frame(property, value, frame)

		elif isinstance(frame, (float, int)):
			# Single frame
			frame = float(frame)
			self._set_property_frame(property, value, frame)

		elif isinstance(frame, list):
			# Custom list of frames
			for fr in frame():
				try:
					fr = float(fr)
				except ValueError:
					raise ValueError(f"Requested frame {fr} was not recognized as a valid float.")

				self._set_property_frame(property, value, fr)

		else:
			raise ValueError(f"Requested frame {frame} could not be parsed. Expected string 'all', a float, or a list of floats.")


	def _set_property_frame(self, property, value, frame):
		"""Changes the value of a DOSFrame property at timestep frame."""
		if frame not in self.dos.keys():
			raise ValueError(f"Tried to set property {property} to invalid frame {frame}.")

		if not hasattr(self, '._variable_list'):
			 self._variable_list =  [var for var in dir(DOSFrame) if not callable(getattr(DOSFrame, var)) and not var.startswith("__")]

		if property not in self._variable_list:
			raise ValueError(f"Tried to set invalid property {property}. Accepted properties {self._variable_list}.")

		setattr(self.dos[frame], property, value)


def yield_potential(fname):
	with open(fname, 'r') as fp:
		# Get number of lines per frame
		fp.readline()
		line = fp.readline()
		nlines = 1
		while '#' not in line:
			nlines += 1
			line = fp.readline()
		fp.seek(0)
		counter = nlines
		# Yield frames
		line = fp.readline()
		counter -= 1
		reset = True
		while line != '':
			if reset:
				pot = []
				reset = False

			if counter < nlines - 1:
				pot.append(line.strip().split())
				if counter == 0:
					counter = nlines
					reset = True
					yield (np.array(pot, dtype=np.float64))

			line = fp.readline()
			counter -= 1

def yield_dos(fname, t):
	try:
		t0 = float(t[0])
		dt = float(t[1])
	except (ValueError, TypeError):
		raise TypeError("Input list to function should contain the initial MD timestep and data sampling stride.")

	with open(fname, 'r') as fp:
		# Get number of MOs per frame
		fp.readline()
		fp.readline()
		line = fp.readline()
		while '#' not in line:
			nmo = int(line.strip().split()[0])
			line = fp.readline()
		nlines = nmo + 2
		fp.seek(0)
		t = t0 - dt
		counter = nlines
		# Yield frames
		line = fp.readline()
		counter -= 1
		reset = True
		while line != '':
			if reset:
				mo = []
				t += dt
				reset = False

			if counter < nmo:

				mo.append(line.strip().split()[1:])
				if counter == 0:
					counter = nlines
					reset = True
					# Convert CP2K output from  Hartree to eV
					dos = np.array(mo, dtype=np.float64)
					dos[:, 0] *= HARTREE_TO_EV
					yield (t, dos, fermi)

			elif counter == nlines - 1:
				fermi =  HARTREE_TO_EV*float(line.strip().split()[-2])


			line = fp.readline()
			counter -= 1

def delta(energies, energy, width):
		""" Return a delta-function centered at energy

		Parameters
		----------
		energies:
			List of energies where DOS is defined
		energy: float
			energy where the gaussian is centered
		width: float
			dispersion parameter

		Return
		------
		delta: numpy array
			array of delta function values
		"""
		x = -0.5*((energies - energy) / width)**2
		return np.exp(x) / (sqrt(2.0*pi) * width)
