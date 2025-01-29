import os, functools, struct
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances, rms
from FoldingAnalysis.clstools import *
from FoldingAnalysis.utilities import *
import json

with open("params.json", 'r') as p:
    json_params = json.load(p)

sysName = json_params['parameters'][0]['sysName']

defaults = dict(native = f'{sysName}.pdb',
				backend='OpenMP')
print("Here")

class Structure:
	"""single-frame structure
	
	Example:
		load the default structure and lazy-compute the distance matrix
		>>> dmap = Structure().dmap
		>>> type(dmap)
		<class 'numpy.ndarray'>
	"""
	def __init__(self, filename_or_universe=None):
		self.u = filename_or_universe
		if not isinstance(filename_or_universe, mda.Universe):
			self.u = mda.Universe(self.u or defaults['native'])
	
	def __getattr__(self, name):
		try:
			return getattr(self.u, name)
		except AttributeError:
			return getattr(self.u.atoms, name)
	
	@lazy_property
	def darr(self):
		"""amino acids distance array with 3 nearest neighbours ignored"""
		return darr(self.u.atoms, selection='all', ignore=3)
	
	@lazy_property
	def sarr(self):
		"""
		(corresponds to cmap in em2cmp.py and FoldingAnalysis):
		distance array of all atoms with 35 neighbours skip 
		and sigmoid_squared function applied
		"""
		return sigmoid_squared(darr(self.u.atoms, selection='all', ignore=35))
	
	def save_sarr(self, filename):
		"""save sarr to binary format for usage with GROMACS ratchet md"""
		save_sarr(self.sarr, filename, n_skip=35, cutoff=7.5)
	
	@lazy_property
	def dmap(self):
		"""amino acids distance matrix"""
		return dmap(self.u.atoms, selection='all', ignore=0)
	
	@lazy_property
	def carr(self):
		"""
		contacts array with 7.5A distance threshold and 3 nearest neighbours ignored
		"""
		return carr(self.u.atoms, selection='all', ignore=3, cutoff=7.5)
	
	@lazy_property
	def carr_soft(self):
		return sigmoid_squared(self.darr)
	
	@lazy_property
	def cmap(self):
		"""protein contact map with 7.5A distance threshold: symmetric square matrix"""
		return cmap(self.u.atoms, selection='all', ignore=0, cutoff=7.5)
	
	@lazy_property
	def cmap_soft(self):
		return sigmoid_squared(self.dmap)


class Frame:
	"""
	helper class for Trajectory to handle the __getitem__ referencing
	"""
	def __init__(self, traj, i):
		self.traj = traj
		self.i = i
	def __getattr__(self, name):
		try:
			return getattr(self.traj, name)[self.i]
		except (AttributeError, TypeError):
			try:
				return getattr(self.traj.trajectory[self.i], name)
			except (AttributeError, TypeError):
				return getattr(self.traj, name)
	def __getitem__(self, i):
		return Frame(self, i)

def _move_to_frame(method):
	@functools.wraps(method)
	def decorated(self, frame):
		self.trajectory[frame]
		return method(Frame(self, frame))
	return decorated

class _call_from_frame(MethodDecorator):
	def __getitem__(self, frame):
		frame = Frame(self._obj, frame)
		return lambda *a, **kw: self._method(frame, *a, **kw)

class TrajectoryFromStructure(type):	# metaclass
	def __new__(cls, clsname, parents, attrs):
		for parent in parents:
			# we exclude parents which were created with this metaclass:
			if parent.__class__ != TrajectoryFromStructure:
				# inheriting and re-decorating already defined lazy methods from parent:
				for attr_name in vars(parent):
					if not attr_name.startswith('_'):	# important because we also inherit from Frame!
						method = getattr(parent, attr_name)
						try:
							method = lazy_array(_move_to_frame(method._method))
						except AttributeError:
							method = _call_from_frame(method)
						# vars()[attr_name] = method
						attrs[attr_name] = method
		return super().__new__(cls, clsname, parents, attrs)

class Trajectory(Frame, Structure, metaclass=TrajectoryFromStructure):
	"""multi-frame trajectory of a structure
	
	Examples:
		Create single-frame default trajectory and compute fraction of native contacts (Q):
		>>> t = Trajectory()
		>>> type(t.q)  # Q not computed yet
		<class 'clstools._LazyArray'>
		>>> t.q[0]     # now Q for frame #0 is computed and cached
		1.0
		>>> t[0].q is t.q[0]  # alternative way to reference a frame
		True
		>>> len(t.q)   # we've got only one frame here, so next step wouldn't dump our array
		1
		>>> t.q[:]     # the whole array of Q-s is computed and cached at this point
		array([1.])
		>>> type(t.q)  # t.q is ndarray now
		<class 'numpy.ndarray'>
		
		Structure methods work for each frame here too:
		>>> type(Trajectory()[0].dmap)
		<class 'numpy.ndarray'>
	"""
	def __init__(self, 
			  filename=None, 
			  ref_filename=None,
			  bias_properties=None):
		
		ref_filename = ref_filename or defaults['native'] 
		filename = filename or ref_filename
		self.u = mda.Universe(ref_filename, filename)
		self.r = mda.Universe(ref_filename)
		self.ref = Structure(self.r)
		self.ref_filename = ref_filename
		self.trajectory = self.u.trajectory
		self.filename = filename
		self.bias_properties = bias_properties
		


		self.folding_time = None
		self.folding_frame = None

	def __getattr__(self, name):
		try:
			return getattr(self.ref, name)
		except AttributeError:
			try:
				arr = np.array([getattr(frame, name) for frame in self.trajectory])
			except AttributeError:
				# raising custom Exception instead of AttributeError (see https://bugs.python.org/issue24983):
				raise Exception(f"{type(self).__name__} object has no attribute '{name}'")
			setattr(self, name, arr)
			return arr
	
	def __len__(self):
		return len(self.trajectory)
	def __iter__(self):
		return (self[i] for i in range(len(self)))

	@lazy_array(dumped=True)
	@_move_to_frame
	def rmsd(self):
		"""protein root mean square distance (RMSD) between all atoms of a frame and the reference"""
		return rmsd(self.r.atoms, self.u.atoms)

	@lazy_array(dumped=True)
	def q(self, frame):
		"""fraction of native contacts between the C-alpha atoms (see `carr` for the definition of a contact)"""
		return self.carr[frame][self.ref.carr].sum() / self.ref.carr.sum()
	
	@lazy_array(dumped=True)
	def q_soft(self, frame):
		return self.carr_soft[frame][self.ref.carr].sum() / self.ref.carr.sum()
	
	@lazy_property
	def folded(self):
		return self.rmsd[-1] < .3

def rmsd(atoms0, atoms1):
	return 0.1 * rms.rmsd(atoms0.positions, atoms1.positions, weights=atoms1.masses, center=True, superposition=True)

def rmsd_traj(ref_filename, traj_filename):
	reference = mda.Universe(ref_filename).atoms
	return np.array([rmsd(frame, reference) for frame in mda.Universe(ref_filename, traj_filename).trajectory])

def rmsd_traj_selection(traj, ref, sel):
	"""RMSD between selecion of atoms
	sel: selecion of atoms in MDAnalysys format"""
	t = Trajectory(traj, ref)
	native = t.r.select_atoms(sel)
	return np.array([rmsd(native, t[f].u.select_atoms(sel)) for f in range(len(t[:].q))])

def ratchet_to_total_force(run):
	"""Ratio between ratchet and total force,
	provide a ratchet.out file"""
	data = np.loadtxt(run)
	ratchet_force = data[:, 4]
	total_force = data[:, 5]
	r_to_t = ratchet_force[1:]/total_force[1:]
	return r_to_t

def get_folding_frame(traj, method='simple', threshold=.3):
	"""
	Returns the frame at which the protein folds
	"""
	return computeFoldingFrame(traj.rmsd[:], method, threshold)

def get_folding_time(traj, dt=None, method='simple', threshold=.3):
	"""
	returns te time at which the protein folds
	in ps from the begnning of the simulation
	"""
	if dt == None:
		dt = traj.dt
	return computeFoldingFrame(traj.rmsd[:], method, threshold) * dt

def get_penalty_at_folding(traj):
	"""
	Value of the cumulatire ratchet force at folding
	"""
	pen = traj.bias_properties.bias_penality
	timesteps = traj.bias_properties.time
	return pen[np.abs(timesteps-get_folding_time(traj)).argmin()]

def get_penalty_at_folding_norm(traj):
	"""
	Value of the cumulatire ratchet force at folding
	normalised to total force
	"""
	pen = traj.bias_properties.bias_penality
	timesteps = traj.bias_properties.time
	tot_forces = traj.bias_properties.cum_tot_force
	return pen[np.abs(timesteps-get_folding_time(traj)).argmin()] / tot_forces[np.abs(timesteps-get_folding_time(traj)).argmin()]


def remove_neighbours_arr(arr, n=3):
	if not n:
		return arr
	mask = np.zeros(len(arr), dtype=bool)
	i, k = 0, round((8*len(arr) + 1)**.5 + 1) // 2
	while i < len(arr):
		mask[i:i+n] = True
		k -= 1
		i += k
	return arr[~mask]

def remove_neighbours_mat_inplace(mat, n=None):
	if n is None:
		return mat
	np.fill_diagonal(mat, 0)
	for i in range(1, n+1):
		r = np.arange(mat.shape[0] - i)
		mat[r, r+i] = 0.
		mat[r+i, r] = 0.
	return mat

def sigmoid_squared(x):	# vectorized version adapted from `utilities.py`
	x = x * x
	cond0 = x <= 151.29
	cond1 = np.abs(x - 56.25) < 10e-5
	x *= ~cond1
	return cond0 * np.where(cond1, 0.6, (1 - (x/56.25)**3) / (1 - (x/56.25)**5))

def get_positions(atoms, selection):
	return atoms.center_of_mass(compound=selection) if selection == 'residues' else atoms.select_atoms(selection).positions

def darr(atoms, selection='all', ignore=3):
	positions = get_positions(atoms, selection)
	return remove_neighbours_arr(distances.self_distance_array(positions, backend=defaults['backend']), n=ignore)

def dmap(atoms, selection='all', ignore=3):
	pos = get_positions(atoms, selection)
	return remove_neighbours_mat_inplace(distances.distance_array(pos, pos, backend=defaults['backend']), n=ignore)

def carr(atoms, selection='all', ignore=3, cutoff=7.5, soft=False):
	_darr = darr(atoms, selection=selection, ignore=ignore)
	if soft:
		return sigmoid_squared(_darr)
	return _darr <= cutoff

def cmap(atoms, selection='all', ignore=3, cutoff=7.5, soft=False):
	if soft:
		return sigmoid_squared(dmap(atoms, selection=selection, ignore=ignore))
	positions = get_positions(atoms, selection)
	return remove_neighbours_mat_inplace(distances.contact_matrix(positions, cutoff=cutoff), n=ignore)

def save_sarr(sarr, filename, n_skip=35, cutoff=7.5):	# adapted from `utilities.py`
	n = round((8*len(sarr) + 1)**.5 + 1) // 2 + n_skip	# recover the square matrix dimension
	mask = sum([[i, n_skip] for i in range(1, n - n_skip)], [630])	# no idea what this is...
	mat = np.empty((n, n))	# probably should rewrite this without intermediate matrix
	mat.T[np.triu_indices(n, n_skip+1)] = sarr
	sarr = mat[np.tril_indices(n, -n_skip-1)]
	# print(sarr[[j-1 + i*(2*(n-n_skip)-j)//2 for j in range(1, n-n_skip) for i in range(j)]])	# almost...
	with open(filename, 'wb') as file:
		file.write(struct.pack('=4sI4sd', b'CMAP', n, b'ABCC', cutoff))
		file.write(np.arange(n, dtype=np.uint32).tobytes())
		file.write(struct.pack('I', len(mask)))
		file.write(np.array(mask, dtype=np.uint16).tobytes())
		file.write(struct.pack('=IIdd', 1, 2, 1, 1/np.inner(sarr, sarr)))
		file.write((sarr*65535 + 0.5).astype(np.uint16).tobytes())

def q(atoms0, atoms1, cutoff=7.5):
	ref_carr = carr(atoms0, cutoff=cutoff)
	positions = atoms1.select_atoms('all')
	darr = remove_neighbours_arr(distances.self_distance_array(positions))
	return (darr <= cutoff)[ref_carr].sum() / ref_carr.sum()
