# -*- coding: utf-8 -*- 
import MDAnalysis as md
from MDAnalysis.analysis import distances
from MDAnalysis.analysis import contacts
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import RMSD
import numpy as np
import itertools 
import warnings
import os
from numba import vectorize
from numba import njit
from numba import prange
import numba
import struct
from FoldingAnalysis.BiasProperties import BiasProperties
import FoldingAnalysis as fa

class _Struct:
    pass

def hardQ(trajectory, ref, use_ref_as_last = False, min_dist = 3, cutoff = 7.5, selection = 'all'):
    
    if use_ref_as_last:
        Q = np.zeros(trajectory.getUniverse().trajectory.n_frames+1)
    else:
        Q = np.zeros(trajectory.getUniverse().trajectory.n_frames)
        
    ref_g1 = md.AtomGroup([i for [i,j] in itertools.combinations(ref.getUniverse().select_atoms(selection),2)])
    ref_g2 = md.AtomGroup([j for [i,j] in itertools.combinations(ref.getUniverse().select_atoms(selection),2)])
    res_id1, res_id2, R0 = distances.dist(ref_g1,ref_g2)
    keep_cutoff = np.abs(res_id1 - res_id2) > min_dist
    keep_native = R0 < cutoff
    R0 = R0[np.logical_and(keep_cutoff, keep_native)]
    
    #print(f'Native contacts: {len(R0)}')
    
    traj_g1 = md.AtomGroup([i for [i,j] in itertools.combinations(trajectory.getUniverse().select_atoms(selection),2)])[np.logical_and(keep_cutoff, keep_native)]
    traj_g2 = md.AtomGroup([j for [i,j] in itertools.combinations(trajectory.getUniverse().select_atoms(selection),2)])[np.logical_and(keep_cutoff, keep_native)]
    
    i = 0
    for ts in trajectory.getUniverse().trajectory:
        R = distances.dist(traj_g1,traj_g2)[2]
        #R = R[np.logical_and(keep_cutoff, keep_native)]
        #Q[i] = contacts.hard_cut_q(R, cutoff=cutoff)
        Q[i] = np.sum(R <= cutoff) / R.shape[0]
        i+=1
    
    if use_ref_as_last:
        Q[i] = 1
    
    return Q

def bestHummerQ(trajectory, ref, use_ref_as_last = False, min_dist = 3, cutoff = 4.5, beta_c = 5, lambda_c = 1.8, selection='all'):
    
    if use_ref_as_last:
        Q = np.zeros(trajectory.getUniverse().trajectory.n_frames+1)
    else:
        Q = np.zeros(trajectory.getUniverse().trajectory.n_frames)
        
    ref_g1 = md.AtomGroup([i for [i,j] in itertools.combinations(ref.getUniverse().select_atoms(selection),2)])
    ref_g2 = md.AtomGroup([j for [i,j] in itertools.combinations(ref.getUniverse().select_atoms(selection),2)])
    res_id1, res_id2, R0 = distances.dist(ref_g1,ref_g2)
    keep_cutoff = np.abs(res_id1 - res_id2) > min_dist
    keep_native = R0 < cutoff
    R0 = R0[np.logical_and(keep_cutoff, keep_native)]
    
    #print(f'Native contacts: {len(R0)}')
    
    traj_g1 = md.AtomGroup([i for [i,j] in itertools.combinations(trajectory.getUniverse().select_atoms(selection),2)])[np.logical_and(keep_cutoff, keep_native)]
    traj_g2 = md.AtomGroup([j for [i,j] in itertools.combinations(trajectory.getUniverse().select_atoms(selection),2)])[np.logical_and(keep_cutoff, keep_native)]
    
    i = 0
    for ts in trajectory.getUniverse().trajectory:
        R = distances.dist(traj_g1,traj_g2)[2]
        #R = R[np.logical_and(keep_cutoff, keep_native)]
        Q[i] = contacts.soft_cut_q(R, R0, beta=beta_c, lambda_constant=lambda_c)
        i+=1
    
    if use_ref_as_last:
        Q[i] = contacts.soft_cut_q(R0, R0, beta=beta_c, lambda_constant=lambda_c)
    
    return Q

def saveCmapOld(f_name, cmap):
    k = int((np.sqrt(1+8*len(cmap))-1)/2)
    f = open(f_name,'w')
    start = 0
    end = k
    for i in range(k,0,-1):
        f.write(' '.join(str(round(n,6) if n != 0 else 0) for n in cmap[start:end]) + '\n')
        start = end
        end = end + i -1
    f.close()

def saveCmapOldLegacy(f_name, cmap):
    k = int((np.sqrt(1+8*len(cmap))-1)/2)
    f = open(f_name,'w')
    start = 0
    end = k - 1
    for i in range(k - 1,0,-1):
        f.write(' '.join(str(round(n,6) if n != 0 else 0) for n in cmap[start:end]) + '\n')
        start = end + 1
        end = end + i
    f.close()

# def computeCmap(trajectory, start=None, end=None, use_ref_as_last = False, ref = None, min_dist = 35, 
#                       verbose = False, map_function = lambda x:x, selection='all'):
    
#     atoms = trajectory.getUniverse().select_atoms(selection)
#     n_atoms = len(atoms)
#     start_f = 0
#     end_f = trajectory.getUniverse().trajectory.n_frames
    
#     if start is not None:
#         start_f = start
#     if end is not None:
#         end_f = end
    
#     if use_ref_as_last and ref is not None:
#         n_frames = end_f - start_f + 1
#     else:
#         n_frames = end_f - start_f
    
    
#     keep_index = np.array([int(((n_atoms*(n_atoms-1))/2) - ((n_atoms-i)*((n_atoms-i)-1))/2 + j - i - 1) for i,j in \
#                   itertools.combinations(range(n_atoms),2) if abs(atoms[i].ix - atoms[j].ix) > min_dist], dtype=np.uint32)

#     trajectory.getUniverse().trajectory[start_f]
#     cmap_traj = np.empty((n_frames, len(map_function(distances.self_distance_array(atoms.positions)[keep_index]))), dtype=np.float32)
    
#     i = 0
#     for ts in trajectory.getUniverse().trajectory[start_f:end_f]:
#         if verbose: print('Cmap '+str(i+1)+' of '+str(n_frames),end='\r',flush=True)
        
#         cmap_traj[i] = map_function(distances.self_distance_array(atoms.positions)[keep_index])
#         i+=1
    
#     if use_ref_as_last and ref is not None:
#         ref_atoms = ref.getUniverse().select_atoms(selection)
#         if verbose: print('Cmap '+str(i+1)+' of '+str(n_frames),end='\r',flush=True)
#         cmap_traj[i] = map_function(distances.self_distance_array(ref_atoms.positions)[keep_index])
    
#     return cmap_traj

def computeRMSD(trajectory, ref, selection='all'):
    res = RMSD(trajectory.getUniverse(), ref.getUniverse(), weights='mass', select=selection).run()
    return res.rmsd[:,2]

def computeFoldingFrame(rmsd, method='simple',threshold = 2, tolerance = 0.5, ignore_last_frames=100):
    if method == 'simple':
        for i in range(len(rmsd)):
            if rmsd[i] < threshold and i < len(rmsd) - ignore_last_frames:
                return i
        return -1
    elif method == 'average':
        for i in range(len(rmsd)):
            if rmsd[i] < threshold and np.mean(rmsd[i:]) < threshold and i < len(rmsd) - ignore_last_frames:
                return i
        return -1
    elif method == 'median':
        for i in range(len(rmsd)):
            if rmsd[i] < threshold and np.median(rmsd[i:]) < threshold and i < len(rmsd) - ignore_last_frames:
                return i
        return -1
    elif method == 'rmse':
        # for each trajectory in the ensable calcolate the rmsd and put it in an arrray. Than calculate the error (residuals predicted - true) and fill another array
        errors = ((rmsd - threshold) * (rmsd > threshold))**2
        print("rmsd", rmsd)
        print("errors", errors)
        # now I want to selelct the fram where all the conditions are respecte. So i iterate over the number of elemnts in the array that i have 
        # created before
        for i in range(len(errors)):
            if rmsd[i] < threshold and np.sqrt(np.sum(errors[i:])/len(errors[i:])) < tolerance and i < len(rmsd) - ignore_last_frames:
                return i
        return -1
    else:
        return computeFoldingFrame(rmsd, threshold=threshold, tolerance=tolerance, ignore_last_frames=ignore_last_frames)

def computeFoldingTime(rmsd, dt, method='simple',threshold = 2, tolerance = 0.5, ignore_last_time=100):
    fold_frame = computeFoldingFrame(rmsd, method=method, threshold=threshold, tolerance=tolerance, 
                                 ignore_last_frames=int(round(ignore_last_time/dt)))
    
    return -1 if fold_frame < 0 else fold_frame*dt

def downsampleOverDistance(cmap_traj, native_cmap, distance=0.1, margin_factor=1, rev=False, return_index=False):
    native_cmap = np.asarray(native_cmap, dtype=np.float32)
    native_cmap = np.atleast_2d(native_cmap)[0]
    cmap_traj = np.asarray(cmap_traj, dtype=np.float32)
    cmap_traj = np.atleast_2d(cmap_traj)
    
    #cmaps = [native_cmap] if rev else [cmap_traj[0]]
    cmaps = [cmap_traj[-1]] if rev else [cmap_traj[0]]
    native_norm = np.linalg.norm(native_cmap)
    warning = False
    i = cmap_traj.shape[0]-1 if rev else 0
    indexes = [i]
    for cmap in (cmap_traj[::-1] if rev else cmap_traj):
        delta_norm = np.linalg.norm(cmap - cmaps[-1])/native_norm
        print('> delta norm: '+str(delta_norm))
        if delta_norm > distance:
            warning = True
            if delta_norm < distance + distance*margin_factor:
                warning = False
                cmaps.append(cmap)
                indexes.append(i)
        i+= -1 if rev else 1
     
    # if not np.allclose(cmaps[-1],native_cmap) and not rev:
    #     cmaps.append(native_cmap)
    
    if warning:
        warnings.warn('resulting trajectory is interrupted, cannot find two consecutive cmaps close enough.')

    if return_index:
        return (np.array(cmaps, dtype=np.float32)[::-1] if rev else np.array(cmaps, dtype=np.float32), np.array(indexes)[::-1] if rev else np.array(indexes))
    else:
        return np.array(cmaps, dtype=np.float32)[::-1] if rev else np.array(cmaps, dtype=np.float32)

def rmd_valsParser(filename):
    rmd_vals = np.loadtxt(filename)
    if rmd_vals.shape[1] == 7:
        return BiasProperties(rmd_vals[:,0],rmd_vals[:,2],rmd_vals[:,4],rmd_vals[:,5],z_coord=rmd_vals[:,6],
                              z_min=rmd_vals[:,1])
    if rmd_vals.shape[1] == 9:
        return BiasProperties(rmd_vals[:,0],rmd_vals[:,3],rmd_vals[:,5],rmd_vals[:,6],s_coord=rmd_vals[:,7],
                              w_coord=rmd_vals[:,8],s_min=rmd_vals[:,1],w_min=rmd_vals[:,2])
    return None

def ratchet_outParser(filename):
    ratchet_out = np.loadtxt(filename)
    return BiasProperties(time=ratchet_out[:,0],progress=ratchet_out[:,1],closeness=ratchet_out[:,2],
                          bias_penality=ratchet_out[:,4], cum_tot_force=ratchet_out[:,5])

def suggestLambda(cmaps, reference):
    norms = np.linalg.norm(cmaps[:-1] - cmaps[1:], axis=1)
    return 1 / np.mean((norms/np.linalg.norm(reference))**2)

@vectorize('float64(float64)',target='parallel')
def sigmoidParallel(x):
    if x > 12.3:
        return 0
    if np.abs(x-7.5) < 10e-5:
        return 0.6
    return (1-(x/7.5)**6)/(1-(x/7.5)**10)

def loadEnsemble(directory, kind = 'legacy', specification = None):
    directory = directory.rstrip('/')
    if kind == 'legacy':
        reference = fa.Trajectory(directory + '/em.pdb', name='reference')
        #reference.configureCmap(use_ref_as_last=False,map_function=sigmoidParallel)
        trajectories = []
        for f in os.listdir(directory):
            if f.startswith('TRIAL_TRAJ'):
                i_traj = directory + '/' + f + '/rmd.trr'
                i_rmd_vals = directory + '/' + f + '/rmd_vals'
                if not os.path.isfile(i_traj) or not os.path.isfile(i_rmd_vals):
                    continue
                bias_sproperties = rmd_valsParser(i_rmd_vals)
                trajectories.append(fa.Trajectory(directory + '/em.pdb', trajectory_file=i_traj, 
                                               name=f,bias_properties=bias_sproperties))
        ensemble = fa.TrajectoryEnsemble(reference=reference, trajectories=trajectories)
        ensemble.setDt(1)
        ensemble.configureFolding(method='rmse',ignore_last_time=100,tolerance=0.3)
        #ensemble.configureCmap(verbose=True, map_function=sigmoidParallel)
        return ensemble
    return None

def saveCmapNew(f_name, cmap_list, indexes, lambdas, colvar):
    cmap_list = np.asarray(cmap_list)
    cmap_list = np.atleast_2d(cmap_list)
    lambdas = np.asarray(lambdas)
    lambdas = np.atleast_1d(lambdas)
    colvar = np.asarray(colvar)
    colvar = np.atleast_1d(colvar)

    resort_indexes = np.lexsort((indexes[:,1], indexes[:,0]))
    indexes = indexes[resort_indexes]
    cmap_list = cmap_list[:,resort_indexes]
    
    n_frames = cmap_list.shape[0]
    ext_indexes = np.unique(np.array(indexes, dtype=np.uint32).flatten())
    
    indexes_t = np.asarray(indexes, dtype=np.uint32)[:,::-1]
    upper2lower = np.lexsort((indexes_t[:,1], indexes_t[:,0]))
    indexes_t = indexes_t[upper2lower,:]
    
    indexes_all = np.array([(j,i) for i, j in itertools.combinations(ext_indexes,2)], dtype=np.uint32)
    indexes_all = indexes_all[np.lexsort((indexes_all[:,1], indexes_all[:,0]))]
    
    cur_pos = 0
    c_0 = True
    mask = []
    n_zeros = 0
    n_ones = 0
    for i in range(len(indexes_all)):
        if c_0:
            if np.array_equal(indexes_all[i], indexes_t[cur_pos]):
                mask.append(n_zeros)
                c_0 = False
                n_zeros = 0
                n_ones += 1
            else:
                if n_zeros == 65535:
                    mask.append(n_zeros)
                    mask.append(0)
                    n_zeros = 0
                n_zeros += 1
        else:
            if cur_pos < len(indexes_t) - 1:
                cur_pos += 1
            if np.array_equal(indexes_all[i], indexes_t[cur_pos]):
                if n_ones == 65535:
                    mask.append(n_ones)
                    mask.append(0)
                    n_ones = 0
                n_ones += 1
            else:
                mask.append(n_ones)
                c_0 = True
                n_ones = 0
                n_zeros += 1
    if c_0:
        mask.append(n_zeros)
    else:
        mask.append(n_ones)
    
    
    with open(f_name, 'wb') as file:
        file.write(struct.pack('4c',b'C',b'M',b'A',b'P'))
        file.write(struct.pack('I',len(ext_indexes)))
        file.write(struct.pack('4c',b'A',b'B',b'C',b'C'))
        file.write(struct.pack('d',7.5))
        file.write((ext_indexes-1).tobytes())
        file.write(struct.pack('I',len(mask)))
        file.write(np.array(mask, dtype=np.uint16).tobytes())
        file.write(struct.pack('I', n_frames))
        file.write(struct.pack('I', 2))
        
        for i in range(n_frames):
            file.write(struct.pack('d', colvar[i]))
            file.write(struct.pack('d', lambdas[i]))
            file.write((cmap_list[i,upper2lower] * 65535 + 0.5).astype(np.uint16).tobytes())

def loadCmapNew(f_name, lower_t = False):
    with open(f_name, 'rb') as file:
        file.read(struct.calcsize('4c'))
        len_ext_ind =  struct.unpack('I', file.read(struct.calcsize('I')))[0]
        c_type = struct.unpack('4c', file.read(struct.calcsize('4c')))
        c_type = ''.join([x.decode('UTF-8') for x in c_type])
        unit =  struct.unpack('d', file.read(struct.calcsize('d')))[0]
        ext_indexes = np.frombuffer(file.read(4 * len_ext_ind), dtype=np.uint32)
        
        indexes_all = np.array([(j,i) for i, j in itertools.combinations(ext_indexes,2)], dtype=np.uint32)
        indexes_all = indexes_all[np.lexsort((indexes_all[:,1], indexes_all[:,0]))]
        
        len_mask = struct.unpack('I', file.read(struct.calcsize('I')))[0]
        mask = np.frombuffer(file.read(2 * len_mask), dtype=np.uint16)
        
        bool_mask = np.full(indexes_all.shape[0], False)
        c_0 = True
        cur_pos = 0
        for i in mask:
            if c_0:
                cur_pos += i
            else:
                bool_mask[cur_pos:cur_pos+i] = True
                cur_pos += i
            c_0 = not c_0
        
        indexes = indexes_all[bool_mask]
        if not lower_t:
            lower2upper = np.lexsort((indexes[:,0], indexes[:,1]))
            indexes = indexes[lower2upper,::-1]
        
        n_frames = struct.unpack('I', file.read(struct.calcsize('I')))[0]
        precision = struct.unpack('I', file.read(struct.calcsize('I')))[0]
        
        colvar = np.empty(n_frames)
        lambdas = np.empty(n_frames)
        cmaps = np.empty((n_frames, indexes.shape[0]))
        
        for i in range(n_frames):
            colvar[i] = struct.unpack('d', file.read(struct.calcsize('d')))[0]
            lambdas[i] = struct.unpack('d', file.read(struct.calcsize('d')))[0]
            cmaps[i] = np.frombuffer(file.read(precision * indexes.shape[0]), dtype=np.uint16) / 65535
            if not lower_t:
                cmaps[i] = cmaps[i,lower2upper]
            
        return (c_type, unit, indexes, colvar, lambdas, cmaps)

def extractReactive(varible, start_tresh, end_tresh):
    reactive_portions = []
    for starting_point in range(len(varible)):
        if start_tresh <= end_tresh and varible[starting_point] < start_tresh:
            break
        if start_tresh > end_tresh and varible[starting_point] > start_tresh:
            break

    cursor = starting_point
    for i in range(starting_point, len(varible)):
        if start_tresh <= end_tresh and varible[i] < start_tresh:
            cursor = i
        if start_tresh > end_tresh and varible[i] > start_tresh:
            cursor = i

        if start_tresh <= end_tresh and varible[i] > end_tresh and (len(reactive_portions) == 0 or cursor > reactive_portions[-1][0]):
            reactive_portions.append((cursor,i))
        if start_tresh > end_tresh and varible[i] < end_tresh and (len(reactive_portions) == 0 or cursor > reactive_portions[-1][0]):
            reactive_portions.append((cursor,i))
    
    return reactive_portions

def downsample(values, criterion='qe', num=10, dist=0.1, margin_factor=1):
        if criterion == 'qe':
            frames = np.zeros(num, dtype=int)
            steps = np.linspace(np.min(values), np.max(values), num)
            for i in range(num):
                frames[i] = np.argmin(np.abs(values-steps[i]))
            
            return frames
            
        elif criterion == 'qp':
            frames = [0]
            warning = False
            i = 0
            for val in values:
                delta_norm = np.abs(val - values[frames[-1]])
                if delta_norm > dist:
                    warning = True
                    if delta_norm < dist + dist*margin_factor:
                        warning = False
                        frames.append(i)
                i+=1

            if warning:
                warnings.warn('Downsampled frames can be interrupted. Distance could be too high.')
            
            return np.array(frames, dtype=int)
            
        else:
            raise TypeError('Parameter "' + str(criterion) + '" does not match any valid criterion')

@njit(parallel=True)
def compute_sigma(cmaps_traj, cmaps_ref_trajs, lamdas, colvars):
    F = cmaps_traj.shape[0]
    R = cmaps_ref_trajs.shape[0]
    N = cmaps_traj.shape[1]
    sigma = np.zeros(cmaps_traj.shape[0], dtype=np.float32)
    for i in prange(F):
        w_tot = 0
        for j in prange(R):
            d = 0
            for k in prange(N):
                tmp = cmaps_traj[i,k] - cmaps_ref_trajs[j,k]
                d += tmp * tmp
            w = np.exp(-d * lamdas[j])
            w_tot += w
            sigma[i] += w * colvars[j]
        sigma[i] /= w_tot
    return sigma

@njit(parallel=True)
def compute_sigma_abcc(cmaps_traj, cmaps_ref_trajs, lamdas, colvars):
    F = cmaps_traj.shape[0]
    R = cmaps_ref_trajs.shape[0]
    N = cmaps_traj.shape[1]
    sigma = np.zeros(cmaps_traj.shape[0], dtype=np.float32)
    for i in prange(F):
        w_tot = 0
        for j in prange(R):
            d = 0
            for k in prange(N):
                if cmaps_traj[i,k] > 0:
                    tmp = cmaps_traj[i,k] - cmaps_ref_trajs[j,k]
                    d += tmp * tmp
            w = np.exp(-d * lamdas[j])
            w_tot += w
            sigma[i] += w * colvars[j]
        sigma[i] /= w_tot
    return sigma

@njit(parallel=True)
def _fast_contacts_dist_int(indexes, min_dist = 35):
    keep_index = []
    for i in range(indexes.shape[0]):
        for j in range(i+1,indexes.shape[0]):
            if np.abs(indexes[i] - indexes[j]) > min_dist:
                keep_index.append((i,j))
    return np.array(keep_index, dtype=np.uint32)

@njit(parallel=True)
def _fast_contacts_dist_ext(indexes, min_dist = 35):
    keep_index = []
    for i in range(indexes.shape[0]):
        for j in range(i+1,indexes.shape[0]):
            if np.abs(indexes[i] - indexes[j]) > min_dist:
                keep_index.append((indexes[i], indexes[j]))
    return np.array(keep_index, dtype=np.uint32)


def internal_signature_dist(atom_group, min_dist = 35):
    return _fast_contacts_dist_int(atom_group.ids, min_dist=min_dist)

def external_signature_dist(atom_group, min_dist = 35):
    return _fast_contacts_dist_ext(atom_group.ids, min_dist=min_dist)

@njit
def identity(x):
    return x

#@njit(parallel=True)
def bare_cmap(X, signature, map_fun, res):
    F = X.shape[0] # numero di frame
    M = signature.shape[0] # numero di contatti trovati dalla funzione signature
    N = X.shape[2] # cordinate di ciascun atomo
    for f in range(F):
        for i in range(M):
            d = 0.0
            for k in range(N):
                # tmp = X(position)[frame, atom involved in contact, coordinate(x,y,z)] -  X(position)[frame, other atom involved in contact, coordinate(x,y,z)]
                tmp = X[f, signature[i,0], k] - X[f, signature[i,1], k]
                # not performed the sqrt because in the next passage we arise at power 3 not 6
                d += tmp * tmp
            res[f, i] += map_fun(d)

@njit
def sigmoid_squared(x):
    if x > 151.29:
        # cut off  rc = 12.3 A so rij > rc was set to zero
        return 0
    if np.abs(x-56.25) < 10e-5:
        return 0.6
    return (1-(x/56.25)**3)/(1-(x/56.25)**5)

def compute_cmap(trajectory, start=None, end=None, use_ref_as_last = False, ref = None, min_dist = 35, 
                 verbose = False, map_function = sigmoid_squared, selection='all', signature_function = internal_signature_dist,
                 res = None, stride=1):

    atoms = trajectory.getUniverse().select_atoms(selection)
    start_f = 0
    end_f = trajectory.getUniverse().trajectory.n_frames
    
    if start is not None:
        start_f = start
    if end is not None:
        end_f = end

    frames = np.arange(start=start_f, stop=end_f, step=stride, dtype=np.int32)

    n_frames = frames.shape[0]
    
    
    signature = signature_function(atoms, min_dist=min_dist)
    # hopfield network understanding
    selected_rows = signature[signature[:, 0] == 0]
    selected_rows_2 = signature[signature[:, 0] == 545]
    print(f"numero di contatti da valutare per l'atomo 0 : {selected_rows.shape}")
    print(f"numero di contatti da valutare per l'atomo 545 : {selected_rows_2.shape}")

    print(f"totale numero di contatti: {signature.shape}")
    
    positions = np.empty((n_frames, atoms.n_atoms, 3), dtype=np.float32)
    print(f"numero di atomi: {atoms.n_atoms}")
    print(f"numero di frame: {n_frames}")

    for i in range(frames.shape[0]):
        trajectory.getUniverse().trajectory[frames[i]]
        positions[i] = atoms.positions

 
    print(f"position shape: {positions.shape}")
    
    if res is not None:
        cmap_traj = res
    else:
        cmap_traj = np.zeros((n_frames, signature.shape[0]), dtype=np.float32)
        
    bare_cmap(positions, signature, map_function, cmap_traj)
    

    if use_ref_as_last and ref is not None:
        ref_pos = ref.getUniverse().select_atoms(selection).positions
        bare_cmap(ref_pos[None,:], signature, map_function, cmap_traj[-1][None,:])
    
    if res is None:
        return cmap_traj

@njit(parallel=True)
def cmap2hard_q(cmap_traj, reference_cmap, cutoff=0.6):
    mask = reference_cmap > cutoff
    n_native = np.sum(mask)
    F = cmap_traj.shape[0]
    q = np.zeros(F)
    for i in prange(F):
        q[i] = np.sum(cmap_traj[i][mask] > cutoff) / n_native
    
    return q
