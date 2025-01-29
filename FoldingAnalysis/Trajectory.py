# -*- coding: utf-8 -*- 
from FoldingAnalysis.utilities import *
from FoldingAnalysis.utilities import _Struct
from FoldingAnalysis.exceptions import *
from FoldingAnalysis.BiasProperties import BiasProperties

import MDAnalysis as md
from MDAnalysis.analysis import distances
from MDAnalysis.analysis import contacts
from MDAnalysis.analysis import align
import numpy as np
import itertools 
import warnings
import os
import matplotlib.pyplot as plt

class Trajectory:
    def __init__(self, topology_file, trajectory_file=None, bias_properties=None, 
                 name=None, reference=None, in_memory=False, dt=None):
        
        self.name = name if name is not None else \
        (os.path.splitext(trajectory_file)[0].replace('/','_').replace('\\','_') if isinstance(trajectory_file, str) else \
        os.path.splitext(topology_file)[0]).replace('/','_').replace('\\','_')
        
        if trajectory_file is not None:
            self.universe = md.Universe(topology_file, trajectory_file, in_memory=in_memory)
        else:
            self.universe = md.Universe(topology_file, in_memory=in_memory)
        
        if reference is not None and not isinstance(reference, Trajectory):
            raise TypeError('reference must be an instance of Trajectory')
            
        if bias_properties is not None and not isinstance(bias_properties, BiasProperties):
            raise TypeError('bias_properties must be an instance of BiasProperties')
        
        self.reference = reference
        self.Q = None
        self.Q_settings = _Struct()
        self.Q_settings.use_ref_as_last = False
        self.Q_settings.min_dist = 3
        self.Q_settings.cutoff = 4.5
        self.Q_settings.beta_c = 5
        self.Q_settings.lambda_c = 1.8
        
        self.Cmap = None
        self.Cmap_settings = _Struct()
        self.Cmap_settings.use_ref_as_last = False
        self.Cmap_settings.min_dist = 35
        self.Cmap_settings.map_function = sigmoid_squared
        self.Cmap_settings.verbose = False
        self.Cmap_settings.start = None
        self.Cmap_settings.end = None
        self.Cmap_settings.selection = 'all'
        
        self.RMSD = None
        self.RMSD_settings = _Struct()
        self.RMSD_settings.selection = 'all'
        
        self.DT = _Struct()
        self.DT.dt = 0 if trajectory_file is None else self.universe.trajectory.dt
        self.DT.guessed = True
        if dt is not None:
            self.DT.dt = float(dt)
            self.DT.guessed = False
           
        
        self.FoldingTime = None
        self.FoldingFrame = None
        self.Folding_settings = _Struct()
        self.Folding_settings.method = 'simple'
        self.Folding_settings.threshold = 2
        self.Folding_settings.tolerance = 0.5
        self.Folding_settings.ignore_last_frames = 0
        self.Folding_settings.ignore_last_time = 0
        
        self.bias_properties = bias_properties
        
    def getUniverse(self):
        return self.universe
    
    def getName(self):
        return self.name
    
    def setName(self, name):
        self.name = name
    
    def getFrameCount(self):
        return self.universe.trajectory.n_frames
    
    def getAtomCount(self):
        return self.universe.trajectory.n_atoms
    
    def getDt(self):
        if self.DT.guessed:
            warnings.warn('Timestep has been guessed from trajectory and might be wrong.\
                           To be sure set it explicitely with setDt(dt)', GuessedTimestep)
        return self.DT.dt
    
    def setDt(self, dt):
        if dt is None:
            raise TypeError('dt cannot be None')
        
        self.DT.dt = float(dt)
        self.DT.guessed = False
        
        self.FoldingTime = None
        self.FoldingFrame = None
        
    def frameToTime(self, frame):
        if self.DT.guessed:
            warnings.warn('Timestep has been guessed from trajectory and might be wrong.\
                           To be sure set it explicitely with setDt(dt)', GuessedTimestep)
            
        return frame*self.DT.dt
            
    def timeToFrame(self, time):
        if self.DT.guessed:
            warnings.warn('Timestep has been guessed from trajectory and might be wrong.\
                           To be sure set it explicitely with setDt(dt)', GuessedTimestep)
        
        return int(round(time/self.DT.dt))
    
    def setReference(self, reference):
        if reference is not None and not isinstance(reference, Trajectory):
            raise TypeError('reference must be an instance of Trajectory')
        
        self.reference = reference
        self.Q = None
        self.Cmap = None
        self.RMSD = None
        self.FoldingTime = None
        self.FoldingFrame = None
        
    def getReference(self):
        return self.reference
        
    def getQ(self):
        if self.reference is None:
            raise MissingReference('Missing reference structure! Add one with setReference(reference)')
            
        if self.Q is None:
            self.Q = bestHummerQ(self, self.reference, use_ref_as_last=self.Q_settings.use_ref_as_last, 
                                min_dist=self.Q_settings.min_dist, cutoff=self.Q_settings.cutoff, 
                                beta_c=self.Q_settings.beta_c, lambda_c=self.Q_settings.lambda_c)
            return self.Q
        else:
            return self.Q
    
    def configureQ(self, use_ref_as_last = False, min_dist = 3, cutoff = 4.5, beta_c = 5, lambda_c = 1.8):
        self.Q_settings.use_ref_as_last = use_ref_as_last
        self.Q_settings.min_dist = min_dist
        self.Q_settings.cutoff = cutoff
        self.Q_settings.beta_c = beta_c
        self.Q_settings.lambda_c = lambda_c
        self.Q = None
        
    def getCmap(self, res = None, cache = False, stride=1):
        if self.Cmap_settings.use_ref_as_last and self.reference is None:
            raise MissingReference('Missing reference structure! Add one with setReference(reference) or \
                                    use configureCmap(use_ref_as_last=False)')
        
            print("##############################")


        if self.Cmap is None:
            if cache:
                self.Cmap = compute_cmap(self, use_ref_as_last=self.Cmap_settings.use_ref_as_last, 
                                        ref=self.reference, min_dist=self.Cmap_settings.min_dist,
                                        start=self.Cmap_settings.start, end=self.Cmap_settings.end,
                                        verbose=self.Cmap_settings.verbose, 
                                        map_function=self.Cmap_settings.map_function, 
                                        selection=self.Cmap_settings.selection,
                                        res=res, stride=stride)
                return self.Cmap
            else:
                return compute_cmap(self, use_ref_as_last=self.Cmap_settings.use_ref_as_last, 
                                   ref=self.reference, min_dist=self.Cmap_settings.min_dist,
                                   start=self.Cmap_settings.start, end=self.Cmap_settings.end,
                                   verbose=self.Cmap_settings.verbose, 
                                   map_function=self.Cmap_settings.map_function,
                                   selection=self.Cmap_settings.selection,
                                   res=res, stride=stride)
        else:
            return self.Cmap
    
    def getCmapAtFrame(self, frame, res = None):
        if frame < 0 or frame >= self.universe.trajectory.n_frames:
            raise FrameOutOfBounds('Selected frame does not exist')
        
        return compute_cmap(self, start=frame, end=frame+1, use_ref_as_last=False, 
                           min_dist=self.Cmap_settings.min_dist, map_function=self.Cmap_settings.map_function,
                           selection=self.Cmap_settings.selection,
                           res=res)
    
    def getCmapDim(self):
        return self.getCmapSignature().shape[0]
    
    def getCmapSignature(self):
        atoms = self.getUniverse().select_atoms(self.Cmap_settings.selection)
        return external_signature_dist(atoms, self.Cmap_settings.min_dist)
    
    def configureCmap(self, start=None, end=None, use_ref_as_last = False, min_dist = 35, verbose=False, 
                      map_function=sigmoid_squared, selection='all'):
        self.Cmap_settings.use_ref_as_last = use_ref_as_last
        self.Cmap_settings.min_dist = min_dist
        self.Cmap_settings.map_function = map_function
        self.Cmap_settings.verbose = verbose
        self.Cmap_settings.start = start
        self.Cmap_settings.end = end
        self.Cmap_settings.selection = selection
        self.Cmap = None
        
        
    def getRMSD(self):
        if self.reference is None:
            raise MissingReference('Missing reference structure! Add one with setReference(reference)')
          
        if self.RMSD is None:
            
            self.RMSD = computeRMSD(self, self.reference, selection=self.RMSD_settings.selection)
            return self.RMSD
        else:
            return self.RMSD
        
    def configureRMSD(self, selection='all'):
        self.RMSD_settings.selection = selection
        self.RMSD = None
        self.FoldingTime = None
        self.FoldingFrame = None
        
        
    def configureFolding(self, method='simple', threshold = 2, tolerance = 0.5, ignore_last_frames=0,
                         ignore_last_time=0):
        self.Folding_settings.method = method
        self.Folding_settings.threshold = threshold
        self.Folding_settings.tolerance = tolerance
        self.Folding_settings.ignore_last_frames = ignore_last_frames
        self.Folding_settings.ignore_last_time = ignore_last_time
        self.FoldingTime = None
        self.FoldingFrame = None
        
        
    
    def getFoldingFrame(self):
        if self.FoldingFrame is None:
            self.FoldingFrame = computeFoldingFrame(self.getRMSD(), method=self.Folding_settings.method, 
                                       threshold=self.Folding_settings.threshold, 
                                       tolerance=self.Folding_settings.tolerance,
                                       ignore_last_frames=self.Folding_settings.ignore_last_frames)
            return self.FoldingFrame
        else:
            return self.FoldingFrame
        
    def getFoldingTime(self):
        if self.FoldingTime is None:
            self.FoldingTime = computeFoldingTime(self.getRMSD(), self.getDt(), method=self.Folding_settings.method, 
                                       threshold=self.Folding_settings.threshold, 
                                       tolerance=self.Folding_settings.tolerance,
                                       ignore_last_time=self.Folding_settings.ignore_last_time)
            return self.FoldingTime
        else:
            return self.FoldingTime
    
    def hasFolded(self):
        return self.getFoldingTime() > 0
    
    def plotRMSD(self, directory, filename=None):
        directory = directory.rstrip('/')
        fig, ax = plt.subplots()
        ax.plot(np.array(range(len(self.getRMSD())))*self.getDt(),self.getRMSD(), linewidth=0.5)
        ax.axhline(y=self.Folding_settings.threshold,linewidth=0.75, linestyle='--', color='gray')
        if self.hasFolded():
            ax.plot(self.getFoldingTime(),self.getRMSD()[self.timeToFrame(self.getFoldingTime())],'o',markersize=5)
            ax.axvline(x=self.getFoldingTime(),linewidth=1,linestyle=':',color='C1')
            
        ax.set_xlabel('Time (ps)')
        ax.set_ylabel('RMSD (Å)')
        ax.minorticks_on()
        ax.grid(color='lightgray', linestyle='-.', linewidth=0.25, which='both')
        fname = filename if filename is not None else (self.name + '_rmsd.pdf')
        fig.savefig(directory + '/' + fname)
        plt.close(fig)
        
    def getBiasProperties(self):
        return self.bias_properties
    
    def setBiasProperites(self, bias_properties):
        if bias_properties is not None and not isinstance(bias_properties, BiasProperties):
            raise TypeError('bias_properties must be an instance of BiasProperties')
        
        self.bias_properties = bias_properties
    
    def getPenalityAtFolding(self):
        try:
            pen = self.bias_properties.bias_penality
            timesteps = self.bias_properties.time
        except:
            raise BiasPropertiesError('Cannot access the needed properties (time and bias_penality) \
                                       of bias_properties. Set it with a correct BiasProperty object')
        
        return pen[np.abs(timesteps-self.getFoldingTime()).argmin()]
    
    def getPenalityNormalAtFolding(self):
        try:
            pen = self.bias_properties.bias_penality
            timesteps = self.bias_properties.time
        except:
            raise BiasPropertiesError('Cannot access the needed properties (time and bias_penality) \
                                       of bias_properties. Set it with a correct BiasProperty object')
        try:
            tot_forces = self.bias_properties.cum_tot_force
            return pen[np.abs(timesteps-self.getFoldingTime()).argmin()] / tot_forces[np.abs(timesteps-self.getFoldingTime()).argmin()]
        except:
            pass
        
        return pen[np.abs(timesteps-self.getFoldingTime()).argmin()] / timesteps[np.abs(timesteps-self.getFoldingTime()).argmin()]
    
    def getPenalityAtTime(self, time):
        try:
            pen = self.bias_properties.bias_penality
            timesteps = self.bias_properties.time
        except:
            raise BiasPropertiesError('Cannot access the needed properties (time and bias_penality) \
                                       of bias_properties. Set it with a correct BiasProperty object')
        
        return pen[np.abs(timesteps-time).argmin()]
    
    def downsample(self, criterion='qe', num=10, dist=0.1, margin_factor=1):
        if criterion == 'qe':
            frames = np.zeros(num, dtype=int)
            steps = np.linspace(np.min(self.getQ()), np.max(self.getQ()), num)
            for i in range(num):
                frames[i] = np.argmin(np.abs(self.getQ()-steps[i]))
            
            return frames
            
        elif criterion == 'qp':
            frames = [0]
            warning = False
            i = 0
            for val in self.getQ():
                delta_norm = np.abs(val - self.getQ()[frames[-1]])
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
