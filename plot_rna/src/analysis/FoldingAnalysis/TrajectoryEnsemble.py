# -*- coding: utf-8 -*- 
from FoldingAnalysis.utilities import *
from FoldingAnalysis.utilities import _Struct
from FoldingAnalysis.exceptions import *
from FoldingAnalysis.Trajectory import Trajectory
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

class TrajectoryEnsemble:
    def __init__(self, reference = None, trajectories = None, dt = None):
        if reference is not None and not isinstance(reference, Trajectory):
            raise TypeError('reference must be an instance of Trajectory')
        
        self.trajectories = None
        self.dt = None
        
        if trajectories is not None:
            if isinstance(trajectories, Trajectory):
                self.trajectories = [trajectories]
            elif isinstance(trajectories, list) and \
            all([isinstance(trajectories[i], Trajectory) for i in range(len(trajectories))]):
                self.trajectories = trajectories
            else:
                raise TypeError('trajectories must be a Trajectory or a list of Trajectory')
        
        if dt is not None:
            self.setDt(dt)
        
        #self.trajectories = trajectories
        self.reference = None
        self.setReference(reference)
        
        self.Q_settings = _Struct()
        self.configureQ()
        
        self.Cmap_settings = _Struct()
        self.configureCmap()
        
        self.Folding_settings = _Struct()
        self.configureFolding()
        
        self.RMSD_settings = _Struct()
        self.configureRMSD()
        
        self.foldingTrajectories = None
    
    def addTrajectory(self, trajectory):
        if isinstance(trajectory, Trajectory):
            if self.trajectories is None:
                self.trajectories = [trajectory]
            else:
                self.trajectories.append(trajectory)
        elif isinstance(trajectory, list) and \
        all([isinstance(trajectory[i], Trajectory) for i in range(len(trajectory))]):
            if self.trajectories is None:
                self.trajectories = trajectory
            else:
                self.trajectories.extend(trajectory)
        else:
            raise TypeError('trajectory must be a Trajectory or a list of Trajectory')
            
        self.setReference(self.reference)
        
        self.configureCmap(use_ref_as_last=self.Cmap_settings.use_ref_as_last, min_dist=self.Cmap_settings.min_dist,
                           verbose=self.Cmap_settings.verbose, start=self.Cmap_settings.start, 
                           end=self.Cmap_settings.end, map_function=self.Cmap_settings.map_function,
                           selection=self.Cmap_settings.selection)
        
        self.configureQ(use_ref_as_last = self.Q_settings.use_ref_as_last, min_dist = self.Q_settings.min_dist, 
                        cutoff = self.Q_settings.cutoff, beta_c = self.Q_settings.beta_c, 
                        lambda_c = self.Q_settings.lambda_c)
        
        self.configureFolding(method=self.Folding_settings.method, threshold = self.Folding_settings.threshold, 
                              tolerance = self.Folding_settings.tolerance, 
                              ignore_last_frames=self.Folding_settings.ignore_last_frames,
                              ignore_last_time=self.Folding_settings.ignore_last_time)
        
        self.configureRMSD(selection=self.RMSD_settings.selection)
        
        self.foldingTrajectories = None
        
        if self.dt is not None:
            self.setDt(self.dt)
                    
    def getTrajectories(self):
        return self.trajectories
    
    def getTrajectoriesCount(self):
        return 0 if self.trajectories is None else len(self.trajectories)
    
    def setReference(self, reference):
        if reference is not None and not isinstance(reference, Trajectory):
            raise TypeError('reference must be an instance of Trajectory')
        
        self.reference = reference
        if self.trajectories is not None:
            for traj in self.trajectories:
                traj.setReference(reference)
                
        self.foldingTrajectories = None
                
    def getReference(self):
        return self.reference
    
    def setDt(self, dt):
        if dt is None:
            raise TypeError('dt cannot be None')
        
        self.dt = dt
        if self.trajectories is not None:
            for traj in self.trajectories:
                traj.setDt(dt)
                
        self.foldingTrajectories = None
        
    def configureRMSD(self, selection='all'):
        self.RMSD_settings.selection = selection
        if self.trajectories is not None:
            for traj in self.trajectories:
                traj.configureRMSD(selection=selection)
                
        self.foldingTrajectories = None
    
    def configureCmap(self, start=None, end=None, use_ref_as_last = True, min_dist = 35, verbose=False, 
                      map_function=sigmoid_squared, selection='all'):
        self.Cmap_settings.use_ref_as_last = use_ref_as_last
        self.Cmap_settings.min_dist = min_dist
        self.Cmap_settings.map_function = map_function
        self.Cmap_settings.verbose = verbose
        self.Cmap_settings.start = start
        self.Cmap_settings.end = end
        self.Cmap_settings.selection = selection
        if self.trajectories is not None:
            for traj in self.trajectories:
                traj.configureCmap(use_ref_as_last = use_ref_as_last, min_dist = min_dist, start=start, end=end,
                                   verbose=verbose, map_function=map_function, selection=selection)
    
    def configureQ(self, use_ref_as_last = True, min_dist = 3, cutoff = 4.5, beta_c = 5, lambda_c = 1.8):
        self.Q_settings.use_ref_as_last = use_ref_as_last
        self.Q_settings.min_dist = min_dist
        self.Q_settings.cutoff = cutoff
        self.Q_settings.beta_c = beta_c
        self.Q_settings.lambda_c = lambda_c
        if self.trajectories is not None:
            for traj in self.trajectories:
                traj.configureQ(use_ref_as_last = use_ref_as_last, min_dist = min_dist, 
                                cutoff = cutoff, beta_c = beta_c, lambda_c = lambda_c)
    
    def configureFolding(self, method='simple', threshold = 2, tolerance = 0.5, ignore_last_frames=0,
                         ignore_last_time=0):
        self.Folding_settings.method = method
        self.Folding_settings.threshold = threshold
        self.Folding_settings.tolerance = tolerance
        self.Folding_settings.ignore_last_frames = ignore_last_frames
        self.Folding_settings.ignore_last_time = ignore_last_time
        if self.trajectories is not None:
            for traj in self.trajectories:
                traj.configureFolding(method=method, threshold = threshold, tolerance = tolerance, 
                                      ignore_last_frames=ignore_last_frames, ignore_last_time=ignore_last_time)
                
        self.foldingTrajectories = None
    
    def getAvgCmapQ(self, n = 40, folded_only = True):
        if self.trajectories is None or len(self.trajectories) == 0:
            raise EmptyEnsemble('No trajectories in the ensemble. Add at least one with addTrajectory(trajectory)')
    
        if self.reference is None:
            raise MissingReference('Missing reference structure! Add one with setReference(reference)')
        
        
        
        if folded_only:
            working_trajectories = self.getFoldedTrajectories()
        else:
            working_trajectories = self.getTrajectories()
        
        if folded_only and len(working_trajectories) == 0:
            raise EmptyEnsemble('No trajectory has folded. Set folded_only to False to compute the average anyway')
        
        start = 0
        end = 1
        for traj in working_trajectories:
            if np.min(traj.getQ()) > start:
                start = np.min(traj.getQ())
            if np.max(traj.getQ()) < end:
                end = np.max(traj.getQ())

        steps = np.linspace(start, end, n)
        average_cmaps = np.zeros((n, working_trajectories[0].getCmapDim()), dtype=np.float32)
        
        for traj in working_trajectories:
            cmaps = traj.getCmap()
            j = 0
            for ts in steps:
                index = np.argmin(np.abs(traj.getQ() - ts))
                average_cmaps[j] += cmaps[index]
                j+=1
            print('> Done with '+traj.getName())
        average_cmaps /= len(working_trajectories)

        return average_cmaps
    
    def getAvgCmapQBins(self, n = 40, folded_only = True):
        if self.trajectories is None or len(self.trajectories) == 0:
            raise EmptyEnsemble('No trajectories in the ensemble. Add at least one with addTrajectory(trajectory)')
    
        if self.reference is None:
            raise MissingReference('Missing reference structure! Add one with setReference(reference)')
            
        if folded_only:
            working_trajectories = self.getFoldedTrajectories()
        else:
            working_trajectories = self.getTrajectories()
        
        if folded_only and len(working_trajectories) == 0:
            raise EmptyEnsemble('No trajectory has folded. Set folded_only to False to compute the average anyway')
        
        start = 1
        end = 0
        for traj in working_trajectories:
            if np.min(traj.getQ()) < start:
                start = np.min(traj.getQ())
            if np.max(traj.getQ()) > end:
                end = np.max(traj.getQ())

        bins = np.linspace(start, end+0.001, n+1)
        bins_count = np.zeros(n)

        average_cmaps = np.zeros((n, working_trajectories[0].getCmapDim()), dtype=np.float32)
        
        for traj in working_trajectories:
            cmaps = traj.getCmap()
            for j in range(n):
                indexs = np.logical_and(traj.getQ() >= bins[j], traj.getQ() < bins[j+1])
                bins_count[j] += np.sum(indexs)
                average_cmaps[j] += np.sum(cmaps[indexs], axis=0)
            print('> Done with '+traj.getName())

        i = 0
        for div in bins_count:
            average_cmaps[i] /= div
            i+=1

        return average_cmaps
    
    def getAverageCmapTime(self, folded_only = True, res=None):
        if self.trajectories is None or len(self.trajectories) == 0:
            raise EmptyEnsemble('No trajectories in the ensemble. Add at least one with addTrajectory(trajectory)')
    
        if self.reference is None:
            raise MissingReference('Missing reference structure! Add one with setReference(reference)')
            
        if folded_only:
            working_trajectories = self.getFoldedTrajectories()
        else:
            working_trajectories = self.getTrajectories()
        
        if folded_only and len(working_trajectories) == 0:
            raise EmptyEnsemble('No trajectory has folded. Set folded_only to False to compute the average anyway')
        
        #n_frames = working_trajectories[0].getFrameCount()

        start_f = 0
        end_f = working_trajectories[0].getFrameCount()
        print('> number of frames in one trajectories: '+str(working_trajectories[0].getCmapDim()))
        
        if self.Cmap_settings.start is not None:
            start_f = self.Cmap_settings.start
        if self.Cmap_settings.end is not None:
            end_f = self.Cmap_settings.end

        n_frames = end_f - start_f

        if self.Cmap_settings.end is not None and not all([traj.getFrameCount() >= end_f for traj in working_trajectories]):
            raise FrameOutOfBounds('Folding trajetotires have a different number of frames')
        
        if not all([traj.getFrameCount() == end_f for traj in working_trajectories]) and self.Cmap_settings.end is None:
            raise FrameOutOfBounds('Folding trajetotires have a different number of frames')
        
        if self.Cmap_settings.use_ref_as_last and self.reference is not None:
            n_frames += 1
        
        if res is not None:
            average_cmaps = res
        else:    
            average_cmaps = np.zeros((n_frames, working_trajectories[0].getCmapDim()), dtype=np.float32)
        
        for traj in working_trajectories:
            traj.getCmap(res=average_cmaps)
            print('> Done with '+traj.getName())
            
        average_cmaps /= len(working_trajectories)
        
        if res is None:
            return average_cmaps
    
    def getAverageRMSDTime(self, folded_only=True):
        if self.trajectories is None or len(self.trajectories) == 0:
            raise EmptyEnsemble('No trajectories in the ensemble. Add at least one with addTrajectory(trajectory)')
    
        if self.reference is None:
            raise MissingReference('Missing reference structure! Add one with setReference(reference)')
        
        if folded_only:
            working_trajectories = self.getFoldedTrajectories()
        else:
            working_trajectories = self.getTrajectories()
        
        if folded_only and len(working_trajectories) == 0:
            raise EmptyEnsemble('No trajectory has folded. Set folded_only to False to compute the average anyway')
            
        n_frames = working_trajectories[0].getFrameCount()
        
        if not all([traj.getFrameCount() == n_frames for traj in working_trajectories]):
            raise FrameOutOfBounds('Folding trajetotires have a different number of frames')
        
        average_RMSD = np.zeros(n_frames)
        for traj in working_trajectories:
            average_RMSD += traj.getRMSD()
            
        average_RMSD /= len(working_trajectories)
        
        return average_RMSD
        
    def getAverageRMSDQ(self, n=40, folded_only=True):
        if self.trajectories is None or len(self.trajectories) == 0:
            raise EmptyEnsemble('No trajectories in the ensemble. Add at least one with addTrajectory(trajectory)')
    
        if self.reference is None:
            raise MissingReference('Missing reference structure! Add one with setReference(reference)')
        
        if folded_only:
            working_trajectories = self.getFoldedTrajectories()
        else:
            working_trajectories = self.getTrajectories()
        
        if folded_only and len(working_trajectories) == 0:
            raise EmptyEnsemble('No trajectory has folded. Set folded_only to False to compute the average anyway')
        
        old_ural = self.Q_settings.use_ref_as_last
        self.configureQ(use_ref_as_last = False, min_dist = self.Q_settings.min_dist, 
                        cutoff = self.Q_settings.cutoff, beta_c = self.Q_settings.beta_c, 
                        lambda_c = self.Q_settings.lambda_c)
        
        start = 0
        end = 1
        for traj in working_trajectories:
            if np.min(traj.getQ()) > start:
                start = np.min(traj.getQ())
            if np.max(traj.getQ()) < end:
                end = np.max(traj.getQ())

        steps = np.linspace(start, end, n)
        average_RMSD = np.zeros(n)
        
        for traj in working_trajectories:
            RMSD = traj.getRMSD()
            j = 0
            for ts in steps:
                index = np.argmin(np.abs(traj.getQ() - ts))
                average_RMSD[j] += RMSD[index]
                j+=1
        
        average_RMSD /= len(working_trajectories)
        
        self.configureQ(use_ref_as_last = old_ural, min_dist = self.Q_settings.min_dist, 
                        cutoff = self.Q_settings.cutoff, beta_c = self.Q_settings.beta_c, 
                        lambda_c = self.Q_settings.lambda_c)
        
        return average_RMSD
    
    def getAveragePenality(self, folded_only=True):
        if self.trajectories is None or len(self.trajectories) == 0:
            raise EmptyEnsemble('No trajectories in the ensemble. Add at least one with addTrajectory(trajectory)')
        
        if folded_only:
            working_trajectories = self.getFoldedTrajectories()
        else:
            working_trajectories = self.getTrajectories()
        
        if folded_only and len(working_trajectories) == 0:
            raise EmptyEnsemble('No trajectory has folded. Set folded_only to False to compute the average anyway')
        
        try:
            pen_dim = len(working_trajectories[0].getBiasProperties().bias_penality)
            
            if not all([len(traj.getBiasProperties().bias_penality) == pen_dim for traj in working_trajectories]):
                raise FrameOutOfBounds('Folding trajetotires have a different length of bias properties')
        except:
            raise BiasPropertiesError('Cannot access the needed properties (time and bias_penality) \
                                       of bias_properties. Set it with a correct BiasProperty object for \
                                       every Trajectory')
        
        average_penality = np.zeros((len(working_trajectories),pen_dim))
        i=0
        for traj in working_trajectories:
            average_penality[i] = traj.getBiasProperties().bias_penality
            i+=1
            
        
        return (np.mean(average_penality,axis=0),np.std(average_penality,axis=0))
        
    
    def getFoldedTrajectories(self):
        if self.trajectories is None or len(self.trajectories) == 0:
            raise EmptyEnsemble('No trajectories in the ensemble. Add at least one with addTrajectory(trajectory)')
        
        if self.foldingTrajectories is None:
            self.foldingTrajectories = [traj for traj in self.trajectories if traj.hasFolded()]
            return self.foldingTrajectories
        else:
            return self.foldingTrajectories
        
    def getFoldedTrajectoriesCount(self):
        if self.trajectories is None or len(self.trajectories) == 0:
            raise EmptyEnsemble('No trajectories in the ensemble. Add at least one with addTrajectory(trajectory)')
        
        if self.foldingTrajectories is None:
            self.foldingTrajectories = [traj for traj in self.trajectories if traj.hasFolded()]
            return len(self.foldingTrajectories)
        else:
            return len(self.foldingTrajectories)
        
    def plotRMSDs(self, directory):
        if self.trajectories is not None:
            for traj in self.trajectories:
                traj.plotRMSD(directory)
                
    def plotMeanPenality(self, directory, filename=None, folded_only=True):
        mean_pen = self.getAveragePenality(folded_only=folded_only)[0]
        std_pen = self.getAveragePenality(folded_only=folded_only)[1]
        directory = directory.rstrip('/')
        plt.plot(mean_pen)
        plt.fill_between(np.arange(len(mean_pen)),mean_pen - std_pen, mean_pen + std_pen,alpha=0.3)
        plt.xlabel('Step')
        plt.ylabel('Penalty')
        plt.minorticks_on()
        plt.grid(color='lightgray', linestyle='-.', linewidth=0.25, which='both')
        fname = filename if filename is not None else ('mean_penalty.pdf')
        plt.savefig(directory + '/' + fname)
        plt.close()
        
    def plotMeanRMSD(self, directory, filename=None, folded_only=True):
        mean_pen = self.getAveragePenality(folded_only=folded_only)[0]
        std_pen = self.getAveragePenality(folded_only=folded_only)[1]
        directory = directory.rstrip('/')
        plt.plot(mean_pen)
        plt.fill_between(np.arange(len(mean_pen)),mean_pen - std_pen, mean_pen + std_pen,alpha=0.3)
        plt.xlabel('Step')
        plt.ylabel('Penalty')
        plt.minorticks_on()
        plt.grid(color='lightgray', linestyle='-.', linewidth=0.25, which='both')
        fname = filename if filename is not None else ('mean_penalty.pdf')
        plt.savefig(directory + '/' + fname)
        plt.close()
        
    
    def getDRP(self):
        if not self.getFoldedTrajectoriesCount() > 0:
            return None
        
        return self.getFoldedTrajectories()[np.argmin([traj.getPenalityAtFolding() for traj in \
                                                       self.getFoldedTrajectories()])]
    def getDRPNormal(self):
        if not self.getFoldedTrajectoriesCount() > 0:
            return None
        
        return self.getFoldedTrajectories()[np.argmin([traj.getPenalityNormalAtFolding() for traj in \
                                                       self.getFoldedTrajectories()])]
    def getDRPComplete(self):
        if not self.getFoldedTrajectoriesCount() > 0:
            return None
        max_ftime = np.max([traj.getFoldingTime() for traj in self.getFoldedTrajectories()])
        return self.getFoldedTrajectories()[np.argmin([traj.getPenalityAtTime(max_ftime) for traj in \
                                                       self.getFoldedTrajectories()])]
    def getMaxFoldingTime(self):
        if not self.getFoldedTrajectoriesCount() > 0:
            return None
        return np.max([traj.getFoldingTime() for traj in self.getFoldedTrajectories()])
    
    def getMinFoldingTime(self):
        if not self.getFoldedTrajectoriesCount() > 0:
            return None
        return np.min([traj.getFoldingTime() for traj in self.getFoldedTrajectories()])
    
    def getMeanFoldingTime(self):
        if not self.getFoldedTrajectoriesCount() > 0:
            return None
        return np.mean([traj.getFoldingTime() for traj in self.getFoldedTrajectories()])
    
    def getMedianFoldingTime(self):
        if not self.getFoldedTrajectoriesCount() > 0:
            return None
        return np.median([traj.getFoldingTime() for traj in self.getFoldedTrajectories()])
    
    def getStdFoldingTime(self):
        if not self.getFoldedTrajectoriesCount() > 0:
            return None
        return np.std([traj.getFoldingTime() for traj in self.getFoldedTrajectories()])
    
    def getSummary(self):
        summary = ''
        summary += '\n--------------------------------------------------------------------------------\n'
        summary += '                         ~ TRAJECTORY ENSEMBLE SUMMARY ~                           '
        summary += '\n--------------------------------------------------------------------------------\n'
        summary += '> Number of trajectories in the ensemble: ' + str(self.getTrajectoriesCount()) + '\n'
        summary += '> Timestep: ' + ('NOT SET' if self.dt is None else str(self.dt) + 'ps\n')
        summary += '> Reference Trajectory: ' + ('NOT SET' if self.reference is None else \
                                                 self.reference.getName()+'\n')
        summary += '\n--------------------------------------------------------------------------------\n'
        if not self.getFoldedTrajectoriesCount() > 0:
            summary += '> No folding trajectories in the ensemble!'
            return summary
        summary += '> A total of ' + str(self.getFoldedTrajectoriesCount()) + ' trajectories folded\n'
        for traj in sorted(self.getFoldedTrajectories(),key=lambda x: x.getFoldingTime()):
            summary += '> ' + traj.getName() + ' folded at ' + str(traj.getFoldingTime())
            
            try:
                summary += ' | Bias: ' + str(traj.getPenalityAtFolding()) + ' | Complete Bias: '+ \
                       str(traj.getPenalityAtTime(self.getMaxFoldingTime())) + ' | Normalised Bias: '+ \
                       str(traj.getPenalityNormalAtFolding()) + '\n'
            except:
                summary += '\n'
        
        summary += '> Max folding time: ' + str(self.getMaxFoldingTime()) + '\n'
        summary += '> Min folding time: ' + str(self.getMinFoldingTime()) + '\n'
        summary += '> Mean folding time: ' + str(self.getMeanFoldingTime()) + '\n'
        summary += '> Median folding time: ' + str(self.getMedianFoldingTime()) + '\n'
        summary += '> Std folding time: ' + str(self.getStdFoldingTime()) + '\n\n'
        try:
            summary += '> DRP is ' + self.getDRP().getName() + ' (' + str(self.getDRP().getPenalityAtFolding()) + ')\n'
            summary += '> DRP (complete bias) is ' + self.getDRPComplete().getName() + \
                       ' (' + str(self.getDRPComplete().getPenalityAtTime(self.getMaxFoldingTime())) + ')\n'
            summary += '> DRP (normalised bias) is ' + self.getDRPNormal().getName() + \
                       ' (' + str(self.getDRPNormal().getPenalityNormalAtFolding()) + ')\n'
        except:
            summary += '> No DRP information available\n'
        
        return summary
