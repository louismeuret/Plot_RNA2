import os
import sys
import numpy as np
import mdtraj as md
import energy_3dplot
from FoldingAnalysis.analysis import *
import plotly.graph_objects as go
from string import ascii_uppercase
import barnaba as bb

ref = "frame_0_SHAW.pdb"
traj = "full_centered_SHAW.xtc"

trajectory = Trajectory(filename=traj, ref_filename=ref)

print(len(trajectory))
test = trajectory.q_soft
values = []
good_rmsd = []
rmsd = bb.rmsd(ref,traj,topology=ref, heavy_atom=True)
print(len(rmsd))

for x in range(len(trajectory)):
    if x % 10 == 0:
        print(x)
        values.append(test[x])
        good_rmsd.append(rmsd[x])

import pandas as pd
df = pd.DataFrame({
    'frame': list(range(0,len(good_rmsd))),
    'Q': values,
    'RMSD': good_rmsd,
    'traj': 'traj_1'
})
print("finished dataframe")
#inputs:
size=65
#selected_regions=[(0.06,0.17,0.95,1.05),(0.15,0.22,0.78,0.92),(0.2,0.27,0.67,0.77),(0.25,0.33,0.49,0.63)] #list of tuples with minQ, maxQ, minRMSD, maxRMSD FROM CLOSEST TO NATIVE TO UNFOLDED!!!
selected_regions=[] #list of tuples with minQ, maxQ, minRMSD, maxRMSD FROM CLOSEST TO NATIVE TO UNFOLDED!!!
#selected_regions=[] #list of tuples with minQ, maxQ, minRMSD, maxRMSD FROM CLOSEST TO NATIVE TO UNFOLDED!!!
#main:
#dataframe, max_RMSD=loaddata('concatenated_Q_RMSD')
dataframe, max_RMSD=df, max(df['RMSD'])
probability_matrix, allframes_matrix, Qbin, RMSDbin=energy_3dplot.make_matrix_probability(dataframe, size, max_RMSD)
energy_3dplot.probability_plot(probability_matrix, Qbin, RMSDbin, max_RMSD)
energy_matrix, real_values=energy_3dplot.make_matrix_energy(probability_matrix, max_RMSD, size)
#energy_3dplot.energy_plot(energy_matrix, Qbin, RMSDbin, max_RMSD, real_values, selected_regions)
energy_3dplot.energy_plot_3d(energy_matrix, Qbin, RMSDbin, max_RMSD, real_values, selected_regions)

energy_3dplot.energy_plot(energy_matrix, Qbin, RMSDbin, max_RMSD, real_values, selected_regions)
