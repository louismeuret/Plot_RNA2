import numpy as np 
import matplotlib.pyplot as plt
import matplotlib as mlp
import pandas as pd
import scipy.stats as stats
import os
import math
import plotly.graph_objects as go
import plotly.express as px
from matplotlib.colors import LinearSegmentedColormap

#NB Kb*T=kT = 4.11×10−21 J according to Wikipedia

#loading data
def loaddata(filename):
	df = pd.read_csv (filename, sep='\t')
	df.columns = ["frame", "Q", 'RMSD', 'traj']
	return (df, max(df['RMSD']))

#makes matrix of empty lists
def make_emptylists_matrix(dimension):
	matrix_empty=np.empty((dimension,dimension), dtype=np.object_) #NB important to specify type of values otherwise numpy can't understand
	for i in range(matrix_empty.shape[0]):
		for j in range(matrix_empty.shape[1]):
			matrix_empty[i,j]=[]
	return(matrix_empty)

#I want a matrix (24x24) to keep bins consistent with previous 2D plots

def make_matrix_probability(df, dimension, maximal_RMSD):
	#generate our matrix
	matrix_counts=np.zeros((dimension,dimension))
	matrix_txt=make_emptylists_matrix((dimension))
	bin_size_Q=float(1/dimension)
	bin_size_RMSD=float(maximal_RMSD/dimension)
	for i in range(len(df)):
		q_normalized=int(df['Q'][i]/bin_size_Q)
		rmsd_normalized=int(df['RMSD'][i]/bin_size_RMSD)
		if q_normalized and rmsd_normalized:	
			matrix_counts[q_normalized-1][rmsd_normalized-1]+=1
			matrix_txt[q_normalized-1][rmsd_normalized-1].append([df['traj'][i],df['frame'][i]])
	matrix_probability=np.divide(matrix_counts, (len(df['frame']))) 
	return (matrix_probability, matrix_txt, bin_size_Q, bin_size_RMSD)

def probability_plot(matrix_probability, bin_size_Q, bin_size_RMSD, maximal_RMSD):
	#initialize...
	fig, ax= plt.subplots(figsize=(10, 10))
	#plan graph...
	ax.set_aspect('equal')
	ax.set_xlim([0, 1])
	ax.set_ylim([0, maximal_RMSD])
	xlocs, xticks=plt.xticks(np.arange(0, 1.1, 0.1), rotation=90)
	xticks[0].set_visible(False) #removes tick at 0 so the plot is less cluttered
	ylocs, yticks=plt.yticks(np.arange(0, maximal_RMSD, 1))
	yticks[0].set_visible(False)
	plt.tick_params (axis = "both", which = "both", bottom = False, top = False, left=False, right=False)
	plt.xlabel('Q', fontsize=14)
	plt.ylabel('RMSD', fontsize=14)
	#generate the plot
	p=plt.imshow(matrix_probability.T, origin = "lower", interpolation = "gaussian", extent=[0,1,0,maximal_RMSD], aspect=float(bin_size_Q/bin_size_RMSD), cmap='Blues') #set aspect to ratio of x unity/y unity to get square plot
	cbar=plt.colorbar(p, ax=ax, ticks=np.arange(0,np.amax(matrix_probability), 0.02), shrink=0.805)
	cbar.set_label('Probability', fontsize=14, rotation=90)
#	plt.show()

def energy_plot_3d(matrix_energy_rescaled, bin_size_Q, bin_size_RMSD, maximal_RMSD, real_values, squares):
    # Calculate the bin centers from edges for plotting
    x_centers = np.linspace(0, 1, matrix_energy_rescaled.shape[0])
    y_centers = np.linspace(0, maximal_RMSD, matrix_energy_rescaled.shape[1])

    # Create a meshgrid for the x and y bin centers
    X, Y = np.meshgrid(x_centers, y_centers)

    # Create the surface plot
    fig = go.Figure(data=[go.Surface(z=matrix_energy_rescaled.T, x=X, y=Y, colorscale='RdYlBu_r', cmin=0, cmax=np.amax(real_values))])

    # Add rectangles for the squares (as 3D boxes)
    colours = dict(zip(squares, px.colors.qualitative.T10[:len(squares)]))
    labels = dict(zip(squares, [ascii_uppercase[i] for i in range(len(squares))]))
    
    for square in squares:
        x0, x1, y0, y1 = square
        fig.add_trace(go.Scatter3d(
            x=[x0, x1, x1, x0, x0],
            y=[y0, y0, y1, y1, y0],
            z=[np.amax(real_values)]*5,
            mode='lines',
            line=dict(color=colours[square], width=10),
            name=labels[square]
        ))

    # Update layout
    fig.update_layout(
        title='3D Energy Plot',
        scene=dict(
            xaxis_title='Q (Froebenius distance)',
            yaxis_title='RMSD (nm)',
            zaxis_title='-ln(p)  (Frequency)',
            zaxis=dict(range=[0, np.amax(real_values)])
        ),
        width=800,
        height=800
    )

    return fig
	

def make_matrix_energy(matrix_probability, maximal, dimension):
	#convert probability into energy
	matrix_energy=np.where(matrix_probability>0, -np.log(matrix_probability), 100) #np.where includes a condition to use first vs second value to build new array in corresponding position
	#rescale so that min value is 0
	matrix_energy_rescaled=np.where(matrix_energy!=100, matrix_energy-np.amin(matrix_energy), 100)
	real_values=np.where(matrix_energy_rescaled!=100, matrix_energy_rescaled, -100)
	matrix_energy_rescaled=np.where(matrix_energy!=100, matrix_energy_rescaled, np.amax(real_values)+1)
	return(matrix_energy_rescaled, real_values)

def energy_plot(matrix_energy_rescaled, bin_size_Q, bin_size_RMSD, maximal_RMSD, real_values, squares, path_landscape):
	#now let's make the figure:
	#initialize...
	colors = ["#E86A44", "#4E4182"]
	custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)
	fig, ax= plt.subplots(figsize=(10, 10))
	ax.set_title('-ln(p) (frequency)')
	colours = dict(zip(squares, plt.cm.tab10.colors[:len(squares)]))
	from string import ascii_uppercase
	labels=dict(zip(squares,[ascii_uppercase[i] for i in range(len(squares))]))
	#plan graph...
	ax.set_aspect('equal')
	ax.set_xlim([0, 1])
	ax.set_ylim([0, maximal_RMSD])
	ax.autoscale(enable=True, axis='y',tight=True)
	xlocs, xticks=plt.xticks(np.arange(0, 1.1, 0.1), rotation=90)
	xticks[0].set_visible(False)
	ylocs, yticks=plt.yticks(np.arange(0, maximal_RMSD, 0.1))
	yticks[0].set_visible(False)
	plt.tick_params (axis = "both", which = "both", bottom = True, top = False, left=True, right=False)
	plt.xlabel('Q (Froebenius distance)', fontsize=14)
	plt.ylabel('RMSD (nm)', fontsize=14)
	#generate the plot
	#cmap = plt.cm.get_cmap('RdYlBu') #grab standard colormap
	cmap=custom_cmap.reversed() #inverts color range to get blue for lower energy as convention!!
	cmap.set_over('white') #needed to set all values above threshold=unreal values to white
	p=plt.imshow(matrix_energy_rescaled.T, origin='lower', extent=[0,1,0.01,maximal_RMSD], interpolation='gaussian', aspect=float(bin_size_Q/bin_size_RMSD), cmap=cmap, vmax=np.amax(real_values)) #set aspect to ratio of x unity/y unity to get square plot
	#cbar=plt.colorbar(p, ax=ax, ticks=np.arange(np.amin(matrix_energy_rescaled.T), np.amax(real_values), 1.0), shrink=0.82)
	cbar=plt.colorbar(p, ax=ax, ticks=np.arange(1, np.amax(real_values), 1.0), shrink=0.805)
	plt.clim(2.5, np.amax(real_values)) #sets range for colorbar
	cbar.set_label('-ln(p)', fontsize=14, rotation=90)
	#optional grid to help select ranges of Q and RMSD:
	plt.grid()
	#optional to show regin where you extracted intermediates:
	for tupl in squares:
		from matplotlib.patches import Rectangle
		ax.add_patch(Rectangle((tupl[0], tupl[2]), float(tupl[1]-tupl[0]), float(tupl[3]-tupl[2]), edgecolor=colours[tupl], facecolor="None" ,label=labels[tupl], linewidth=1.5 ))
	ax.legend(loc='upper right', framealpha=1, edgecolor='white')
	#plt.show()
	plt.savefig(path_landscape)


#inputs:
size=65
#selected_regions=[(0.06,0.17,0.95,1.05),(0.15,0.22,0.78,0.92),(0.2,0.27,0.67,0.77),(0.25,0.33,0.49,0.63)] #list of tuples with minQ, maxQ, minRMSD, maxRMSD FROM CLOSEST TO NATIVE TO UNFOLDED!!!
selected_regions=[] #list of tuples with minQ, maxQ, minRMSD, maxRMSD FROM CLOSEST TO NATIVE TO UNFOLDED!!!
#selected_regions=[] #list of tuples with minQ, maxQ, minRMSD, maxRMSD FROM CLOSEST TO NATIVE TO UNFOLDED!!!
#main:
#dataframe, max_RMSD=loaddata('concatenated_Q_RMSD')
#dataframe, max_RMSD=loaddata('concat_Q_RMSD_rescaled')
#probability_matrix, allframes_matrix, Qbin, RMSDbin=make_matrix_probability(dataframe, size, max_RMSD)
#probability_plot(probability_matrix, Qbin, RMSDbin, max_RMSD)
#energy_matrix, real_values=make_matrix_energy(probability_matrix, max_RMSD, size)
#energy_plot(energy_matrix, Qbin, RMSDbin, max_RMSD, real_values, selected_regions)







