import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os

def plot_column(data, time_column, column_to_plot, output_folder):
    plt.figure(figsize=(10, 6))
    plt.plot(data[time_column], data[column_to_plot])
    plt.title(f'{column_to_plot} over time')
    plt.xlabel('Time')
    plt.ylabel(column_to_plot)
    plt.savefig(os.path.join(output_folder, f"{column_to_plot}_over_time.png"))
    plt.close()

def read_biased(csv_file):
    """
    Reads biased data from a CSV file and generates plots for all numerical columns using Matplotlib.

    Parameters:
    - csv_file: Path to the "ratchet.out" file.
    """
    # Load data using NumPy
    data = np.genfromtxt(csv_file, delimiter='\t', names=True, dtype=None, encoding=None)
    
    # Extract column names
    column_names = data.dtype.names
    
    # Exclude the first column (assuming it's time) and plot each of the other columns
    time_column = column_names[0]

    # Create a folder for saving plots in the same directory as ratchet.out
    output_folder = os.path.join(os.path.dirname(csv_file), "plots")
    os.makedirs(output_folder, exist_ok=True)

    for column_to_plot in column_names[1:]:
        plot_column(data, time_column, column_to_plot, output_folder)

# Example usage:
read_biased("/home/cya/Téléchargements/Projet_STAGE/Analysis_rna_GAAA/Compair/ratchet_43_rna/t/ratchet.out")

