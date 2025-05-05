from create_plots import *
import json
import orjson
import time
from celery import Celery
import os
import pandas as pd
import barnaba as bb
from FoldingAnalysis.analysis import Trajectory
import numpy as np
import mdtraj as md
from itertools import combinations
import energy_3dplot
from plotly.io import to_json
import pickle

# Configure Celery
app2 = Celery('tasks')

# Celery configuration
app2.conf.update(
    broker_url='redis://localhost:6379/0',
    result_backend='redis://localhost:6379/0',
    task_serializer='json',
    result_serializer='json',
    accept_content=['json'],
    enable_utc=True,
    task_track_started=True,
    task_ignore_result=False,
    result_expires=None
)

#def plotly_to_json(fig):
#    return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)

def plotly_to_json(fig):
    return to_json(fig, validate=False, engine="orjson")

def best_hummer_q(traj, native):
    """Compute the fraction of native contacts according the definition from
    Best, Hummer and Eaton [1]

    Parameters
    ----------
    traj : md.Trajectory
        The trajectory to do the computation for
    native : md.Trajectory
        The 'native state'. This can be an entire trajecory, or just a single frame.
        Only the first conformation is used

    Returns
    -------
    q : np.array, shape=(len(traj),)
        The fraction of native contacts in each frame of `traj`

    References
    ----------
    ..[1] Best, Hummer, and Eaton, "Native contacts determine protein folding
          mechanisms in atomistic simulations" PNAS (2013)
    """

    BETA_CONST = 50  # 1/nm
    LAMBDA_CONST = 1.8
    NATIVE_CUTOFF = 0.45  # nanometers

    # get the indices of all of the heavy atoms
    heavy = native.topology.select('not element H')
    # get the pairs of heavy atoms which are farther than 3
    # residues apart
    heavy_pairs = np.array(
        [(i,j) for (i,j) in combinations(heavy, 2)
            if abs(native.topology.atom(i).residue.index - \
                   native.topology.atom(j).residue.index) > 3])

    # compute the distances between these pairs in the native state
    heavy_pairs_distances = md.compute_distances(native[0], heavy_pairs)[0]
    # and get the pairs s.t. the distance is less than NATIVE_CUTOFF
    native_contacts = heavy_pairs[heavy_pairs_distances < NATIVE_CUTOFF]
    print("Number of native contacts", len(native_contacts))

    # now compute these distances for the whole trajectory
    r = md.compute_distances(traj, native_contacts)
    # and recompute them for just the native state
    r0 = md.compute_distances(native[0], native_contacts)

    q = np.mean(1.0 / (1 + np.exp(BETA_CONST * (r - LAMBDA_CONST * r0))), axis=1)
    return q

@app2.task(bind=True, max_retries=3)
def generate_contact_maps_plot(self, native_pdb, traj_xtc, download_path, plot_path, session_id):
    try:
        create_contact_maps(traj_xtc, native_pdb, 1, session_id)
        path_save_figure = os.path.join('static', session_id, 'contact_map_plotly.png')
        return json.dumps({'path': path_save_figure})
    except Exception as exc:
        self.retry(exc=exc, countdown=60)  # Retry after 60 seconds

@app2.task(bind=True, max_retries=3)
def generate_rmsd_plot(self, native_pdb_path, traj_xtc_path, download_path, plot_path, session_id):
    try:
        start_time = time.time()

        # Calculate RMSD
        rmsd_start = time.time()
        rmsd = bb.rmsd(
            native_pdb_path,
            traj_xtc_path,
            topology=native_pdb_path,
            heavy_atom=True,
        )
        print(f"RMSD calculation time: {time.time() - rmsd_start:.2f}s")

        # Create plot
        plot_start = time.time()
        fig = plot_rmsd(rmsd)
        print(f"Plot creation time: {time.time() - plot_start:.2f}s")

        # Save data and plot
        save_start = time.time()
        rmsd_df = pd.DataFrame({"RMSD": rmsd})
        rmsd_df.to_csv(os.path.join(download_path, "rmsd_values.csv"), index=False)
        fig.write_html(os.path.join(plot_path, "rmsd_plot.html"))
        print(f"Save time: {time.time() - save_start:.2f}s")

        # Convert to JSON
        json_start = time.time()
        plotly_data = plotly_to_json(fig)
        print(f"JSON conversion time: {time.time() - json_start:.2f}s")

        print(f"Total RMSD plot generation time: {time.time() - start_time:.2f}s")
        return plotly_data
    except Exception as exc:
        self.retry(exc=exc, countdown=60)

@app2.task(bind=True, max_retries=3)
def generate_ermsd_plot(self, native_pdb_path, traj_xtc_path, download_path, plot_path, session_id):
    try:
        start_time = time.time()

        # Calculate ERMSD
        ermsd_start = time.time()
        ermsd = bb.ermsd(native_pdb_path, traj_xtc_path, topology=native_pdb_path)
        print(f"ERMSD calculation time: {time.time() - ermsd_start:.2f}s")

        # Create plot
        plot_start = time.time()
        fig = plot_ermsd(ermsd)
        print(f"Plot creation time: {time.time() - plot_start:.2f}s")

        # Save data and plot
        save_start = time.time()
        ermsd_df = pd.DataFrame({"ERMSD": ermsd})
        ermsd_df.to_csv(os.path.join(download_path, "ermsd_values.csv"), index=False)
        fig.write_html(os.path.join(plot_path, "ermsd_plot.html"))
        print(f"Save time: {time.time() - save_start:.2f}s")

        # Convert to JSON
        json_start = time.time()
        plotly_data = plotly_to_json(fig)
        print(f"JSON conversion time: {time.time() - json_start:.2f}s")

        print(f"Total ERMSD plot generation time: {time.time() - start_time:.2f}s")
        return plotly_data
    except Exception as exc:
        self.retry(exc=exc, countdown=60)

@app2.task(bind=True, max_retries=3)
def generate_torsion_plot(self, traj_xtc_path, native_pdb_path, download_path, plot_path, session_id):
    try:
        angles, res = bb.backbone_angles(traj_xtc_path, topology=native_pdb_path)
        print(res)
        fig = plot_torsion(angles, res, session["torsionResidue"])
        fig.write_html(os.path.join(plot_path, "torsion_plot.html"))
        plotly_data = plotly_to_json(fig)
        return plotly_data
    except Exception as exc:
        self.retry(exc=exc, countdown=60)

@app2.task(bind=True, max_retries=3)
def generate_sec_structure_plot(self, native_pdb_path,traj_xtc_path, download_path, plot_path, session_id):
    try:
        stackings, pairings, res = bb.annotate(
            traj_xtc_path, topology=native_pdb_path
        )
        dotbracket_data, res2 = bb.dot_bracket(pairings, res)
        return [dotbracket_data, res2.strip()]
    except Exception as exc:
        self.retry(exc=exc, countdown=60)

@app2.task(bind=True, max_retries=3)
def generate_2Dpairing_plot(self, native_pdb_path, traj_xtc_path, download_path, plot_path, session_id):
    try:
        rvecs_traj, res_traj = bb.dump_rvec(traj_xtc_path, topology=native_pdb_path, cutoff=100.0)
        nonzero = np.where(np.sum(rvecs_traj**2, axis=3) > 0.01)
        rr = rvecs_traj[nonzero]
        fig = base_pairs_visualisation(rr)
        fig.write_html(os.path.join(plot_path, "2D_pairing.html"))
        plotly_data = plotly_to_json(fig)
        return plotly_data

    except Exception as exc:
        self.retry(exc=exc, countdown=60)

@app2.task(bind=True, max_retries=3)
def generate_annotate_plot(self, traj_xtc_path, native_pdb_path, download_path, plot_path, session_id):
    try:
        stackings, pairings, res = bb.annotate(
            traj_xtc_path, topology=native_pdb_path
        )
        return [plot, "annotate", stackings, pairings, res]
    except Exception as exc:
        self.retry(exc=exc, countdown=60)

@app2.task(bind=True, max_retries=3)
def generate_dotbracket_plot(self, native_pdb_path, traj_xtc_path, download_path, plot_path, session_id):
    try:
        stackings, pairings, res = bb.annotate(
            traj_xtc_path, topology=native_pdb_path
        )
        dotbracket_data = bb.dot_bracket(pairings, res)[0]
        fig = plot_dotbracket(dotbracket_data)
        fig.write_html(os.path.join(plot_path, "dotbracket_timeline_plot.html"))
        plotly_data = plotly_to_json(fig)
        return plotly_data
    except Exception as exc:
        self.retry(exc=exc, countdown=60)

@app2.task(bind=True, max_retries=3)
def generate_arc_plot(self, native_pdb_path, traj_xtc_path, download_path, plot_path, session_id):
    try:
        print("ARC")
        stackings, pairings, res = bb.annotate(
            traj_xtc_path, topology=native_pdb_path
        )
        dotbracket_data = bb.dot_bracket(pairings, res)[0]
        #Dotbracket for the native state:
        stackings, pairings, res = bb.annotate(native_pdb_path)
        dotbracket_native = bb.dot_bracket(pairings, res)[0]
        print(f"DOTBRACKET NATIVE = {str(dotbracket_native[0])}")

        sequence = ''.join([item[0] for item in res])
        resids = [item.split("_")[1] for item in res]
        resids_json = orjson.dumps(resids)
        print(f"ARC SEQUENCE = {str(sequence)}")
        dotbracket_df = pd.DataFrame(dotbracket_data, columns=["DotBracket"])
        dotbracket_df.to_csv(os.path.join(download_path, "dotbracket_data.csv"), index=False)
        fig = plot_diagram_frequency(sequence, dotbracket_data, dotbracket_native)
        print(f"ARC resids = {resids_json}")
        fig.write_html(os.path.join(plot_path, "arc_diagram_plot.html"))
        plotly_data = plotly_to_json(fig)
        return [plotly_data, resids]
    except Exception as exc:
        self.retry(exc=exc, countdown=60)

@app2.task(bind=True, max_retries=3)
def generate_ds_motif_plot(self, traj_xtc_path, native_pdb_path, download_path, plot_path, session_id):
    # Implement DS_MOTIF plot
    pass

@app2.task(bind=True, max_retries=3)
def generate_ss_motif_plot(self, traj_xtc_path, native_pdb_path, download_path, plot_path, session_id):
    # Implement SS_MOTIF plot
    pass

@app2.task(bind=True, max_retries=3)
def generate_jcoupling_plot(self, traj_xtc_path, native_pdb_path, download_path, plot_path, session_id):
    try:
        couplings, res = bb.jcouplings(traj_xtc_path, topology=native_pdb_path)
        return [plot, "jcoupling", couplings]
    except Exception as exc:
        self.retry(exc=exc, countdown=60)

@app2.task(bind=True, max_retries=3)
def generate_escore_plot(self, traj_xtc_path, native_pdb_path, download_path, plot_path, session_id):
    # Implement ESCORE plot
    pass

@app2.task(bind=True, max_retries=3)
def update_contact_map_plot(self, generate_data_path, plot_path,  frame_number, session_id):

    with open(os.path.join(generate_data_path, "contact_map_data.pkl"), 'rb') as f:
        loaded_data = pickle.load(f)

    frames_dict = loaded_data['frames_dict']
    sequence = loaded_data['sequence']

    frame_num, frame_data = sorted(frames_dict.items())[frame_number]
    print(f"EXPECTED DIRECTORY TO SAVE PLOTS = {os.path.join(plot_path, 'contact_map.html')}")
    fig = plot_rna_contact_map(frame_data, sequence, output_file=os.path.join(plot_path, "contact_map.html"), frame_number=frame_num)
    plotly_data = plotly_to_json(fig)
    return plotly_data


@app2.task(bind=True, max_retries=3)
def update_landscape_frame(self, generate_data_path, coordinates):
    with open(os.path.join(generate_data_path, "dataframe.pkl"), 'rb') as f:
        loaded_data = pickle.load(f)
    df = loaded_data
    target_Q = coordinates['Q']
    target_RMSD = coordinates['RMSD']

    # Calculate Euclidean distance between target coordinates and all points
    distances = np.sqrt((df['Q'] - target_Q)**2 + (df['RMSD'] - target_RMSD)**2)
    closest_frame_idx = distances.idxmin()
    closest_frame = int(df.iloc[closest_frame_idx]['frame'])  # Convert to Python int
    return closest_frame



@app2.task(bind=True, max_retries=3)
def generate_contact_map_plot(self, native_pdb_path, traj_xtc_path, download_path, plot_path, generate_data_path, session_id):
    def process_barnaba_pairings(pairings, res):
        """
        Process barnaba pairings and residue information into frames dictionary

        Parameters:
        -----------
        pairings : list
            List of pairs for each frame from barnaba.annotate
        res : list
            List of residue information

        Returns:
        --------
        frames_dict : dict
            Dictionary with frame numbers as keys and DataFrames with base pair information as values
        sequence : list
            RNA sequence
        """
        frames_dict = {}

        # Create sequence from residue information
        sequence = [r[0] for r in res]  # Assuming res contains nucleotide types

        # Process each frame
        for frame_num, frame_data in enumerate(pairings):
            base_pairs = []

            # Each frame contains a list of pairs and their annotations
            if len(frame_data) == 2:
                pair_indices = frame_data[0]
                annotations = frame_data[1]

                for pair_idx, pair in enumerate(pair_indices):
                    if not pair:
                        continue

                    res_i = pair[0] + 1  # Convert to 1-based indexing
                    res_j = pair[1] + 1

                    # Get residue names from the sequence
                    if 0 <= res_i - 1 < len(sequence) and 0 <= res_j - 1 < len(sequence):
                        res_i_name = f"{sequence[res_i - 1]}{res_i}"
                        res_j_name = f"{sequence[res_j - 1]}{res_j}"
                    else:
                        res_i_name = f"N{res_i}"
                        res_j_name = f"N{res_j}"

                    anno = annotations[pair_idx] if pair_idx < len(annotations) else 'XXX'

                    base_pairs.append({
                        'res_i': res_i,
                        'res_j': res_j,
                        'res_i_name': res_i_name,
                        'res_j_name': res_j_name,
                        'anno': anno
                    })

            # Always assign a DataFrame (even if empty)
            frames_dict[frame_num] = pd.DataFrame(base_pairs)

        return frames_dict, sequence

    # Process barnaba results into our format
    stackings, pairings, res = stackings, pairings, res = bb.annotate(traj_xtc_path, topology=native_pdb_path)
    frames_dict, sequence = process_barnaba_pairings(pairings, res)
    print(len(frames_dict))
    print(f"RNA length: {len(sequence)} nucleotides")
    print(f"Found {len(frames_dict)} frames in the data")
    frame_num, frame_data = sorted(frames_dict.items())[0]

    # Save pairings to CSV
    pairings_df = pd.DataFrame(pairings)
    pairings_df.to_csv(os.path.join(download_path, "pairings.csv"), index=False)

    # Save frames_dict and sequence to pickle
    data_to_save = {
        'frames_dict': frames_dict,
        'sequence': sequence
    }
    with open(os.path.join(generate_data_path, "contact_map_data.pkl"), 'wb') as f:
        pickle.dump(data_to_save, f)

    fig = plot_rna_contact_map(frame_data, sequence, output_file=os.path.join(plot_path, "contact_map.html"), frame_number=frame_num)
    plotly_data = plotly_to_json(fig)
    return plotly_data


@app2.task(bind=True, max_retries=3)
def generate_landscape_plot(self, native_pdb_path,traj_xtc_path, download_path, plot_path, session_id, parameters, generate_data_path):
    def calculate_structural_metrics(parameters, traj_load, native_load, native_pdb_path, traj_xtc_path):
        """
        Calculate structural metrics for trajectory analysis based on specified parameters.

        Args:
            parameters (list): List of metrics to calculate
            traj_load: Loaded trajectory data
            native_load: Loaded native structure data
            native_pdb_path (str): Path to native PDB file
            traj_xtc_path (str): Path to trajectory XTC file

        Returns:
            dict: Dictionary containing calculated metrics
        """
        # Define a dictionary mapping calculation types to their functions
        calculation_functions = {
            "Q": lambda: best_hummer_q(traj_load, native_load),
            "RMSD": lambda: bb.rmsd(
                native_pdb_path,
                traj_xtc_path,
                topology=native_pdb_path,
                heavy_atom=True,
            ),
            "eRMSD": lambda: bb.ermsd(
                native_pdb_path,
                traj_xtc_path,
                topology=native_pdb_path
            ),
            "Torsion": lambda: print('torsion')
        }

        # Initialize results dictionary
        results = {}

        # Calculate requested metrics
        for metric in parameters:
            if metric in calculation_functions:
                results[metric] = calculation_functions[metric]()
            else:
                results[metric] = None
                print(f"Warning: Unknown metric '{metric}'. Skipping calculation.")

        return results

    try:
        traj_load = md.load_xtc(traj_xtc_path,native_pdb_path)
        native_load = md.load_pdb(native_pdb_path)
        print(parameters)
        metrics_to_calculate = [parameters[1],parameters[2]]

        results = calculate_structural_metrics(
                metrics_to_calculate,
                traj_load,
                native_load,
                native_pdb_path,
                traj_xtc_path
        )
        print(parameters[0])
        stride = int(parameters[0])
        component1 = [results[parameters[1]][x] for x in range(0, len(traj_load), stride)]
        component2 = [results[parameters[2]][x] for x in range(0, len(traj_load), stride)]
        print(results)
        #component1 = results[parameters[1]][::stride]
        #component2 = results[parameters[2]][::stride]


        df = pd.DataFrame(
            {
                "frame": list(range(0, len(component2))),
                "Q": component1,
                "RMSD": component2,
                "traj": "traj_1",
            }
        )
        with open(os.path.join(generate_data_path, "dataframe.pkl"), 'wb') as f:
            pickle.dump(df, f)
        print("finished dataframe")
        size = 65
        selected_regions = []

        dataframe, max_RMSD = df, max(df["RMSD"])
        (
            probability_matrix,
            allframes_matrix,
            Qbin,
            RMSDbin,
        ) = energy_3dplot.make_matrix_probability(dataframe, size, max_RMSD)
        energy_matrix, real_values = energy_3dplot.make_matrix_energy(
            probability_matrix, max_RMSD, size
        )
        fig = plot_landscapes_3D(
            energy_matrix, Qbin, RMSDbin, max_RMSD, real_values, selected_regions
        )
        path_landscape_3d = f"static/{session_id}/download_plot/LANDSCAPE/landscape.html"
        fig.write_html(path_landscape_3d)
        path_landscape = f"static/{session_id}/download_plot/LANDSCAPE/landscape.png"
        energy_3dplot.energy_plot(
            energy_matrix,
            Qbin,
            RMSDbin,
            max_RMSD,
            real_values,
            selected_regions,
            path_landscape,
        )
        plotly_data = plotly_to_json(fig)
        return plotly_data
    except Exception as exc:
        self.retry(exc=exc, countdown=60)
