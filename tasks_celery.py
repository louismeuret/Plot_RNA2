from create_plots import *
import json
import kaleido
import plotly
import time
from celery import Celery
import os
import pandas as pd
import barnaba as bb
from FoldingAnalysis.analysis import Trajectory
import energy_3dplot

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

def plotly_to_json(fig):
    return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)

import numpy as np
import mdtraj as md
from itertools import combinations

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
        dotbracket_df = pd.DataFrame(dotbracket_data, columns=["DotBracket"])
        dotbracket_df.to_csv(os.path.join(download_path, "dotbracket_data.csv"), index=False)
        fig = plot_diagram_frequency(sequence, dotbracket_data, dotbracket_native)
        fig.write_html(os.path.join(plot_path, "arc_diagram_plot.html"))
        plotly_data = plotly_to_json(fig)
        return plotly_data
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
def generate_landscape_plot(self, native_pdb_path,traj_xtc_path, download_path, plot_path, session_id, landscape_stride):
    try:
        traj_load = md.load_xtc(traj_xtc_path,native_pdb_path)
        native_load = md.load_pdb(native_pdb_path)
        q_values = best_hummer_q(traj_load, native_load)
        rmsd = bb.rmsd(
            native_pdb_path,
            traj_xtc_path,
            topology=native_pdb_path,
            heavy_atom=True,
        )
        print(len(rmsd))

        stride = int(landscape_stride)
        values = [q_values[x] for x in range(0, len(traj_load), stride)]
        good_rmsd = [rmsd[x] for x in range(0, len(traj_load), stride)]

        df = pd.DataFrame(
            {
                "frame": list(range(0, len(good_rmsd))),
                "Q": values,
                "RMSD": good_rmsd,
                "traj": "traj_1",
            }
        )
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

        path_landscape = f"static/{session_id}/landscape.png"
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
