from create_plots import *
import json
import plotly
import time

def plotly_to_json(fig):
    return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)

def generate_contacts_maps_plot(native_pdb, traj_xtc, download_path, plot_path, session_id):
    create_contact_maps(traj_xtc, native_pdb, 1, session_id)
    return plotly

def generate_rmsd_plot(native_pdb_path, traj_xtc_path, download_path, plot_path):
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

def generate_ermsd_plot(native_pdb_path, traj_xtc_path, download_path, plot_path):
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

def generate_torsion_plot(traj_xtc_path, native_pdb_path, download_path, plot_path, session_id):
    angles, res = bb.backbone_angles(traj_xtc_path, topology=native_pdb_path)
    print(res)
    fig = plot_torsion(angles, res, session["torsionResidue"])
    fig.write_html(os.path.join(plot_path, "torsion_plot.html"))
    plotly_data = plotly_to_json(fig)
    return plotly_data

def generate_sec_structure_plot(traj_xtc_path, native_pdb_path, download_path, plot_path, session_id):
    sec_structure = bb.secondary_structure(
        traj_xtc_path, topology=native_pdb_path
    )
    fig = plot_sec_structure(sec_structure)
    fig.write_html(os.path.join(plot_path, "secondary_structure_plot.html"))
    plotly_data = plotly_to_json(fig)
    return plotly_data

def generate_annotate_plot(traj_xtc_path, native_pdb_path, download_path, plot_path, session_id):
    stackings, pairings, res = bb.annotate(
        traj_xtc_path, topology=native_pdb_path
    )
    plot_data.append([plot, "annotate", stackings, pairings, res])

def generate_dotbracket_plot(dotbracket_data, download_path, plot_path, session_id):
    stackings, pairings, res = bb.annotate(
        traj_xtc_path, topology=native_pdb_path
    )
    dotbracket_data = bb.dot_bracket(pairings, res)[0]
    fig = plot_dotbracket(dotbracket_data)
    fig.write_html(os.path.join(plot_path, "dotbracket_timeline_plot.html"))
    plotly_data = plotly_to_json(fig)
    return plotly_data

def generate_arc_plot(traj_xtc_path, native_pdb_path, download_path, plot_path, session_id):
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



def generate_ds_motif_plot(traj_xtc_path, native_pdb_path, download_path, plot_path, session_id):
    # Implement DS_MOTIF plot
    pass
def generate_ss_motif_plot(traj_xtc_path, native_pdb_path, download_path, plot_path, session_id):
    # Implement SS_MOTIF plot
    pass
def generate_jcoupling_plot(traj_xtc_path, native_pdb_path, download_path, plot_path, session_id):
    couplings, res = bb.jcouplings(traj_xtc_path, topology=native_pdb_path)
    plot_data.append([plot, "jcoupling", couplings])
def generate_escore_plot(traj_xtc_path, native_pdb_path, download_path, plot_path, session_id):
    # Implement ESCORE plot
    pass

def generate_landscape_plot(traj_xtc_path, native_pdb_path, download_path, plot_path, session_id):
    trajectory = Trajectory(
        filename=traj_xtc_path, ref_filename=native_pdb_path
    )
    test = trajectory.q_soft
    values = []
    good_rmsd = []
    rmsd = bb.rmsd(
        native_pdb_path,
        traj_xtc_path,
        topology=native_pdb_path,
        heavy_atom=True,
    )
    print(len(rmsd))

    for x in range(len(trajectory)):
        if x % int(session["landscape_stride"]) == 0:
            print(x)
            values.append(test[x])
            good_rmsd.append(rmsd[x])

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
    # selected_regions=[(0.06,0.17,0.95,1.05),(0.15,0.22,0.78,0.92),(0.2,0.27,0.67,0.77),(0.25,0.33,0.49,0.63)] #list of tuples with minQ, maxQ, minRMSD, maxRMSD FROM CLOSEST TO NATIVE TO UNFOLDED!!!
    selected_regions = (
        []
    )  # list of tuples with minQ, maxQ, minRMSD, maxRMSD FROM CLOSEST TO NATIVE TO UNFOLDED!!!

    dataframe, max_RMSD = df, max(df["RMSD"])
    (
        probability_matrix,
        allframes_matrix,
        Qbin,
        RMSDbin,
    ) = energy_3dplot.make_matrix_probability(dataframe, size, max_RMSD)
    # energy_3dplot.probability_plot(probability_matrix, Qbin, RMSDbin, max_RMSD)
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