import os
import mdtraj as md
import plotly.graph_objects as go
import plotly.express as px
import numpy as np
import pandas as pd
import barnaba as bb
import matplotlib.pyplot as plt
from plotly.subplots import make_subplots
import energy_3dplot
from collections import defaultdict
import MDAnalysis as mda
import generate_contact_maps
import logging

# Configure logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def plot_ermsd(ermsd):
    fig = go.Figure(
        data=go.Scattergl(y=ermsd, mode="markers", legendrank=1),
        layout=dict(
            title="ERMSD Scatter Plot",
            xaxis_title="Frame",
            yaxis_title="ERMSD Value",
            hovermode="closest"
        )
    )
    return fig

def plot_rmsd(rmsd):
    fig = go.Figure(
        data=go.Scattergl(y=rmsd, mode="markers", legendrank=1),
        layout=dict(
            title="RMSD Scatter Plot",
            xaxis_title="Frame",
            yaxis_title="RMSD Value",
            hovermode="closest"
        )
    )
    return fig

def plot_dotbracket(dotbracket_data):
    reverse_mapping = {}
    for frame, dotbracket in enumerate(dotbracket_data, start=1):
        reverse_mapping.setdefault(dotbracket, []).append(frame)

    structures_data = []
    for dotbracket, frames in reverse_mapping.items():
        for frame in frames:
            structures_data.append({"Frame": frame, "DotBracket": dotbracket})

    structures_df = pd.DataFrame(structures_data)
    structures_df.sort_values(by=["DotBracket", "Frame"], inplace=True)

    unique_structures = structures_df["DotBracket"].unique()
    structure_colors = [
        f"rgb{tuple(int(255 * x) for x in plt.cm.tab20b(i)[:3])}"
        for i in np.linspace(0, 1, len(unique_structures))
    ]
    color_dict = dict(zip(unique_structures, structure_colors))

    traces = []
    for i, dotbracket in enumerate(unique_structures):
        structure_df = structures_df[structures_df["DotBracket"] == dotbracket]
        x_values = []
        y_values = []
        prev_frame = None
        for _, row in structure_df.iterrows():
            frame = row["Frame"]
            if prev_frame is not None and frame != prev_frame + 1:
                x_values.append(None)
                y_values.append(None)
            x_values.append(frame)
            y_values.append(i + 1 - 0.2)
            prev_frame = frame
        trace = go.Scattergl(
            x=x_values,
            y=y_values,
            mode="lines",
            line=dict(color=color_dict[dotbracket], width=8),
            name=dotbracket,
        )
        traces.append(trace)

    layout = go.Layout(
        title="Timeline of RNA Structures",
        xaxis=dict(title="Frame", showgrid=False, zeroline=False, showline=False),
        yaxis=dict(
            title="Dot-Bracket Structure",
            tickmode="array",
            tickvals=list(range(1, len(unique_structures) + 1)),
            ticktext=unique_structures,
            showgrid=False,
            zeroline=False,
            showline=False,
            ticklen=1,
        ),
        showlegend=False,
        plot_bgcolor="rgba(255, 255, 255, 0)",
    )

    fig = go.Figure(data=traces, layout=layout)
    return fig

def plot_torsion(angles, res, torsionResidue):
    if torsionResidue.isdigit():
        residue_index = int(torsionResidue)
    else:
        residue_index = res.index(torsionResidue)

    residue_index = 2
    specific_residue_angles = angles[:, residue_index, :]  # All frames for residue 2
    angles_names = ["alpha", "beta", "gamma", "delta", "epsilon", "zeta", "chi"]

    # Define colors for each angle in hexadecimal
    colors = [
        "#636efa",
        "#ef553b",
        "#00cc96",
        "#a45bf9",
        "#fa9f5e",
        "#1ad0f2",
        "#ff6692",
    ]

    # Create subplots with 7 rows and 2 columns
    names = [
        f"Angle: {angles_names[i]} - {plot_type}"
        for i in range(7)
        for plot_type in ["Line Plot", "Histogram"]
    ]
    fig = make_subplots(
        rows=7,
        cols=2,
        shared_xaxes=False,
        subplot_titles=names,
        horizontal_spacing=0.1,  # Adjust spacing between plots
        vertical_spacing=0.1,
    )

    # Add traces for each angle (line plots in the first column and histograms in the second column)
    for i in range(7):
        y_values = specific_residue_angles[:, i]
        x_values = np.arange(len(y_values))
        valid_indices = ~np.isnan(y_values)
        color = colors[i]

        # Line plot
        fig.add_trace(
            go.Scattergl(
                x=x_values[valid_indices],
                y=y_values[valid_indices],
                mode="lines",
                name=f"{angles_names[i]} - Line Plot",
                line=dict(color=color),
            ),
            row=i + 1,
            col=1,
        )

        # Histogram
        fig.add_trace(
            go.Histogram(
                x=y_values[valid_indices],
                name=f"{angles_names[i]} - Histogram",
                marker=dict(color=color),
            ),
            row=i + 1,
            col=2,
        )

    # Update layout
    fig.update_layout(
        height=1400,
        width=1600,
        title_text=f"Angles for Residue {res[residue_index]} Across All Frames",
        showlegend=False,
    )

    return fig

def plot_sec_structure(sec_structure):
    fig = go.Figure(data=go.Scattergl(y=sec_structure, mode="markers"))
    fig.update_layout(
        title="Secondary Structure Plot",
        xaxis_title="Frame",
        yaxis_title="Secondary Structure",
    )
    return fig

def plot_landscapes_3D(
    energy_matrix, Qbin, RMSDbin, max_RMSD, real_values, selected_regions
):
    fig = energy_3dplot.energy_plot_3d(
        energy_matrix, Qbin, RMSDbin, max_RMSD, real_values, selected_regions
    )
    return fig

def plot_diagram_frequency(sequence, dot_bracket_list, dotbracket_native):
    def parse_dot_bracket(dot_bracket):
        stack = []
        pairs = []
        logging.debug(f"Parsing dot-bracket: {dot_bracket}")

        for i, char in enumerate(dot_bracket):
            if char == '(':
                stack.append(i)
            elif char == ')':
                if stack:
                    pairs.append((stack.pop(), i))

        logging.debug(f"Pairs found by parser: {pairs}")
        return pairs

    # Function to calculate pair frequencies from a list of dot-bracket structures
    def calculate_pair_frequencies(dot_bracket_list):
        pair_counts = defaultdict(int)
        total_pairs = 0

        for dot_bracket in dot_bracket_list:
            pairs = parse_dot_bracket(dot_bracket)
            for pair in pairs:
                pair_counts[pair] += 1
                total_pairs += 1

        pair_frequencies = {pair: count / total_pairs for pair, count in pair_counts.items()}
        logging.debug(f"Pair frequencies: {pair_frequencies}")
        return pair_frequencies

    # Function to map nucleotide to color
    def nucleotide_color(nucleotide):
        color_map = {
            'A': 'red',
            'U': 'green',
            'G': 'blue',
            'C': 'orange'
        }
        return color_map.get(nucleotide, 'black')  # Default to black if unknown nucleotide

    # Function to plot arc diagram
    def plot_arc_diagram(sequence, pair_frequencies, pairs):
        logging.debug(f"Pairs for arc diagram: {pairs}")

        fig = go.Figure()

        # Plot nodes with nucleotide-specific colors
        for i, nucleotide in enumerate(sequence):
            color = nucleotide_color(nucleotide)
            fig.add_trace(go.Scattergl(
                x=[i], y=[0], mode='markers+text', text=nucleotide, textposition="bottom center",
                marker=dict(size=10, color=color),
                showlegend=False
            ))

        # Get the frequency for each pair
        pairs_unique = list(pair_frequencies.keys())
        frequencies = [pair_frequencies.get(pair, 0) for pair in pairs_unique]

        min_frequency = min(frequencies)
        max_frequency = max(frequencies)

        if frequencies:
            norm_frequencies = [(f - min(frequencies)) / (max(frequencies) - min(frequencies)) for f in frequencies]
        else:
            norm_frequencies = [0] * len(pairs)

        # Custom color scale from blue (least frequent) to red (most frequent)
        custom_colorscale = [
            [0.0, 'rgb(0, 0, 255)'],   # Blue
            [1.0, 'rgb(255, 0, 0)']    # Red
        ]

        # Plot arcs for edges with frequency-based colors
        for (start, end), norm_frequency, frequency in zip(pairs_unique, norm_frequencies, frequencies):
            center = (start + end) / 2
            radius = (end - start) / 2
            theta = np.linspace(0, np.pi, 100)
            x = center + radius * np.cos(theta)
            y = radius * np.sin(theta)
            color = px.colors.sample_colorscale(custom_colorscale, norm_frequency)
            fig.add_trace(go.Scattergl(
                x=x, y=y, mode='lines',
                line=dict(color=color[0], width=2),
                showlegend=False,
                hoverinfo='text',
                text=f'Frequency: {frequency:.2f}'
            ))

        fig.add_trace(go.Scattergl(
            x=[None], y=[None],
            mode='markers',
            marker=dict(
                colorscale=custom_colorscale,
                cmin=min_frequency, cmax=max_frequency,
                colorbar=dict(
                    title="Frequency",
                    tickvals=np.linspace(min_frequency, max_frequency, num=5),
                    ticktext=[f'{v:.2f}' for v in np.linspace(min_frequency, max_frequency, num=5)],
                    x=1.1,  # Position color bar to the right of the plot
                    y=0.5,
                    len=0.75,
                    thickness=20
                ),
                showscale=True
            ),
            hoverinfo='none'
        ))

        # Update layout
        fig.update_layout(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            showlegend=False,
            width=800,
            height=400,
            plot_bgcolor='white',
            title="RNA Secondary Structure Arc Diagram with Frequency-Based Coloring"
        )
        return fig

    logging.debug(f"Dot-bracket list: {dot_bracket_list}")
    # Calculate pair frequencies from the dot-bracket list
    pair_frequencies = calculate_pair_frequencies(dot_bracket_list)
    logging.debug(f"Final pair frequencies: {pair_frequencies}")

    pairs = parse_dot_bracket(dotbracket_native[0])
    logging.debug(f"Pairs found in native structure: {pairs}")
    # Plot arc diagram
    fig = plot_arc_diagram(sequence, pair_frequencies, pairs)
    return fig

def base_pairs_visualisation(rr):
    z = rr[:, 2]
    rho = np.sqrt(rr[:, 0]**2 + rr[:, 1]**2)
    pairs = rr[np.where(np.abs(z) < 0.18)]
    # in case it's having to much issues chuggling the data
    sampled_pairs = pairs[::]

    # Define x and y for the contour plot
    x = sampled_pairs[:, 0]
    y = sampled_pairs[:, 1]

    # Create grid and histogram
    grid_size = 100
    x_edges = np.linspace(-1.1, 1.1, grid_size)
    y_edges = np.linspace(-1.1, 1.1, grid_size)
    heatmap, x_edges, y_edges = np.histogram2d(x, y, bins=[x_edges, y_edges])
    logging.debug(f"Size x and y = {len(x), len(y), x, y}")
    if len(x) > 100000:
        max_threshold = np.percentile(heatmap, 98)  # Clip at the 98th percentile to highlight middle/high-density
        min_threshold = 5
    else:
        min_threshold = 0
        max_threshold = np.percentile(heatmap, 100)  # Clip at the 100th percentile to highlight all values
    heatmap = np.clip(heatmap, min_threshold, max_threshold)  # Clip the histogram data

    # Create the contour plot using the clipped data
    # Ensure that size is positive and within an appropriate range
    size_value = (max_threshold - min_threshold) / 20
    if size_value <= 0:
        size_value = 0.1  # Set a default size if the calculated value is not positive

    fig = go.Figure(data=go.Contour(
        z=heatmap.T,
        x=x_edges,
        y=y_edges,
        colorscale="OrRd",
        contours=dict(
            start=min_threshold,  # Start at the minimum threshold
            end=max_threshold,  # End at the maximum threshold
            size=size_value,  # Use validated size
            coloring="heatmap",
            showlabels=True,
            labelfont=dict(size=12, color="white")
        ),
        colorbar=dict(
            title="Density",
            titleside="right",
            titlefont=dict(size=14, family="Arial, sans-serif")
        )
    ))

    # Add concentric circles and angular lines
    circle_radii = [0.5, 0.75, 1.0]
    for r in circle_radii:
        fig.add_shape(type="circle",
                      xref="x", yref="y",
                      x0=-r, y0=-r, x1=r, y1=r,
                      line=dict(color="black", width=1, dash="dash"))

    for r, label in zip(circle_radii, ["r=0.5 nm", "r=0.75 nm", "r=1.0 nm"]):
        fig.add_annotation(x=-r, y=0, text=label, showarrow=False, font=dict(size=13), xshift=-10)

    # Add region labels
    region_labels = {
        (0.35, 0.45): "Watson-Crick",
        (0.1, 0.6): "GU",
        (0.5, 0.3): "GU",
        (-0.6, 0.4): "Hoogsteen",
        (0.6, -0.4): "Sugar"
    }

    for (x, y), text in region_labels.items():
        fig.add_annotation(x=x, y=y, text=text, showarrow=False, font=dict(size=17), font_color="black")

    def create_polygon_path(center, radius, num_vertices, orientation=0):
        points = []
        for i in range(num_vertices):
            angle = orientation + i * (2 * np.pi / num_vertices)
            x = center[0] + radius * np.cos(angle)
            y = center[1] + radius * np.sin(angle)
            points.append(f"L {x} {y}")
        path = f"M {points[0][2:]} " + " ".join(points[1:]) + " Z"
        return path

    # Hexagon (centered at [0, 0], radius = 0.28, orientation = π/2)
    hexagon_path = create_polygon_path(center=[0, 0], radius=0.28, num_vertices=6, orientation=0)
    fig.add_shape(type="path", path=hexagon_path, line_color="black")

    # Pentagon (centered at [-0.375, -0.225], radius = 0.24, orientation = -0.42 radians)
    pentagon_path = create_polygon_path(center=[-0.375, -0.225], radius=0.24, num_vertices=5, orientation=-0.08)
    fig.add_shape(type="path", path=pentagon_path, line_color="black")

    # Customize layout
    fig.update_layout(
        xaxis=dict(range=[-1.1, 1.1], title="x (nm)", zeroline=False),
        yaxis=dict(range=[-1.1, 1.1], title="y (nm)", zeroline=False),
        height=1000,
        title="Visualisation of 2d base pairing",
        showlegend=False
    )

    return fig

def save_nth_frame(traj_file, top_file, index, session):
    # Load the trajectory and topology
    u = mda.Universe(top_file, traj_file)

    # Select the nth frame
    u.trajectory[index]

    # Define the path to save the frame
    path_save = os.path.join('static', session, 'temp_frame.pdb')

    # Save the frame
    with mda.Writer(path_save, multiframe=False) as W:
        W.write(u)

    return path_save

def create_contact_maps(traj_file, top_file, index, session):
    path_save_contact_maps = os.path.join('static', session, 'target.cmp')
    path_save_figure = os.path.join('static', session, 'contact_map_plotly.png')

    path_traj_file = os.path.join('static', session, traj_file)
    path_top_file = os.path.join('static', session, top_file)

    path_save_nth_frame = save_nth_frame(path_traj_file, path_top_file, index, session)
    contact_map, signature = generate_contact_maps.get_contact_map_and_signature(path_save_nth_frame)
    generate_contact_maps.write_contact_map(path_save_contact_maps, contact_map, signature)
    generate_contact_maps.read_contact_map(path_save_contact_maps, path_save_figure)
