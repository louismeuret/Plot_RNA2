import os
import sys

# Set environment variables to prevent conflicts before importing eventlet
os.environ['EVENTLET_NO_GREENDNS'] = '1'

# Try eventlet import with better error handling
try:
    import eventlet
    # Conservative monkey patching to avoid conflicts
    eventlet.monkey_patch(socket=True, dns=False, time=False, select=True, thread=False, os=False)
    EVENTLET_AVAILABLE = True
except ImportError as e:
    print(f"Warning: eventlet import failed: {e}")
    EVENTLET_AVAILABLE = False
except Exception as e:
    print(f"Warning: eventlet configuration failed: {e}")
    EVENTLET_AVAILABLE = False

import time
from flask import (
    Flask,
    request,
    render_template,
    send_from_directory,
    redirect,
    url_for,
    session,
    send_file,
    jsonify,
)
from werkzeug.utils import secure_filename
import uuid
import os
import glob
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.coordinates import PDB
from MDAnalysis.analysis import rms
from MDAnalysis.topology.guessers import guess_types
from MDAnalysis.analysis.align import alignto
from rq import Queue
from redis import Redis
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import numpy as np
import json
import plotly
from plotly.tools import mpl_to_plotly

from flask_socketio import SocketIO, join_room, leave_room, emit
from io import BytesIO
from PIL import Image, ImageDraw
import sys
import numpy as np

from FoldingAnalysis.analysis import *
import plotly.graph_objects as go
from string import ascii_uppercase
import barnaba as bb
import pickle
import io
import zipfile

from tasks_celery import *
from celery.result import AsyncResult
import threading
import logging
from dataclasses import dataclass
from typing import Dict, List, Optional
from datetime import datetime

app = Flask(__name__)
app.secret_key = "pi"
app.debug = True

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

@dataclass
class CalculationRequirement:
    """Represents a calculation requirement for a plot"""
    calculation_type: str
    estimated_time: float  # in seconds
    memory_requirement: float  # in MB
    cpu_intensive: bool
    dependencies: List[str]
    
class CalculationPlanner:
    """Plans and manages calculation resources"""
    
    CALCULATION_SPECS = {
        'RMSD': CalculationRequirement('rmsd', 30.0, 100, True, []),
        'ERMSD': CalculationRequirement('ermsd', 45.0, 150, True, []),
        'TORSION': CalculationRequirement('torsion', 60.0, 200, True, []),
        'SEC_STRUCTURE': CalculationRequirement('annotate', 120.0, 300, True, []),
        'DOTBRACKET': CalculationRequirement('annotate', 120.0, 300, True, []),
        'ARC': CalculationRequirement('annotate', 130.0, 350, True, []),
        'CONTACT_MAPS': CalculationRequirement('annotate', 180.0, 500, True, []),
        'ANNOTATE': CalculationRequirement('annotate', 120.0, 300, True, []),
        'DS_MOTIF': CalculationRequirement('motif', 90.0, 250, False, []),
        'SS_MOTIF': CalculationRequirement('motif', 90.0, 250, False, []),
        'JCOUPLING': CalculationRequirement('jcoupling', 75.0, 200, True, []),
        'ESCORE': CalculationRequirement('escore', 60.0, 150, False, []),
        'LANDSCAPE': CalculationRequirement('landscape', 300.0, 800, True, ['RMSD']),
        'BASE_PAIRING': CalculationRequirement('base_pairing', 150.0, 400, True, [])
    }
    
    def __init__(self, session_id: str):
        self.session_id = session_id
        self.logger = logging.getLogger(f'{__name__}.{session_id}')
        
    def plan_calculations(self, selected_plots: List[str]) -> Dict[str, CalculationRequirement]:
        """Plan and log all calculations needed"""
        self.logger.info(f"Planning calculations for plots: {selected_plots}")
        
        planned_calculations = {}
        total_time = 0
        total_memory = 0
        
        # Group calculations by type to avoid duplicates
        calculation_groups = {}
        
        for plot in selected_plots:
            if plot in self.CALCULATION_SPECS:
                spec = self.CALCULATION_SPECS[plot]
                calc_type = spec.calculation_type
                
                if calc_type not in calculation_groups:
                    calculation_groups[calc_type] = {
                        'plots': [plot],
                        'spec': spec
                    }
                else:
                    calculation_groups[calc_type]['plots'].append(plot)
                    
                planned_calculations[plot] = spec
                total_time += spec.estimated_time
                total_memory = max(total_memory, spec.memory_requirement)
        
        self.logger.info(f"Calculation plan:")
        for calc_type, group in calculation_groups.items():
            self.logger.info(f"  {calc_type}: {group['plots']} - {group['spec'].estimated_time}s, {group['spec'].memory_requirement}MB")
            
        self.logger.info(f"Total estimated time: {total_time}s ({total_time/60:.1f}min)")
        self.logger.info(f"Peak memory requirement: {total_memory}MB")
        
        return planned_calculations
        
    def optimize_calculation_order(self, planned_calculations: Dict[str, CalculationRequirement]) -> List[str]:
        """Optimize the order of calculations for efficiency"""
        # Group by calculation type to share results
        calc_groups = {}
        for plot, spec in planned_calculations.items():
            calc_type = spec.calculation_type
            if calc_type not in calc_groups:
                calc_groups[calc_type] = []
            calc_groups[calc_type].append(plot)
            
        # Order: CPU-intensive first, then memory-intensive
        ordered_plots = []
        for calc_type, plots in calc_groups.items():
            ordered_plots.extend(plots)
            
        self.logger.info(f"Optimized calculation order: {ordered_plots}")
        return ordered_plots
socketio = SocketIO(
    app,
    logger=True,
    engineio_logger=True,
    async_mode='eventlet' if EVENTLET_AVAILABLE else 'threading',
    cors_allowed_origins=[
        "https://arny-plotter.rpbs.univ-paris-diderot.fr",
        "http://arny-plotter.rpbs.univ-paris-diderot.fr",
        "http://localhost:4242",
        "http://127.0.0.1:4242",
        "http://172.27.7.130:4242"
    ],

    cors_credentials=True,     # Allow credentials in CORS requests
    allow_upgrades=True,       # Allow protocol upgrades
    transports=['websocket', 'polling'],  # Support both transport methods

    # Additional proxy configurations
    ping_timeout=60,           # Increase timeout for proxy delays
    ping_interval=25,          # Regular ping to keep connection alive

    # Handle proxy headers
    engineio_options={
        'ping_timeout': 60,
        'ping_interval': 25,
        'upgrade_timeout': 30,
        'max_http_buffer_size': 1000000,
        # Additional engineio CORS options
        'cors_allowed_origins': [
            "https://arny-plotter.rpbs.univ-paris-diderot.fr",
            "http://arny-plotter.rpbs.univ-paris-diderot.fr",
            "http://localhost:4242",
            "http://127.0.0.1:4242",
            "http://172.27.7.130:4242"
        ],
        'cors_credentials': True,
    }
)
redis_conn = Redis()
plot_queue = Queue('plot_queue', connection=redis_conn)

@socketio.on("connect")
def handle_connect():
    session_id = request.args.get("session_id")
    if session_id:
        join_room(session_id)
        socketio.emit(
            "command_output",
            {"output": "Connected to room: " + session_id},
            to=session_id,
        )
    print("Client connected")

def create_session_id():
    return uuid.uuid4().hex


@app.route("/", methods=["GET", "POST"])
def index():
    if request.method == "POST":
        directory = request.form["directory"]
        session_id = request.form["session_id"]
        # Perform actions based on the form input, e.g., list files
        files = list_files(directory)
        return render_template("results.html", files=files, session_id=session_id)
    session_id = uuid.uuid4().hex
    session["session_id"] = session_id
    return render_template("index.html", session_id=session_id)

@app.route("/get_session", methods=["GET"])
def get_session():
    session["session_id"] = uuid.uuid4().hex

    session_id = session["session_id"]
    session_path = os.path.join(app.static_folder, session_id)

    # Check if the session_id exists as a directory in /static
    session_exists = os.path.isdir(session_path)

    return jsonify({"session_id": session_id, "exists": session_exists})

@app.route("/cgu")
def cgu():
    return render_template("cgu.html")

@app.route("/authors")
def authors():
    return render_template("authors.html")

@app.route("/documentation")
def documentation():
    return render_template("documentation.html")

@app.route("/documentationll")
def documentation2():
    return render_template("documentationll.html")

def list_files(directory):
    # Your implementation of listing files
    return os.listdir(directory)  # Simplified example




def plotly_to_json(fig):
    return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)

@app.route("/upload-files", methods=["POST"])
def upload_files():
    print(request.files)
    print(request.form)
    if "nativePdb" not in request.files or "trajXtc" not in request.files:
        return redirect(request.url)

    native_pdb = request.files["nativePdb"]
    traj_xtc = request.files["trajXtc"]

    if native_pdb.filename == "" or traj_xtc.filename == "":
        return redirect(request.url)

    # Validate file extensions
    """
    if not (native_pdb.filename.endswith('.pdb') and traj_xtc.filename.endswith('.xtc')):
        return "Invalid file type. Please upload a PDB and XTC file.", 400
    """
    print(request.form)
    session_id = request.form.get("session_id")
    if not session_id:
        session_id = session.get("session_id")

    print(f"Session_id = {session_id}")
    os.makedirs(os.path.join(app.static_folder, session_id), exist_ok=True)

    native_pdb_path = os.path.join(app.static_folder, session_id, secure_filename(native_pdb.filename))
    traj_xtc_path = os.path.join(app.static_folder, session_id, secure_filename(traj_xtc.filename))
    print(traj_xtc_path)

    try:
        native_pdb.save(native_pdb_path)
        traj_xtc.save(traj_xtc_path)
    except Exception as e:
        app.logger.error(f"Error saving files: {e}")
        return "Error saving files.", 500

    selected_plots = []
    n_frames = int(request.form.get("n_frames", 1))  # Default to 1 if not specified
    first_frame = request.form.get("firstFrame", "")
    last_frame = request.form.get("lastFrame", "")
    stride = request.form.get("stride", "")
    frame_range = "all" if not (first_frame and last_frame and stride) else f"{first_frame}-{last_frame}:{stride}"

    for plot in [
        "RMSD",
        "ERMSD",
        "CONTACT_MAPS",
        "TORSION",
        "SEC_STRUCTURE",
        "DOTBRACKET",
        "ARC",
        "ANNOTATE",
        "DS_MOTIF",
        "SS_MOTIF",
        "JCOUPLING",
        "ESCORE",
        "LANDSCAPE",
        "BASE_PAIRING",
    ]:  # Adjust based on your available plots
        app.logger.info(plot)
        print(str(plot.lower()))
        if plot.lower() in request.form:
            app.logger.info(request.form)
            app.logger.info("In the Request part of the code")
            print(plot)
            selected_plots.append(plot)

    session_data = {
        "selected_plots": selected_plots,
        "n_frames": n_frames,
        "frame_range": frame_range,
        "torsionResidue": request.form.get("torsionResidue", 0),
        "landscape_stride": request.form.get("landscape_stride", 0),
        "landscape_first_component": request.form.get("firstDimension", 0),
        "landscape_second_component": request.form.get("secondDimension", 0),
        "form": list(request.form),
        "files": {
            "nativePdb": str(request.files["nativePdb"].filename),
            "trajXtc": str(request.files["trajXtc"].filename),
            "trajpdb": str(request.files["trajXtc"].filename)
        }
    }

    session.update(session_data)

    session_file_path = os.path.join(app.static_folder, session_id, "session_data.json")
    with open(session_file_path, "w") as session_file:
        json.dump(session_data, session_file)

    app.logger.info(selected_plots)

    return redirect(
        url_for(
            "view_trajectory",
            session_id=session_id,
            native_pdb=native_pdb.filename,
            traj_xtc=traj_xtc.filename,
        )
    )


@app.route("/retrieve-results", methods=["GET"])
def retrieve_results():
    print("RETRIVE RESULTS")
    session_id = request.args.get("session_id")
    print(session_id)

    directory_path = os.path.join(app.static_folder, session_id)
    pickle_file_path = os.path.join(directory_path, "plot_data.pkl")

    with open(os.path.join(directory_path, "session_data.json"), "r") as file:
        session_data = json.load(file)
        print(session_data)

    native_pdb = session_data['files']['nativePdb']
    traj_xtc = session_data['files']['trajXtc']


    native_pdb_path = os.path.join(directory_path, native_pdb)
    traj_xtc_path = os.path.join(directory_path, traj_xtc)
    print(f"Native PDB: {native_pdb_path}")
    print(f"Trajectory XTC: {traj_xtc_path}")

    u = mda.Universe(native_pdb_path, traj_xtc_path)

    if not os.path.exists(pickle_file_path):
        return "Session results not found", 404

    with open(pickle_file_path, "rb") as f:
        plot_data = pickle.load(f)

    with open(os.path.join(app.static_folder,'explanations.json')) as f:
        explanations = json.load(f)

    print(explanations)

    return render_template(
        "view_trajectory.html",
        session_id=session_id,
        native_pdb=native_pdb,
        traj_xtc=traj_xtc,
        plot_data=plot_data,
        trajectory_length=len(u.trajectory),
        explainations=explanations
    )


@app.route("/download/plot_data/<session_id>/<plot_id>")
def download_plot_data(plot_id, session_id):
    directory_path = os.path.join(app.static_folder, session_id)
    download_path = os.path.join(directory_path, "download_data")
    files_path = os.path.join(download_path, plot_id)
    print(f"Download path: {files_path}")
    file_paths = glob.glob(f"{files_path}/*")
    print(file_paths)

    if len(file_paths) > 1:
        memory_file = io.BytesIO()
        download_filename = f"{plot_id}_data.zip"
        with zipfile.ZipFile(memory_file, "w", zipfile.ZIP_DEFLATED) as zf:
            for file_path in file_paths:
                with open(file_path, "rb") as f:
                    # Use the file name as the arcname
                    zf.writestr(file_path.split("/")[-1], f.read())
    else:
        memory_file = io.BytesIO()
        with open(file_paths[0], "rb") as f:
            memory_file.write(f.read())
        download_filename = os.path.basename(file_paths[0])

    memory_file.seek(0)

    # Send the zip file as a response
    return send_file(memory_file, download_name=download_filename, as_attachment=True)


@app.route("/download/plot/<session_id>/<plot_id>")
def download_plot(plot_id, session_id):
    directory_path = os.path.join(app.static_folder, session_id)
    download_path = os.path.join(directory_path, "download_plot")
    files_path = os.path.join(download_path, plot_id)
    print(f"Download path: {files_path}")
    file_paths = glob.glob(f"{files_path}/*")
    print(file_paths)

    if len(file_paths) > 1:
        memory_file = io.BytesIO()
        download_filename = f"{plot_id}_plots.zip"
        with zipfile.ZipFile(memory_file, "w", zipfile.ZIP_DEFLATED) as zf:
            for file_path in file_paths:
                with open(file_path, "rb") as f:
                    # Use the file name as the arcname
                    zf.writestr(file_path.split("/")[-1], f.read())
    else:
        memory_file = io.BytesIO()
        with open(file_paths[0], "rb") as f:
            memory_file.write(f.read())
        download_filename = os.path.basename(file_paths[0])

    memory_file.seek(0)

    # Send the zip file as a response
    return send_file(memory_file, download_name=download_filename, as_attachment=True)

@app.route('/view-trajectory/<session_id>/<native_pdb>/<traj_xtc>')
def view_trajectory(session_id, native_pdb, traj_xtc):
    start_time = time.time()
    socketio.emit('update_progress', {"progress": 0, "message": "Initializing trajectory analysis..."}, to=session_id)
    socketio.sleep(0.1)  # Ensure the message is sent

    # Setup paths and load explanations
    directory_path = os.path.join(app.static_folder, session_id)
    download_path = os.path.join(directory_path, "download_data")
    download_plot = os.path.join(directory_path, "download_plot")
    generate_data_path = os.path.join(directory_path, "generated_data")

    # Create directories if they don't exist
    os.makedirs(download_path, exist_ok=True)
    os.makedirs(download_plot, exist_ok=True)
    os.makedirs(generate_data_path, exist_ok=True)

    with open(os.path.join(app.static_folder, 'explanations.json')) as f:
        explanations = json.load(f)

    socketio.emit('update_progress', {"progress": 10, "message": "Paths set up and explanations loaded."}, to=session_id)
    socketio.sleep(0.1)  # Ensure the message is sent

    # Validate paths and session data
    native_pdb_path = os.path.join(directory_path, native_pdb)
    traj_xtc_path = os.path.join(directory_path, traj_xtc)
    if not os.path.exists(native_pdb_path) or not os.path.exists(traj_xtc_path):
        socketio.emit('update_progress', {"progress": 100, "message": "Error: File not found."}, to=session_id)
        return

    selected_plots = session.get("selected_plots", [])
    n_frames = session.get("n_frames", 1)
    frame_range = session.get("frame_range", "all")

    socketio.emit('update_progress', {"progress": 20, "message": "Session data validated."}, to=session_id)
    socketio.sleep(0.1)  # Ensure the message is sent

    # Load the trajectory
    u = mda.Universe(native_pdb_path, traj_xtc_path)
    ref = mda.Universe(native_pdb_path)
    residue_names = [residue.resname for residue in u.residues]
    print(residue_names)
    print(len(u.select_atoms('nucleicbackbone').positions))
    print(len(ref.select_atoms('nucleicbackbone').positions))
    r_matrix = mda.analysis.align.rotation_matrix(u.select_atoms('nucleicbackbone').positions, ref.select_atoms('nucleicbackbone').positions)[0]
    r_matrix_str = str(r_matrix.tolist())
    print(r_matrix_str)
    print(r_matrix)
    if frame_range != "all":
        start, end_stride = frame_range.split("-")
        end, stride = end_stride.split(":")
        start, end, stride = int(start), int(end), int(stride)
        u.trajectory[start:end:stride]
        traj_xtc_path = os.path.join(directory_path, f"traj_{start}_{end}.xtc")
        with mda.Writer(traj_xtc_path, n_atoms=u.atoms.n_atoms) as W:
            for ts in u.trajectory[start:end:stride]:
                W.write(u)
    else:
        traj_xtc_path = os.path.join(directory_path, f"traj_uploaded.xtc")
        with mda.Writer(traj_xtc_path, n_atoms=u.atoms.n_atoms) as W:
            for ts in u.trajectory[::]:
                W.write(u)


    socketio.emit('update_progress', {"progress": 40, "message": "Trajectory loaded."}, to=session_id)
    socketio.sleep(0.1)
    
    # Plan calculations with proper resource management
    planner = CalculationPlanner(session_id)
    planned_calculations = planner.plan_calculations(selected_plots)
    optimized_order = planner.optimize_calculation_order(planned_calculations)
    
    socketio.emit('update_progress', {"progress": 45, "message": "Calculations planned and optimized."}, to=session_id)
    socketio.sleep(0.1)
    
    plot_data = []
    for plot in optimized_order:
        files_path = os.path.join(download_path, plot)
        plot_dir = os.path.join(download_plot, plot)
        os.makedirs(files_path, exist_ok=True)
        os.makedirs(plot_dir, exist_ok=True)
        
        logger.info(f"Starting calculation for {plot} - estimated time: {planned_calculations[plot].estimated_time}s")

        try:
            job = None
            plot_style = "default"
            
            if plot == "RMSD":
                job = generate_rmsd_plot.apply_async(args=[native_pdb_path, traj_xtc_path, files_path, plot_dir, session_id])
                plot_style = "scatter2"
            elif plot == "ERMSD":
                job = generate_ermsd_plot.apply_async(args=[native_pdb_path, traj_xtc_path, files_path, plot_dir, session_id])
                plot_style = "scatter"
            elif plot == "TORSION":
                torsion_residue = session.get("torsionResidue", 0)
                job = generate_torsion_plot.apply_async(args=[native_pdb_path, traj_xtc_path, files_path, plot_dir, session_id, torsion_residue])
                plot_style = "torsion"
            elif plot == "SEC_STRUCTURE":
                job = generate_sec_structure_plot.apply_async(args=[native_pdb_path, traj_xtc_path, files_path, plot_dir, session_id])
                plot_style = "bar"
            elif plot == "DOTBRACKET":
                job = generate_dotbracket_plot.apply_async(args=[native_pdb_path, traj_xtc_path, files_path, plot_dir, session_id])
                plot_style = "dotbracket"
            elif plot == "ARC":
                job = generate_arc_plot.apply_async(args=[native_pdb_path, traj_xtc_path, files_path, plot_dir, session_id])
                plot_style = "arc"
            elif plot == "CONTACT_MAPS":
                job = generate_contact_map_plot.apply_async(args=[native_pdb_path, traj_xtc_path, files_path, plot_dir, generate_data_path, session_id])
                plot_style = "CONTACT_MAPS"
            elif plot == "ANNOTATE":
                job = generate_annotate_plot.apply_async(args=[native_pdb_path, traj_xtc_path, files_path, plot_dir, session_id])
                plot_style = "annotate"
            elif plot == "DS_MOTIF":
                job = generate_ds_motif_plot.apply_async(args=[native_pdb_path, traj_xtc_path, files_path, plot_dir, session_id])
                plot_style = "motif"
            elif plot == "SS_MOTIF":
                job = generate_ss_motif_plot.apply_async(args=[native_pdb_path, traj_xtc_path, files_path, plot_dir, session_id])
                plot_style = "motif"
            elif plot == "JCOUPLING":
                job = generate_jcoupling_plot.apply_async(args=[native_pdb_path, traj_xtc_path, files_path, plot_dir, session_id])
                plot_style = "scatter"
            elif plot == "ESCORE":
                job = generate_escore_plot.apply_async(args=[native_pdb_path, traj_xtc_path, files_path, plot_dir, session_id])
                plot_style = "scatter"
            elif plot == "LANDSCAPE":
                landscape_params = [
                    session.get("landscape_stride", 1),
                    session.get("landscape_first_component", 1),
                    session.get("landscape_second_component", 1)
                ]
                job = generate_landscape_plot.apply_async(args=[native_pdb_path, traj_xtc_path, files_path, plot_dir, session_id, landscape_params, generate_data_path])
                plot_style = "surface"
            elif plot == "BASE_PAIRING":
                job = generate_2Dpairing_plot.apply_async(args=[native_pdb_path, traj_xtc_path, files_path, plot_dir, session_id])
                plot_style = "2Dpairing"
            else:
                logger.warning(f"Unknown plot type: {plot}")
                continue
                
            if job:
                plot_data.append([plot, plot_style, job.id])
                logger.info(f"Enqueued {plot} calculation with job ID: {job.id}")
                
        except Exception as e:
            logger.error(f"Failed to enqueue {plot} calculation: {str(e)}")
            socketio.emit('update_progress', {"progress": 60, "message": f"Error enqueueing {plot}: {str(e)}"}, to=session_id)

        socketio.emit('update_progress', {"progress": 60, "message": f"{plot} plot enqueued."}, to=session_id)
        socketio.sleep(0.1)  # Ensure the message is sent

    socketio.emit('update_progress', {"progress": 70, "message": "Waiting for plot jobs to complete..."}, to=session_id)
    socketio.sleep(0.1)  # Ensure the message is sent

    # Wait for and process results with better error handling
    completed_plot_data = []
    total_jobs = len(plot_data)
    completed_jobs = 0
    
    for plot_type, plot_style, job_id in plot_data:
        try:
            job = app2.AsyncResult(job_id)
            start_time = time.time()
            
            # Wait with timeout and progress updates
            while not job.ready():
                elapsed = time.time() - start_time
                if elapsed > 600:  # 10 minute timeout
                    logger.error(f"Timeout waiting for {plot_type} calculation")
                    socketio.emit('update_progress', {"progress": 80, "message": f"Timeout: {plot_type} calculation taking too long"}, to=session_id)
                    break
                time.sleep(0.1)

            if job.successful():
                completed_plot_data.append([plot_type, plot_style, job.result])
                completed_jobs += 1
                progress = 70 + (completed_jobs / total_jobs) * 20
                logger.info(f"Completed {plot_type} calculation successfully")
                socketio.emit('update_progress', {"progress": progress, "message": f"{plot_type} plot completed ({completed_jobs}/{total_jobs})."}, to=session_id)
            else:
                error_msg = str(job.result) if job.result else "Unknown error"
                logger.error(f"Error in {plot_type} calculation: {error_msg}")
                socketio.emit('update_progress', {"progress": 80, "message": f"Error with {plot_type} plot: {error_msg}"}, to=session_id)
                
        except Exception as e:
            logger.error(f"Exception processing {plot_type} result: {str(e)}")
            socketio.emit('update_progress', {"progress": 80, "message": f"Exception with {plot_type}: {str(e)}"}, to=session_id)

    pickle_file_path = os.path.join(directory_path, "plot_data.pkl")
    print(f"Dir Path = {directory_path}")
    print(f"Completed Plot Data: {completed_plot_data}")  # Debugging line to check the data
    with open(pickle_file_path, "wb") as f:
        pickle.dump(completed_plot_data, f, protocol=pickle.HIGHEST_PROTOCOL)

    socketio.emit('update_progress', {"progress": 100, "message": "Trajectory analysis complete."}, to=session_id)

    return render_template(
            "view_trajectory.html",
            session_id=session_id,
            native_pdb=native_pdb,
            traj_xtc="traj_uploaded.xtc",
            plot_data=completed_plot_data,
            trajectory_length=len(u.trajectory),
            explainations=explanations,
            rotation_matrix=r_matrix_str,
        )



@app.route("/plot-biased", methods=["POST"])
def plot_biased():
    directory = request.form["directory"]
    df = read_biased(os.path.join(directory, "ratchet.out"))
    plot_filenames = save_all_plots(df, directory)
    return render_template("biased_plots.html", plot_filenames=plot_filenames)


def read_biased(file_path):
    df = pd.read_csv(file_path, sep="\t", header=None)
    return df


@socketio.on('update_frame_displayed')
def seek_right_frame_landscape(data, second_arg):
    #As for input coordinates on the energy landscape, and return the correspondin closest frame
    print(f"Data received from update frame displayed: {data}")
    session = data['session_id']
    directory_path = os.path.join(app.static_folder, session)
    generate_data_path = os.path.join(directory_path, "generated_data")
    coordinates = {
        'Q': float(data['value'][0]),  # Convert to float to ensure JSON serialization
        'RMSD': float(data['value'][1])  # Convert to float to ensure JSON serialization
    }
    print(f"Coordinates received: Q={coordinates['Q']}, RMSD={coordinates['RMSD']}")
    job = update_landscape_frame.apply_async(args=[generate_data_path, coordinates])
    result = job.get()
    print(result)
    emit('update_frame_landscape_click', {'frame': int(result)}, room=session)  # Convert result to int for JSON serialization



@socketio.on('slider_value')
def handle_slider_value(data):
    slider_value = int(data['value'])
    traj = data['traj']
    native_state = data['native_state']
    session = data['session_id']
    print(f"slider value = {slider_value}")
    path_contactmap_figure = create_contact_maps(traj, native_state, slider_value, session)


    path_save_2 = 'contact_map_plotly.png'

    emit('image_update', {'image_url': path_save_2})

@socketio.on('update_contact_map')
def update_contact_map_to_right_frame(data, second_arg):
    #print(second_arg)
    print('Trying to update contact map')
    session = data['session_id']
    slider_value = int(data['value'])
    directory_path = os.path.join(app.static_folder, session)
    generate_data_path = os.path.join(directory_path, "generated_data")
    download_plot = os.path.join(directory_path, "download_plot", "CONTACT_MAPS")
    print(f"SESSION = {session}")
    job = update_contact_map_plot.apply_async(args=[generate_data_path, download_plot, slider_value, session])
    result = job.get()
    emit('contact_map_plot_update', {'plotData': result}, room=session)


def generate_image(value):
    # Create a simple image based on the slider value
    print("generating image")
    width, height = 200, 100
    image = Image.new('RGB', (width, height), color=(255, 255, 255))
    draw = ImageDraw.Draw(image)
    draw.rectangle([10, 10, 10 + value, 90], fill=(255, 0, 0))
    return image


def save_plot(df, time_column, column_to_plot, output_folder):
    plt.figure(figsize=(10, 6))
    plt.plot(df[time_column], df[column_to_plot])
    plt.title(f"{column_to_plot} over time")
    plt.xlabel("Time")
    plt.ylabel(column_to_plot)
    plot_filename = f"{column_to_plot}_over_time.png"
    plt.savefig(os.path.join(output_folder, plot_filename))
    plt.close()
    return plot_filename


def save_all_plots(df, directory):
    output_folder = os.path.join(directory, "plots")
    os.makedirs(output_folder, exist_ok=True)
    time_column = df.columns[0]
    plot_filenames = []
    for column_to_plot in df.columns[1:]:
        plot_filename = save_plot(df, time_column, column_to_plot, output_folder)
        plot_filenames.append(plot_filename)
    return plot_filenames

if __name__ == "__main__":
    # eventlet.monkey_patch() already done at top
    with app.app_context():
        # app.run(debug=True)
        socketio.run(app)
