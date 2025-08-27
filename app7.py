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
import pickle
import hashlib
import concurrent.futures
import multiprocessing
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, Any, Optional, Tuple
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

from task_celery import *
from task_celery import app as celery_app
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


class OptimizedTrajectoryManager:
    """Optimized trajectory loading and management with parallel processing"""
    
    def __init__(self):
        self._trajectory_cache = {}
        self._metadata_cache = {}
        self._lock = threading.Lock()
    
    def preload_trajectory_metadata(self, native_pdb_path: str, traj_xtc_path: str) -> Dict[str, Any]:
        """Preload trajectory metadata for optimization"""
        cache_key = f"{native_pdb_path}:{traj_xtc_path}"
        
        with self._lock:
            if cache_key in self._metadata_cache:
                return self._metadata_cache[cache_key]
        
        try:
            # Load minimal trajectory info without full data
            import mdtraj as md
            traj = md.load_frame(traj_xtc_path, 0, top=native_pdb_path)
            
            metadata = {
                'n_frames': md.load(traj_xtc_path, top=native_pdb_path).n_frames,
                'n_atoms': traj.n_atoms,
                'n_residues': traj.n_residues,
                'topology': traj.topology,
                'box_vectors': traj.unitcell_vectors[0] if traj.unitcell_vectors is not None else None
            }
            
            with self._lock:
                self._metadata_cache[cache_key] = metadata
            
            logger.info(f"Preloaded trajectory metadata: {metadata['n_frames']} frames, {metadata['n_atoms']} atoms")
            return metadata
            
        except Exception as e:
            logger.warning(f"Failed to preload trajectory metadata: {e}")
            return {}
    
    def load_trajectory_parallel(self, native_pdb_path: str, traj_xtc_path: str, stride: int = 1) -> Tuple[Any, Any]:
        """Load trajectory with parallel optimization"""
        cache_key = f"{native_pdb_path}:{traj_xtc_path}:{stride}"
        
        with self._lock:
            if cache_key in self._trajectory_cache:
                logger.info("Trajectory loaded from memory cache")
                return self._trajectory_cache[cache_key]
        
        def load_native():
            import mdtraj as md
            return md.load(native_pdb_path)
        
        def load_trajectory():
            import mdtraj as md
            if stride > 1:
                return md.load(traj_xtc_path, top=native_pdb_path, stride=stride)
            else:
                return md.load(traj_xtc_path, top=native_pdb_path)
        
        # Parallel loading
        with ThreadPoolExecutor(max_workers=2) as executor:
            native_future = executor.submit(load_native)
            traj_future = executor.submit(load_trajectory)
            
            native_load = native_future.result()
            traj_load = traj_future.result()
        
        # Cache for reuse
        with self._lock:
            self._trajectory_cache[cache_key] = (traj_load, native_load)
        
        logger.info(f"Loaded trajectory with {traj_load.n_frames} frames in parallel")
        return traj_load, native_load

trajectory_manager = OptimizedTrajectoryManager()

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
        'LANDSCAPE': CalculationRequirement('landscape', 300.0, 800, True, ['RMSD', 'ERMSD']),
        'BASE_PAIRING': CalculationRequirement('base_pairing', 150.0, 400, True, [])
    }

    def __init__(self, session_id: str):
        self.session_id = session_id
        self.logger = logging.getLogger(f'{__name__}.{session_id}')

    def get_shared_computations(self, selected_plots: List[str]) -> Dict[str, List[str]]:
        """Identify computations that can be shared between plots"""
        shared_computations = {}

        # Map computation types to plots that need them
        computation_to_plots = {}
        for plot in selected_plots:
            if plot in self.CALCULATION_SPECS:
                spec = self.CALCULATION_SPECS[plot]
                calc_type = spec.calculation_type

                if calc_type not in computation_to_plots:
                    computation_to_plots[calc_type] = []
                computation_to_plots[calc_type].append(plot)

                # Add dependencies
                for dep in spec.dependencies:
                    if dep in self.CALCULATION_SPECS:
                        dep_calc_type = self.CALCULATION_SPECS[dep].calculation_type
                        if dep_calc_type not in computation_to_plots:
                            computation_to_plots[dep_calc_type] = []
                        if plot not in computation_to_plots[dep_calc_type]:
                            computation_to_plots[dep_calc_type].append(plot)

        # Identify shared computations (needed by multiple plots)
        for calc_type, plots in computation_to_plots.items():
            if len(plots) > 1:
                shared_computations[calc_type] = plots

        return shared_computations

    def plan_calculations(self, selected_plots: List[str]) -> Dict[str, CalculationRequirement]:
        """Plan and log all calculations needed with optimization for shared computations"""
        self.logger.info(f"Planning calculations for plots: {selected_plots}")

        # Identify shared computations
        shared_computations = self.get_shared_computations(selected_plots)
        if shared_computations:
            self.logger.info(f"Shared computations identified: {shared_computations}")

        planned_calculations = {}
        total_time = 0
        total_memory = 0
        optimized_time = 0

        # Group calculations by type to avoid duplicates
        calculation_groups = {}
        unique_computations = set()

        for plot in selected_plots:
            if plot in self.CALCULATION_SPECS:
                spec = self.CALCULATION_SPECS[plot]
                calc_type = spec.calculation_type

                if calc_type not in calculation_groups:
                    calculation_groups[calc_type] = {
                        'plots': [plot],
                        'spec': spec
                    }
                    unique_computations.add(calc_type)
                    optimized_time += spec.estimated_time
                else:
                    calculation_groups[calc_type]['plots'].append(plot)

                planned_calculations[plot] = spec
                total_time += spec.estimated_time
                total_memory = max(total_memory, spec.memory_requirement)

        # Calculate time savings from optimization
        time_saved = total_time - optimized_time

        self.logger.info(f"Calculation plan:")
        for calc_type, group in calculation_groups.items():
            is_shared = calc_type in [ct for ct, plots in shared_computations.items()]
            shared_indicator = " [SHARED]" if is_shared else ""
            self.logger.info(f"  {calc_type}{shared_indicator}: {group['plots']} - {group['spec'].estimated_time}s, {group['spec'].memory_requirement}MB")

        self.logger.info(f"Total estimated time (unoptimized): {total_time}s ({total_time/60:.1f}min)")
        self.logger.info(f"Optimized time (shared computations): {optimized_time}s ({optimized_time/60:.1f}min)")
        self.logger.info(f"Time savings: {time_saved}s ({time_saved/60:.1f}min, {(time_saved/total_time)*100:.1f}%)")
        self.logger.info(f"Peak memory requirement: {total_memory}MB")

        return planned_calculations

    def optimize_calculation_order(self, planned_calculations: Dict[str, CalculationRequirement]) -> List[str]:
        """Optimize calculation order using dependency graph and parallel execution potential"""
        from collections import defaultdict, deque
        
        # Build dependency graph
        dependencies = defaultdict(set)
        dependents = defaultdict(set)
        
        for plot, spec in planned_calculations.items():
            for dep in spec.dependencies:
                if dep in planned_calculations:
                    dependencies[plot].add(dep)
                    dependents[dep].add(plot)
        
        # Topological sort with optimization for parallel execution
        in_degree = {plot: len(dependencies[plot]) for plot in planned_calculations}
        queue = deque([plot for plot in planned_calculations if in_degree[plot] == 0])
        ordered_plots = []
        
        # Priority: dependencies first, then by computation type optimization
        priority_weights = {
            'rmsd': 1,     # High priority - needed by landscape
            'ermsd': 2,
            'landscape': 10,  # Lower priority - depends on others
            'annotate': 5,
            'motif': 3,
            'jcoupling': 4,
            'escore': 4,
            'base_pairing': 6
        }
        
        while queue:
            # Sort queue by priority and parallel execution potential
            current_batch = list(queue)
            queue.clear()
            
            # Sort by priority weight and CPU intensity
            current_batch.sort(key=lambda plot: (
                priority_weights.get(planned_calculations[plot].calculation_type, 5),
                -planned_calculations[plot].estimated_time,  # Longer tasks first
                planned_calculations[plot].cpu_intensive
            ))
            
            for plot in current_batch:
                ordered_plots.append(plot)
                
                # Update dependencies
                for dependent in dependents[plot]:
                    in_degree[dependent] -= 1
                    if in_degree[dependent] == 0:
                        queue.append(dependent)
        
        # Group plots that can run in parallel (no dependencies between them)
        parallel_groups = []
        current_group = []
        processed_deps = set()
        
        for plot in ordered_plots:
            plot_deps = dependencies[plot]
            if plot_deps.issubset(processed_deps):
                current_group.append(plot)
            else:
                if current_group:
                    parallel_groups.append(current_group)
                current_group = [plot]
            processed_deps.add(plot)
        
        if current_group:
            parallel_groups.append(current_group)
        
        self.logger.info(f"Optimized calculation order with parallel groups: {parallel_groups}")
        self.logger.info(f"Final order: {ordered_plots}")
        
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

def parse_topology_file(file_path):
    """Parse topology file and extract residue information for torsion angle calculation"""
    try:
        u = mda.Universe(file_path)
        residues_info = []
        
        for residue in u.residues:
            res_info = {
                'index': int(residue.resindex),  # Convert numpy int64 to Python int
                'name': str(residue.resname),    # Ensure string type
                'id': int(residue.resid),        # Convert numpy int64 to Python int
                'full_name': f"{residue.resname}{residue.resid}"
            }
            residues_info.append(res_info)
        
        # Standard RNA torsion angles
        torsion_angles = [
            {'name': 'alpha', 'description': 'Alpha (O3\'-P-O5\'-C5\')'},
            {'name': 'beta', 'description': 'Beta (P-O5\'-C5\'-C4\')'},
            {'name': 'gamma', 'description': 'Gamma (O5\'-C5\'-C4\'-C3\')'},
            {'name': 'delta', 'description': 'Delta (C5\'-C4\'-C3\'-O3\')'},
            {'name': 'epsilon', 'description': 'Epsilon (C4\'-C3\'-O3\'-P)'},
            {'name': 'zeta', 'description': 'Zeta (C3\'-O3\'-P-O5\')'},
            {'name': 'chi', 'description': 'Chi (O4\'-C1\'-N1-C2/C4)'}
        ]
        
        return {
            'residues': residues_info,
            'torsion_angles': torsion_angles,
            'num_residues': len(residues_info),
            'sequence': [str(res.resname) for res in u.residues]  # Ensure strings
        }
    except Exception as e:
        logger.error(f"Error parsing topology file {file_path}: {str(e)}")
        return None


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

@app.route("/parse-topology", methods=["POST"])
def parse_topology():
    """Parse uploaded topology file and return residue information"""
    if "topologyFile" not in request.files:
        return jsonify({"error": "No topology file provided"}), 400
    
    topology_file = request.files["topologyFile"]
    if topology_file.filename == "":
        return jsonify({"error": "No file selected"}), 400
    
    session_id = request.form.get("session_id")
    if not session_id:
        return jsonify({"error": "No session ID provided"}), 400
    
    # Create session directory if it doesn't exist
    session_dir = os.path.join(app.static_folder, session_id)
    os.makedirs(session_dir, exist_ok=True)
    
    # Save the file temporarily for parsing
    temp_filename = secure_filename(topology_file.filename)
    temp_path = os.path.join(session_dir, f"temp_{temp_filename}")
    
    try:
        topology_file.save(temp_path)
        
        # Parse the file
        topology_info = parse_topology_file(temp_path)
        
        if topology_info is None:
            return jsonify({"error": "Failed to parse topology file"}), 500
        
        return jsonify({
            "success": True,
            "residues": topology_info["residues"],
            "torsion_angles": topology_info["torsion_angles"],
            "num_residues": topology_info["num_residues"],
            "sequence": topology_info["sequence"]
        })
        
    except Exception as e:
        logger.error(f"Error parsing topology file: {str(e)}")
        return jsonify({"error": f"Error parsing file: {str(e)}"}), 500
    finally:
        # Clean up temporary file
        if os.path.exists(temp_path):
            try:
                os.remove(temp_path)
            except:
                pass

@app.route("/cgu")
def cgu():
    return render_template("cgu.html")

@app.route("/authors")
def authors():
    return render_template("authors.html")

@app.route("/documentation")
def documentation():
    return render_template("documentation.html")

@app.route("/cache-debug")
def cache_debug():
    return render_template("cache-debug.html")

@app.route("/simple-test")
def simple_test():
    return render_template("simple-test.html")


@app.route("/benchmark-performance/<session_id>", methods=['POST'])
def benchmark_performance_route(session_id):
    """Run performance benchmark for different optimization strategies"""
    try:
        # Get file paths from session
        directory_path = os.path.join(app.static_folder, session_id)
        
        with open(os.path.join(directory_path, "session_data.json"), "r") as file:
            session_data = json.load(file)
        
        native_pdb = session_data['files']['nativePdb']
        traj_xtc = session_data['files']['trajXtc']
        
        native_pdb_path = os.path.join(directory_path, native_pdb)
        traj_xtc_path = os.path.join(directory_path, traj_xtc)
        
        # Start benchmark task
        benchmark_job = performance_benchmark.apply_async(
            args=[native_pdb_path, traj_xtc_path, session_id]
        )
        
        # Wait for results (with timeout)
        try:
            results = benchmark_job.get(timeout=300)  # 5 minute timeout
            return jsonify({
                'success': True,
                'benchmark_results': results,
                'session_id': session_id
            })
        except Exception as e:
            return jsonify({
                'success': False,
                'error': f'Benchmark timeout or failed: {str(e)}',
                'task_id': benchmark_job.id
            }), 500
        
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

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

    # Handle torsion angle parameters
    torsion_residues = request.form.getlist("torsionResidues")
    torsion_angles = request.form.getlist("torsionAngles")
    torsion_mode = request.form.get("torsionMode", "single")  # single, multiple, all
    
    session_data = {
        "selected_plots": selected_plots,
        "n_frames": n_frames,
        "frame_range": frame_range,
        "torsionResidue": request.form.get("torsionResidue", 0),
        "torsionResidues": torsion_residues,
        "torsionAngles": torsion_angles,
        "torsionMode": torsion_mode,
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

    # Optimized trajectory loading with parallel preprocessing
    socketio.emit('update_progress', {"progress": 25, "message": "Loading trajectory with optimizations..."}, to=session_id)
    
    def load_and_preprocess_trajectory():
        """Load trajectory with parallel optimization"""
        # Preload metadata first
        metadata = trajectory_manager.preload_trajectory_metadata(native_pdb_path, traj_xtc_path)
        
        # Load trajectory efficiently
        u = mda.Universe(native_pdb_path, traj_xtc_path)
        ref = mda.Universe(native_pdb_path)
        
        # Parallel computation of trajectory properties
        def compute_alignment():
            nucleic_backbone = u.select_atoms('nucleicbackbone').positions
            ref_backbone = ref.select_atoms('nucleicbackbone').positions
            return mda.analysis.align.rotation_matrix(nucleic_backbone, ref_backbone)[0]
        
        def extract_residue_info():
            return [residue.resname for residue in u.residues]
        
        # Execute in parallel
        with ThreadPoolExecutor(max_workers=2) as executor:
            r_matrix_future = executor.submit(compute_alignment)
            residue_future = executor.submit(extract_residue_info)
            
            r_matrix = r_matrix_future.result()
            residue_names = residue_future.result()
        
        return u, ref, r_matrix, residue_names
    
    # Load with optimization
    u, ref, r_matrix, residue_names = load_and_preprocess_trajectory()
    logger.info(f"Computed the rotation matrix: {r_matrix}")
    try:
        r_matrix_str = str(r_matrix.tolist())
    except:
        logger.warning("Got an issue with r_matrix computation")     
    logger.info(f"Loaded trajectory: {residue_names}")
    logger.info(f"Nucleic backbone atoms: {len(u.select_atoms('nucleicbackbone').positions)}")
    logger.info(f"Reference backbone atoms: {len(ref.select_atoms('nucleicbackbone').positions)}")
    
    # Optimized trajectory writing with parallel processing
    def write_trajectory_range():
        if frame_range != "all":
            start, end_stride = frame_range.split("-")
            end, stride = end_stride.split(":")
            start, end, stride = int(start), int(end), int(stride)
            u.trajectory[start:end:stride]
            output_path = os.path.join(directory_path, f"traj_{start}_{end}.xtc")
        else:
            output_path = os.path.join(directory_path, f"traj_uploaded.xtc")
            start, end, stride = 0, len(u.trajectory), 1
        
        # Parallel trajectory writing
        if frame_range != "all":
            with mda.Writer(output_path, n_atoms=u.atoms.n_atoms) as W:
                for ts in u.trajectory[start:end:stride]:
                    W.write(u)
        else:
            # If frame_range is "all", don't rewrite - just move/rename if needed
            if traj_xtc_path != output_path:
                #import shutil
                # trying symlink
                os.symlink(traj_xtc_path, output_path)
                #shutil.move(traj_xtc_path, output_path)
        
        return output_path
    
    # Use ThreadPoolExecutor for non-blocking trajectory writing
    with ThreadPoolExecutor(max_workers=1) as executor:
        traj_future = executor.submit(write_trajectory_range)
        traj_xtc_path = traj_future.result()

    socketio.emit('update_progress', {"progress": 40, "message": "Trajectory loaded and optimized."}, to=session_id)
    socketio.sleep(0.1)

    # Plan calculations with proper resource management
    planner = CalculationPlanner(session_id)
    planned_calculations = planner.plan_calculations(selected_plots)
    optimized_order = planner.optimize_calculation_order(planned_calculations)

    socketio.emit('update_progress', {"progress": 45, "message": "Calculations planned and optimized."}, to=session_id)
    socketio.sleep(0.1)

    # Phase 1: Compute metrics needed by plots
    metrics_needed = []
    if "RMSD" in selected_plots:
        metrics_needed.append("rmsd")
    if "ERMSD" in selected_plots:
        metrics_needed.append("ermsd")
    if "ANNOTATE" in selected_plots:
        metrics_needed.append("annotate")
    
    if metrics_needed:
        socketio.emit('update_progress', {"progress": 50, "message": "Computing metrics..."}, to=session_id)
        
        # Create Celery chain for metrics computation
        from celery import chain
        metrics_jobs = []
        
        if "rmsd" in metrics_needed:
            metrics_jobs.append(compute_rmsd.s(native_pdb_path, traj_xtc_path, session_id))
        if "ermsd" in metrics_needed:
            metrics_jobs.append(compute_ermsd.s(native_pdb_path, traj_xtc_path, session_id))
        if "annotate" in metrics_needed:
            metrics_jobs.append(compute_annotate.s(native_pdb_path, traj_xtc_path, session_id))
        
        # Execute metrics computation chain
        if metrics_jobs:
            metrics_chain = chain(*metrics_jobs)
            metrics_result = metrics_chain.apply_async()
            
            # Wait for metrics computation to complete
            while not metrics_result.ready():
                socketio.sleep(0.1)
            
            if metrics_result.successful():
                logger.info("All metrics computed successfully")
                socketio.emit('update_progress', {"progress": 55, "message": "Metrics computed."}, to=session_id)
            else:
                logger.error(f"Metrics computation failed: {metrics_result.result}")
                socketio.emit('update_progress', {"progress": 100, "message": "Error: Metrics computation failed"}, to=session_id)
                return

    # Phase 2: Generate plots using computed metrics  
    socketio.emit('update_progress', {"progress": 60, "message": "Generating plots..."}, to=session_id)
    plot_data = []
    active_jobs = {}  # Track running jobs
    max_parallel_jobs = min(multiprocessing.cpu_count(), 4)  # Limit concurrent jobs
    
    def submit_plot_job(plot):
        """Submit a single plot job and return job info"""
        files_path = os.path.join(download_path, plot)
        plot_dir = os.path.join(download_plot, plot)
        os.makedirs(files_path, exist_ok=True)
        os.makedirs(plot_dir, exist_ok=True)

        logger.info(f"Starting calculation for {plot} - estimated time: {planned_calculations[plot].estimated_time}s")

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
            return None

        if job:
            logger.info(f"Enqueued {plot} calculation with job ID: {job.id}")
            return [plot, plot_style, job.id]
        return None
    
    # Parallel job submission with dependency management
    completed_dependencies = set()
    pending_plots = list(optimized_order)
    
    while pending_plots or active_jobs:
        # Submit new jobs for plots whose dependencies are met
        submittable_plots = []
        for plot in pending_plots[:]:
            plot_spec = planned_calculations[plot]
            dependencies_met = all(dep in completed_dependencies for dep in plot_spec.dependencies)
            
            if dependencies_met and len(active_jobs) < max_parallel_jobs:
                submittable_plots.append(plot)
                pending_plots.remove(plot)
        
        # Submit jobs in parallel
        for plot in submittable_plots:
            try:
                job_info = submit_plot_job(plot)
                if job_info:
                    plot_data.append(job_info)
                    active_jobs[plot] = job_info[2]  # Store the job ID
                    logger.info(f"Submitted {plot} for parallel execution")
                    
                    socketio.emit('update_progress', {"progress": 60, "message": f"{plot} plot enqueued for parallel execution."}, to=session_id)
                    socketio.sleep(0.1)
                    
            except Exception as e:
                logger.error(f"Failed to enqueue {plot} calculation: {str(e)}")
                socketio.emit('update_progress', {"progress": 60, "message": f"Error enqueueing {plot}: {str(e)}"}, to=session_id)
        
        # Check for completed jobs
        completed_jobs = []
        for plot, job_id in active_jobs.items():
            job = celery_app.AsyncResult(job_id)
            if job.ready():
                completed_jobs.append(plot)
                completed_dependencies.add(plot)
                logger.info(f"Completed {plot} calculation in parallel")
        
        # Remove completed jobs
        for plot in completed_jobs:
            del active_jobs[plot]
        
        # Brief sleep to prevent busy waiting
        if active_jobs:
            time.sleep(0.1)

    socketio.emit('update_progress', {"progress": 70, "message": "Waiting for plot jobs to complete..."}, to=session_id)
    socketio.sleep(0.1)  # Ensure the message is sent

    # Wait for and process results with better error handling
    completed_plot_data = []
    total_jobs = len(plot_data)
    completed_jobs = 0

    for plot_type, plot_style, job_id in plot_data:
        try:
            job = celery_app.AsyncResult(job_id)
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
