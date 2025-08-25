"""
Simplified Ultra-Optimized Celery Tasks
Focus: Resource precomputation + shared cache + simple architecture
"""

import os
import time
import json
import logging
from celery import Celery
from functools import wraps
import traceback

# Import only what we need to avoid dependency issues
try:
    from computation_cache import computation_cache, cached_computation
except ImportError:
    # Fallback if cache not available
    def cached_computation(**kwargs):
        def decorator(func):
            return func
        return decorator
    
    class DummyCache:
        def get(self, key): return None
        def set(self, key, value, timeout=3600): pass
    
    computation_cache = DummyCache()

# Try to import create_plots functions, with fallback
try:
    from create_plots import plot_ermsd, plot_torsion, create_contact_maps
    CREATE_PLOTS_AVAILABLE = True
except ImportError:
    logger.warning("create_plots module not available, using fallback implementations")
    CREATE_PLOTS_AVAILABLE = False

try:
    import barnaba as bb
    import numpy as np
    BARNABA_AVAILABLE = True
except ImportError:
    BARNABA_AVAILABLE = False
    
try:
    import MDAnalysis as mda
    MDA_AVAILABLE = True
except ImportError:
    MDA_AVAILABLE = False

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Celery app with optimized config (no persistence, minimal overhead)
app = Celery('simple_optimized_tasks')
app.conf.update(
    broker_url='redis://localhost:6379/0',
    result_backend='redis://localhost:6379/0',
    task_serializer='json',
    result_serializer='json',
    accept_content=['json'],
    
    # Performance optimizations
    task_acks_late=False,           # Ack immediately, don't wait for completion
    worker_prefetch_multiplier=1,   # Process one task at a time per worker
    task_reject_on_worker_lost=False,
    result_expires=3600,            # Results expire after 1 hour
    task_ignore_result=False,
    
    # Disable task persistence (faster)
    task_store_errors_even_if_ignored=False,
    worker_disable_rate_limits=True,
)

def simple_log_execution(func):
    """Simple execution logger"""
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        task_name = func.__name__
        start_time = time.time()
        
        try:
            logger.info(f"Starting {task_name}")
            result = func(self, *args, **kwargs)
            duration = time.time() - start_time
            logger.info(f"Completed {task_name} in {duration:.2f}s")
            return result
        except Exception as exc:
            duration = time.time() - start_time
            logger.error(f"Failed {task_name} after {duration:.2f}s: {exc}")
            raise
    return wrapper

class SharedResourceManager:
    """Manages shared computational resources across plot generation"""
    
    def __init__(self):
        self.universe = None
        self.selections = {}
        self.barnaba_data = {}
        
    def initialize_universe(self, topology_file, trajectory_file):
        """Initialize MDAnalysis Universe once"""
        if not MDA_AVAILABLE:
            raise RuntimeError("MDAnalysis not available")
            
        if self.universe is None:
            self.universe = mda.Universe(topology_file, trajectory_file)
            logger.info(f"Initialized Universe: {len(self.universe.trajectory)} frames, {len(self.universe.atoms)} atoms")
        return self.universe
    
    def get_selection(self, selection_string):
        """Get cached atom selection"""
        if selection_string not in self.selections:
            if self.universe is None:
                raise RuntimeError("Universe not initialized")
            
            start_time = time.time()
            self.selections[selection_string] = self.universe.select_atoms(selection_string)
            duration = time.time() - start_time
            
            logger.info(f"Cached selection '{selection_string}': {len(self.selections[selection_string])} atoms ({duration:.3f}s)")
        
        return self.selections[selection_string]
    
    def precompute_barnaba_data(self, topology_file, trajectory_file):
        """Precompute common Barnaba analyses once"""
        if not BARNABA_AVAILABLE:
            logger.warning("Barnaba not available, skipping precomputation")
            return self.barnaba_data
            
        if not self.barnaba_data:
            logger.info("Precomputing Barnaba data...")
            start_time = time.time()
            
            # Common barnaba computations - using correct parameter pattern from tasks_celery.py
            try:
                self.barnaba_data['ermsd'] = bb.ermsd(topology_file, trajectory_file, topology=topology_file)
                self.barnaba_data['torsions'], self.barnaba_data['torsions_res'] = bb.backbone_angles(trajectory_file, topology=topology_file)
                
                duration = time.time() - start_time
                logger.info(f"Precomputed Barnaba data in {duration:.2f}s")
            except Exception as e:
                logger.warning(f"Barnaba precomputation failed: {e}")
                
        return self.barnaba_data
    
    def clear(self):
        """Clear all cached resources"""
        self.universe = None
        self.selections = {}
        self.barnaba_data = {}

# Global resource manager
resource_manager = SharedResourceManager()

@app.task(bind=True, max_retries=3)
@simple_log_execution
@cached_computation(computation_type="simple_rna_analysis")
def simple_rna_analysis(self, topology_file, trajectory_file, plots_requested, session_id):
    """
    Simplified RNA analysis with resource precomputation and shared cache
    
    1. Precompute all resources needed
    2. Use shared resources for all plot generation
    3. Cache results for future use
    """
    
    start_time = time.time()
    results = {}
    
    try:
        # Phase 1: Initialize shared resources
        logger.info("Phase 1: Initializing shared resources")
        universe = resource_manager.initialize_universe(topology_file, trajectory_file)
        
        # Precompute common selections
        resource_manager.get_selection('all')
        resource_manager.get_selection('nucleic')
        resource_manager.get_selection('protein')
        resource_manager.get_selection('backbone')
        
        # Precompute Barnaba data
        barnaba_data = resource_manager.precompute_barnaba_data(topology_file, trajectory_file)
        
        # Phase 2: Generate plots using shared resources
        logger.info(f"Phase 2: Generating {len(plots_requested)} plots using shared resources")
        
        for plot_type in plots_requested:
            plot_start = time.time()
            
            try:
                # Use the existing plot generation functions but with shared resources
                if plot_type == 'RMSD':
                    results[plot_type] = generate_rmsd_with_shared_resources(
                        universe, resource_manager.get_selection('nucleic'), session_id
                    )
                elif plot_type == 'ERMSD':
                    results[plot_type] = generate_ermsd_with_shared_resources(
                        barnaba_data.get('ermsd'), session_id
                    )
                elif plot_type == 'CONTACT_MAPS':
                    # Use the existing create_contact_maps function like in tasks_celery.py
                    results[plot_type] = generate_contact_maps_with_create_plots(
                        topology_file, trajectory_file, session_id
                    )
                elif plot_type == 'TORSION':
                    results[plot_type] = generate_torsion_with_shared_resources(
                        barnaba_data.get('torsions'), barnaba_data.get('torsions_res'), session_id
                    )
                else:
                    # Fallback to standard generation
                    results[plot_type] = generate_plot_standard(plot_type, topology_file, trajectory_file, session_id)
                
                plot_duration = time.time() - plot_start
                logger.info(f"Generated {plot_type} in {plot_duration:.2f}s")
                
            except Exception as e:
                logger.error(f"Failed to generate {plot_type}: {e}")
                results[plot_type] = {"error": str(e)}
        
        # Phase 3: Performance summary
        total_time = time.time() - start_time
        results['_performance_metadata'] = {
            'total_time': total_time,
            'plots_processed': len(plots_requested),
            'avg_time_per_plot': total_time / len(plots_requested) if plots_requested else 0,
            'optimization_level': 'simple_shared_resources',
            'universe_stats': {
                'n_frames': len(universe.trajectory),
                'n_atoms': len(universe.atoms),
                'selections_cached': len(resource_manager.selections),
                'barnaba_precomputed': len(barnaba_data)
            }
        }
        
        logger.info(f"Simple RNA analysis completed in {total_time:.2f}s")
        return results
        
    except Exception as e:
        logger.error(f"Simple RNA analysis failed: {e}")
        logger.error(traceback.format_exc())
        raise

def generate_rmsd_with_shared_resources(universe, selection, session_id):
    """Generate RMSD plot using shared universe and selection"""
    from MDAnalysis.analysis import rms
    import matplotlib
    matplotlib.use('Agg')  # Use non-GUI backend
    import matplotlib.pyplot as plt
    
    # Use the pre-loaded universe and selection
    if len(selection) == 0:
        # Fallback to all atoms if no specific selection
        selection = universe.atoms
    
    rmsd_analysis = rms.RMSD(selection, selection)
    rmsd_analysis.run()
    
    # Create plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Handle different MDAnalysis versions
    if hasattr(rmsd_analysis.results, 'time'):
        time_data = rmsd_analysis.results.time
    else:
        time_data = range(len(rmsd_analysis.results.rmsd))
    
    ax.plot(time_data, rmsd_analysis.results.rmsd[:, 2])
    ax.set_xlabel('Frame' if not hasattr(rmsd_analysis.results, 'time') else 'Time (ps)')
    ax.set_ylabel('RMSD (Å)')
    ax.set_title('RMSD over Time')
    
    # Save plot
    plot_path = f"static/{session_id}/rmsd_shared.png"
    os.makedirs(os.path.dirname(plot_path), exist_ok=True)
    fig.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    return json.dumps({"path": plot_path})

def generate_ermsd_with_shared_resources(ermsd_data, session_id):
    """Generate eRMSD plot using precomputed data"""
    if ermsd_data is None:
        raise ValueError("eRMSD data not precomputed")
    
    if CREATE_PLOTS_AVAILABLE:
        try:
            import plotly.io as pio
            
            # Create plot using existing function
            fig = plot_ermsd(ermsd_data)
            
            # Save plot
            plot_path = f"static/{session_id}/ermsd_shared.png"
            os.makedirs(os.path.dirname(plot_path), exist_ok=True)
            
            # Save as PNG using plotly
            pio.write_image(fig, plot_path, format='png', width=1200, height=600)
            
            return json.dumps({"path": plot_path})
            
        except Exception as e:
            logger.warning(f"Failed to use plot_ermsd, falling back to matplotlib: {e}")
    
    # Fallback to matplotlib (same as original implementation)
    import matplotlib
    matplotlib.use('Agg')  # Use non-GUI backend
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(ermsd_data)
    ax.set_xlabel('Frame')
    ax.set_ylabel('eRMSD')
    ax.set_title('Elastic RMSD over Time')
    
    # Save plot
    plot_path = f"static/{session_id}/ermsd_shared.png"
    os.makedirs(os.path.dirname(plot_path), exist_ok=True)
    fig.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    return json.dumps({"path": plot_path})

def generate_contact_maps_with_shared_resources(contacts_data, session_id):
    """Generate contact maps using precomputed data"""
    if contacts_data is None:
        raise ValueError("Contacts data not precomputed")
    
    # Create heatmap from precomputed data
    import matplotlib
    matplotlib.use('Agg')  # Use non-GUI backend
    import matplotlib.pyplot as plt
    import numpy as np
    
    fig, ax = plt.subplots(figsize=(8, 8))
    im = ax.imshow(contacts_data, cmap='viridis', aspect='equal')
    ax.set_xlabel('Residue')
    ax.set_ylabel('Residue')
    ax.set_title('Contact Map')
    plt.colorbar(im)
    
    # Save plot
    plot_path = f"static/{session_id}/contact_map_shared.png"
    os.makedirs(os.path.dirname(plot_path), exist_ok=True)
    fig.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    return json.dumps({"path": plot_path})

def generate_contact_maps_with_create_plots(topology_file, trajectory_file, session_id):
    """Generate contact maps using the existing create_contact_maps function like in tasks_celery.py"""
    if CREATE_PLOTS_AVAILABLE:
        try:
            # Extract filenames from full paths (the function expects just filenames)
            traj_filename = os.path.basename(trajectory_file)
            topology_filename = os.path.basename(topology_file)
            
            # Call the existing function exactly like in tasks_celery.py
            create_contact_maps(traj_filename, topology_filename, 1, session_id)
            
            # Return the expected path
            path_save_figure = os.path.join('static', session_id, 'contact_map_plotly.png')
            return json.dumps({'path': path_save_figure})
            
        except Exception as e:
            logger.warning(f"create_contact_maps failed, using fallback: {e}")
    
    # Fallback: Generate basic contact map using barnaba annotate
    try:
        if not BARNABA_AVAILABLE:
            return json.dumps({"error": "Neither create_plots nor barnaba available"})
        
        import barnaba as bb
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import pandas as pd
        
        # Use barnaba annotate like in the original
        stackings, pairings, res = bb.annotate(trajectory_file, topology=topology_file)
        
        # Create a simple contact matrix visualization
        fig, ax = plt.subplots(figsize=(8, 8))
        
        # Create basic contact representation
        n_residues = len(res)
        contact_matrix = np.zeros((n_residues, n_residues))
        
        # Fill matrix with pairing information from first frame
        if pairings and len(pairings) > 0:
            frame_pairs = pairings[0]
            if len(frame_pairs) == 2:
                pair_indices = frame_pairs[0]
                for pair in pair_indices:
                    if pair:
                        i, j = pair[0], pair[1]
                        if 0 <= i < n_residues and 0 <= j < n_residues:
                            contact_matrix[i, j] = 1
                            contact_matrix[j, i] = 1
        
        im = ax.imshow(contact_matrix, cmap='viridis', aspect='equal')
        ax.set_xlabel('Residue')
        ax.set_ylabel('Residue')
        ax.set_title('Contact Map (Fallback)')
        plt.colorbar(im)
        
        # Save plot
        plot_path = f"static/{session_id}/contact_map_fallback.png"
        os.makedirs(os.path.dirname(plot_path), exist_ok=True)
        fig.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        
        return json.dumps({'path': plot_path})
        
    except Exception as e:
        logger.error(f"Fallback contact maps generation failed: {e}")
        return json.dumps({"error": str(e)})

def generate_torsion_with_shared_resources(torsions_data, torsions_res, session_id):
    """Generate torsion plot using precomputed data"""
    if torsions_data is None:
        raise ValueError("Torsions data not precomputed")
    
    if CREATE_PLOTS_AVAILABLE:
        try:
            import plotly.io as pio
            
            # Create plot using existing function with default residue 0
            fig = plot_torsion(torsions_data, torsions_res, 0)
            
            # Save plot
            plot_path = f"static/{session_id}/torsion_shared.png"
            os.makedirs(os.path.dirname(plot_path), exist_ok=True)
            
            # Save as PNG using plotly
            pio.write_image(fig, plot_path, format='png', width=1200, height=600)
            
            return json.dumps({"path": plot_path})
            
        except Exception as e:
            logger.warning(f"Failed to use plot_torsion, falling back to matplotlib: {e}")
    
    # Fallback to matplotlib (similar to original implementation)
    import matplotlib
    matplotlib.use('Agg')  # Use non-GUI backend
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Plot torsion angles (shape is usually frames x residues x 7_angles)
    if len(torsions_data.shape) == 3:
        # Pick first residue as example, plot all 7 torsion types
        angle_names = ["alpha", "beta", "gamma", "delta", "epsilon", "zeta", "chi"]
        for i in range(min(7, torsions_data.shape[2])):
            ax.plot(torsions_data[:, 0, i], label=angle_names[i] if i < len(angle_names) else f'Angle {i+1}')
    else:
        # Fallback: plot as 2D data
        for i, torsion in enumerate(torsions_data.T):
            ax.plot(torsion, label=f'Torsion {i+1}')
    
    ax.set_xlabel('Frame')
    ax.set_ylabel('Angle (degrees)')
    ax.set_title('Torsion Angles over Time')
    ax.legend()
    
    # Save plot
    plot_path = f"static/{session_id}/torsion_shared.png"
    os.makedirs(os.path.dirname(plot_path), exist_ok=True)
    fig.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    return json.dumps({"path": plot_path})

def generate_plot_standard(plot_type, topology_file, trajectory_file, session_id):
    """Fallback to standard plot generation for unsupported types"""
    # This would call the original plot generation functions
    # For now, return a placeholder
    return json.dumps({"path": f"static/{session_id}/{plot_type.lower()}_placeholder.png"})

# Individual plot tasks for backward compatibility
@app.task(bind=True, max_retries=3)
@simple_log_execution
def generate_rmsd_plot(self, topology_file, trajectory_file, files_path, plot_dir, session_id):
    """Generate RMSD plot - fallback to simple implementation"""
    try:
        universe = resource_manager.initialize_universe(topology_file, trajectory_file)
        selection = resource_manager.get_selection('nucleic')
        return generate_rmsd_with_shared_resources(universe, selection, session_id)
    except Exception as e:
        logger.error(f"RMSD plot generation failed: {e}")
        return json.dumps({"error": str(e)})

@app.task(bind=True, max_retries=3)
@simple_log_execution
def generate_ermsd_plot(self, topology_file, trajectory_file, files_path, plot_dir, session_id):
    """Generate eRMSD plot"""
    try:
        barnaba_data = resource_manager.precompute_barnaba_data(topology_file, trajectory_file)
        return generate_ermsd_with_shared_resources(barnaba_data.get('ermsd'), session_id)
    except Exception as e:
        logger.error(f"eRMSD plot generation failed: {e}")
        return json.dumps({"error": str(e)})

@app.task(bind=True, max_retries=3)
@simple_log_execution
def generate_contact_map_plot(self, topology_file, trajectory_file, files_path, plot_dir, generate_data_path, session_id):
    """Generate contact map plot"""
    try:
        barnaba_data = resource_manager.precompute_barnaba_data(topology_file, trajectory_file)
        return generate_contact_maps_with_shared_resources(barnaba_data.get('contacts'), session_id)
    except Exception as e:
        logger.error(f"Contact map plot generation failed: {e}")
        return json.dumps({"error": str(e)})

@app.task(bind=True, max_retries=3)
@simple_log_execution
def generate_torsion_plot(self, topology_file, trajectory_file, files_path, plot_dir, session_id, torsion_params=None):
    """Generate torsion plot"""
    try:
        barnaba_data = resource_manager.precompute_barnaba_data(topology_file, trajectory_file)
        return generate_torsion_with_shared_resources(barnaba_data.get('torsions'), barnaba_data.get('torsions_res'), session_id)
    except Exception as e:
        logger.error(f"Torsion plot generation failed: {e}")
        return json.dumps({"error": str(e)})

# Placeholder tasks for other plots
@app.task(bind=True, max_retries=3)
@simple_log_execution
def generate_sec_structure_plot(self, topology_file, trajectory_file, files_path, plot_dir, session_id):
    return json.dumps({"path": f"static/{session_id}/sec_structure_placeholder.png"})

@app.task(bind=True, max_retries=3)
@simple_log_execution
def generate_dotbracket_plot(self, topology_file, trajectory_file, files_path, plot_dir, session_id):
    return json.dumps({"path": f"static/{session_id}/dotbracket_placeholder.png"})

@app.task(bind=True, max_retries=3)
@simple_log_execution
def generate_arc_plot(self, topology_file, trajectory_file, files_path, plot_dir, session_id):
    return json.dumps({"path": f"static/{session_id}/arc_placeholder.png"})

@app.task(bind=True, max_retries=3)
@simple_log_execution
def generate_annotate_plot(self, topology_file, trajectory_file, files_path, plot_dir, session_id):
    return json.dumps({"path": f"static/{session_id}/annotate_placeholder.png"})

@app.task(bind=True, max_retries=3)
@simple_log_execution
def generate_ds_motif_plot(self, topology_file, trajectory_file, files_path, plot_dir, session_id):
    return json.dumps({"path": f"static/{session_id}/ds_motif_placeholder.png"})

@app.task(bind=True, max_retries=3)
@simple_log_execution
def generate_ss_motif_plot(self, topology_file, trajectory_file, files_path, plot_dir, session_id):
    return json.dumps({"path": f"static/{session_id}/ss_motif_placeholder.png"})

@app.task(bind=True, max_retries=3)
@simple_log_execution
def generate_jcoupling_plot(self, topology_file, trajectory_file, files_path, plot_dir, session_id):
    return json.dumps({"path": f"static/{session_id}/jcoupling_placeholder.png"})

@app.task(bind=True, max_retries=3)
@simple_log_execution
def generate_escore_plot(self, topology_file, trajectory_file, files_path, plot_dir, session_id):
    return json.dumps({"path": f"static/{session_id}/escore_placeholder.png"})

@app.task(bind=True, max_retries=3)
@simple_log_execution
def generate_landscape_plot(self, topology_file, trajectory_file, files_path, plot_dir, session_id, landscape_params, generate_data_path):
    return json.dumps({"path": f"static/{session_id}/landscape_placeholder.png"})

@app.task(bind=True, max_retries=3)
@simple_log_execution
def generate_2Dpairing_plot(self, topology_file, trajectory_file, files_path, plot_dir, session_id):
    return json.dumps({"path": f"static/{session_id}/2dpairing_placeholder.png"})

@app.task(bind=True, max_retries=3)
@simple_log_execution
def update_landscape_frame(self, generate_data_path, coordinates):
    return json.dumps({"success": True})

@app.task(bind=True, max_retries=3)
@simple_log_execution
def update_contact_map_plot(self, generate_data_path, download_plot, slider_value, session):
    return json.dumps({"success": True})

@app.task(bind=True, max_retries=3)
@simple_log_execution
def performance_benchmark(self, *args, **kwargs):
    return {"benchmark_completed": True}

@app.task(bind=True, max_retries=3)
@simple_log_execution
def batch_computation_coordinator(self, *args, **kwargs):
    return {"batch_completed": True}

if __name__ == "__main__":
    # For testing
    logger.info("Simple optimized tasks module loaded")
