"""
Clean Celery Tasks for RNA Analysis
Two-phase approach: Metrics computation then Plot generation
"""

import os
import time
import pickle
import logging
from celery import Celery
from functools import wraps
from create_plots import *
from plotly.io import to_json

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Celery app configuration
app = Celery('rna_tasks')
app.conf.update(
    broker_url='redis://localhost:6379/0',
    result_backend='redis://localhost:6379/0',
    task_serializer='json',
    result_serializer='json',
    accept_content=['json'],
    task_acks_late=False,
    worker_prefetch_multiplier=1,
    result_expires=3600,
)

# Check dependencies
try:
    import barnaba as bb
    BARNABA_AVAILABLE = True
except ImportError:
    BARNABA_AVAILABLE = False
    logger.warning("Barnaba not available")

def plotly_to_json(fig):
    return to_json(fig, validate=False, engine="orjson")
    
def log_task(func):
    """Simple task logging decorator"""
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

# Phase 1: Metric Computation Tasks
@app.task(bind=True, max_retries=3)
@log_task
def compute_rmsd(self, *args):
    """Compute RMSD metric"""
    # Handle chain arguments
    if len(args) == 4:  # previous_result, topology_file, trajectory_file, session_id
        _, topology_file, trajectory_file, session_id = args
    else:  # topology_file, trajectory_file, session_id
        topology_file, trajectory_file, session_id = args
        
    if not BARNABA_AVAILABLE:
        raise ImportError("Barnaba not available")
    
    import barnaba as bb
    rmsd_result = bb.rmsd(topology_file, trajectory_file, topology=topology_file, heavy_atom=True)
    
    # Save to session directory
    session_dir = os.path.join("static", session_id)
    os.makedirs(session_dir, exist_ok=True)
    
    rmsd_path = os.path.join(session_dir, "rmsd_data.pkl")
    with open(rmsd_path, 'wb') as f:
        pickle.dump(rmsd_result, f)
    
    logger.info(f"RMSD computed and saved to {rmsd_path}")
    return {"metric": "rmsd", "status": "success", "path": rmsd_path}

@app.task(bind=True, max_retries=3)
@log_task
def compute_ermsd(self, *args):
    """Compute eRMSD metric"""
    # Handle chain arguments
    if len(args) == 4:  # previous_result, topology_file, trajectory_file, session_id
        _, topology_file, trajectory_file, session_id = args
    else:  # topology_file, trajectory_file, session_id
        topology_file, trajectory_file, session_id = args
        
    if not BARNABA_AVAILABLE:
        raise ImportError("Barnaba not available")
    
    import barnaba as bb
    ermsd_result = bb.ermsd(topology_file, trajectory_file, topology=topology_file)
    
    # Save to session directory
    session_dir = os.path.join("static", session_id)
    os.makedirs(session_dir, exist_ok=True)
    
    ermsd_path = os.path.join(session_dir, "ermsd_data.pkl")
    with open(ermsd_path, 'wb') as f:
        pickle.dump(ermsd_result, f)
    
    logger.info(f"eRMSD computed and saved to {ermsd_path}")
    return {"metric": "ermsd", "status": "success", "path": ermsd_path}

@app.task(bind=True, max_retries=3)
@log_task
def compute_annotate(self, *args):
    """Compute annotate metric"""
    # Handle chain arguments
    if len(args) == 4:  # previous_result, topology_file, trajectory_file, session_id
        _, topology_file, trajectory_file, session_id = args
    else:  # topology_file, trajectory_file, session_id
        topology_file, trajectory_file, session_id = args
        
    if not BARNABA_AVAILABLE:
        raise ImportError("Barnaba not available")
    
    import barnaba as bb
    annotate_result = bb.annotate(trajectory_file, topology=topology_file)
    
    # Save to session directory
    session_dir = os.path.join("static", session_id)
    os.makedirs(session_dir, exist_ok=True)
    
    annotate_path = os.path.join(session_dir, "annotate_data.pkl")
    with open(annotate_path, 'wb') as f:
        pickle.dump(annotate_result, f)
    
    logger.info(f"Annotate computed and saved to {annotate_path}")
    return {"metric": "annotate", "status": "success", "path": annotate_path}

# Phase 2: Plot Generation Tasks
@app.task(bind=True, max_retries=3)
@log_task
def generate_rmsd_plot(self, topology_file, trajectory_file, files_path, plot_dir, session_id):
    """Generate RMSD plot using pre-computed data"""
    try:
        # Load pre-computed RMSD data
        rmsd_path = os.path.join("static", session_id, "rmsd_data.pkl")
        if os.path.exists(rmsd_path):
            with open(rmsd_path, 'rb') as f:
                rmsd = pickle.load(f)
            print(f"LOADED RMSD FROM SAVED DATA")
        else:
            # Fallback: compute if not available
            print(f"USED FALLBACK")
            if not BARNABA_AVAILABLE:
                raise ImportError("Barnaba not available and no pre-computed data")
            import barnaba as bb
            rmsd = bb.rmsd(topology_file, trajectory_file, topology=topology_file, heavy_atom=True)

        # Create plot using create_plots functions
        #try:
        logging.info("Default plots are used")
        fig = plot_rmsd(rmsd)
        print(fig)
    
        # Save data and plot
        import pandas as pd
        rmsd_df = pd.DataFrame({"RMSD": rmsd})
        rmsd_df.to_csv(os.path.join(files_path, "rmsd_values.csv"), index=False)
        fig.write_html(os.path.join(plot_dir, "rmsd_plot.html"))

        # Convert to JSON
        plotly_data = plotly_to_json(fig)
        return plotly_data
        """
        except ImportError:
            # Fallback to simple matplotlib plot
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.plot(rmsd)
            ax.set_xlabel('Frame')
            ax.set_ylabel('RMSD (Å)')
            ax.set_title('RMSD Analysis')
            
            plot_path = os.path.join(plot_dir, "rmsd_plot.png")
            os.makedirs(plot_dir, exist_ok=True)
            plt.savefig(plot_path, dpi=150, bbox_inches='tight')
            plt.close()
            
            return {"path": plot_path, "status": "success"}
        """
        
    except Exception as e:
        logger.error(f"RMSD plot generation failed: {e}")
        raise

@app.task(bind=True, max_retries=3)
@log_task
def generate_ermsd_plot(self, topology_file, trajectory_file, files_path, plot_dir, session_id):
    """Generate eRMSD plot using pre-computed data"""
    try:
        # Load pre-computed eRMSD data
        ermsd_path = os.path.join("static", session_id, "ermsd_data.pkl")
        if os.path.exists(ermsd_path):
            with open(ermsd_path, 'rb') as f:
                ermsd = pickle.load(f)
            print(f"LOADED eRMSD FROM SAVED DATA")
        else:
            # Fallback: compute if not available
            print(f"USED FALLBACK")
            if not BARNABA_AVAILABLE:
                raise ImportError("Barnaba not available and no pre-computed data")
            import barnaba as bb
            ermsd = bb.ermsd(topology_file, trajectory_file, topology=topology_file)

        # Create plot using create_plots functions
        try:
            fig = plot_ermsd(ermsd)
            
            # Save data and plot
            import pandas as pd
            ermsd_df = pd.DataFrame({"ERMSD": ermsd})
            ermsd_df.to_csv(os.path.join(files_path, "ermsd_values.csv"), index=False)
            fig.write_html(os.path.join(plot_dir, "ermsd_plot.html"))

            # Convert to JSON
            plotly_data = plotly_to_json(fig)
            return plotly_data
            
        except ImportError:
            # Fallback to simple matplotlib plot
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.plot(ermsd)
            ax.set_xlabel('Frame')
            ax.set_ylabel('eRMSD')
            ax.set_title('eRMSD Analysis')
            
            plot_path = os.path.join(plot_dir, "ermsd_plot.png")
            os.makedirs(plot_dir, exist_ok=True)
            plt.savefig(plot_path, dpi=150, bbox_inches='tight')
            plt.close()
            
            return {"path": plot_path, "status": "success"}
        
    except Exception as e:
        logger.error(f"eRMSD plot generation failed: {e}")
        raise

# Placeholder tasks for other plots that don't need metrics
@app.task(bind=True, max_retries=3)
@log_task  
def generate_torsion_plot(self, topology_file, trajectory_file, files_path, plot_dir, session_id, torsion_residue=0):
    """Generate torsion plot"""
    try:
        if not BARNABA_AVAILABLE:
            raise ImportError("Barnaba not available")
            
        import barnaba as bb
        angles, res = bb.backbone_angles(trajectory_file, topology=topology_file)
        logger.info(f"Calculated torsion angles for {len(res)} residues")
        
        try:
            from create_plots import plot_torsion, plot_torsion_enhanced, plotly_to_json
            
            # Handle backward compatibility
            if isinstance(torsion_residue, dict):
                fig = plot_torsion_enhanced(angles, res, torsion_residue)
            else:
                # Old single residue format
                fig = plot_torsion(angles, res, torsion_residue)
                
            fig.write_html(os.path.join(plot_dir, "torsion_plot.html"))
            
            # Save data
            import pandas as pd
            angles_df = pd.DataFrame(angles.reshape(-1, angles.shape[-1]), 
                                   columns=["alpha", "beta", "gamma", "delta", "epsilon", "zeta", "chi"])
            angles_df.to_csv(os.path.join(files_path, "torsion_angles.csv"), index=False)
            
            plotly_data = plotly_to_json(fig)
            return plotly_data
            
        except ImportError:
            return {"path": f"static/{session_id}/torsion_plot.png", "status": "fallback"}
            
    except Exception as exc:
        logger.error(f"Torsion calculation failed: {str(exc)}")
        raise exc

@app.task(bind=True, max_retries=3)
@log_task
def generate_sec_structure_plot(self, topology_file, trajectory_file, files_path, plot_dir, session_id):
    """Generate secondary structure plot"""
    return {"path": f"static/{session_id}/sec_structure_plot.png", "status": "placeholder"}

@app.task(bind=True, max_retries=3)
@log_task
def generate_dotbracket_plot(self, topology_file, trajectory_file, files_path, plot_dir, session_id):
    """Generate dot-bracket plot"""
    return {"path": f"static/{session_id}/dotbracket_plot.png", "status": "placeholder"}

@app.task(bind=True, max_retries=3)
@log_task
def generate_arc_plot(self, topology_file, trajectory_file, files_path, plot_dir, session_id):
    """Generate arc plot"""
    return {"path": f"static/{session_id}/arc_plot.png", "status": "placeholder"}

@app.task(bind=True, max_retries=3)
@log_task
def generate_contact_map_plot(self, topology_file, trajectory_file, files_path, plot_dir, generate_data_path, session_id):
    """Generate contact map plot"""
    return {"path": f"static/{session_id}/contact_map_plot.png", "status": "placeholder"}

@app.task(bind=True, max_retries=3)
@log_task
def generate_annotate_plot(self, topology_file, trajectory_file, files_path, plot_dir, session_id):
    """Generate annotate plot using pre-computed data"""
    try:
        # Load pre-computed annotate data
        annotate_path = os.path.join("static", session_id, "annotate_data.pkl")
        if os.path.exists(annotate_path):
            with open(annotate_path, 'rb') as f:
                annotate_data = pickle.load(f)
        else:
            # Fallback: compute if not available
            if not BARNABA_AVAILABLE:
                raise ImportError("Barnaba not available and no pre-computed data")
            import barnaba as bb
            annotate_data = bb.annotate(trajectory_file, topology=topology_file)
        
        # Simple placeholder for annotate plot
        return {"path": f"static/{session_id}/annotate_plot.png", "status": "success"}
        
    except Exception as e:
        logger.error(f"Annotate plot generation failed: {e}")
        raise

@app.task(bind=True, max_retries=3)
@log_task
def generate_ds_motif_plot(self, topology_file, trajectory_file, files_path, plot_dir, session_id):
    """Generate double-strand motif plot"""
    return {"path": f"static/{session_id}/ds_motif_plot.png", "status": "placeholder"}

@app.task(bind=True, max_retries=3)
@log_task
def generate_ss_motif_plot(self, topology_file, trajectory_file, files_path, plot_dir, session_id):
    """Generate single-strand motif plot"""
    return {"path": f"static/{session_id}/ss_motif_plot.png", "status": "placeholder"}

@app.task(bind=True, max_retries=3)
@log_task
def generate_jcoupling_plot(self, topology_file, trajectory_file, files_path, plot_dir, session_id):
    """Generate J-coupling plot"""
    return {"path": f"static/{session_id}/jcoupling_plot.png", "status": "placeholder"}

@app.task(bind=True, max_retries=3)
@log_task
def generate_escore_plot(self, topology_file, trajectory_file, files_path, plot_dir, session_id):
    """Generate e-score plot"""
    return {"path": f"static/{session_id}/escore_plot.png", "status": "placeholder"}

@app.task(bind=True, max_retries=3)
@log_task
def generate_landscape_plot(self, topology_file, trajectory_file, files_path, plot_dir, session_id, landscape_params, generate_data_path):
    """Generate landscape plot"""
    return {"path": f"static/{session_id}/landscape_plot.png", "status": "placeholder"}

@app.task(bind=True, max_retries=3)
@log_task
def generate_2Dpairing_plot(self, topology_file, trajectory_file, files_path, plot_dir, session_id):
    """Generate 2D pairing plot"""
    return {"path": f"static/{session_id}/2dpairing_plot.png", "status": "placeholder"}

if __name__ == "__main__":
    logger.info("Clean RNA analysis tasks module loaded")
