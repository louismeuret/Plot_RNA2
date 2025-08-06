"""
Ultra-Optimized Celery tasks using advanced MDAnalysis techniques
Implements 2024 cutting-edge optimizations for maximum performance
"""

from create_plots import *
import json
import time
from celery import Celery
import os
import pandas as pd
import barnaba as bb
import numpy as np
import logging
from functools import wraps
from typing import Optional, Dict, Any, List
import traceback
import pickle
from itertools import combinations
import mdtraj as md
import energy_3dplot

# Import our optimization systems
from computation_cache import computation_cache, cached_computation
from mdanalysis_optimizations import OptimizedTrajectoryProcessor, RNAAnalysisOptimizer, OptimizationConfig

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Configure Celery
app3 = Celery('ultra_optimized_tasks')
app3.conf.update(
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

def ultra_log_task_execution(func):
    """Enhanced task execution logger with optimization tracking"""
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        task_name = func.__name__
        session_id = args[-1] if args else 'unknown'
        
        logger.info(f"Starting ULTRA-optimized task {task_name} for session {session_id}")
        start_time = time.time()
        
        try:
            result = func(self, *args, **kwargs)
            execution_time = time.time() - start_time
            logger.info(f"ULTRA-optimized task {task_name} completed for session {session_id} in {execution_time:.2f}s")
            return result
        except Exception as exc:
            execution_time = time.time() - start_time
            logger.error(f"ULTRA-optimized task {task_name} failed for session {session_id} after {execution_time:.2f}s: {str(exc)}")
            logger.error(f"Traceback: {traceback.format_exc()}")
            raise
    return wrapper

def plotly_to_json(fig):
    from plotly.io import to_json
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
    logger.info(f"Number of native contacts: {len(native_contacts)}")

    # now compute these distances for the whole trajectory
    r = md.compute_distances(traj, native_contacts)
    # and recompute them for just the native state
    r0 = md.compute_distances(native[0], native_contacts)

    q = np.mean(1.0 / (1 + np.exp(BETA_CONST * (r - LAMBDA_CONST * r0))), axis=1)
    return q

@app3.task(bind=True, max_retries=3)
@ultra_log_task_execution
def ultra_optimized_rna_analysis(self, native_pdb_path, traj_xtc_path, download_path, plot_path, 
                                session_id, selected_plots, optimization_level='maximum'):
    """
    Ultra-optimized RNA analysis using advanced MDAnalysis techniques
    Processes all requested plots with maximum efficiency
    """
    try:
        total_start_time = time.time()
        
        # Configure optimization level
        config = OptimizationConfig(
            use_parallel_processing=True,
            n_workers=min(8, len(selected_plots)),  # Adaptive worker count
            chunk_size=50,  # Smaller chunks for better parallelization
            preload_selections=True,
            enable_caching=True
        )
        
        if optimization_level == 'conservative':
            config.use_parallel_processing = False
            config.n_workers = 2
        
        logger.info(f"Using optimization config: {config}")
        
        # Initialize optimized processor
        processor = OptimizedTrajectoryProcessor(native_pdb_path, traj_xtc_path, config)
        
        # Get optimization stats
        stats = processor.get_optimization_stats()
        logger.info(f"Processor stats: {stats}")
        
        results = {}
        
        # Phase 1: Precompute common selections
        logger.info("Phase 1: Precomputing atom selections")
        selections = {
            'protein': 'protein',
            'backbone': 'backbone',
            'nucleic': 'nucleic',
            'all': 'all'
        }
        processor.precompute_selections(selections)
        
        # Phase 2: Use RNA-specific optimizer for Barnaba calculations
        logger.info("Phase 2: Optimized Barnaba calculations")
        barnaba_plots = [plot for plot in selected_plots if plot in 
                        {'RMSD', 'ERMSD', 'TORSION', 'SEC_STRUCTURE', 'DOTBRACKET', 'ARC', 'CONTACT_MAPS'}]
        
        if barnaba_plots:
            barnaba_results = RNAAnalysisOptimizer.optimize_barnaba_calculations(
                processor, barnaba_plots
            )
            
            # Process Barnaba results into plot format
            for plot_type in barnaba_plots:
                plot_start = time.time()
                
                if plot_type == 'RMSD' and 'rmsd' in barnaba_results:
                    rmsd_data = barnaba_results['rmsd']
                    fig = plot_rmsd(rmsd_data)
                    fig.write_html(os.path.join(plot_path, "rmsd_plot.html"))
                    
                    # Save data
                    rmsd_df = pd.DataFrame({"RMSD": rmsd_data})
                    rmsd_df.to_csv(os.path.join(download_path, "rmsd_values.csv"), index=False)
                    
                    results['RMSD'] = plotly_to_json(fig)
                
                elif plot_type == 'ERMSD' and 'ermsd' in barnaba_results:
                    ermsd_data = barnaba_results['ermsd']
                    fig = plot_ermsd(ermsd_data)
                    fig.write_html(os.path.join(plot_path, "ermsd_plot.html"))
                    
                    # Save data
                    ermsd_df = pd.DataFrame({"ERMSD": ermsd_data})
                    ermsd_df.to_csv(os.path.join(download_path, "ermsd_values.csv"), index=False)
                    
                    results['ERMSD'] = plotly_to_json(fig)
                
                elif plot_type == 'TORSION' and 'torsion' in barnaba_results:
                    angles, res = barnaba_results['torsion']
                    
                    # Get torsion parameters from cache or defaults
                    torsion_params = computation_cache.get(session_id, 'torsion_params', ()) or {
                        'torsionMode': 'single',
                        'torsionResidue': 0,
                        'torsionResidues': [],
                        'torsionAngles': ['alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'chi']
                    }
                    
                    if isinstance(torsion_params, dict):
                        fig = plot_torsion_enhanced(angles, res, torsion_params)
                    else:
                        fig = plot_torsion(angles, res, 0)
                    
                    fig.write_html(os.path.join(plot_path, "torsion_plot.html"))
                    
                    # Save data
                    angles_df = pd.DataFrame(
                        angles.reshape(-1, angles.shape[-1]), 
                        columns=["alpha", "beta", "gamma", "delta", "epsilon", "zeta", "chi"]
                    )
                    angles_df.to_csv(os.path.join(download_path, "torsion_angles.csv"), index=False)
                    
                    results['TORSION'] = plotly_to_json(fig)
                
                elif plot_type in ['SEC_STRUCTURE', 'DOTBRACKET', 'ARC'] and 'annotate_data' in barnaba_results:
                    stackings, pairings, res = barnaba_results['annotate_data']
                    
                    if plot_type == 'SEC_STRUCTURE':
                        dotbracket_data = bb.dot_bracket(pairings, res)
                        results['SEC_STRUCTURE'] = [dotbracket_data[0], res.strip()]
                    
                    elif plot_type == 'DOTBRACKET':
                        if 'dotbracket' in barnaba_results:
                            dotbracket_data = barnaba_results['dotbracket'][0]
                        else:
                            dotbracket_data = bb.dot_bracket(pairings, res)[0]
                        
                        fig = plot_dotbracket(dotbracket_data)
                        fig.write_html(os.path.join(plot_path, "dotbracket_timeline_plot.html"))
                        results['DOTBRACKET'] = plotly_to_json(fig)
                    
                    elif plot_type == 'ARC':
                        # Get trajectory dotbracket
                        if 'dotbracket' in barnaba_results:
                            dotbracket_data = barnaba_results['dotbracket'][0]
                        else:
                            dotbracket_data = bb.dot_bracket(pairings, res)[0]
                        
                        # Get native dotbracket
                        native_cached = computation_cache.get(session_id, 'native_dotbracket', (native_pdb_path,))
                        if native_cached is not None:
                            dotbracket_native = native_cached
                        else:
                            stackings_native, pairings_native, res_native = bb.annotate(native_pdb_path)
                            dotbracket_native = bb.dot_bracket(pairings_native, res_native)[0]
                            computation_cache.set(session_id, 'native_dotbracket', (native_pdb_path,), dotbracket_native)
                        
                        fig = plot_arc(dotbracket_data, dotbracket_native)
                        fig.write_html(os.path.join(plot_path, "arc_plot.html"))
                        results['ARC'] = plotly_to_json(fig)
                
                logger.info(f"Ultra-optimized {plot_type} completed in {time.time() - plot_start:.2f}s")
        
        # Phase 3: Handle remaining plots with standard optimization
        remaining_plots = [plot for plot in selected_plots if plot not in results]
        
        if remaining_plots:
            logger.info(f"Phase 3: Processing remaining plots with standard optimization: {remaining_plots}")
            
            for plot_type in remaining_plots:
                if plot_type == 'CONTACT_MAPS':
                    # Contact maps still require separate processing
                    from create_plots import create_contact_maps
                    create_contact_maps(traj_xtc_path, native_pdb_path, 1, session_id)
                    path_save_figure = os.path.join('static', session_id, 'contact_map_plotly.png')
                    results['CONTACT_MAPS'] = json.dumps({'path': path_save_figure})
        
        # Phase 4: Cache optimization results
        logger.info("Phase 4: Caching optimization results")
        optimization_cache_key = f"ultra_optimization_{session_id}"
        computation_cache.set(session_id, optimization_cache_key, 
                            (native_pdb_path, traj_xtc_path), 
                            processor.get_optimization_stats())
        
        total_time = time.time() - total_start_time
        avg_time_per_plot = total_time / len(selected_plots) if selected_plots else 0
        
        logger.info(f"ULTRA-OPTIMIZED analysis completed:")
        logger.info(f"  Total time: {total_time:.2f}s")
        logger.info(f"  Plots processed: {len(selected_plots)}")
        logger.info(f"  Average time per plot: {avg_time_per_plot:.2f}s")
        logger.info(f"  Optimization level: {optimization_level}")
        
        # Add performance metadata to results (ensure JSON serializable)
        optimization_stats = processor.get_optimization_stats()
        serializable_stats = {
            k: v for k, v in optimization_stats.items() 
            if isinstance(v, (int, float, str, bool, list, dict, type(None)))
        }
        
        results['_performance_metadata'] = {
            'total_time': total_time,
            'avg_time_per_plot': avg_time_per_plot,
            'optimization_level': optimization_level,
            'plots_processed': len(selected_plots),
            'optimization_stats': serializable_stats
        }
        
        return results
        
    except Exception as exc:
        logger.error(f"Ultra-optimized RNA analysis failed: {str(exc)}")
        logger.error(f"Traceback: {traceback.format_exc()}")
        self.retry(exc=exc, countdown=60)

@app3.task(bind=True, max_retries=3)
@ultra_log_task_execution
def parallel_rmsd_analysis(self, native_pdb_path, traj_xtc_path, download_path, plot_path, 
                          session_id, reference_frame=0, selection='protein'):
    """
    Parallel RMSD calculation using MDAnalysis 2024 techniques
    """
    try:
        config = OptimizationConfig(use_parallel_processing=True, n_workers=4)
        processor = OptimizedTrajectoryProcessor(native_pdb_path, traj_xtc_path, config)
        
        # Use optimized parallel RMSD calculation
        rmsd_values = processor.efficient_rmsd_calculation(reference_frame, selection)
        
        # Create plot
        fig = plot_rmsd(rmsd_values)
        fig.write_html(os.path.join(plot_path, "rmsd_plot.html"))
        
        # Save data
        rmsd_df = pd.DataFrame({"RMSD": rmsd_values})
        rmsd_df.to_csv(os.path.join(download_path, "rmsd_values.csv"), index=False)
        
        return plotly_to_json(fig)
        
    except Exception as exc:
        logger.error(f"Parallel RMSD analysis failed: {str(exc)}")
        self.retry(exc=exc, countdown=60)

@app3.task(bind=True, max_retries=3)
@ultra_log_task_execution
def batch_properties_analysis_task(self, native_pdb_path, traj_xtc_path, download_path, plot_path, 
                                  session_id, analysis_type='batch_properties'):
    """
    Batch property analysis for multiple calculations in single trajectory pass
    """
    try:
        config = OptimizationConfig(preload_selections=True, enable_caching=True)
        processor = OptimizedTrajectoryProcessor(native_pdb_path, traj_xtc_path, config)
        
        if analysis_type == 'batch_properties':
            # Define multiple properties to calculate in single pass
            properties = {
                'radius_of_gyration': lambda u: u.select_atoms('protein').radius_of_gyration() if u.select_atoms('protein') else 0,
                'center_of_mass': lambda u: u.select_atoms('protein').center_of_mass()[2] if u.select_atoms('protein') else 0,  # Z-coordinate
            }
            
            results = processor.batch_property_calculation(properties)
            
            # Create plots for each property
            plot_results = {}
            for prop_name, values in results.items():
                import plotly.graph_objects as go
                fig = go.Figure()
                fig.add_trace(go.Scattergl(y=values, mode='markers', name=prop_name))
                fig.update_layout(title=f"{prop_name} vs Frame", xaxis_title="Frame", yaxis_title=prop_name)
                
                # Save plot
                fig.write_html(os.path.join(plot_path, f"{prop_name}_plot.html"))
                plot_results[prop_name] = plotly_to_json(fig)
                
                # Save data
                df = pd.DataFrame({prop_name: values})
                df.to_csv(os.path.join(download_path, f"{prop_name}_values.csv"), index=False)
            
            return plot_results
        
    except Exception as exc:
        logger.error(f"Batch properties analysis failed: {str(exc)}")
        self.retry(exc=exc, countdown=60)

# Performance monitoring task
# Individual cached task implementations for backward compatibility

@app3.task(bind=True, max_retries=3)
@ultra_log_task_execution
def generate_contact_maps_plot(self, native_pdb, traj_xtc, download_path, plot_path, session_id):
    try:
        from create_plots import create_contact_maps
        create_contact_maps(traj_xtc, native_pdb, 1, session_id)
        path_save_figure = os.path.join('static', session_id, 'contact_map_plotly.png')
        return json.dumps({'path': path_save_figure})
    except Exception as exc:
        logger.error(f"Contact maps generation failed: {str(exc)}")
        self.retry(exc=exc, countdown=60)

@app3.task(bind=True, max_retries=3)
@ultra_log_task_execution
def generate_rmsd_plot(self, native_pdb_path, traj_xtc_path, download_path, plot_path, session_id):
    try:
        # Check cache first
        cached_rmsd = computation_cache.get(session_id, 'bb_rmsd', (native_pdb_path, traj_xtc_path))
        
        if cached_rmsd is not None:
            logger.info("Using cached RMSD data")
            rmsd = cached_rmsd
        else:
            logger.info("Computing RMSD data")
            rmsd = bb.rmsd(
                native_pdb_path,
                traj_xtc_path,
                topology=native_pdb_path,
                heavy_atom=True,
            )
            # Cache the result
            computation_cache.set(session_id, 'bb_rmsd', (native_pdb_path, traj_xtc_path), rmsd)

        # Create plot
        fig = plot_rmsd(rmsd)
        
        # Save data and plot
        rmsd_df = pd.DataFrame({"RMSD": rmsd})
        rmsd_df.to_csv(os.path.join(download_path, "rmsd_values.csv"), index=False)
        fig.write_html(os.path.join(plot_path, "rmsd_plot.html"))

        return plotly_to_json(fig)
    except Exception as exc:
        logger.error(f"RMSD plot generation failed: {str(exc)}")
        self.retry(exc=exc, countdown=60)

@app3.task(bind=True, max_retries=3)
@ultra_log_task_execution
def generate_ermsd_plot(self, native_pdb_path, traj_xtc_path, download_path, plot_path, session_id):
    try:
        # Check cache first
        cached_ermsd = computation_cache.get(session_id, 'bb_ermsd', (native_pdb_path, traj_xtc_path))
        
        if cached_ermsd is not None:
            logger.info("Using cached ERMSD data")
            ermsd = cached_ermsd
        else:
            logger.info("Computing ERMSD data")
            ermsd = bb.ermsd(native_pdb_path, traj_xtc_path, topology=native_pdb_path)
            # Cache the result
            computation_cache.set(session_id, 'bb_ermsd', (native_pdb_path, traj_xtc_path), ermsd)

        # Create plot
        fig = plot_ermsd(ermsd)
        
        # Save data and plot
        ermsd_df = pd.DataFrame({"ERMSD": ermsd})
        ermsd_df.to_csv(os.path.join(download_path, "ermsd_values.csv"), index=False)
        fig.write_html(os.path.join(plot_path, "ermsd_plot.html"))

        return plotly_to_json(fig)
    except Exception as exc:
        logger.error(f"ERMSD plot generation failed: {str(exc)}")
        self.retry(exc=exc, countdown=60)

@app3.task(bind=True, max_retries=3)
@ultra_log_task_execution
def generate_torsion_plot(self, native_pdb_path, traj_xtc_path, download_path, plot_path, session_id, torsion_params=None):
    try:
        # Check cache first
        cached_angles = computation_cache.get(session_id, 'bb_backbone_angles', (native_pdb_path, traj_xtc_path))
        
        if cached_angles is not None:
            logger.info("Using cached torsion angles data")
            angles, res = cached_angles
        else:
            logger.info("Computing torsion angles data")
            angles, res = bb.backbone_angles(traj_xtc_path, topology=native_pdb_path)
            # Cache the result
            computation_cache.set(session_id, 'bb_backbone_angles', (native_pdb_path, traj_xtc_path), (angles, res))
        
        logger.info(f"Calculated torsion angles for {len(res)} residues")
        
        # Handle backward compatibility
        if isinstance(torsion_params, dict):
            fig = plot_torsion_enhanced(angles, res, torsion_params)
        else:
            # Old single residue format
            torsion_residue = torsion_params if torsion_params is not None else 0
            fig = plot_torsion(angles, res, torsion_residue)
            
        fig.write_html(os.path.join(plot_path, "torsion_plot.html"))
        
        # Save data
        angles_df = pd.DataFrame(angles.reshape(-1, angles.shape[-1]), 
                               columns=["alpha", "beta", "gamma", "delta", "epsilon", "zeta", "chi"])
        angles_df.to_csv(os.path.join(download_path, "torsion_angles.csv"), index=False)
        
        return plotly_to_json(fig)
    except Exception as exc:
        logger.error(f"Torsion calculation failed: {str(exc)}")
        if self.request.retries < self.max_retries:
            logger.info(f"Retrying torsion calculation (attempt {self.request.retries + 1}/{self.max_retries})")
            raise self.retry(exc=exc, countdown=60)
        else:
            logger.error(f"Torsion calculation failed permanently after {self.max_retries} retries")
            raise exc

@app3.task(bind=True, max_retries=3)
@ultra_log_task_execution
def generate_sec_structure_plot(self, native_pdb_path, traj_xtc_path, download_path, plot_path, session_id):
    try:
        # Check cache first
        cached_annotate = computation_cache.get(session_id, 'bb_annotate', (native_pdb_path, traj_xtc_path))
        
        if cached_annotate is not None:
            logger.info("Using cached annotate data")
            stackings, pairings, res = cached_annotate
        else:
            logger.info("Computing annotate data")
            stackings, pairings, res = bb.annotate(traj_xtc_path, topology=native_pdb_path)
            # Cache the result
            computation_cache.set(session_id, 'bb_annotate', (native_pdb_path, traj_xtc_path), (stackings, pairings, res))
        
        dotbracket_data, res2 = bb.dot_bracket(pairings, res)
        return [dotbracket_data, res2.strip()]
    except Exception as exc:
        logger.error(f"Secondary structure generation failed: {str(exc)}")
        self.retry(exc=exc, countdown=60)

@app3.task(bind=True, max_retries=3)
@ultra_log_task_execution
def generate_2Dpairing_plot(self, native_pdb_path, traj_xtc_path, download_path, plot_path, session_id):
    try:
        # Check cache first
        cache_key = 'bb_dump_rvec'
        cached_rvecs = computation_cache.get(session_id, cache_key, (native_pdb_path, traj_xtc_path))
        
        if cached_rvecs is not None:
            logger.info("Using cached rvec data")
            rvecs_traj, res_traj = cached_rvecs
        else:
            logger.info("Computing rvec data")
            rvecs_traj, res_traj = bb.dump_rvec(traj_xtc_path, topology=native_pdb_path, cutoff=100.0)
            # Cache the result
            computation_cache.set(session_id, cache_key, (native_pdb_path, traj_xtc_path), (rvecs_traj, res_traj))
        
        nonzero = np.where(np.sum(rvecs_traj**2, axis=3) > 0.01)
        rr = rvecs_traj[nonzero]
        fig = base_pairs_visualisation(rr)
        fig.write_html(os.path.join(plot_path, "2D_pairing.html"))
        return plotly_to_json(fig)
    except Exception as exc:
        logger.error(f"2D pairing plot generation failed: {str(exc)}")
        self.retry(exc=exc, countdown=60)

@app3.task(bind=True, max_retries=3)
@ultra_log_task_execution
def generate_annotate_plot(self, traj_xtc_path, native_pdb_path, download_path, plot_path, session_id):
    try:
        # Check cache first
        cached_annotate = computation_cache.get(session_id, 'bb_annotate', (native_pdb_path, traj_xtc_path))
        
        if cached_annotate is not None:
            logger.info("Using cached annotate data")
            stackings, pairings, res = cached_annotate
        else:
            logger.info("Computing annotate data")
            stackings, pairings, res = bb.annotate(traj_xtc_path, topology=native_pdb_path)
            # Cache the result
            computation_cache.set(session_id, 'bb_annotate', (native_pdb_path, traj_xtc_path), (stackings, pairings, res))
        
        return ["plot", "annotate", stackings, pairings, res]
    except Exception as exc:
        logger.error(f"Annotate plot generation failed: {str(exc)}")
        self.retry(exc=exc, countdown=60)

@app3.task(bind=True, max_retries=3)
@ultra_log_task_execution
def generate_dotbracket_plot(self, native_pdb_path, traj_xtc_path, download_path, plot_path, session_id):
    try:
        # Check cache first
        cached_dotbracket = computation_cache.get(session_id, 'dotbracket', (native_pdb_path, traj_xtc_path))
        
        if cached_dotbracket is not None:
            logger.info("Using cached dotbracket data")
            dotbracket_data = cached_dotbracket
        else:
            # Try to get annotate data from cache
            cached_annotate = computation_cache.get(session_id, 'bb_annotate', (native_pdb_path, traj_xtc_path))
            
            if cached_annotate is not None:
                logger.info("Using cached annotate data for dotbracket")
                stackings, pairings, res = cached_annotate
            else:
                logger.info("Computing annotate data for dotbracket")
                stackings, pairings, res = bb.annotate(traj_xtc_path, topology=native_pdb_path)
                # Cache the annotate result
                computation_cache.set(session_id, 'bb_annotate', (native_pdb_path, traj_xtc_path), (stackings, pairings, res))
            
            dotbracket_data = bb.dot_bracket(pairings, res)[0]
            # Cache the dotbracket result
            computation_cache.set(session_id, 'dotbracket', (native_pdb_path, traj_xtc_path), dotbracket_data)
        
        fig = plot_dotbracket(dotbracket_data)
        fig.write_html(os.path.join(plot_path, "dotbracket_timeline_plot.html"))
        return plotly_to_json(fig)
    except Exception as exc:
        logger.error(f"Dotbracket plot generation failed: {str(exc)}")
        self.retry(exc=exc, countdown=60)

@app3.task(bind=True, max_retries=3)
@ultra_log_task_execution
def generate_arc_plot(self, native_pdb_path, traj_xtc_path, download_path, plot_path, session_id):
    try:
        # Check cache for trajectory dotbracket
        cached_dotbracket = computation_cache.get(session_id, 'dotbracket', (native_pdb_path, traj_xtc_path))
        
        if cached_dotbracket is not None:
            logger.info("Using cached trajectory dotbracket data")
            dotbracket_data = cached_dotbracket
            # We still need res for sequence info, so get annotate data
            cached_annotate = computation_cache.get(session_id, 'bb_annotate', (native_pdb_path, traj_xtc_path))
            if cached_annotate is not None:
                stackings, pairings, res = cached_annotate
            else:
                logger.info("Computing annotate data for sequence info")
                stackings, pairings, res = bb.annotate(traj_xtc_path, topology=native_pdb_path)
                computation_cache.set(session_id, 'bb_annotate', (native_pdb_path, traj_xtc_path), (stackings, pairings, res))
        else:
            # Need to compute both
            cached_annotate = computation_cache.get(session_id, 'bb_annotate', (native_pdb_path, traj_xtc_path))
            if cached_annotate is not None:
                logger.info("Using cached annotate data")
                stackings, pairings, res = cached_annotate
            else:
                logger.info("Computing trajectory annotate data")
                stackings, pairings, res = bb.annotate(traj_xtc_path, topology=native_pdb_path)
                computation_cache.set(session_id, 'bb_annotate', (native_pdb_path, traj_xtc_path), (stackings, pairings, res))
            
            dotbracket_data = bb.dot_bracket(pairings, res)[0]
            computation_cache.set(session_id, 'dotbracket', (native_pdb_path, traj_xtc_path), dotbracket_data)
        
        # Get native dotbracket (cached separately)
        cached_native_dotbracket = computation_cache.get(session_id, 'native_dotbracket', (native_pdb_path,))
        
        if cached_native_dotbracket is not None:
            logger.info("Using cached native dotbracket data")
            dotbracket_native = cached_native_dotbracket
        else:
            logger.info("Computing native dotbracket data")
            stackings_native, pairings_native, res_native = bb.annotate(native_pdb_path)
            dotbracket_native = bb.dot_bracket(pairings_native, res_native)[0]
            computation_cache.set(session_id, 'native_dotbracket', (native_pdb_path,), dotbracket_native)
        
        sequence = ''.join([item[0] for item in res])
        resids = [item.split("_")[1] for item in res]
        
        dotbracket_df = pd.DataFrame(dotbracket_data, columns=["DotBracket"])
        dotbracket_df.to_csv(os.path.join(download_path, "dotbracket_data.csv"), index=False)
        
        fig = plot_diagram_frequency(sequence, dotbracket_data, dotbracket_native)
        fig.write_html(os.path.join(plot_path, "arc_diagram_plot.html"))
        
        return [plotly_to_json(fig), resids]
    except Exception as exc:
        logger.error(f"Arc plot generation failed: {str(exc)}")
        self.retry(exc=exc, countdown=60)

@app3.task(bind=True, max_retries=3)
@ultra_log_task_execution
def generate_ds_motif_plot(self, traj_xtc_path, native_pdb_path, download_path, plot_path, session_id):
    # Placeholder implementation - to be completed based on requirements
    logger.info("DS_MOTIF plot generation - placeholder implementation")
    return {"status": "not_implemented"}

@app3.task(bind=True, max_retries=3)
@ultra_log_task_execution
def generate_ss_motif_plot(self, traj_xtc_path, native_pdb_path, download_path, plot_path, session_id):
    # Placeholder implementation - to be completed based on requirements
    logger.info("SS_MOTIF plot generation - placeholder implementation")
    return {"status": "not_implemented"}

@app3.task(bind=True, max_retries=3)
@ultra_log_task_execution
def generate_jcoupling_plot(self, traj_xtc_path, native_pdb_path, download_path, plot_path, session_id):
    try:
        # Check cache first
        cached_jcouplings = computation_cache.get(session_id, 'bb_jcouplings', (native_pdb_path, traj_xtc_path))
        
        if cached_jcouplings is not None:
            logger.info("Using cached J-coupling data")
            couplings, res = cached_jcouplings
        else:
            logger.info("Computing J-coupling data")
            couplings, res = bb.jcouplings(traj_xtc_path, topology=native_pdb_path)
            # Cache the result
            computation_cache.set(session_id, 'bb_jcouplings', (native_pdb_path, traj_xtc_path), (couplings, res))
        
        return ["plot", "jcoupling", couplings]
    except Exception as exc:
        logger.error(f"J-coupling plot generation failed: {str(exc)}")
        self.retry(exc=exc, countdown=60)

@app3.task(bind=True, max_retries=3)
@ultra_log_task_execution
def generate_escore_plot(self, traj_xtc_path, native_pdb_path, download_path, plot_path, session_id):
    # Placeholder implementation - to be completed based on requirements
    logger.info("ESCORE plot generation - placeholder implementation")
    return {"status": "not_implemented"}

@app3.task(bind=True, max_retries=3)
@ultra_log_task_execution
def update_contact_map_plot(self, generate_data_path, plot_path, frame_number, session_id):
    try:
        import pickle
        with open(os.path.join(generate_data_path, "contact_map_data.pkl"), 'rb') as f:
            loaded_data = pickle.load(f)

        frames_dict = loaded_data['frames_dict']
        sequence = loaded_data['sequence']

        frame_num, frame_data = sorted(frames_dict.items())[frame_number]
        logger.info(f"Updating contact map for frame {frame_num}")
        
        fig = plot_rna_contact_map(frame_data, sequence, output_file=os.path.join(plot_path, "contact_map.html"), frame_number=frame_num)
        return plotly_to_json(fig)
    except Exception as exc:
        logger.error(f"Contact map update failed: {str(exc)}")
        self.retry(exc=exc, countdown=60)

@app3.task(bind=True, max_retries=3)
@ultra_log_task_execution
def update_landscape_frame(self, generate_data_path, coordinates):
    try:
        import pickle
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
    except Exception as exc:
        logger.error(f"Landscape frame update failed: {str(exc)}")
        self.retry(exc=exc, countdown=60)

@app3.task(bind=True, max_retries=3)
@ultra_log_task_execution
def generate_contact_map_plot(self, native_pdb_path, traj_xtc_path, download_path, plot_path, generate_data_path, session_id):
    try:
        def process_barnaba_pairings(pairings, res):
            """Process barnaba pairings and residue information into frames dictionary"""
            frames_dict = {}
            sequence = [r[0] for r in res]
            
            for frame_num, frame_data in enumerate(pairings):
                base_pairs = []
                
                if len(frame_data) == 2:
                    pair_indices = frame_data[0]
                    annotations = frame_data[1]
                    
                    for pair_idx, pair in enumerate(pair_indices):
                        if not pair:
                            continue
                        
                        res_i = pair[0] + 1
                        res_j = pair[1] + 1
                        
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
                
                frames_dict[frame_num] = pd.DataFrame(base_pairs)
            
            return frames_dict, sequence
        
        # Check cache first
        cached_annotate = computation_cache.get(session_id, 'bb_annotate', (native_pdb_path, traj_xtc_path))
        
        if cached_annotate is not None:
            logger.info("Using cached annotate data for contact map")
            stackings, pairings, res = cached_annotate
        else:
            logger.info("Computing annotate data for contact map")
            stackings, pairings, res = bb.annotate(traj_xtc_path, topology=native_pdb_path)
            # Cache the result
            computation_cache.set(session_id, 'bb_annotate', (native_pdb_path, traj_xtc_path), (stackings, pairings, res))
        
        # Process results
        frames_dict, sequence = process_barnaba_pairings(pairings, res)
        logger.info(f"RNA length: {len(sequence)} nucleotides")
        logger.info(f"Found {len(frames_dict)} frames in the data")
        
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
        return plotly_to_json(fig)
    except Exception as exc:
        logger.error(f"Contact map generation failed: {str(exc)}")
        self.retry(exc=exc, countdown=60)

@app3.task(bind=True, max_retries=3)
@ultra_log_task_execution
def generate_landscape_plot(self, native_pdb_path, traj_xtc_path, download_path, plot_path, session_id, parameters, generate_data_path):
    try:
        import mdtraj as md
        import energy_3dplot
        
        def calculate_structural_metrics(parameters, traj_load, native_load, native_pdb_path, traj_xtc_path):
            """Calculate structural metrics for trajectory analysis based on specified parameters."""
            calculation_functions = {
                "Fraction of Contact Formed": lambda: best_hummer_q(traj_load, native_load),
                "RMSD": lambda: (
                    computation_cache.get(session_id, 'bb_rmsd', (native_pdb_path, traj_xtc_path))
                    if computation_cache.get(session_id, 'bb_rmsd', (native_pdb_path, traj_xtc_path)) is not None else
                    bb.rmsd(native_pdb_path, traj_xtc_path, topology=native_pdb_path, heavy_atom=True)
                ),
                "eRMSD": lambda: (
                    computation_cache.get(session_id, 'bb_ermsd', (native_pdb_path, traj_xtc_path))
                    if computation_cache.get(session_id, 'bb_ermsd', (native_pdb_path, traj_xtc_path)) is not None else
                    bb.ermsd(native_pdb_path, traj_xtc_path, topology=native_pdb_path)
                ),
                "Torsion": lambda: logger.info('torsion')
            }

            results = {}
            for metric in parameters:
                if metric in calculation_functions:
                    result = calculation_functions[metric]()
                    results[metric] = result
                    # Cache RMSD/eRMSD results if they weren't cached
                    if metric == "RMSD" and not computation_cache.get(session_id, 'bb_rmsd', (native_pdb_path, traj_xtc_path)):
                        computation_cache.set(session_id, 'bb_rmsd', (native_pdb_path, traj_xtc_path), result)
                    elif metric == "eRMSD" and not computation_cache.get(session_id, 'bb_ermsd', (native_pdb_path, traj_xtc_path)):
                        computation_cache.set(session_id, 'bb_ermsd', (native_pdb_path, traj_xtc_path), result)
                else:
                    results[metric] = None
                    logger.warning(f"Unknown metric '{metric}'. Skipping calculation.")
            return results
        
        traj_load = md.load_xtc(traj_xtc_path, native_pdb_path)
        native_load = md.load(native_pdb_path)
        metrics_to_calculate = [parameters[1], parameters[2]]

        results = calculate_structural_metrics(
            metrics_to_calculate,
            traj_load,
            native_load,
            native_pdb_path,
            traj_xtc_path
        )
        
        stride = int(parameters[0])
        component1 = [results[parameters[1]][x] for x in range(0, len(traj_load), stride)]
        component2 = [results[parameters[2]][x] for x in range(0, len(traj_load), stride)]

        df = pd.DataFrame({
            "frame": list(range(0, len(component2))),
            "Q": component1,
            "RMSD": component2,
            "traj": "traj_1",
        })
        
        with open(os.path.join(generate_data_path, "dataframe.pkl"), 'wb') as f:
            pickle.dump(df, f)
        
        logger.info("Finished dataframe creation")
        
        size = 65
        selected_regions = []
        dataframe, max_RMSD, max_Q = df, max(df["RMSD"]), max(df['Q'])
        
        (
            probability_matrix,
            allframes_matrix,
            Qbin,
            RMSDbin,
        ) = energy_3dplot.make_matrix_probability(dataframe, size, max_RMSD, max_Q)
        
        energy_matrix, real_values = energy_3dplot.make_matrix_energy(
            probability_matrix, max_RMSD, size
        )
        
        fig = plot_landscapes_3D(
            energy_matrix, Qbin, RMSDbin, max_RMSD, max_Q, real_values, selected_regions, metrics_to_calculate
        )
        fig2 = plot_landscapes_2D(
            energy_matrix, Qbin, RMSDbin, max_RMSD, max_Q, real_values, selected_regions, metrics_to_calculate
        )
        
        path_landscape_3d = f"static/{session_id}/download_plot/LANDSCAPE/landscape.html"
        fig.write_html(path_landscape_3d)
        
        return [plotly_to_json(fig), plotly_to_json(fig2)]
    except Exception as exc:
        logger.error(f"Landscape plot generation failed: {str(exc)}")
        self.retry(exc=exc, countdown=60)

@app3.task(bind=True, max_retries=3)
@ultra_log_task_execution
def batch_computation_coordinator(self, native_pdb_path, traj_xtc_path, download_path, plot_path, session_id, selected_plots):
    """
    Ultra-optimized batch computation coordinator
    Uses the consolidated ultra_optimized_rna_analysis for maximum efficiency
    """
    try:
        logger.info(f"Starting ultra-optimized batch coordinator for {len(selected_plots)} plots")
        
        # Use the ultra-optimized consolidated analysis
        results = ultra_optimized_rna_analysis(
            native_pdb_path, traj_xtc_path, download_path, plot_path, session_id, selected_plots
        )
        
        logger.info(f"Batch coordination completed with {len(results)} results")
        return results
        
    except Exception as exc:
        logger.error(f"Batch coordination failed: {str(exc)}")
        self.retry(exc=exc, countdown=60)

@app3.task(bind=True)
def performance_benchmark(self, native_pdb_path, traj_xtc_path, session_id):
    """
    Benchmark different optimization strategies
    """
    benchmark_results = {}
    
    # Test different configurations
    configs = {
        'basic': OptimizationConfig(use_parallel_processing=False, n_workers=1),
        'parallel': OptimizationConfig(use_parallel_processing=True, n_workers=4),
        'maximum': OptimizationConfig(
            use_parallel_processing=True, n_workers=8, 
            preload_selections=True, enable_caching=True
        )
    }
    
    for config_name, config in configs.items():
        start_time = time.time()
        
        try:
            processor = OptimizedTrajectoryProcessor(native_pdb_path, traj_xtc_path, config)
            
            # Run a simple RMSD calculation as benchmark
            rmsd_values = processor.efficient_rmsd_calculation()
            
            benchmark_results[config_name] = {
                'time': time.time() - start_time,
                'success': True,
                'n_frames': len(rmsd_values)
            }
            
        except Exception as e:
            benchmark_results[config_name] = {
                'time': time.time() - start_time,
                'success': False,
                'error': str(e)
            }
    
    return benchmark_results