"""
Optimized Celery tasks with computation caching and batch processing
Reduces redundant file I/O and computation time
"""

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
import logging
from functools import wraps
from typing import Optional, Dict, Any, List
import traceback
from celery.exceptions import Retry

# Import our computation cache system
from computation_cache import computation_cache, cached_computation, batch_manager

# Configure logging for Celery tasks
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Configure Celery
app2 = Celery('tasks')

def log_task_execution(func):
    """Decorator to log task execution and handle errors"""
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        task_name = func.__name__
        session_id = args[-1] if args else 'unknown'
        
        logger.info(f"Starting optimized task {task_name} for session {session_id}")
        start_time = time.time()
        
        try:
            result = func(self, *args, **kwargs)
            execution_time = time.time() - start_time
            logger.info(f"Completed optimized task {task_name} for session {session_id} in {execution_time:.2f}s")
            return result
        except Exception as exc:
            execution_time = time.time() - start_time
            logger.error(f"Optimized task {task_name} failed for session {session_id} after {execution_time:.2f}s: {str(exc)}")
            logger.error(f"Traceback: {traceback.format_exc()}")
            raise
    return wrapper

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
    return to_json(fig, validate=False, engine="orjson")

@app2.task(bind=True, max_retries=3)
@log_task_execution
def generate_rmsd_plot_optimized(self, native_pdb_path, traj_xtc_path, download_path, plot_path, session_id):
    """Optimized RMSD plot generation with caching"""
    try:
        start_time = time.time()

        # Try to get cached result
        cached_rmsd = computation_cache.get(session_id, 'bb_rmsd', (native_pdb_path, traj_xtc_path))
        
        if cached_rmsd is not None:
            logger.info("Using cached RMSD data")
            rmsd = cached_rmsd
        else:
            # Calculate RMSD
            rmsd_start = time.time()
            rmsd = bb.rmsd(
                native_pdb_path,
                traj_xtc_path,
                topology=native_pdb_path,
                heavy_atom=True,
            )
            logger.info(f"RMSD calculation time: {time.time() - rmsd_start:.2f}s")
            
            # Cache the result
            computation_cache.set(session_id, 'bb_rmsd', (native_pdb_path, traj_xtc_path), rmsd)

        # Create plot
        plot_start = time.time()
        fig = plot_rmsd(rmsd)
        logger.info(f"Plot creation time: {time.time() - plot_start:.2f}s")

        # Save data and plot
        save_start = time.time()
        rmsd_df = pd.DataFrame({"RMSD": rmsd})
        rmsd_df.to_csv(os.path.join(download_path, "rmsd_values.csv"), index=False)
        fig.write_html(os.path.join(plot_path, "rmsd_plot.html"))
        logger.info(f"Save time: {time.time() - save_start:.2f}s")

        # Convert to JSON
        json_start = time.time()
        plotly_data = plotly_to_json(fig)
        logger.info(f"JSON conversion time: {time.time() - json_start:.2f}s")

        logger.info(f"Total optimized RMSD plot generation time: {time.time() - start_time:.2f}s")
        return plotly_data
    except Exception as exc:
        self.retry(exc=exc, countdown=60)

@app2.task(bind=True, max_retries=3)
@log_task_execution
def generate_ermsd_plot_optimized(self, native_pdb_path, traj_xtc_path, download_path, plot_path, session_id):
    """Optimized ERMSD plot generation with caching"""
    try:
        start_time = time.time()

        # Try to get cached result
        cached_ermsd = computation_cache.get(session_id, 'bb_ermsd', (native_pdb_path, traj_xtc_path))
        
        if cached_ermsd is not None:
            logger.info("Using cached ERMSD data")
            ermsd = cached_ermsd
        else:
            # Calculate ERMSD
            ermsd_start = time.time()
            ermsd = bb.ermsd(native_pdb_path, traj_xtc_path, topology=native_pdb_path)
            logger.info(f"ERMSD calculation time: {time.time() - ermsd_start:.2f}s")
            
            # Cache the result
            computation_cache.set(session_id, 'bb_ermsd', (native_pdb_path, traj_xtc_path), ermsd)

        # Create plot
        plot_start = time.time()
        fig = plot_ermsd(ermsd)
        logger.info(f"Plot creation time: {time.time() - plot_start:.2f}s")

        # Save data and plot
        save_start = time.time()
        ermsd_df = pd.DataFrame({"ERMSD": ermsd})
        ermsd_df.to_csv(os.path.join(download_path, "ermsd_values.csv"), index=False)
        fig.write_html(os.path.join(plot_path, "ermsd_plot.html"))
        logger.info(f"Save time: {time.time() - save_start:.2f}s")

        # Convert to JSON
        json_start = time.time()
        plotly_data = plotly_to_json(fig)
        logger.info(f"JSON conversion time: {time.time() - json_start:.2f}s")

        logger.info(f"Total optimized ERMSD plot generation time: {time.time() - start_time:.2f}s")
        return plotly_data
    except Exception as exc:
        self.retry(exc=exc, countdown=60)

@app2.task(bind=True, max_retries=3)
@log_task_execution
def generate_torsion_plot_optimized(self, native_pdb_path, traj_xtc_path, download_path, plot_path, session_id, torsion_params=None):
    """Optimized torsion plot generation with caching"""
    try:
        start_time = time.time()
        
        # Try to get cached result
        cached_angles = computation_cache.get(session_id, 'bb_backbone_angles', (native_pdb_path, traj_xtc_path))
        
        if cached_angles is not None:
            logger.info("Using cached backbone angles data")
            angles, res = cached_angles
        else:
            # Calculate angles
            angles_start = time.time()
            angles, res = bb.backbone_angles(traj_xtc_path, topology=native_pdb_path)
            logger.info(f"Backbone angles calculation time: {time.time() - angles_start:.2f}s")
            logger.info(f"Calculated torsion angles for {len(res)} residues")
            
            # Cache the result
            computation_cache.set(session_id, 'bb_backbone_angles', (native_pdb_path, traj_xtc_path), (angles, res))
        
        # Handle backward compatibility
        plot_start = time.time()
        if isinstance(torsion_params, dict):
            fig = plot_torsion_enhanced(angles, res, torsion_params)
        else:
            # Old single residue format
            torsion_residue = torsion_params if torsion_params is not None else 0
            fig = plot_torsion(angles, res, torsion_residue)
        logger.info(f"Plot creation time: {time.time() - plot_start:.2f}s")
            
        fig.write_html(os.path.join(plot_path, "torsion_plot.html"))
        
        # Save data
        save_start = time.time()
        angles_df = pd.DataFrame(angles.reshape(-1, angles.shape[-1]), 
                               columns=["alpha", "beta", "gamma", "delta", "epsilon", "zeta", "chi"])
        angles_df.to_csv(os.path.join(download_path, "torsion_angles.csv"), index=False)
        logger.info(f"Save time: {time.time() - save_start:.2f}s")
        
        plotly_data = plotly_to_json(fig)
        logger.info(f"Total optimized torsion plot generation time: {time.time() - start_time:.2f}s")
        return plotly_data
    except Exception as exc:
        logger.error(f"Torsion calculation failed: {str(exc)}")
        if self.request.retries < self.max_retries:
            logger.info(f"Retrying torsion calculation (attempt {self.request.retries + 1}/{self.max_retries})")
            raise self.retry(exc=exc, countdown=60)
        else:
            logger.error(f"Torsion calculation failed permanently after {self.max_retries} retries")
            raise exc

@app2.task(bind=True, max_retries=3)
@log_task_execution
def generate_annotate_plots_batch(self, native_pdb_path, traj_xtc_path, download_path, plot_path, session_id, plot_types):
    """
    Batch processing for plots that use bb.annotate():
    - SEC_STRUCTURE
    - DOTBRACKET  
    - ARC
    - CONTACT_MAPS
    """
    try:
        start_time = time.time()
        results = {}
        
        # Try to get cached annotate result
        cached_annotate = computation_cache.get(session_id, 'bb_annotate', (native_pdb_path, traj_xtc_path))
        
        if cached_annotate is not None:
            logger.info("Using cached bb.annotate data for batch processing")
            stackings, pairings, res = cached_annotate
        else:
            # Calculate annotate once for all plots
            annotate_start = time.time()
            stackings, pairings, res = bb.annotate(traj_xtc_path, topology=native_pdb_path)
            logger.info(f"bb.annotate calculation time: {time.time() - annotate_start:.2f}s")
            
            # Cache the result
            computation_cache.set(session_id, 'bb_annotate', (native_pdb_path, traj_xtc_path), (stackings, pairings, res))
        
        # Get or compute dotbracket data if needed
        dotbracket_data = None
        if any(plot_type in ['DOTBRACKET', 'ARC'] for plot_type in plot_types):
            cached_dotbracket = computation_cache.get(session_id, 'dotbracket_data', (native_pdb_path, traj_xtc_path))
            
            if cached_dotbracket is not None:
                logger.info("Using cached dotbracket data")
                dotbracket_data = cached_dotbracket
            else:
                dotbracket_start = time.time()
                dotbracket_data = bb.dot_bracket(pairings, res)[0]
                logger.info(f"Dotbracket calculation time: {time.time() - dotbracket_start:.2f}s")
                
                # Cache the result
                computation_cache.set(session_id, 'dotbracket_data', (native_pdb_path, traj_xtc_path), dotbracket_data)
        
        # Process each requested plot type
        for plot_type in plot_types:
            plot_start = time.time()
            
            if plot_type == 'SEC_STRUCTURE':
                results['SEC_STRUCTURE'] = [dotbracket_data if dotbracket_data else bb.dot_bracket(pairings, res)[0], res.strip()]
            
            elif plot_type == 'DOTBRACKET':
                fig = plot_dotbracket(dotbracket_data)
                fig.write_html(os.path.join(plot_path, "dotbracket_timeline_plot.html"))
                results['DOTBRACKET'] = plotly_to_json(fig)
            
            elif plot_type == 'ARC':
                # Get native dotbracket for ARC plot
                native_cached = computation_cache.get(session_id, 'native_dotbracket', (native_pdb_path,))
                if native_cached is not None:
                    dotbracket_native = native_cached
                else:
                    stackings_native, pairings_native, res_native = bb.annotate(native_pdb_path)
                    dotbracket_native = bb.dot_bracket(pairings_native, res_native)[0]
                    computation_cache.set(session_id, 'native_dotbracket', (native_pdb_path,), dotbracket_native)
                
                fig = plot_arc(dotbracket_data, dotbracket_native[0])
                fig.write_html(os.path.join(plot_path, "arc_plot.html"))
                results['ARC'] = plotly_to_json(fig)
            
            elif plot_type == 'CONTACT_MAPS':
                # This still requires separate computation
                from create_plots import create_contact_maps
                create_contact_maps(traj_xtc_path, native_pdb_path, 1, session_id)
                path_save_figure = os.path.join('static', session_id, 'contact_map_plotly.png')
                results['CONTACT_MAPS'] = json.dumps({'path': path_save_figure})
            
            logger.info(f"{plot_type} processing time: {time.time() - plot_start:.2f}s")
        
        logger.info(f"Total batch annotate processing time: {time.time() - start_time:.2f}s")
        logger.info(f"Processed {len(plot_types)} plots in single batch operation")
        
        return results
        
    except Exception as exc:
        logger.error(f"Batch annotate processing failed: {str(exc)}")
        self.retry(exc=exc, countdown=60)

@app2.task(bind=True, max_retries=3)
@log_task_execution  
def batch_computation_coordinator(self, native_pdb_path, traj_xtc_path, download_path, plot_path, session_id, selected_plots):
    """
    Coordinates batch processing of all computations to minimize file I/O
    Groups related computations and executes them efficiently
    """
    try:
        start_time = time.time()
        logger.info(f"Starting batch computation coordinator for {len(selected_plots)} plots")
        
        results = {}
        
        # Group plots by computation type
        individual_plots = []
        annotate_plots = []
        
        for plot_type in selected_plots:
            if plot_type in ['SEC_STRUCTURE', 'DOTBRACKET', 'ARC', 'CONTACT_MAPS']:
                annotate_plots.append(plot_type)
            else:
                individual_plots.append(plot_type)
        
        # Process annotate-based plots in batch
        if annotate_plots:
            logger.info(f"Processing annotate-based plots in batch: {annotate_plots}")
            annotate_results = generate_annotate_plots_batch.apply_async(
                args=[native_pdb_path, traj_xtc_path, download_path, plot_path, session_id, annotate_plots]
            ).get()
            results.update(annotate_results)
        
        # Process individual plots with caching
        for plot_type in individual_plots:
            if plot_type == 'RMSD':
                result = generate_rmsd_plot_optimized.apply_async(
                    args=[native_pdb_path, traj_xtc_path, download_path, plot_path, session_id]
                ).get()
                results['RMSD'] = result
            
            elif plot_type == 'ERMSD':
                result = generate_ermsd_plot_optimized.apply_async(
                    args=[native_pdb_path, traj_xtc_path, download_path, plot_path, session_id]
                ).get()
                results['ERMSD'] = result
            
            elif plot_type == 'TORSION':
                # Get torsion parameters from session data if available
                torsion_params = None  # TODO: Extract from session data
                result = generate_torsion_plot_optimized.apply_async(
                    args=[native_pdb_path, traj_xtc_path, download_path, plot_path, session_id, torsion_params]
                ).get()
                results['TORSION'] = result
        
        total_time = time.time() - start_time
        logger.info(f"Batch computation coordinator completed in {total_time:.2f}s")
        logger.info(f"Average time per plot: {total_time/len(selected_plots):.2f}s")
        
        return results
        
    except Exception as exc:
        logger.error(f"Batch computation coordinator failed: {str(exc)}")
        self.retry(exc=exc, countdown=60)

# Keep original tasks for backward compatibility
from tasks_celery import *