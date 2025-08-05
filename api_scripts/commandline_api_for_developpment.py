import requests
import time
import json
import os
import zipfile
import io
import argparse
import sys
from urllib.parse import urljoin
import warnings
import webbrowser
from datetime import datetime
import subprocess
import shutil

class RNATrajectoryAnalysis:
    """
    Client for interacting with the RNA trajectory analysis web service.
    """

    def __init__(self, base_url="http://localhost:5000", verify_ssl=True):
        """
        Initialize the RNA trajectory analysis client.

        Args:
            base_url (str): The base URL of the RNA trajectory analysis service
            verify_ssl (bool): Whether to verify SSL certificates
        """
        self.base_url = base_url
        self.verify_ssl = verify_ssl
        self.session_id = None
        self.native_pdb = None
        self.traj_xtc = None
        self.timing_data = {}

        if not verify_ssl:
            warnings.warn("SSL verification is disabled. This is not recommended for production use.")

    def start_session(self):
        """
        Start a new analysis session.

        Returns:
            str: The session ID for the new session
        """
        response = requests.get(urljoin(self.base_url, "/get_session"), verify=self.verify_ssl)
        # Extract session ID from JSON response
        session_id = response.json().get("session_id")
        if not session_id:
            # Alternative approach - make a POST request to get a session ID
            response = requests.post(urljoin(self.base_url, "/"), verify=self.verify_ssl)
            session_id = response.cookies.get("session")

        self.session_id = session_id
        print(session_id)
        return session_id

    def upload_files(self, native_pdb_path, traj_xtc_path, analysis_options=None, frame_options=None, landscape_params=None):
        """
        Upload PDB and XTC files for analysis.

        Args:
            native_pdb_path (str): Path to the native PDB file
            traj_xtc_path (str): Path to the trajectory XTC file
            analysis_options (list): List of analysis types to perform (e.g., ["RMSD", "ERMSD", "TORSION"])
            frame_options (dict): Dict with frame selection options (e.g., {"n_frames": 10, "first_frame": 1, "last_frame": 100, "stride": 1})
            landscape_params (dict): Parameters for LANDSCAPE plot (e.g., {"x_param": "ERMSD", "y_param": "Q"})

        Returns:
            dict: Response from the server with the status of the upload
        """
        upload_start_time = time.time()
        print(f"[TIMING] Starting file upload at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        
        if not self.session_id:
            session_start_time = time.time()
            self.start_session()
            session_duration = time.time() - session_start_time
            self.timing_data['session_creation'] = {
                'duration_seconds': round(session_duration, 3),
                'description': 'Time taken to create a new analysis session with the server',
                'timestamp': datetime.now().isoformat()
            }
            print(f"[TIMING] Session creation completed in {session_duration:.3f} seconds")

        # Set default options if not provided
        if analysis_options is None:
            analysis_options = ["RMSD", "ERMSD"]

        if frame_options is None:
            frame_options = {"n_frames": 1, "first_frame": "", "last_frame": "", "stride": ""}

        if landscape_params is None:
            landscape_params = {"x_param": "ERMSD", "y_param": "Q"}

        # Prepare form data
        form_data = {
            "n_frames": frame_options.get("n_frames", 1),
            "firstFrame": frame_options.get("first_frame", ""),
            "lastFrame": frame_options.get("last_frame", ""),
            "stride": frame_options.get("stride", ""),
            "session_id": self.session_id
        }

        # Add selected plots to form data
        for plot in analysis_options:
            form_data[plot.lower()] = "on"

        # Add LANDSCAPE parameters if LANDSCAPE is selected
        if "LANDSCAPE" in analysis_options:
            x_param = landscape_params.get("x_param", "ERMSD")
            y_param = landscape_params.get("y_param", "Q")


            form_data["firstDimension"] = landscape_params.get("x_param", "ERMSD")
            form_data["secondDimension"] = landscape_params.get("y_param", "Q")
            form_data["landscape_stride"] = 1

            # Replace ERMSD with eRMSD if needed
            if x_param == "ERMSD":
                form_data["firstDimension"] = "eRMSD"
            if y_param == "ERMSD":
                form_data["secondDimension"] = "eRMSD"


        # Prepare files
        files = {
            "nativePdb": (os.path.basename(native_pdb_path), open(native_pdb_path, "rb"), "chemical/x-pdb"),
            "trajXtc": (os.path.basename(traj_xtc_path), open(traj_xtc_path, "rb"), "application/octet-stream")
        }

        # Store filenames for later use
        self.native_pdb = os.path.basename(native_pdb_path)
        self.traj_xtc = os.path.basename(traj_xtc_path)
        print(form_data)

        try:
            # Make the upload request
            request_start_time = time.time()
            print(f"[TIMING] Starting HTTP upload request at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
            
            response = requests.post(
                urljoin(self.base_url, "/upload-files"),
                files=files,
                data=form_data,
                verify=self.verify_ssl
            )
            
            request_duration = time.time() - request_start_time
            upload_duration = time.time() - upload_start_time
            
            # Store detailed timing information
            self.timing_data['file_upload'] = {
                'total_duration_seconds': round(upload_duration, 3),
                'http_request_duration_seconds': round(request_duration, 3),
                'description': f'Time taken to upload files ({os.path.basename(native_pdb_path)} and {os.path.basename(traj_xtc_path)}) to the server and receive response',
                'file_sizes': {
                    'pdb_file_mb': round(os.path.getsize(native_pdb_path) / 1024 / 1024, 3),
                    'xtc_file_mb': round(os.path.getsize(traj_xtc_path) / 1024 / 1024, 3)
                },
                'analysis_options': analysis_options or [],
                'frame_options': frame_options or {},
                'timestamp': datetime.now().isoformat()
            }
            
            print(f"[TIMING] File upload completed in {upload_duration:.3f} seconds (HTTP request: {request_duration:.3f}s)")

            # Print detailed response information for debugging
            #print(f"Response status code: {response.status_code}")
            #print(f"Response headers: {response.headers}")
            #print(f"Response text: {response.text}")

            # Close file handles
            for file_obj in files.values():
                file_obj[1].close()

            if response.status_code == 200:
                return {"status": "success", "session_id": self.session_id}
            else:
                return {"status": "error", "message": f"Upload failed with status code {response.status_code}", "response": response.text}

        except Exception as e:
            # Close file handles on error
            for file_obj in files.values():
                file_obj[1].close()
            return {"status": "error", "message": str(e)}

    def get_analysis_status(self, max_retries=30, retry_interval=2, show_progress=True):
        """
        Check the status of the analysis and wait until it's complete.

        Args:
            max_retries (int): Maximum number of retries
            retry_interval (int): Interval between retries in seconds
            show_progress (bool): Whether to show a progress animation

        Returns:
            dict: Analysis results or error information
        """
        if not self.session_id or not self.native_pdb or not self.traj_xtc:
            return {"status": "error", "message": "No active analysis session"}
            
        analysis_start_time = time.time()
        print(f"[TIMING] Starting analysis status monitoring at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

        # URL for retrieving results
        results_url = urljoin(
            self.base_url,
            f"/retrieve-results?session_id={self.session_id}"
        )

        # ASCII loading animations
        loading_frames = [
            "O       o O       o O       o\n| O   o | | O   o | | O   o |\n| | O | | | | O | | | | O | |\n| o   O | | o   O | | o   O |\no       O o       O o       O",
            "o       O o       O o       O\n| o   O | | o   O | | o   O |\n| | O | | | | O | | | | O | |\n| O   o | | O   o | | O   o |\nO       o O       o O       o",
            "O       o O       o O       o\n| O   o | | O   o | | O   o |\n| | O | | | | O | | | | O | |\n| o   O | | o   O | | o   O |\no       O o       O o       O",
            "o       O o       O o       O\n| o   O | | o   O | | o   O |\n| | O | | | | O | | | | O | |\n| O   o | | O   o | | O   o |\nO       o O       o O       o"
        ]

        # Then poll for results
        for attempt in range(max_retries):
            try:
                response = requests.get(results_url, verify=self.verify_ssl)

                if response.status_code == 200:
                    # Successful result
                    analysis_duration = time.time() - analysis_start_time
                    
                    # Store analysis completion timing
                    self.timing_data['analysis_completion'] = {
                        'duration_seconds': round(analysis_duration, 3),
                        'description': f'Time taken for server to complete analysis and become ready for result retrieval (waited {attempt + 1} attempts)',
                        'attempts_made': attempt + 1,
                        'retry_interval_seconds': retry_interval,
                        'timestamp': datetime.now().isoformat()
                    }
                    
                    print(f"[TIMING] Analysis completed in {analysis_duration:.3f} seconds after {attempt + 1} attempts")
                    
                    if show_progress:
                        # Clear the loading animation
                        sys.stdout.write("\033[F" * 5)  # Move cursor up 5 lines
                        sys.stdout.write("\033[K" * 5)  # Clear 5 lines
                        sys.stdout.flush()
                    return {
                        "status": "success",
                        "message": "Analysis complete",
                        "results_url": results_url
                    }
                elif response.status_code == 404:
                    # Results not ready yet
                    if show_progress:
                        # Show loading animation
                        frame_idx = attempt % len(loading_frames)
                        # Clear previous frame if not the first attempt
                        if attempt > 0:
                            sys.stdout.write("\033[F" * 5)  # Move cursor up 5 lines
                        print(loading_frames[frame_idx])
                        sys.stdout.flush()
                    else:
                        print(f"Analysis in progress... (attempt {attempt+1}/{max_retries})")
                    time.sleep(retry_interval)
                else:
                    # Unexpected status code
                    if show_progress:
                        # Clear the loading animation
                        sys.stdout.write("\033[F" * 5)  # Move cursor up 5 lines
                        sys.stdout.write("\033[K" * 5)  # Clear 5 lines
                        sys.stdout.flush()
                    return {
                        "status": "error",
                        "message": f"Got unexpected status code: {response.status_code}"
                    }

            except Exception as e:
                if show_progress:
                    # Clear the loading animation
                    sys.stdout.write("\033[F" * 5)  # Move cursor up 5 lines
                    sys.stdout.write("\033[K" * 5)  # Clear 5 lines
                    sys.stdout.flush()
                return {"status": "error", "message": str(e)}

        if show_progress:
            # Clear the loading animation
            sys.stdout.write("\033[F" * 5)  # Move cursor up 5 lines
            sys.stdout.write("\033[K" * 5)  # Clear 5 lines
            sys.stdout.flush()
        analysis_timeout_duration = time.time() - analysis_start_time
        self.timing_data['analysis_timeout'] = {
            'duration_seconds': round(analysis_timeout_duration, 3),
            'description': f'Analysis timed out after {max_retries} attempts over {analysis_timeout_duration:.3f} seconds',
            'max_retries': max_retries,
            'retry_interval_seconds': retry_interval,
            'timestamp': datetime.now().isoformat()
        }
        print(f"[TIMING] Analysis timed out after {analysis_timeout_duration:.3f} seconds")
        return {"status": "error", "message": "Analysis timed out"}

    def download_plot_data(self, plot_id):
        """
        Download data for a specific plot.

        Args:
            plot_id (str): ID of the plot (e.g., "RMSD", "ERMSD")

        Returns:
            dict: Dictionary with file content or error message
        """
        if not self.session_id:
            return {"status": "error", "message": "No active session"}

        download_url = urljoin(
            self.base_url,
            f"/download/plot_data/{self.session_id}/{plot_id}"
        )

        try:
            response = requests.get(download_url, verify=self.verify_ssl)

            if response.status_code == 200:
                filename = plot_id + "_data.zip"
                content_type = response.headers.get('Content-Type', '')

                if 'zip' in content_type:
                    # Return as zip file bytes
                    return {
                        "status": "success",
                        "filename": filename,
                        "content": response.content,
                        "content_type": content_type
                    }
                else:
                    # Return as raw data
                    return {
                        "status": "success",
                        "filename": plot_id + "_data.csv",
                        "content": response.content,
                        "content_type": content_type
                    }
            else:
                return {"status": "error", "message": f"Download failed with status code {response.status_code}"}

        except Exception as e:
            return {"status": "error", "message": str(e)}

    def download_plot(self, plot_id):
        """
        Download plot visualization for a specific analysis.

        Args:
            plot_id (str): ID of the plot (e.g., "RMSD", "ERMSD")

        Returns:
            dict: Dictionary with file content or error message
        """
        if not self.session_id:
            return {"status": "error", "message": "No active session"}

        download_url = urljoin(
            self.base_url,
            f"/download/plot/{self.session_id}/{plot_id}"
        )

        try:
            response = requests.get(download_url, verify=self.verify_ssl)

            if response.status_code == 200:
                filename = plot_id + "_plots.zip"
                content_type = response.headers.get('Content-Type', '')

                return {
                    "status": "success",
                    "filename": filename,
                    "content": response.content,
                    "content_type": content_type
                }
            else:
                return {"status": "error", "message": f"Download failed with status code {response.status_code}"}

        except Exception as e:
            return {"status": "error", "message": str(e)}

    def save_file(self, file_data, output_dir="./outputs"):
        """
        Save downloaded file data to disk.

        Args:
            file_data (dict): File data returned from download_plot or download_plot_data
            output_dir (str): Directory to save the file

        Returns:
            str: Path to the saved file or error message
        """
        if file_data.get("status") != "success":
            return file_data.get("message", "Unknown error")

        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)

        # Save the file
        output_path = os.path.join(output_dir, file_data["filename"])
        with open(output_path, "wb") as f:
            f.write(file_data["content"])

        return output_path

    def extract_zip(self, zip_data, output_dir="./outputs"):
        """
        Extract contents of a zip file returned from download_plot or download_plot_data.

        Args:
            zip_data (dict): Zip file data returned from download functions
            output_dir (str): Directory to extract the files to

        Returns:
            list: List of extracted file paths or error message
        """
        if zip_data.get("status") != "success":
            return zip_data.get("message", "Unknown error")

        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)

        # Extract the zip file
        try:
            extracted_files = []
            with zipfile.ZipFile(io.BytesIO(zip_data["content"])) as z:
                z.extractall(output_dir)
                extracted_files = [os.path.join(output_dir, name) for name in z.namelist()]

            return extracted_files
        except Exception as e:
            return [f"Error extracting zip: {str(e)}"]

    def get_available_plots(self):
        """
        Get a list of available plot types for RNA analysis.

        Returns:
            list: List of available plot types
        """
        return [
            "RMSD",
            "ERMSD",
            "CONTACT_MAPS",
            "SEC_STRUCTURE",
            "DOTBRACKET",
            "ARC",
            "ANNOTATE",
            "DS_MOTIF",
            "SS_MOTIF",
            "JCOUPLING",
            "ESCORE",
            "LANDSCAPE",
            "BASE_PAIRING"
        ]

    def get_landscape_parameters(self):
        """
        Get available parameters for LANDSCAPE plots.

        Returns:
            dict: Dictionary with available X and Y parameters
        """
        return {
            "x_params": ["RMSD", "ERMSD", "TORSION", "Q"],
            "y_params": ["RMSD", "ERMSD", "TORSION", "Q"]
        }
        
    def save_timing_data(self, output_file=None):
        """
        Save timing data to a JSON file with comprehensive information about the analysis procedure.
        
        Args:
            output_file (str): Path to the output JSON file. If None, generates timestamped filename.
            
        Returns:
            str: Path to the saved timing file
        """
        # Generate timestamped filename if not provided
        if output_file is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_file = f"analysis_timing_{timestamp}.json"
        # Calculate total procedure time by summing all relevant timing phases
        total_time = 0
        phase_times = {}
        
        # Define all possible timing phases and their duration keys
        timing_phases = {
            'session_creation': 'duration_seconds',
            'file_upload': 'total_duration_seconds', 
            'analysis_completion': 'duration_seconds',
            'analysis_timeout': 'duration_seconds',
            'selenium_page_load': 'total_duration_seconds',
            'browser_operation_total': 'total_duration_seconds',
            'webbrowser_fallback': 'total_duration_seconds',
            'overall_procedure': 'duration_seconds'
        }
        
        # Sum up all available timing data
        for phase, duration_key in timing_phases.items():
            if phase in self.timing_data and duration_key in self.timing_data[phase]:
                phase_time = self.timing_data[phase][duration_key]
                phase_times[phase] = phase_time
                # Only add to total if it's not already included in overall_procedure
                if phase != 'overall_procedure' and phase != 'browser_operation_total':
                    total_time += phase_time
        
        # Use overall_procedure time if available (most accurate), otherwise use calculated sum
        if 'overall_procedure' in phase_times:
            total_time = phase_times['overall_procedure']
            
        # Create comprehensive timing report
        timing_report = {
            'procedure_summary': {
                'total_duration_seconds': round(total_time, 3),
                'total_duration_formatted': f"{int(total_time // 60)}m {total_time % 60:.3f}s" if total_time >= 60 else f"{total_time:.3f}s",
                'description': 'Complete time from session start to analysis completion including browser opening',
                'session_id': self.session_id,
                'files_analyzed': {
                    'pdb_file': self.native_pdb,
                    'xtc_file': self.traj_xtc
                },
                'procedure_completed_at': datetime.now().isoformat(),
                'phases_included': list(phase_times.keys()),
                'timing_method': 'selenium_webdriver' if 'selenium_page_load' in self.timing_data else 'standard'
            },
            'detailed_timings': self.timing_data,
            'performance_breakdown': {},
            'comparison_metrics': {
                'file_sizes_mb': {},
                'browser_performance': {},
                'analysis_efficiency': {}
            }
        }
        
        # Add performance breakdown percentages and comparison metrics
        if total_time > 0:
            for phase, data in self.timing_data.items():
                # Get duration from the appropriate field
                duration = None
                if 'duration_seconds' in data:
                    duration = data['duration_seconds']
                elif 'total_duration_seconds' in data:
                    duration = data['total_duration_seconds']
                
                if duration is not None:
                    percentage = (duration / total_time) * 100
                    timing_report['performance_breakdown'][phase] = {
                        'percentage_of_total': round(percentage, 2),
                        'duration_seconds': duration,
                        'description': data.get('description', f'Time spent in {phase} phase')
                    }
        
        # Add comparison metrics for analysis
        if 'file_upload' in self.timing_data and 'file_sizes' in self.timing_data['file_upload']:
            file_sizes = self.timing_data['file_upload']['file_sizes']
            timing_report['comparison_metrics']['file_sizes_mb'] = file_sizes
            
            # Calculate upload speed if available
            if 'total_duration_seconds' in self.timing_data['file_upload']:
                upload_time = self.timing_data['file_upload']['total_duration_seconds']
                total_size_mb = file_sizes.get('pdb_file_mb', 0) + file_sizes.get('xtc_file_mb', 0)
                if upload_time > 0 and total_size_mb > 0:
                    timing_report['comparison_metrics']['upload_speed_mbps'] = round(total_size_mb / upload_time, 3)
        
        # Add browser performance metrics if available
        if 'selenium_page_load' in self.timing_data:
            browser_data = self.timing_data['selenium_page_load']
            timing_report['comparison_metrics']['browser_performance'] = {
                'browser_used': browser_data.get('browser_used', 'unknown'),
                'page_load_speed_seconds': browser_data.get('full_load_duration_seconds', 0),
                'page_size_kb': browser_data.get('page_source_size_kb', 0),
                'navigation_speed_seconds': browser_data.get('navigation_duration_seconds', 0)
            }
        elif 'webbrowser_fallback' in self.timing_data:
            timing_report['comparison_metrics']['browser_performance'] = {
                'browser_used': 'system_default',
                'method': 'webbrowser_fallback',
                'launch_time_seconds': self.timing_data['webbrowser_fallback'].get('browser_launch_duration_seconds', 0)
            }
        
        # Add analysis efficiency metrics
        if 'analysis_completion' in self.timing_data:
            analysis_data = self.timing_data['analysis_completion']
            timing_report['comparison_metrics']['analysis_efficiency'] = {
                'completion_time_seconds': analysis_data.get('duration_seconds', 0),
                'attempts_required': analysis_data.get('attempts_made', 0),
                'retry_interval_seconds': analysis_data.get('retry_interval_seconds', 0),
                'efficiency_score': round(1.0 / max(analysis_data.get('attempts_made', 1), 1), 3)  # Lower attempts = higher efficiency
            }
        
        # Save to file
        with open(output_file, 'w') as f:
            json.dump(timing_report, f, indent=2)
            
        print(f"[TIMING] Comprehensive timing data saved to {output_file}")
        print(f"[TIMING] Total procedure time: {timing_report['procedure_summary']['total_duration_formatted']}")
        
        return output_file
        
    def _check_selenium_available(self):
        """
        Check if Selenium is available and suggest installation if not.
        
        Returns:
            bool: True if Selenium is available, False otherwise
        """
        try:
            import selenium
            return True
        except ImportError:
            print("[SELENIUM] Selenium not available. Install with: pip install selenium")
            return False
    
    def _setup_webdriver(self, visible=False):
        """
        Set up WebDriver with automatic driver management.
        Tries Chrome first, then Firefox, with fallback options.
        
        Args:
            visible (bool): Whether to run browser in visible mode (not headless)
        
        Returns:
            tuple: (driver, browser_name) or (None, None) if setup fails
        """
        try:
            from selenium import webdriver
            from selenium.webdriver.chrome.service import Service as ChromeService
            from selenium.webdriver.firefox.service import Service as FirefoxService
            from selenium.webdriver.chrome.options import Options as ChromeOptions
            from selenium.webdriver.firefox.options import Options as FirefoxOptions
        except ImportError:
            return None, None
            
        # Try Chrome first
        try:
            # Check if chromedriver is available
            chrome_options = ChromeOptions()
            if not visible:
                chrome_options.add_argument('--headless=new')  # Run in background
            chrome_options.add_argument('--no-sandbox')
            chrome_options.add_argument('--disable-dev-shm-usage')
            if not visible:
                chrome_options.add_argument('--disable-gpu')
            chrome_options.add_argument('--disable-web-security')
            
            # Try to use chromedriver from PATH first
            driver = webdriver.Chrome(options=chrome_options)
            return driver, "Chrome"
            
        except Exception as chrome_error:
            print(f"[SELENIUM] Chrome setup failed: {chrome_error}")
            
        # Try Firefox as fallback
        try:
            firefox_options = FirefoxOptions()
            if not visible:
                firefox_options.add_argument('--headless')
            
            driver = webdriver.Firefox(options=firefox_options)
            return driver, "Firefox"
            
        except Exception as firefox_error:
            print(f"[SELENIUM] Firefox setup failed: {firefox_error}")
            
        return None, None
    
    def open_results_in_browser(self, results_url, use_selenium=True, fallback_to_webbrowser=True, visible_browser=False):
        """
        Open the analysis results page and measure actual page load completion time.
        Uses Selenium WebDriver to detect when the page finishes loading completely.
        
        Args:
            results_url (str): URL to the results page
            use_selenium (bool): Whether to try Selenium for accurate timing
            fallback_to_webbrowser (bool): Whether to fallback to webbrowser if Selenium fails
            visible_browser (bool): Whether to run browser in visible mode (not headless)
            
        Returns:
            bool: True if browser opened successfully, False otherwise
        """
        total_start_time = time.time()
        print(f"[BROWSER] Opening results page with full load timing: {results_url}")
        
        # Try Selenium first for accurate page load timing
        if use_selenium and self._check_selenium_available():
            return self._open_with_selenium(results_url, total_start_time, fallback_to_webbrowser, visible_browser)
        
        # Fallback to traditional method
        if fallback_to_webbrowser:
            print("[BROWSER] Using fallback method (webbrowser module)")
            return self._open_with_webbrowser(results_url, total_start_time)
        
        return False
    
    def _open_with_selenium(self, results_url, total_start_time, fallback_to_webbrowser, visible_browser=False):
        """
        Open browser using Selenium WebDriver for accurate page load timing.
        """
        try:
            from selenium.webdriver.support.ui import WebDriverWait
            from selenium.webdriver.support import expected_conditions as EC
            from selenium.webdriver.common.by import By
            from selenium.common.exceptions import TimeoutException, WebDriverException
            
            # Set up WebDriver
            setup_start = time.time()
            driver, browser_name = self._setup_webdriver(visible=visible_browser)
            
            if not driver:
                print("[SELENIUM] WebDriver setup failed")
                if fallback_to_webbrowser:
                    return self._open_with_webbrowser(results_url, total_start_time)
                return False
                
            setup_duration = time.time() - setup_start
            mode = "visible" if visible_browser else "headless"
            print(f"[SELENIUM] {browser_name} WebDriver initialized in {mode} mode in {setup_duration:.3f} seconds")
            
            try:
                # Navigate to page and measure load time
                page_start = time.time()
                print(f"[TIMING] Starting page navigation at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
                
                driver.get(results_url)
                navigation_duration = time.time() - page_start
                
                # Wait for page to be fully loaded
                load_start = time.time()
                print("[TIMING] Waiting for page load completion...")
                
                # Wait for document ready state
                WebDriverWait(driver, 30).until(
                    lambda driver: driver.execute_script("return document.readyState") == "complete"
                )
                
                # Additional wait for any dynamic content/JavaScript
                # Wait for body to be present and any loading indicators to disappear
                WebDriverWait(driver, 10).until(
                    EC.presence_of_element_located((By.TAG_NAME, "body"))
                )
                
                # Try to wait for common loading indicators to disappear
                try:
                    # Wait a bit more for any AJAX/dynamic content to load
                    time.sleep(1)  # Small buffer for final rendering
                except:
                    pass
                
                full_load_duration = time.time() - load_start
                total_page_duration = time.time() - page_start
                total_browser_duration = time.time() - total_start_time
                
                # Get page information
                page_title = driver.title
                page_source_size = len(driver.page_source.encode('utf-8'))
                
                # Store comprehensive timing data
                self.timing_data['selenium_page_load'] = {
                    'total_duration_seconds': round(total_page_duration, 3),
                    'navigation_duration_seconds': round(navigation_duration, 3),
                    'full_load_duration_seconds': round(full_load_duration, 3),
                    'webdriver_setup_duration_seconds': round(setup_duration, 3),
                    'browser_used': browser_name,
                    'page_title': page_title,
                    'page_source_size_kb': round(page_source_size / 1024, 2),
                    'description': f'Complete page load timing using {browser_name} WebDriver - measures from navigation start to full page load completion',
                    'url': results_url,
                    'timestamp': datetime.now().isoformat()
                }
                
                self.timing_data['browser_operation_total'] = {
                    'total_duration_seconds': round(total_browser_duration, 3),
                    'description': 'Total time for complete browser operation including WebDriver setup, navigation, and full page load',
                    'method': 'selenium_webdriver',
                    'timestamp': datetime.now().isoformat()
                }
                
                print(f"[SELENIUM] Page navigation completed in {navigation_duration:.3f} seconds")
                print(f"[SELENIUM] Full page load completed in {full_load_duration:.3f} seconds")
                print(f"[SELENIUM] Total page operation: {total_page_duration:.3f} seconds")
                print(f"[SELENIUM] Page title: {page_title}")
                print(f"[TIMING] Complete browser operation finished in {total_browser_duration:.3f} seconds")
                
                # Keep browser open for user interaction
                print("[SELENIUM] Browser opened successfully. Page is fully loaded!")
                print("[SELENIUM] Browser will remain open for viewing results.")
                
                # Don't close the driver immediately - let user interact with the page
                # driver.quit() will be called when the script ends
                
                return True
                
            except TimeoutException:
                load_duration = time.time() - total_start_time
                self.timing_data['selenium_page_load'] = {
                    'total_duration_seconds': round(load_duration, 3),
                    'error': 'Page load timed out after 30 seconds',
                    'browser_used': browser_name,
                    'description': f'Page load attempt using {browser_name} WebDriver timed out',
                    'url': results_url,
                    'timestamp': datetime.now().isoformat()
                }
                print("[SELENIUM] Page load timed out, but browser window is open")
                return True  # Browser is still open even if load timed out
                
            except Exception as selenium_error:
                print(f"[SELENIUM] Error during page load: {selenium_error}")
                if fallback_to_webbrowser:
                    return self._open_with_webbrowser(results_url, total_start_time)
                return False
            
        except ImportError as import_error:
            print(f"[SELENIUM] Import error: {import_error}")
            if fallback_to_webbrowser:
                return self._open_with_webbrowser(results_url, total_start_time)
            return False
        
        except Exception as e:
            print(f"[SELENIUM] Unexpected error: {e}")
            if fallback_to_webbrowser:
                return self._open_with_webbrowser(results_url, total_start_time)
            return False
    
    def _open_with_webbrowser(self, results_url, total_start_time):
        """
        Fallback method using webbrowser module.
        """
        try:
            print("[BROWSER] Using webbrowser module (fallback method)")
            browser_start = time.time()
            webbrowser.open(results_url)
            browser_duration = time.time() - browser_start
            total_duration = time.time() - total_start_time
            
            self.timing_data['webbrowser_fallback'] = {
                'browser_launch_duration_seconds': round(browser_duration, 3),
                'total_duration_seconds': round(total_duration, 3),
                'description': 'Fallback browser opening using webbrowser module - cannot measure actual page load completion',
                'url': results_url,
                'timestamp': datetime.now().isoformat()
            }
            
            print(f"[BROWSER] Browser launched in {browser_duration:.3f} seconds (fallback method)")
            print("[BROWSER] Note: Actual page load completion time not measured with fallback method")
            
            return True
            
        except Exception as e:
            total_duration = time.time() - total_start_time
            self.timing_data['browser_failure'] = {
                'total_duration_seconds': round(total_duration, 3),
                'error': str(e),
                'description': 'Failed to open browser with any method',
                'url': results_url,
                'timestamp': datetime.now().isoformat()
            }
            print(f"[BROWSER] Failed to open browser: {str(e)}")
            return False

def main():
    overall_start_time = time.time()
    print(f"[TIMING] Starting RNA trajectory analysis procedure at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Set up argument parser
    parser = argparse.ArgumentParser(description='RNA Trajectory Analysis Tool')
    parser.add_argument('trajectory', help='Path to the trajectory XTC file')
    parser.add_argument('topology', help='Path to the topology PDB file')
    parser.add_argument('--analyses', nargs='+', default=['RMSD', 'ERMSD'],
                      choices=['RMSD', 'ERMSD', 'CONTACT_MAPS', 'TORSION', 'SEC_STRUCTURE',
                              'DOTBRACKET', 'ARC', 'LANDSCAPE', 'BASE_PAIRING', 'ALL'],
                      help='List of analyses to run (or "ALL" to run everything)')
    parser.add_argument('--n_frames', type=int, default=1,
                       help='Number of frames to analyze')
    parser.add_argument('--first_frame', type=int, default=0,
                       help='First frame to analyze')
    parser.add_argument('--last_frame', type=int, default=None,
                       help='Last frame to analyze')
    parser.add_argument('--stride', type=int, default=1,
                       help='Stride for frame selection')
    parser.add_argument('--server', default='http://localhost:4242',
                       help='Server URL for the analysis service')
    parser.add_argument('--no-verify-ssl', action='store_true',
                       help='Disable SSL certificate verification (not recommended for production)')
    parser.add_argument('--landscape-x', choices=['RMSD', 'ERMSD', 'TORSION', 'Q'], default='ERMSD',
                       help='X-axis parameter for LANDSCAPE plot')
    parser.add_argument('--landscape-y', choices=['RMSD', 'ERMSD', 'TORSION', 'Q'], default='Q',
                       help='Y-axis parameter for LANDSCAPE plot')
    parser.add_argument('--no-progress', action='store_true',
                       help='Disable progress animation')
    parser.add_argument('--timing-file', default=None,
                       help='File to save timing data (default: analysis_timing_YYYYMMDD_HHMMSS.json)')
    parser.add_argument('--no-browser', action='store_true',
                       help='Do not automatically open results in browser')
    parser.add_argument('--no-selenium', action='store_true',
                       help='Disable Selenium WebDriver and use fallback browser opening')
    parser.add_argument('--visible-browser', action='store_true',
                       help='Run browser in visible mode instead of headless (when using Selenium)')

    args = parser.parse_args()

    all_analyses = ['RMSD', 'ERMSD', 'CONTACT_MAPS',
                    'SEC_STRUCTURE', 'DOTBRACKET', 'ARC',
                    'LANDSCAPE', 'BASE_PAIRING']

    if 'ALL' in args.analyses:
        args.analyses = all_analyses

    # Create frame options dictionary
    frame_options = {
        'n_frames': args.n_frames,
        'first_frame': args.first_frame,
        'last_frame': args.last_frame,
        'stride': args.stride
    }

    # Create landscape parameters dictionary
    landscape_params = {
        'x_param': args.landscape_x,
        'y_param': args.landscape_y
    }

    # Initialize the client
    client = RNATrajectoryAnalysis(base_url=args.server, verify_ssl=not args.no_verify_ssl)
    
    # Store overall procedure start time
    client.timing_data['overall_procedure'] = {
        'start_time': datetime.now().isoformat(),
        'description': 'Complete RNA trajectory analysis procedure from initialization to results availability'
    }

    # Upload files and start analysis
    print("Uploading files and starting analysis...")
    upload_result = client.upload_files(
        native_pdb_path=args.topology,
        traj_xtc_path=args.trajectory,
        analysis_options=args.analyses,
        frame_options=frame_options,
        landscape_params=landscape_params
    )

    if upload_result['status'] != 'success':
        print(f"Error uploading files: {upload_result['message']}")
        # Save timing data even on failure
        client.save_timing_data(args.timing_file)
        return

    # Wait for analysis to complete
    print("Waiting for analysis to complete...")
    status_result = client.get_analysis_status(
        max_retries=60,
        retry_interval=5,
        show_progress=not args.no_progress
    )

    if status_result['status'] != 'success':
        print(f"Error during analysis: {status_result['message']}")
        # Still save timing data even on failure
        overall_duration = time.time() - overall_start_time
        client.timing_data['overall_procedure']['duration_seconds'] = round(overall_duration, 3)
        client.timing_data['overall_procedure']['end_time'] = datetime.now().isoformat()
        timing_file_path = client.save_timing_data(args.timing_file)
        print(f"Timing data saved to: {timing_file_path}")
        return

    # Print the results URL
    print("\n\n\nAnalysis complete!")
    print(f"Access raw data at: {status_result['results_url']}")
    
    # Open results in browser unless disabled (BEFORE saving timing data)
    if not args.no_browser:
        print("[TIMING] Starting browser operation...")
        client.open_results_in_browser(
            status_result['results_url'], 
            use_selenium=not args.no_selenium,
            fallback_to_webbrowser=True,
            visible_browser=args.visible_browser
        )
    
    # Calculate and store overall procedure timing (AFTER browser operation)
    overall_duration = time.time() - overall_start_time
    client.timing_data['overall_procedure']['duration_seconds'] = round(overall_duration, 3)
    client.timing_data['overall_procedure']['end_time'] = datetime.now().isoformat()
    
    # Save timing data to JSON file (AFTER all operations complete)
    timing_file_path = client.save_timing_data(args.timing_file)
    print(f"Timing data saved to: {timing_file_path}")
    
    # Print final summary
    print(f"[TIMING] Complete procedure finished in {overall_duration:.3f} seconds")

if __name__ == "__main__":
    main()
