import unittest
import os
import tempfile
import shutil
import numpy as np
import pandas as pd
from unittest.mock import Mock, patch, MagicMock
import logging

# Import the modules we want to test
from app7 import CalculationPlanner, CalculationRequirement
from create_plots import plot_rmsd, plot_ermsd, plot_dotbracket, plot_torsion, plot_rna_contact_map
import tasks_celery

# Configure logging for tests
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

class TestCalculationPlanner(unittest.TestCase):
    """Test the CalculationPlanner class"""
    
    def setUp(self):
        self.session_id = "test_session_123"
        self.planner = CalculationPlanner(self.session_id)
        
    def test_calculation_specs_exist(self):
        """Test that all expected calculation specs are defined"""
        expected_plots = [
            'RMSD', 'ERMSD', 'TORSION', 'SEC_STRUCTURE', 'DOTBRACKET',
            'ARC', 'CONTACT_MAPS', 'ANNOTATE', 'DS_MOTIF', 'SS_MOTIF',
            'JCOUPLING', 'ESCORE', 'LANDSCAPE', 'BASE_PAIRING'
        ]
        
        for plot in expected_plots:
            self.assertIn(plot, self.planner.CALCULATION_SPECS)
            spec = self.planner.CALCULATION_SPECS[plot]
            self.assertIsInstance(spec, CalculationRequirement)
            self.assertGreater(spec.estimated_time, 0)
            self.assertGreater(spec.memory_requirement, 0)
            self.assertIsInstance(spec.cpu_intensive, bool)
            self.assertIsInstance(spec.dependencies, list)
    
    def test_plan_calculations_basic(self):
        """Test basic calculation planning"""
        selected_plots = ['RMSD', 'ERMSD']
        planned = self.planner.plan_calculations(selected_plots)
        
        self.assertEqual(len(planned), 2)
        self.assertIn('RMSD', planned)
        self.assertIn('ERMSD', planned)
        
    def test_plan_calculations_with_dependencies(self):
        """Test calculation planning with dependencies"""
        selected_plots = ['LANDSCAPE']  # Has RMSD dependency
        planned = self.planner.plan_calculations(selected_plots)
        
        self.assertEqual(len(planned), 1)
        self.assertIn('LANDSCAPE', planned)
        
    def test_optimize_calculation_order(self):
        """Test calculation order optimization"""
        selected_plots = ['RMSD', 'ERMSD', 'TORSION', 'LANDSCAPE']
        planned = self.planner.plan_calculations(selected_plots)
        optimized_order = self.planner.optimize_calculation_order(planned)
        
        self.assertEqual(len(optimized_order), len(selected_plots))
        self.assertEqual(set(optimized_order), set(selected_plots))
        
    def test_unknown_plot_type(self):
        """Test handling of unknown plot types"""
        selected_plots = ['UNKNOWN_PLOT']
        planned = self.planner.plan_calculations(selected_plots)
        
        self.assertEqual(len(planned), 0)


class TestPlottingFunctions(unittest.TestCase):
    """Test the plotting functions in create_plots.py"""
    
    def setUp(self):
        # Create mock data for testing
        self.rmsd_data = np.random.rand(100) * 5  # 100 frames, RMSD values 0-5
        self.ermsd_data = np.random.rand(100) * 3  # 100 frames, ERMSD values 0-3
        self.sequence = ['A', 'U', 'G', 'C'] * 5  # 20 nucleotides
        self.angles = np.random.rand(100, 20, 7)  # 100 frames, 20 residues, 7 angles
        self.residues = [f"A{i+1}" for i in range(20)]
        
    def test_plot_rmsd_basic(self):
        """Test basic RMSD plotting"""
        fig = plot_rmsd(self.rmsd_data)
        self.assertIsNotNone(fig)
        self.assertEqual(len(fig.data), 1)
        self.assertEqual(len(fig.data[0].y), len(self.rmsd_data))
        
    def test_plot_ermsd_basic(self):
        """Test basic ERMSD plotting"""
        fig = plot_ermsd(self.ermsd_data)
        self.assertIsNotNone(fig)
        self.assertEqual(len(fig.data), 1)
        self.assertEqual(len(fig.data[0].y), len(self.ermsd_data))
        
    def test_plot_dotbracket_basic(self):
        """Test basic dot-bracket plotting"""
        dotbracket_data = ['((()))', '((.))', '....', '((()))', '((.))']
        fig = plot_dotbracket(dotbracket_data)
        self.assertIsNotNone(fig)
        self.assertGreater(len(fig.data), 0)
        
    def test_plot_torsion_basic(self):
        """Test basic torsion plotting"""
        fig = plot_torsion(self.angles, self.residues, 2)
        self.assertIsNotNone(fig)
        # Should have 14 traces (7 angles × 2 plot types)
        self.assertEqual(len(fig.data), 14)
        
    def test_plot_rna_contact_map_empty(self):
        """Test RNA contact map with empty data"""
        empty_df = pd.DataFrame()
        fig = plot_rna_contact_map(empty_df, self.sequence)
        self.assertIsNotNone(fig)
        
    def test_plot_rna_contact_map_with_data(self):
        """Test RNA contact map with sample data"""
        base_pairs_df = pd.DataFrame({
            'res_i': [1, 2, 3],
            'res_j': [10, 9, 8],
            'res_i_name': ['A1', 'U2', 'G3'],
            'res_j_name': ['U10', 'A9', 'C8'],
            'anno': ['WCc', 'GUc', 'WHc']
        })
        fig = plot_rna_contact_map(base_pairs_df, self.sequence)
        self.assertIsNotNone(fig)
        self.assertGreater(len(fig.data), 0)
        
    def test_plot_with_invalid_data(self):
        """Test plotting functions with invalid data"""
        # Test with None data
        with self.assertRaises(Exception):
            plot_rmsd(None)
            
        # Test with empty array
        with self.assertRaises(Exception):
            plot_rmsd(np.array([]))


class TestCeleryTasks(unittest.TestCase):
    """Test Celery task functions"""
    
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.session_id = "test_session"
        self.native_pdb = os.path.join(self.temp_dir, "native.pdb")
        self.traj_xtc = os.path.join(self.temp_dir, "traj.xtc")
        self.download_path = os.path.join(self.temp_dir, "download")
        self.plot_path = os.path.join(self.temp_dir, "plots")
        
        # Create directories
        os.makedirs(self.download_path, exist_ok=True)
        os.makedirs(self.plot_path, exist_ok=True)
        
        # Create dummy files
        open(self.native_pdb, 'a').close()
        open(self.traj_xtc, 'a').close()
        
    def tearDown(self):
        shutil.rmtree(self.temp_dir)
        
    @patch('tasks_celery.bb.rmsd')
    @patch('tasks_celery.plot_rmsd')
    @patch('tasks_celery.plotly_to_json')
    def test_generate_rmsd_plot_success(self, mock_plotly_to_json, mock_plot_rmsd, mock_bb_rmsd):
        """Test successful RMSD plot generation"""
        # Mock the barnaba calculation
        mock_bb_rmsd.return_value = np.random.rand(100)
        
        # Mock the plotting function
        mock_fig = Mock()
        mock_plot_rmsd.return_value = mock_fig
        
        # Mock the JSON conversion
        mock_plotly_to_json.return_value = '{"data": []}'
        
        # Create a mock task instance
        mock_task = Mock()
        mock_task.request.retries = 0
        mock_task.max_retries = 3
        
        # Call the function
        result = tasks_celery.generate_rmsd_plot(
            mock_task, self.native_pdb, self.traj_xtc, 
            self.download_path, self.plot_path, self.session_id
        )
        
        # Verify the result
        self.assertEqual(result, '{"data": []}')
        mock_bb_rmsd.assert_called_once()
        mock_plot_rmsd.assert_called_once()
        mock_plotly_to_json.assert_called_once()
        
    @patch('tasks_celery.bb.ermsd')
    @patch('tasks_celery.plot_ermsd')
    @patch('tasks_celery.plotly_to_json')
    def test_generate_ermsd_plot_success(self, mock_plotly_to_json, mock_plot_ermsd, mock_bb_ermsd):
        """Test successful ERMSD plot generation"""
        # Mock the barnaba calculation
        mock_bb_ermsd.return_value = np.random.rand(100)
        
        # Mock the plotting function
        mock_fig = Mock()
        mock_plot_ermsd.return_value = mock_fig
        
        # Mock the JSON conversion
        mock_plotly_to_json.return_value = '{"data": []}'
        
        # Create a mock task instance
        mock_task = Mock()
        mock_task.request.retries = 0
        mock_task.max_retries = 3
        
        # Call the function
        result = tasks_celery.generate_ermsd_plot(
            mock_task, self.native_pdb, self.traj_xtc,
            self.download_path, self.plot_path, self.session_id
        )
        
        # Verify the result
        self.assertEqual(result, '{"data": []}')
        mock_bb_ermsd.assert_called_once()
        mock_plot_ermsd.assert_called_once()
        mock_plotly_to_json.assert_called_once()


class TestIntegration(unittest.TestCase):
    """Integration tests for the calculation workflow"""
    
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.session_id = "integration_test"
        
    def tearDown(self):
        shutil.rmtree(self.temp_dir)
        
    def test_calculation_workflow_planning(self):
        """Test the complete calculation workflow planning"""
        planner = CalculationPlanner(self.session_id)
        
        # Test with a realistic set of plots
        selected_plots = ['RMSD', 'ERMSD', 'TORSION', 'CONTACT_MAPS', 'LANDSCAPE']
        
        # Plan calculations
        planned = planner.plan_calculations(selected_plots)
        self.assertEqual(len(planned), len(selected_plots))
        
        # Optimize order
        optimized_order = planner.optimize_calculation_order(planned)
        self.assertEqual(len(optimized_order), len(selected_plots))
        self.assertEqual(set(optimized_order), set(selected_plots))
        
    def test_resource_estimation(self):
        """Test resource estimation accuracy"""
        planner = CalculationPlanner(self.session_id)
        
        # Test resource calculations
        cpu_intensive_plots = ['RMSD', 'ERMSD', 'TORSION']
        planned = planner.plan_calculations(cpu_intensive_plots)
        
        total_cpu_time = sum(spec.estimated_time for spec in planned.values() if spec.cpu_intensive)
        total_memory = max(spec.memory_requirement for spec in planned.values())
        
        self.assertGreater(total_cpu_time, 0)
        self.assertGreater(total_memory, 0)


class TestErrorHandling(unittest.TestCase):
    """Test error handling in various scenarios"""
    
    def test_missing_files(self):
        """Test handling of missing input files"""
        planner = CalculationPlanner("test_session")
        
        # This should not crash, just return empty planning
        selected_plots = ['RMSD']
        planned = planner.plan_calculations(selected_plots)
        self.assertEqual(len(planned), 1)
        
    def test_invalid_session_data(self):
        """Test handling of invalid session data"""
        planner = CalculationPlanner("")  # Empty session ID
        
        selected_plots = ['RMSD', 'ERMSD']
        planned = planner.plan_calculations(selected_plots)
        self.assertEqual(len(planned), 2)
        
    def test_plot_function_error_handling(self):
        """Test error handling in plot functions"""
        # Test with malformed data
        with self.assertRaises(Exception):
            plot_rmsd("not_an_array")
            
        with self.assertRaises(Exception):
            plot_ermsd([])  # Empty list instead of numpy array


if __name__ == '__main__':
    # Create a test suite
    suite = unittest.TestSuite()
    
    # Add test classes
    suite.addTest(unittest.makeSuite(TestCalculationPlanner))
    suite.addTest(unittest.makeSuite(TestPlottingFunctions))
    suite.addTest(unittest.makeSuite(TestCeleryTasks))
    suite.addTest(unittest.makeSuite(TestIntegration))
    suite.addTest(unittest.makeSuite(TestErrorHandling))
    
    # Run the tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    # Print summary
    if result.wasSuccessful():
        print(f"\n✅ All {result.testsRun} tests passed!")
    else:
        print(f"\n❌ {len(result.failures)} failures, {len(result.errors)} errors out of {result.testsRun} tests")
        
    # Exit with appropriate code
    exit(0 if result.wasSuccessful() else 1)