"""
WSGI entry point for Plot RNA application
This file is used by production servers like Gunicorn
"""

import os
import sys

# Add the project directory to Python path
project_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, project_dir)
# Add parent directory for FoldingAnalysis module
sys.path.insert(0, os.path.join(project_dir, '..'))

from src.core.app import app, socketio

if __name__ == "__main__":
    socketio.run(app, host='0.0.0.0', port=4242, allow_unsafe_werkzeug=True)
