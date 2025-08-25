#!/usr/bin/env python3
"""
Start Celery worker with proper path configuration
"""
import os
import sys

# Add parent directories to path for imports
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.join(current_dir, '..')
sys.path.insert(0, parent_dir)

# Import and start Celery
from tasks_celery import app2

if __name__ == '__main__':
    # Start the Celery worker
    app2.start([
        'worker',
        '--loglevel=info',
        '--concurrency=2',
        '--prefetch-multiplier=1'
    ])