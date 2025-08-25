#!/bin/bash
"""
Setup script for Plot RNA application
"""

set -e

echo "Setting up Plot RNA application..."

# Create virtual environment if it doesn't exist
if [ ! -d "venv" ]; then
    echo "Creating virtual environment..."
    python3 -m venv venv
fi

# Activate virtual environment
source venv/bin/activate

# Upgrade pip
pip install --upgrade pip

# Install requirements
echo "Installing Python dependencies..."
pip install -r requirements.txt

# Create necessary directories
echo "Creating directories..."
mkdir -p logs data/cache data/temp data/uploads static/uploads

# Set permissions
chmod 755 scripts/*.sh

# Install Redis if not available
if ! command -v redis-server &> /dev/null; then
    echo "Redis not found. Please install Redis server:"
    echo "Ubuntu/Debian: sudo apt-get install redis-server"
    echo "CentOS/RHEL: sudo yum install redis"
    echo "macOS: brew install redis"
fi

echo "Setup complete!"
echo ""
echo "To start the application:"
echo "1. Activate virtual environment: source venv/bin/activate"
echo "2. Start services: ./scripts/start.sh"
echo ""
echo "For development:"
echo "export FLASK_ENV=development"
echo "python wsgi.py"