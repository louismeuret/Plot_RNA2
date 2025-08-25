#!/bin/bash
"""
Production startup script for Plot RNA
"""

set -e

# Configuration
WORKERS=${WORKERS:-4}
PORT=${PORT:-8000}
HOST=${HOST:-0.0.0.0}
TIMEOUT=${TIMEOUT:-120}
KEEPALIVE=${KEEPALIVE:-5}

# Environment
export FLASK_ENV=production
export PYTHONPATH="${PYTHONPATH}:$(pwd)"

echo "Starting Plot RNA application..."
echo "Workers: $WORKERS"
echo "Port: $PORT"
echo "Host: $HOST"

# Start Redis if not running
if ! pgrep redis-server > /dev/null; then
    echo "Starting Redis server..."
    redis-server --daemonize yes
fi

# Start Celery worker in background
echo "Starting Celery worker..."
celery -A src.core.tasks worker --loglevel=info --detach

# Start Celery flower monitoring (optional)
if [ "$ENABLE_FLOWER" = "true" ]; then
    echo "Starting Celery Flower monitoring..."
    celery -A src.core.tasks flower --detach
fi

# Start Gunicorn
echo "Starting Gunicorn server..."
exec gunicorn \
    --bind $HOST:$PORT \
    --workers $WORKERS \
    --worker-class eventlet \
    --timeout $TIMEOUT \
    --keepalive $KEEPALIVE \
    --max-requests 1000 \
    --max-requests-jitter 100 \
    --preload \
    --access-logfile - \
    --error-logfile - \
    wsgi:app