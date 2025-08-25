"""
Logging configuration for Plot RNA application
"""

import os
import logging
import logging.config
from datetime import datetime

def setup_logging(config):
    """Setup logging configuration"""
    
    # Ensure log directory exists
    log_dir = os.path.dirname(config.LOG_FILE)
    os.makedirs(log_dir, exist_ok=True)
    
    # Logging configuration
    LOGGING_CONFIG = {
        'version': 1,
        'disable_existing_loggers': False,
        'formatters': {
            'default': {
                'format': '[%(asctime)s] %(levelname)s in %(module)s: %(message)s',
                'datefmt': '%Y-%m-%d %H:%M:%S'
            },
            'detailed': {
                'format': '[%(asctime)s] %(levelname)s [%(name)s:%(lineno)d] %(message)s',
                'datefmt': '%Y-%m-%d %H:%M:%S'
            },
            'performance': {
                'format': '[%(asctime)s] PERF [%(name)s] %(message)s',
                'datefmt': '%Y-%m-%d %H:%M:%S'
            }
        },
        'handlers': {
            'console': {
                'class': 'logging.StreamHandler',
                'level': config.LOG_LEVEL,
                'formatter': 'default',
                'stream': 'ext://sys.stdout'
            },
            'file': {
                'class': 'logging.handlers.RotatingFileHandler',
                'level': config.LOG_LEVEL,
                'formatter': 'detailed',
                'filename': config.LOG_FILE,
                'maxBytes': 10485760,  # 10MB
                'backupCount': 5,
                'encoding': 'utf8'
            },
            'performance': {
                'class': 'logging.handlers.RotatingFileHandler',
                'level': 'INFO',
                'formatter': 'performance',
                'filename': os.path.join(log_dir, 'performance.log'),
                'maxBytes': 10485760,  # 10MB
                'backupCount': 3,
                'encoding': 'utf8'
            },
            'error': {
                'class': 'logging.handlers.RotatingFileHandler',
                'level': 'ERROR',
                'formatter': 'detailed',
                'filename': os.path.join(log_dir, 'errors.log'),
                'maxBytes': 10485760,  # 10MB
                'backupCount': 5,
                'encoding': 'utf8'
            }
        },
        'loggers': {
            '': {  # Root logger
                'handlers': ['console', 'file', 'error'],
                'level': config.LOG_LEVEL,
                'propagate': False
            },
            'plot_rna': {
                'handlers': ['console', 'file', 'error'],
                'level': config.LOG_LEVEL,
                'propagate': False
            },
            'plot_rna.performance': {
                'handlers': ['performance'],
                'level': 'INFO',
                'propagate': False
            },
            'celery': {
                'handlers': ['file', 'error'],
                'level': 'INFO',
                'propagate': False
            },
            'werkzeug': {
                'handlers': ['file'],
                'level': 'WARNING',
                'propagate': False
            }
        }
    }
    
    # Apply logging configuration
    logging.config.dictConfig(LOGGING_CONFIG)
    
    # Log startup message
    logger = logging.getLogger('plot_rna')
    logger.info(f"Plot RNA logging initialized at {datetime.now()}")
    logger.info(f"Log level: {config.LOG_LEVEL}")
    logger.info(f"Log file: {config.LOG_FILE}")

def get_performance_logger():
    """Get performance logger for timing measurements"""
    return logging.getLogger('plot_rna.performance')