"""
Plot RNA Configuration Settings
Centralized configuration management for the Plot RNA application
"""

import os
from typing import Dict, Any

# Base directory
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

class Config:
    """Base configuration class"""
    
    # Flask settings
    SECRET_KEY = os.environ.get('SECRET_KEY', 'your-secret-key-change-in-production')
    DEBUG = False
    
    # Celery settings
    CELERY_BROKER_URL = os.environ.get('CELERY_BROKER_URL', 'redis://localhost:6379/0')
    CELERY_RESULT_BACKEND = os.environ.get('CELERY_RESULT_BACKEND', 'redis://localhost:6379/0')
    CELERY_TASK_SERIALIZER = 'json'
    CELERY_RESULT_SERIALIZER = 'json'
    CELERY_ACCEPT_CONTENT = ['json']
    CELERY_TIMEZONE = 'UTC'
    CELERY_ENABLE_UTC = True
    
    # Application paths
    DATA_DIR = os.path.join(BASE_DIR, 'data')
    CACHE_DIR = os.path.join(DATA_DIR, 'cache')
    TEMP_DIR = os.path.join(DATA_DIR, 'temp')
    STATIC_DIR = os.path.join(BASE_DIR, 'static')
    TEMPLATES_DIR = os.path.join(BASE_DIR, 'templates')
    
    # Upload settings
    MAX_CONTENT_LENGTH = 500 * 1024 * 1024  # 500MB max file size
    UPLOAD_FOLDER = os.path.join(DATA_DIR, 'uploads')
    ALLOWED_EXTENSIONS = {'pdb', 'xtc', 'gro', 'dcd', 'trr'}
    
    # Computation settings
    MAX_PARALLEL_JOBS = int(os.environ.get('MAX_PARALLEL_JOBS', '4'))
    COMPUTATION_TIMEOUT = int(os.environ.get('COMPUTATION_TIMEOUT', '3600'))  # 1 hour
    CACHE_ENABLED = os.environ.get('CACHE_ENABLED', 'true').lower() == 'true'
    COMPRESSION_ENABLED = os.environ.get('COMPRESSION_ENABLED', 'true').lower() == 'true'
    
    # Performance settings
    ENABLE_PARALLEL_PROCESSING = True
    ENABLE_NUMBA_JIT = True
    MEMORY_LIMIT_GB = int(os.environ.get('MEMORY_LIMIT_GB', '8'))
    
    # Logging settings
    LOG_LEVEL = os.environ.get('LOG_LEVEL', 'INFO')
    LOG_FILE = os.path.join(BASE_DIR, 'logs', 'plot_rna.log')
    
    @staticmethod
    def init_app(app):
        """Initialize application with config"""
        # Create necessary directories
        for directory in [Config.DATA_DIR, Config.CACHE_DIR, Config.TEMP_DIR, 
                         Config.UPLOAD_FOLDER, os.path.dirname(Config.LOG_FILE)]:
            os.makedirs(directory, exist_ok=True)

class DevelopmentConfig(Config):
    """Development configuration"""
    DEBUG = True
    LOG_LEVEL = 'DEBUG'

class ProductionConfig(Config):
    """Production configuration"""
    DEBUG = False
    SECRET_KEY = os.environ.get('SECRET_KEY')
    
    # Security settings
    SESSION_COOKIE_SECURE = True
    SESSION_COOKIE_HTTPONLY = True
    SESSION_COOKIE_SAMESITE = 'Lax'
    
    # Performance settings
    SEND_FILE_MAX_AGE_DEFAULT = 31536000  # 1 year cache for static files

class TestingConfig(Config):
    """Testing configuration"""
    TESTING = True
    WTF_CSRF_ENABLED = False

# Configuration mapping
config = {
    'development': DevelopmentConfig,
    'production': ProductionConfig,
    'testing': TestingConfig,
    'default': DevelopmentConfig
}

def get_config() -> Config:
    """Get configuration based on environment"""
    env = os.environ.get('FLASK_ENV', 'development')
    return config.get(env, config['default'])