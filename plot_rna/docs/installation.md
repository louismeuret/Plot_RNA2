# Installation Guide

## System Requirements

### Minimum Requirements
- Python 3.10 or higher
- 4 GB RAM
- 2 CPU cores
- 10 GB disk space

### Recommended Requirements
- Python 3.11
- 16 GB RAM
- 8 CPU cores
- 50 GB disk space (for large trajectories)

### Dependencies
- Redis server
- Git
- C/C++ compiler (for NumPy/SciPy compilation)

## Installation Methods

### 1. Docker Installation (Recommended for Production)

```bash
# Clone repository
git clone <repository-url>
cd plot_rna

# Start with Docker Compose
docker-compose up -d

# Check status
docker-compose ps

# View logs
docker-compose logs -f plot-rna
```

### 2. Manual Installation

#### Ubuntu/Debian
```bash
# Install system dependencies
sudo apt-get update
sudo apt-get install python3.10 python3.10-venv python3.10-dev
sudo apt-get install redis-server gcc g++ gfortran
sudo apt-get install libhdf5-dev libnetcdf-dev

# Clone and setup
git clone <repository-url>
cd plot_rna
./scripts/setup.sh
```

#### CentOS/RHEL
```bash
# Install system dependencies
sudo yum install python3 python3-pip python3-devel
sudo yum install redis gcc gcc-c++ gcc-gfortran
sudo yum install hdf5-devel netcdf-devel

# Clone and setup
git clone <repository-url>
cd plot_rna
./scripts/setup.sh
```

#### macOS
```bash
# Install Homebrew (if not installed)
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install dependencies
brew install python@3.10 redis hdf5 netcdf

# Clone and setup
git clone <repository-url>
cd plot_rna
./scripts/setup.sh
```

### 3. Development Installation

```bash
# Clone repository
git clone <repository-url>
cd plot_rna

# Create virtual environment
python3 -m venv venv
source venv/bin/activate

# Install development dependencies
pip install -r requirements.txt
pip install -r requirements-dev.txt

# Setup environment
export FLASK_ENV=development
export PYTHONPATH=$(pwd)

# Run development server
python wsgi.py
```

## Configuration

### Environment Variables

Create a `.env` file:
```bash
FLASK_ENV=production
SECRET_KEY=your-very-secret-key-here
CELERY_BROKER_URL=redis://localhost:6379/0
CELERY_RESULT_BACKEND=redis://localhost:6379/0
MAX_PARALLEL_JOBS=4
LOG_LEVEL=INFO
```

### Redis Configuration

Default Redis configuration should work for most cases. For production:

```bash
# Edit Redis configuration
sudo nano /etc/redis/redis.conf

# Recommended settings:
# maxmemory 2gb
# maxmemory-policy allkeys-lru
# save 900 1
# save 300 10
# save 60 10000

# Restart Redis
sudo systemctl restart redis
```

## Verification

### Test Installation
```bash
# Activate environment
source venv/bin/activate

# Run tests
python -m pytest tests/

# Check services
curl http://localhost:8000/health
```

### Performance Test
```bash
# Run with example data
python -c "
from src.core.app import create_app
from config import get_config

app = create_app(get_config())
print('✅ Application loaded successfully')
"
```

## Troubleshooting

### Common Issues

**1. Redis Connection Error**
```bash
# Check Redis status
sudo systemctl status redis

# Start Redis
sudo systemctl start redis
```

**2. Permission Errors**
```bash
# Fix permissions
chmod +x scripts/*.sh
chown -R $USER:$USER logs/ data/
```

**3. Import Errors**
```bash
# Ensure Python path is set
export PYTHONPATH=$(pwd)

# Reinstall dependencies
pip install -r requirements.txt --force-reinstall
```

**4. Memory Issues**
```bash
# Increase memory limits
export MEMORY_LIMIT_GB=16

# Use swap if needed
sudo swapon /swapfile
```

## Production Deployment

### Nginx Setup
```bash
# Install Nginx
sudo apt-get install nginx

# Copy configuration
sudo cp deployment/nginx.conf /etc/nginx/sites-available/plot-rna
sudo ln -s /etc/nginx/sites-available/plot-rna /etc/nginx/sites-enabled/
sudo nginx -t
sudo systemctl reload nginx
```

### SSL Certificate
```bash
# Using Let's Encrypt
sudo apt-get install certbot python3-certbot-nginx
sudo certbot --nginx -d yourdomain.com
```

### Process Management
```bash
# Using systemd
sudo cp deployment/plot-rna.service /etc/systemd/system/
sudo systemctl enable plot-rna
sudo systemctl start plot-rna
```

## Updates

### Updating the Application
```bash
# Pull latest changes
git pull origin main

# Update dependencies
pip install -r requirements.txt --upgrade

# Restart services
docker-compose restart  # For Docker
# OR
sudo systemctl restart plot-rna  # For systemd
```

### Database Migrations
```bash
# Backup data
cp -r data/ data_backup/

# Apply updates
python scripts/migrate.py

# Verify
python scripts/verify_installation.py
```