# Plot RNA - Advanced RNA Trajectory Analysis Platform

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Docker](https://img.shields.io/badge/docker-supported-blue.svg)](https://www.docker.com/)

Plot RNA is a comprehensive web-based platform for analyzing and visualizing RNA molecular dynamics trajectories. It provides advanced computational tools for RMSD analysis, energy landscapes, contact maps, and structural motif analysis with high-performance optimization.

## ✨ Features

### Core Analysis Tools
- **RMSD Analysis** - Root Mean Square Deviation calculations with caching
- **eRMSD Analysis** - Enhanced RMSD with structural alignment
- **Energy Landscapes** - 3D/2D free energy surface visualization
- **Contact Maps** - Residue-residue contact analysis
- **Secondary Structure** - RNA base pairing and motif analysis
- **Torsion Analysis** - Backbone and side-chain dihedral angles

### Performance Optimizations
- **Parallel Processing** - Multi-core CPU utilization
- **Smart Caching** - Eliminates redundant computations (30-75% speedup)
- **Vectorized Computing** - NumPy/Numba optimization (3-10x faster)
- **Memory Efficiency** - Streaming and compression for large trajectories
- **Dependency Management** - Optimal task execution order

### Web Interface
- **Real-time Progress** - WebSocket-based status updates
- **Interactive Plots** - Plotly.js visualizations
- **Responsive Design** - Mobile and desktop compatible
- **File Management** - Drag-and-drop trajectory uploads
- **Export Options** - Multiple output formats (HTML, PNG, CSV)

## 🚀 Quick Start

### Docker Deployment (Recommended)

```bash
# Clone the repository
git clone <repository-url>
cd plot_rna

# Start with Docker Compose
docker-compose up -d

# Access the application
open http://localhost
```

### Manual Installation

```bash
# Setup environment
./scripts/setup.sh

# Activate virtual environment
source venv/bin/activate

# Start services
./scripts/start.sh
```

### Development Setup

```bash
# Install dependencies
pip install -r requirements.txt

# Set development environment
export FLASK_ENV=development

# Start development server
python wsgi.py
```

## 📁 Project Structure

```
plot_rna/
├── src/                    # Source code
│   ├── core/              # Core application logic
│   │   ├── app.py         # Main Flask application
│   │   ├── tasks.py       # Celery background tasks
│   │   └── utils.py       # Utility functions
│   ├── analysis/          # Analysis modules
│   │   └── FoldingAnalysis/  # Trajectory analysis tools
│   ├── plotting/          # Visualization functions
│   │   ├── create_plots.py   # Plot generation
│   │   └── energy_3dplot.py  # 3D energy plots
│   └── api/               # API endpoints
├── config/                # Configuration files
│   ├── settings.py        # Application settings
│   ├── logging.py         # Logging configuration
│   └── explanations.json  # Plot explanations
├── templates/             # HTML templates
├── static/                # Static assets (CSS, JS, images)
├── data/                  # Data directory
│   ├── examples/          # Example datasets
│   ├── cache/             # Computation cache
│   └── uploads/           # User uploads
├── deployment/            # Deployment configurations
├── scripts/               # Utility scripts
├── docs/                  # Documentation
└── tests/                 # Test suite
```

## 🔧 Configuration

### Environment Variables

```bash
# Application
FLASK_ENV=production
SECRET_KEY=your-secret-key

# Celery
CELERY_BROKER_URL=redis://localhost:6379/0
CELERY_RESULT_BACKEND=redis://localhost:6379/0

# Performance
MAX_PARALLEL_JOBS=4
MEMORY_LIMIT_GB=8
CACHE_ENABLED=true

# Logging
LOG_LEVEL=INFO
```

### Configuration Files

- `config/settings.py` - Main application configuration
- `config/logging.py` - Logging setup
- `deployment/nginx.conf` - Nginx configuration
- `docker-compose.yml` - Docker services

## 📊 Performance

### Optimization Results
- **RMSD + Landscape**: ~70% time reduction
- **Multiple plots**: ~55-85% time reduction
- **Large trajectories**: ~65% time reduction
- **CPU utilization**: Up to 4x improvement
- **Memory usage**: 30-50% reduction

### Supported File Formats
- **Trajectories**: XTC, DCD, TRR, NetCDF
- **Structures**: PDB, GRO
- **Outputs**: HTML, PNG, CSV, JSON

## 🧪 API Usage

### REST Endpoints

```python
# Submit analysis job
POST /api/submit
{
    "plots": ["RMSD", "LANDSCAPE"],
    "files": {"native": "structure.pdb", "trajectory": "traj.xtc"}
}

# Check job status
GET /api/status/<job_id>

# Download results
GET /api/results/<job_id>
```

### Python Client

```python
from plot_rna_client import PlotRNAClient

client = PlotRNAClient("http://localhost:8000")
job = client.submit_analysis("structure.pdb", "trajectory.xtc", ["RMSD", "LANDSCAPE"])
results = client.wait_for_completion(job.id)
```

## 🐳 Docker Deployment

### Single Container
```bash
docker build -t plot-rna .
docker run -p 8000:8000 plot-rna
```

### Production Stack
```bash
docker-compose -f docker-compose.prod.yml up -d
```

### Kubernetes
```bash
kubectl apply -f deployment/k8s/
```

## 🔬 Scientific Background

Plot RNA implements state-of-the-art algorithms for RNA structural analysis:

- **RMSD Calculations**: Best-fit alignment using Kabsch algorithm
- **Q-factor Analysis**: Native contact fraction (Best-Hummer method)
- **Free Energy Landscapes**: Boltzmann-weighted probability distributions
- **Contact Analysis**: Distance-based residue interaction networks

## 📚 Documentation

- [Installation Guide](docs/installation.md)
- [User Manual](docs/user_guide.md)
- [API Reference](docs/api.md)
- [Performance Tuning](docs/performance.md)
- [Contributing](docs/contributing.md)

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit changes (`git commit -m 'Add amazing feature'`)
4. Push to branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🏆 Citation

If you use Plot RNA in your research, please cite:

```bibtex
@software{plot_rna_2024,
  title={Plot RNA: Advanced RNA Trajectory Analysis Platform},
  author={Plot RNA Development Team},
  year={2024},
  version={2.0.0},
  url={https://github.com/your-repo/plot-rna}
}
```

## 🆘 Support

- 📧 Email: support@plotrna.org
- 🐛 Issues: [GitHub Issues](https://github.com/your-repo/plot-rna/issues)
- 💬 Discussions: [GitHub Discussions](https://github.com/your-repo/plot-rna/discussions)
- 📖 Wiki: [GitHub Wiki](https://github.com/your-repo/plot-rna/wiki)

## 🎯 Roadmap

- [ ] GPU acceleration with CUDA
- [ ] Machine learning-based structure prediction
- [ ] Cloud deployment on AWS/GCP/Azure
- [ ] Integration with molecular visualization tools
- [ ] Advanced statistical analysis modules

---

**Plot RNA** - Empowering RNA research through advanced computational analysis.