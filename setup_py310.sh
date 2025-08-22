#!/bin/bash

# GraphConv Multi-Class App - Python 3.10 Setup Script
echo "🧬 Setting up GraphConv Multi-Class Classifier for Python 3.10..."

# Check Python version
python_version=$(python3 --version 2>&1 | grep -oE '[0-9]+\.[0-9]+')
echo "📋 Detected Python version: $python_version"

if [[ "$python_version" != "3.10" ]]; then
    echo "⚠️  Warning: This app is optimized for Python 3.10"
    echo "   Current version: $python_version"
    echo "   Consider using Python 3.10 for best compatibility"
fi

# Create virtual environment
echo "🔧 Creating virtual environment..."
python3 -m venv graphconv_py310_env

# Activate virtual environment
echo "🔄 Activating virtual environment..."
source graphconv_py310_env/bin/activate

# Upgrade pip
echo "⬆️  Upgrading pip..."
pip install --upgrade pip

# Install compatible packages
echo "📦 Installing Python 3.10 compatible packages..."
pip install -r requirements_py310.txt

echo "✅ Setup complete!"
echo ""
echo "🚀 To run the app:"
echo "   source graphconv_py310_env/bin/activate"
echo "   streamlit run app_graph_multiclass.py"
echo ""
echo "🔗 Or run directly:"
echo "   ./run_app.sh"
