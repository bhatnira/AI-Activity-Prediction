#!/bin/bash

echo "🔧 Fixing NumPy/TensorFlow compatibility for Python 3.10..."
echo "📋 This will install compatible package versions"

# Show current versions
echo "Current package versions:"
pip3 list | grep -E "(numpy|tensorflow|deepchem)" || echo "Some packages not found"

echo ""
echo "🚀 Installing compatible versions..."

# Install compatible NumPy first
echo "📦 Installing NumPy 1.24.3..."
pip3 install numpy==1.24.3 --force-reinstall --no-deps

# Install compatible TensorFlow
echo "📦 Installing TensorFlow 2.13.0..."
pip3 install tensorflow==2.13.0

# Install compatible DeepChem
echo "📦 Installing DeepChem 2.7.1..."
pip3 install deepchem==2.7.1

# Install other required packages
echo "📦 Installing other dependencies..."
pip3 install pandas==2.0.3 matplotlib==3.7.2 seaborn==0.12.2
pip3 install rdkit==2023.3.2 scikit-learn==1.3.0
pip3 install streamlit==1.25.0 openpyxl==3.1.2

echo ""
echo "✅ Installation complete!"
echo ""
echo "🧪 Testing imports..."

# Test critical imports
python3 -c "
import numpy as np
print(f'✅ NumPy {np.__version__} imported successfully')

try:
    import tensorflow as tf
    print(f'✅ TensorFlow {tf.__version__} imported successfully')
except Exception as e:
    print(f'❌ TensorFlow import failed: {e}')

try:
    import deepchem as dc
    print(f'✅ DeepChem imported successfully')
except Exception as e:
    print(f'❌ DeepChem import failed: {e}')

try:
    from rdkit import Chem
    print(f'✅ RDKit imported successfully')
except Exception as e:
    print(f'❌ RDKit import failed: {e}')
"

echo ""
echo "🎉 Python 3.10 compatibility fix complete!"
echo "🚀 You can now run: streamlit run app_graph_multiclass.py"
