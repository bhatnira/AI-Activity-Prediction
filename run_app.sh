#!/bin/bash

# GraphConv Multi-Class App Runner
echo "🧬 Starting GraphConv Multi-Class Classifier..."

# Check if virtual environment exists
if [ -d "graphconv_py310_env" ]; then
    echo "🔄 Activating virtual environment..."
    source graphconv_py310_env/bin/activate
else
    echo "⚠️  Virtual environment not found. Running setup first..."
    ./setup_py310.sh
    source graphconv_py310_env/bin/activate
fi

# Run the app
echo "🚀 Launching Streamlit app..."
streamlit run app_graph_multiclass.py --server.port 8503 --server.address 0.0.0.0
