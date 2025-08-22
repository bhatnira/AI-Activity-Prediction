import streamlit as st
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
from rdkit import Chem
from rdkit.Chem import Draw

st.set_page_config(
    page_title="GraphConv Multi-Class Classifier",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

def create_ios_header(title, subtitle=""):
    """Create iOS-style header with clean design"""
    st.markdown(f"""
    <div style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
                padding: 2rem 1rem; 
                border-radius: 15px; 
                margin-bottom: 2rem;
                box-shadow: 0 8px 32px rgba(0,0,0,0.1);">
        <h1 style="color: white; 
                   margin: 0; 
                   font-size: 2.5rem; 
                   font-weight: 700;
                   text-align: center;
                   text-shadow: 2px 2px 4px rgba(0,0,0,0.3);">{title}</h1>
        {f'<p style="color: rgba(255,255,255,0.9); margin: 0.5rem 0 0 0; font-size: 1.1rem; text-align: center;">{subtitle}</p>' if subtitle else ''}
    </div>
    """, unsafe_allow_html=True)

def check_dependencies():
    """Check if all required dependencies are available"""
    missing_deps = []
    
    try:
        import deepchem
        st.success("✅ DeepChem is available")
    except Exception as e:
        missing_deps.append("deepchem")
        st.error(f"❌ DeepChem: {str(e)}")
    
    try:
        import tensorflow
        st.success("✅ TensorFlow is available")
    except Exception as e:
        missing_deps.append("tensorflow")
        st.error(f"❌ TensorFlow: {str(e)}")
    
    try:
        import sklearn
        st.success("✅ Scikit-learn is available")
    except Exception as e:
        missing_deps.append("scikit-learn")
        st.error(f"❌ Scikit-learn: {str(e)}")
    
    return missing_deps

def main():
    create_ios_header("🧬 GraphConv Multi-Class Classifier", "Python 3.10 Compatibility Check")
    
    st.markdown("## 🔧 Environment Diagnostics")
    
    # Python version check
    python_version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
    
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Python Version", python_version)
    with col2:
        st.metric("NumPy Version", np.__version__)
    with col3:
        st.metric("Platform", sys.platform)
    
    st.markdown("---")
    
    # Dependency check
    st.markdown("### 📦 Dependencies Status")
    missing_deps = check_dependencies()
    
    if missing_deps:
        st.markdown("---")
        st.error("🚨 **Compatibility Issues Detected**")
        
        st.markdown("""
        ### 🔧 **Python 3.10 Fix for NumPy/TensorFlow Compatibility**
        
        The error you're seeing is a common issue with Python 3.10 where NumPy and TensorFlow have binary incompatibilities.
        
        #### **Quick Fix (Recommended):**
        ```bash
        # Install compatible versions
        pip3 install numpy==1.24.3 tensorflow==2.13.0 deepchem==2.7.1 --force-reinstall
        
        # Or use the provided requirements file
        pip3 install -r requirements_py310.txt
        ```
        
        #### **Alternative: Create Fresh Environment**
        ```bash
        # Run the setup script
        chmod +x setup_py310.sh
        ./setup_py310.sh
        
        # Or manually:
        python3 -m venv graphconv_env
        source graphconv_env/bin/activate
        pip install -r requirements_py310.txt
        ```
        
        #### **The Issue:**
        - NumPy 1.26+ has breaking changes in dtype structures
        - TensorFlow compiled against older NumPy versions
        - DeepChem depends on specific TensorFlow versions
        
        #### **Compatible Versions for Python 3.10:**
        - NumPy: 1.24.3
        - TensorFlow: 2.13.0
        - DeepChem: 2.7.1
        """)
        
        # Show the requirements file content
        st.markdown("### 📋 **Compatible Requirements (requirements_py310.txt):**")
        try:
            with open("requirements_py310.txt", "r") as f:
                requirements_content = f.read()
            st.code(requirements_content, language="text")
        except FileNotFoundError:
            st.warning("requirements_py310.txt not found in current directory")
        
        # Provide download for requirements
        requirements_content = """# Python 3.10 Compatible Requirements
numpy==1.24.3
pandas==2.0.3
matplotlib==3.7.2
seaborn==0.12.2
streamlit==1.25.0
rdkit==2023.3.2
scikit-learn==1.3.0
tensorflow==2.13.0
protobuf==3.20.3
deepchem==2.7.1
Pillow==10.0.0
openpyxl==3.1.2
xlrd==2.0.1
h5py==3.9.0
joblib==1.3.2"""
        
        st.download_button(
            label="📥 Download Compatible Requirements",
            data=requirements_content,
            file_name="requirements_py310_fixed.txt",
            mime="text/plain"
        )
        
    else:
        st.success("🎉 **All dependencies are working correctly!**")
        st.info("You can now use the full GraphConv Multi-Class Classifier application.")
        
        if st.button("🚀 Launch Full Application", type="primary"):
            st.info("Restart the app with the main application file to access all features.")
    
    st.markdown("---")
    st.markdown("### 📊 **Simple SMILES Validator (Always Available)**")
    
    smiles_input = st.text_input("🧬 Test SMILES validation:", placeholder="CCO")
    
    if smiles_input:
        try:
            mol = Chem.MolFromSmiles(smiles_input)
            if mol:
                st.success(f"✅ Valid SMILES: {smiles_input}")
                st.info(f"📊 Molecular formula: {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
                st.info(f"⚖️ Molecular weight: {Chem.rdMolDescriptors.CalcExactMolWt(mol):.2f}")
                
                # Draw molecule
                img = Draw.MolToImage(mol, size=(300, 300))
                st.image(img, caption=f"Structure: {smiles_input}")
            else:
                st.error("❌ Invalid SMILES string")
        except Exception as e:
            st.error(f"❌ Error processing SMILES: {str(e)}")

if __name__ == "__main__":
    main()
