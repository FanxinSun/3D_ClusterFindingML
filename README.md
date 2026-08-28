# 3D_ClusterFindingML

> **New here?** Start with **[FIRST_RUN.md](FIRST_RUN.md)** — environment setup, a
> two-minute demo of the whole chain, and the full production run. The notes below
> describe the original container-based approach, which the current detached
> pipeline (v6.1) no longer requires.

Use Singularity for using official sPHENIX geometry and conditioned DB: 
https://github.com/sPHENIX-Collaboration/Singularity

Without a bnl account, the only way is to use a 2022 ana based on gcc=8.3, as the README.md says in the above repo.
Some open-sourced CDB are collected targeting a combination of the ideal geometry for simulation usage.

# Simulation with frozen ana.331 combining CDB_offlines
The geometry (ideal one) is frozen and open-sourced in 2022, while 177 open-sourced Conditional DBs are downloaded and can be loaded to the ideal geometry, together with beam and magnetic field configurations.
The CDB mapping is internal and cannot be accessed without BNL account, but a random combination is considered sufficient for ML training.

# Add offline CDBs as a combination of offline cvmfs

### Geant4 Particle Types                
 
  Common Hadrons:
  
  •  pi-  /  pi+  /  pi0  (Pions)
  •  proton  /  anti_proton 
  •  neutron  /  anti_neutron 
  •  kaon+  /  kaon- 
  
  Leptons & Photons:
  
  •  e-  /  e+  (Electron / Positron)
  •  mu-  /  mu+  (Muon / Anti-muon)
  •  gamma  (Photon)
  
  Heavy Ions / Nuclei:
  
  •  alpha 
  •  deuteron 
  •  triton 
  
  Debugging:
  
  •  geantino

A Machine Learning project focused on identifying and reconstructing 3D clusters within the **sPHENIX TPC** detector using real data. This project aims to bridge the gap between traditional reconstruction and modern deep learning by creating a manually labeled Ground Truth dataset.

## 🔗 Project Resources
*   **Shared Project Documentation:** [Google Docs Link - Placeholder]

## 🔬 Dataset Properties & Geometry (Sanity Checks)
Through empirical analysis of the `ntp_hit` ROOT tree, the following sPHENIX detector properties have been established for the Machine Learning pipeline:

### 1. Detector Subsystem Mapping
The dataset contains hits from multiple sub-detectors. They are mapped by `layer`:
*   **MVTX (Inner Silicon Pixels):** Layers 0, 1, 2
*   **INTT (Silicon Strips):** Layers 3, 4, 5, 6
*   **TPC (Gas Drift Volume):** Layers 7 through 54
*   **TPOT (Micromegas):** Layers 55, 56

### 2. TPC Coordinate Relationships
For the **TPC (Layers 7-54)**, the empirical formulas used to generate the spatial coordinates are:
*   **Radius (`r`):** Strictly fixed per layer (Standard Deviation = 0.0 cm). The `layer` integer acts as a perfect 1:1 lookup for the physical radius.
*   **X & Y:** Calculated from the continuous physical radius `r` (not the discrete `layer` index).
    *   `x = r * cos(phi)`
    *   `y = r * sin(phi)`
*   **Z:** Calculated linearly from drift time (`tbin`), dependent on the `zelem` branch (which indicates the side of the central membrane):
    *   If `zelem == 0`: `z = 0.42341 * tbin + 0.28000`
    *   If `zelem == 1`: `z = -0.42341 * tbin + 102.60500`

### 3. Bulletproof Data Cleaning Cuts
To isolate pure, physical TPC hits for ML training, the following cuts must be applied sequentially to the `ntp_hit` DataFrame:
1.  `layer >= 7 & layer <= 54`: Isolates the TPC.
2.  `hitID > 0`: Removes background/noise hits.
3.  `adc > 0`: Removes zero-energy ghost hits/artifacts.
4.  `dropna()`: Removes NaNs (which mostly populate the outer TPOT layers).
5.  `zelem.isin([0, 1])`: Enforces valid TPC geometry bounds.

## 🚀 ML Training Pipeline
The project follows a three-stage strategy to ensure robust cluster identification:

1.  **Numerical Classification:** Utilize traditional ML models (Random Forest, XGBoost, or MLP) on numerical hit data to perform initial categorization and feature importance analysis.
2.  **Geometric Deep Learning:** Implement Graph Neural Networks (GNNs) to exploit the spatial and topological relationships of hits, predicting clusters and track segments directly from the hit-graph.
3.  **Ensemble & Validation:** Merge and cross-verify results from both the numerical and GNN models to improve purity and efficiency, using the manually labeled Ground Truth for final performance metrics.

---

## 💻 Environment Setup

### Mac (Apple Silicon M-series)
1.  **Xcode Command Line Tools:** `xcode-select --install`
2.  **Homebrew:** Install [Homebrew](https://brew.sh/) if not already present.
3.  **ROOT CERN:** `brew install root`
4.  **Python Environment:**
    ```bash
    conda create -n cluster_ml python=3.10
    conda activate cluster_ml
    pip install numpy matplotlib pandas scikit-learn torch torch-geometric
    ```

### WSL (Windows Subsystem for Linux - Ubuntu)
1.  **System Dependencies:**
    ```bash
    sudo apt-get update
    sudo apt-get install dpkg-dev cmake g++ gcc binutils libx11-dev libxpm-dev libxft-dev libxext-dev python3 libssl-dev
    ```
2.  **ROOT CERN:** It is recommended to install ROOT via `conda` for easiest dependency management in WSL:
    ```bash
    conda create -n cluster_ml -c conda-forge root python=3.10
    conda activate cluster_ml
    pip install numpy matplotlib pandas scikit-learn torch torch-geometric
    ```

---

## 📊 Data Visualization
To inspect the data and verify labeling:
1.  **Extract Data:**
    ```bash
    root -l -b -q Visualizing/src/GenerateEvent74CSVs.C
    ```
2.  **Run Plotter:**
    ```bash
    python Visualizing/src/plot_event74_hits.py
    ```

---

## 📂 Project Structure
*   `Visualizing/`: Scripts for data extraction and 3D hit inspection.
*   `Articles/`: Reference papers for TPC reconstruction and ML methodologies.
*   `*.ipynb`: Research and training notebooks.
*   `memory/`: (Private) Project development notes.