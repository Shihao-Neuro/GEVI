# GEVI Voltage Imaging Analysis

**Associated with the manuscript: _Bright and Photostable Voltage Sensors Derived from mBaoJin_**

This repository provides code for the analysis of voltage imaging data acquired from 1-photon cultured neurons and in vivo experiments using Genetically Encoded Voltage Indicators (GEVIs). The provided scripts cover the entire workflow from motion correction, ROI segmentation, trace extraction, to spike waveform analysis.

---

## 📁 Repository Structure

```

GEVI
├── data/                   # Placeholder for data storage
├── LICENSE
├── README.md
└── src/
    ├── 1p_culture/         # Python scripts for 1-photon culture cell data analysis
    └── 1p_in_vivo/         # MATLAB scripts for 1-photon in vivo data analysis

```

---

## 🧪 1-Photon Culture Imaging Analysis (`src/1p_culture/`)

These Python scripts implement a complete analysis pipeline based on [VolPy](https://github.com/flatironinstitute/CaImAn/tree/master/caiman/source_extraction/volpy). The pipeline is tailored for GEVI imaging in cultured neurons.

Runtime: < 60 minutes per 60 sessions on a standard desktop.

### 📜 Scripts

- **`run_volpy_pipeline.py`**  
  Preprocessing and signal extraction pipeline:  
  - Motion correction  
  - ROI segmentation (manual segmentation based)  
  - VolPy spike and subthreshold activity inference  

- **`analyze_volpy_results.py`**  
  Post-analysis of extracted traces:  
  - SNR computation  
  - Spike waveform extraction  
  - Segment-wise spike comparison  
  - .mat export for downstream statistics or visualization

---

## 🧬 1-Photon In Vivo Imaging Analysis (`src/1p_in_vivo/`)

These MATLAB scripts were adapted from [WKLabION/Confocal-Light-Field-Microscopy](https://github.com/WKLabION/Confocal-Light-Field-Microscopy) and modified for longitudinal in vivo GEVI recordings.

Runtime: ~10–15 minutes per image on a standard desktop.

### 📜 Scripts

- **`Folder_generator.m`**  
  Creates structured directories (`001`, `002`, ...) and organizes `.nd2` image sequences accordingly.

- **`motion_pca_correction_Electra_batch.m`**  
  Batch processing of `.nd2` files for motion correction, global background subtraction, and PCA-based dimensionality reduction. Edited based on Bai L et al. Nat Methods 2024.  
  **Note:** Requires the bioformat package for function calls and loading nd2 images.

- **`get_spikeinfo_bkg_1fov_Electra_batch.m`**  
  Subtracts PCs, extracts ROI-based traces, detects spikes, and stores all results in a `Spike_Info` folder. Edited based on Bai L et al. Nat Methods 2024.  
  **Note:** Also requires Bio-Formats for image I/O.

- **`rename_spike_info.m`**  
  Optional utility to prepend identifiers to spike info files for batch traceability.

- **`Spike_info_extraction.m`**  
  Aggregates and exports spike trace statistics (SNR, spike times, etc.) across all ROIs in the working directory.  

- **`subfunctions/`**  
  Auxiliary I/O utilities:
  - `imstackread.m`, `imstackwrite.m`: ND image reading/writing
  - `spike_denoise.m`: spike train denoising utility

---

## 📊 Demo Dataset

A demo dataset for two intermitent recordings of ElectraOFF were provided for benchmarking and testing.  
**[Download link – https://westlakeu-my.sharepoint.com/:f:/g/personal/kiryl_piatkevich_westlake_edu_cn/EqVKdQ7gf0hHjB2rmoKNLuwB4L_q74VChhpzV1l1yQ6xCQ?e=TzmtEp]**

---

## 📦 Dependencies

### Python
- Python ≥ 3.8  
- [CaImAn (VolPy module)](https://github.com/flatironinstitute/CaImAn)
- NumPy, SciPy, Matplotlib, h5py

### MATLAB
- MATLAB ≥ R2020
- bioformat package

---

## 📄 License

This repository is distributed under the MIT License. See [LICENSE](LICENSE) for details.
