# Syncytial coupling of mid-capillary pericytes underlies seizure-associated electro-metabolic signaling
This repository contains the Python and R analysis scripts used to generate the figures and statistical results for the publication "Syncytial coupling of mid-capillary pericytes underlies seizure-associated electro-metabolic signaling". It includes pipelines for morphological clustering, electrophysiology data visualization, and calcium imaging analysis.

Raw data (ABF files, CSVs, imaging data) are available on Figshare: [Figshare DOI link here]

Repository Structure
pericyte-syncytial-coupling-seizure-signaling/
├── python/
│   ├── morphology_pipeline.py      # Morphological clustering pipeline
│   ├── cluster_pipeline.py         # UMAP/PCA clustering analysis
│   ├── fig_utils.py                # Figure utility functions
│   └── abf_export.py               # ABF electrophysiology data export
├── r/
│   └── [your R scripts]            # Statistical analysis scripts
├── requirements.txt                # Python dependencies
├── LICENSE
└── README.md

Requirements
**Python**
Python 3.x
numpy
pandas
matplotlib
pyabf
umap-learn
scikit-learn

**R**
R 4.x
nparcomp

**Installation**
bashgit clone https://github.com/mgrotelambers/pericyte-syncytial-coupling-seizure-signaling.git
cd pericyte-syncytial-coupling-seizure-signaling
pip install -r requirements.txt

Data Availability
Raw and processed data are deposited on Figshare and can be downloaded separately: [Insert Figshare DOI and link here]

Place downloaded data files in a data/ folder in the root of the repository before running the scripts.

Usage
Download the data from Figshare, then run the analysis scripts from the python/ or r/ directories. Refer to the comments within each script for details on inputs and outputs.

**Authors**
Mirja grote Lambers, Charité – Universitätsmedizin Berlin, Department of Neurophysiology
Majed Kikhia, Charité – Universitätsmedizin Berlin, Department of Experimental Neurology  


**Citation**
If you use this code, please cite:
M. grote Lambers et al. (2026). Syncytial coupling of mid-capillary pericytes underlies seizure-associated electro-metabolic signaling. [Journal name]. [DOI]


**License**
This project is licensed under the BSD 3-Clause License. See the LICENSE file for details.
