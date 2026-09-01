# Pericyte Paper 2026

This repository hosts the code used in the manuscript: *"[MANUSCRIPT TITLE]"* [DOI / journal link].

## Content

The repository is organised into three folders by language/tool:

- **`python/`** — Code for reproducing the manuscript figures, the electrophysiological and morphological clustering, and the two scripts used for analysing the ABF recordings:
    - `fig_utils.py` — Shared plotting utilities (style configuration, scalebars, ABF/CSV trace plotting, statistical plots, and fluorescence/calcium imaging plots). All figure notebooks import from this module.
    - `cluster_pipeline.py` — Electrophysiological clustering pipeline (IV-curve processing, polynomial fitting, morphology matching, and `KmeansWithNulls` clustering).
    - `morphology_pipeline.py` — Morphological clustering pipeline (feature analysis, clustering with NaN support, silhouette evaluation, feature selection, and UMAP/PCA visualisation).
    - Figure notebooks (e.g. `Figure_2.ipynb` … `Figure_5.ipynb`) — Configuration and execution for each manuscript figure. Note that the folder `abf_files` and the file `fig_utils.py` must reside in the same directory as the notebook.
    - Resting-membrane-potential and 10 mV voltage-step scripts — Standalone scripts (pyabf) for extracting resting potential, membrane resistance ($R_m$), and series resistance ($R_s$) from current- and voltage-clamp recordings.
- **`r/`** — R scripts for the statistical analysis.
- **`fiji/`** — Fiji/ImageJ macro by Kikhia et al. (a coauthor of the paper) for single-cell morphology and skeleton analysis of dye-filled pericytes after whole-cell patch clamp.

## Data

To run this code, please download the associated data available here: [DATA REPOSITORY URL]. The data need to be stored in the paths expected by the notebooks and scripts; in particular, the `abf_files` folder and `fig_utils.py` must sit in the same directory as the notebook that uses them.

The recordings are provided as `.abf` files, or as `.npz`/`.csv` files each with a companion `.json` metadata file. The loader (`load_abf_or_csv` in `fig_utils.py`) auto-detects the format by extension and returns a `pyabf.ABF`-compatible object, so the notebooks run unchanged regardless of which format is provided.

## System requirements

- For the Python code (figures, ABF analysis, and clustering):
    - Python 3.9 with the following libraries:
        - Jupyter Notebook
        - NumPy
        - Pandas
        - Matplotlib
        - seaborn
        - SciPy
        - pyabf
        - scikit-learn
        - umap-learn
    - The clustering pipelines additionally depend on the `KmeansWithNulls` package, which provides a K-means implementation supporting missing values (imported as `from KmeansWithNulls import KmeansWithNulls`).
- For the statistical analysis:
    - R with the packages used in the provided R scripts.
- For the morphological image analysis:
    - Fiji (ImageJ) with the following plugin:
        - MorphoLibJ. Website: https://imagej.net/plugins/morpholibj. Plugin name in the update list: IJPB-plugins.

The `arial` font is used for all figures; if it is not installed, Matplotlib will fall back to a default font (this affects appearance only, not results).

## Hardware

The code files were run on a standard laptop; no GPU is required. Runtimes below are indicative for such a machine.

## Installation guide

1. Download the source code together with the associated data.
2. Create a virtual environment with Python and the libraries mentioned above.
3. Install R and the packages used in the R scripts.
4. Install Fiji and add MorphoLibJ from the update list.

### Installation time

- Downloading the code and data requires several minutes, depending on internet speed.
- Installing all necessary software and packages might require up to an hour. However, the code can be run immediately if all required software is available.

## Instructions for analysis

After installing the dependencies and downloading the data, the following analyses can be run.

### Python (`python/`)

- **Figure notebooks (Figures 1–5):**
    - Open the notebook for the figure of interest.
    - Ensure the `abf_files` folder and `fig_utils.py` are in the same directory as the notebook.
    - Run the notebook. Figures are displayed inline.
    - To save the figures to disk, uncomment `enable_saving()` at the top of the notebook; all figures with a defined `save_path` are then written out.

- **Resting membrane potential:**
    - Resting membrane potential is determined from I = 0 current-clamp recordings. For each `.abf` file, data are pooled across all sweeps and the median membrane potential is computed per channel (median is preferred over mean to reduce sensitivity to transient noise artifacts). The script detects single- and dual-channel recordings automatically and exports the results to a spreadsheet.
    - Set the input directory, then run the script.

- **10 mV voltage-step analysis:**
    - Membrane resistance ($R_m$) and series resistance ($R_s$) are extracted from voltage-clamp recordings (19 sweeps, 300 ms, 10 kHz, holding potential −100 mV). A 50 ms pre-step epoch at 0 mV serves as baseline, followed by a 50 ms +10 mV step. Baseline-corrected peak-transient and steady-state currents are averaged across sweeps and used to derive $R_s$ and $R_m$ via Ohm's law:

        $$R_s = \frac{10}{(I_{peak} - I_{ss}) + I_{ss}} \times 1000 \qquad R_{m+s} = \frac{10}{I_{ss}} \times 1000 \qquad R_m = R_{m+s} - R_s$$

        where current is in pA and resistance in MΩ.
    - Set the input directory, then run the script.

- **Electrophysiological pericyte cluster analysis:**
    - Pericytes are classified by their electrophysiological profile using IV-curve processing, polynomial fitting, morphology matching, and `KmeansWithNulls` clustering. Pipeline logic lives in `cluster_pipeline.py`; the notebook handles configuration and execution.

- **Morphological pericyte cluster analysis (OHSC):**
    - OHSC pericytes are clustered by morphological features using `KmeansWithNulls`. Pipeline logic lives in `morphology_pipeline.py`; the notebook handles configuration and execution. UMAP reproducibility is enforced via a fixed random seed and `n_jobs=1`.

### R (`r/`)

- The R scripts contain the statistical analysis reported in the manuscript. Open the relevant script and run it, adjusting input paths as needed.

### Fiji (`fiji/`)

- **Single-cell morphology (Fiji macro):**
    - Open Fiji, then open the macro in the `fiji/` folder.
    - Set the `output_folder` at the top of the macro.
    - Open an image and run the macro. The macro is interactive: it prompts the user to mark where the pipette ends in Z, to draw ROIs around the pipette and its shadow, to enter the maximum intensity of the pipette shadow, and to confirm or adjust the automated threshold.
    - Output: 6 TIFF files (MIP of the original image, 2D and 3D morphology, and three skeleton-analysis images), 3 CSV files (morphology parameters, skeleton parameters, and the manually entered values), and 2 ROIs (pipette and pipette shadow, saved for reference and reproducibility).

## Licenses

- Unless stated otherwise, the code in this repository is licensed under the BSD-style license found in the LICENSE file in the root directory of this source tree.
- The scalebar functions in `fig_utils.py` are adapted from pyABF (Scott W. Harden) and published under the MIT license. See the source references within the code.
- The Fiji macro is by Kikhia et al. and licensed under BSD-3 (see the header in the macro file). It requires the MorphoLibJ plugin.

## Citation

[AUTHORS]. *[MANUSCRIPT TITLE].* [Journal] (2026). [DOI]

## Author

[YOUR NAME], [AFFILIATION].
Code adapted from other authors is indicated within the code files.
