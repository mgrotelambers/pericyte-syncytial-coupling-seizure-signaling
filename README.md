This repository hosts the code used in the manuscript: *"Syncytial coupling of mid-capillary pericytes underlies seizure-associated electro-metabolic signaling"* doi: https://doi.org/10.64898/2026.03.16.711912.

## Content
The repository is organised into three folders:

- **`python/`** — Code for reproducing the manuscript figures, the electrophysiological and morphological clustering, and the two scripts used for analysing the ABF recordings:
    - `fig_utils.py` — Shared plotting utilities (style configuration, scalebars, ABF/CSV trace plotting, statistical plots, and fluorescence/calcium imaging plots). All figure notebooks import from this module.
    - `cluster_pipeline.py` — Electrophysiological clustering pipeline (IV-curve processing, polynomial fitting, morphology matching, and `KmeansWithNulls` clustering).
    - `morphology_pipeline.py` — Morphological clustering pipeline (feature analysis, clustering with NaN support, silhouette evaluation, feature selection, and UMAP/PCA visualisation).
    - Figure notebooks — Configuration and execution for each manuscript figure. The folder `abf_files` and the file `fig_utils.py` must reside in the same directory as the notebook to be able to run the code.
    - Resting-membrane-potential, 10 mV voltage-step script — These scripts were used to extract the resting potential, membrane resistance ($R_m$) and series resistance ($R_s$) from abf files of current- and voltage-clamp recordings.
- **`r/`** — R scripts used for statistical analysis.
- **`fiji/`** — Fiji/ImageJ macro for single-cell morphology and skeleton analysis of dye-filled pericytes after whole-cell patch clamp.

## Data
The following file types were used: `.abf` files or `.npz` files paired with a corresponding `.json` metadata file, as well as `.csv` files.

## System requirements
- For the Python code (figures, ABF analysis, and clustering):
    - Python 3.12 and the following libraries:
        - Jupyter Notebook
        - NumPy
        - Pandas
        - Matplotlib
        - seaborn
        - SciPy
        - pyabf
        - scikit-learn
        - umap-learn
        - kmeanswithnulls
- For the statistical analysis:
    - R with the packages used in the provided R scripts.
- For the morphological image analysis:
    - Fiji (ImageJ) with the following plugin:
        - MorphoLibJ. Website: https://imagej.net/plugins/morpholibj. Plugin name in the update list: IJPB-plugins.


To run the code installation of the dependencies and downloading the data is necessary. The data that support the findings of this study are available on request from the corresponding author.

## Licenses
- The code in this repository is licensed under the BSD-style license (see LICENSE file).
- The scalebar function used in `fig_utils.py` is adapted from pyABF (Scott W. Harden) and published under the MIT license. Please see the source references within the code.
- The Fiji macro made by M. Kikhia is licensed under BSD-3 (see the header in the macro file). It requires the MorphoLibJ plugin.

## Citation
Syncytial coupling of mid-capillary pericytes underlies seizure-associated electro-metabolic signaling
Mirja grote Lambers, Majed Kikhia, Agustin Liotta, Han Wang, Henrike Planert, Thilo Kalbhenn, Ran Xu, Julia Onken, Thomas Sauvigny, Ulrich-W Thomale, Angela M Kaindl, Martin Holtkamp, Pawel Fidzinski, Matthias Simon, Henrik Alle, Jörg RP Geiger, Christian Madry, Richard Kovács
bioRxiv 2026.03.16.711912; doi: https://doi.org/10.64898/2026.03.16.711912 
