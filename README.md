# Persistent Cohomological Cycle Matching

## Outline 

This repository is forked from the code for the "Persistent Cohomological Cycle Matching" approach developed in the paper ["Fast Topological Signal Identification and Persistent Cohomological Cycle Matching" (García-Redondo, Monod, and Song 2022)](https://arxiv.org/abs/2209.15446). The original codes were co-written by [Inés García-Redondo](https://sites.google.com/view/ines-garcia-redondo/home) and [Anna Song](https://sites.google.com/view/annasong), performing cycle matching while incorporating the computational advantages given by Ripser [2] and Ripser-image [3].

### Fork Contributions
This fork of the original repository contributes new functionality and makes significant changes to the code structure.
#### Functionality
- generalized to accept arbitrary (pre-computed) distance metrics in the bootstrapping case (see XX.py)
- additional utilities for cycle registration over data bootstraps (see `create_ldm_images.py` and `generate_tag_subindex.py` in `utils_match`)
- integration of cycle matching with statistical permutation testing (see `comp_permtest_dists.py`, `bootstrap_dists.py`, `submit_bsdist/permdist_sbatch.sh`, and `permtest/bootstrap_distances.sh`)
- new modules for computation of Wasserstein variants (see `diagram_distances.py`)
- new visualizations of aggregated statistcs and distributions of bootstrapped cycles (see, e.g., `compare_topostats.py` and `distributional_summaries.py`)
- new toy examples (and corresponding manipulations) for validation and exploration (see `toy_models.py`)
#### Structure
- split cycle-matching implementation into modular structure (see `compute.py`, `extract.py`, and `match.py` in `utils_match`)
- re-factored script-style code into Pythonic functinoal programming style
- optimized for integration with high-performance computing (HPC) environment (including bash scripts with example calls: see `submit_*.sh` scripts in `slurm_bash` below)
- re-organized directory structure

### Structure of the repository

This repository is organised as follows.
#### `modified_ripser`
Garcia-Redondo and Song's modifications to [Ripser](https://github.com/Ripser/ripser/tree/image-persistence-simple) [2] and [Ripser-image](https://github.com/Ripser/ripser/tree/tight-representative-cycles) [3] needed to implement cycle matching. The changes constitute a single inserted line of code, which extracts a lexicographical refinement of indices correposding to simplices of persistence pairs. Garcia-Redondo and Song cede all credit for the files in these folders to the authors in [2] and [3]. 
- `ripser-image-persistence-simple`: Line 474 is altered in `ripser-image-persistence-simple/ripser.cpp`. Forked from: [Ripser-image](https://github.com/Ripser/ripser/tree/tight-representative-cycles)
- `ripser-tight-representative-cycles`: Line 829 is altered in `ripseer-tight-representative-cycles/ripser.cpp`. Forked from: [Ripser](https://github.com/Ripser/ripser/tree/image-persistence-simple)
#### `utils_match`
Prepare Ripser input, perform post-Ripser interval matching, and compute immediate derivative values from match data.
- `compute.py`: Jaccard index and (currently deprecated) Ripser calls
- `create_ldm_images.py`: build and write out Ripser inputs
- `extract.py`: collect and reformat persistence data from Ripser output
- `generate_subindex.py`: generate and bit-encode bootstrapping indices
- `matching.py`: perform cycle-matching operations
- `prevalence.py`: define and read/write prevalence scoring
#### `visualization`
Visualize persistence and cycle-match output data, define and calculate quantities of post-hoc analysis.
- `bootstrap_dists.py`: computes distance between diagrams each containing bootstrap-matched cycles
- `compare_topostats.py`: compares distributions of topological features
- `comp_permtest_dists.py`: compute distribution of distances between data and its null-permuted counterparts
- `diagram_distances.py `: define and compute distance metrics (i.e., Wasserstein variants) between persitence diagrams
- `distributional_summaries.py`: create graphical summaries of topological features
- `prevwt_PD.py`: prevalence-weighted persistence diagram
- `toy_models.py`: module for creation of several toy-model validation spaces (concentric circles, S1 wegde S2 wedge S1, S2 with a diameter, and the torus)
#### `slurm_bash`
SLURM scripting at problem scale, including calls to Ripser.
- `bootstrap_distances.sh`: distributes computation of distances between diagrams with matched cycles
- `calc_dists.sh`: distributes computation of metric (Gram) matrices over dataset(s)
- `calc_dists_w_perms.sh`: distributes computation of metric (Gram) matrices over permuted dataset(s)
- `match_bootstraps.sh`: distributes cycle-matching computation
- `phom_bootstraps.sh`: distrubtes calls to Ripser/Ripser-image
- `prevwt_PDs.sh`: distributes visualization for prevalence-weighted persistence diagrams
- `permtest_distances.sh`: distributes computation of diagram distances between permuted-null persistence diagrams
- `prevalence.sh`: distributes computation of prevalence (i.e., cycle-match quality) information
- `submit_bsdist_sbatch.sh`: SLURM-submits computation of distance between diagrams of matched cycles
- `submit_distmtx_sbatch.sh`: SLURM-submits computation of metric (Gram) matrix
- `submit_permdist_sbatch.sh`: SLURM-submits computation of distances between diagrams and their null-permutation generated counterparts
- `submit_match_sbatch.sh`: SLURM-submits computation of matched cycles
- `submit_ripser_sbatch.sh`: SLURM-submits Ripser (or Ripser-image) call

## Preparations

### About C++

C++ is a general-purpose programming language which has object-oriented, generic, and functional features in addition to facilities for low-level memory manipulation. It is the language chosen for the codes Ripser [2] and Ripser-image [3]. The C++ files for those, with a slight modification needed to implement cycle matching, can be found in the folder `modified ripser`. 

### Compiling the C++ programmes
Before running the code to perform cycle matching in this repository, one needs to compile the C++ files in the `modified ripser` folder. For that
- Install a C++ compiler in your computer. We recommend getting the compiler [GCC](https://gcc.gnu.org/).
	- *For Linux*: the default Ubuntu repositories contain a meta-package named build-essential that contains the GCC compiler and a lot of libraries and other utilities required for compiling software. You only need to run `sudo apt install build-essential` in a terminal to install it.
	- *For Windows*: you can install [Mingw-w64](https://www.mingw-w64.org/) which supports the GCC compiler on Windows systems. You can get this through the installation packages for [MSYS2](https://www.msys2.org/).
	- *For MacOS*: see this [link](https://macappstore.org/gcc/) to install GCC.
- The relevant Makefiles are included in the corresponding folders, so the compilation can be done by running the command line `make` in a terminal opened in the folder. 
- The compiled files should be in the same directory than the python scripts/notebooks in which the cycle matching code is invoked.

### About Python
Python is a high-level, general-purpose programming language. It is the language we use for our code for cycle matching: the 'match' directory defines a small cycle-matching python module.

## ------------------------- UPDATE DEPENDENCIES DISUCSSION NEXT ------------------------- 

### Installing Python libraries
The implementation of cycle matching requires the installation of Python on your computer. 

Additionally, the Python code in `match` requires the installation of the following libraries (follow the corresponding link to find the documentation and installation guidelines):
- [numpy](https://numpy.org/)
- [matplotlib](https://matplotlib.org/stable/index.html)

For the notebook on `applications` one must also install the library
- [scikit-image](https://scikit-image.org/)

We recommend installing Python and these libraries through [Anaconda](https://www.anaconda.com/) and [conda-forge](https://conda-forge.org/).

## Academic use

This code is available and is fully adaptable for individual user customization. If you use the our methods, please cite as the following:

```tex
@misc{garcia-redondo_fast_2022,
	title = {Fast {Topological} {Signal} {Identification} and {Persistent} {Cohomological} {Cycle} {Matching}},
	url = {http://arxiv.org/abs/2209.15446},
	urldate = {2022-10-03},
	publisher = {arXiv},
	author = {García-Redondo, Inés and Monod, Anthea and Song, Anna},
	month = sep,
	year = {2022},
	note = {arXiv:2209.15446 [math, stat]},
	keywords = {Mathematics - Algebraic Topology, Statistics - Machine Learning},
}
```

## References
[1] Reani, Yohai, and Omer Bobrowski. 2021. ‘Cycle Registration in Persistent Homology with Applications in Topological Bootstrap’, January. https://arxiv.org/abs/2101.00698v1.

[2] Bauer, Ulrich. 2021. ‘Ripser: Efficient Computation of Vietoris-Rips Persistence Barcodes’. Journal of Applied and Computational Topology 5 (3): 391–423. https://doi.org/10.1007/s41468-021-00071-5.

[3] Bauer, Ulrich, and Maximilian Schmahl. 2022. ‘Efficient Computation of Image Persistence’. ArXiv:2201.04170 [Cs, Math], January. http://arxiv.org/abs/2201.04170.

[4] Čufar, Matija, and Žiga Virk. 2021. ‘Fast Computation of Persistent Homology Representatives with Involuted Persistent Homology’. ArXiv:2105.03629 [Math], May. http://arxiv.org/abs/2105.03629.
