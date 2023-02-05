# Internship-IBIVU

This repository is a collection of the scripts that I made during my internship at IBIVU. Included in this repository are python workflows for umbrella sampling and metadynamics simulations of amyloid fibril elongation. For umbrella sampling, the files can be found in `UmbrellaSampling/noSASA` and `UmbrellaSampling/SASA`, which should be able to reproduce the umbrella sampling simulations. Similarly, the files required for metadynamics can be found in `MetaDynamics/noSASA` and `MetaDynamics/SASA`. Keep in mind that these workflows have not been tested on a different system, so some tinkering may be necessary to get it to work on compute clusters other than the VU BAZIS cluster.

Finally, there is a folder called `archive`, which is purely just meant as that, an archive. Many of these scripts have hardcoded variables or need a specific working directory. For me personally this serves as 1) an archive in case I need anything later, and 2) a way to record my current coding skills which I can compare myself to later. Having said that, feel free to look around and use any scripts or code snippets to your liking!

## Setting up
In order to run the workflow scripts there are few prerequisites that need to be installed first. These are:
- [Amber20](http://ambermd.org/)
- [PLUMED 2.8.1](https://github.com/plumed/plumed2/releases/tag/v2.8.1)
- [WHAM 2.0.11](http://membrane.urmc.rochester.edu/?page_id=126)
- [Conda 4.12.0](https://docs.conda.io/en/latest/miniconda.html) with a specific environment (see below)

### Amber20
The Amber20 simulation suite should be compiled and installed using GPU support. All my scripts were used with a NVIDIA RTX 2070 and [CUDA Toolkit version 10.1.243](https://developer.nvidia.com/cuda-10.1-download-archive-base). Compilation was done using [GCC 8.3.0](https://gcc.gnu.org/gcc-8/).

### PLUMED 2.8.1
PLUMED 2.8.1 needs to be compiled and installed including the `sasa` module. When compiling, make sure to include `--install-modules=sasa` directive. After installation, make sure to inlcude the PLUMED binaries in your `$PATH` variable and dynamically link the `$PLUMED_KERNEL` environment variable to `<your install prefix>/lib/libplumedKernel.so`. I recommend placing these commands in your `.bash_profile` file.

### WHAM
This program runs the weighted histogram analysis method. To fully reproduce my results for umbrella sampling, compile the program using the kJ/mol units as described in [their documentation](http://membrane.urmc.rochester.edu/sites/default/files/wham/doc.pdf)

### Conda environment
Make sure conda 4.12.0 is installed first. Then clone the repository to a local folder, `cd` there and install the required packages using
```bash
conda env -n AnalysisTools -f requirements.yaml
```
Then activate it using
```bash
conda activate AnalysisTools
```

## Running the workflow
The entire workflow is run seperately for both simulation schemes. The `noSASA` scheme is the regular generalized borne implicit water model, while the `SASA` scheme includes a temperature dependent potential. 

### Umbrella Sampling
Run the following two notebooks in order:
1. `run_US_all_<scheme>.ipynb` - This prepares all the input files, parameters, PLUMED files and structures for Umbrella Sampling. If running this on the VU BAZIS cluster, it can be submitted directly to the job scheduler. If not, please modify `batch_window.sh` to suit your needs.
2. `WHAM_all_<scheme>.ipynb` - This runs the WHAM program to produce the free energy files, and also obtains some of the plots.

### Metadynamics
Run the following notebooks in order:
1. `run_metad_<scheme>.ipynb` - This prepares all the input files, parameters, PLUMED files and the structure for metadynamics. If running this on the VU BAZIS cluster, it can be submitted directly to the job scheduler. If not, please modify `batch_window-noscratch.sh` or `batch_window.sh` to suit your needs.
2. `FESfigures.ipynb` - This runs the postprocessing of the gaussian hills and produces the free energy landscape data. It also produces the plots. For one of the simulations (300K SASA), it also produces a convergence figure. Note that this requires all simulations to be finished.
