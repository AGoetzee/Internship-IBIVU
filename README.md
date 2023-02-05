# Internship-IBIVU

This repository is a collection of the scripts that I made during my internship at IBIVU. It is purely meant as an archive, many of these scripts have hardcoded variables or need a specific working directory. For me personally this repo serves as 1) an archive in case I need anything later, and 2) a way to record my current coding skills which I can compare myself to later.
Having said that, feel free to look around and use any scripts or code snippets to your liking!
Included here is a python notebook `run_workflow.ipynb`, which should be able to reproduce the umbrella sampling simulations. It has not been tested on a different system, so some tinkering may be necessary to get it to work.

## Setting up
In order to run the workflow scripts there are few prerequisites that need to be installed first. These are:
- [Amber20](http://ambermd.org/)
- [PLUMED 2.8.1](https://github.com/plumed/plumed2/releases/tag/v2.8.1)
- [WHAM 2.0.11](http://membrane.urmc.rochester.edu/?page_id=126)
- [Conda 4.12.0](https://docs.conda.io/en/latest/miniconda.html) with a specific environment (see below)

### Amber20
The Amber20 simulation suite should be compiled and installed using GPU support. All my scripts were used with a GTX Titan and [CUDA 10.1.243](https://developer.nvidia.com/cuda-10.1-download-archive-base). Compilation was done using [GCC 8.3.0](https://gcc.gnu.org/gcc-8/).

### PLUMED 2.8.1
PLUMED 2.8.1 needs to be compiled and installed including the `sasa` module. When compiling, make sure to include `--install-modules=sasa` directive. After installation, make sure to inlcude the PLUMED binaries in your `$PATH` variable and dynamically link the `$PLUMED_KERNEL` environment variable to `<your install prefix>/lib/libplumedKernel.so`. I recommend placing these commands in your `.bash_profile` file.

### WHAM
This program runs the weighted histogram analysis method. To fully reproduce my results, compile the program using the kJ/mol units as described in [their documentation](http://membrane.urmc.rochester.edu/sites/default/files/wham/doc.pdf)

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
The entire workflow is run seperately for both simulation schemes. The `noSASA` scheme is the regular generalized borne implicit water model, while the `SASA` scheme includes a temperature dependent potential. Run the following two notebooks in order:
1. `run_US_all_<scheme>.ipynb` - This prepares all the input files, parameters, PLUMED files and structures for Umbrella Sampling. If on running this on the VU BAZIS cluster, it can be submitted directly to the job scheduler. If not, please modify `batch_window.sh` to suit your needs.
2. `WHAM_all_<scheme>.ipynb` - This runs the WHAM program to produce the free energy files, and also obtains some of the plots.
