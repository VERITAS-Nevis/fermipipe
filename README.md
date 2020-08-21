# fermipipe
Fermi analysis pipeline built on FermiPy. Provides:
- An installation script suitable for the Columbia/Barnard cluster at Nevis
- A turnkey configuration file and analysis script with default settings suitable for most typical analyses
- Functionality to calculate light curves more robustly

## [Documentation](https://twiki.nevis.columbia.edu/twiki/bin/view/Veritas/FermiPyAnalysis)

The complete instructions for setup and use, as well as links to more information about Fermi analysis in general and this pipeline in particular, can be found on the [Nevis wiki](https://twiki.nevis.columbia.edu/twiki/bin/view/Veritas/FermiPyAnalysis).

## Quickstart

### Installation

Note: Only steps 2 and 3 are necessary on the Nevis cluster, as the pipeline has already been installed.

1. `git clone https://github.com/aribrill/fermipipe.git`
2. `cat fermipipe/bash_setup.txt >> ~/.myprofile`
3. `source ~/.bashrc`
4. `$FERMIPIPE/setup_pipeline.sh -p <CONDA_PATH>`

### Running an analysis

To perform an analysis, you will have to set a minimum of three parameters: the start time, the stop time, and the name of the target to analyze (if your target is not a known Fermi source, you can define the target using coordinates instead). The analysis start and stop times must be provided in Fermi mission elapsed time (MET).

Set up your environment and activate the conda environment: `fermisetup`

Run a base analysis: `fermianalysis` (alias for `nohup python $FERMIPIPE/run_analysis.py $FERMI_ANALYSIS_DIR/pipeline_config.yml &`)

Once the base analysis is complete, run the analysis for the first section of the light curve:

`python $FERMIPIPE/run_analysis.py $FERMI_ANALYSIS_DIR/pipeline_config.yml --lightcurve --section 0`

After running the analyses for all light curve sections, combine the output files:

`python $FERMIPIPE/combine_lightcurve_results.py $FERMI_ANALYSIS_DIR/pipeline_config.yml`


