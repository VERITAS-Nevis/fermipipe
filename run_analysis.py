import argparse
import yaml

from fermipy.gtanalysis import GTAnalysis

def run_analysis(fermipy_config, prefix):

    # Load FermiPy config file
    with open(fermipy_config, 'r') as config_file:
        config = yaml.safe_load(config_file)
    source_name = config['selection']['target']

    # Set the outdir to be the same as the prefix
    if 'fileio' not in config:
        config['fileio'] = {}
    if 'outdir' not in config['fileio']:
        config['fileio']['outdir'] = prefix

    gta = GTAnalysis(config, logging={'verbosity': 3})

    # Prepare the data and perform calculations needed for the analysis
    gta.setup()

    # Fit the normalization and spectral shape parameters of all model components
    # and compute the TS of all sources in the ROI
    gta.optimize()

    # Remove undetected sources (TS < 1) to simplify the model
    gta.delete_sources(minmax_ts=[None, 1], exclude=['galdiff', 'isodiff'])

    # Free the normalization of high TS sources and those close to the ROI center
    gta.free_sources(minmax_ts=[10, None], pars='norm')
    gta.free_sources(distance=3.0, pars='norm')

    # Free all parameters of the galactic and isotropic diffuse components
    gta.free_source('galdiff')
    gta.free_source('isodiff')

    # Perform the model fit, and fix all sources again when done
    gta.fit()
    gta.free_sources(free=False)

    # Save the results so they can be loaded later
    gta.write_roi(prefix, make_plots=True)

    # Plot TS and residual maps
    model = {
        'SpatialModel': 'PointSource',
        'Index': 2.5,
        'SpectrumType': 'PowerLaw'
        }
    tsmap = gta.tsmap(prefix, model=model, make_plots=True)
    resmap = gta.residmap(prefix, model=model, make_plots=True)

    # Calculate the spectral energy distribution
    sed = gta.sed(source_name, make_plots=True)

if __name__ == "__main__":

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Run a Fermipy analysis.")
    parser.add_argument('config', help="path to pipeline configuration file")
    args = parser.parse_args()

    with open(args.config, 'r') as config_file:
        pipeline_config = yaml.safe_load(config_file)

    run_analysis(**pipeline_config)
