import argparse
import yaml

from fermipy.gtanalysis import GTAnalysis

def run_analysis(fermipy_config, prefix):
    """
    Run a Fermipy analysis.

    Arguments:
        fermipy_config -- path to Fermipy config file
        prefix -- prefix for output files and output directory
    """

    gta = GTAnalysis(fermipy_config, logging={'verbosity': 3})

    # Prepare the data and perform calculations needed for the analysis
    gta.setup()

    # Fit the normalization and spectral shape parameters
    # of all model components and compute the TS of all sources in the ROI
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

    # Plot TS and residual maps
    model = {
        'SpatialModel': 'PointSource',
        'Index': 2.5,
        'SpectrumType': 'PowerLaw'
        }
    tsmap = gta.tsmap(prefix, model=model, make_plots=True)
    resmap = gta.residmap(prefix, model=model, make_plots=True)

    # Calculate the spectral energy distribution
    sed = gta.sed(fermipy_config['selection']['target'], make_plots=True)

    # Save the results so they can be loaded later
    gta.write_roi(prefix, make_plots=True)


def run_lightcurve(fermipy_config, prefix, **kwargs):
    """
    Run a Fermipy lightcurve analysis.

    Arguments:
        fermipy_config -- path to Fermipy config file
        prefix -- prefix for output files and output directory
    """
    gta = GTAnalysis(fermipy_config, logging={'verbosity': 3})
    gta.load_roi('{}.npy'.format(prefix))
    gta.lightcurve(fermipy_config['selection']['target'],
                   make_plots=True)


if __name__ == "__main__":

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Run a Fermipy analysis.")
    parser.add_argument('config', help="path to pipeline configuration file")
    parser.add_argument('--lightcurve', action='store_true',
                        help="perform a lightcurve analysis")
    args = parser.parse_args()

    with open(args.config, 'r') as config_file:
        pipeline_config = yaml.safe_load(config_file)

    # Load pipeline configuration
    with open(pipeline_config['fermipy_config'], 'r') as config_file:
        fermipy_config = yaml.safe_load(config_file)
    prefix = pipeline_config['prefix']

    # Set the outdir to be the same as the prefix
    if 'fileio' not in fermipy_config:
        fermipy_config['fileio'] = {}
    if 'outdir' not in fermipy_config['fileio']:
        fermipy_config['fileio']['outdir'] = prefix

    run_analysis(fermipy_config, prefix)

    if args.lightcurve:
        run_lightcurve(fermipy_config, prefix **pipeline_config['lightcurve'])
