import argparse
import glob
import os

from fermipy.gtanalysis import GTAnalysis
import numpy as np
import yaml

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


def get_times(fermipy_config):
    """
    Calculate sequence of time bins as defined in Fermipy.

    https://github.com/fermiPy/fermipy/blob/master/fermipy/lightcurve.py
    """

    tmin = fermipy_config['selection']['tmin']
    tmax = fermipy_config['selection']['tmax']
    lc_config = fermipy_config['lightcurve']

    # Make array of time values in MET
    if lc_config.get('time_bins'):
        times = np.array(lc_config['time_bins'])
    elif lc_config.get('nbins'):
        times = np.linspace(tmin, tmax, lc_config['nbins'] + 1)
    elif lc_config.get('binsz'):
        times = np.arange(tmin, tmax, lc_config['binsz'])
    else:
        raise ValueError("None of time_bins, nbins, or binsz specified")
    return times


def run_lightcurve(fermipy_config, prefix, num_sections=None, section=0):
    """
    Run a Fermipy lightcurve analysis.

    Arguments:
        fermipy_config -- path to Fermipy config file
        prefix -- prefix for output files and output directory
        num_sections -- number of sections to split the lightcurve into
                        to compartmentalize the calculation when it crashes
                        if None, do not section the analysis
        section -- the section index to perform the analysis for
    """

    # Split the lightcurve into sections
    if num_sections is not None:
        times = get_times(fermipy_config)
        split_times = np.array_split(times, num_sections)
        selected_times = split_times[section]
        # Prevent the last bin from being lost due to the split
        if section < num_sections - 1:
            next_section_times = split_times[section + 1]
            if next_section_times.size != 0:
                last_bin = next_section_times[0]
                selected_times = np.append(selected_times, last_bin)
        selected_times = selected_times.tolist()
        fermipy_config['lightcurve']['time_bins'] = selected_times

    # Perform the analysis (on the selected section)
    gta = GTAnalysis(fermipy_config, logging={'verbosity': 3})
    gta.load_roi('{}.npy'.format(prefix))
    gta.lightcurve(fermipy_config['selection']['target'], make_plots=True)

    # Rename output file to prevent sections from overwriting each other
    outdir = fermipy_config['fileio']['outdir']
    outfiles = glob.glob(os.path.join(outdir, "*_lightcurve.*"))
    for outfile in outfiles:
        ext = outfile.rsplit('.', 1)[1]
        os.rename(outfile, os.path.join(outdir,
            "{}_lightcurve_{}.{}".format(prefix, section, ext)))


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

    if args.lightcurve:
        num_sections = pipeline_config['num_sections']
        section = pipeline_config['section']
        run_lightcurve(fermipy_config, prefix, num_sections, section)
    else:
        run_analysis(fermipy_config, prefix)

