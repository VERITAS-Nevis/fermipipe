import argparse
import copy
import glob
import os
import shutil
import time

from astropy.io.fits.verify import VerifyError
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

    gta.print_roi()

    # Remove undetected sources to simplify the model
    gta.delete_sources(minmax_ts=[None, 4], exclude=['galdiff', 'isodiff'])
    gta.delete_sources(minmax_npred=[None, 1], exclude=['galdiff', 'isodiff'])

    gta.print_roi()

    # Free high TS sources and those close to the ROI center
    gta.free_sources(minmax_ts=[25, None])
    gta.free_sources(distance=5.0)

    # Free all parameters of the galactic and isotropic diffuse components
    gta.free_source('galdiff')
    gta.free_source('isodiff')

    # Perform the model fit, and fix all sources again when done
    gta.fit()
    gta.free_sources(free=False)

    gta.print_roi()

    # Plot TS and residual maps
    model = {
        'SpatialModel': 'PointSource',
        'Index': 2.0,
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
    outdir = os.path.join(os.getcwd(), fermipy_config['fileio']['outdir'])
    original_config = copy.deepcopy(fermipy_config)

    times = get_times(fermipy_config)
    if num_sections is None:
        selected_times = times
    else:  # Split the lightcurve into sections
        split_times = np.array_split(times, num_sections)
        selected_times = split_times[section]
        # Prevent the last bin from being lost due to the split
        if section < num_sections - 1:
            next_section_times = split_times[section + 1]
            if next_section_times.size != 0:
                last_bin = next_section_times[0]
                selected_times = np.append(selected_times, last_bin)

    # Consolidate bins encompassed by periods with no data
    # such as non-science periods or ToOs
    # Otherwise, analysis will encounter a divide by zero error and segfault

    # https://fermi.gsfc.nasa.gov/ssc/observations/timeline/posting/cal/
    no_data_periods = {
        (586362562, 586479265),  # Non-science Sunpoint
        (560543825, 562120025),  # Extended south (+50) rock
        (557887205, 558482705),  # Extended south (+50) rock
        (556044185, 556770005),  # Extended south (+50) rock
        (554261765, 554843765),  # Extended south (+50) rock
        (551667365, 552502385),  # Extended south (+50) rock
        (550290005, 551107805),  # Extended south (+50) rock
        (547923665, 548557805),  # Extended south (+50) rock
        (542869922, 546791705),  # Spacecraft operational anomaly, recovery,
                                 # and extended south (+50) rock
        (491961604, 492389764),  # ToO PSR J1119
        (415063383, 415324023),  # Solar ToO
        (407981463, 408414003),  # ToO Nova Cen
        (258507868, 258671110),  # Non-science Unknown
        }

    for period_start, period_end in no_data_periods:
        # If there's only a few hours of data, the fit will likely fail
        # Add a margin of error to avoid this
        margin_of_error = 10800  # 3 hours
        period_start -= margin_of_error
        period_end += margin_of_error
        # Get the bins encompassed by each no-data period
        period_times = selected_times[(selected_times >= period_start)
                                      & (selected_times <= period_end)]
        if period_times.size == 0:
            continue
        # Remove the encompassed bins
        # Note that this doesn't guarantee that the resulting merged bin
        # has sufficient exposure, and that data may be thrown out
        # if the non-science period was at the start or end of the lightcurve
        selected_times = np.setdiff1d(selected_times, period_times,
                                      assume_unique=True)

    selected_times = selected_times.tolist()
    print("Time bins:", selected_times)
    fermipy_config['lightcurve']['time_bins'] = selected_times

    # Perform the analysis (on the selected section)
    gta = GTAnalysis(fermipy_config, logging={'verbosity': 3})

    for __ in range(3):
        try:
            gta.load_roi('{}.npy'.format(prefix))
            break
        except RuntimeError:
            # Maybe the file is in use? Wait, then try again
            time.sleep(5)

    try:
        gta.lightcurve(fermipy_config['selection']['target'], make_plots=True)
    except (AttributeError, IOError, OSError, RuntimeError, TypeError,
            VerifyError) as err:

        if type(err) == RuntimeError:
            print(str(err))
            allowed_runtime_error_starts = [
                "File not in FITS or Root format",
                "Requested energy",
                "Could not open FITS extension",
                "Failed to converge after",
                "File not found",
                ]
            allowed_start = False
            for start in allowed_runtime_error_starts:
                if str(err).startswith(start):
                    allowed_start = True
            if not allowed_start:
                raise err

        # There is a missing or corrupt file in the directory for a bin
        # Delete the directory of the last bin to enable a clean restart
        bin_dirs = glob.glob(os.path.join(outdir, "lightcurve_*"))
        print("selected_times:", selected_times)

        def selected(bin_dir):
            bin_edges = bin_dir.split('/')[-1].split('_')[1:]
            for edge in bin_edges:
                if int(edge) not in selected_times:
                    return False
            return True

        selected_bin_dirs= [bin_dir for bin_dir in bin_dirs
                            if selected(bin_dir)]
        # Remove the last modified directory
        last_bin_dir = max(selected_bin_dirs, key=os.path.getmtime)
        print("Removing corrupted directory {} ...".format(last_bin_dir))
        shutil.rmtree(last_bin_dir)

        # Restart the lightcurve analysis
        print("Restarting the lightcurve analysis...")
        return run_lightcurve(original_config, prefix, num_sections, section)

    # Rename output file to prevent sections from overwriting each other
    if num_sections is not None:
        outfiles = glob.glob(os.path.join(outdir, "*_lightcurve.*"))
        for outfile in outfiles:
            ext = outfile.rsplit('.', 1)[1]
            os.rename(outfile, os.path.join(outdir,
                "{}_lightcurve_{}.{}".format(prefix, section, ext)))


if __name__ == "__main__":

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Run a Fermipy analysis.")
    parser.add_argument('config', help="path to pipeline configuration file")
    parser.add_argument('-p', '--prefix',
                        help="analysis prefix (overrides config)")
    parser.add_argument('-l', '--lightcurve', action='store_true',
                        help="perform a lightcurve analysis")
    parser.add_argument('-s', '--section', type=int,
                        help="lightcurve section (overrides config)")
    args = parser.parse_args()

    with open(args.config, 'r') as config_file:
        pipeline_config = yaml.safe_load(config_file)
    prefix = args.prefix if args.prefix else pipeline_config['prefix']

    # Set the Fermipy config filename to a default if not otherwise specified
    if not pipeline_config.get('fermipy_config'):
        pipeline_config['fermipy_config'] = prefix + '_config.yml'

    # Load Fermipy configuration
    with open(pipeline_config['fermipy_config'], 'r') as config_file:
        fermipy_config = yaml.safe_load(config_file)

    # Set the outdir to be the same as the prefix
    if 'fileio' not in fermipy_config:
        fermipy_config['fileio'] = {}
    if 'outdir' not in fermipy_config['fileio']:
        fermipy_config['fileio']['outdir'] = prefix

    if args.lightcurve:
        num_sections = pipeline_config['num_sections']
        section = (args.section if args.section is not None
                   else pipeline_config['section'])
        run_lightcurve(fermipy_config, prefix, num_sections, section)
    else:
        run_analysis(fermipy_config, prefix)

