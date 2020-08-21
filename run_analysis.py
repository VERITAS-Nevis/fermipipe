import argparse
import copy
import glob
import os
import shutil
import time

from astropy.io import fits
from fermipy.gtanalysis import GTAnalysis
import numpy as np
import yaml


def run_analysis(fermipy_config, prefix,
                 delete_source=None, delete_sources=None,
                 free_source=None, free_sources=None):
    """
    Run a Fermipy analysis.

    Arguments:
        fermipy_config -- path to Fermipy config file
        prefix -- prefix for output files and output directory
        delete_source --
            list of names of sources to delete from the model
        delete_sources --
            list of dictionaries containing parameters
            defining sources to delete from the model,
            with isodiff, galdiff, and all custom sources
            automatically excluded
        free_source -- list of names of sources to free in the model fitting
        free_sources --
            list of dictionaries containing parameters
            defining sources to free in the model fitting
    """

    # Set default values
    def default(arg):
        return [] if arg is None else arg

    delete_source = default(delete_source)
    delete_sources = default(delete_sources)
    free_source = default(free_source)
    free_sources = default(free_sources)

    # Prepare the analysis
    gta = GTAnalysis(fermipy_config, logging={'verbosity': 3})

    # Prepare the data and perform calculations needed for the analysis
    gta.setup()

    # Fit the normalization and spectral shape parameters
    # of all model components and compute the TS of all sources in the ROI
    gta.optimize()

    gta.print_roi()

    # Remove sources to simplify the model
    custom_sources = fermipy_config.get('model', {}).get('sources', [])
    exclude = [source['name'] for source in custom_sources]
    exclude.extend(['galdiff', 'isodiff'])
    for source in delete_source:
        gta.delete_source(source)
    for args in delete_sources:
        args['exclude'] = exclude
        gta.delete_sources(**args)

    gta.print_roi()

    # Free the specified sources (all others fixed)
    for source in free_source:
        gta.free_source(source)
    for args in free_sources:
        gta.free_sources(**args)

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

def run_lightcurve(fermipy_config, prefix, num_sections=None, section=0,
                   no_data_periods=None):
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
    def get_no_data_periods():

        DEG_TO_RAD = np.pi / 180
        RAD_TO_DEG = 180 / np.pi

        # Output of gtmktime from base analysis
        events_file = os.path.join(outdir, "ft1_00.fits")
        with fits.open(events_file) as events:
            data = events[1].data
            time_met = data.TIME
            event_ra = data.RA * DEG_TO_RAD
            event_dec = data.DEC * DEG_TO_RAD

        with fits.open(fermipy_config['data']['scfile']) as spacecraft:
            data = spacecraft[1].data
            pointing_ra = data.RA_SCZ * DEG_TO_RAD
            pointing_dec = data.DEC_SCZ * DEG_TO_RAD
            edges = np.append(data.START, data.STOP[-1])

        containdata = np.digitize(time_met, edges) - 1
        pointing_ra = pointing_ra[containdata]
        pointing_dec = pointing_dec[containdata]

        def angular_distance(ra1, dec1, ra2, dec2):
            distance = np.arccos(np.sin(dec1)*np.sin(dec2)
                                 + np.cos(dec1)*np.cos(dec2)*np.cos(ra1 - ra2))
            return distance

        # Distance between the telescope pointing z-direction and each event
        dist_events = RAD_TO_DEG * angular_distance(pointing_ra, pointing_dec,
                                                    event_ra, event_dec)

        # Consider only events within 50 degrees of the telescope Z-direction
        # Chosen to roughly match time bins previously chosen manually,
        # through trial and error
        time_met = time_met[dist_events < 50]
        starts = time_met[:-1]
        stops = time_met[1:]

        min_bin_size = np.amin(times[1:] - times[:-1])
        no_data = (stops - starts) > min_bin_size
        no_data_periods = list(zip(starts[no_data], stops[no_data]))
        return no_data_periods

    if no_data_periods is None:
        no_data_periods = get_no_data_periods()

    for period_start, period_end in no_data_periods:
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
            fits.VerifyError) as err:

        if type(err) == RuntimeError:
            print(str(err))
            allowed_runtime_error_starts = [
                "File not in FITS or Root format",
                "Requested energy",
                "Could not open FITS extension",
                "Failed to converge after",
                "File not found",
                "mage in extension",
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
        return run_lightcurve(original_config, prefix, num_sections, section,
                              no_data_periods)

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
        run_analysis(fermipy_config, prefix,
                     delete_source=pipeline_config.get('delete_source'),
                     delete_sources=pipeline_config.get('delete_sources'),
                     free_source=pipeline_config.get('free_source'),
                     free_sources=pipeline_config.get('free_sources'))
