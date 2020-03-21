import argparse
import os
import sys

import numpy as np
import yaml


def combine_lightcurve_results(outdir, prefix, num_sections):
    """
    Combine lightcurve output files from the Fermipy pipelie.

    Arguments:
        outdir -- analysis output directory
        prefix -- prefix for output files and output directory
        num_sections -- number of lightcurve sections
    """
    print("Combining lightcurve output files...")
    combined_results = {}

    for i in range(num_sections):
        path = os.path.join(outdir, "{}_lightcurve_{}.npy".format(prefix, i))
        try:
            results = np.load(path, encoding='latin1').item()
        except IOError:
            print("Section {}: no file found!".format(i))
            continue
        num_bins = len(results['tmin'])
        print("Section {}: {} bins".format(i, num_bins))
        non_array_keys = ['name', 'file', 'ts_var', 'config']
        for key, value in results.items():
            if key in non_array_keys:
                if key not in combined_results:
                    combined_results[key] = [value]
                else:
                    combined_results[key].append(value)
            else:
                if key not in combined_results:
                    combined_results[key] = value
                else:
                    arrays = [combined_results[key], value]
                    combined_results[key] = np.concatenate(arrays)

    filename = os.path.join(outdir, "{}_lightcurve_combined.npy".format(prefix))
    np.save(filename, combined_results)
    print("Combined results saved to: {}".format(filename))


if __name__ == "__main__":

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Combine lightcurve output files from the pipeline.")
    parser.add_argument('config', help="path to pipeline configuration file")
    parser.add_argument('-p', '--prefix',
                        help="analysis prefix (overrides config)")
    parser.add_argument('-n', '--num_sections', type=int,
                        help="number of sections (overrides config)")
    args = parser.parse_args()

    with open(args.config, 'r') as config_file:
        pipeline_config = yaml.safe_load(config_file)
    prefix = args.prefix if args.prefix else pipeline_config['prefix']
    num_sections = (args.num_sections if args.num_sections
                    else pipeline_config['num_sections'])

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
    outdir = os.path.join(os.getcwd(), fermipy_config['fileio']['outdir'])

    combine_lightcurve_results(outdir, prefix, num_sections)
