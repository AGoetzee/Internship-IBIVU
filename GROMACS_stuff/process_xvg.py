#!/usr/bin/python3
"""
process_xvg.py: Processes GROMACS .xvg data files into plots.
"""

__author__ = 'Arthur Goetzee'
__version__ = '1.0'
__email__ = 'a.g.goetzee@student.vu.nl'

import pandas as pd
import matplotlib.pyplot as plt
import os
import argparse


def load_data(file, pandas=True):
    """
    Loads .xvg files and returns it in a python-friendly format

    Arguments:
        file (str): GROMACS .xvg file.
        pandas (logic): If set to True returns a pd dataframe (default). If false returns a python dict.
    """

    # Open file
    f = open(file)

    # Move data into lists
    cols = []
    rows = []  # This will be a nested list

    for line in f.readlines():
        # Skip comment lines and add x column name
        if line.startswith('#') or line.startswith('@') and not line.startswith('@ s'):
            if line.startswith('@    x'):
                cols.append(line[17:].strip(' "\n'))
            continue
            # Add y column names to list
        elif line.startswith('@ s'):
            cols.append(line.split('"')[-2]
                        .replace('#', '')
                        .replace('*', '')
                        .replace('.', ''))
            # Add row to list of rows
        else:
            rows.append(line.split())

    f.close()

    # Initialize df
    df = {}
    for col in cols:
        df[col] = []

    # Transpose data into dict
    for row in rows:
        for i, col in enumerate(cols):
            df[col].append(float(row[i]))

    # Load as pandas df
    if pandas:
        df = pd.DataFrame.from_dict(df).set_index(cols[0], drop=True)

    return df


def check_name(path, prefix, plotname):
    """
    Ensures that the name is not already present in working dir.
    Note: This function ensures plots won't be overwritten.

    Arguments:
        name: str. Relative path to check in.
        prefix: Prefix name for files
    """

    # Check if dir was given
    if path == '':
        files = os.listdir()
    else:
        try:  # Create dir if it doesnt exist
            files = os.listdir(path)
        except:
            os.mkdir(path)
            print(path, 'does not exist. Creating...')
            files = os.listdir(path)

    name = prefix + plotname
    i = -1  # Suffix to be given

    # Check if file exists
    while name + '.png' in files:
        name = prefix + plotname + str(i)
        i -= 1
    return path + name


def make_plots(df, name='', checkname=True):
    """
    Creates plots from all columns in the dataframe.

    Arguments:
        df: pandas dataframe or python dict
        name: str. optional. can help to specify a folder.
        checkname: bool. When True invokes check_name() function.
    """

    # Determine if input is python dict or pd dataframe
    pandas = type(df) == pd.core.frame.DataFrame

    if checkname:
        path = name.split('/')[:-1]
        path.append('')
        path = '/'.join(path)
        prefix = name.split('/')[-1]

    if pandas:
        for col in df.columns:

            if checkname:
                savename = check_name(path, prefix, col) + '.png'
            else:
                savename = name + col + '.png'

            plt.plot(df.index, df[col])
            plt.xlabel(df.index.name)
            plt.ylabel(col)

            plt.savefig(savename)
            plt.clf()

    else:
        df_cols = list(df.keys())
        for col in df_cols:

            if checkname:
                savename = check_name(path, prefix, col) + '.png'
            else:
                savename = name + col + '.png'

            if col == df_cols[0]:  # Skip X axis
                continue

            plt.plot(df[df_cols[0]], df[col])
            plt.xlabel(df_cols[0])
            plt.ylabel(col)

            plt.savefig(savename)
            plt.clf()


def parse_args():
    """
    Parses command line arguments. All args are optional.
    Arguments:
        -n --name: Argument for usage in make_plots()
        -v --verbose: Print verbose
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--name', action='store', type=str, help='Name/path for plot files', required=False,
                        default='')
    parser.add_argument('-v', '--verbose', action='store_true', help='Print verbose', required=False, default=False)
    args = parser.parse_args()

    return args


def main():
    # Get arguments
    args = parse_args()

    if args.verbose:
        print('GROMACS automatic XVG processor.\nArthur Goetzee, 2021\n\n')
        print('Flags given:\n-v')
        if args.name != '':
            print('-n', args.name, '\n')
        print('Looking for xvg files...')

    # Look for xvg files in working dir
    files = []
    for filename in os.listdir():

        if filename.endswith(".xvg"):
            files.append(filename)
            if args.verbose:
                print('Found', filename)

    if args.verbose:
        if not files:
            print('No files found. Exiting')
        print('Files found.\n')

    # Perform functions with found xvg files
    for file in files:
        if args.verbose:
            print('Processing', file)
        data = load_data(file)
        make_plots(data, args.name)

    print('All .xvg files have been processed')


main()
