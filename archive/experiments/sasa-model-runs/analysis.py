import pytraj as pt  # For Amber trajectories
import pandas as pd  # Data analysis
import numpy as np  # Data analysis
import matplotlib.pyplot as plt  # Plotting
import plumed  # For loading PLUMED files
import os
from celluloid import Camera  # For creating gifs
from subprocess import call, DEVNULL, Popen
import re


# call(["source", "~/plumed/source.sh"])
# call(["source", "~/amber20/amber.sh"])

def check_dirs(path):
    if not os.path.exists(path + '/figs'):
        print(f'{path}/figs does not exist, creating...')
        os.makedirs(f'{path}/figs')

    if not os.path.exists(path + '/fesdata'):
        print(f'{path}/fesdata does not exist, creating...')
        os.makedirs(f'{path}/fesdata')


def get_fes_files(path, stride=None, out="fesdata/fes", verbose=False):
    """Calls PLUMED sum_hills utility to obtain the FES data. 
    When stride is specified, hills are summed over the stride interval.
    
    :param path: str - Specifies the path to the HILLS file.
    :param stride: int - Specifies size of stride.
    :param out: str - Specifies output directory. Default 'fesdata/'.
    """

    if verbose:
        verbosity = None
    else:
        verbosity = DEVNULL

    if stride is not None:
        assert isinstance(stride, int), f"stride should be int, got {type(stride)} instead."
        # call(map(str, ["plumed", "sum_hills", "--hills", path, "--stride", stride, "--outfile", out]),stdout=verbosity,env=os.environ)
        call(f"plumed sum_hills --hills {path} --stride {stride} --outfile {out}_", shell=True)
    else:
        # call(map(str, ["plumed", "sum_hills", "--hills", path, "--outfile", out+'_']),stdout=verbosity,env=os.environ)
        call(f"plumed sum_hills --hills {path} --outfile {out}", shell=True)
    return


def plot_fesdata(data, figname='fes', title='Free Energy Landscape', suffix=''):
    """Plots FES data using matplotlib and saves it.
    
    :param data: DataFrame or str - Either path to FES data or data object itself
    :param figname: str - Name (and path) to give to the plot.
    :param title: str - Title to give to the plot.
    :param suffix: str - Suffix string to appended to the title."""

    # Check if this is a str. if so convert to df
    if isinstance(data, str):
        data = plumed.read_as_pandas(data)

    plt.figure(dpi=150)
    plt.plot(data['t1'], data['file.free'])
    plt.xlabel('RSD')
    plt.ylabel('Relative FE (kJ/mol)')

    title += suffix
    plt.title(title + suffix)

    if not figname.endswith('.png'):
        figname += '.png'

    plt.savefig(figname)

    return


def plot_convergence(path, file_suffix='fes_', step=1, unit='', figname='/figs/conv', title='Free Energy Landscape',
                     suffix=''):
    """Plots convergence data.
    
    :param step: Stepsize for strides. Default 1.
    :param file_suffix: str - Suffix title for fes files.
    :param path: str - Path to fesdata directory.
    :param step: list - List containing names of data to be plotted.
    :param unit: str - Unit to be given in legend.
    :param figname: str - Name (and path) to give to the plot.
    :param title: str - Title to give to the plot.
    :param suffix: str - Suffix string to appended to the title."""

    plt.figure(dpi=150)

    if not path.endswith('/'):
        path += '/'

    if not figname.endswith('.png'):
        figname += '.png'

    # Get number of strides
    with os.scandir(path) as files:
        max = 0
        for f in files:
            try:
                filenum = re.sub('\D', '', f.name)
                if int(filenum) > max:
                    max = int(filenum)
            except:
                print(f'{f} is not a number.')

    if isinstance(step, list):  # Only plot steps in list
        for i in step:
            temp = plumed.read_as_pandas(f'{path}{file_suffix}{i}.dat')
            plt.plot(temp['t1'], temp['file.free'], label=f'{i + 1} {unit}', linewidth=0.5)
    elif isinstance(step, int):
        assert step >= 1, f'ERROR: Step cannot be smaller than 1, got {step} instead.'
        for i in range(0, max + 1, step):
            temp = plumed.read_as_pandas(f'{path}{file_suffix}{i}.dat')
            plt.plot(temp['t1'], temp['file.free'], label=f'{i + 1} {unit}', linewidth=0.5)
    elif step == 'all':
        for i in range(0, max + 1):
            temp = plumed.read_as_pandas(f'{path}{file_suffix}{i}.dat')
            plt.plot(temp['t1'], temp['file.free'], label=f'{i + 1} {unit}', linewidth=0.5)
    else:
        print('Something went wrong, that\'s all we know.')

    # plt.plot(fesdata['t1'], fesdata['file.free'], linestyle='dashed', linewidth=2.5, label='Final', color='black')
    plt.legend()

    plt.xlabel('RSD')
    plt.ylabel('Relative FE (kJ/mol)')

    plt.title(title + suffix)
    plt.savefig(figname)

    return


def merge_colvar_hills(colvar, fesdata, convert_time=True):
    # Check if this is a str. if so convert to df
    if isinstance(colvar, str):
        colvar = plumed.read_as_pandas(colvar)

    # Check if this is a str. if so convert to df
    if isinstance(fesdata, str):
        fesdata = plumed.read_as_pandas(fesdata)

    # Fix time to be in line with trajectory data from Amber
    if convert_time:
        colvar = colvar[colvar['time'] % 4 == 0]
        colvar['time'] = colvar['time'] + 20

    # Join df's on RMSD
    colvar['RMSD_bin'] = np.digitize(colvar['t1'], fesdata['t1'])
    fesdata['RMSD_bin'] = np.digitize(fesdata['t1'], fesdata['t1'])
    df = pd.merge(colvar, fesdata, how='inner', on='RMSD_bin')

    return df


def create_features(df):
    # Create combi function and translate values
    df['combi'] = np.sqrt((df['c2'].values - 1.01) ** 2 + (df['c27'].values - 1.09) ** 2)
    df['c15_translated'] = df['c15'] - 1.01

    # Remove data involving barrier
    df = df[df['t1_y'] < 6]

    return df


def create_projection(df, figname='/figs/projection', title='Projection of FE to individual distances', suffix=''):
    assert 'file.free' in df.columns, "Data requires FES data. Run merge_colvar_hills() first."
    assert 'c15_translated' in df.columns, "Data requires engineered features. Run create_features() first."

    plt.figure(dpi=150)
    plt.scatter(df['combi'], df['c15_translated'], c=df['file.free'].values)

    plt.colorbar(label='Relative FE (kJ/mol)')
    plt.xlabel('c2 and c27')
    plt.ylabel('c15 (nm)')
    plt.title(title + suffix)

    if not figname.endswith('.png'):
        figname += '.png'

    plt.savefig(figname)

    return


def create_conv_animation(path, file_suffix='fes_', unit='', figname='/figs/conv_time', title='FES through time',
                          suffix=''):
    if not path.endswith('/'):
        path += '/'

    if not figname.endswith('.gif'):
        figname += '.gif'

    # Get number of strides
    with os.scandir(path) as files:
        max = 0
        for f in files:
            try:
                filenum = re.sub('\D', '', f.name)
                if int(filenum) > max:
                    max = int(filenum)
            except:
                print(f'{f} is not a number.')

    fig = plt.figure(dpi=150)
    camera = Camera(fig)
    for i in range(0, max + 1):
        temp = plumed.read_as_pandas(f'{path}{file_suffix}{i}.dat')
        t = plt.plot(temp['t1'], temp['file.free'], label=f'{i + 1} {unit}', linewidth=0.5)
        plt.title(title + suffix)
        plt.legend(t, [f'{i + 1}'])

        camera.snap()

    animation = camera.animate()
    animation.save(figname, writer='imagemagick')


def create_projection_animation(df, stride=250, figname='/figs/proj_time'):
    assert 'file.free' in df.columns, "Data requires FES data. Run merge_colvar_hills() first."
    assert 'c15_translated' in df.columns, "Data requires engineered features. Run create_features() first."

    fig = plt.figure()
    camera = Camera(fig)
    for i in range(0, df.shape[0] + 1, stride):
        t = plt.scatter(df.loc[:i, ['combi']], df.loc[:i, ['c15_translated']], c=df.loc[:i, ['file.free']].values,
                        label=f'{i} ns')
        plt.xlim([0, 4])
        plt.ylim([-0.5, 1.5])
        camera.snap()

    if not figname.endswith('.gif'):
        figname += '.gif'

    animation = camera.animate()
    animation.save(figname, writer='imagemagick')

def generate_difference_plot(path, file_suffix='fes_', step=1, unit='', figname='/figs/diff', title='Free Energy Landscape',
                     suffix=''):

    plt.figure(dpi=150)

    diffs = []
    mins = []
    means = []

    # Get number of strides
    with os.scandir(path) as files:
        max = 0
        for f in files:
            try:
                filenum = re.sub('\D', '', f.name)
                if int(filenum) > max:
                    max = int(filenum)
            except:
                print(f'{f} is not a number.')


    if isinstance(step, list):  # Only plot steps in list
        iterator = step
    elif isinstance(step, int):
        assert step >= 1, f'ERROR: Step cannot be smaller than 1, got {step} instead.'
        iterator = range(0, max + 1, step)
    elif step == 'all':
        iterator = range(0, max + 1)
    else:
        print('Something went wrong, that\'s all we know.')

    for i in iterator:
        temp = plumed.read_as_pandas(f'{path}{file_suffix}{i}.dat')
        mean = temp[temp['t1'].between(2, 6)]['file.free'].mean()
        minim = temp[temp['t1'].between(0, 1)]['file.free'].min()
        diff = temp['file.free'].min() - mean

        diffs.append(diff)
        mins.append(minim)
        means.append(mean)

    plt.plot(range(len(diffs)), diffs, linestyle='dashed', linewidth=2.5, label='Final', color='black')
    plt.legend()

    plt.xlabel('Time (ns)')
    plt.ylabel('Difference in relative FE (kJ/mol)')

    plt.title(title + suffix)
    plt.savefig(figname)

    return



if __name__ == '__main__':

    master_dir = os.getcwd()
    print('Current dir:', master_dir)

    if not master_dir.endswith('/'):
        master_dir += '/'

    sim_dirs = []
    with os.scandir(master_dir) as dirs:
        for entry in dirs:
            if entry.name.isnumeric():
                sim_dirs.append(entry.name)

    print(f'Found sim dirs {sim_dirs}')

    for sim in sim_dirs:
        dirname = master_dir + sim
        print(f'\t****Processing {sim}*****')

        check_dirs(dirname)

        # Currently I can't get the terminal calls to work properly. This is done in a separate script.
        # print('Getting FES files...')
        # get_fes_files(dirname + f'/HILLS_{sim}',verbose=True)
        # get_fes_files(dirname + f'/HILLS_{sim}', stride=2000,verbose=True)
        # print('FES files obtained!')

        print('Plotting FES data...')
        plot_fesdata(dirname + '/fesdata/fes.dat', figname=dirname + '/figs/fes')

        print('Plotting convergence...')
        plot_convergence(dirname + '/fesdata/', step=10, figname=dirname + '/figs/convergence')

        print('Getting difference through time...')
        generate_difference_plot(dirname + '/fesdata/', step='all', figname=dirname + '/figs/difference')

        print('Creating convergence animation...')
        create_conv_animation(dirname + '/fesdata/', figname=dirname + '/figs/convergence_time')


        print('Processing FES data for projection...')
        df = merge_colvar_hills(colvar=dirname + f'/COLVAR_{sim}', fesdata=dirname + '/fesdata/fes.dat')
        df = create_features(df)
        print('Plotting projection...')
        create_projection(df, figname=dirname + '/figs/projection')
        create_projection_animation(df, figname=dirname + '/figs/proj_time.gif')
        print(f'\t****Finished {sim}*****')

