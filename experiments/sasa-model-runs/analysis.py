import pytraj as pt # For Amber trajectories
import pandas as pd # Data analysis
import numpy as np # Data analysis
import matplotlib.pyplot as plt # Plotting
import plumed # For loading PLUMED files
import os
from celluloid import Camera # For creating gifs
from subprocess import call
call(["source", "~/plumed/source.sh"])
call(["source", "~/amber20/amber.sh"])


def get_fes_files(path,stride,out="fesdata/"):
    """Calls PLUMED sum_hills utility to obtain the FES data. 
    When stride is specified, hills are summed over the stride interval.
    
    :param stride: int - Specifies size of stride.
    :param out: str - Specifies output directory. Default 'fesdata/'
    """
    
    if stride is not None:
        assert type(stride) == int(), f"stride should be int, got {type(stride)} instead."
        call(["plumed", "sum_hills", "--hills", path, "--stride", stride,"--outfile",out])
    else:
        call(["plumed", "sum_hills", "--hills", path, "--outfile", out])

    return

def plot_fesdata(data,figname='figs/fes',title='Free Energy Landscape',suffix=''):
    """Plots FES data using matplotlib and saves it.
    
    :param data: DataFrame or str - Either path to FES data or data object itself
    :param figname: str - Name (and path) to give to the plot.
    :param title: str - Title to give to the plot.
    :param suffix: str - Suffix string to appended to the title."""
    
    
    # Check if this is a str. if so convert to df
    if type(data) is str():
        data = plumed.read_as_pandas(data)
    
    plt.figure(dpi=150)
    plt.plot(fesdata['t1'],fesdata['file.free'])
    plt.xlabel('RSD')
    plt.ylabel('Relative FE (kJ/mol)')
    
    title += suffix
    plt.title(title+suffix)
    
    if not figname.endswith('.png'):
        figname += '.png'
    
    plt.savefig(figname)
    
    return

def plot_convergence(path, conv_list, unit='', figname='figs/fes',title='Free Energy Landscape',suffix=''):
    """Plots convergence data
    
    :param path: str - Path to fesdata directory.
    :param conv_list: list - List containing names of data to be plotted.
    :param unit: str - Unit to be given in legend.
    :param figname: str - Name (and path) to give to the plot.
    :param title: str - Title to give to the plot.
    :param suffix: str - Suffix string to appended to the title."""
    
    plt.figure(dpi=150)
    
    if not path.endswith('/'):
        path += '/'
    
    if not figname.endswith('.png'):
        figname += '.png'
    
    for i in times:
        temp = plumed.read_as_pandas(f'{TEMP}{i}.dat')
        plt.plot(temp['t1'],temp['file.free'],label=f'{i+1} {unit}',linewidth=0.5)
    
    plt.plot(fesdata['t1'],fesdata['file.free'],linestyle='dashed',linewidth=2.5, label='Final', color='black')
    plt.legend()

    plt.xlabel('RSD')
    plt.ylabel('Relative FE (kJ/mol)')

    plt.title(title+suffix)
    plt.savefig(figname)

    return

def merge_colvar_hills(colvar,fesdata,convert_time=True, create_bins=True)    
        
   
    # Check if this is a str. if so convert to df
    if type(colvar) is str():
        colvar = plumed.read_as_pandas(colvar)
    
    # Check if this is a str. if so convert to df
    if type(fesdata) is str():
        fesdata = plumed.read_as_pandas(fesdata)
    
    # Fix time to be in line with trajectory data from Amber
    if convert_time:
        colvar = colvar[colvar['time']%4==0]
        colvar['time'] = colvar['time'] + 20
    
    # Join df's on 
    colvar['RMSD_bin'] = np.digitize(colvar['t1'],fesdata['t1'])
    fesdata['RMSD_bin'] = np.digitize(fesdata['t1'],fesdata['t1'])
    df = pd.merge(colvar,fesdata,how='inner',on='RMSD_bin')
    
    return df

def create_features(df):
    
    assert type(df) is pd.DataFrame(), f"Expected DataFrame, got {type(df)} instead."
    
    # Create combi function and translate values
    df['combi'] = np.sqrt((df['c2'].values - 1.01 )**2 + (df['c27'].values - 1.09) ** 2)
    df['c15_translated'] = df['c15'] - 1.01
    
    # Remove data involving barrier
    df=df[df['t1_y']<6]
    
    return df
    
def create_projection(df, figname='figs/projection',title='Projection of FE to individual distances',suffix=''):
    
    assert type(df) is pd.DataFrame(), f"Expected DataFrame, got {type(df)} instead."
    assert 'file.free' in df.columns, "Data requires FES data. Run merge_colvar_hills() first.'
    assert 'c15_translated' in df.columns, "Data requires engineered features. Run create_features() first."
    
    plt.figure(dpi=150)
    plt.scatter(df['combi'],df['c15_translated'],c=df['file.free'].values)
    
    plt.colorbar(label='Relative FE (kJ/mol)')
    plt.xlabel('c2 and c27')
    plt.ylabel('c15 (nm)')
    plt.title(title+suffix)
    
    if not figname.endswith('.png'):
        figname += '.png'
    
    plt.savefig(figname)
    
    return

if __name__ == '__main__':

    master_dir = '/home/goetzee/simulations/sasa-model-runs/'
    
    sim_dirs = []
    for root, dirs, files in os.walk('/home/goetzee/simulations/sasa-model-runs/'):
      for f in files:
        if f.isnumeric():
          sim_dirs.append(f)
    
    print(f'Found sim dirs {sim_dirs}')
    
    for sim in sim_dirs 
        dirname = master_dir + sim
        
        get_fes_files(dirname+f'/HILLS_{sim}')
        get_fes_files(dirname+f'/HILLS_{sim}',stride=2000)
        
        plot_fesdata()