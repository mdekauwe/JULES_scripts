#!/usr/bin/env python

"""
Plot visual benchmark (average seasonal cycle) of old vs new model runs.

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (18.10.2017)"
__email__ = "mdekauwe@gmail.com"

import netCDF4 as nc
import matplotlib.pyplot as plt
import sys
import datetime as dt
import pandas as pd
import numpy as np
from matplotlib.ticker import FixedLocator
import xarray as xr

def main(fname, plot_fname=None):

    df = read_jules_file(fname)
    df = resample_to_seasonal_cycle(df)

    fig = plt.figure(figsize=(6,9))
    fig.subplots_adjust(hspace=0.3)
    fig.subplots_adjust(wspace=0.2)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['font.size'] = 12
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12

    ax1 = fig.add_subplot(3,2,1)
    ax2 = fig.add_subplot(3,2,2)
    ax3 = fig.add_subplot(3,2,3)
    ax4 = fig.add_subplot(3,2,4)
    ax5 = fig.add_subplot(3,2,5)
    ax6 = fig.add_subplot(3,2,6)

    axes = [ax1, ax2, ax3, ax4, ax5, ax6]
    vars = ["GPP", "Qh","Qle", "LAI", "TVeg", "Evap"]

    for a, v in zip(axes, vars):
        a.plot(df.month, df[v], c="black", lw=2.0, ls="-")

    labels = ["GPP (g C m$^{-2}$ d$^{-1}$)","Qh (W m$^{-2}$)",\
              "Qle (W m$^{-2}$)", "LAI (m$^{2}$ m$^{-2}$)",\
              "TVeg (mm d$^{-1}$)", "Evap (mm d$^{-1}$)"]
    for a, l in zip(axes, labels):
        a.set_title(l, fontsize=12)

    xtickagaes_minor = FixedLocator([2, 3, 4, 5, 7, 8, 9, 10, 11])
    for i,a in enumerate(axes):
        a.set_xticks([1, 6, 12])
        if i != 1:
            a.set_ylim(ymin=0)
        a.xaxis.set_minor_locator(xtickagaes_minor)
        a.set_xticklabels(['Jan', 'Jun', 'Dec'])
        if i < 4:
            plt.setp(a.get_xticklabels(), visible=False)

    if plot_fname is None:
        plt.show()
    else:
        fig.savefig(plot_fname, bbox_inches='tight', pad_inches=0.1)


def read_jules_file(fname, type=None):

    vars_to_keep = ['GPP','Qle','Qh','LAI','TVeg','Evap']

    ds = xr.open_dataset(fname, decode_times=False)
    time_jump = int(ds.time[1].values) - int(ds.time[0].values)

    if time_jump == 3600:
        freq = "H"
    elif time_jump == 1800:
        freq = "30M"
    elif time_jump == 10800:
        freq = "3H"
    else:
        raise("Time problem")

    units, reference_date = ds.time.attrs['units'].split('since')
    #df = ds[vars_to_keep].squeeze(dim=["x","y","soil"], drop=True).to_dataframe()

    ds = ds[vars_to_keep].squeeze(dim=["x","y"], drop=True)
    df = pd.DataFrame(ds['GPP'].values[:], columns=['GPP'])
    df['Qle'] = ds['Qle'].values[:]
    df['Qh'] = ds['Qh'].values[:]
    df['LAI'] = ds['LAI'].values[:]
    df['TVeg'] = ds['TVeg'].values[:]
    df['Evap'] = ds['Evap'].values[:]

    start = reference_date.strip().split(" ")[0].replace("-","/")
    df['dates'] = pd.date_range(start=start, periods=len(df), freq=freq)
    df = df.set_index('dates')

    return df


def resample_to_seasonal_cycle(df, OBS=False):

    UMOL_TO_MOL = 1E-6
    MOL_C_TO_GRAMS_C = 12.0
    SEC_2_HOUR = 3600.
    HOUR_2_DAY = 24.
    KG_2_G = 1000.

    # umol/m2/s -> g/C/d
    #df['GPP'] *= UMOL_TO_MOL * MOL_C_TO_GRAMS_C * SEC_2_HOUR * HOUR_2_DAY

    # kg/m2/s -> g/C/d
    df['GPP'] *= KG_2_G * SEC_2_HOUR * HOUR_2_DAY

    # kg/m2/s -> mm/d
    df['TVeg'] *= SEC_2_HOUR * HOUR_2_DAY
    df['Evap'] *= SEC_2_HOUR * HOUR_2_DAY

    method = {'GPP':'mean', 'Qle':'mean', 'Qh':'mean', 'LAI':'mean',
              'TVeg':'mean', 'Evap':'mean'}
    df = df.resample("M").agg(method).groupby(lambda x: x.month).mean()
    df['month'] = np.arange(1,13)

    return df

if __name__ == "__main__":

    """
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("-f", "--fname", dest="fname",
                      action="store", help="filename",
                      type="string")
    parser.add_option("-p", "--plot_fname", dest="plot_fname", action="store",
                      help="Benchmark plot filename", type="string")
    (options, args) = parser.parse_args()
    main(options.fname, options.plot_fname)
    """

    fname = "/Users/xj21307/research/JULES/runs/roses/AU_Tum/outputs/AU_Tum/local_AU-Tum_fluxnet2015_AU-Tum.AU-Tum.nc"

    main(fname)
