#!/usr/bin/env python

"""
Plot GPP and LE, showing standard CABLE, CABLE + hydraulics model and Ozflux
OBS.

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (23.03.2022)"
__email__ = "mdekauwe@gmail.com"

import netCDF4 as nc
import matplotlib.pyplot as plt
import sys
import datetime as dt
import pandas as pd
import numpy as np
from matplotlib.ticker import FixedLocator
import datetime
import os
import glob
from optparse import OptionParser
import string
from summary_stats import rmse, bias, nash_sutcliffe, willmott_agreement_indx
from scipy.stats import linregress
import xarray as xr

def main(site, met_fname, flux_fname, fname1, fname2, plot_fname=None):

    df1 = read_jules_file(fname1)
    df1 = resample_timestep(df1, type="JULES")

    df2 = read_cable_file(fname2, type="CABLE")
    df2 = resample_timestep(df2, type="CABLE")

    df_flx = read_cable_file(flux_fname, type="FLUX")
    df_flx = resample_timestep(df_flx, type="FLUX")

    df_met = read_cable_file(met_fname, type="MET")
    df_met = resample_timestep(df_met, type="MET")

    df1_drt = df1[(df1.index > '2006-11-1') & (df1.index <= '2007-4-1')]
    df2_drt = df2[(df2.index > '2006-11-1') & (df2.index <= '2007-4-1')]
    df_flx_drt = df_flx[(df_flx.index > '2006-11-1') & (df_flx.index <= '2007-4-1')]


    print("LE - JULES")
    m = df1_drt.Qle.values
    o = df_flx_drt.Qle.values
    print("RMSE = %.2f" % rmse(m, o))
    print("Nash-Sutcliffe Coefficient = %.2f" % nash_sutcliffe(m, o))
    slope, intercept, r_value, p_value, std_err = linregress(m, o)
    print("Pearson's r = %.2f" % (r_value))

    print("\n")
    print("LE - CABLE")
    m = df2_drt.Qle.values
    print("RMSE = %.2f" % rmse(m, o))
    print("Nash-Sutcliffe Coefficient = %.2f" % nash_sutcliffe(m, o))
    slope, intercept, r_value, p_value, std_err = linregress(m, o)
    print("Pearson's r = %.2f" % (r_value))

    print("\nJULES min LAI: ", df1.LAI.min())


    fig = plt.figure(figsize=(9,9))
    fig.subplots_adjust(hspace=0.1)
    fig.subplots_adjust(wspace=0.2)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['font.size'] = 12
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12

    labels_gen = label_generator('lower', start="(", end=")")
    colours = plt.cm.Set2(np.linspace(0, 1, 7))

    ax1 = fig.add_subplot(4,1,1)
    ax2 = fig.add_subplot(4,1,2)
    ax3 = fig.add_subplot(4,1,3)
    ax4 = fig.add_subplot(4,1,4)
    axx1 = ax1.twinx()
    axx2 = ax2.twinx()
    axx3 = ax3.twinx()
    axx4 = ax4.twinx()

    axes = [ax1, ax2, ax3, ax4]
    axes2 = [axx1, axx2, axx3, axx4]
    vars = ["GPP", "Qle", "LAI", "theta_weight"]

    props = dict(boxstyle='round', facecolor='white', alpha=1.0,
                 ec="white")
    for a, x, v in zip(axes, axes2, vars):

        if v != "LAI" and v != "theta_weight":
            a.plot(df_flx[v].index.to_pydatetime(),
                df_flx[v].rolling(window=5).mean(), c=colours[1], lw=2.0,
                ls="-", label="Observations")
        a.plot(df1[v].index.to_pydatetime(), df1[v].rolling(window=5).mean(),
               c=colours[0], lw=1.5, ls="-", label="JULES")
        a.plot(df2[v].index.to_pydatetime(), df2[v].rolling(window=5).mean(),
               c=colours[2], lw=1.5, ls="-", label="CABLE")

        x.bar(df_met.index, df_met["Rainf"], alpha=0.3, color="black")



        fig_label = "%s" % (next(labels_gen))
        a.text(0.02, 0.95, fig_label,
                transform=a.transAxes, fontsize=14, verticalalignment='top',
                bbox=props)

    ax1.set_ylim(0, 17)
    ax2.set_ylim(0, 150)
    ax3.set_ylim(0, 4)
    ax4.set_ylim(0, 0.4)
    axx1.set_yticks([0, 15, 30])
    axx2.set_yticks([0, 15, 30])
    axx3.set_yticks([0, 15, 30])
    axx4.set_yticks([0, 15, 30])
    labels = ["GPP (g C m$^{-2}$ d$^{-1}$)", "LE (W m$^{-2}$)", \
              "LAI (m$^{2}$ m$^{-2}$)", "$\Theta$ (m$^{3}$ m$^{-3}$)"]
    for a, l in zip(axes, labels):
        a.set_ylabel(l, fontsize=12)

    axx1.set_ylabel("Rainfall (mm d$^{-1}$)", fontsize=12, position=(0.5, 0.0))


    from matplotlib.ticker import MaxNLocator
    ax1.yaxis.set_major_locator(MaxNLocator(5))
    ax2.yaxis.set_major_locator(MaxNLocator(5))
    ax3.yaxis.set_major_locator(MaxNLocator(5))
    ax4.yaxis.set_major_locator(MaxNLocator(5))
    #axx1.yaxis.set_major_locator(MaxNLocator(3))
    #axx2.yaxis.set_major_locator(MaxNLocator(3))

    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax3.get_xticklabels(), visible=False)
    ax1.legend(numpoints=1, loc="best", frameon=False)


    for a in axes:
        #a.set_xlim([datetime.date(2002,10,1), datetime.date(2003, 4, 1)])
        #a.set_xlim([datetime.date(2002,12,1), datetime.date(2003, 5, 1)])
        a.set_xlim([datetime.date(2006,11,1), datetime.date(2007, 4, 1)])
        #a.set_xlim([datetime.date(2012,7,1), datetime.date(2013, 8, 1)])
        #a.set_xlim([datetime.date(2006,11,1), datetime.date(2007, 4, 1)])

    if plot_fname is None:
        plt.show()
    else:
        #fig.autofmt_xdate()
        fig.savefig(plot_fname, bbox_inches='tight', pad_inches=0.1)

def read_cable_file(fname, type=None):
    import xarray as xr

    if type == "CABLE":
        #vars_to_keep = ['weighted_psi_soil','psi_leaf','psi_soil','SoilMoist',\
        #                'LAI','GPP','Qle','TVeg','NEE']
        vars_to_keep = ['LAI','GPP','Qle','TVeg','NEE','SoilMoist']
    elif type == "FLUX":
        vars_to_keep = ['GPP','Qle']
    elif type == "MET":
        vars_to_keep = ['Rainf']

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

    if type == "CABLE":
        ds = ds[vars_to_keep].squeeze(dim=["x","y"], drop=True)
        df = pd.DataFrame(ds['GPP'].values[:], columns=['GPP'])
        df['LAI'] = ds['LAI'].values[:]
        df['NEE'] = ds['NEE'].values[:]
        df['Qle'] = ds['Qle'].values[:]
        df['TVeg'] = ds['TVeg'].values[:]


        zse = np.array([.022, .058, .154, .409, 1.085, 2.872])

        frac1 = zse[0] / sum(zse)
        frac2 = zse[1] / sum(zse)
        frac3 = zse[2] / sum(zse)
        frac4 = zse[3] / sum(zse)
        frac5 = zse[4] / sum(zse)
        frac6 = zse[5] / sum(zse)

        df['theta_weight'] = (ds['SoilMoist'][:,0] * frac1) + \
                             (ds['SoilMoist'][:,1] * frac2) + \
                             (ds['SoilMoist'][:,2] * frac3) + \
                             (ds['SoilMoist'][:,3] * frac4) + \
                             (ds['SoilMoist'][:,4] * frac5) + \
                             (ds['SoilMoist'][:,5] * frac6)


    elif type == "MET":
        ds = ds[vars_to_keep].squeeze(dim=["x","y"], drop=True)
        df = pd.DataFrame(ds['Rainf'], columns=['Rainf'])
    elif type == "FLUX":
        ds = ds[vars_to_keep].squeeze(dim=["x","y"], drop=True)
        df = pd.DataFrame(ds['Qle'].values[:], columns=['Qle'])
        df['GPP'] = ds['GPP'].values[:]

    start = reference_date.strip().split(" ")[0].replace("-","/")
    df['dates'] = pd.date_range(start=start, periods=len(df), freq=freq)
    df = df.set_index('dates')

    return df

def resample_timestep(df, type=None):

    UMOL_TO_MOL = 1E-6
    MOL_C_TO_GRAMS_C = 12.0
    SEC_2_HLFHOUR = 1800.
    SEC_2_HOUR = 3600.
    HOUR_2_DAY = 24.
    KG_2_G = 1000.

    if type == "JULES":
        # kg/m2/s -> g/C/60min
        df['GPP'] *= KG_2_G * SEC_2_HOUR

        # kg/m2/s -> mm/60min
        df['TVeg'] *= SEC_2_HOUR
        df['Evap'] *= SEC_2_HOUR

        method = {'GPP':'sum', 'TVeg':'sum', "Qle":"mean", "LAI":"mean",
                  "theta_weight":"mean"}

    elif type == "CABLE":
        # umol/m2/s -> g/C/60min
        #df['GPP'] *= UMOL_TO_MOL * MOL_C_TO_GRAMS_C * SEC_2_HLFHOUR
        df['GPP'] *= UMOL_TO_MOL * MOL_C_TO_GRAMS_C * SEC_2_HOUR

        # kg/m2/s -> mm/60min
        #df['TVeg'] *= SEC_2_HLFHOUR
        df['TVeg'] *= SEC_2_HOUR

        method = {'GPP':'sum', 'TVeg':'sum', "Qle":"mean", "LAI":"mean",
                 "theta_weight":"mean"}

    elif type == "FLUX":
        # umol/m2/s -> g/C/60min
        #df['GPP'] *= UMOL_TO_MOL * MOL_C_TO_GRAMS_C * SEC_2_HLFHOUR
        df['GPP'] *= UMOL_TO_MOL * MOL_C_TO_GRAMS_C * SEC_2_HOUR

        method = {'GPP':'sum', "Qle":"mean"}

    elif type == "MET":
        # kg/m2/s -> mm/60min
        df['Rainf'] *= SEC_2_HOUR

        method = {'Rainf':'sum'}

    df = df.resample("D").agg(method)

    return df

def read_jules_file(fname, type=None):

    vars_to_keep = ['GPP','Qle','Qh','LAI','TVeg','Evap','SoilMoist']

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

    zse = np.array([.1, 0.25, 0.65, 2.0]) * 1000 # m to mm

    frac1 = zse[0] / sum(zse)
    frac2 = zse[1] / sum(zse)
    frac3 = zse[2] / sum(zse)
    frac4 = zse[3] / sum(zse)


    df['theta_weight'] = (ds['SoilMoist'][:,0] / zse[0] * frac1) + \
                         (ds['SoilMoist'][:,1] / zse[1] * frac2) + \
                         (ds['SoilMoist'][:,2] / zse[2] * frac3) + \
                         (ds['SoilMoist'][:,3] / zse[3] * frac4)

    start = reference_date.strip().split(" ")[0].replace("-","/")
    df['dates'] = pd.date_range(start=start, periods=len(df), freq=freq)
    df = df.set_index('dates')

    return df

def label_generator(case='lower', start='', end=''):
    choose_type = {'lower': string.ascii_lowercase,
                   'upper': string.ascii_uppercase}
    generator = ('%s%s%s' %(start, letter, end) for letter in choose_type[case])

    return generator

if __name__ == "__main__":

    site = "TumbaFluxnet"
    output_dir = "outputs"
    met_dir = "/Users/xj21307/research/OzFlux"
    flux_dir = "/Users/xj21307/research/OzFlux"

    parser = OptionParser()
    parser.add_option("-a", "--fname1", dest="fname1",
                      action="store", help="filename",
                      type="string",
                      default=os.path.join(output_dir, "original_out.nc"))
    parser.add_option("-b", "--fname2", dest="fname2",
                      action="store", help="filename",
                      type="string",
                      default=os.path.join(output_dir, "%s_out.nc" % site))
    parser.add_option("-p", "--plot_fname", dest="plot_fname", action="store",
                      help="Benchmark plot filename", type="string")
    (options, args) = parser.parse_args()


    flux_fname = os.path.join(flux_dir, "TumbarumbaOzFlux2.0_flux.nc")
    met_fname = os.path.join(met_dir, "TumbarumbaOzFlux2.0_met.nc")

    fname1 = "/Users/xj21307/research/JULES/runs/roses/AU_Tum/outputs/AU_Tum/local_AU-Tum_fluxnet2015_AU-Tum.AU-Tum.nc"
    #fname2 = "outputs/profitmax_tumba.nc"
    #fname2 = "/Users/xj21307/research/CABLE/runs/ozflux_sites_profitmax_hydraulics_paper/outputs/standard_tumba.nc"
    fname2 = "/Users/xj21307/research/CABLE/runs/quick_jules_comp_tumba/outputs/standard_tumba.nc"
    #fname2 = "/Users/xj21307/research/CABLE/runs/ozflux_sites_profitmax_hydraulics_paper/outputs/profitmax_tumba.nc"

    #odir = "/Users/mdekauwe/Dropbox/Documents/papers/Future_euc_drought_paper/figures/figs/"
    odir = "/Users/xj21307/Desktop/"
    plot_fname = os.path.join(odir, "tumba.pdf")
    main(site, met_fname, flux_fname, fname1, fname2, plot_fname=plot_fname)
