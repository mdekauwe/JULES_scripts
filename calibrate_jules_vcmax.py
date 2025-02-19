#!/usr/bin/env python

"""
Calibrate the JULES physiology params

NB. you need to change the vcmax flag (l_trait_phys=False) by hand yourself,
but we should change this

We're doing this for the first half of the timeseries 2000-8 and between May
and Sept

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (18.02.2025)"
__email__ = "mdekauwe@gmail.com"

import sys
import pandas as pd
import numpy as np
import datetime
import os
import glob
import subprocess
import shutil
import xarray as xr
import re
import os
import shutil
import subprocess
import numpy as np
from summary_stats import rmse, bias, nash_sutcliffe, joint_rmse
from scipy.stats import linregress

def calibrate_variable(fname_flux, output_nc_fname, output_dir,
                       jules_runme_path, var_string, var_name, var_value,
                       kmax, neff, zse, move_file=False):
    #print(f"{var_string}: {var_value}")
    if var_string == "kmax_pft_io":
        print(f"{var_string}: {var_value*1e9}")
    else:
        print(f"{var_string}: {var_value}")

    update_variable(output_file, var_string, var_value)

    # Call the bash script
    try:
        subprocess.run(['bash', jules_runme_path], check=True, \
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running the bash script: {e}")
        print(f"Error message:\n {e.stderr.decode()}")


    # Plot the outputs
    #script_to_run = "plot_timeseries_comparison_CTL_vs_hydraulics.py"
    script_to_run = "plot_calibration_result.py"

    if var_string == "kmax_pft_io":
        argument = f"FR_Pue_{var_name}_equals_{var_value*1e9}.png"
    else:
        argument = f"FR_Pue_{var_name}_equals_{var_value}.png"

    command = ['python', script_to_run, argument]

    try:
        subprocess.run(command, check=True)
        #print("Script executed successfully.")
    except subprocess.CalledProcessError as e:
        print("Error while executing script:", e.stderr.decode())

    (rmse_gpp, rmse_qle, rmse_joint,
     nash_sutcliffe_gpp, nash_sutcliffe_qle,
     r_value_gpp, r_value_qle,
     bias_gpp, bias_qle) = calc_stats(fname_flux, output_nc_fname)

    print("RMSE GPP, RMSE Qle, RMSE jt")
    print(rmse_gpp, rmse_qle, rmse_joint)
    print("NS_GPP, NS_Qle, r_GPP, r_Qle, bias_GPP, bias_Qle")
    print(nash_sutcliffe_gpp,nash_sutcliffe_qle, r_value_gpp, r_value_qle,
          bias_gpp, bias_qle)
    print("\n")

    if move_file:
        # Move the output file
        out_fname = os.path.join(output_dir,
                                 f"FR_Pue_{var_name}_equals_{var_value:.6f}.nc")
        try:
            shutil.move(output_nc_fname, out_fname)
            print(f"File '{output_nc_fname}' has been moved to '{out_fname}'")
        except Exception as e:
            print(f"Error while moving the file: {e}")

    return (rmse_gpp, rmse_qle, rmse_joint)

def calc_stats(flux_fname, model_fname):

    df_mod,freq = read_jules_file(model_fname, zse, type="hydraulics")
    df_obs = read_flux_file(fname_flux)

    df_mod = resample_timestep(df_mod, type="JULES", time_step=freq)
    df_obs = resample_timestep(df_obs, type="FLUX", time_step=freq)

    # Define date range - use first half of the timeseries
    start_date = datetime.date(2000, 1, 1)
    end_date = datetime.date(2008, 1, 1)

    #df_mod = df_mod[(df_mod.index >= pd.Timestamp(start_date)) & \
    #                (df_mod.index < pd.Timestamp(end_date))]

    #df_obs = df_obs[(df_obs.index >= pd.Timestamp(start_date)) & \
    #                (df_obs.index < pd.Timestamp(end_date))]

    # Filter only for May to Sept
    df_mod = df_mod[(df_mod.index.month >= 5) & (df_mod.index.month < 10)]
    df_obs = df_obs[(df_obs.index.month >= 5) & (df_obs.index.month < 10)]

    m1 = df_mod["GPP"].values
    m2 = df_mod["Qle"].values
    o1 = df_obs["GPP"].values
    o2 = df_obs["Qle"].values

    mask = ~np.isnan(m1) & ~np.isnan(o1)
    m1 = m1[mask]
    o1 = o1[mask]

    mask = ~np.isnan(m2) & ~np.isnan(o2)
    m2 = m2[mask]
    o2 = o2[mask]

    rmse_gpp = rmse(m1, o1)
    rmse_qle = rmse(m2, o2)
    rmse_joint = joint_rmse(rmse_gpp, len(m1), rmse_qle, len(m2))
    nash_sutcliffe_gpp = nash_sutcliffe(m1, o1)
    nash_sutcliffe_qle = nash_sutcliffe(m2, o2)
    (slope_gpp, intercept_gpp, r_value_gpp,
     p_value_gpp, std_err_gpp) = linregress(m1, o1)
    (slope_qlq, intercept_qle, r_value_qle,
     p_value_qle, std_err_qle) = linregress(m2, o2)
    bias_gpp = bias(m1, o1)
    bias_qle = bias(m2, o2)

    return (rmse_gpp, rmse_qle, rmse_joint, nash_sutcliffe_gpp,
            nash_sutcliffe_qle, r_value_gpp, r_value_qle, bias_gpp, bias_qle)


def read_jules_file(fname, zse, type=None):
    # Define variables to keep based on the type
    vars_to_keep = ['GPP', 'Qle', 'Qh', 'LAI', 'TVeg', 'Evap',\
                    'SoilMoist', 'fsmc']
    if type == "hydraulics":
        vars_to_keep += ['psi_root_zone_pft', 'psi_leaf_pft']

    # Load dataset
    ds = xr.open_dataset(fname, decode_times=False)

    # Determine frequency based on time step
    time_jump = int(ds.time[1].values) - int(ds.time[0].values)
    freq_map = {3600: "H", 1800: "30min", 10800: "3H"}
    freq = freq_map.get(time_jump)

    if not freq:
        raise ValueError("Unsupported time step detected in the dataset.")

    # Extract units and reference date from time attributes
    _, reference_date = ds.time.attrs['units'].split('since')

    ds = ds[vars_to_keep].squeeze(dim=["x","y"], drop=True)
    df = pd.DataFrame(ds['GPP'].values[:], columns=['GPP'])

    df['Qle'] = ds['Qle'].values[:]
    df['Qh'] = ds['Qh'].values[:]
    df['LAI'] = ds['LAI'].values[:]
    df['TVeg'] = ds['TVeg'].values[:]
    df['Evap'] = ds['Evap'].values[:]
    df['beta'] = ds['fsmc'][:,0] # just EBF (But only using one pft)
    if type == "hydraulics":
        df['psi_rootzone'] = ds['psi_root_zone_pft'][:,0] # just EBF (But only using one pft)
        df['psi_can'] = ds['psi_leaf_pft'][:,0] # just EBF (But only using one pft)


    df['SoilMoist1'] = ds['SoilMoist'].values[:,0] / (zse[0] * 1000) # divided by soil thickness (m) x water density (kg/m3 = 1000)
    df['SoilMoist2'] = ds['SoilMoist'].values[:,1] / (zse[1] * 1000)
    df['SoilMoist3'] = ds['SoilMoist'].values[:,2] / (zse[2] * 1000)
    df['SoilMoist4'] = ds['SoilMoist'].values[:,3] / (zse[3] * 1000)

    #def compute_weighted_avg(soil_moist, weights):
    #    return (soil_moist / zse[:len(weights)] * weights).sum(axis=1)
    #total_thickness = zse.sum()
    #df['theta_weight'] = compute_weighted_avg(ds['SoilMoist'], zse / (total_thickness * 1000)


    # Add date index
    start = reference_date.strip().split()[0].replace("-", "/")
    df['dates'] = pd.date_range(start=start, periods=len(df), freq=freq)
    df.set_index('dates', inplace=True)

    return df, freq

def read_flux_file(fname):

    # Read the CSV file
    df = pd.read_csv(fname)

    # Convert timestamp to datetime and set as index
    df['date'] = pd.to_datetime(df['TIMESTAMP_START'], format='%Y%m%d%H%M')
    df.set_index('date', inplace=True)

    # Drop unnecessary columns
    df.drop(columns=['TIMESTAMP_START', 'TIMESTAMP_END'], inplace=True)

    # Add time-based columns
    df['year'] = df.index.year
    df['doy'] = df.index.dayofyear
    df['hod'] = df.index.hour

    # Rename columns for consistency
    rename_map = {
        "GPP_NT_VUT_REF": "GPP",
        "LE_F_MDS": "Qle",
        "H_F_MDS": "Qh",
        "P_ERA": "Precip",
    }
    df.rename(columns=rename_map, inplace=True)

    # Keep only relevant variables
    vars_to_keep = list(rename_map.values())
    df = df[vars_to_keep]

    return df


def resample_timestep(df, type=None, time_step='30min'):

    # Constants
    UMOL_TO_MOL = 1E-6
    MOL_C_TO_GRAMS_C = 12.0
    KG_TO_G = 1000.0

    # Time step conversion factor
    if time_step == '30min':
        SEC_CONVERSION = 1800.0
    elif time_step == 'hour':
        SEC_CONVERSION = 3600.0
    else:
        raise ValueError("Invalid time_step. Choose '30min' or 'hour'.")

    # Transformation functions for common operations
    def apply_conversion(data, conversion_factor):
        return data * conversion_factor

    def handle_negative_values(data, column):
        """Replace negative values in a column with NaN."""
        data.loc[data[column] < 0, column] = np.nan

    # Define transformations and aggregation methods based on the type
    if type == "JULES":
        df['GPP'] = apply_conversion(df['GPP'], KG_TO_G * SEC_CONVERSION)
        df['TVeg'] = apply_conversion(df['TVeg'], SEC_CONVERSION)
        df['Evap'] = apply_conversion(df['Evap'], SEC_CONVERSION)
        handle_negative_values(df, "Qle")
        aggregation_methods = {
            'GPP': 'sum', 'TVeg': 'sum', 'Qle': 'mean', 'LAI': 'mean',
            'SoilMoist1': 'mean', 'SoilMoist2': 'mean', 'SoilMoist3': 'mean',
            'SoilMoist4': 'mean', 'beta': 'min'
        }

    elif type == "CABLE":
        df['GPP'] = apply_conversion(df['GPP'], UMOL_TO_MOL * \
                                     MOL_C_TO_GRAMS_C * SEC_CONVERSION)
        df['TVeg'] = apply_conversion(df['TVeg'], SEC_CONVERSION)
        df['ESoil'] = apply_conversion(df['ESoil'], SEC_CONVERSION)
        df['ECanop'] = apply_conversion(df['ECanop'], SEC_CONVERSION)
        handle_negative_values(df, "Qle")
        aggregation_methods = {
            'GPP': 'sum', 'TVeg': 'sum', 'Qle': 'mean', 'LAI': 'mean',
            'theta_weight': 'mean', 'theta_20': 'mean', 'theta_50': 'mean',
            'beta': 'min', 'ESoil': 'sum', 'ECanop': 'sum',
            'psi_can': 'min', 'psi_rootzone': 'min'
        }

    elif type == "FLUX":
        df['GPP'] = apply_conversion(df['GPP'], UMOL_TO_MOL * \
                                     MOL_C_TO_GRAMS_C * SEC_CONVERSION)
        handle_negative_values(df, "Qle")
        aggregation_methods = {
            'Qle': 'mean', 'GPP': 'sum', 'Precip': 'sum'
        }

    elif type == "MET":
        df['Precip'] = apply_conversion(df['Precip'], SEC_CONVERSION)
        aggregation_methods = {
            'Precip': 'sum', 'Tair': 'max', 'Qair': 'mean',
            'Psurf': 'mean', 'VPD': 'max'
        }

    else:
        raise ValueError("Invalid type provided. Choose from 'JULES', \
                         'CABLE', 'FLUX', or 'MET'.")

    # Resample and aggregate data
    resampled_df = df.resample("D").agg(aggregation_methods)

    return resampled_df


def update_variable(output_file, var_string, value):
    new_values = f"5*{value},"

    with open(output_file, 'r') as file:
        lines = file.readlines()

    with open(output_file, 'w') as file:
        for line in lines:
            if line.strip().startswith(var_string):
                file.write(f"{var_string} = {new_values}\n")  # Replace with new values
            else:
                file.write(line)

def update_timestamps(file_path, new_timestamp):
    with open(file_path, 'r') as file:
        content = file.read()

    # Define regex pattern to match the timestamps
    pattern = r"(main_run_end|spinup_end)='\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}'"

    # Replace timestamps with the new value
    updated_content = re.sub(pattern,
                             lambda m: f"{m.group(1)}='{new_timestamp}'",
                             content)

    with open(file_path, 'w') as file:
        file.write(updated_content)

def update_flag(file_path, parameter, set_value=True):

    new_lines = []
    target_line = f"{parameter}="
    new_value = ".true." if set_value else ".false."

    with open(file_path, "r") as file:
        for line in file:
            # Ignore commented-out lines
            if target_line in line and "!!" not in line:
                line = f"{parameter}={new_value},\n"
            new_lines.append(line)

    with open(file_path, "w") as file:
        file.writelines(new_lines)


if __name__ == "__main__":

    # Calirbation script currently only works with this option, make sure it is
    # set
    flags_fname = 'namelists/jules_vegetation.nml'
    update_flag(flags_fname, "l_trait_phys", set_value=False)

    neff_range = np.linspace(0.3e-3, 1.8e-3, 10)

    fname_flux = "FLX_FR-Pue_FLUXNET2015_FULLSET_HH_2000-2014_2-4.csv"
    output_file = 'namelists/pft_params.nml'
    output_nc_fname = "outputs/FR_Pue/local_FR_Pue_fluxnet2015_FR_Pue.FR_Pue.nc"
    output_dir = "outputs/FR_Pue/"
    jules_runme_path = os.path.join(os.getcwd(), "runMe.sh")


    move_file = False

    # Soil layer thickness (mm)
    #zse = np.array([0.1000,0.2500,0.6500,2.0000]) # 3.0 m
    zse = np.array([0.1,0.15,0.2,0.25]) # 0.7m FC=(0.39-WP=0.17)*700=154 mm

    # Reset initial defaults
    neff = 0.8e-3
    #kmax = 0.6e-9 # default
    kmax = 0.357e-9 # calibrated
    update_variable(output_file, "neff_io", neff)
    update_variable(output_file, "kmax_pft_io", kmax)

    # Shorten the run to 2000-2008 for the calibration
    time_fname = "namelists/timesteps.nml"
    update_timestamps(time_fname, '2009-01-01 00:00:00')

    # Specify the file name
    rmse_file_name = "rmse_neff.txt"

    # Open the file in append mode
    fp = open(rmse_file_name, "a")

    var_string = "neff_io"

    for neff in neff_range:
        (rmse_gpp, rmse_qle,
         rmse_joint) = calibrate_variable(fname_flux, output_nc_fname,
                                          output_dir, jules_runme_path,
                                          var_string, "neff", neff, kmax,
                                          neff, zse, move_file)
        fp.write("%f %f %f %f\n" % (neff, rmse_gpp, rmse_qle, rmse_joint))
    fp.close()

    # Return timestamp file to normal
    time_fname = "namelists/timesteps.nml"
    update_timestamps(time_fname, '2015-01-01 00:00:00')
