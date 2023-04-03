#!/usr/bin/env python

"""
A series of summary stats comparing model and obs

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (04.07.2014)"
__email__ = "mdekauwe@gmail.com"

import numpy as np

def filter_nan(m,o):
    """
    Remove the data from model and obs arrays whereever the observed data
    contains nan's. Used by all functions.
    """
    data = np.array([m,o])
    data = np.transpose(data)
    data = data[~np.isnan(data).any(1)]
    return data[:,0],data[:,1]

def pc_bias(m,o):
    """ percentage bias """
    #m,o = filter_nan(m,o)
    return 100.0*sum(m-o)/sum(o)

def apb(m,o):
    """absolute percent bias"""
    #m,o = filter_nan(m,o)
    return 100.0*sum(abs(m-o))/sum(o)

def rmse(m,o):
    """
    -same units as forecast, obs
    -sensitive to large errors)
    """
    #m,o = filter_nan(m,o)
    return np.sqrt(np.mean((m-o)**2))

def mae(m,o):
    """ mean absolute error """
    #m,o = filter_nan(m,o)
    return np.mean(abs(m-o))

def bias(m,o):
    """
    Mean error or bias

    positive: forecast is on average above the obs (too warm)
    negative: forecast is on average below the obs (too cold)
    """
    #m,o = filter_nan(m,o)
    return np.mean(m-o)

def nash_sutcliffe(m,o):
    """ Nash-Sutcliffe model efficiency (ME) coefficient
    ME = 1 : perfect match
    ME = 0 : model predictions are as accurate as the mean of the observed data
    ME < 0 : the observed mean is a better predictor than the model
    """
    #m,o = filter_nan(m,o)
    return 1 - sum((m-o)**2)/sum((o-np.mean(o))**2)

def correlation(m,o):
    #m,o = filter_nan(m,o)
    return np.corrcoef(o, s)[0,1]

def willmott_agreement_indx(m,o):
    """
    computes the Index of Agreement between sim and obs, with treatment of
    missing values.

    The Index of Agreement (d) developed by Willmott (1981) as a standardized
    measure of the degree of model prediction error and varies between 0 and 1.
    A value of 1 indicates a perfect match, and 0 indicates no agreement at
    all (Willmott, 1981).

    The index of agreement can detect additive and proportional differences in
    the observed and simulated means and variances; however, it is overly
    sensitive to extreme values due to the squared differences
    (Legates and McCabe, 1999).

    """
    #m,o = filter_nan(m,o)
    d = (1 - ((np.sum((m-o)**2)) / np.sum((np.abs(m-np.mean(o)) +
               np.abs(m-np.mean(o)))**2)))
    return d
