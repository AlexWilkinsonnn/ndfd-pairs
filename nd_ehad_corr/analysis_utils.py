"""
Author: Colin Weber (webe1077@umn.edu)
Date: 16 September 2025
Purpose: To collect the various functions used in developing and testing 
the ad hoc ND Ehad correction.
"""

# Imports
import argparse
import os
import h5py
import shutil

import polars as pl
import copy as cp
import numpy as np
import matplotlib.pyplot as plt

from h5_utils import *
from math import sqrt, floor
from lmfit import Parameters, Minimizer, report_fit
from sklearn.linear_model import TheilSenRegressor, RANSACRegressor
from scipy.stats import truncnorm, norm, gaussian_kde
from scipy.optimize import minimize_scalar


def apply_cuts(df, cut_col, cut_min, cut_max):
  ''' Apply cuts to the input dataframe such that events are kept if their 
  value is between cut_min and cut_max in the cut column.
  
  Inputs
  ------
  df : Polars DataFrame
    The dataframe to apply cuts to
    
  cut_col : str
    The column of df to apply the cuts to
    
  cut_min, cut_max : float
    The extrema that set the limits of the cut. Only events with a value 
    between cut_min and cut_max in cut_col are kept.
    
  Returns
  -------
  return_df : Polars DataFrame
    The dataframe with the cuts applied
  '''
  return_df = cp.deepcopy(df.filter(
    (pl.col(cut_col) > cut_min) & 
    (pl.col(cut_col) < cut_max)))
  
  return return_df


def get_binning_scheme(df, column, bin_exponent, cond_column="", 
  cond_column_min=0, cond_column_max=0):
  '''Calculate the binning scheme for a column of a dataframe with n 
  events in it such that each bin has the same number of events. 
  Optionally, only look at events within a certain range based on a 
  second column.
  
  Inputs
  ------
  df : Polars DataFrame object
    The dataframe that contains the data to bin
    
  column : str
    The column in the dataframe to bin with respect to
    
  bin_exponent : float
    The number of events to the bin exponent determines the number of 
    bins. Should be less than 1/2.
    
  cond_column : str (optional, defaults to "")
    A second column that can be used to restrict the dataframe rows to 
    include in the binning.
    
  cond_column_min : float (optional, defaults to 0)
    The minimum value of the conditional column that events should have
    
  cond_column_max : float (optional, defaults to 0)
    The maximum value of the conditional column that events should have
    
  Returns
  -------
  edges : array-like
    A list containing the bin edges for the column that was binned
  '''
  # Optionally apply the conditions to the dataframe
  if cond_column != "":
    scope_df = apply_cuts(df, cond_column, cond_column_min, 
      cond_column_max)
  else:
    scope_df = cp.deepcopy(df)
  
  # Extract the data to bin from the new dataframes
  data = scope_df.collect()[column].sort().fill_null(0)
  
  # Determine the number of bins to use
  n = len(data)
  n_bins = floor(n ** (bin_exponent))
  
  '''# Determine the edges for same stats in each bin
  edges = []
  if(n == 0 or n_bins == 0):
    # If there are no events or bins, set the edges to be 0-120, the 
    # widest possible range of energies.
    edges.append(0)
    edges.append(120.0)
  elif(n_bins == 1):
    # If there is 1 bin, that's straightforward. Use 0 - the max entry
    edges.append(0)
    edges.append(data[-1])  
  else:
    n_events_per_bin = n // n_bins
    # The bottom edge is the lowest value
    edges.append(data[0])
    # The middle edges are set so that each bin has the same number of 
    # events
    for i in range(n_bins - 1):
      last_in = data[n_events_per_bin * (i + 1) - 1]
      first_out = data[n_events_per_bin * (i + 1) + 1]
      edge = (first_out + last_in) / 2
      edges.append(edge)
    # The last edge is the highest value. The remainder will go here
    edges.append(data[-1])'''
  
  # Determine the edges for even binning
  edges = []
  if(n == 0 or n_bins == 0):
    # If there are no events or bins, set the edges to be 0-120, the 
    # widest possible range of energies.
    edges.append(0)
    edges.append(120.0)
  elif(n_bins == 1):
    # If there is 1 bin, that's straightforward. Use 0 - the max entry
    edges.append(0)
    edges.append(data[-1])
  else:
    edges = np.linspace(had_min, had_max, n_bins + 1)
  
  return edges

  
def fill_reco_edges(h5_df, caf_df, n_h5, n_caf, reco_col, bin_exponent, 
  true_col, true_edges):
  '''Calculate the bin edges for the reco column within each bin of the 
  true column.
  
  Inputs
  ------
  h5_df : Polars DataFrame
    The dataframe containing the events to be shifted
    
  caf_df : Polars DataFrame
    The dataframe containing the events that will determine the shift
    
  n_hf, caf_df : int
    The number of events in each dataframe
    
  reco_col : str
    The dataframe column we want to calculate the binning for
    
  bin_exponent : float
    The number of events to the bin exponent determines the number of 
    bins. Should be less than 1/2.
    
  true_col : str
    The dataframe column that we've already binned according to the edges 
    given by 'true_edges'.
    
  true_edges : array-like
    The already-binned edges of the true column
    
  Returns
  -------
  reco_edges : array-like of len(true_edges) - 1
    A list of lists, in which each component list are the edges for a bin 
    of the true column.
  '''
  reco_edges = []
  print_modulo = math.ceil((len(true_edges) - 1) / 100)
  for i in range(len(true_edges) - 1): # One more edge than bins
    if(i % print_modulo == 0):
      print(
        "Binning for true bin " + str(i) + " / " + 
        str(len(true_edges) - 1))
  
    if n_h5 < n_caf:
      RecoBinScheme = get_binning_scheme(h5_df, reco_col, bin_exponent, 
        true_col, true_edges[i], true_edges[i + 1])
    else:
      RecoBinScheme = get_binning_scheme(caf_df, reco_col, bin_exponent, 
        true_col, true_edges[i], true_edges[i + 1])
    reco_edges.append(RecoBinScheme)
    
  return reco_edges


def get_histogram(df, col, edges, cond_col="", cond_min=0, cond_max=0):
  '''Calculate the histogram bin heights for a column from the input 
  dataframe, with the given edges. Optionally, only look at events within 
  a certain range based on a second column. Return the statistical 
  uncertainties for each bin as well
  
  Inputs
  ------
  df : Polars DataFrame object
    The dataframe that contains the data to put in a histogram
    
  col : str
    The column in the dataframe to that contains the data
    
  edges : array-like
    The edges to use when making the histogram
    
  cond_col: str (optional, defaults to "")
    A second column that can be used to restrict the dataframe rows to 
    include in the histogram.
    
  cond_min : float (optional, defaults to 0)
    The minimum value of the conditional column that events should have
    
  cond_max : float (optional, defaults to 0)
    The maximum value of the conditional column that events should have
    
  Returns
  -------
  histo : array-like of len(edges) - 1
    A list containing the histo bin heights for the column that was binned
    
  unc : array-like of len(edges) - 1
  '''
  # Optionally apply the conditions to the dataframe
  if cond_col != "":
    scope_df = apply_cuts(df, cond_col, cond_min, cond_max)
  else:
    scope_df = cp.deepcopy(df)
  
  # Get the histogram  
  histo = np.histogram(scope_df.collect()[col], bins=edges)[0]
  unc = [sqrt(histo[i]) for i in range(len(histo))]
  
  return histo, unc
  
def get_stats(df, col, cond_col="", cond_min=0, cond_max=0):
  '''Calculate the mean and standard deviation of the df in col. 
  Optionally, only look at events within a certain range based on a s
  econd column. Return the statistical uncertainties for each bin as well
  
  Inputs
  ------
  df : Polars DataFrame object
    The dataframe that contains the data to put in a histogram
    
  col : str
    The column in the dataframe to that contains the data
    
  cond_col: str (optional, defaults to "")
    A second column that can be used to restrict the dataframe rows to 
    include in the histogram.
    
  cond_min : float (optional, defaults to 0)
    The minimum value of the conditional column that events should have
    
  cond_max : float (optional, defaults to 0)
    The maximum value of the conditional column that events should have
    
  Returns
  -------
  cv : float
    Some metric for the central value of the data (mean, median, etc.)
    
  stddev : float
    The standard deviation of the data
  '''
  # Optionally apply the conditions to the dataframe
  if cond_col != "":
    scope_df = apply_cuts(df, cond_col, cond_min, cond_max)
  else:
    scope_df = cp.deepcopy(df)
  
  # Get the stats
  data = np.array(scope_df.collect()[col])
  '''kde = gaussian_kde(data)
  objective = lambda x: -kde(x)
  bounds = (min(data), max(data))
  result = minimize_scalar(objective, bounds=bounds)
  if result.success:
    mean = result.x
  else:
    mean = np.median(data)'''
  cv = np.median(data)
  stddev = np.std(data)

  return cv, stddev


def truncgauss_objective_w_xerr(params, x, y, x_unc, y_unc):
  ''' Define the quantity to minimize to be the difference between a 
  truncated Gaussian model given by params and the y data, 
  weighted according to the uncertainties given by x_unc and y_unc.
  
  Inputs
  ------
  params : lmfit Parameters object
    The object holding the truncated Gaussian parameters defining the 
    model
    
  x, y : array-like
    The data forming the ordered pairs to perform the fit to
    
  x_unc, y_unc : array-like
    The uncertainties in x and y, respectively
    
  Returns
  -------
  loss : array-like
    The differences between the model and the data at each point
  '''
  model = truncnorm.pdf(x, -params['cen'] / params['sig'], (max(y) - 
    params['cen']) / params['sig'], loc=params['cen'], 
    scale=params['sig'])
  dmodel_dx = np.gradient(model) / np.gradient(x)
  dmodel = np.sqrt(y_unc**2 + (x_unc*dmodel_dx)**2)
  
  loss = (y - model) / dmodel
  
  return loss
  
  
def norm_objective_w_xerr(params, x, y, x_unc, y_unc):
  ''' Define the quantity to minimize to be the difference between a 
  Gaussian model given by params and the y data, weighted according to 
  the uncertainties given by x_unc and y_unc.
  
  Inputs
  ------
  params : lmfit Parameters object
    The object holding the truncated Gaussian parameters defining the 
    model
    
  x, y : array-like
    The data forming the ordered pairs to perform the fit to
    
  x_unc, y_unc : array-like
    The uncertainties in x and y, respectively
    
  Returns
  -------
  loss : array-like
    The differences between the model and the data at each point
  '''
  model = norm.pdf(x, loc=params['cen'], scale=params['sig'])
  dmodel_dx = np.gradient(model) / np.gradient(x)
  dmodel = np.sqrt(y_unc**2 + (x_unc*dmodel_dx)**2)
  
  loss = (y - model) / dmodel
  
  return loss


def lin_objective_w_xerr(params, x, y, x_unc, y_unc):
  ''' Define the quantity to minimize to be the difference between a 
  linear model given by params and the y data, 
  weighted according to the uncertainties given by x_unc and y_unc.
  
  Inputs
  ------
  params : lmfit Parameters object
    The object holding the linear parameters defining the model
    
  x, y : array-like
    The data forming the ordered pairs to perform the fit to
    
  x_unc, y_unc : array-like
    The uncertainties in x and y, respectively
    
  Returns
  -------
  loss : array-like
    The differences between the model and the data at each point
  '''
  model = params['slope'] * x + params['int']
  dmodel_dx = np.gradient(model) / np.gradient(x)
  dmodel = np.sqrt(y_unc**2 + (x_unc*dmodel_dx)**2)
  
  loss = (model - y) / dmodel
  
  return loss
  
  
def quad_objective_w_xerr(params, x, y, x_unc, y_unc):
  ''' Define the quantity to minimize to be the difference between a 
  quadratic model given by params and the y data, 
  weighted according to the uncertainties given by x_unc and y_unc.
  
  Inputs
  ------
  params : lmfit Parameters object
    The object holding the linear parameters defining the model
    
  x, y : array-like
    The data forming the ordered pairs to perform the fit to
    
  x_unc, y_unc : array-like
    The uncertainties in x and y, respectively
    
  Returns
  -------
  loss : array-like
    The differences between the model and the data at each point
  '''
  model = params['a'] * x**2 + params['b'] * x + params['c']
  dmodel_dx = np.gradient(model) / np.gradient(x)
  dmodel = np.sqrt(y_unc**2 + (x_unc*dmodel_dx)**2)
  
  loss = (model - y) / dmodel
  
  return loss


def fit_trunc_gaussian(x, y, x_unc, y_unc, true_min, true_max, data_type,
  var, debug=False):
  ''' Fit a Gaussian to the data given by x and y, with uncertainties.
  
  Inputs
  ------
  x, y : array-like
    The data forming the ordered pairs to perform the fit to
    
  x_unc, y_unc : array-like
    The uncertainties in x and y, respectively
    
  true_min, true_max : float
    The extrema of the true bin that is being fit. Used only for plot 
    labels.
    
  data_type, var : str
    The data type and variable that is being plotted. Used only for plot 
    labelling.
    
  debug : bool, defaults to False
    If true, show a plot with the fit results overlaid on the data
    
  Returns
  -------
  mean, stddev : float
    The parameters of the truncated Gaussian fit
    
  valid : bool
    If the fit failed, return false. We define a failed fit as one in 
    which the uncertainty on the center or sigma is greater than 100%, if 
    the best fit value for sigma is negative, if the uncertainties 
    could not be estimated, or if the fitter's 'success' parameter is 
    false
  '''
  # create a set of Parameters
  params = Parameters()
  params.add('cen', value=np.median(x), min=0.0)
  params.add('sig', value=np.std(x), min=0.0)
  
  # We are fitting to a PDF, so scale the y and y_unc variables
  integral = sum(np.multiply(y, np.multiply(x_unc, 2.0)))
  y_norm = np.divide(y, integral)
  y_unc_norm = np.divide(y_unc, integral)
  n_data = len(y_unc_norm)
  # Use the following block to either remove points with 0 events, or set 
  # their uncertainties to 1 (default is 0). This must be done because 
  # the fitter calculates an undefined number if consecutive points both 
  # have 0 events and 0 uncertainty.
  '''for i in range(n_data):
    if y_unc_norm[n_data - i - 1] == 0:
      x = np.delete(x, n_data - i - 1)
      y_norm = np.delete(y_norm, n_data - i - 1)
      x_unc = np.delete(x_unc, n_data - i - 1)
      y_unc_norm = np.delete(y_unc_norm, n_data - i - 1)
      y_unc_norm[n_data - i - 1] = 1'''

  # This chunk for diagnosing why fitter was failing
  '''data = [x, y_norm, x_unc, y_unc_norm]
  columns = ["x", "y", "x_unc", "y_unc"]
  save_df = pl.DataFrame(data, schema=columns)
  save_df.write_csv("{}_{:.4f}_{:.4f}.csv".format(data_type, true_min, true_max))'''
  
  # Perform the fit
  mini = Minimizer(truncgauss_objective_w_xerr, params, 
    fcn_args=(x, y_norm, x_unc, y_unc_norm))
  result = mini.minimize()
  
  # Assess the validity of the fit
  # The fit fails if errors can't be computed, the result was not deemed 
  # a success by the fitter, or the reduced chi squared metric is too far 
  # from 1.
  valid = True
  if((not result.errorbars) or (not result.success) or 
    (result.params['cen'].stderr is None)):
    valid = False
  elif(result.redchi > 1.99 or result.redchi < 0.01):
    valid = False
    
  # Debug the fit
  if(debug):
    report_fit(result)
    x_fine = np.linspace(min(x), max(x), num=500)
    bestfit = truncnorm.pdf(x_fine, 
      -result.params['cen'] / result.params['sig'], 
      (max(y_norm) - result.params['cen']) / result.params['sig'], 
      loc=result.params['cen'], scale=result.params['sig'])
    plt.errorbar(x, y_norm, yerr=y_unc_norm, xerr=x_unc, 
      linestyle='None', mec='blue', mfc='blue')
    plt.plot(x_fine, bestfit, color='tab:red')
    plt.text(max(x) * 0.8, max(y), 'Fit Results')
    if(valid):
      plt.text(max(x) * 0.8, max(y) * 0.9, 'Center = {:.2f} +/- {:.2f}'.format(result.params['cen'].value, result.params['cen'].stderr))
      plt.text(max(x) * 0.8, max(y) * 0.8, 'Width = {:.2f} +/- {:.2f}'.format(result.params['sig'].value, result.params['sig'].stderr))
    else:
      plt.text(max(x) * 0.8, max(y) * 0.9, 'Center = {:.2f} +/- ?'.format(result.params['cen'].value))
      plt.text(max(x) * 0.8, max(y) * 0.8, 'Width = {:.2f} +/- ?'.format(result.params['sig'].value))
    plt.title('{:.5f} < Etrue < {:.5f}'.format(true_min, true_max))
    plt.savefig('{}_{}_{:.5f}_{:.5f}_truncGaussfit.png'.format(
      data_type, var, true_min, true_max), dpi=500)
    plt.close()  
    
  return result.params['cen'].value, result.params['sig'].value, valid
  
  
def fit_peak(x, y, x_unc, y_unc, true_min, true_max, data_type,
  var, debug=False):
  ''' Fit a Gaussian to the peak of the data given by x and y, with 
  uncertainties.
  
  Inputs
  ------
  x, y : array-like
    The data forming the ordered pairs to perform the fit to
    
  x_unc, y_unc : array-like
    The uncertainties in x and y, respectively
    
  true_min, true_max : float
    The extrema of the true bin that is being fit. Used only for plot 
    labels.
    
  data_type, var : str
    The data type and variable that is being plotted. Used only for plot 
    labelling.
    
  debug : bool, defaults to False
    If true, show a plot with the fit results overlaid on the data
    
  Returns
  -------
  mean, stddev : float
    The parameters of the truncated Gaussian fit
    
  valid : bool
    If the fit failed, return false. We define a failed fit as one in 
    which the uncertainty on the center or sigma is greater than 100%, if 
    the best fit value for sigma is negative, if the uncertainties 
    could not be estimated, or if the fitter's 'success' parameter is 
    false
  '''  
  # We are fitting to a PDF, so scale the y and y_unc variables
  integral = sum(np.multiply(y, np.multiply(x_unc, 2.0)))
  y_norm = np.divide(y, integral)
  y_unc_norm = np.divide(y_unc, integral)
  n_data = len(y_unc_norm)
  
  # Find the peak
  max_y = max(y_norm)
  peak_bin_idx = np.where(y_norm == max_y)[0][0]
  
  # We only want to fit to the max and the few bins in its vicinity
  pct_keep_had = 0.15
  n_keep = floor(n_data * pct_keep_had)
  min_bin_idx = max(peak_bin_idx - floor(n_keep / 2), 0)
  max_bin_idx = min_bin_idx + n_keep + 1
  x_trim = x[min_bin_idx:max_bin_idx]
  y_norm_trim = y_norm[min_bin_idx:max_bin_idx]
  x_unc_trim = x_unc[min_bin_idx:max_bin_idx]
  y_unc_norm_trim = y_unc_norm[min_bin_idx:max_bin_idx]
      
  # create a set of Parameters
  params = Parameters()
  params.add('cen', value=np.median(x_trim), min=0.0)
  params.add('sig', value=np.std(x_trim), min=0.0)
  
  # Perform the fit
  mini = Minimizer(norm_objective_w_xerr, params, 
    fcn_args=(x_trim, y_norm_trim, x_unc_trim, y_unc_norm_trim))
  result = mini.minimize()
  
  # Assess the validity of the fit
  valid = True
  if((not result.errorbars) or (not result.success) or 
    (result.params['cen'].stderr is None)):
    valid = False
  elif(result.redchi > 1.99 or result.redchi < 0.01):
    valid = False
    
  # Debug the fit
  if(debug):
    report_fit(result)
    x_fine = np.linspace(min(x), max(x), num=500)
    bestfit = norm.pdf(x_fine, loc=result.params['cen'], 
      scale=result.params['sig'])
    plt.errorbar(x, y_norm, yerr=y_unc_norm, xerr=x_unc, 
      linestyle='None', mec='blue', mfc='blue')
    plt.plot(x_fine, bestfit, color='tab:red')
    plt.text(max(x) * 0.8, max(y), 'Fit Results')
    if(valid):
      plt.text(max(x) * 0.8, max(y) * 0.9, 'Center = {:.2f} +/- {:.2f}'.format(result.params['cen'].value, result.params['cen'].stderr))
      plt.text(max(x) * 0.8, max(y) * 0.8, 'Width = {:.2f} +/- {:.2f}'.format(result.params['sig'].value, result.params['sig'].stderr))
    else:
      plt.text(max(x) * 0.8, max(y) * 0.9, 'Center = {:.2f} +/- ?'.format(result.params['cen'].value))
      plt.text(max(x) * 0.8, max(y) * 0.8, 'Width = {:.2f} +/- ?'.format(result.params['sig'].value))
    plt.title('{:.5f} < Etrue < {:.5f}'.format(true_min, true_max))
    plt.savefig('{}_{}_{:.5f}_{:.5f}_truncGaussfit.png'.format(
      data_type, var, true_min, true_max), dpi=500)
    plt.close()  
    
  return result.params['cen'].value, result.params['cen'].stderr, valid
  
  
def fit_line(x, y, x_unc, y_unc, data_type, var, vary_y_int=False, 
  debug=False):
  ''' Fit a line to the data.
  
  Inputs
  ------
  x, y : array-like
    The data forming the ordered pairs to perform the fit to
    
  x_unc, y_unc : array-like
    The uncertainties in x and y, respectively
    
  data_type, var : str
    The data type and variable that's being fit. Used as a label for the 
    plot.  
    
  vary_y_int : bool, defaults to False (optional)
    If False, hold the y-intercept fixed at 0. This is appropriate for p
    rotons and neutrons. If True, allow it to vary. this is appropriate 
    for pions. 
    
  debug : bool
    If true, show a plot with the fit results overlaid on the data
    
  Results
  -------
  slope, int : float
    The results of the fit
  '''
  # create a set of Parameters
  params = Parameters()
  params.add('slope', value=1.0)
  if(vary_y_int):
    params.add('int', value=0.14, vary=True)
  else:
    params.add('int', value=0.0, vary=False)
  
  # Perform the fit
  mini = Minimizer(lin_objective_w_xerr, params, 
    fcn_args=(x, y, x_unc, y_unc))
  result = mini.minimize()
  
  # Debug the fit
  if(debug):
    report_fit(result)
    x_fine = np.linspace(min(x), max(x), 500)
    bestfit = x_fine * result.params['slope'] + result.params['int'] 
    plt.errorbar(x, y, yerr=y_unc, xerr=x_unc, linestyle='None', 
      mfc='blue', mec='blue')
    plt.plot(x_fine, bestfit, color='tab:red')
    plt.text(max(x) * 0.8, max(y) * 0.4, 'Fit Results')
    plt.text(max(x) * 0.8, max(y) * 0.3, 'Slope = {:.2f} +/- {:.2f}'.format(result.params['slope'].value, result.params['slope'].stderr))
    plt.text(max(x) * 0.8, max(y) * 0.2, 'Int = {:.2f} +/- {:.2f}'.format(result.params['int'].value, result.params['int'].stderr))
    plt.title('{} {} linear fit'.format(data_type, var))
    plt.savefig('{}_{}_linearfit.png'.format(data_type, var), dpi=500)
    plt.close()
  
  return result.params['slope'].value, result.params['int'].value
  
  
def perform_lsq_fit(df, col1, col2):
  ''' Fit a line to unbinned data using the least squares method
  
  Inputs
  ------
  df : Polars DataFrame
    The DataFrame containing the data to fit to. Must contain columns 
    given by 'col1' and 'col2'.
    
  col1 : str
    The column in 'df' that contains the x-axis data.
    
  col2 : str
    The column in 'df' that contains the y-axis data.
    
  Returns
  -------
  m, b : float
    The slope and intercept, respectively, of the line of best fit.
  '''
  # Grab the data
  x = np.array(df.collect()[col1])
  y = np.array(df.collect()[col2])
  
  # Section for least-squares (sensitive to outliers, which we have lots 
  # of)
  '''n = np.size(x)
  
  m_x = np.mean(x)
  m_y = np.mean(y)
  
  SS_xy = np.sum(y * x) - n * m_x * m_y # this is n * the cross term
  SS_xx = np.sum(x * x) - n * m_x * m_x # This is n * var(x)
  
  m = SS_xy / SS_xx
  b = m_y - m * m_x'''
  
  # Section for RANSAC Regressor
  ransac = RANSACRegressor()
  ransac.fit(x.reshape(-1, 1), y.reshape(-1, 1))
  line_x = np.array([x[0], x[1]]).reshape(-1, 1)
  line_y = ransac.predict(line_x)
  m = (line_y[1, 0] - line_y[0, 0]) / (line_x[1, 0] - line_x[0, 0])
  b = line_y[0, 0] - m * line_x[0, 0]
  
  return m, b
  
  
def fit_quad(x, y, x_unc, y_unc, data_type, var, debug=False):
  ''' Fit a quadratic to the data.
  
  Inputs
  ------
  x, y : array-like
    The data forming the ordered pairs to perform the fit to
    
  x_unc, y_unc : array-like
    The uncertainties in x and y, respectively
    
  data_type, var : str
    The data type and variable that's being fit. Used as a label for the 
    plot.  

    
  debug : bool
    If true, show a plot with the fit results overlaid on the data
    
  Results
  -------
  a, b, c : float
    The results of the fit
  '''
  # create a set of Parameters
  params = Parameters()
  params.add('a', value=0.0)
  params.add('b', value=1.0, vary=True)
  params.add('c', value=0.0, vary=True)
  
  # Perform the fit
  mini = Minimizer(quad_objective_w_xerr, params, 
    fcn_args=(x, y, x_unc, y_unc))
  result = mini.minimize()
  
  # Debug the fit
  if(debug):
    report_fit(result)
    x_fine = np.linspace(min(x), max(x), 500)
    bestfit = x_fine**2 * result.params['a'] + \
      x_fine * result.params['b'] + result.params['c'] 
    plt.errorbar(x, y, yerr=y_unc, xerr=x_unc, linestyle='None', 
      mfc='blue', mec='blue')
    plt.plot(x_fine, bestfit, color='tab:red')
    plt.text(max(x) * 0.8, max(y) * 0.4, 'Fit Results')
    plt.text(max(x) * 0.8, max(y) * 0.3, 'a = {:.2f} +/- {:.2f}'.format(result.params['a'].value, result.params['a'].stderr))
    plt.text(max(x) * 0.8, max(y) * 0.2, 'b = {:.2f} +/- {:.2f}'.format(result.params['b'].value, result.params['b'].stderr))
    plt.text(max(x) * 0.8, max(y) * 0.1, 'c = {:.2f} +/- {:.2f}'.format(result.params['c'].value, result.params['c'].stderr))
    plt.title('{} {} quadratic fit'.format(data_type, var))
    plt.savefig('{}_{}_quadraticfit.png'.format(data_type, var), dpi=500)
    plt.close()
  
  return result.params['a'].value, result.params['b'].value, \
    result.params['c'].value


def plot_2D_hist(df, col1, col2, edges, data_type, var, line=False, 
  slope=0, yint=0):
  ''' Simple plotting utility
  
  Inputs
  ------
  df : Polars Dataframe
    The dataframe containing the data to use for the pairs
    
  col1, col2 : str
    The dataframe columns to plot. col1 corresponds to the horizontal 
    axis.
    
  edges : array-like
    The edges to use for both axes of the histogram
    
  data_type, var : str
    Strings used for plot labelling
    
  line : bool, defaults to False
    If True, also plot a line on top of the 2D histogram with parameters 
    given by 'slope' and 'yint'.
  '''
  fig, axs = plt.subplots(1, 1, figsize=(8, 8), layout='constrained')
  h, xedges, yedges, img = plt.hist2d(df.collect()[col1], 
    df.collect()[col2], bins=edges, range=[[0, 5], [0, 5]])
  axs.set_xlim(0, 5)
  axs.set_ylim(0, 5)
  plt.colorbar(img, ax=axs)  
  plt.title('{} {}'.format(data_type, var))
  
  # Optionally plot a line
  if(line):
    plt.axline(xy1=(0, yint), slope=slope, color='red', 
    linestyle='solid', marker='')
  
  plt.savefig('{}_{}_2D_hist'.format(data_type, var), dpi=500)    
  plt.close()


def plot_scatter(df, col1, col2, n_entries, data_type, var, line=False, 
  slope=0, yint=0, min_val=0.0, max_val=0.0):
  ''' Simple plotting utility
  
  Inputs
  ------
  df : Polars Dataframe
    The dataframe containing the data to use for the pairs
    
  col1, col2 : str
    The dataframe columns to plot. col1 corresponds to the horizontal 
    axis.
    
  n_entries : int
    The number of dataframe rows to plot.
    
  data_type, var : str
    Strings used for plot labelling
    
  line : bool, defaults to False
    If True, also plot a line on top of the 2D histogram with parameters 
    given by 'slope' and 'yint'.
    
  min_val : float, defaults to 0.0
    The minimum value to show on the axes.
    
  max_val : float, defaults to 0.0
    The maximum value to show on the axes. The default is a placeholder. 
    If it is set to 0, then the max value is calculated to be the max of 
    the data from either column.
  '''
  fig, axs = plt.subplots(1, 1, figsize=(8, 8), layout='constrained')
  plt.scatter(df.collect()[col1][:n_entries], 
    df.collect()[col2][:n_entries], marker='.', linewidths=0.1)
  if(max_val == 0.0):
    max_val = max((max(df.collect()[col1]), max(df.collect()[col2])))
  plt.xlim(min_val, max_val)
  plt.ylim(min_val, max_val)
  plt.title('{} {}'.format(data_type, var))
  
  # Optionally plot a line
  if(line):
    plt.axline(xy1=(0, yint), slope=slope, color='red', 
    linestyle='solid', marker='')
  
  plt.savefig('{}_{}_scatter'.format(data_type, var), dpi=500)    
  plt.close()


def get_pairs(h5_df, caf_df, true_col, reco_col, true_edges, reco_edges, 
  var, debug=False):
  ''' For each bin of true_col as defined by true edges, use the reco 
  edges for the reco col to create a 1D histogram. Then, fit a truncated 
  Gaussian to the histogram and extract the center. The center of the 
  Gaussian for that bin of true_col constitue a pair. The set of these, 
  with their uncertainties, are returned.
  
  Inputs
  ------
  h5_df, caf_df : Polars Dataframe
    The dataframes containing the data to use for the pairs
    
  true_col : str
    The column name to use for the 0th element of the pairs
    
  reco_col : str
    The column name to use for the 1st element of the pairs
    
  true_edges : array-like
    The binning edges for true_col
    
  reco_edges : array-like
    A list of lists of bin edges. There is one list for each bin of 
    true_col
  
  var : str
    The variable being paired. Used only for plotting purposes.
  
  debug : bool, defaults to False
    If true, make plots of the fits
  
  Returns
  -------
  pairs : NumPy array, shape(len(true_edges) - 1, 4)
    The pairs obtained from the fitting procedure
  '''
  TrueCenters = [(true_edges[i + 1] + true_edges[i]) / 2 
    for i in range(len(true_edges) - 1)]
  TrueWidths = [true_edges[i + 1] - true_edges[i] 
    for i in range(len(true_edges) - 1)]
  
  # If debug, print out the 2D histogram we're about to fit to
  if(debug):
    plot_2D_hist(h5_df, true_col, reco_col, true_edges, "LArbath", var)
    plot_2D_hist(caf_df, true_col, reco_col, true_edges, "NDHall", var)
    
  h5_pairs = []
  caf_pairs = []
  # For each of the true bins:
  for i in range(len(true_edges) - 1):
    # Get the reco bin widths and centers
    RecoWidths = np.array(
      [reco_edges[i][j + 1] - reco_edges[i][j] 
      for j in range(len(reco_edges[i]) - 1)])
    RecoCenters = np.array(
      [reco_edges[i][j] + (RecoWidths[j] / 2) 
      for j in range(len(RecoWidths))])
      
    # Get the histogram
    reco_hist_h5, reco_unc_h5 = get_histogram(h5_df, reco_col, 
      reco_edges[i], true_col, true_edges[i], true_edges[i + 1])
    reco_hist_caf, reco_unc_caf = get_histogram(caf_df, reco_col, 
      reco_edges[i], true_col, true_edges[i], true_edges[i + 1])
      
    # Divide by bin widths
    reco_hist_h5 = np.divide(reco_hist_h5, RecoWidths)
    reco_unc_h5 = np.divide(reco_unc_h5, RecoWidths)
    reco_hist_caf = np.divide(reco_hist_caf, RecoWidths)
    reco_unc_caf = np.divide(reco_unc_caf, RecoWidths)
    
    # Perform either truncated Gaussian fits to the entire range, or 
    # Gaussian fits to just the peak region.
    # Section for the truncated Gaussian fit
    h5_mean, h5_unc, h5_valid = fit_trunc_gaussian(RecoCenters, 
      reco_hist_h5, np.array(RecoWidths) / 2.0, reco_unc_h5, 
      true_edges[i], true_edges[i+1], "LArbath", var, debug)
    caf_mean, caf_unc, caf_valid = fit_trunc_gaussian(RecoCenters, 
      reco_hist_caf, np.array(RecoWidths) / 2.0, reco_unc_caf, 
      true_edges[i], true_edges[i+1], "NDHall", var, debug)
    
    # Section for fitting to just the peak  
    '''h5_mean, h5_unc, h5_valid = fit_peak(RecoCenters, 
      reco_hist_h5, np.array(RecoWidths) / 2.0, reco_unc_h5, 
      true_edges[i], true_edges[i+1], "LArbath", var, debug)
    caf_mean, caf_unc, caf_valid = fit_peak(RecoCenters, 
      reco_hist_caf, np.array(RecoWidths) / 2.0, reco_unc_caf, 
      true_edges[i], true_edges[i+1], "NDHall", var, debug)'''
    
    # Append the pair, if the fit succeeded
    if(h5_valid and caf_valid):
      h5_pairs.append([TrueCenters[i], h5_mean, 
        TrueWidths[i] / 2.0, h5_unc])
      caf_pairs.append([TrueCenters[i], caf_mean, 
        TrueWidths[i] / 2.0, caf_unc])
    # Optionally, if we want each point to have some value, can compute a 
    # central value and width in other ways.
    '''else:
      # If one of the fits fails, use the computed mean and standard
      # deviation
      h5_mean, h5_stddev = get_stats(h5_df, reco_col, 
        true_col, true_edges[i], true_edges[i + 1])
      caf_mean, caf_stddev = get_stats(caf_df, reco_col, 
        true_col, true_edges[i], true_edges[i + 1])
      h5_pairs.append([TrueCenters[i], h5_mean, TrueWidths[i] / 2.0, 
        h5_stddev])
      caf_pairs.append([TrueCenters[i], caf_mean, TrueWidths[i] / 2.0, 
        caf_stddev])'''
        
  return np.array(h5_pairs), np.array(caf_pairs)
  
  
def get_shifts(h5_df, caf_df, true_col, reco_col, true_edges, reco_edges, 
  var, debug=False):
  ''' For each bin of true_col as defined by true edges, use the reco 
  edges for the reco col to create a 1D histogram. Then, fit a truncated 
  Gaussian to the histogram and extract the center and width. For each 
  dataset, return the center of the true bin, the center of the Gaussian, 
  the half-width of the true bin, and the width of the Gaussian. If the 
  Gaussian fit fails, return the mean and standard deviation in place of 
  the fit parameters. This therefore calculates a shift for each bin.
  
  Inputs
  ------
  h5_df, caf_df : Polars Dataframe
    The dataframes containing the data to use for the pairs
    
  true_col : str
    The column name to use for the 0th element of the pairs
    
  reco_col : str
    The column name to use for the 1st element of the pairs
    
  true_edges : array-like
    The binning edges for true_col
    
  reco_edges : array-like
    A list of lists of bin edges. There is one list for each bin of 
    true_col
  
  var : str
    The variable being paired. Used only for plotting purposes.
  
  debug : bool, defaults to False
    If true, make plots of the fits
  
  Returns
  -------
  pairs : NumPy array, shape(len(true_edges) - 1, 4)
    The pairs obtained from the fitting procedure
  '''
  TrueCenters = [(true_edges[i + 1] + true_edges[i]) / 2 
    for i in range(len(true_edges) - 1)]
  TrueWidths = [true_edges[i + 1] - true_edges[i] 
    for i in range(len(true_edges) - 1)]
  
  # If debug, print out the 2D histogram we're about to fit to
  if(debug):
    plot_2D_hist(h5_df, true_col, reco_col, true_edges, "LArbath", var)
    plot_2D_hist(caf_df, true_col, reco_col, true_edges, "NDHall", var)
    
  h5_pairs = []
  caf_pairs = []
  # For each of the true bins:
  for i in range(len(true_edges) - 1):
    # Get the reco bin widths and centers
    RecoWidths = np.array(
      [reco_edges[i][j + 1] - reco_edges[i][j] 
      for j in range(len(reco_edges[i]) - 1)])
    RecoCenters = np.array(
      [reco_edges[i][j] + (RecoWidths[j] / 2) 
      for j in range(len(RecoWidths))])
      
    # Get the histogram
    reco_hist_h5, reco_unc_h5 = get_histogram(h5_df, reco_col, 
      reco_edges[i], true_col, true_edges[i], true_edges[i + 1])
    reco_hist_caf, reco_unc_caf = get_histogram(caf_df, reco_col, 
      reco_edges[i], true_col, true_edges[i], true_edges[i + 1])
      
    # Divide by bin widths
    reco_hist_h5 = np.divide(reco_hist_h5, RecoWidths)
    reco_unc_h5 = np.divide(reco_unc_h5, RecoWidths)
    reco_hist_caf = np.divide(reco_hist_caf, RecoWidths)
    reco_unc_caf = np.divide(reco_unc_caf, RecoWidths)
    
    # Fit the data either with a truncated Gaussian across the entire 
    # range, or a Gaussian to just the peak region.
    # Perform the truncated Gaussian fit
    '''h5_mean, h5_width, h5_valid = fit_trunc_gaussian(RecoCenters, 
      reco_hist_h5, np.array(RecoWidths) / 2.0, reco_unc_h5, 
      true_edges[i], true_edges[i+1], "LArbath", var, debug)
    caf_mean, caf_width, caf_valid = fit_trunc_gaussian(RecoCenters, 
      reco_hist_caf, np.array(RecoWidths) / 2.0, reco_unc_caf, 
      true_edges[i], true_edges[i+1], "NDHall", var, debug)'''
    
    # Perform the Gaussian fit to the peak region.  
    h5_mean, h5_width, h5_valid = fit_peak(RecoCenters, 
      reco_hist_h5, np.array(RecoWidths) / 2.0, reco_unc_h5, 
      true_edges[i], true_edges[i+1], "LArbath", var, debug)
    caf_mean, caf_width, caf_valid = fit_peak(RecoCenters, 
      reco_hist_caf, np.array(RecoWidths) / 2.0, reco_unc_caf, 
      true_edges[i], true_edges[i+1], "NDHall", var, debug)
    
    # Append the pair, if the fit succeeded
    if(h5_valid and caf_valid):
      h5_pairs.append([TrueCenters[i], h5_mean, 
        TrueWidths[i] / 2.0, h5_width])
      caf_pairs.append([TrueCenters[i], caf_mean, 
        TrueWidths[i] / 2.0, caf_width])
    else:
      # If one of the fits fails, use the computed mean and standard
      # deviation
      h5_mean, h5_stddev = get_stats(h5_df, reco_col, 
        true_col, true_edges[i], true_edges[i + 1])
      caf_mean, caf_stddev = get_stats(caf_df, reco_col, 
        true_col, true_edges[i], true_edges[i + 1])
      h5_pairs.append([TrueCenters[i], h5_mean, TrueWidths[i] / 2.0, 
        h5_stddev])
      caf_pairs.append([TrueCenters[i], caf_mean, TrueWidths[i] / 2.0, 
        caf_stddev])
        
  return np.array(h5_pairs), np.array(caf_pairs)


def compare_lin_fit(h5_slope, h5_int, caf_slope, caf_int, var, 
  plot_data=True, h5_pairs=[], caf_pairs=[], lim=(0.0, 0.0)):
  ''' Plotting utility for comparing the results of both fits
  
  Inputs
  ------    
  h5_slope, caf_slope, h5_int, caf_int : float
    The best fit parameters for a linear fit to each of the pairs
    
  var : str
    The variable being plotted. Used only for plot labelling
    
  plot_pairs : bool, defaults to True
    If True, plot the data given by 'h5_pairs' and 'caf_pairs'
    
  h5_pairs, caf_pairs : optional, array-like, shape(:, 4)
    The ordered pairs, with uncertainty, for each dataset. The order is 
    (x, y, x_unc, y_unc). If given, these will be plotted.
    
  lim : tuple, defaults to (0.0, 0.0)
    If plot_data is True and lim[1] == 0.0, the max value is determined 
    automatically. Otherwise, use the user-defined limits.
  '''
  # Calculate the best fit data
  if(plot_data):
    bestfit_h5 = h5_pairs[:, 0] * h5_slope + h5_int 
    bestfit_caf = caf_pairs[:, 0] * caf_slope + caf_int
  
    # Make the plots. First plot the data
    f, axs = plt.subplots(1, 1, figsize=(10, 8))
    axs.errorbar(h5_pairs[:, 0], h5_pairs[:, 1], yerr=h5_pairs[:, 3], 
      xerr=h5_pairs[:, 2], linestyle='None', mfc='tab:blue', 
      mec='tab:blue')
    axs.errorbar(caf_pairs[:, 0], caf_pairs[:, 1], yerr=caf_pairs[:, 3], 
      xerr=caf_pairs[:, 2], linestyle='None', mfc='tab:orange', 
      mec='tab:orang')
    # Plot the fits
    axs.plot(h5_pairs[:, 0], bestfit_h5, color='tab:blue', 
      label='LArbath fit')
    axs.plot(caf_pairs[:, 0], bestfit_caf, color='tab:orange', 
      label='NDHall fit')
    if(lim[1] != 0):
      plt.xlim(lim[0], lim[1])
      plt.ylim(lim[0], lim[1])
  else:
    f, axs = plt.subplots(1, 1, figsize=(10, 8))
    # Plot the fits
    axs.axline(xy1=(0, h5_int), slope=h5_slope, color='tab:blue', 
      linestyle='solid', marker='', label='LArbath fit')
    axs.axline(xy1=(0, caf_int), slope=caf_slope, color='tab:orange', 
      linestyle='solid', marker='', label='NDHall fit')
    plt.xlim(lim[0], lim[1])
    plt.ylim(lim[0], lim[1])
  
  # Labels
  axs.set_title('LArbath v. NDHall: {}'.format(var))
  axs.set_xlabel('True')
  axs.set_ylabel('Reco')
  
  # Also show the fit parameters
  axs.text(1.02, 0.5, 
    'LArbath slope = {:.2f}\nLArbath int = {:.2f}\nND Hall slope = {:.2f}\nND Hall int = {:.2f}'.format(
    h5_slope, h5_int, caf_slope, caf_int), transform=axs.transAxes, 
    ha='left', va='center', fontsize=16)
  
  # Legend, save, and close
  plt.legend(loc='upper left')
  plt.savefig('LArbathvNDHall_lin{}'.format(var), dpi=500, 
    bbox_inches='tight')
  plt.close()
  
  
def compare_quad_fit(h5_pairs, h5_a, h5_b, h5_c, caf_pairs, caf_a, caf_b,
  caf_c, var):
  ''' Plotting utility for comparing the results of both fits
  
  Inputs
  ------
  h5_pairs, caf_pairs : array-like, shape(:, 4)
    The ordered pairs, with uncertainty, for each dataset. The order is 
    (x, y, x_unc, y_unc)
    
  h5_a, caf_a, h5_b, caf_b, h5_c, caf_c : float
    The best fit parameters for a quadratic fit to each of the pairs. The 
    fit takes the form ax**2 + bx + c
    
  var : str
    The variable being plotted. Used only for plot labelling
  '''
  x_fine = np.linspace(min(h5_pairs[:, 0]), max(h5_pairs[:, 0]), 500)
  # Calculate the best fit data
  bestfit_h5 = x_fine**2 * h5_a + x_fine * h5_b + h5_c 
  bestfit_caf = x_fine**2 * caf_a + x_fine * caf_b + \
    caf_c
  
  # Make the plots. First plot the data
  f, axs = plt.subplots(1, 1, figsize=(10, 8))
  axs.errorbar(h5_pairs[:, 0], h5_pairs[:, 1], yerr=h5_pairs[:, 3], 
    xerr=h5_pairs[:, 2], linestyle='None', mfc='tab:blue', mec='tab:blue')
  axs.errorbar(caf_pairs[:, 0], caf_pairs[:, 1], yerr=caf_pairs[:, 3], 
    xerr=caf_pairs[:, 2], linestyle='None', mfc='tab:orange', 
    mec='tab:orang')
  # Plot the fits
  axs.plot(x_fine, bestfit_h5, color='tab:blue', 
    label='LArbath fit')
  axs.plot(x_fine, bestfit_caf, color='tab:orange', 
    label='NDHall fit')
  
  # Labels
  axs.set_title('LArbath v. NDHall: {}'.format(var))
  axs.set_xlabel('True')
  axs.set_ylabel('Reco')
  
  # Also show the fit parameters
  axs.text(1.02, 0.5, 
    'LArbath a = {:.2f}\nLArbath b = {:.2f}\nLArbath c = {:.2f}\nND Hall a = {:.2f}\nND Hall b = {:.2f}\nND Hall c = {:.2f}'.format(
    h5_a, h5_b, h5_c, caf_a, caf_b, caf_c), transform=axs.transAxes, 
    ha='left', va='center', fontsize=16)
  
  # Legend, save, and close
  plt.legend(loc='upper left')
  plt.savefig('LArbathvNDHall_quad{}'.format(var), dpi=500, 
    bbox_inches='tight')
  plt.close()


def create_csv(input_h5_df, h5_params, caf_params, outname):
  ''' Write a CSV file using the dataframe, but with selected columns 
  overwritten using the parameters from the HDF5 and CAF file fits. The 
  other hadrons category is fixed to be the corrected reconstructed 
  hadronic energy, minus the corrected reco energies from all the 
  selected particle types.
  
  Inputs
  ------
  input_h5_df : Polars DataFrame
    The data to be overwritten
    
  h5_params, caf_params : dict
    A dictionary where the key is a variable alias, and the value is a 
    list of length 2 in which the 0th item is the slope and the 1st item 
    is the y-intercept for that variable's fit
    
  outname : str
    The name of the output file
  '''
  # Top correction is for if each bin of true gets its own shift.
  '''# Calculate the corrections
  old_Ehad_reco_data = input_h5_df.collect()["ND_Ehad_reco"].to_list()
  Ehad_true_data = input_h5_df.collect()["Ehad_true"].to_list()
  new_Ehad_reco_data = []
  n_data = len(old_Ehad_reco_data)
  print_modulo = math.ceil(n_data / 100)
  for i in range(n_data):
    if(i % print_modulo == 0):
      print("Correcting CSV data point " + str(i) + " / " + str(n_data)) 
    # Find the shift parameters
    h5_Ehad_shifts = h5_params['Ehad']
    caf_Ehad_shifts = caf_params['Ehad']
    n_shifts = h5_Ehad_shifts.shape[0]
    found = False
    for j in range(n_shifts):
      if found == False:
        h5_Ehad_shift = h5_Ehad_shifts[j, :]
        caf_Ehad_shift = caf_Ehad_shifts[j, :]
        if(abs(Ehad_true_data[i] - h5_Ehad_shift[0]) < h5_Ehad_shift[2]):
          h5_mean = h5_Ehad_shift[1]
          h5_width = h5_Ehad_shift[3]
          caf_mean = caf_Ehad_shift[1]
          caf_width = caf_Ehad_shift[3]
          new_val = (caf_width * \
              (old_Ehad_reco_data[i] - h5_mean) / h5_width) + caf_mean
          if(new_val > 0):
            new_Ehad_reco_data.append(new_val)
          else:
            new_Ehad_reco_data.append(old_Ehad_reco_data[i])
          found = True
          continue
        if j == (n_shifts) - 1:
          new_Ehad_reco_data.append(old_Ehad_reco_data[i])  

  output_df = input_h5_df.with_columns([
    pl.lit(pl.Series("ND_Ehad_reco_corr", new_Ehad_reco_data))])'''

  # Correction based on linear fits. The mean shifts as expected for each 
  # bin of true, but the width is rescaled.
  output_df = input_h5_df.with_columns([
    pl.when((pl.col("ND_Ehad_reco") - h5_params['Ehad'][1]) * 
    caf_params['Ehad'][0] / h5_params['Ehad'][0] + 
    caf_params['Ehad'][1] > 0 ).then((pl.col("ND_Ehad_reco") - h5_params['Ehad'][1]) * 
    caf_params['Ehad'][0] / h5_params['Ehad'][0] + 
    caf_params['Ehad'][1]).otherwise(pl.col("ND_Ehad_reco")).alias("ND_Ehad_reco_corr"),
    
    pl.when((pl.col("ND_proton_reco_E") - h5_params['EP'][1]) * 
    caf_params['EP'][0] / h5_params['EP'][0] + 
    caf_params['EP'][1] > 0).then((pl.col("ND_proton_reco_E") - h5_params['EP'][1]) * 
    caf_params['EP'][0] / h5_params['EP'][0] + 
    caf_params['EP'][1]).otherwise(pl.col("ND_proton_reco_E")).alias("ND_proton_reco_E_corr"),
    
    pl.when((pl.col("ND_neutron_reco_E") - h5_params['EN'][1]) * 
    caf_params['EN'][0] / h5_params['EN'][0] + 
    caf_params['EN'][1] > 0).then((pl.col("ND_neutron_reco_E") - h5_params['EN'][1]) * 
    caf_params['EN'][0] / h5_params['EN'][0] + 
    caf_params['EN'][1]).otherwise(pl.col("ND_neutron_reco_E")).alias("ND_neutron_reco_E_corr"),
    
    ((pl.col("ND_pip_reco_E") - h5_params['EPip'][1]) * 
    caf_params['EPip'][0] / h5_params['EPip'][0] + 
    caf_params['EPip'][1]).alias("ND_pip_reco_E_corr"),
    
    ((pl.col("ND_pim_reco_E") - h5_params['EPim'][1]) * 
    caf_params['EPim'][0] / h5_params['EPim'][0] + 
    caf_params['EPim'][1]).alias("ND_pim_reco_E_corr"),
    
    ((pl.col("ND_pi0_reco_E") - h5_params['EPi0'][1]) * 
    caf_params['EPi0'][0] / h5_params['EPi0'][0] + 
    caf_params['EPi0'][1]).alias("ND_pi0_reco_E_corr")])
  # The Other column is set so that the the sum of all the exclusive 
  # hadron reco variables equals the inclusive variable.
  output_df = output_df.with_columns([
    pl.max_horizontal((
      pl.col("ND_Ehad_reco_corr") - pl.col("ND_proton_reco_E_corr") - 
      pl.col("ND_neutron_reco_E_corr") - pl.col("ND_pip_reco_E_corr") - 
      pl.col("ND_pim_reco_E_corr") - pl.col("ND_pi0_reco_E_corr")), 
      pl.lit(0)).alias("ND_other_reco_E_corr"),
      (pl.col("ND_Ehad_reco_corr") + pl.col("ND_Elep_reco")).alias(
      "ND_Ev_reco_corr")])  
  
  # Remove the uncorrected columns  
  output_df_dropped = output_df.drop(["ND_Ehad_reco", "ND_Ev_reco", 
    "ND_proton_reco_E", "ND_neutron_reco_E", "ND_pip_reco_E", 
    "ND_pim_reco_E", "ND_pi0_reco_E", "ND_other_reco_E"])
    
  # Rename the corrected columns, effectively replacing the uncorrected 
  # columns
  output_df_renamed = output_df_dropped.rename({
    "ND_Ev_reco_corr": "ND_Ev_reco",
    "ND_Ehad_reco_corr": "ND_Ehad_reco",
    "ND_proton_reco_E_corr": "ND_proton_reco_E",
    "ND_neutron_reco_E_corr": "ND_neutron_reco_E",
    "ND_pip_reco_E_corr": "ND_pip_reco_E",
    "ND_pim_reco_E_corr": "ND_pim_reco_E",
    "ND_pi0_reco_E_corr": "ND_pi0_reco_E",
    "ND_other_reco_E_corr": "ND_other_reco_E"
    })
  '''output_df_renamed = output_df_dropped.rename({
    "ND_Ev_reco_corr": "ND_Ev_reco",
    "ND_Ehad_reco_corr": "ND_Ehad_reco"
    })'''
  
  output_df_renamed.sink_csv(outname, include_header=True, separator=',')


def correct_Ev_reco(output_h5_file, subfile, h5_params, caf_params):
  ''' Correct down the reconstructed neutrino energy in a subfile of the 
  output HDF5 file so that it looks like a CAF file.
  
  Inputs
  ------
  output_h5_file : HDF5 file
    The file to correct
    
  subfile : str
    The subfile withing output_h5_file to correct
    
  h5_params, caf_params : dict
    A dictionary where the key is a variable alias, and the value is a 
    list of length 2 in which the 0th item is the slope and the 1st item 
    is the y-intercept for that variable's fit
  '''
  # Extract the old data. Hadronic energy is not calculated in the .h5 
  # files, so we calculated it here as the difference between Ev_reco and 
  # Elep_reco
  old_Ev_reco_data = output_h5_file[subfile]['nd_paramreco']['Ev_reco'][:]
  old_Elep_reco_data = output_h5_file[subfile]['nd_paramreco']['Elep_reco'][:]
  old_Ehad_reco_data = old_Ev_reco_data - old_Elep_reco_data
  Ev_true_data = output_h5_file[subfile]['nd_paramreco']['Ev'][:]
  Elep_true_data = output_h5_file[subfile]['nd_paramreco']['LepE'][:]
  Ehad_true_data = Ev_true_data - Elep_true_data
  
  # Apply the correction
  new_Ehad_reco_data = (old_Ehad_reco_data - h5_params['Ehad'][1]) * caf_params['Ehad'][0] / h5_params['Ehad'][0] + caf_params['Ehad'][1]
  new_Ev_reco_data = new_Ehad_reco_data + old_Elep_reco_data
  
  # For if each bin of true gets its own shift
  '''new_Ehad_reco_data = old_Ehad_reco_data - h5_params['Ehad'][1] - \
    caf_params['Ehad'][1] - Ehad_true_data * \
    (h5_params['Ehad'][0] - caf_params['Ehad'][0])
  new_Ehad_reco_data = [max(new_Ehad_reco_data[i], 0.0) 
    for i in range(len(new_Ehad_reco_data))]'''
  '''new_Ehad_reco_data = []
  for i in range(len(old_Ehad_reco_data)):
    # Find the shift parameters
    h5_Ehad_shifts = h5_params['Ehad']
    caf_Ehad_shifts = caf_params['Ehad']
    n_shifts = h5_Ehad_shifts.shape[0]
    found = False
    for j in range(n_shifts):
      if found == False:
        h5_Ehad_shift = h5_Ehad_shifts[j, :]
        caf_Ehad_shift = caf_Ehad_shifts[j, :]
        if(abs(Ehad_true_data[i] - h5_Ehad_shift[0]) < h5_Ehad_shift[2]):
          h5_mean = h5_Ehad_shift[1]
          h5_width = h5_Ehad_shift[3]
          caf_mean = caf_Ehad_shift[1]
          caf_width = caf_Ehad_shift[3]
          new_val = (caf_width * \
              (old_Ehad_reco_data[i] - h5_mean) / h5_width) + caf_mean
          if(new_val > 0):
            new_Ehad_reco_data.append(new_val)
          else:
            new_Ehad_reco_data.append(old_Ehad_reco_data[i])
          found = True
          continue
        if j == (n_shifts) - 1:
          new_Ehad_reco_data.append(old_Ehad_reco_data[i])'''
          
  output_h5_file[subfile]['nd_paramreco']['Ev_reco'][:] = new_Ev_reco_data
  
  
def correct_Epart_reco(output_h5_file, subfile, h5_params, caf_params, 
  Epart):
  ''' Correct down the reconstructed neutrino energy in a subfile of the 
  output HDF5 file so that it looks like a CAF file.
  
  Inputs
  ------
  output_h5_file : HDF5 file
    The file to correct
    
  subfile : str
    The subfile withing output_h5_file to correct
    
  h5_params, caf_params : dict
    A dictionary where the key is a variable alias, and the value is a 
    list of length 2 in which the 0th item is the slope and the 1st item 
    is the y-intercept for that variable's fit
    
  Epart : str
    The energy variable to correct. Must be one of 
    {'Ehad', 'EP', 'EN', 'EPip', 'EPim', 'EPi0', 'EOther'}
  '''
  # Get the old data
  h5_var_str = epart_to_h5_str_dict[Epart]
  old_Epart_reco_data = output_h5_file[subfile]['nd_paramreco'][h5_var_str][:]
  
  # Apply the correction
  new_Epart_reco_data = (old_Epart_reco_data - h5_params[Epart][1]) * caf_params[Epart][0] / h5_params[Epart][0] + caf_params[Epart][1]

  output_h5_file[subfile]['nd_paramreco'][h5_var_str][:] = new_Epart_reco_data  


def correct_Eother_reco(output_h5_file, subfile, h5_params, caf_params):
  ''' Correct down the reconstructed energy associated with other 
  particles (excluding protons, neutrons, pions) in a subfile of the 
  output HDF5 file so that it looks like a CAF file.
  
  Inputs
  ------
  output_h5_file : HDF5 file
    The file to correct
    
  subfile : str
    The subfile withing output_h5_file to correct
    
  h5_params, caf_params : dict
    A dictionary where the key is a variable alias, and the value is a 
    list of length 2 in which the 0th item is the slope and the 1st item 
    is the y-intercept for that variable's fit
  '''
  # Set the reconstructed energy attributed to "other" particles to be 
  # the difference between the reconstructed hadronic energy and the 
  # reconstructed energies for the different particle types.
  # Calculate the reco. hadronic energy first. This has already been 
  # corrected
  corr_Ev_reco_data = output_h5_file[subfile]['nd_paramreco']['Ev_reco'][:]
  corr_Elep_reco_data = output_h5_file[subfile]['nd_paramreco']['Elep_reco'][:]
  corr_Ehad_reco_data = np.array(corr_Ev_reco_data - corr_Elep_reco_data)
  
  # Calculate the reco. energy for the different particle types
  corr_Epart_reco_data = np.array(
    output_h5_file[subfile]['nd_paramreco']['eRecoP'][:])
  for part in ['eRecoN', 'eRecoPip', 'eRecoPim', 'eRecoPi0']:
    corr_Epart_reco_data += np.array(
      output_h5_file[subfile]['nd_paramreco'][part][:])
      
  # Perform the subtraction
  new_Eother_reco_data = corr_Ehad_reco_data - corr_Epart_reco_data
  
  output_h5_file[subfile]['nd_paramreco']['eRecoOther'][:] = new_Eother_reco_data


