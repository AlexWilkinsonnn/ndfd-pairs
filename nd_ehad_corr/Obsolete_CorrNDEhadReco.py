"""
Author: Colin Weber (webe1077@umn.edu)
Date: 18 September 2025
Purpose: This correction is obsolete, but implements the version where 
each bin of true has a Gaussian fit to it, and then the centers of these 
are used to fit a line to the reco v. true data. The lines are then used \
to calculate the correction. By default, bin edges are determined so that 
each bin contains roughly the same number of events.

Command: python CorrNDEhadReco.py /exp/dune/data/users/colweber/larbath_ndfd_pairs/tdr_sample/1e6_FHC_numu_HitEhadCorr.h5 /exp/dune/data/users/colweber/TDR_CAFs/AlexBooth_ND_/PRISMCut_1e7_ND_FHC_TDR_CAF.root /exp/dune/data/users/colweber/larbath_ndfd_pairs/tdr_sample/1e6_FHC_numu_AdHocCorr_BruteForce.h5 --debug

Inputs
------
  input_h5 : the filename of the input paired training dataset. Must be 
  an HDF5 file.

  input_CAF: the filename of the input CAF reference file. Must be a ROOT 
  file.

  output_h5 : the filename fo the corrected output file. Must be an HDF5 
  file.

Outputs
-------
  output_h5 : an HDF5 file with the ND Ehad reco quantities shifted to 
  statistically look like the reference CAF file.	
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
from analysis_utils import *

from math import sqrt, floor
from lmfit import Parameters, Minimizer, report_fit
from sklearn.linear_model import TheilSenRegressor, RANSACRegressor
from scipy.stats import truncnorm, norm, gaussian_kde
from scipy.optimize import minimize_scalar


''' Global variables. The exponents are the powers that the number of 
events is taken to in order to determine the number of bins. Lower values 
means coarser binning. The dictionary is just to help with bookkeeping.'''
had_exp = 1/5
p_exp = 1/7
n_exp = 1/7
pip_exp = 1/7
pim_exp = 1/7
pi0_exp = 1/6

epart_to_h5_str_dict = {
  'EP': 'eP',
  'EN': 'eN',
  'EPip': 'ePip',
  'EPim': 'ePim',
  'EPi0': 'ePi0',
  'EOther': 'eOther'
}
  
  
def main(args):
  '''Open the associated CSV's for each file. If the CSV's don't exist, 
  make them. CSV's are much faster to work with in python than the HDF5 
  files. Similarly, turn the CAF into a CSV for uniformity.'''
  print("Opening the files...")
  # Get the filenames of the (to-be-created) csv's
  input_h5_csv_str = args.input_h5.replace("h5", "csv")
  input_caf_csv_str = args.input_caf.replace("root", "csv")
  
  # Check to see if they exist. If not, create them
  if not os.path.isfile(input_h5_csv_str):
    with h5py.File(args.input_h5) as input_h5_file:
      input_h5_csv_str = make_csv_from_paired(input_h5_file)
    
  if not os.path.isfile(input_caf_csv_str):
    raise NameError(input_caf_csv_str + " has not been created yet."\
      "\n Please run the function\n"\
      "root -q \"caf_utils.C(\\\"make_csv_from_nd_caf\\\", \\\"" + \
      args.input_caf + "\\\", \\\"caf\\\")\"\n"\
      "before performing this analysis.")
      
  # Open the CSVs as Polars DataFrames
  input_h5_df = pl.scan_csv(input_h5_csv_str)
  # input_caf_df = pl.scan_csv(input_caf_csv_str)
  
  # Make the output CSV string
  output_csv = args.output_h5.replace("h5", "csv")
  
  '''The quantities we will shift are EhadReco, eRecoP, eRecoN, eRecoPip, 
  eRecoPim, and eRecoPi0. Each of these needs its own binning scheme. We 
  will use a procedure where the number of bins along the true axis will 
  be the 1/6 root of the number of events (rounded down), and then the 
  bin edges will be set so that each bin of true energy has the same 
  number of events. We will then repeat this procedure to get the reco 
  axis binning. In this next chunk, we'll get the binning scheme for all 
  the true quantities'''
  print("Copying dataframes for each variable to smear...")
  # Make copies of the dataframes with the particle filters applied.
  input_h5_df_P = cp.deepcopy(
    input_h5_df.filter(pl.col('ND_n_proton_true') > 0))
  input_h5_df_N = cp.deepcopy(
    input_h5_df.filter(pl.col('ND_n_neutron_true') > 0))
  input_h5_df_Pip = cp.deepcopy(
    input_h5_df.filter(pl.col('ND_n_pip_true') > 0))
  input_h5_df_Pim = cp.deepcopy(
    input_h5_df.filter(pl.col('ND_n_pim_true') > 0))
  input_h5_df_Pi0 = cp.deepcopy(
    input_h5_df.filter(pl.col('ND_n_pi0_true') > 0))
    
  input_caf_df_P = cp.deepcopy(
    input_caf_df.filter(pl.col('ND_n_proton_true') > 0))
  input_caf_df_N = cp.deepcopy(
    input_caf_df.filter(pl.col('ND_n_neutron_true') > 0))
  input_caf_df_Pip = cp.deepcopy(
    input_caf_df.filter(pl.col('ND_n_pip_true') > 0))
  input_caf_df_Pim = cp.deepcopy(
    input_caf_df.filter(pl.col('ND_n_pim_true') > 0))
  input_caf_df_Pi0 = cp.deepcopy(
    input_caf_df.filter(pl.col('ND_n_pi0_true') > 0))
  
  # Count the number of events in each of the 7 categories for both input 
  # files. Assume PRISM cuts have already been applied
  n_h5 = input_h5_df.collect().select(pl.len()).item()
  n_P_h5 = input_h5_df_P.collect().select(pl.len()).item()
  n_N_h5 = input_h5_df_N.collect().select(pl.len()).item()
  n_Pip_h5 = input_h5_df_Pip.collect().select(pl.len()).item()
  n_Pim_h5 = input_h5_df_Pim.collect().select(pl.len()).item()
  n_Pi0_h5 = input_h5_df_Pi0.collect().select(pl.len()).item()
  
  n_caf = input_caf_df.collect().select(pl.len()).item()
  n_P_caf = input_caf_df_P.collect().select(pl.len()).item()
  n_N_caf = input_caf_df_N.collect().select(pl.len()).item()
  n_Pip_caf = input_caf_df_Pip.collect().select(pl.len()).item()
  n_Pim_caf = input_caf_df_Pim.collect().select(pl.len()).item()
  n_Pi0_caf = input_caf_df_Pi0.collect().select(pl.len()).item()
  
  # Get the binning schemes for the true column (bin edges)
  # Use the file with the lowest number of events
  print("Getting the bin edges of the true axis...")
  EhadTrueEdges = get_binning_scheme(input_h5_df, "Ehad_true", had_exp) \
    if n_h5 < n_caf else \
    get_binning_scheme(input_caf_df, "Ehad_true", had_exp)
  ePEdges = get_binning_scheme(input_h5_df_P, "ND_proton_true_E", 
    p_exp) if n_P_h5 < n_P_caf else \
    get_binning_scheme(input_caf_df_P, "ND_proton_true_E", p_exp)
  eNEdges = get_binning_scheme(input_h5_df_N, "ND_neutron_true_E", 
    n_exp) if n_N_h5 < n_N_caf else \
    get_binning_scheme(input_caf_df_N, "ND_neutron_true_E", n_exp)
  ePipEdges = get_binning_scheme(input_h5_df_Pip, "ND_pip_true_E", 
    pip_exp) if n_Pip_h5 < n_Pip_caf else \
    get_binning_scheme(input_caf_df_Pip, "ND_pip_true_E", pip_exp)
  ePimEdges = get_binning_scheme(input_h5_df_Pim, "ND_pim_true_E", 
    pim_exp) if n_Pim_h5 < n_Pim_caf else \
    get_binning_scheme(input_caf_df_Pim, "ND_pim_true_E", pim_exp)
  ePi0Edges = get_binning_scheme(input_h5_df_Pi0, "ND_pi0_true_E", 
    pi0_exp) if n_Pi0_h5 < n_Pi0_caf else \
    get_binning_scheme(input_caf_df_Pi0, "ND_pi0_true_E", pi0_exp)
    
  # Get the binning schemes for the reco columns
  # There will be one scheme for each gin of truth
  # Again, use the file with the fewest events
  print("Getting the bin edges of the reco axis...")
  # The following line can be used to just duplicate the binning scheme
  # EhadRecoEdges = [EhadTrueEdges for i in range(len(EhadTrueEdges) - 1)]
  EhadRecoEdges = fill_reco_edges(input_h5_df, input_caf_df, n_h5, n_caf,
    "ND_Ehad_reco", had_exp, "Ehad_true", EhadTrueEdges)
  eRecoPEdges = fill_reco_edges(input_h5_df_P, input_caf_df_P, n_h5, 
    n_caf, "ND_proton_reco_E", p_exp, "ND_proton_true_E", ePEdges)
  eRecoNEdges = fill_reco_edges(input_h5_df_N, input_caf_df_N, n_h5, 
    n_caf, "ND_neutron_reco_E", n_exp, "ND_neutron_true_E", eNEdges)
  eRecoPipEdges = fill_reco_edges(input_h5_df_Pip, input_caf_df_Pip, 
    n_h5, n_caf, "ND_pip_reco_E", pip_exp, "ND_pip_true_E", ePipEdges)
  eRecoPimEdges = fill_reco_edges(input_h5_df_Pim, input_caf_df_Pim, 
    n_h5, n_caf, "ND_pim_reco_E", pim_exp, "ND_pim_true_E", ePimEdges)
  eRecoPi0Edges = fill_reco_edges(input_h5_df_Pi0, input_caf_df_Pi0, 
    n_h5, n_caf, "ND_pi0_reco_E", pi0_exp, "ND_pi0_true_E", ePi0Edges)
  
  '''Now we will get the (Etrue, Ereco) ordered pairs for each dataset. 
  Using the reco binning scheme for each bin of true, we create a 
  histogram in reco space. After dividing each bin by its width, we fit a 
  Gaussian to the histogram. The mean is taken to be the Ereco value for 
  that bin of Etrue.'''
  print("Calculating the (true, reco_cv, true_unc, reco_width) values...")
  h5_Ehad_pairs, caf_Ehad_pairs = get_pairs(input_h5_df, input_caf_df, 
    "Ehad_true", "ND_Ehad_reco", EhadTrueEdges, EhadRecoEdges, "Ehad", 
    debug=args.debug)
  h5_EP_pairs, caf_EP_pairs = get_pairs(input_h5_df_P, input_caf_df_P, 
    "ND_proton_true_E", "ND_proton_reco_E", ePEdges, eRecoPEdges, 
    "Eproton", debug=args.debug)
  h5_EN_pairs, caf_EN_pairs = get_pairs(input_h5_df_N, input_caf_df_N, 
    "ND_neutron_true_E", "ND_neutron_reco_E", eNEdges, eRecoNEdges, 
    "Eneutron", debug=args.debug)
  h5_EPip_pairs, caf_EPip_pairs = get_pairs(input_h5_df_Pip, 
    input_caf_df_Pip, "ND_pip_true_E", "ND_pip_reco_E", ePipEdges, 
    eRecoPipEdges, "Epip", debug=args.debug)
  h5_EPim_pairs, caf_EPim_pairs = get_pairs(input_h5_df_Pim, 
    input_caf_df_Pim, "ND_pim_true_E", "ND_pim_reco_E", ePimEdges, 
    eRecoPimEdges, "Epim", debug=args.debug)
  h5_EPi0_pairs, caf_EPi0_pairs = get_pairs(input_h5_df_Pi0, 
    input_caf_df_Pi0, "ND_pi0_true_E", "ND_pi0_reco_E", ePi0Edges, 
    eRecoPi0Edges, "Epi0", debug=args.debug)

  # Fit lines.
  print("Fitting lines to the pairs...")
  h5_Ehad_slope, h5_Ehad_int = fit_line(h5_Ehad_pairs[:, 0], 
    h5_Ehad_pairs[:, 1], h5_Ehad_pairs[:, 2], h5_Ehad_pairs[:, 3], 
    "LArbath", "Ehad", debug=args.debug, vary_y_int=True)
  caf_Ehad_slope, caf_Ehad_int = fit_line(caf_Ehad_pairs[:, 0], 
    caf_Ehad_pairs[:, 1], caf_Ehad_pairs[:, 2], caf_Ehad_pairs[:, 3], 
    "NDHall", "Ehad", debug=args.debug, vary_y_int=True)
  
  h5_EP_slope, h5_EP_int = fit_line(h5_EP_pairs[:, 0], 
    h5_EP_pairs[:, 1], h5_EP_pairs[:, 2], h5_EP_pairs[:, 3], 
    "LArbath", "Eproton", False, debug=args.debug)
  caf_EP_slope, caf_EP_int = fit_line(caf_EP_pairs[:, 0], 
    caf_EP_pairs[:, 1], caf_EP_pairs[:, 2], caf_EP_pairs[:, 3], 
    "NDHall", "Eproton", debug=args.debug)
  
  h5_EN_slope, h5_EN_int = fit_line(h5_EN_pairs[:, 0], 
    h5_EN_pairs[:, 1], h5_EN_pairs[:, 2], h5_EN_pairs[:, 3], 
    "LArbath", "Eneutron", debug=args.debug)
  caf_EN_slope, caf_EN_int = fit_line(caf_EN_pairs[:, 0], 
    caf_EN_pairs[:, 1], caf_EN_pairs[:, 2], caf_EN_pairs[:, 3], 
    "NDHall", "Eneutron", debug=args.debug)
  
  h5_EPip_slope, h5_EPip_int = fit_line(h5_EPip_pairs[:, 0], 
    h5_EPip_pairs[:, 1], h5_EPip_pairs[:, 2], h5_EPip_pairs[:, 3], 
    "LArbath", "Epip", vary_y_int=True, debug=args.debug)
  caf_EPip_slope, caf_EPip_int = fit_line(caf_EPip_pairs[:, 0], 
    caf_EPip_pairs[:, 1], caf_EPip_pairs[:, 2], caf_EPip_pairs[:, 3], 
    "NDHall", "Epip", vary_y_int=True, debug=args.debug)
  
  h5_EPim_slope, h5_EPim_int = fit_line(h5_EPim_pairs[:, 0], 
    h5_EPim_pairs[:, 1], h5_EPim_pairs[:, 2], h5_EPim_pairs[:, 3], 
    "LArbath", "Epim", vary_y_int=True, debug=args.debug)
  caf_EPim_slope, caf_EPim_int = fit_line(caf_EPim_pairs[:, 0], 
    caf_EPim_pairs[:, 1], caf_EPim_pairs[:, 2], caf_EPim_pairs[:, 3], 
    "NDHall", "Epim", vary_y_int=True, debug=args.debug)
  
  h5_EPi0_slope, h5_EPi0_int = fit_line(h5_EPi0_pairs[:, 0], 
    h5_EPi0_pairs[:, 1], h5_EPi0_pairs[:, 2], h5_EPi0_pairs[:, 3], 
    "LArbath", "Epi0", vary_y_int=True, debug=args.debug)
  caf_EPi0_slope, caf_EPi0_int = fit_line(caf_EPi0_pairs[:, 0], 
    caf_EPi0_pairs[:, 1], caf_EPi0_pairs[:, 2], caf_EPi0_pairs[:, 3], 
    "NDHall", "Epi0", vary_y_int=True, debug=args.debug)
    
  # Use the least squares method to fit a line to the data, and debug
  '''h5_slopes = []
  h5_ints = []
  caf_slopes = []
  caf_ints = []
  for i in range(100):
    print(i)
    h5_Ehad_slope, h5_Ehad_int = perform_lsq_fit(input_h5_df, 
      "Ehad_true", "ND_Ehad_reco")
    caf_Ehad_slope, caf_Ehad_int = perform_lsq_fit(input_caf_df, 
      "Ehad_true", "ND_Ehad_reco")
    h5_slopes.append(h5_Ehad_slope)
    h5_ints.append(h5_Ehad_int)
    caf_slopes.append(caf_Ehad_slope)
    caf_ints.append(caf_Ehad_int)
  plt.hist(h5_slopes)
  print(np.mean(h5_slopes), np.std(h5_slopes))
  plt.show()
  plt.hist(h5_ints)
  print(np.mean(h5_ints), np.std(h5_ints))
  plt.show()
  plt.hist(caf_slopes)
  print(np.mean(caf_slopes), np.std(caf_slopes))
  plt.show()
  plt.hist(caf_ints)
  print(np.mean(caf_ints), np.std(caf_ints))
  plt.show()
    
  plot_scatter(input_h5_df, 'Ehad_true', 'ND_Ehad_reco', 10000, 
    'LArbath', 'Ehad', line=True, slope=h5_Ehad_slope, yint=h5_Ehad_int, 
    min_val=0.0, max_val=1.5)
  plot_scatter(input_caf_df, 'Ehad_true', 'ND_Ehad_reco', 10000, 
    'NDHall', 'Ehad', line=True, slope=caf_Ehad_slope, yint=caf_Ehad_int, 
    min_val=0.0, max_val=1.5)'''
  
  if(args.debug):
    compare_lin_fit(h5_Ehad_slope, h5_Ehad_int, caf_Ehad_slope, 
      caf_Ehad_int, "Ehad", plot_data=True, h5_pairs=h5_Ehad_pairs, 
      caf_pairs=caf_Ehad_pairs)
    compare_lin_fit(h5_EP_pairs, h5_EP_slope, h5_EP_int, 
      caf_EP_pairs, caf_EP_slope, caf_EP_int, "Eproton")
    compare_lin_fit(h5_EN_pairs, h5_EN_slope, h5_EN_int, 
      caf_EN_pairs, caf_EN_slope, caf_EN_int, "Eneutron")
    compare_lin_fit(h5_EPip_pairs, h5_EPip_slope, h5_EPip_int, 
      caf_EPip_pairs, caf_EPip_slope, caf_EPip_int, "Epip")
    compare_lin_fit(h5_EPim_pairs, h5_EPim_slope, h5_EPim_int, 
      caf_EPim_pairs, caf_EPim_slope, caf_EPim_int, "Epim")
    compare_lin_fit(h5_EPi0_pairs, h5_EPi0_slope, h5_EPi0_int, 
      caf_EPi0_pairs, caf_EPi0_slope, caf_EPi0_int, "Epi0")
  
  ''' Calculate correction and save to output files
  Using the two equations for the lines above, we solve each for the 
  true value, and then set those equal to each other. Using that 
  equation, we solve for the reco CAF value in terms of the reco HDF5 
  value. This tells us how to adjust the HDF5 value to look like a CAF 
  file. We perform this for Ehad reco and all of the particle type energy 
  variables except "Other", which is set so that Ehad is the sum of all 
  the particle type energies.'''
  print("Calculating correction and applying it to the data...")
  h5_params = {
    'Ehad': [h5_Ehad_slope, h5_Ehad_int],
    'EP': [h5_EP_slope, h5_EP_int],
    'EN': [h5_EN_slope, h5_EN_int],
    'EPip': [h5_EPip_slope, h5_EPip_int],
    'EPim': [h5_EPim_slope, h5_EPim_int],
    'EPi0': [h5_EPi0_slope, h5_EPi0_int]
  }
  caf_params = {
    'Ehad': [caf_Ehad_slope, caf_Ehad_int],
    'EP': [caf_EP_slope, caf_EP_int],
    'EN': [caf_EN_slope, caf_EN_int],
    'EPip': [caf_EPip_slope, caf_EPip_int],
    'EPim': [caf_EPim_slope, caf_EPim_int],
    'EPi0': [caf_EPi0_slope, caf_EPi0_int]
  }
  
  create_csv(input_h5_df, h5_params, caf_params, output_csv)
  
  shutil.copyfile(args.input_h5, args.output_h5)
  with h5py.File(args.output_h5, 'a') as output_h5_file:
    n_subfiles = len(output_h5_file.keys())
    print_modulo = math.ceil(n_subfiles / 100)
    for i, subfile in enumerate(output_h5_file.keys()):
      if(i % print_modulo == 0):
        print("Correcting subfile " + str(i) + " / " + str(n_subfiles))
      correct_Ev_reco(output_h5_file, subfile, h5_params, caf_params)
      correct_Epart_reco(output_h5_file, subfile, h5_params, caf_params, 
        'EP')
      correct_Epart_reco(output_h5_file, subfile, h5_params, caf_params, 
        'EN')
      correct_Epart_reco(output_h5_file, subfile, h5_params, caf_params, 
        'EPip')
      correct_Epart_reco(output_h5_file, subfile, h5_params, caf_params, 
        'EPim')
      correct_Epart_reco(output_h5_file, subfile, h5_params, caf_params, 
        'EPi0')
      correct_Eother_reco(output_h5_file, subfile, h5_params, 
        caf_params)

  
def parse_arguments():
  parser = argparse.ArgumentParser()
  
  parser.add_argument(
    "input_h5", type=str, 
    help="""Should be a single HDF5 file that has multiple HDF5 files 
    merged together inside. The inner HDF5 files should all be of the 
    same format."""
  )
  
  parser.add_argument(
    "input_caf", type=str, help="""Should be an ND CAF file."""
  )
  
  parser.add_argument(
    "output_h5", type=str, 
    help="""The name of the output file with the ND Ehad reco variables 
    shifted. It will be saved in the same directory as input_h5."""
  )
  
  parser.add_argument(
    "--debug", action='store_true',
    help="""If this option is used, print debug plots."""
  )
  
  args = parser.parse_args()
  
  return args
  

if __name__ == "__main__":
  args = parse_arguments()
  main(args)
