"""
Author: Colin Weber (webe1077@umn.edu)
Date: 18 September 2025
Purpose: Using linear models of the training dataset's Ereco v. Etrue 
variables for hadrons at the ND, and the corresponding relationships for 
CAF files, correct the training dataset's Ereco distributions for each 
bin of Etrue. The correction is of the form 
E_corr = m_ND / m_LAr * (E - b_LAr) + b_ND, and each parameter's value 
was set using a guess-and-check method to find values that minimized a 
reduced chi squared metric between the training data and CAF file.

Command: python CorrNDEhadReco.py /exp/dune/data/users/colweber/larbath_ndfd_pairs/tdr_sample/1e6_FHC_numu_HitEhadCorr.h5 /exp/dune/data/users/colweber/larbath_ndfd_pairs/tdr_sample/1e6_FHC_numu_AdHocCorr_BruteForce.h5

Inputs
------
  input_h5 : the filename of the input paired training dataset. Must be 
  an HDF5 file.

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


epart_to_h5_str_dict = {
  'EP': 'eP',
  'EN': 'eN',
  'EPip': 'ePip',
  'EPim': 'ePim',
  'EPi0': 'ePi0',
  'EOther': 'eOther'
}
  
  
def main(args):
  '''Open the associated CSV for the input. If the CSV don't exist, 
  make it. CSVs are much faster to work with in python than the HDF5 
  files.'''
  print("Opening the file...")
  # Get the filenames of the (to-be-created) csv
  input_h5_csv_str = args.input_h5.replace("h5", "csv")
  
  # Check to see if they exist. If not, create them
  if not os.path.isfile(input_h5_csv_str):
    with h5py.File(args.input_h5) as input_h5_file:
      input_h5_csv_str = make_csv_from_paired(input_h5_file)
      
  # Open the CSV as Polars DataFrames
  input_h5_df = pl.scan_csv(input_h5_csv_str)
  
  # Make the output CSV string
  output_csv = args.output_h5.replace("h5", "csv")
  
  ''' Calculate correction and save to output files
  The sets of parameters below describe the [slope, int] parameters that 
  can calculate a reco value, given a true value. These parameters are set 
  by hand. The relationship that they describe can be rearranged to give a 
  true value, given a reco value. Doing this for both datasets and setting 
  them equal to each other gives a means of going from a reco LArbath 
  value to a reco CAF value, which is what we want to reproduce. This 
  correction is performed for all exclusive types of hadronic energy 
  except "Other", which instead is set so that the final Ehad reco 
  variable is the some of the Ehad particle variables.'''
  # print("Calculating correction and applying it to the data...")
  h5_params = {
    'Ehad': [0.67, 0.001],
    'EP': [0.86, 0.0001],
    'EN': [0.42, 0.001],
    'EPip': [0.79, 0.06],
    'EPim': [0.72, 0.082],
    'EPi0': [0.86, 0.13]
  }
  caf_params = {
    'Ehad': [0.665, 0.0],
    'EP': [0.8599, 0.0],
    'EN': [0.40, 0.0],
    'EPip': [0.78, 0.05],
    'EPim': [0.71, 0.08],
    'EPi0': [0.858, 0.13]
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
    "output_h5", type=str, 
    help="""The name of the output file with the ND Ehad reco variables 
    shifted. It will be saved in the same directory as input_h5."""
  )
  
  args = parser.parse_args()
  
  return args
  

if __name__ == "__main__":
  args = parse_arguments()
  main(args)
