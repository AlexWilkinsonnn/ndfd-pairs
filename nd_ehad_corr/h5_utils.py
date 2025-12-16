"""Functions used to extract information from .h5 files"""

# Author: Colin Weber (webe1077@umn.edu)
# Date: 28 August 2024


import os

import h5py
import pandas as pd
import math


def make_csv_from_paired(f):
    """Checks to see if the .h5 file has been turned into a .csv file 
    yet by searching for the associated .csv file in the same directory. 
    If not, this function extracts the data from the input .h5 file 
    and organizes it into a Pandas DataFrame, which is much easier to 
    work with. Either way, it returns the name of the .csv file. 
    
    Parameters
    __________
    f : h5py File object
        The .h5 file to extract information from.
        f should be the output of 
        ndfd-pairs/fd_detsim_reco/larsoft_area_dunetpc/add_h5_ndfdreco.py
        
    Returns
    _______
    csv_filename : str
        The full path and name of the dataframe saved as a .csv file. 
    """
    # Construct the name of the .csv file
    csv_filename = f.filename.replace("h5", "csv")
    
    # Check to see if the .csv file exists. If not create it.
    if not os.path.isfile(csv_filename):
        '''Loop over the groups in the h5 file, extract the desired 
        variables, and append to the data array (which will be converted 
        into the DataFrame).'''
        event_counter = 0
        data_array = []
        for file in f.keys():
            n_events = len(f[file]["fd_reco"]["eventID"])
            for event in range(n_events):
                print(".h5 event: %i" % event_counter)
                eventID = event_counter
                Ev_true = f[file][
                    "nd_paramreco"]["Ev"][event]
                ND_Ev_reco = f[file][
                    "nd_paramreco"]["Ev_reco"][event]
                FD_numu_nu_E = f[file][
                    "fd_reco"]["numu_nu_E"][event]
                Elep_true = f[file][
                    "nd_paramreco"]["LepE"][event]
                ND_Elep_reco = f[file][
                    "nd_paramreco"]["Elep_reco"][event]
                FD_numu_lep_E = f[file][
                    "fd_reco"]["numu_lep_E"][event]
                    
                Ehad_true = Ev_true - Elep_true
                
                ND_Ehad_reco = ND_Ev_reco - ND_Elep_reco
                
                FD_Ehad_reco = f[file][
                    "fd_reco"]["numu_had_E"][event]
                    
                ND_n_proton_true = f[file][
                    "nd_paramreco"]["nP"][event]
                ND_proton_reco_E = f[file][
                    "nd_paramreco"]["eRecoP"][event]
                ND_proton_true_E = f[file][
                    "nd_paramreco"]["eP"][event]
                ND_n_neutron_true = f[file][
                    "nd_paramreco"]["nN"][event]
                ND_neutron_reco_E = f[file][
                    "nd_paramreco"]["eRecoN"][event]
                ND_neutron_true_E = f[file][
                    "nd_paramreco"]["eN"][event]
                ND_n_pip_true = f[file][
                    "nd_paramreco"]["nipip"][event]
                ND_pip_reco_E = f[file][
                    "nd_paramreco"]["eRecoPip"][event]
                ND_pip_true_E = f[file][
                    "nd_paramreco"]["ePip"][event]
                ND_n_pim_true = f[file][
                    "nd_paramreco"]["nipim"][event]
                ND_pim_reco_E = f[file][
                    "nd_paramreco"]["eRecoPim"][event]
                ND_pim_true_E = f[file][
                    "nd_paramreco"]["ePim"][event]
                ND_n_pi0_true = f[file][
                    "nd_paramreco"]["nipi0"][event]
                ND_pi0_reco_E = f[file][
                    "nd_paramreco"]["eRecoPi0"][event]
                ND_pi0_true_E = f[file][
                    "nd_paramreco"]["ePi0"][event]
                ND_n_other_true = \
                    f[file]["nd_paramreco"]["nikp"][event] + \
                    f[file]["nd_paramreco"]["nikm"][event] + \
                    f[file]["nd_paramreco"]["nik0"][event] + \
                    f[file]["nd_paramreco"]["niem"][event] + \
                    f[file]["nd_paramreco"]["niother"][event] + \
                    f[file]["nd_paramreco"]["nNucleus"][event] + \
                    f[file]["nd_paramreco"]["nUNKNOWN"][event]
                ND_other_reco_E = f[file][
                    "nd_paramreco"]["eRecoOther"][event]
                ND_other_true_E = f[file][
                    "nd_paramreco"]["eOther"][event]
                ND_reco_theta = (f[file][
                    "nd_paramreco"]["theta_reco"][event]) * 180 / math.pi
                ND_true_theta = (f[file][
                    "nd_paramreco"]["LepNuAngle"][event]) * 180 / math.pi
                FD_numu_score = f[file][
                    "fd_reco"]["numu_score"][event]
                FD_nue_score = f[file][
                    "fd_reco"]["nue_score"][event]
                FD_nutau_score = f[file][
                    "fd_reco"]["nutau_score"][event]
                FD_nc_score = f[file][
                    "fd_reco"]["nc_score"][event]
                ND_reco_numu = f[file][
                    "nd_paramreco"]["reco_numu"][event]
                ND_reco_nue = f[file][
                    "nd_paramreco"]["reco_nue"][event]
                ND_reco_nc = f[file][
                    "nd_paramreco"]["reco_nc"][event]
                ND_reco_q = f[file][
                    "nd_paramreco"]["reco_q"][event]
                TrueNuPDG = f[file][
                    "nd_paramreco"]["nuPDG"][event]
                CC_flag = f[file][
                    "nd_paramreco"]["isCC"][event]
                ND_vtx_x = f[file][
                    "nd_paramreco"]["vtx_x"][event]
                ND_vtx_y = f[file][
                    "nd_paramreco"]["vtx_y"][event]
                ND_vtx_z = f[file][
                    "nd_paramreco"]["vtx_z"][event]
                FD_vtx_x = f[file][
                    "fd_vertices"]["x_vert"][event]
                FD_vtx_y = f[file][
                    "fd_vertices"]["y_vert"][event]
                FD_vtx_z = f[file][
                    "fd_vertices"]["z_vert"][event]
                ND_muon_contained = f[file][
                    "nd_paramreco"]["muon_contained"][event]
                ND_muon_tracker = f[file][
                    "nd_paramreco"]["muon_tracker"][event]
                ND_muon_ecal = f[file][
                    "nd_paramreco"]["muon_ecal"][event]
                ND_muon_exit = f[file][
                    "nd_paramreco"]["muon_exit"][event]
                ND_reco_lepton_pdg = f[file][
                    "nd_paramreco"]["reco_lepton_pdg"][event]
                ND_reco_q = f[file][
                    "nd_paramreco"]["reco_q"][event]
                ND_muon_end_x = f[file][
                    "nd_paramreco"]["muon_endpoint_x"][event]
                ND_muon_end_y = f[file][
                    "nd_paramreco"]["muon_endpoint_y"][event]
                ND_muon_end_z = f[file][
                    "nd_paramreco"]["muon_endpoint_z"][event]
                ND_Ehad_veto = f[file][
                    "nd_paramreco"]["Ehad_veto"][event]
                ND_lep_momx = f[file][
                    "nd_paramreco"]["LepMomX"][event]
                ND_lep_momy = f[file][
                    "nd_paramreco"]["LepMomY"][event]
                ND_lep_momz = f[file][
                    "nd_paramreco"]["LepMomZ"][event]
                ND_nu_momx = f[file][
                    "nd_paramreco"]["NuMomX"][event]
                ND_nu_momy = f[file][
                    "nd_paramreco"]["NuMomY"][event]
                ND_nu_momz = f[file][
                    "nd_paramreco"]["NuMomZ"][event]
                data_row = [eventID, Ev_true, ND_Ev_reco, FD_numu_nu_E,
                        Elep_true, ND_Elep_reco, FD_numu_lep_E, Ehad_true,
                        ND_Ehad_reco, FD_Ehad_reco, ND_n_proton_true,
                        ND_proton_reco_E, ND_proton_true_E,
                        ND_n_neutron_true, ND_neutron_reco_E,
                        ND_neutron_true_E, ND_n_pip_true, ND_pip_reco_E,
                        ND_pip_true_E, ND_n_pim_true, ND_pim_reco_E,
                        ND_pim_true_E, ND_n_pi0_true, ND_pi0_reco_E,
                        ND_pi0_true_E, ND_n_other_true, ND_other_reco_E,
                        ND_other_true_E, ND_reco_theta, ND_true_theta,
                        FD_numu_score, FD_nue_score, FD_nutau_score,
                        FD_nc_score, ND_reco_numu, ND_reco_nue, ND_reco_nc,
                        ND_reco_q, TrueNuPDG, CC_flag, ND_vtx_x, ND_vtx_y, 
                        ND_vtx_z, FD_vtx_x, FD_vtx_y, FD_vtx_z, 
                        ND_muon_contained, ND_muon_tracker, ND_muon_ecal, 
                        ND_muon_exit, ND_reco_lepton_pdg, ND_reco_q, 
                        ND_muon_end_x, ND_muon_end_y, ND_muon_end_z, 
                        ND_Ehad_veto, ND_lep_momx, ND_lep_momy, 
                        ND_lep_momz, ND_nu_momx, ND_nu_momy, ND_nu_momz]
            
                data_array.append(data_row)
                event_counter += 1
                
        # Create the pandas dataframe
        dataframe = pd.DataFrame(data=data_array, columns=[
            "eventID", "Ev_true", "ND_Ev_reco", "FD_numu_nu_E",
            "Elep_true", "ND_Elep_reco", "FD_numu_lep_E", "Ehad_true",
            "ND_Ehad_reco", "FD_Ehad_reco", "ND_n_proton_true",
            "ND_proton_reco_E", "ND_proton_true_E", "ND_n_neutron_true",
            "ND_neutron_reco_E", "ND_neutron_true_E", "ND_n_pip_true",
            "ND_pip_reco_E", "ND_pip_true_E", "ND_n_pim_true",
            "ND_pim_reco_E", "ND_pim_true_E", "ND_n_pi0_true",
            "ND_pi0_reco_E", "ND_pi0_true_E", "ND_n_other_true",
            "ND_other_reco_E", "ND_other_true_E", "ND_reco_theta",
            "ND_true_theta", "FD_numu_score", "FD_nue_score",
            "FD_nutau_score", "FD_nc_score", "ND_reco_numu",
            "ND_reco_nue", "ND_reco_nc", "ND_reco_q", "TrueNuPDG", 
            "CC_flag", "ND_vtx_x", "ND_vtx_y", "ND_vtx_z", "FD_vtx_x", 
            "FD_vtx_y", "FD_vtx_z", "ND_muon_contained", 
            "ND_muon_tracker", "ND_muon_ecal", "ND_muon_exit", 
            "ND_reco_lepton_pdg", "ND_reco_q", "ND_muon_end_x", 
            "ND_muon_end_y", "ND_muon_end_z", "ND_Ehad_veto",
            "ND_lep_momx", "ND_lep_momy", "ND_lep_momz", 
            "ND_nu_momx","ND_nu_momy", "ND_nu_momz"])
    
        # Save the dataframe
        dataframe.to_csv(csv_filename)
    
    return csv_filename
    print("converted pairs to csv")
