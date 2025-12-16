// Author: Colin Weber (webe1077@umn.edu)
// Date: 5 September 2024
// Purpose: to define useful functions for working with CAF .root files, 
// such as turning ,root files into .csv files, which are much easier to 
// work with in python.
//
// Command: root -q "caf_utils.C(\"/exp/dune/data/users/colweber/TDR_CAFs/AlexBooth_ND_/PRISMCut_1e7_FD_FHC_TDR_CAF.root\", \"caf\")"

// Includes
// ROOT includes
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

// C++ includes
#include <fstream>
#include <sys/stat.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <sstream>

std::string make_csv_from_nd_caf(std::string full_filename, \
	std::string tree_name) {
/*	Checks to see if the .root file has been turned into a .csv file 
    yet by searching for the associated .csv file in the same directory. 
    If not, this function extracts the data from the input .root file 
    and organizes it into a .csv, which is much easier to 
    work with. Either way, it returns the name of the .csv file. 
    
    Parameters
    __________
    full_filename : std::string
        The full filename of the .root file to extract information from.
        f should be a DUNE ND CAF file from the directory 
        /exp/dune/data/users/colweber/TDR_CAFs/{AlexBooth_ND_, CiarinHasnip_ND_}
            
    tree_name : std::string
        The name of the tree in f to extract information from. Should be 
        one of {caf [for Alex Booth's CAFs], cafTree (for Ciarin Hasnip's 
        CAFs]}.
                
    Returns
    _______
    csv_filename : str
            The full path and name of the dataframe saved as a .csv file. 
*/	
	// Construct the name of the .csv file
	std::string csv_filename = full_filename;
	csv_filename.replace(csv_filename.find("root"), 4, "csv");
	
	// If the file doesn't exist, create it and fill it up
	struct stat buffer;
	bool exists = (stat(csv_filename.c_str(), &buffer) == 0);
	if(! exists) {
		// Open the .root file
		TFile* root_file = TFile::Open(full_filename.c_str());

		// Declare the tree reader and the values that it will read
		TTreeReader reader(tree_name.c_str(), root_file);
		TTreeReaderValue<double> Ev_true(reader, "Ev");
		TTreeReaderValue<double> ND_Ev_reco(reader, "Ev_reco");
		TTreeReaderValue<double> Elep_true(reader, "LepE");
		TTreeReaderValue<double> ND_Elep_reco(reader, "Elep_reco");
		TTreeReaderValue<int> n_proton_true(reader, "nP");
		TTreeReaderValue<double> ND_proton_reco_E(reader, "eRecoP");
		TTreeReaderValue<double> ND_proton_true_E(reader, "eP");
		TTreeReaderValue<int> n_neutron_true(reader, "nN");
		TTreeReaderValue<double> ND_neutron_reco_E(reader, "eRecoN");
		TTreeReaderValue<double> ND_neutron_true_E(reader, "eN");
		TTreeReaderValue<int> n_pip_true(reader, "nipip");
		TTreeReaderValue<double> ND_pip_reco_E(reader, "eRecoPip");
		TTreeReaderValue<double> ND_pip_true_E(reader, "ePip");
		TTreeReaderValue<int> n_pim_true(reader, "nipim");
		TTreeReaderValue<double> ND_pim_reco_E(reader, "eRecoPim");
		TTreeReaderValue<double> ND_pim_true_E(reader, "ePim");
		TTreeReaderValue<int> n_pi0_true(reader, "nipi0");
		TTreeReaderValue<double> ND_pi0_reco_E(reader, "eRecoPi0");
		TTreeReaderValue<double> ND_pi0_true_E(reader, "ePi0");
		TTreeReaderValue<int> n_kp_true(reader, "nikp");
		TTreeReaderValue<int> n_km_true(reader, "nikm");
		TTreeReaderValue<int> n_k0_true(reader, "nik0");
		TTreeReaderValue<int> n_em_true(reader, "niem");
		TTreeReaderValue<int> n_other_true(reader, "niother");
		TTreeReaderValue<int> n_nucleus_true(reader, "nNucleus");
		TTreeReaderValue<int> n_UNKNOWN_true(reader, "nUNKNOWN");
		TTreeReaderValue<double> ND_other_reco_E(reader, "eRecoOther");
		TTreeReaderValue<double> ND_other_true_E(reader, "eOther");
		TTreeReaderValue<double> ND_reco_theta_rad(reader, "theta_reco");
		TTreeReaderValue<double> ND_true_theta_rad(reader, "LepNuAngle");
		TTreeReaderValue<int> ND_numu_reco(reader, "reco_numu");
		TTreeReaderValue<int> ND_nue_reco(reader, "reco_nue");
		TTreeReaderValue<int> ND_nc_reco(reader, "reco_nc");
		TTreeReaderValue<int> ND_reco_q(reader, "reco_q");
		TTreeReaderValue<int> TrueNuPDG(reader, "nuPDG");
		TTreeReaderValue<int> CC_flag(reader, "isCC");
		TTreeReaderValue<double> ND_vtx_x(reader, "vtx_x");
		TTreeReaderValue<double> ND_vtx_y(reader, "vtx_y");
		TTreeReaderValue<double> ND_vtx_z(reader, "vtx_z");
		TTreeReaderValue<int> ND_muon_contained(reader, "muon_contained");
		TTreeReaderValue<int> ND_muon_tracker(reader, "muon_tracker");
		TTreeReaderValue<int> ND_muon_ecal(reader, "muon_ecal");
		TTreeReaderValue<int> ND_muon_exit(reader, "muon_exit");
		TTreeReaderValue<int> ND_reco_lepton_pdg(reader, "reco_lepton_pdg");
		TTreeReaderArray<float> ND_muon_end(reader, "muon_endpoint");
		TTreeReaderValue<double> ND_Ehad_veto(reader, "Ehad_veto");
		TTreeReaderValue<double> ND_lep_momx(reader, "LepMomX");
		TTreeReaderValue<double> ND_lep_momy(reader, "LepMomY");
		TTreeReaderValue<double> ND_lep_momz(reader, "LepMomZ");
		TTreeReaderValue<double> ND_nu_momx(reader, "NuMomX");
		TTreeReaderValue<double> ND_nu_momy(reader, "NuMomY");
		TTreeReaderValue<double> ND_nu_momz(reader, "NuMomZ");
		
		// Create the .csv file
		std::ofstream outfile(csv_filename);
		
		// Write the column headers
		outfile << ",eventID,Ev_true,ND_Ev_reco,Elep_true,ND_Elep_reco,Ehad_true,ND_Ehad_reco,ND_n_proton_true,ND_proton_reco_E,ND_proton_true_E,ND_n_neutron_true,ND_neutron_reco_E,ND_neutron_true_E,ND_n_pip_true,ND_pip_reco_E,ND_pip_true_E,ND_n_pim_true,ND_pim_reco_E,ND_pim_true_E,ND_n_pi0_true,ND_pi0_reco_E,ND_pi0_true_E,ND_n_other_true,ND_other_reco_E,ND_other_true_E,ND_reco_theta,ND_true_theta,ND_reco_numu,ND_reco_nue,ND_reco_nc,ND_reco_q,TrueNuPDG,CC_flag,ND_vtx_x,ND_vtx_y,ND_vtx_z,ND_muon_contained,ND_muon_tracker,ND_muon_ecal,ND_muon_exit,ND_reco_lepton_pdg,ND_muon_end_x,ND_muon_end_y,ND_muon_end_z,ND_Ehad_veto,ND_lep_momx,ND_lep_momy,ND_lep_momz,ND_nu_momx,ND_nu_momy,ND_nu_momz" << std::endl;

		// Loop over the tree
		int event_counter = 0;
		while(reader.Next()) {
			std::cout << event_counter << std::endl;
			// Perform some calculations
			double Ehad_true = *Ev_true - *Elep_true;
			double ND_Ehad_reco = *ND_Ev_reco - *ND_Elep_reco;
			int n_other_true = *n_kp_true + *n_km_true + *n_k0_true + \
				*n_em_true + n_other_true + *n_nucleus_true + \
				*n_UNKNOWN_true;
			double ND_reco_theta = *ND_reco_theta_rad * 180. / M_PI;
			double ND_true_theta = *ND_true_theta_rad * 180. / M_PI;
			outfile << event_counter << "," << event_counter << "," << *Ev_true << "," << *ND_Ev_reco << "," << *Elep_true << "," << *ND_Elep_reco << "," << Ehad_true << "," << ND_Ehad_reco << "," << *n_proton_true << "," << *ND_proton_reco_E << "," << *ND_proton_true_E << "," << *n_neutron_true << "," << *ND_neutron_reco_E << "," << *ND_neutron_true_E << "," << *n_pip_true << "," << *ND_pip_reco_E << "," << *ND_pip_true_E << "," << *n_pim_true << "," << *ND_pim_reco_E << "," << *ND_pim_true_E << "," << *n_pi0_true << "," << *ND_pi0_reco_E << "," << *ND_pi0_true_E << "," << n_other_true << "," << *ND_other_reco_E << "," << *ND_other_true_E << "," << ND_reco_theta << "," << ND_true_theta << "," << *ND_numu_reco << "," << *ND_nue_reco << "," << *ND_nc_reco << "," << *ND_reco_q << "," << *TrueNuPDG << "," << *CC_flag << "," << *ND_vtx_x << "," << *ND_vtx_y << "," << *ND_vtx_z << "," << *ND_muon_contained << "," << *ND_muon_tracker << "," << *ND_muon_ecal << "," << *ND_muon_exit << "," << *ND_reco_lepton_pdg << "," << ND_muon_end[0] << "," << ND_muon_end[1] << "," << ND_muon_end[2] << "," << *ND_Ehad_veto << "," << *ND_lep_momx << "," << *ND_lep_momy << "," << *ND_lep_momz << "," << *ND_nu_momx << "," << *ND_nu_momy << "," << *ND_nu_momz << std::endl;
			event_counter++;
		}
		root_file->Close();
		outfile.close();
	}
	return csv_filename;
}

void caf_utils(std::string util, std::string full_filename, \
	std::string tree_name) {
	if(util == "make_csv_from_nd_caf") {
		std::string nd_str = make_csv_from_nd_caf(full_filename, \
			tree_name);
	}
	else {
		std::cout << "Unknown value for util" << util << std::endl;
	}
}

