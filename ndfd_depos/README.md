# ndfd_depos

Generates a neutrino event (genie with ND flux) and runs particle propagation (edep-sim) in a large liquid argon volume (LArBath). Event goes through translation and rotation (around beam axis) throws in ND to find a realisation that is selected. This event is moved to FD (Earth's curvature correction applied) and translation throws applied until a selected realisation is found. The result is the truth level information required to produce ND-FD pairs.

## System

Recommended to do this on a dunegpvm

## Setup

```
# Get code on ndfd_pairs branch of nd-sim-tools
git submodule update --remote nd-sim-tools
# Prep for submitting jobs
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup duneutil v09_78_06d00 -q e20:prof
setup jobsub_client v_lite
setup_fnal_security
```

## Instructions

1. Create job tarball
  ```
  cd nd-sim-tools/inputs
  tar -czvf jobdata.tar.gz ND_CAFMaker/ DUNE_ND_GeoEff/ sim_inputs_larbath_selected_ndfd_pairs/
  ```

2. Edit `nd-sim-tools/produce_scripts/produce_edep_larbath_transrots.sh` to set directories for
   outputted files, whether to save intermediary files, and flux options

3. Submit to grid
  ```
  jobsub_submit -G dune -N 200 --disk=60Gb --memory=3000MB --expected-lifetime=6h --cpu=1 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC,OFFSITE --tar_file_name=dropbox:///<path_to_repo>/ndfd_depos/nd-sim-tools/inputs/jobdata.tar.gz --use-cvmfs-dropbox -l '+SingularityImage=\"/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest\"' --append_condor_requirements='(TARGET.HAS_Singularity==true&&TARGET.HAS_CVMFS_dune_opensciencegrid_org==true&&TARGET.HAS_CVMFS_larsoft_opensciencegrid_org==true&&TARGET.CVMFS_dune_opensciencegrid_org_REVISION>=1105&&TARGET.HAS_CVMFS_fifeuser1_opensciencegrid_org==true&&TARGET.HAS_CVMFS_fifeuser2_opensciencegrid_org==true&&TARGET.HAS_CVMFS_fifeuser3_opensciencegrid_org==true&&TARGET.HAS_CVMFS_fifeuser4_opensciencegrid_org==true)' file:///dune/app/users/awilkins/nd-sim-tools/produce_scripts/produce_edep_larbath_transrots.sh 0 1e15
  ```
  `1e15` POT generates files with ~300 neutrinos, if you change this you will also need to scale
  `--disk` and `--expected_lifetime`. The first argument to the script (`0` in above) is just the
  offset for output file numbering, if you want to produce extra files to the same output directory after a previous job just increase this.

## Output

The outputted file is a HDF5 file with the following keys:

* `fd_deps` (selected FD event energy depositions)
  ```
  dtype([('eventID', '<u4'), ('uniqID', '<u4'), ('x_start', '<f4'), ('x_end', '<f4'), ('x', '<f4'), ('y_start', '<f4'), ('y_end', '<f4'), ('y', '<f4'), ('z_start', '<f4'), ('z_end', '<f4'), ('z', '<f4'), ('t0_start', '<f8'), ('t0_end', '<f8'), ('t0', '<f8'), ('dx', '<f4'), ('dEdx', '<f4'), ('dE', '<f4'), ('outsideNDLAr', '<i4')])
  ```

* `fd_vertices` (selected FD event vertex)
  ```
  dtype([('eventID', '<u4'), ('x_vert', '<f4'), ('y_vert', '<f4'), ('z_vert', '<f4')])
  ```

* `lepton` (selected ND event exiting muon truth information)
  ```
  dtype([('eventID', '<u4'), ('pdg', '<u4'), ('p0_MeV', '<f4'), ('p1_MeV', '<f4'), ('p2_MeV', '<f4'), ('p3_MeV', '<f4'), ('nn_lep_contained_prob', '<f4'), ('nn_lep_tracker_prob', '<f4'), ('lep_ke_outside_ndlar_MeV', '<f4')])
  ```

* `primaries` (4-vector of interaction primary particles)
  ```
  dtype([('eventID', '<u4'), ('pdg', '<u4'), ('p0_MeV', '<f4'), ('p1_MeV', '<f4'), ('p2_MeV', '<f4'), ('p3_MeV', '<f4')])
  ```

* `segments` (selected ND event energy depositions)
  ```
  dtype([('eventID', '<u4'), ('trackID', '<u4'), ('uniqID', '<u4'), ('pdgID', '<i4'), ('x_start', '<f4'), ('x_end', '<f4'), ('x', '<f4'), ('y_start', '<f4'), ('y_end', '<f4'), ('y', '<f4'), ('z_start', '<f4'), ('z_end', '<f4'), ('z', '<f4'), ('t_start', '<f4'), ('t_end', '<f4'), ('t', '<f4'), ('t0_start', '<f4'), ('t0_end', '<f4'), ('t0', '<f4'), ('n_electrons', '<u4'), ('n_photons', '<u4'), ('tran_diff', '<f4'), ('long_diff', '<f4'), ('dx', '<f4'), ('dEdx', '<f4'), ('dE', '<f4'), ('pixel_plane', '<i4')])
  ```

* `vertices` (selected ND event vertex)
  ```
  dtype([('eventID', '<u4'), ('x_vert', '<f4'), ('y_vert', '<f4'), ('z_vert', '<f4')])
  ```

## Notes

* `eventID` is unique within each HDF5 file
* `segments` and `vertices` are the dataset names required by `larnd-sim`
* `lepton` dataset is included to provide information required to make estimate the TMS reco
* The selection throws do not account for the fate of the muon in the ND but only the hadronic veto. The output of the geometric efficiency muon fate prediction neural network can be found in the `nn_lep_contained_prob` and `nn_lep_tracker_prob` columns of the `lepton` dataset. May want to apply some cut on this.
* There is also a `produce_edep-paramreco_larbath_transrots.sh` script which also runs the old ND parameterised reconstruction and writes this the HDF5. The reco is not very useful since it uses the ND event before the throws to find a selected realisation. It does provide true neutrino information summarised neatly on the HDF5 though (`Ev`, `lepE`, `isCC`, ...)

