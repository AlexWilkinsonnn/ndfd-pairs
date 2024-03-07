# fd\_detsim\_reco

Loads FD energy depositions into an art-root file and runs FD simulation to generate detector
simulation and reconstruction. ND detector simulation is also loaded so it can be projected onto FD
wires to make ND-FD detector response pairs.

## System

Recommended to do this on a dunegpvm

## Setup

```
# Get code on ndfd_pairs branch of duneextrapolation
git submodule update --init --remote fd_detsim_reco/duneextrapolation

# Build larsoft module
cd fd_detsim_reco/
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
mkdir larsoft_area
cd larsoft_area/
mrb newDev -v v09_78_04 -q e20:prof
source localProducts_larsoft_v09_78_04_e20_prof/setup
cd srcs/
cp -r ../../duneextrapolation .
mrb g --tag v09_78_03d01 dunereco
mrb uc
cd ../
cd $MRB_BUILDDIR
mrbsetenv && mrb i --generator=ninja && mrbslp

# Prep for submitting jobs
setup duneutil v09_78_03d01 -q e20:prof
setup jobsub_client v_lite
setup_fnal_security
```

## Input Data

Requires HDF5 files from previous `nd_detsim` stage. Should be in a directory accessible from
jobsub (`/pnfs/...`)

## Instructions

1. Create job tarball
  ```
  cd larsoft_area/
  tar -czvf jobdata.tar.gz srcs/ build_slf7.x86_64/ localProducts_larsoft_v09_78_04_e20_prof/
  ```

2. Edit the file `larsoft_area/srcs/duneextrapolation/scripts/jobs/produce_fd_pair_reco_resp.sh` to
   set output paths

3. Submit to grid
  ```
  jobsub_submit -G dune -N 100 --disk=60Gb --memory=6000MB --expected-lifetime=54h --cpu=1 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC,OFFSITE --tar_file_name=dropbox:///<path_to_repo>/fd_detsim_reco/larsoft_area/jobdata.tar.gz --use-cvmfs-dropbox -l '+SingularityImage=\"/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest\"' --append_condor_requirements='(TARGET.HAS_Singularity==true&&TARGET.HAS_CVMFS_dune_opensciencegrid_org==true&&TARGET.HAS_CVMFS_larsoft_opensciencegrid_org==true&&TARGET.CVMFS_dune_opensciencegrid_org_REVISION>=1105&&TARGET.HAS_CVMFS_fifeuser1_opensciencegrid_org==true&&TARGET.HAS_CVMFS_fifeuser2_opensciencegrid_org==true&&TARGET.HAS_CVMFS_fifeuser3_opensciencegrid_org==true&&TARGET.HAS_CVMFS_fifeuser4_opensciencegrid_org==true)' file:///<path_to_repo>/fd_detsim_reco/larsoft_area/srcs/duneextrapolation/scripts/jobs/produce_fd_pair_reco_resp.sh <path_to_nd_detsim_output>
  ```
  There are two `<path_to_repo>` a one `<path_to_nd_detsim_output>` that need to be substituted for. `-N` should be set to the number of files in the input directory (ND detsim HDF5 files). Disk usage and lifetime is for files with ~300 events, scale accordingly.

## Output

The job will add the following keys to the HDF5 files:

* `fd_reco` (FD reconstruction of full LArBath event)
  ```
  dtype([('eventID', '<i4'), ('numu_score', '<f4'), ('nue_score', '<f4'), ('nuc_score', '<f4'), ('nutau_score', '<f4'), ('antinu_score', '<f4'), ('0_p_score', '<f4'), ('1_p_score', '<f4'), ('2_p_score', '<f4'), ('N_p_score', '<f4'), ('0_pi_score', '<f4'), ('1_pi_score', '<f4'), ('2_pi_score', '<f4'), ('N_pi_score', '<f4'), ('0_pi0_score', '<f4'), ('1_pi0_score', '<f4'), ('2_pi0_score', '<f4'), ('N_pi0_score', '<f4'), ('0_n_score', '<f4'), ('1_n_score', '<f4'), ('2_n_score', '<f4'), ('N_n_score', '<f4'), ('numu_nu_E', '<f4'), ('numu_had_E', '<f4'), ('numu_lep_E', '<f4'), ('numu_reco_method', '<i4'), ('numu_longest_track_contained', '<i4'), ('numu_longest_track_mom_method', '<i4'), ('nue_nu_E', '<f4'), ('nue_had_E', '<f4'), ('nue_lep_E', '<f4'), ('nue_reco_method', '<i4'), ('nc_nu_E', '<f4'), ('nc_had_E', '<f4'), ('nc_lep_E', '<f4'), ('nc_reco_method', '<i4')])
  ```

* `fd_resp` (FD detsim of event with inside ND-LAr mask applied)
  * `eventID`
    * `tpc_set` ([0, 11])
      * `rop_id` ([0,3] 0 is U (induction) 1 is V (induction) 2&3 are Z (collection))
        Induction: `<HDF5 dataset "0": shape (800, 4492), type "<i2">`
        Collection: `<HDF5 dataset "0": shape (480, 4492), type "<i2">`
        The first dimension is channel number and the second is time tick.
        Note that the induction channels wrap around the APA

* `nd_packet_wire_projs` (ND 3d packet positions rotated and translated to align with FD energy
depositions and projected onto wire and tick)
  * `eventID`
    * `tpc_set` ([0, 11])
      * `rop_id` ([0,3] 0 is U (induction) 1 is V (induction) 2&3 are Z (collection))
        ```
        dtype([('adc' '<u4'), ('local_ch', '<i4'), ('tick', '<i4'), ('nd_drift_dist', '<f8'), ('fd_drift_dist', '<f8'), ('nd_x_module', '<f8'), ('wire_dist', '<f8'), ('forward_facing_anode', '<f4')])
        ```

* `ndfd_vertices_projs` (Wire and tick projection of shared vertex after aligning ND packets with FD energy depositions)
  ```
  dtype([('eventID', '<i4'), ('local_ch', '<i4'), ('tick', '<i4'), ('tpc_set', '<i4'), ('readout', '<i4')])
  ```

## Notes

* The drift coordinate reconstruction in the `nd_detsim` stage was a little short and I could not
figure out why. This lead to positions in forward facing anodes being smaller and backward facing
anode being larger. I correct this empirically with the `NDProjForwardAnodeXShift` and
`NDProjBackwardAnodeXShift`.
* The default time readout window is 6000 ticks but 4492 is enough to have the full FD volume in
the trigger

## TDR Pairs

Loads FD energy depositions into an art-root file with duneextrapolation. Runs ionisation and electron drift to make SimChannels. Writes SimChannels to original HDF5 file. SimChannels read into art-root file with dunetpc. Run TDR-era FD simulation to generate reconstruction. FD reconstruction written back to HDF5 file to finish ND-FD reconstruction pairs.

### Setup

Not bothering with git submodules for the duneextrapolation part of this build because it is just going to get messy

```
# Get main branch of dunetpc
git submodule update --init --remote fd_detsim_reco/dunetpc

# Build duneextrapolation larsoft module
cd fd_detsim_reco/
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
mkdir larsoft_area_duneextrapolation
cd larsoft_area_duneextrapolation/
mrb newDev -v v09_78_03 -q e20:prof
cd srcs/
git clone git@github.com:AlexWilkinsonnn/duneextrapolation.git
cd duneextrapolation/
git checkout simchannels_for_v07_06_02_dunetpc
cd ../
git clone git@github.com:AlexWilkinsonnn/larsim.git
cd larsim/
git checkout v09_38_01-elecdrift_fix
cd ../../
cp srcs/duneextrapolation/scripts/dev/* .
# go into setup.sh and update paths
source setup.sh
mrb g --tag v09_78_03d01 dunecore
mrb g --tag v09_78_03 larsoft
mrb g --tag v09_08_03 larexamples
mrb uc
source build_ninja.sh

# Open a fresh shell...

# Build dunetpc larsoft module
cd fd_detsim_reco/
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
mkdir larsoft_area_dunetpc
cd larsoft_area_dunetpc/
mrb newDev -v v07_06_02 -q e17:prof
cd srcs/
cp -r ../../dunetpc .
cd ../
cp srcs/dunetpc/scripts/dev/* .
# go into setup.sh and update paths
source setup.sh
mrb uc
source build_ninja.sh # this will fail and is a required step
source setup_old_cmake.sh
source build_ninja.sh # this should now build

# Open a fresh shell...

# Prep for submitting jobs
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup duneutil v09_78_03d01 -q e20:prof
setup jobsub_client v_lite
setup_fnal_security
```

### Instructions

1. Create job tarballs
  ```
  cd larsoft_area_duneextrapolation/
  tar -czvf jobdata_duneextrapolation.tar.gz srcs/ build_slf7.x86_64/ localProducts_larsoft_v09_78_03_e20_prof/
  cd ../larsoft_area_dunetpc/
  tar -czvf jobdata_dunetpc.tar.gz srcs/ build_slf7.x86_64/ localProducts_larsoft_v07_06_02_e17_prof/
  ```

2. Edit the files `larsoft_area_duneextrapolation/srcs/duneextrapolation/scripts/jobs/produce_fd_pair_simchannels.sh` and `larsoft_area_dunetpc/srcs/dunetpc/scripts/jobs/produce_fd_pair_reco_from_simchannels.sh` to set output paths

3. Submit to grid the duneextrapolation SimChannel jobs
  ```
  jobsub_submit -G dune -N 100 --disk=30Gb --memory=6000MB --expected-lifetime=2h --cpu=1 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC,OFFSITE --tar_file_name=dropbox://<path_to_repo>/fd_detsim_reco/larsoft_area_duneextrapolation/jobdata_duneextrapolation.tar.gz --use-cvmfs-dropbox -l '+SingularityImage=\"/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest\"' --append_condor_requirements='(TARGET.HAS_Singularity==true&&TARGET.HAS_CVMFS_dune_opensciencegrid_org==true&&TARGET.HAS_CVMFS_larsoft_opensciencegrid_org==true&&TARGET.CVMFS_dune_opensciencegrid_org_REVISION>=1105&&TARGET.HAS_CVMFS_fifeuser1_opensciencegrid_org==true&&TARGET.HAS_CVMFS_fifeuser2_opensciencegrid_org==true&&TARGET.HAS_CVMFS_fifeuser3_opensciencegrid_org==true&&TARGET.HAS_CVMFS_fifeuser4_opensciencegrid_org==true)' file://<path_to_repo>/fd_detsim_reco/larsoft_area_duneextrapolation/srcs/duneextrapolation/scripts/jobs/produce_fd_pair_simchannels.sh <path_to_ndfd_depos>
  ```
  In this step and the next there are two `<path_to_repo>` a one `<path_to_ndfd_depos>/<path_to_ndfd_fdsimchannels>` that need to be substituted for. `-N` should be set to the number of files in the input directory (HDF5 files). Disk usage and lifetime is for files with ~300 events, scale accordingly.

4. Submit to grid the dunetpc reco jobs
  ```
  jobsub_submit -G dune -N 100 --disk=30Gb --memory=6000MB --expected-lifetime=12h --cpu=1 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC,OFFSITE --tar_file_name=dropbox://<path_to_repo>/fd_detsim_reco/larsoft_area_dunetpc/jobdata_dunetpc.tar.gz --use-cvmfs-dropbox -l '+SingularityImage=\"/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest\"' --append_condor_requirements='(TARGET.HAS_Singularity==true&&TARGET.HAS_CVMFS_dune_opensciencegrid_org==true&&TARGET.HAS_CVMFS_larsoft_opensciencegrid_org==true&&TARGET.CVMFS_dune_opensciencegrid_org_REVISION>=1105&&TARGET.HAS_CVMFS_fifeuser1_opensciencegrid_org==true&&TARGET.HAS_CVMFS_fifeuser2_opensciencegrid_org==true&&TARGET.HAS_CVMFS_fifeuser3_opensciencegrid_org==true&&TARGET.HAS_CVMFS_fifeuser4_opensciencegrid_org==true)' file://<path_to_repo>/fd_detsim_reco/larsoft_area_dunetpc/srcs/dunetpc/scripts/jobs/produce_fd_pair_reco_from_simchannels.sh <path_to_ndfd_fdsimchannels>
  ```

