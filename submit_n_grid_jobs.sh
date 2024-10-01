 # Author: Colin Weber (webe1077@umn.edu)
# Date: 3 June 2024
# Purpose: To send off n grid jobs in a for loop

# Command: bash submit_n_grid_jobs.sh \
#	[-i input script] \
#	[-t input tarbell] \
#	[-n number of times to run jobsub_submit.py]

# She-bang!
#!/bin/bash

# Read the inputs
input_script=$2
input_tar=$4
n_submits=$6

# Define constants
submissions_per_job=200

# Execute the for loop
for (( i=1; i<=$n_submits; i++ ));
do
	offset=$(($submissions_per_job * ($i - 1)))

	jobsub_submit -G dune -N $submissions_per_job --disk=10Gb --memory=3000MB --expected-lifetime=10h --cpu=1 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC,OFFSITE --tar_file_name=dropbox://$input_tar --use-cvmfs-dropbox -l '+SingularityImage=\"/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest\"' --append_condor_requirements='(TARGET.HAS_Singularity==true&&TARGET.HAS_CVMFS_dune_opensciencegrid_org==true&&TARGET.HAS_CVMFS_larsoft_opensciencegrid_org==true&&TARGET.CVMFS_dune_opensciencegrid_org_REVISION>=1105&&TARGET.HAS_CVMFS_fifeuser1_opensciencegrid_org==true&&TARGET.HAS_CVMFS_fifeuser2_opensciencegrid_org==true&&TARGET.HAS_CVMFS_fifeuser3_opensciencegrid_org==true&&TARGET.HAS_CVMFS_fifeuser4_opensciencegrid_org==true)' file://$input_script $offset 1e15
	
done
