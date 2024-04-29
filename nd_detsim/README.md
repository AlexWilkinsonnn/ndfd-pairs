# nd_detsim

Runs ND-LAr detector simulation (larnd-sim) and uses packet timestamps and triggers to reconstruct drift coordinate (larpixsoft).

## System

Requires GPU resources to larnd-sim. Job submission scripts and instructions on how to use them are for clusters using Slurm. For other workload managers these scripts will need to be changed.

## Setup

```
# Get code on v0.3.4_ndfd_pairs branch of fork of larnd-sim
git submodule update --init --remote nd_detsim/larnd-sim
# Get code on ndfd_pairs branch of larpixsoft
git submodule update --init --remote nd_detsim/larpixsoft
cd nd_detsim
```

## Input Data

Requires HDF5 files from the previous `ndfd_depos` stage

## Instructions

1. Setup a python environment using the `requirements.txt` (venv or conda). Do this on the GPU node
   you will run on so that all the correct cuda stuff is there. I found it easier to use conda to
   get packages since it installs `cudatoolkit` for you. Just make sure you have the conda-forge
   channel setup: `conda config --add channels conda-forge`.

2. Make a file `larnd-sim/setup.sh` that activates this python environment and adds larnd-sim to
   the python path:
  ```
  export PYTHONPATH="/<base_path>/ndfd-pairs/ndfd_detsim/larnd-sim:$PYTHONPATH"
  ```

3. Edit `run_larndsim_ndfdpairs.sh` to configure the slurm job (`#SBATCH` options) and set the
   directory paths. The array Slurm parameter should be configured to like `1-<number_of_files_process>`.  `SCRATCH_DIR` should point to the path of the local disk for the GPU worker
   node.

4. Submit the job with:
  ```
  sbatch run_larndsim_ndfdpairs.sh /<path_to_ndfd_depos_data_directory> /<path_to_output_directory>
  ```

## Output

The job will add the following keys to the HDF5 files:

* `3d_packets` (ND-LAr detsim with reconstructed drift coordinate)
  ```
  dtype([('eventID', '<u4'), ('adc', '<f4'), ('x', '<f4'), ('x_module', '<f4'), ('y', '<f4'), ('z', '<f4'), ('z_module', '<f4'), ('forward_facing_anode', '<u4')])
  ```

* `mc_packets_assn` (contribution from to ND-LAr packets from true tracks)
  ```
  dtype([('track_ids', '<i8', (5,)), ('fraction', '<f8', (5,))])
  ```

* `packets` (raw ND-LAr detsim)
  ```
  dtype([('io_group', 'u1'), ('io_channel', 'u1'), ('chip_id', 'u1'), ('packet_type', 'u1'), ('downstream_marker', 'u1'), ('parity', 'u1'), ('valid_parity', 'u1'), ('channel_id', 'u1'), ('timestamp', '<u8'), ('dataword', 'u1'), ('trigger_type', 'u1'), ('local_fifo', 'u1'), ('shared_fifo', 'u1'), ('register_address', 'u1'), ('register_data', 'u1'), ('direction', 'u1'), ('local_fifo_events', 'u1'), ('shared_fifo_events', '<u2'), ('counter', '<u4'), ('fifo_diagnostics_enabled', 'u1'), ('first_packet', 'u1'), ('receipt_timestamp', '<u4')])
  ```

* `tracks` (ND event energy depositions, renamed from `segments` in previous stage)

* `_header`, `configs`, `messages` (unused, added by larnd-sim)

## Notes

* The `run_larndsim_ndfdpairs_specfile.sh` script can be used to run a job for a single input file

