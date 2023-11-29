# ndfd-pairs

Code for ND-FD event pair making. The process is split into 3 stages:

1. Neutrino event generation, particle propagation, and selection (`ndfd_depos/`)
2. Near detector detector simulation and 3d position reconstruction (`nd_detsim/`)
3. Far detector detector simulation and reconstructions (`fd_detsim_reco/`)

To produce pairs, visit the directories above in order and follow instructions. Due to differences
in resource and software requirements this workflow is not expected to be done on one machine. This
repository will need to be cloned at a few different places and the relevant submodules for each
step downloaded.

The outputted HDF5 files contain all the required information to make paired data of (ND detector
response, FD reconstruction) and (ND detector response, FD detector response).

