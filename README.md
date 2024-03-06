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

## TDR-era simulation ndfd pairs

There is code for making ND-FD reco-reco pairs with as close to the TDR software as I can get.
This is for training models that will be applicable to analyses that use MCC11.

TDR era means tag `last_unstructured_cafs` of `DUNE/ND_CAFMaker.git` for ND, tag `v07_06_02` of
`dunetpc` for FD simulation and reconstruction, and tag `v07_09_00` of `dunetpc` for FD
CAFMaker. These versions come from asking, snooping, and guessing.

Caveats of this dataset due to software restrictions:
- Cannot load SimEnergyDeposit into old `dunetpc`. Load them into newer larsoft with
  `duneextrapolation` and run the drifting module to get SimChannel. Load these into `dunetpc`.
- Diffucult to load selection throws output into ND CAFMaker. Don't do selection throws at the ND,
  FD event still has translations applied.

This workflow does not require the `nd_detsim/` step. Instructions can be found at the bottom of
the READMEs.

## Moving data around

To make use of both Fermigrid and GPU resources I had to copy files to and from dCache to an
external cluster. To make these copies safely we need to use `ifdh cp`. Doing this nto from a
dunegpvm can be tricky:

### Setup

```
git clone git@github.com:fermitools/cigetcert.git
cd cigetcert
git checkout tags/1.21
cd ../
export SINGULARITY_BIND="/cvmfs,/home,/tmp,..."
singularity shell /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7\:latest
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup python v3_9_2
python -m venv .venv.3_9_2_cigetcert
source .venv.3_9_2_cigetcert/bin/activate
pip install swig lxml M2Crypto pykerberos pyOpenSSL
```

### Copying

Put your krb5_fnal.conf on your home directory

```
export SINGULARITY_BIND="/cvmfs,/home,/tmp,..."
singularity shell /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7\:latest
export KRB5_CONFIG=/home/awilkins/krb5_fnal.conf
kinit
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
source .venv.3_9_2_cigetcert/bin/activate
python cigetcert/cigetcert -i 'Fermi National Accelerator Laboratory' -n
voms-proxy-init -noregen -rfc -voms dune:/dune/Role=Analysis
setup ifdhc
ifdh cp /pnfs/dune/persistent/...
```

## Using Output Data

The final output comes from the `fd_detsim_reco` stage and can be used to train macgine learning algorithms that map from ND detector response to FD reconstrution or to FD detector response. To prepare a dataset, iterate over the hd5f files and write the desired data to the desired format for input to a model. Note that when iterating the hdf5 files it is possible (though rare) that `fd_resp` or `fd_reco` will be missing for an event, in which case just skip this event.

