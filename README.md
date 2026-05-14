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

## Recent developments: 'DiffEdep' and ND sim. improvements

The main branch of this repository has not changed meaningfully in the recent past, as we've 
considered different improvements and directions the ndfd-transformer program will go. These
improvements and developments have taken place in the two other branches, first in `lep_swapper`,
and then in `DiffEdep`.

### `lep_swapper` branch

In Spring 2024, we thought we were ready to develop a training dataset with oscillations, where 
a muon neutrino event at the ND would be paired with an oscillated electron neutrino event at 
the FD. This feature was going to be developed in the `lep_swapper` branch. However, we soon 
discovered that we were not ready for this step because we saw that we were not simulating 
the ND properly. Because we had been simulating the entire ND hall as an infinite vat of LAr, 
the simulated muonic and hadronic components were both wrong, leading to the wrong 
relationship between true and reco. energies at the ND. For muons, the error was due to the 
gap between TMS/ND-GAr and ND-LAr, which when simulated as LAr instead of air leads to way 
more energy deposited by the muon in that region, an overall shorter muon track, and not 
enough energy reconstructed for tracker-matched muons. The fix here was to simulate just the muon 
in a realistic ND hall, and use those reco. quantities stitched together with the hadronic 
system as a complete ND event. The code to do this is in 
`ndfd_depos/nd-sim-tools/inputs/ND_CAFMaker/makeCAF_resim-muon.cxx`, which is called by an 
edited version of 
`/ndfd_depos/nd-sim-tools/produce_scripts/produce_edep-paramreco_larbath_transrots_tdr.sh`. 
The hadronic problem is similar, but the fix is more complicated. If we simulate the ND as 
a vat of LAr, then the dead region is not dense enough, so on average the hadrons travel 
too far in simulation, depositing too much energy in active regions, resulting in a reco. 
energy that's too high. While we don't have a full fix of this, we have a two-step partial 
fix. The first step is the Ehad density correction, where we add in the energy of the inactive 
hits, but with weights applied to account for the density differences. This is done in 
`ndfd_depos/nd-sim-tools/inputs/sim_inputs_larbath_selected_ndfd_pairs/dumpTree_tdr_nogeoeff_larbath.py`, 
lines 253-634. This doesn't totally take care of the unrealistic reconstruction, so we also 
add an *ad hoc* correction. The code and instructions for applying this are in the `nd_ehad_corr` 
directory of `lep_swapper`.

### `DiffEdep` branch

After failing to perfectly recover the ND performance using the hadronic energy corrections, 
we decided instead to start experimenting with a dataset that breaks the symmetry in energy 
deposits between the two detectors, hence the name 'Diff[erent]Edep'. In this training dataset, 
we take the same GENIE event, and place it once in the ND and once in the FD, ensuring that 
the energy deposits in each case are realistic and that we get the right reconstruction 
performance. This branch also includes an additional branch in the paired dataset, 
`nd_paramreco_part`, which contains true FS particle information from each event for the 
purpose of benchmarking the NDFD transformer.

## Moving data around

To make use of both Fermigrid and GPU resources I had to copy files to and from dCache to an
external cluster. To make these copies safely we need to use `ifdh cp`. Doing this not from a
dunegpvm can be tricky, I have found two ways to do this. The "Alternative Method" is probably
best to try first.

### Setup

```
git clone git@github.com:fermitools/cigetcert.git
cd cigetcert
git checkout tags/1.21
cd ../
export APPTAINER_BIND="/cvmfs,/home,/tmp,..."
apptainer shell /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7\:latest
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup python v3_9_2
python -m venv .venv.3_9_2_cigetcert
source .venv.3_9_2_cigetcert/bin/activate
pip install swig lxml M2Crypto pykerberos pyOpenSSL
```

### Copying

Put your krb5_fnal.conf on your home directory

```
export APPTAINER_BIND="/cvmfs,/home,/tmp"
apptainer shell /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7\:latest
export KRB5_CONFIG=/home/awilkins/krb5_fnal.conf
kinit
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
source .venv.3_9_2_cigetcert/bin/activate
python cigetcert/cigetcert -i 'Fermi National Accelerator Laboratory' -n
voms-proxy-init -noregen -rfc -voms dune:/dune/Role=Analysis
setup ifdhc
ifdh cp /pnfs/dune/persistent/...
```

### Alternative Method

The method above just stopped working for me at some point so here is an alternate and maybe better
one:

```
export APPTAINER_BIND="/cvmfs,/home,/tmp,/run"
apptainer shell /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7\:latest
export KRB5_CONFIG=/home/awilkins/krb5_fnal.conf
kinit
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
# May ask you to paste a link into browser to authenticate
htgettoken -a htvaultprod.fnal.gov -i dune
# Your NUID will be printed to terminal from the above command in this format
# Storing vault token in /tmp/vt_uNUID
# Storing bearer token in /run/user/NUID/bt_uNUID
export BEARER_TOKEN_FILE=/run/user/<NUID>/bt_u<NUID>
setup ifdhc
export IFDH_PROXY_ENABLE=0
export IFDH_TOKEN_ENABLE=1
ifdh cp /pnfs/dune/persistent/...
```

## Using Output Data

The final output comes from the `fd_detsim_reco` stage and can be used to train machine learning algorithms that map from ND detector response to FD reconstrution or to FD detector response. To prepare a dataset, iterate over the hd5f files and write the desired data to the desired format for input to a model. Note that when iterating the hdf5 files it is possible (though rare) that `fd_resp` or `fd_reco` will be missing for an event, in which case just skip this event.

