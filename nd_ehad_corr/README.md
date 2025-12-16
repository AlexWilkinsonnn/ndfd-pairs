# nd_ehad_corr

This is the final step in the NDFD Transformer training dataset generation procedure. It takes the paired dataset output from the `fd_detsim_reco` module, which already has a density-based correction to the ND hadronic energy reconstruction, and applies one additional correction of the form $$E_\text{corr.} = \frac{m_\text{ND}}{m_\text{LAr}} * (E - b_\text{LAr}) + b_\text{ND}$$. A correction of this form is applied to each of the ND reco. hadronic energy variables in the training dataset except "Other", which is set so that the reconstructed hadronic energy at the ND in each event is equal to the sum of the reconstructed hadronic energy for each particle type. A correction is not applied if it would make the reco. energy negative. Each particle type gets its own correction parameters, which are set using a guess-and-check method to approximately minimize the $\chi^2_\text{red.}$ between the training dataset and a reference CAF file.

## System

Recommended to do this on a dunegpvm within an SL7 container. DUNE software should not yet be set up.

## Setup

```
cd nd_ehad_corr
source setup_nd_ehad_corr_venv.sh
```

## Running
```
python CorrNDEhadReco.py <initial paired dataset> <output paired dataset>

```
