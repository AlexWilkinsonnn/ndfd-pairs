# nd_ehad_corr

This is the final step in the NDFD Transformer training dataset generation procedure. It takes the paired .h5 dataset output from the `fd_detsim_reco` module, which already has a density-based correction to the ND hadronic energy reconstruction, and applies one additional correction of the form $$E_\text{corr.} = \frac{m_\text{ND}}{m_\text{LAr}} * (E - b_\text{LAr}) + b_\text{ND}$$. A correction of this form is applied to each of the ND reco. hadronic energy variables in the training dataset except "Other", which is set so that the reconstructed hadronic energy at the ND in each event is equal to the sum of the reconstructed hadronic energy for each particle type. A correction is not applied if it would make the reco. energy negative. Each particle type gets its own correction parameters, which are set using a guess-and-check method to approximately minimize the $\chi^2_\text{red.}$ between the training dataset and a reference CAF file.

For the FHC on-axis orientation, the following parameters are used for the shift (format: $m_\text{LAr}$, $b_\text{LAr}$, $m_\text{ND}$, $b_\text{ND}$):
- Hadronic energy (which affects `Ev_reco`): [0.67, 0.001, 0.665, 0.0]
- `eRecoP`: [0.86, 0.0001, 0.8599, 0.0]
- `eRecoN`: [0.42, 0.002, 0.4, 0.0]
- `eRecoPip`: [0.79, 0.05, 0.78, 0.05]
- `eRecoPim`: [0.72, 0.082, 0.71, 0.08]
- `eRecoPi0`: [0.86, 0.13, 0.858, 0.13]

The progression of $\chi^2_\text{red/}$ from no correction, to the density correction, to the *ad hoc* correction, is:
- `Ev`: $6.818 \rightarrow 2.297 \rightarrow 1.758$
- Hadronic energy: $7.395 \rightarrow 1.954 \rightarrow 1.521$
- `eP`: $18.82 \rightarrow 2.03 \rightarrow 1.915$
- `eN`: $21.35 \rightarrow 12.18 \rightarrow 2.028$
- `ePip`: $6.732 \rightarrow 1.926 \rightarrow 1.574$
- `ePim`: $2.632 \rightarrow 1.62 \rightarrow 1.487$
- `ePi0`: $15.79 \rightarrow 2.459 \rightarrow 2.283$

For reference, $\chi^2_\text{red.}$ comparing the FD reco to true variables, and the ND leptonic energy, are
- FD `Ev`: 1.091
- FD `Elep`: 1.371
- FD `Ehad`: 1.72
- ND `Elep: 1.172

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
