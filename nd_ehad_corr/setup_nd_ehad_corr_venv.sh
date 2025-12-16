venv_dir_name=".venv_3_9_2_nd_ehad_corr"

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup python v3_9_2
setup ifdhc
mkdir -p ${venv_dir_name}
python -m venv ${venv_dir_name}
source ${venv_dir_name}/bin/activate
pip install argparse h5py polars numpy matplotlib lmfit pandas scikit-learn
