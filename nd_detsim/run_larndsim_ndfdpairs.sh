#!/bin/bash
#SBATCH -p GPU
#SBATCH -N 1
#SBATCH -c 2
#SBATCH -J larnd-sim
#SBATCH -t 360
#SBATCH --array=1-183
#SBATCH --gres=gpu:1
#SBATCH --error=/home/awilkins/ndfd_pairs/logs/err/%x.%A_%a.err
#SBATCH --output=/home/awilkins/ndfd_pairs/logs/out/%x.%A_%a.out

################################################################################
# Options

LARNDSIM_WORK_DIR="/home/awilkins/ndfd_pairs/ndfd-pairs/nd_detsim/larnd-sim"
LARPIXSOFT_WORK_DIR="/home/awilkins/ndfd_pairs/ndfd-pairs/nd_detsim/larpixsoft"
SCRATCH_DIR="/state/partition1/awilkins/scratch/${SLURM_JOB_ID}"

INPUT_DIR=$1
OUTPUT_DIR=$2

################################################################################

mkdir -p ${SCRATCH_DIR}

input_name=$(ls $INPUT_DIR | head -n $SLURM_ARRAY_TASK_ID | tail -n -1)
input_file=${INPUT_DIR}/${input_name}
output_name_larndsim=${input_name%.*}_larndsim.h5
output_name_final=${input_name%.*}_larndsim_3dpackets.h5
output_file_larndsim=${SCRATCH_DIR}/${output_name_larndsim}
output_file_final=${SCRATCH_DIR}/${output_name_final}

echo "Job id ${SLURM_JOB_ID}"
echo "Job array task id ${SLURM_ARRAY_TASK_ID}"
echo "Running on ${SLURM_JOB_NODELIST}"
echo "With cuda device ${CUDA_VISIBLE_DEVICES}"
echo "Input file is ${input_file}"
echo "Output file will be ${output_file_final}"

cd $LARNDSIM_WORK_DIR
source setup.sh

python cli/simulate_pixels.py --input_filename $input_file \
                              --detector_properties larndsim/detector_properties/ndlar-module.yaml \
                              --pixel_layout larndsim/pixel_layouts/multi_tile_layout-3.0.40.yaml \
                              --output_filename $output_file_larndsim \
                              --response_file larndsim/bin/response_38.npy
# If larnd-sim command crashed during execution there will be missing data down the line,
# just throw away these event if this happens
exit_code=$?
if [ $exit_code -ne 0 ]; then
  rm -r ${SCRATCH_DIR}
  exit 1
fi

cd $LARPIXSOFT_WORK_DIR

python add_3d_toh5.py $output_file_larndsim \
                      $output_file_final \
                      ${LARNDSIM_WORK_DIR}/larndsim/detector_properties/ndlar-module.yaml \
                      ${LARNDSIM_WORK_DIR}/larndsim/pixel_layouts/multi_tile_layout-3.0.40.yaml

cp ${SCRATCH_DIR}/${output_name_final} ${OUTPUT_DIR}/${output_name_final}

rm -r ${SCRATCH_DIR}

