#!/bin/bash
#SBATCH --job-name=dynamics
#SBATCH --output=dynamics_%a.out
#SBATCH --partition=gpuA5500
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --time=3-00:00:00
#SBATCH --mem=16GB
#SBATCH --array=0-3

# Activate conda environment
source ~/.bashrc
micromamba activate anc_acyl
module load cgenff/2023.1
module load mmtsb/rockyos/mmtsb/mmtsb

# Pick protein based on array idx
proteins=("caze" "anc451" "anc452" "anc454")
protein=${proteins[$SLURM_ARRAY_TASK_ID]}
echo "Protein: $protein"

# Format arguments
repo_dir='/home/azamh/workspace/anc_acyl'
solv_dir="$repo_dir/r_substrate/dyn/output/solvate/${protein}"
dyn_dir="$repo_dir/r_substrate/dyn/output/${protein}"
mkdir -p "$dyn_dir"

input_pdb="$solv_dir/solvated.pdb"
input_psf="$solv_dir/solvated.psf"
c5_str="$repo_dir/add_5cthioester/5cthioester.str"
ligand_str="$repo_dir/r_substrate/substrate/r_substrate.str"
workdir="$dyn_dir"
ion_type="CLA"
nequil=250000
nprod=15000000
nsavc=30000
toppar_dir="$repo_dir/toppar"
CHARMM_LIB_DIR='/home/azamh/charmm/pycharmm_2024/pycharmm_install/lib'
output_log="$workdir/dyn.log"

# Run dyn.py
script_path="$repo_dir/dyn/dyn.py"
python "$script_path" \
    --input_pdb "$input_pdb" \
    --input_psf "$input_psf" \
    --c5_str "$c5_str" \
    --ligand_str "$ligand_str" \
    --workdir "$workdir" \
    --ion_type "$ion_type" \
    --nequil "$nequil" \
    --nprod "$nprod" \
    --nsavc "$nsavc" \
    --toppar_dir "$toppar_dir" \
    --CHARMM_LIB_DIR "$CHARMM_LIB_DIR" \
    > "$output_log" 2>&1

# Dynamics simulation parameters in manuscript:
# --nequil 250000 --nprod 15000000 --nsavc 30000