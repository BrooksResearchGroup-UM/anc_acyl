#!/bin/bash
#SBATCH --job-name=solvate
#SBATCH --output=solvate_%a.out
#SBATCH --partition=gpuA5500
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --time=0-01:00:00
#SBATCH --mem=16GB
#SBATCH --array=2

# Activate conda environment
source ~/.bashrc
micromamba activate anc_acyl
module load cgenff/2023.1
module load mmtsb/rockyos/mmtsb/mmtsb

# Pick protein based on array idx
proteins=("caze" "anc451" "anc452" "anc454")
cleave_resis=(26 20 23 20)
protein=${proteins[$SLURM_ARRAY_TASK_ID]}
cleave_resi=${cleave_resis[$SLURM_ARRAY_TASK_ID]}
echo "Protein: $protein"
echo "Cleave residue: $cleave_resi"

repo_dir='/home/azamh/workspace/anc_acyl'
input_pdb="$repo_dir/r_substrate/add_substrate/output/${protein}_5cthioester_r_substrate.pdb" 
input_psf="$repo_dir/r_substrate/add_substrate/output/${protein}_5cthioester_r_substrate.psf"
C5_str="$repo_dir/add_5cthioester/5cthioester.str"
ligand_str="$repo_dir/r_substrate/substrate/r_substrate.str"
workdir="$repo_dir/r_substrate/dyn/output/solvate/${protein}"
mkdir -p "$workdir"
toppar_dir="$repo_dir/toppar"
CHARMM_LIB_DIR='/home/azamh/charmm/pycharmm_2024/pycharmm_install/lib'
output_log="$workdir/solvate.log"

# Run solvate.py
script_path="$repo_dir/dyn/solvate.py"
python "$script_path" \
    --input_pdb "$input_pdb" \
    --input_psf "$input_psf" \
    --C5_str "$C5_str" \
    --ligand_str "$ligand_str" \
    --workdir "$workdir" \
    --cleave_resi "$cleave_resi" \
    --toppar_dir "$toppar_dir" \
    --CHARMM_LIB_DIR "$CHARMM_LIB_DIR" \
    > "$output_log" 2>&1