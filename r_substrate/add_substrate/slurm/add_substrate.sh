#!/bin/bash
#SBATCH --job-name=add_substrate
#SBATCH --output=add_substrate_%a.out
#SBATCH --partition=gpuA5500
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:0
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
protein=${proteins[$SLURM_ARRAY_TASK_ID]}
echo "Protein: $protein"

# Set arguments
repo_dir='/home/azamh/workspace/anc_acyl'
substrate_pdb="$repo_dir/r_substrate/substrate/r_substrate.pdb"
substrate_str="$repo_dir/r_substrate/substrate/r_substrate.str"
protein_pdb="$repo_dir/add_5cthioester/5cthioester/${protein}_5cthioester.pdb"
protein_psf="$repo_dir/add_5cthioester/5cthioester/${protein}_5cthioester.psf"
C5_str="$repo_dir/add_5cthioester/5cthioester.str"
C5_segname='5CT1'
output_pdb="$repo_dir/r_substrate/add_substrate/output/${protein}_5cthioester_r_substrate.pdb"
output_psf="$repo_dir/r_substrate/add_substrate/output/${protein}_5cthioester_r_substrate.psf"
CHARMM_LIB_DIR='/home/azamh/charmm/pycharmm_2024/pycharmm_install/lib'
toppar_dir="$repo_dir/toppar"
output_log="$repo_dir/r_substrate/add_substrate/output/${protein}_5cthioester_r_substrate.log"

# Run add_substrate.py 
script_path="$repo_dir/add_substrate/add_substrate.py"
python $script_path \
    --substrate_pdb "$substrate_pdb" \
    --substrate_str "$substrate_str" \
    --protein_pdb "$protein_pdb" \
    --protein_psf "$protein_psf" \
    --C5_str "$C5_str" \
    --C5_segname "$C5_segname" \
    --output_pdb "$output_pdb" \
    --output_psf "$output_psf" \
    --CHARMM_LIB_DIR "$CHARMM_LIB_DIR" \
    --toppar_dir "$toppar_dir" \
    > "$output_log" 2>&1