output_dir=./output_dyn
[ ! -d "$output_dir" ] && mkdir -p "$output_dir"
python ../dyn.py \
    --input_pdb ../solvate/caze_5cthioester_anc_substrate/solvated.pdb \
    --input_psf ../solvate/caze_5cthioester_anc_substrate/solvated.psf \
    --c5_str ../5cthioester.str \
    --ligand_str ../anc_substrate.str \
    --workdir $output_dir \
    --ion_type CLA \
    --nequil 2500 \
    --nprod 150000 \
    --nsavc 300 \
    --toppar_dir ../../toppar \
    --CHARMM_LIB_DIR /home/azamh/charmm/08302023/build_pyCHARMM/install/lib \
    > $output_dir/dyn.log 2>&1

# Dynamics simulation parameters in manuscript:
# --nequil 250000 --nprod 15000000 --nsavc 30000