output_dir=./output_solvate
[ ! -d "$output_dir" ] && mkdir -p "$output_dir"
python ../solvate.py \
    --input_pdb ../../add_substrate/substrate/caze_5cthioester_anc_substrate.pdb \
    --input_psf ../../add_substrate/substrate/caze_5cthioester_anc_substrate.psf \
    --C5_str ../5cthioester.str \
    --ligand_str ../anc_substrate.str \
    --workdir $output_dir \
    --cleave_resi 26 \
    --toppar_dir ../../toppar \
    --CHARMM_LIB_DIR /home/azamh/charmm/08302023/build_pyCHARMM/install/lib \
    > $output_dir/solvate.log 2>&1

# Residues to cleave at by protein (to remove long N-terminal tail)
# anc451 20
# anc452 23
# anc454 20
# caze 26