output_dir=./output
[ ! -d "$output_dir" ] && mkdir -p "$output_dir"
python ../add_substrate.py \
    --substrate_pdb ../anc451_top_pose.pdb \
    --substrate_str ../anc_substrate.str \
    --protein_pdb ../../add_5cthioester/5cthioester/caze_5cthioester.pdb \
    --protein_psf ../../add_5cthioester/5cthioester/caze_5cthioester.psf \
    --C5_str ../5cthioester.str \
    --C5_segname 5CT1 \
    --output_pdb $output_dir/caze_5cthioster_anc_substrate.pdb \
    --output_psf $output_dir/caze_5cthioster_anc_substrate.psf \
    --CHARMM_LIB_DIR /home/azamh/charmm/08302023/build_pyCHARMM/install/lib \
    --toppar_dir ../../toppar \
    > $output_dir/add_substrate.log 2>&1