output_dir=./output
[ ! -d "$output_dir" ] && mkdir -p "$output_dir"
python ../add_5cthioester.py \
    --protein_pdb ../../vacuum_min/pdb_min/caze.pdb \
    --protein_psf ../../vacuum_min/pdb_min/caze.psf \
    --c5_pdb ../5cthioester.pdb \
    --c5_str ../5cthioester.str \
    --c5_segname 5CT1 \
    --output_pdb $output_dir/caze_5cthioester.pdb \
    --output_psf $output_dir/caze_5cthioester.psf \
    --output_energy $output_dir/caze_5cthioester_energy.txt \
    --charmm_lib_dir /home/azamh/charmm/08302023/build_pyCHARMM/install/lib \
    --toppar_dir ../../toppar \
    > $output_dir/add_5cthioester.log 2>&1