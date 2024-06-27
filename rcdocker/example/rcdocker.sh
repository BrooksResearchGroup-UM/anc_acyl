output_dir=./output
[ ! -d "$output_dir" ] && mkdir -p "$output_dir"
python ../rcdocker.py \
    --receptor_pdb ../../add_5cthioester/5cthioester/caze_5cthioester.pdb \
    --receptor_psf ../../add_5cthioester/5cthioester/caze_5cthioester.psf \
    --c5_str ../../add_5cthioester/5cthioester.str \
    --substrate_pdb ../anc_substrate.pdb \
    --substrate_str ../anc_substrate.str \
    --substrate_segname ANC1 \
    --xcen 26.711349487304688 \
    --ycen 4.538999557495117 \
    --zcen 65.70699691772461 \
    --work_dir $output_dir \
    --charmm_lib_dir /home/azamh/charmm/08302023/build_pyCHARMM/install/lib \
    --toppar_dir ../../toppar \
    > $output_dir/rcdocker.log 2>&1