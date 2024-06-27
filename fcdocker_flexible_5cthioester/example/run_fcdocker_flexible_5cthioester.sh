output_dir=./output
[ ! -d "$output_dir" ] && mkdir -p "$output_dir"
python ../run_fcdocker_flexible_5cthioester.py \
    --ligPDB ../anc_substrate.pdb \
    --ligandrtf ../anc_substrate.str \
    --receptorPDB ../../add_5cthioester/5cthioester/caze_5cthioester.pdb \
    --receptorPSF ../../add_5cthioester/5cthioester/caze_5cthioester.psf \
    --flexchain_csv ../flexres/caze_5cthioester_anc_substrate.csv \
    --C5_rtf ../../add_5cthioester/5cthioester.str \
    --C5_segname 5CT1 \
    --workdir $output_dir \
    --num 5 \
    --copy 1 \
    --xcen 26.711349487304688 \
    --ycen 4.538999557495117 \
    --zcen 65.70699691772461 \
    --maxlen 40 \
    --topdir ../../toppar \
    --CHARMM_LIB_DIR /home/azamh/charmm/08302023/build_pyCHARMM/install/lib \
    > "$output_dir"/run_fcdocker_flexible_5cthioester.log 2>&1

# Original values for num, copy, and maxlen:
# --num 25 \
# --copy 20 \
# --maxlen 40 \
