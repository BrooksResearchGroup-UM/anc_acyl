"""Solvate an input structure with water molecules."""

import os
import argparse
import numpy as np

# region Parse arguments
parser = argparse.ArgumentParser(
    description="Solvate an input structure with water molecules."
)


# Helper function to get absolute path
def abspath(path):
    return os.path.abspath(path)


# Required arguments
parser.add_argument(
    "--input_pdb",
    type=abspath,
    help="Receptor PDB file",
    required=True,
)
parser.add_argument(
    "--input_psf",
    type=abspath,
    help="Receptor PSF file",
    required=True,
)
parser.add_argument(
    "--C5_str",
    type=abspath,
    help="5CThioester RTF file",
    required=True,
)
parser.add_argument(
    "--ligand_str",
    type=abspath,
    help="Ligand RTF file",
    required=True,
)
parser.add_argument(
    "--workdir",
    type=abspath,
    help="Working directory",
)
parser.add_argument(
    "--cleave_resi",
    type=str,
    help="Residue to cleave up to from the N-terminus",
    required=True,
)

parser.add_argument(
    "--toppar_dir",
    type=abspath,
    help="Directory containing topology and parameter files",
    required=True,
)
parser.add_argument(
    "--CHARMM_LIB_DIR",
    type=abspath,
    help="CHARMM library directory",
    required=True,
)

args = parser.parse_args()

# Print arguments
print("Arguments:")
for arg in vars(args):
    print(arg, getattr(args, arg))

# endregion Parse arguments

# Create working directory
os.makedirs(args.workdir, exist_ok=True)

# Import pyCHARMM
os.environ["CHARMM_LIB_DIR"] = args.CHARMM_LIB_DIR
import pycharmm
import pycharmm.generate as gen
import pycharmm.ic as ic
import pycharmm.coor as coor
import pycharmm.energy as energy
import pycharmm.dynamics as dyn
import pycharmm.nbonds as nbonds
import pycharmm.minimize as minimize
import pycharmm.crystal as crystal
import pycharmm.image as image
import pycharmm.psf as psf
import pycharmm.read as read
import pycharmm.write as write
import pycharmm.settings as settings
import pycharmm.cons_harm as cons_harm
import pycharmm.cons_fix as cons_fix
import pycharmm.select as select
import pycharmm.shake as shake
import pycharmm.settings as settings
from pycharmm.select_atoms import SelectAtoms
from pycharmm.lingo import charmm_script
from pycharmm.lib import charmm as libcharmm


def setup_pbc(boxhalf, solvated=True, ion_type=None, blade=False):
    """Apply periodic boundary conditions to a system.

    Args:
        boxhalf (float): Half the box size.
        solvated (bool, optional): Whether the system is solvated. Defaults to True.
        ion_type (str, optional): The ion type. Defaults to None.

    Assumes solute is segid LIGA and DUM and solvent is resname TIP3.
    Protein is segid PROA, PROB, etc.
    """
    crystal.define_cubic(boxhalf * 2)
    crystal.build(boxhalf)

    if blade:
        boxhalf = 0.0  # center at origin for blade

    # Turn on image centering - bysegment for protein, by residue for solvent and ions
    for segid in psf.get_segid():
        segid = segid.upper()
        if segid[0:3] in ["LIG", "DUM", "PRO", "5CT"]:
            print(f"Setting up segment {segid}")
            image.setup_segment(boxhalf, boxhalf, boxhalf, segid)

    if ion_type is not None:
        print(f"Setting up residue {ion_type}")
        image.setup_residue(boxhalf, boxhalf, boxhalf, ion_type)

    if solvated:
        print("Setting up residue TIP3")
        image.setup_residue(boxhalf, boxhalf, boxhalf, "TIP3")

    return


# Read in topology and parameter files
settings.set_bomb_level(-1)
read.rtf(os.path.join(args.toppar_dir, "top_all36_prot.rtf"))
read.rtf(os.path.join(args.toppar_dir, "top_all36_cgenff.rtf"), append=True)
read.prm(os.path.join(args.toppar_dir, "par_all36m_prot.prm"), flex=True)
read.prm(os.path.join(args.toppar_dir, "par_all36_cgenff.prm"), flex=True, append=True)
charmm_script(f"stream {args.C5_str}")
charmm_script(f"stream {args.ligand_str}")
charmm_script(f'stream {os.path.join(args.toppar_dir, "toppar_water_ions.str")}')
settings.set_bomb_level(0)


# Read in input structure
read.psf_card(args.input_psf)
read.pdb(args.input_pdb, resid=True)

# Delete first n flexible residues
delete_sel = None
protein_sel = SelectAtoms().by_seg_id("PROA")
for i in range(1, int(args.cleave_resi)):
    res_sel = SelectAtoms().by_res_id(str(i)) & protein_sel
    if delete_sel is None:
        delete_sel = res_sel
    else:
        delete_sel = delete_sel | res_sel

psf.delete_atoms(delete_sel)
charmm_script(f"patch ACE PROA {args.cleave_resi}")
charmm_script("IC GENErate")
ic.prm_fill(replace_all=True)
ic.build()

# Orient system
coor.orient(by_rms=False, by_mass=False, by_noro=False)

# Write pdb to working directory
input_pdb = os.path.join(args.workdir, "input.pdb")
write.coor_pdb(input_pdb)

# Setup nonbonds
my_nbonds = pycharmm.NonBondedScript(
    cutnb=12,
    ctonnb=8,
    ctofnb=10,
    eps=1.0,
    cdie=True,
    atom=True,
    vatom=True,
    switch=True,
    vfswitch=False,
    vswitch=False,
    noEwald=True,
)
my_nbonds.run()

# Calculate total charge to add neutralizing ions
q = psf.get_charges()
Ntot = round((np.sum(q)))
if Ntot > 0:
    ion_type = "CLA"
if Ntot < 0:
    ion_type = "SOD"
if Ntot != 0:
    ions = "-ions {}:{}".format(ion_type, np.abs(Ntot))
if np.abs(Ntot) < 1e-2:
    ions = ""
    ion_type = None
print(f"Adding {ions} to neutralize the system.")


# Solvate system
water_ion_pdb = os.path.join(args.workdir, "water_ion.pdb")
water_pdb = os.path.join(args.workdir, "wt00.pdb")
ion_pdb = os.path.join(args.workdir, "ion.pdb")
solvate_command = f"convpdb.pl -solvate -cutoff 10 {ions} -cubic -out charmm22 {input_pdb} > {water_ion_pdb};"
solvate_command += f"convpdb.pl -segnames -nsel TIP3 {water_ion_pdb} > {water_pdb};"
solvate_command += f"convpdb.pl -segnames -nsel ion {water_ion_pdb} > {ion_pdb}"
os.system(solvate_command)

# replace HETATM by ATOM in ions
ion_segid = "IONA"
new_ion_pdb = ""
with open(ion_pdb, "r") as f:
    for line in f.read().splitlines():
        line = line.replace("HETATM", "ATOM  ")
        if line.startswith("ATOM"):
            line = list(line)
            line[-4:] = ion_segid.split()
            line = "".join(line)
        new_ion_pdb += line + "\n"
with open(ion_pdb, "w") as f:
    f.write(new_ion_pdb)

# Replace water segid with WT00
water_segid = "WT00"
new_water_pdb = ""
with open(water_pdb, "r") as f:
    for line in f.read().splitlines():
        if line.startswith("ATOM"):
            line = list(line)
            line[-4:] = water_segid.split()
            line = "".join(line)
        new_water_pdb += line + "\n"
with open(water_pdb, "w") as f:
    f.write(new_water_pdb)

# Read in water sequence
read.sequence_pdb(water_pdb)
gen.new_segment(water_segid, angle=False, dihedral=False)
read.pdb(water_pdb, resid=True)

# Read in ion sequence
if ions:
    read.sequence_pdb(ion_pdb)
    gen.new_segment(ion_segid, angle=False, dihedral=False)
    read.pdb(ion_pdb, resid=True)

# get the coor statistics to construct boxlengths
# coor stat
stats = coor.stat()

# boxsize
xsize = stats["xmax"] - stats["xmin"]
ysize = stats["ymax"] - stats["ymin"]
zsize = stats["zmax"] - stats["zmin"]
boxsize = max(xsize, ysize, zsize)

# half box size
boxhalf = boxsize / 2.0

print("Boxsize: {}".format(boxsize))
print("Boxhalf: {}".format(boxhalf))

# Get segids
segids = psf.get_segid()
print(f"PSF SEGIDS: {segids}")

# translate coordinates by boxhalf to be consistent with openmm
coor.set_positions(coor.get_positions() + boxhalf)

# # Set up periodic boundary conditions
# setup_pbc(boxhalf, ion_type=ion_type, solvated=True)

# Run Minimization
print("Minimizing Solvent...")
cons_fix.setup(
    selection=~(
        pycharmm.SelectAtoms(seg_id="WT00") | pycharmm.SelectAtoms(seg_id="IONA")
    )
)
minimize.run_sd(nstep=500, tolenr=1e-3, tolgrd=1e-3)
cons_fix.turn_off()

# Write boxhalf to file
with open(os.path.join(args.workdir, "boxhalf.txt"), "w") as f:
    f.write(str(boxhalf))

# Write ion type to file
with open(os.path.join(args.workdir, "ion_type.txt"), "w") as f:
    f.write(str(ion_type))

# Write psf and pdb
output_pdb = os.path.join(args.workdir, "solvated.pdb")
output_psf = os.path.join(args.workdir, "solvated.psf")
write.coor_pdb(output_pdb)
write.psf_card(output_psf)
