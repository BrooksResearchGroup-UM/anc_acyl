"""Minimize 5cthioester group in vacuum using pyCHARMM"""

import os
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser()

# Required arguments

# Input protein pdb
parser.add_argument(
    "--protein_pdb",
    help="Input protein pdb",
    required=True,
)

# Input protein psf
parser.add_argument(
    "--protein_psf",
    help="Input protein psf",
    required=True,
)

# Input prosthetic pdb
parser.add_argument(
    "--c5_pdb",
    help="Input 5cthioester pdb",
    required=True,
)

# Input prosthetic str
parser.add_argument(
    "--c5_str",
    help="Input 5cthioester str",
    required=True,
)

# Input prosthetic segname
parser.add_argument(
    "--c5_segname",
    help="Input 5cthioester segname",
    required=True,
)

# Output pdb file
parser.add_argument(
    "--output_pdb",
    help="Output pdb file",
    required=True,
)

# Output psf file
parser.add_argument(
    "--output_psf",
    help="Output psf file",
    required=True,
)

# Output energy file
parser.add_argument(
    "--output_energy",
    help="Output energy file",
    required=True,
)

# Path to CHARMM_LIB_DIR
parser.add_argument(
    "--charmm_lib_dir",
    help="Path to CHARMM_LIB_DIR",
    required=True,
)

# Path to toppar directory
parser.add_argument(
    "--toppar_dir",
    help="Path to toppar directory",
    required=True,
)

# Optional arguments

# Protein segname
parser.add_argument(
    "--protein_segname",
    help="Protein segname",
    default="PROA",
)

# Number of SD steps
parser.add_argument(
    "--nstep_sd",
    help="Number of SD steps",
    default=200,
)

# Number of ABNR steps
parser.add_argument(
    "--nstep_abnr",
    help="Number of ABNR steps",
    default=1000,
)

# Parse arguments
args = parser.parse_args()

# Print arguments
print("Arguments:")
for arg in vars(args):
    print(arg, getattr(args, arg))

# Import pyCHARMM
os.environ["CHARMM_LIB_DIR"] = args.charmm_lib_dir
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

# Read in topology and parameter files
settings.set_bomb_level(-1)
read.rtf(os.path.join(args.toppar_dir, "top_all36_prot.rtf"))
read.rtf(os.path.join(args.toppar_dir, "top_all36_cgenff.rtf"), append=True)
read.prm(os.path.join(args.toppar_dir, "par_all36m_prot.prm"), flex=True)
read.prm(os.path.join(args.toppar_dir, "par_all36_cgenff.prm"), flex=True, append=True)
charmm_script(f"stream {args.c5_str}")
settings.set_bomb_level(0)

# Read in protein
read.psf_card(args.protein_psf)
read.pdb(args.protein_pdb, resid=True)

# Generate prosthetic psf
read.sequence_pdb(args.c5_pdb)
gen.new_segment(seg_name=args.c5_segname, setup_ic=True)
read.pdb(args.c5_pdb, resid=True)
charmm_script("auto angle dihe")
write.psf_card(args.output_psf)

# Setup nonbonds
my_nbonds = pycharmm.NonBondedScript(
    cutnb=10.0,
    ctonnb=8.0,
    ctofnb=7.0,
    eps=1.0,
    cdie=True,
    atom=True,
    vatom=True,
    fswitch=True,
    vfswitch=True,
)
# Implement these non-bonded parameters by "running" them.
my_nbonds.run()
energy.show()

# Select backbone atoms and flexible side chains around prosthetic group
protein_sel = SelectAtoms(seg_id=args.protein_segname) & ~SelectAtoms(hydrogens=True)
prosthetic_sel = SelectAtoms(seg_id=args.c5_segname)
all_sel = SelectAtoms(select_all=True)

backbone_sel = protein_sel & (
    SelectAtoms(atom_type="CA")
    | SelectAtoms(atom_type="N")
    | SelectAtoms(atom_type="C")
    | SelectAtoms(atom_type="O")
)
flexres_sel = (
    pycharmm.SelectAtoms(seg_id=args.protein_segname)
    & pycharmm.SelectAtoms(seg_id=args.c5_segname).around(5)
).whole_residues()
flexres_sc_sel = flexres_sel & (~backbone_sel)

print("all atoms:", all_sel.get_n_selected())
print("protein atoms:", protein_sel.get_n_selected())
print("prosthetic atoms:", prosthetic_sel.get_n_selected())
print("backbone atoms", backbone_sel.get_n_selected())
print("flexible residues around prosthetic_sel atoms:", flexres_sel.get_n_selected())
print(
    "flexible residues sidechains around prosthetic_sel atoms:",
    flexres_sc_sel.get_n_selected(),
)

# restrain backbone atoms first letting sidechains and prosthetic groups move
cons_harm.setup_absolute(selection=backbone_sel, force_constant=50)
minimize.run_sd(nstep=args.nstep_sd, tolenr=1e-3, tolgrd=1e-4)
cons_harm.turn_off()

# restrain all atoms except nearby side chains and prosthetic group
cons_harm.setup_absolute(selection=((~flexres_sc_sel) & protein_sel), force_constant=50)
minimize.run_sd(nstep=args.nstep_sd, tolenr=1e-3, tolgrd=1e-4)
cons_harm.turn_off()

# restrain everything except prosthetic group
cons_harm.setup_absolute(selection=protein_sel, force_constant=50)
minimize.run_sd(nstep=args.nstep_sd, tolenr=1e-3, tolgrd=1e-4)
minimize.run_abnr(nstep=args.nstep_abnr, tolenr=1e-3, tolgrd=1e-4)
cons_harm.turn_off()

# print energy
energy.show()

# Compute interaction energy
charmm_script(
    f"INTE sele segid {args.c5_segname} end sele segid {args.protein_segname} end"
)

# Write energy file
charmm_script(
    f"""
open unit 1 write form name {args.output_energy}
echu 1
echo TOT: ?ener
echo VDW: ?vdw
echo ELEC: ?elec
"""
)

# Save pdb
write.coor_pdb(args.output_pdb)
