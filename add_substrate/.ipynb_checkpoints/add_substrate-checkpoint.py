"""Minimize substrate in vacuum using pyCHARMM"""

import os
import argparse
import pandas as pd
import numpy as np

# region Parse arguments
parser = argparse.ArgumentParser(
    description="Minimize substrate in vacuum using pyCHARMM"
)


# Helper function to get absolute path
def abspath(path):
    return os.path.abspath(path)


# Required arguments
parser.add_argument(
    "--substrate_pdb",
    type=abspath,
    help="Ligand PDB file",
    required=True,
)
parser.add_argument(
    "--substrate_str",
    type=abspath,
    help="Ligand RTF file",
    required=True,
)
parser.add_argument(
    "--protein_pdb",
    type=abspath,
    help="Receptor PDB file",
    required=True,
)
parser.add_argument(
    "--protein_psf",
    type=abspath,
    help="Receptor PSF file",
    required=True,
)
parser.add_argument(
    "--C5_str",
    type=abspath,
    help="C5_ RTF file",
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

# Input 5cthioester segname
parser.add_argument(
    "--C5_segname",
    help="Input 5cthioester segname",
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

# Protein segname
parser.add_argument(
    "--protein_segname",
    help="Protein segname",
    default="PROA",
)

# Input substrate segname
parser.add_argument(
    "--substrate_segname",
    help="Input substrate segname",
    default="LIGA",
)

args = parser.parse_args()

# Print arguments
print("Arguments:")
for arg in vars(args):
    print(arg, getattr(args, arg))

# endregion Parse arguments

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

# Read in topology and parameter files
settings.set_bomb_level(-1)
read.rtf(os.path.join(args.toppar_dir, "top_all36_prot.rtf"))
read.rtf(os.path.join(args.toppar_dir, "top_all36_cgenff.rtf"), append=True)
read.prm(os.path.join(args.toppar_dir, "par_all36m_prot.prm"), flex=True)
read.prm(os.path.join(args.toppar_dir, "par_all36_cgenff.prm"), flex=True, append=True)
charmm_script(f"stream {args.C5_str}")
charmm_script(f"stream {args.substrate_str}")
settings.set_bomb_level(0)

# # Read in protein
read.psf_card(args.protein_psf)
read.pdb(args.protein_pdb, resid=True)

# Generate substrate psf
read.sequence_pdb(args.substrate_pdb)
gen.new_segment(seg_name=args.substrate_segname, setup_ic=True)
read.pdb(args.substrate_pdb, resid=True)
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

# Select backbone atoms and flexible side chains around 5cthioester group
protein_sel = SelectAtoms(seg_id=args.protein_segname) & ~SelectAtoms(hydrogens=True)
C5_sel = SelectAtoms(seg_id=args.C5_segname)
substrate_sel = SelectAtoms(seg_id=args.substrate_segname)
all_sel = SelectAtoms(select_all=True)

backbone_sel = protein_sel & (
    SelectAtoms(atom_type="CA")
    | SelectAtoms(atom_type="N")
    | SelectAtoms(atom_type="C")
    | SelectAtoms(atom_type="O")
)

flexres_sel = (
    (protein_sel | C5_sel) & pycharmm.SelectAtoms(seg_id=args.C5_segname).around(5)
).whole_residues()
flexres_sc_sel = flexres_sel & (~backbone_sel)

print("all atoms:", all_sel.get_n_selected())
print("protein atoms:", protein_sel.get_n_selected())
print("5cthioester atoms:", C5_sel.get_n_selected())
print("substrate atoms:", substrate_sel.get_n_selected())
print("backbone atoms", backbone_sel.get_n_selected())
print("flexible residues around substrate_sel atoms:", flexres_sel.get_n_selected())
print(
    "flexible residues sidechains around substrate_sel atoms:",
    flexres_sc_sel.get_n_selected(),
)


# first restrain backbone atoms and substrate and letting sidechains and C5_ groups move
cons_harm.setup_absolute(selection=backbone_sel | substrate_sel, force_constant=50)
minimize.run_sd(nstep=args.nstep_sd, tolenr=1e-3, tolgrd=1e-4)
cons_harm.turn_off()

# restrain all atoms except nearby side chains and C5_ group
cons_harm.setup_absolute(selection=((~flexres_sc_sel) & protein_sel), force_constant=50)
minimize.run_sd(nstep=args.nstep_sd, tolenr=1e-3, tolgrd=1e-4)
cons_harm.turn_off()

# restrain everything except substrate
cons_harm.setup_absolute(selection=~substrate_sel, force_constant=50)
minimize.run_sd(nstep=args.nstep_sd, tolenr=1e-3, tolgrd=1e-4)
minimize.run_abnr(nstep=args.nstep_abnr, tolenr=1e-3, tolgrd=1e-4)
cons_harm.turn_off()

# print energy
energy.show()


# Save pdb
write.coor_pdb(args.output_pdb)
