"""Minimize input pdb file in vacuum using pyCHARMM"""

import os
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser()

# Required arguments

# Input pdb file
parser.add_argument("-i", "--input", help="Input pdb file", required=True)

# Output pdb file
parser.add_argument("-o", "--output", help="Output pdb file", required=True)

# Output psf file
parser.add_argument("-p", "--psf", help="Output psf file", required=True)

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

# Parse arguments
args = parser.parse_args()

# Print arguments
print("Arguments:")
for arg in vars(args):
    print(arg, getattr(args, arg))

# Validate input arguments
assert args.input.endswith(".pdb"), "Input file must be a pdb file"
assert os.path.exists(args.input), "Input file does not exist"
assert args.output.endswith(".pdb"), "Output file must be a pdb file"
assert args.psf.endswith(".psf"), "Output file must be a psf file"
assert os.path.exists(args.charmm_lib_dir), "CHARMM_LIB_DIR does not exist"
assert os.path.exists(args.toppar_dir), "toppar directory does not exist"

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
read.rtf(os.path.join(args.toppar_dir, "top_all36_prot.rtf"))
read.prm(os.path.join(args.toppar_dir, "par_all36m_prot.prm"), flex=True)

# Read in protein sequence and generate psf
read.sequence_pdb(args.input)
gen.new_segment(seg_name="PROA", first_patch="ACE", last_patch="CT3", setup_ic=True)
read.pdb(args.input, resid=True)
ic.prm_fill(replace_all=False)
ic.build()
write.psf_card(args.psf)

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

# Impose restraints on heavy atoms
cons_harm.setup_absolute(
    selection=~pycharmm.SelectAtoms(hydrogens=True), force_constant=50
)
# minimize with steepest descent
minimize.run_sd(nstep=1000, tolenr=1e-3, tolgrd=1e-4)
cons_harm.turn_off()

# Save minimized pdb
write.coor_pdb(args.output)
