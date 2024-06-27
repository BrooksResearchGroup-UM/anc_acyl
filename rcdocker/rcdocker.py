"""Dock substrate to receptors using simulated annealing"""

import os
import shutil
import argparse


# Path type
def path(x):
    return os.path.abspath(x)


# Parse arguments
parser = argparse.ArgumentParser(
    description="Dock substrate to receptors using simulated annealing"
)

# Required arguments

# Receptor pdb
parser.add_argument(
    "--receptor_pdb",
    type=path,
    required=True,
)

# Receptor psf
parser.add_argument(
    "--receptor_psf",
    type=path,
    required=True,
)

# 5cthioester group str
parser.add_argument(
    "--c5_str",
    type=path,
    required=True,
)

# Substrate pdb
parser.add_argument(
    "--substrate_pdb",
    type=path,
    required=True,
)

# Substrate str
parser.add_argument(
    "--substrate_str",
    type=path,
    required=True,
)

# Substrate segname
parser.add_argument(
    "--substrate_segname",
    type=str,
    required=True,
)

# Xcen of binding site
parser.add_argument(
    "--xcen",
    type=float,
    required=True,
)

# Ycen of binding site
parser.add_argument(
    "--ycen",
    type=float,
    required=True,
)

# Zcen of binding site
parser.add_argument(
    "--zcen",
    type=float,
    required=True,
)

# Work directory for docking
parser.add_argument(
    "--work_dir",
    type=path,
    required=True,
)

# Path to CHARMM_LIB_DIR
parser.add_argument(
    "--charmm_lib_dir",
    help="Path to CHARMM_LIB_DIR",
    type=path,
    required=True,
)

# Path to toppar directory
parser.add_argument(
    "--toppar_dir",
    type=path,
    help="Path to toppar directory",
    required=True,
)

# Optional arguments

# Max length of grid
parser.add_argument(
    "--maxlen",
    type=int,
    default=20,
)

# Total number of conformers
parser.add_argument(
    "--num_conformer",
    type=int,
    default=10,
)

# Number of placements
parser.add_argument(
    "--num_place",
    type=int,
    default=100,
)

# Number of copies
parser.add_argument(
    "--num_copy",
    type=int,
    default=1000,
)

# Path to obrotamer executable
parser.add_argument(
    "--obrotamer",
    help="Path to obrotamer executable",
    type=path,
    default=shutil.which("obrotamer"),
)

# Parse arguments
args = parser.parse_args()

# Print arguments
print("Arguments:")
for arg in vars(args):
    print(arg, getattr(args, arg))


# Generate workdir
if not os.path.exists(args.work_dir):
    os.makedirs(args.work_dir)

# Change directory to workdir
os.chdir(args.work_dir)

# Directory for grids
grid_dir = os.path.join(args.work_dir, "grid/")
if not os.path.exists(grid_dir):
    os.makedirs(grid_dir)

# Directory for conformers
conformer_dir = os.path.join(args.work_dir, "conformer/")
if not os.path.exists(conformer_dir):
    os.makedirs(conformer_dir)

# Generate conformers if they don't exist
if not os.path.exists(os.path.join(conformer_dir, "1.pdb")):
    print("Generating conformers...")
    for r in range(1, args.num_conformer + 1):
        os.system(
            f'{args.obrotamer} {args.substrate_pdb} > {os.path.join(conformer_dir, str(r) + ".pdb")}'
        )

# Directory for placements
placement_dir = os.path.join(args.work_dir, "placement/")
if not os.path.exists(placement_dir):
    os.makedirs(placement_dir)

# Directory for dockresult
dockresult_dir = os.path.join(args.work_dir, "dockresult/")
if not os.path.exists(dockresult_dir):
    os.makedirs(dockresult_dir)

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
from pycharmm.cdocker import Rigid_CDOCKER

# Read in topology and parameter files
settings.set_bomb_level(-1)
read.rtf(os.path.join(args.toppar_dir, "top_all36_prot.rtf"))
read.rtf(os.path.join(args.toppar_dir, "top_all36_cgenff.rtf"), append=True)
read.prm(os.path.join(args.toppar_dir, "par_all36m_prot.prm"), flex=True)
read.prm(os.path.join(args.toppar_dir, "par_all36_cgenff.prm"), flex=True, append=True)
charmm_script(f"stream {args.c5_str}")
charmm_script(f"stream {args.substrate_str}")
settings.set_bomb_level(0)

# Format arguments for rigid cdocker
rigid_args = {
    "xcen": args.xcen,
    "ycen": args.ycen,
    "zcen": args.zcen,
    "maxlen": args.maxlen,
    "dielec": 3,
    "rcta": 0,
    "rctb": 0,
    "hmax": 0,
    "flag_grid": False,
    "flag_rdie": True,
    "flag_form": False,
    "flag_delete_grid": False,
    "probeFile": os.path.join(args.toppar_dir, "fftdock_c36prot_cgenff_probes.txt"),
    "softGridFile": os.path.join(grid_dir, "grid-emax-0.6-mine--0.4-maxe-0.4.bin"),
    "hardGridFile": os.path.join(grid_dir, "grid-emax-3-mine--30-maxe-30.bin"),
    "nativeGridFile": os.path.join(grid_dir, "grid-emax-100-mine--100-maxe-100.bin"),
    "receptorPDB": args.receptor_pdb,
    "receptorPSF": args.receptor_psf,
    "ligPDB": args.substrate_pdb,
    "ligSeg": args.substrate_segname,
    "confDir": conformer_dir,
    "placementDir": placement_dir,
    "exhaustiveness": "high",
    "numPlace": args.num_place,
    "numCopy": args.num_copy,
    "flag_delete_conformer": False,
    "flag_delete_placement": False,
    "flag_save_all": True,
    "flag_save_cluster": True,
    "flag_save_top": True,
    "flag_suppress_print": True,
    "flag_center_ligand": True,
    "flag_fast_grid": False,
    "flag_use_hbond": False,
    "flag_fast_placement": True,
    "threshold": 2500,
    "sort_energy": "total_energy",
    "saveDir": dockresult_dir,
}

print("Rigid docking arguments:")
for arg in rigid_args:
    print(arg, rigid_args[arg])

# print('Running rigid docking...')
Rigid_CDOCKER(**rigid_args)

# Documentation for Rigid_CDOCKER
# Rigid_CDOCKER(xcen = 0, ycen = 0, zcen = 0, maxlen = 10, dielec = 3,
#                   rcta = 0, rctb = 0, hmax = 0, flag_grid = False,
#                   flag_rdie = True, flag_form = False, flag_delete_grid = True,
#                   probeFile = '"../Toppar/fftdock_c36prot_cgenff_probes.txt"',
#                   softGridFile = 'grid-emax-0.6-mine--0.4-maxe-0.4.bin',
#                   hardGridFile = 'grid-emax-3-mine--30-maxe-30.bin',
#                   nativeGridFile = 'grid-emax-100-mine--100-maxe-100.bin',
#                   receptorPDB = './protein.pdb', receptorPSF = './protein.psf',
#                   ligPDB = './ligand.pdb', ligSeg = 'LIGA', confDir = './conformer/',
#                   placementDir = './placement/',  exhaustiveness = 'high',
#                   numPlace = 100, numCopy = 1000, flag_delete_conformer= True,
#                   flag_delete_placement = True, flag_save_all = True,
#                   flag_save_cluster= True, flag_save_top = True,
#                   flag_suppress_print = True, flag_center_ligand = True,
#                   flag_fast_grid = False, flag_use_hbond = False,
#                   flag_fast_placement = True, threshold = 2500,
#                   sort_energy = 'total_energy', saveDir = './dockresult/'):
