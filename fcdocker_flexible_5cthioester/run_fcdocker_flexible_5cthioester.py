"""Run FCDOCKER for a given ligand and receptor"""

import os
import argparse

# Backbone atoms for 5C thioester group
backbone_5cthioester = [
    "C6",
    "O6",
    "C3",
    "O5",
    "H25",
    "H18",
    "C2",
    "C5",
    "H24",
    "H22",
    "H23",
    "H21",
    "C4",
    "H19",
    "H20",
    "C1",
    "H17",
    "H16",
    "O4",
    "P1",
    "O1",
    "H14",
    "O2",
    "O3",
]


# region Parse arguments
parser = argparse.ArgumentParser(
    description="Run FCDOCKER for a given ligand and receptor"
)


# Helper function to get absolute path
def abspath(path):
    return os.path.abspath(path)


# Required arguments
parser.add_argument(
    "--ligPDB",
    type=abspath,
    help="Ligand PDB file",
)
parser.add_argument(
    "--ligandrtf",
    type=abspath,
    help="Ligand RTF file",
)
parser.add_argument(
    "--receptorPDB",
    type=abspath,
    help="Receptor PDB file",
)
parser.add_argument(
    "--receptorPSF",
    type=abspath,
    help="Receptor PSF file",
)
parser.add_argument(
    "--flexchain_csv",
    type=abspath,
    help="Flexible chain CSV file",
)
parser.add_argument(
    "--C5_rtf",
    type=abspath,
    help="5cthioester RTF file",
)
parser.add_argument(
    "--C5_segname",
    type=str,
    help="Segment name for 5cthioester",
    required=True,
)
parser.add_argument(
    "--workdir",
    type=abspath,
    help="Working directory",
)
parser.add_argument(
    "--num",
    type=int,
    help="Number of conformations to generate",
)
parser.add_argument(
    "--copy",
    type=int,
    help="Number of copies of each conformation",
)
parser.add_argument(
    "--xcen",
    type=float,
    help="X center of grid box",
)
parser.add_argument(
    "--ycen",
    type=float,
    help="Y center of grid box",
)
parser.add_argument(
    "--zcen",
    type=float,
    help="Z center of grid box",
)
parser.add_argument(
    "--maxlen",
    type=float,
    help="Maximum length of grid box",
)

parser.add_argument(
    "--topdir",
    type=abspath,
    default="/home/azamh/acyl/toppar",
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

# region Import pyCHARMM
os.environ["CHARMM_LIB_DIR"] = "/home/azamh/charmm/08302023/build_pyCHARMM/install/lib"
import pandas as pd
import pycharmm
import pycharmm.lib as lib
import pycharmm.read as read
import pycharmm.lingo as lingo
import pycharmm.settings as settings
from modified_cdocker import FCDOCKER_flexible_5cthioester

# endregion Import pyCHARMM

# region Topology and parameter files
settings.set_bomb_level(-1)
read.rtf(os.path.join(args.topdir, "top_all36_prot.rtf"))
read.rtf(os.path.join(args.topdir, "top_all36_cgenff.rtf"), append=True)
read.prm(os.path.join(args.topdir, "par_all36m_prot.prm"), flex=True)
read.prm(os.path.join(args.topdir, "par_all36_cgenff.prm"), append=True, flex=True)
settings.set_bomb_level(0)
lingo.charmm_script("stream " + args.ligandrtf)
lingo.charmm_script("stream " + args.C5_rtf)

# endregion Topology and parameter files

# Read in the flexible chain information
flexchain = pd.read_csv(args.flexchain_csv, sep=",", index_col=0)


# Add in 5cthioester group to flexchain
flexchain.loc[len(flexchain)] = [
    "1",
    args.C5_segname,
]
print(flexchain)

# Workdir
if not os.path.exists(args.workdir):
    os.makedirs(args.workdir)

# Change to workdir
os.chdir(args.workdir)

clusterResult, dockResult = FCDOCKER_flexible_5cthioester(
    C5_backbone=backbone_5cthioester,
    C5_segid=args.C5_segname,
    xcen=args.xcen,
    ycen=args.ycen,
    zcen=args.zcen,
    maxlen=args.maxlen,
    ligPDB=args.ligPDB,
    receptorPDB=args.receptorPDB,
    receptorPSF=args.receptorPSF,
    num=args.num,
    copy=args.copy,
    flexchain=flexchain,
    probeFile=os.path.join(args.topdir, "fftdock_c36prot_cgenff_probes.txt"),
)
