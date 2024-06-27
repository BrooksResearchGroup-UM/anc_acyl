"""Run molecular dynamics with OpenMM"""

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
    "--c5_str",
    type=abspath,
    help="5cthioester RTF file",
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
    required=True,
)
parser.add_argument(
    "--ion_type",
    type=str,
    help="Ion type",
    required=True,
)

parser.add_argument(
    "--nequil",
    type=int,
    required=True,
    help="Number of equilibration steps",
)

parser.add_argument(
    "--nprod",
    type=int,
    required=True,
    help="Number of production steps",
)

parser.add_argument(
    "--nsavc",
    type=int,
    required=True,
    help="Number of steps between saving coordinates",
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


def is_factor(n):
    """Ensure that FFT grid is product of small primes 2, 3, 5"""
    if n % 2 != 0:
        return False  # favors even number
    while n:
        flag = False
        for x in (2, 3, 5):
            if n % x == 0:
                n = n / x
                flag = True
                break

        if flag:
            continue
        break

    if n == 1:
        return True
    return False


def checkfft(n, margin=5):
    """Check FFT grid size"""
    n = int(n) + margin
    while 1:
        if is_factor(n):
            break
        else:
            n += 1
    return n


# Convert time steps to nanoseconds and print
# 2 fs time step
dt = 2e-6  # ns
print(f"Equilibration steps: {args.nequil}")
print(f"Production steps: {args.nprod}")
print(f"Equilibration time: {args.nequil * dt} ns")
print(f"Production time: {args.nprod * dt} ns")

# Read in topology and parameter files
settings.set_bomb_level(-1)
read.rtf(os.path.join(args.toppar_dir, "top_all36_prot.rtf"))
read.rtf(os.path.join(args.toppar_dir, "top_all36_cgenff.rtf"), append=True)
read.prm(os.path.join(args.toppar_dir, "par_all36m_prot.prm"), flex=True)
read.prm(os.path.join(args.toppar_dir, "par_all36_cgenff.prm"), flex=True, append=True)
charmm_script(f"stream {args.c5_str}")
charmm_script(f"stream {args.ligand_str}")
charmm_script(f'stream {os.path.join(args.toppar_dir, "toppar_water_ions.str")}')
settings.set_bomb_level(0)


# Read in input structure
read.psf_card(args.input_psf)
read.pdb(args.input_pdb, resid=True)


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

# Setup periodic boundary conditions
print(f"Setting up periodic boundary conditions: boxhalf: {boxhalf}")
setup_pbc(boxhalf, solvated=True, ion_type=args.ion_type)

# determine fft grid size
fft = checkfft(n=np.ceil(boxhalf) * 2, margin=0)
print(f"FFT grid size: {fft}")

# Setup nonbonds with particle mesh ewald
nb_wPME = pycharmm.NonBondedScript(
    cutnb=12,
    cutim=12,
    ctonnb=8,
    ctofnb=10,
    eps=1.0,
    cdie=True,
    atom=True,
    vatom=True,
    switch=True,
    vfswitch=False,
    vswitch=False,
    inbfrq=-1,
    imgfrq=-1,
    ewald=True,
    pmewald=True,
    kappa=0.32,
    fftx=fft,
    ffty=fft,
    fftz=fft,
    order=4,
)
nb_wPME.run()

# Turn on shake
shake.on(bonh=True, fast=True, tol=1e-7)

# Turn on openmm
charmm_script("omm on platform cuda deviceid 0")
energy.show()

imgfrq = -1
dyn.set_fbetas(np.full((psf.get_natom()), 1.0, dtype=float))

# Paths to outputs
restart_file_path = os.path.join(args.workdir, "dyn.res")
lambda_file_path = os.path.join(args.workdir, "dyn.lam")
dcd_file_path = os.path.join(args.workdir, "dyn.dcd")
pdb_file_path = os.path.join(args.workdir, "dyn.pdb")

# CHARMM files
res_file = pycharmm.CharmmFile(
    file_name=restart_file_path,
    file_unit=2,
    formatted=True,
    read_only=False,
)
lam_file = pycharmm.CharmmFile(
    file_name=lambda_file_path,
    file_unit=3,
    formatted=False,
    read_only=False,
)

# Equilibrate
equil_script = pycharmm.DynamicsScript(
    leap=False,
    lang=True,
    start=True,
    nstep=args.nequil,
    timest=0.002,
    firstt=298.0,
    finalt=298.0,
    tbath=298.0,
    tstruc=298.0,
    teminc=0.0,
    twindh=0.0,
    twindl=0.0,
    iunwri=res_file.file_unit,
    iunlam=lam_file.file_unit,
    inbfrq=-1,
    imgfrq=imgfrq,
    iasors=0,
    iasvel=1,
    ichecw=0,
    iscale=0,
    iscvel=0,
    echeck=-1.0,
    nsavc=0,
    nsavv=0,
    nsavl=0,
    ntrfrq=0,
    isvfrq=args.nsavc,
    iprfrq=2 * args.nsavc,
    nprint=args.nsavc,
    ihtfrq=0,
    ieqfrq=0,
    ilbfrq=0,
    ihbfrq=0,
    blade=False,
    omm="gamma 2 prmc pref 1 iprsfrq 100",
)
equil_script.run()
res_file.close()
lam_file.close()

# Run production
res_file = pycharmm.CharmmFile(
    file_name=restart_file_path,
    file_unit=2,
    formatted=True,
    read_only=False,
)
lam_file = pycharmm.CharmmFile(
    file_name=lambda_file_path,
    file_unit=3,
    formatted=False,
    read_only=False,
)
dcd_file = pycharmm.CharmmFile(
    file_name=dcd_file_path,
    file_unit=1,
    formatted=False,
    read_only=False,
)

pycharmm.DynamicsScript(
    leap=False,
    lang=True,
    restart=True,
    nstep=args.nprod,
    timest=0.002,
    firstt=298.0,
    finalt=298.0,
    tbath=298.0,
    tstruc=298.0,
    teminc=0.0,
    twindh=0.0,
    twindl=0.0,
    iunwri=res_file.file_unit,  # restart file to write to
    iunrea=res_file.file_unit,  # this is the restart file to read from
    iunlam=lam_file.file_unit,  # this is lambda file to write to (for lambda dynamics?)
    inbfrq=-1,
    imgfrq=imgfrq,
    iuncrd=dcd_file.file_unit,  # this is the trajectory file to write to
    iasors=0,
    iasvel=1,
    ichecw=0,
    iscale=0,
    iscvel=0,
    echeck=-1.0,
    nsavc=args.nsavc,
    nsavv=0,
    nsavl=0,
    ntrfrq=0,
    isvfrq=args.nsavc,
    iprfrq=5 * args.nsavc,
    nprint=args.nsavc,
    ihtfrq=0,
    ieqfrq=0,
    ilbfrq=0,
    ihbfrq=0,
    blade=False,
    omm="gamma 2 prmc pref 1 iprsfrq 100",
).run()
res_file.close()
lam_file.close()
dcd_file.close()
write.coor_pdb(pdb_file_path)
