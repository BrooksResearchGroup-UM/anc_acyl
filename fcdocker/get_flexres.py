"""From Rigid CDOCKER results, define set of flexible residues in contact with ligand poses"""

import os
import argparse
import numpy as np
import biotite.structure as struc
import biotite.structure.io as strucio
import pandas as pd

# region Parse arguments
parser = argparse.ArgumentParser(
    description="Analzye Rigid CDOCKER results to define set of flexible residues in contact with ligand poses"
)


# Helper function to get absolute path
def abspath(path):
    return os.path.abspath(path)


parser.add_argument(
    "--receptor_pdb",
    type=abspath,
    help="Receptor PDB file",
    required=True,
)
parser.add_argument(
    "--rcdocker_results_dir",
    type=abspath,
    help="Directory with Rigid CDOCKER results",
    required=True,
)

# Optional arguments
parser.add_argument(
    "--contact_dist",
    type=float,
    help="Contact distance to define contacts",
    default=4.5,
)

args = parser.parse_args()

# Print arguments
print("Arguments:")
for arg in vars(args):
    print(arg, getattr(args, arg))

# endregion Parse arguments


# Helper function to get top ener poses
def get_top_ener(rcdocker_results_dir, n=10):
    dockresult_dir = os.path.join(rcdocker_results_dir, "dockresult")
    top_ener_dir = os.path.join(dockresult_dir, "top_ener")
    return [os.path.join(top_ener_dir, f"top_{i}.pdb") for i in range(1, n + 1)]


# Helper function to get top cluster poses
def get_top_cluster(rcdocker_results_dir, n=10):
    dockresult_dir = os.path.join(rcdocker_results_dir, "dockresult")
    top_cluster_dir = os.path.join(dockresult_dir, "cluster")
    return [os.path.join(top_cluster_dir, f"top_{i}.pdb") for i in range(1, n + 1)]


# Helper function to get cluster energy
def get_top_cluster_ener(rcdocker_results_dir):
    dockresult_dir = os.path.join(rcdocker_results_dir, "dockresult")
    top_cluster_ener_path = os.path.join(dockresult_dir, "clusterResult.tsv")
    df = pd.read_csv(top_cluster_ener_path, sep="\t")
    return df["total_energy"].tolist()


# Get top cluster poses
top_cluster = get_top_cluster(args.rcdocker_results_dir)

# Parse structures with biotite
receptor = strucio.load_structure(args.receptor_pdb)
poses = [strucio.load_structure(pose) for pose in top_cluster]

# Filter out hydrogens
receptor = receptor[receptor.element != "H"]
poses = [pose[pose.element != "H"] for pose in poses]

# Filter out 5cthioester
prosthetic = receptor[receptor.res_name == "5CT"]
receptor = receptor[receptor.res_name != "5CT"]

# Find contacts with each pose
resi_contacts = set()
for pose in poses:
    cell_list = struc.CellList(pose, args.contact_dist)
    contacts = cell_list.get_atoms(receptor.coord, radius=args.contact_dist)
    contact_indices = np.where((contacts != -1).any(axis=1))[0]
    contact_res_ids = receptor.res_id[contact_indices]
    resi_contacts.update(contact_res_ids)

print("Residues in contact with ligand poses:")
resi_contacts = sorted(resi_contacts)
print(resi_contacts)
print(f"Number of residues: {len(resi_contacts)}")
