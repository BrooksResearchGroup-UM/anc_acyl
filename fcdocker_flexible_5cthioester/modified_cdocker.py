"""Run FCDOCKER with flexible 5cthioester group"""

## Import module
import pycharmm
import pycharmm.ic as ic
import pycharmm.psf as psf
import pycharmm.coor as coor
import pycharmm.read as read
import pycharmm.grid as grid
import pycharmm.write as write
import pycharmm.lingo as lingo
import pycharmm.generate as gen
import pycharmm.energy as energy
import pycharmm.cons_fix as cons_fix
import pycharmm.minimize as minimize
import pycharmm.settings as settings
from pycharmm.implicit_solvent import FACTS
from pycharmm.cdocker import (
    FCDOCKER_init_place,
    FCDOCKER_fast_init_place,
    FCDOCKER_crossover,
    FCDOCKER_mutation,
    FCDOCKER_calc_dihedral,
    default_ommd_sian,
    cluster_mmtsb,
    top_N_cluster,
    protein_dihedral,
    _rm_,
    _mkdir_,
    _mv_,
    _cp_,
    _fcdocker_cluster_check_,
)

import numpy as np
import pandas as pd
from os import listdir, system
from math import pi, exp, atan2
from random import random, choices
from scipy.spatial.transform import Rotation as R


def FCDOCKER_flexible_5cthioester(
    C5_backbone,
    C5_segid,
    xcen=0,
    ycen=0,
    zcen=0,
    maxlen=10,
    num=20,
    copy=25,
    generation=2,
    threshold_init=2500,
    threshold_mutate=100,
    flag_grid=False,
    flag_form=False,
    flag_delete_grid=True,
    probeFile='"../Toppar/fftdock_c36prot_cgenff_probes.txt"',
    softGridFile="grid-emax-0.6-mine--0.4-maxe-0.4.bin",
    hardGridFile="grid-emax-3-mine--30-maxe-30.bin",
    nativeGridFile="grid-emax-100-mine--100-maxe-100.bin",
    receptorPDB="./protein.pdb",
    receptorPSF="./protein.psf",
    receptorCard='"./flexchain.crd"',
    saveLig="./ligand/",
    saveProt="./protein/",
    crossoverLig="./crossover_ligand/",
    crossoverProt="./crossover_protein/",
    saveLigFinal="./ligand_final/",
    saveProtFinal="./protein_final/",
    ligPDB="./ligand.pdb",
    ligSeg="LIGA",
    flexchain=None,
    placementDir="./placement/",
    flag_save_all=False,
    flag_save_cluster=True,
    flag_save_placement=False,
    flag_save_crossover=False,
    flag_suppress_print=True,
    flag_center_ligand=True,
    flag_fast_grid=False,
    flag_fast_placement=False,
    exhaustiveness="high",
    top_N_result=10,
    sort_energy="total_energy",
    saveDir="./dockresult/",
):
    """Flexible CDOCKER method for 5cthioester group

    Parameters
    ----------
    xcen: float
           center of the docking pocket
    ycen: float
           center of the docking pocket
    zcen: float
           center of the docking pocket
    maxlen: float
             size of the grid box
    num: float
          number of conformer
    copy: float
           number of copies for each conformer
    threshold_init: float
                     energy threshold for initial placement
    threshold_mutate: float
                       energy threshold for mutation
    flag_grid: bool
                whether or not grid need to be generated
    flag_form: bool
                whether or not grid form is formatted
    flag_delete_grid: bool
                       whether or not delete grid after calculation
    probeFile: str
                gpu probe file names
    softGridFile: str
                   soft grid file name
    hardGridFile: str
                   hard grid file name
    nativeGridFile: str
                     native grid file name
    receptorPDB: str
              protein pdb file name
    receptorPSF: str
              protein psf file name
    receptorCard: str
                   receptor coordinate card name
    saveLig: str
              ligand docked result saved folder name
    saveProt: str
               protein docked result saved folder name
    crossoverLig: str
                   crossover ligand folder before mutation
    crossoverProt: str
                    crossover protein folder before mutation
    saveLigFinal: str
                   ligand final docked result saved folder name
    saveProtFinal: str
                    protein final docked result saved folder name
    ligPDB: str
             ligand pdb file name
    ligSeg: str
             ligand segment ID
    flexchain: pd.DataFrame
                dataframe of the flexchain selection
    placementDir: str
                   placement folder name
    flag_save_all: bool
                    whether or not save all docked pose
    flag_save_cluster: bool
                        whether or not save clustered results
    flag_save_placement: bool
                          whether or not save inital placement after docking
    flag_save_crossover: bool
                          whether or not save crossover after docking
    flag_suppress_print: bool
                          whether or not suppress printing
    flag_center_ligand: bool
                         whether or not ligand need to be centered
    flag_fast_grid: bool
                     whether or not just use grid minimize result i.e., skip all atom mini
    flag_fast_placement : bool
                          whether or not use fast initial placement
    exhaustiveness : str
                     exhaustiveness for fast placement, high, medium, low
    top_N_result: int
                   number of top N clusters, final generation uses top N + 5 clusters
    sort_energy: str
                  sorting method
    saveDir: str
              folder name for docked result

    Returns
    -------
    clusterResult: pd.DataFrame
          clustering result
    dockResult : pd.DataFrame
          docking result
    """
    ## Center ligand if flag_center_ligand = True
    if flag_center_ligand:
        read.sequence_pdb(ligPDB)
        settings.set_bomb_level(-1)  ## for 3-mem ring
        gen.new_segment(seg_name=ligSeg)
        settings.set_bomb_level(0)
        read.pdb(ligPDB, resid=True)
        ligand = pycharmm.SelectAtoms().by_seg_id(ligSeg)

        xyz = coor.get_positions().to_numpy()
        gridCenter = np.array([xcen, ycen, zcen])
        ligCenter = (np.amin(xyz, axis=0) + np.amax(xyz, axis=0)) / 2
        xyz = xyz - ligCenter + gridCenter
        new_xyz = pd.DataFrame(xyz, columns=["x", "y", "z"])
        coor.set_positions(new_xyz)
        write.coor_pdb(
            ligPDB,
            selection=ligand,
            title="""Title
                       * Move ligand to the center of binding pocket""",
        )

    ## Generate grids for docking if no grid files generated before
    if not flag_grid:
        ## Read in protein and prepare flex side chain selection
        if psf.get_natom() > 0:
            psf.delete_atoms()
        read.psf_card(receptorPSF, append=True)
        read.pdb(receptorPDB, resid=True)

        for idx, row in flexchain.iterrows():
            resid = str(row["res_id"])
            segid = row["seg_id"]
            tmpI = pycharmm.SelectAtoms().by_res_id(res_id=resid)
            tmpJ = pycharmm.SelectAtoms().by_seg_id(seg_id=segid)
            if idx == 0:
                flexSideChain = tmpI & tmpJ
            else:
                flexSideChain = flexSideChain | (tmpI & tmpJ)

        settings.set_bomb_level(-1)
        psf.delete_atoms(flexSideChain)
        settings.set_bomb_level(0)

        # TODO: Figure out why you can't delete atoms and generate grid consecutively
        # Write psf and pdb
        write.psf_card("receptor_grid.psf")
        write.coor_pdb("receptor_grid.pdb")
        psf.delete_atoms()
        settings.set_bomb_level(0)
        # Read in pdb
        read.psf_card("receptor_grid.psf")
        read.pdb("receptor_grid.pdb", resid=True)
        ic.prm_fill(replace_all=False)
        ic.build()

        settings.set_bomb_level(0)

        ## Generate grids for docking
        genGrid = grid.Grid()
        gridSet = {
            "xCen": xcen,
            "yCen": ycen,
            "zCen": zcen,
            "xMax": maxlen,
            "yMax": maxlen,
            "zMax": maxlen,
            "emax": 0.6,
            "maxe": 0.4,
            "mine": -0.4,
            "gridFile": softGridFile,
            "flag_form": flag_form,
            "probeFile": probeFile,
        }
        genGrid.setVar(gridSet)
        status = genGrid.generate()
        print("Grid generation for " + softGridFile + " is ", status)

        gridSet = {"emax": 3, "maxe": 30, "mine": -30, "gridFile": hardGridFile}
        genGrid.setVar(gridSet)
        status = genGrid.generate()
        print("Grid generation for " + hardGridFile + " is ", status)

        gridSet = {"emax": 100, "maxe": 100, "mine": -100, "gridFile": nativeGridFile}
        genGrid.setVar(gridSet)
        status = genGrid.generate()
        print("Grid generation for " + nativeGridFile + " is ", status)

    ## Fast initial placement
    if flag_fast_placement:
        FCDOCKER_fast_init_place(
            receptorPDB=receptorPDB,
            receptorPSF=receptorPSF,
            ligPDB=ligPDB,
            ligSeg=ligSeg,
            placementDir=placementDir,
            exhaustiveness=exhaustiveness,
            num=num,
            copy=copy,
        )

    ## Prepare protein explicit atoms
    if psf.get_natom() > 0:
        psf.delete_atoms()
    read.psf_card(receptorPSF, append=True)
    read.pdb(receptorPDB, resid=True)

    for idx, row in flexchain.iterrows():
        resid = str(row["res_id"])
        segid = row["seg_id"]
        tmpI = pycharmm.SelectAtoms().by_res_id(res_id=resid)
        tmpJ = pycharmm.SelectAtoms().by_seg_id(seg_id=segid)
        if idx == 0:
            flexSideChain = tmpI & tmpJ
        else:
            flexSideChain = flexSideChain | (tmpI & tmpJ)
    settings.set_bomb_level(-1)
    psf.delete_atoms(flexSideChain.__invert__())
    settings.set_bomb_level(0)

    ## Read in ligand
    read.sequence_pdb(ligPDB)
    settings.set_bomb_level(-1)  ## for 3-mem ring
    gen.new_segment(seg_name=ligSeg)
    settings.set_bomb_level(0)
    read.pdb(ligPDB, resid=True)
    ligand = pycharmm.SelectAtoms().by_seg_id(ligSeg)

    ## Declare receptor(protein) backbone
    backboneAtom = [
        "N",
        "CA",
        "C",
        "O",
        "HA",
        "HN",
        "HA1",
        "HA2",
        "HT1",
        "HT2",
        "HT3",
        "OT1",
        "OT2",
    ]
    backboneAtom = np.asarray(backboneAtom)
    for idx in np.arange(len(backboneAtom)):
        atom = backboneAtom[idx]
        if idx == 0:
            backbone = pycharmm.SelectAtoms().by_atom_type(atom)
        backbone = backbone | pycharmm.SelectAtoms().by_atom_type(atom)

    print(f"Number of protein backbone atoms: {backbone.get_n_selected()}")

    # Select prosthetic group backbone
    prosthetic_atoms = pycharmm.SelectAtoms().by_seg_id(C5_segid)
    print(f"Number of prosthetic group atoms: {prosthetic_atoms.get_n_selected()}")
    for idx, atom in enumerate(C5_backbone):
        if idx == 0:
            C5_backbone = pycharmm.SelectAtoms().by_atom_type(atom) & prosthetic_atoms
        else:
            C5_backbone = C5_backbone | (
                pycharmm.SelectAtoms().by_atom_type(atom) & prosthetic_atoms
            )

    print(f"Number of prosthetic group backbone atoms: {C5_backbone.get_n_selected()}")

    ## Add to backbone
    backbone = backbone | C5_backbone

    print(f"Total number of backbone atoms: {backbone.get_n_selected()}")

    ## Define fix part and flex part of the system
    fixsc = ligand.__invert__() & backbone
    flexsc = fixsc.__invert__()
    receptor = ligand.__invert__()
    cons_fix.setup(fixsc)
    write.coor_card(
        receptorCard,
        selection=receptor,
        title="""Title
                    * Receptor Side Chain Initial Position""",
    )

    ## Update nonbond interactions
    update_nonbond = pycharmm.UpdateNonBondedScript(
        atom=True,
        switch=True,
        vswitch=True,
        vdwe=True,
        elec=True,
        rdie=True,
        cutnb=12,
        ctofnb=10,
        ctonnb=8,
        emax=10000,
        maxe=10000,
        mine=-10000,
        epsilon=3,
    ).run()

    ## Coor stats
    allCoor = coor.get_positions()
    allCoor["atomIndex"] = np.arange(np.shape(list(ligand))[0])
    ligand_xyz = allCoor.loc[np.asarray(list(ligand)), ["x", "y", "z"]]
    receptor_xyz = allCoor.loc[np.asarray(list(receptor)), ["x", "y", "z", "atomIndex"]]

    ligStat = coor.stat(selection=ligand)
    xcen = ligStat["xave"]
    ycen = ligStat["yave"]
    zcen = ligStat["zave"]

    ## Set up the OpenMM system
    ommd = grid.OMMD()
    ommdSet = {
        "softGridFile": softGridFile,
        "hardGridFile": hardGridFile,
        "fix_select": fixsc,
        "flex_select": flexsc,
        "flag_form": flag_form,
        "numCopy": num * copy,
    }
    ommd.setVar(ommdSet)
    ommd.create()
    print("OMMD set up for flexible docking")
    print("Flexible docking for the first generation")

    ## OpenMM docking for initial generation
    if flag_suppress_print:
        settings.set_verbosity(1)
        settings.set_warn_level(1)
        settings.set_bomb_level(-1)
    if not flag_fast_placement:
        FCDOCKER_init_place(
            receptorCard=receptorCard,
            receptorSel=receptor,
            flexSel=flexsc,
            ligSel=ligand,
            ligPDB=ligPDB,
            ligSeg=ligSeg,
            hardGridFile=hardGridFile,
            nativeGridFile=nativeGridFile,
            placementDir=placementDir,
            num=num,
            copy=copy,
            threshold=threshold_init,
            flag_form=flag_form,
        )

    idxCopy = 1
    while idxCopy <= num * copy:
        read.pdb(placementDir + str(idxCopy) + ".pdb", resid=True)
        ommd.set_coor(idxCopy=idxCopy)
        print("Set OMMD ligand copy " + str(idxCopy))
        idxCopy += 1

    default_ommd_sian(ommd)
    _rm_("[0-9]*")

    hardSet = {"selection": flexsc, "flag_form": flag_form, "gridFile": hardGridFile}
    hardGrid = grid.CDOCKER()
    hardGrid.setVar(hardSet)
    hardGrid.read()
    energy.show()
    idxCopy = 1
    _mkdir_(saveLig + " " + saveProt)
    dockEner = []
    while idxCopy <= num * copy:
        ommd.copy_coor(idxCopy=idxCopy)
        print("Grid minimize OMMD ligand copy " + str(idxCopy))
        minimize.run_sd(nstep=50)
        minimize.run_abnr(nstep=1000, tolenr=1e-3)
        energy.show()
        totalener = energy.get_total()
        write.coor_pdb(str(idxCopy), selection=ligand)
        write.coor_pdb(
            saveLig + str(idxCopy) + ".pdb",
            selection=ligand,
            title="""Title
                       * The docked pose ID is %d
                       * The total energy is %8.6f """
            % (idxCopy, totalener),
        )
        write.coor_pdb(
            saveProt + str(idxCopy) + ".pdb",
            selection=receptor,
            title="""Title
                       * The docked pose ID is %d
                       * The total energy is %8.6f """
            % (idxCopy, totalener),
        )
        dockEner.append(totalener)
        idxCopy += 1

    hardGrid.off()
    hardGrid.clear()
    dock_result = np.zeros((num * copy, 2))
    dock_result[:, 0] = np.asarray(dockEner)
    dock_result[:, 1] = np.arange(num * copy) + 1
    settings.set_verbosity(5)
    settings.set_warn_level(5)
    settings.set_bomb_level(0)

    ## Loop through following generations
    idxGen = 2
    radius = 1
    n_cluster = 0
    while n_cluster <= 0:
        cluster_mmtsb(radius=radius, name="[0-9]*")
        n_cluster = _fcdocker_cluster_check_(total=num * copy)
        radius += 0.5
    cluster_result = top_N_cluster(N=top_N_result, total=num * copy)

    while idxGen <= generation:
        print("Flexible docking for generation: " + str(idxGen))
        if flag_suppress_print:
            settings.set_verbosity(1)
            settings.set_warn_level(1)
            settings.set_bomb_level(-1)

        ## Crossover
        final_result, pair_list = FCDOCKER_crossover(
            cluster_result, dock_result, num=num, copy=copy
        )
        _mkdir_(crossoverLig + " " + crossoverProt)
        for idx in np.arange(num * copy) + 1:
            ## Copy ligand
            ligandPose = saveLig + str(pair_list[idx - 1, 0]) + ".pdb"
            placement = crossoverLig + str(idx) + ".pdb"
            _cp_(ligandPose, placement)

            ## Copy protein
            proteinPose = saveProt + str(pair_list[idx - 1, 1]) + ".pdb"
            placement = crossoverProt + str(idx) + ".pdb"
            _cp_(proteinPose, placement)

        ## OpenMM docking for second generation
        FCDOCKER_mutation(
            final_result,
            pair_list,
            receptorSel=receptor,
            flexSel=flexsc,
            ligSel=ligand,
            ligPDB=ligPDB,
            ligSeg=ligSeg,
            saveLig=saveLig,
            saveProt=saveProt,
            crossoverLig=crossoverLig,
            crossoverProt=crossoverProt,
            hardGridFile=hardGridFile,
            nativeGridFile=nativeGridFile,
            placementDir=placementDir,
            num=num,
            copy=copy,
            threshold=threshold_mutate,
            flag_form=flag_form,
        )

        idxCopy = 1
        while idxCopy <= num * copy:
            read.pdb(placementDir + str(idxCopy) + ".pdb", resid=True)
            ommd.set_coor(idxCopy=idxCopy)
            print("Set OMMD ligand copy " + str(idxCopy))
            idxCopy += 1
        default_ommd_sian(ommd)

        _rm_("[0-9]* cluster.log cluster_list")

        hardSet = {
            "selection": flexsc,
            "flag_form": flag_form,
            "gridFile": hardGridFile,
        }
        hardGrid = grid.CDOCKER()
        hardGrid.setVar(hardSet)
        hardGrid.read()
        energy.show()
        idxCopy = 1
        _mkdir_(saveLig + " " + saveProt)
        dockEner = []
        while idxCopy <= num * copy:
            ommd.copy_coor(idxCopy=idxCopy)
            print("Grid minimize OMMD ligand copy " + str(idxCopy))
            minimize.run_sd(nstep=50)
            minimize.run_abnr(nstep=1000, tolenr=1e-3)
            energy.show()
            totalener = energy.get_total()
            write.coor_pdb(str(idxCopy), selection=ligand)
            write.coor_pdb(
                saveLig + str(idxCopy) + ".pdb",
                selection=ligand,
                title="""Title
                           * The docked pose ID is %d
                           * The total energy is %8.6f """
                % (idxCopy, totalener),
            )
            write.coor_pdb(
                saveProt + str(idxCopy) + ".pdb",
                selection=receptor,
                title="""Title
                           * The docked pose ID is %d
                           * The total energy is %8.6f """
                % (idxCopy, totalener),
            )
            dockEner.append(totalener)
            idxCopy += 1

        hardGrid.off()
        hardGrid.clear()
        dock_result = np.zeros((num * copy, 2))
        dock_result[:, 0] = np.asarray(dockEner)
        dock_result[:, 1] = np.arange(num * copy) + 1
        settings.set_verbosity(5)
        settings.set_warn_level(5)
        settings.set_bomb_level(0)

        ## Cluster docked poses
        radius = 1
        n_cluster = 0
        while n_cluster <= 0:
            cluster_mmtsb(radius=radius, name="[0-9]*")
            n_cluster = _fcdocker_cluster_check_(total=num * copy)
            radius += 0.5
        if idxGen < generation:
            cluster_result = top_N_cluster(N=top_N_result, total=num * copy)
        else:
            cluster_result = top_N_cluster(N=top_N_result + 5, total=num * copy)

        idxGen += 1

    ## Clear OMMD system
    ommd.clear()

    ## Analyze clustering result after final generation is done
    settings.set_verbosity(1)
    settings.set_warn_level(1)
    settings.set_bomb_level(-1)

    # Create copy of flexchain without prosthetic group
    flexchain_no_prosthetic = flexchain.copy(deep=True)
    flexchain_no_prosthetic = flexchain_no_prosthetic[
        flexchain_no_prosthetic["seg_id"] != C5_segid
    ]
    assert len(flexchain) == len(flexchain_no_prosthetic) + 1

    entropy, clusterID, sc_entropy = FCDOCKER_calc_dihedral(
        flexchain_no_prosthetic, cluster_result, saveProt, protein_dihedral
    )
    settings.set_verbosity(5)
    settings.set_warn_level(5)
    settings.set_bomb_level(0)
    cluster_size = []
    for cluster in np.unique(cluster_result[:, 1]):
        cluster_size.append(len(cluster_result[cluster_result[:, 1] == cluster]))

    ## Get list of pdb needed for all atom explicit minimization
    if flag_fast_grid:
        pdbID = np.ones(len(clusterID))
        for cluster in clusterID:
            tmpEner = []
            pdbIdx = cluster_result[cluster_result[:, 1] == cluster][:, 0]
            for pdb in pdbIdx:
                tmpEner.append(dock_result[dock_result[:, 1] == pdb][:, 0])
            cResult = pd.DataFrame()
            cResult["PDB"] = np.asarray(pdbIdx)
            cResult["Energy"] = np.asarray(tmpEner)
            cResult = cResult.sort_values(by=["Energy"], ignore_index=True)
            pdbID[cluster - 1] = cResult.iloc[0]["PDB"]
        pdbID = pdbID.astype(int)
    else:
        tmp_entropy = []
        for cluster in clusterID:
            tmp = np.ones(cluster_size[cluster - 1]) * entropy[cluster - 1]
            tmp_entropy = tmp_entropy + tmp.tolist()
        pdbID = cluster_result[:, 0]
        clusterID = cluster_result[:, 1]
        pdbID = pdbID.astype(int)
        entropy = np.asarray(tmp_entropy)

    ## Read in the explicit protein and ligand
    cons_fix.turn_off()
    psf.delete_atoms()
    read.psf_card(receptorPSF, append=True)
    read.pdb(receptorPDB, resid=True)
    read.sequence_pdb(ligPDB)
    settings.set_bomb_level(-1)  ## for 3-mem ring
    gen.new_segment(seg_name=ligSeg)
    settings.set_bomb_level(0)
    read.pdb(ligPDB, resid=True)
    for idx, row in flexchain.iterrows():
        resid = str(row["res_id"])
        segid = row["seg_id"]
        tmpI = pycharmm.SelectAtoms().by_res_id(res_id=resid)
        tmpJ = pycharmm.SelectAtoms().by_seg_id(seg_id=segid)
        if idx == 0:
            flexSideChain = tmpI & tmpJ
        else:
            flexSideChain = flexSideChain | (tmpI & tmpJ)
    ligand = pycharmm.SelectAtoms().by_seg_id(ligSeg)
    flexPart = flexSideChain | ligand
    fixPart = flexPart.__invert__()

    ## Minimize docked pose with explicit protein atoms
    vdw = []
    elec = []
    totalEner = []
    _mkdir_("tmpligand/ tmpprot/")
    cons_fix.setup(fixPart)
    if flag_suppress_print:
        settings.set_verbosity(1)
        settings.set_warn_level(1)
        settings.set_bomb_level(-1)

    for idx in np.arange(len(pdbID)):
        read.pdb(saveLig + str(pdbID[idx]) + ".pdb", resid=True)
        read.pdb(saveProt + str(pdbID[idx]) + ".pdb", resid=True)
        minimize.run_sd(nstep=50)
        minimize.run_abnr(nstep=1000, tolenr=0.001)
        tmptotal = energy.get_total() + entropy[idx]
        tmpelec = energy.get_elec()
        tmpvdw = energy.get_vdw()
        totalEner.append(tmptotal)
        elec.append(tmpelec)
        vdw.append(tmpvdw)

        write.coor_pdb(
            "tmpligand/" + str(pdbID[idx]) + ".pdb",
            selection=ligand,
            title="""Title
                       * The docked pose ID is %d
                       * The cluster ID is %d
                       * The total energy is %8.6f
                       * The vdw energy is %8.6f
                       * The elec energy is %8.6f
                       * The entropy is %8.6f """
            % (pdbID[idx], clusterID[idx], tmptotal, tmpvdw, tmpelec, entropy[idx]),
        )

        write.coor_pdb(
            "tmpprot/" + str(pdbID[idx]) + ".pdb",
            selection=flexSideChain,
            title="""Title
                       * The docked pose ID is %d
                       * The cluster ID is %d
                       * The total energy is %8.6f
                       * The vdw energy is %8.6f
                       * The elec energy is %8.6f
                       * The entropy is %8.6f """
            % (pdbID[idx], clusterID[idx], tmptotal, tmpvdw, tmpelec, entropy[idx]),
        )

    settings.set_verbosity(5)
    settings.set_warn_level(5)
    settings.set_bomb_level(0)

    ## Save docking result
    clusterResult = pd.DataFrame()
    if flag_fast_grid:
        clusterResult["total_energy"] = totalEner
        clusterResult["enthalpy"] = np.asarray(totalEner) - entropy
        clusterResult["vdw"] = vdw
        clusterResult["elec"] = elec
        clusterResult["entropy"] = entropy
        clusterResult["cluster_size"] = cluster_size
        clusterResult["PDB_name"] = pdbID
        clusterResult["cluster_id"] = clusterID
    else:
        explicitResult = pd.DataFrame()
        explicitResult["total_energy"] = totalEner
        explicitResult["enthalpy"] = np.asarray(totalEner) - entropy
        explicitResult["vdw"] = vdw
        explicitResult["elec"] = elec
        explicitResult["entropy"] = entropy
        explicitResult["PDB_name"] = pdbID
        explicitResult["cluster_id"] = cluster_result[:, 1]

        vdw = []
        elec = []
        pdbID = []
        entropy = []
        totalEner = []
        clusterID = np.unique(clusterID)
        for cluster in clusterID:
            tmp_result = explicitResult.loc[explicitResult["cluster_id"] == cluster]
            tmp_result = tmp_result.sort_values(by=["enthalpy"], ignore_index=True)
            pdbID.append(tmp_result.iloc[0]["PDB_name"])
            tmp_result = tmp_result.to_numpy()

            ## Get probability and ensemble average
            enthalpy = tmp_result[:, 1]
            enthalpy = enthalpy - np.amin(enthalpy)
            prob = -(1 / 0.593) * enthalpy
            prob = np.exp(prob)
            prob = prob / np.sum(prob)
            tmptotal = np.sum(tmp_result[:, 1] * prob)
            tmpvdw = np.sum(tmp_result[:, 2] * prob)
            tmpelec = np.sum(tmp_result[:, 3] * prob)
            tmpentropy = np.sum(tmp_result[:, 4] * prob)
            tmpcluster = np.sum(tmp_result[:, 6] * prob)

            vdw.append(tmpvdw)
            elec.append(tmpelec)
            entropy.append(tmpentropy)
            totalEner.append(tmptotal + tmpentropy)

        pdbID = np.asarray(pdbID)
        pdbID = pdbID.astype(int)
        clusterResult["total_energy"] = totalEner
        clusterResult["enthalpy"] = np.asarray(totalEner) - np.asarray(entropy)
        clusterResult["vdw"] = vdw
        clusterResult["elec"] = elec
        clusterResult["entropy"] = entropy
        clusterResult["cluster_size"] = cluster_size
        clusterResult["PDB_name"] = pdbID
        clusterResult["cluster_id"] = clusterID

    dockResult = pd.DataFrame()
    dockResult["enthalpy"] = dock_result[:, 0]
    dockResult["PDB_name"] = dock_result[:, 1].astype(int)
    dockResult = dockResult.sort_values(by=["enthalpy"], ignore_index=True)

    ## Sort dataframe and save results
    clusterResult = clusterResult.sort_values(by=[sort_energy], ignore_index=True)
    pdbList = clusterResult["PDB_name"].tolist()
    idx = 1
    _mkdir_(saveLigFinal + " " + saveProtFinal)
    for pdb in pdbList:
        ## Copy Ligand
        source = "tmpligand/" + str(pdb) + ".pdb"
        target = saveLigFinal + "top_" + str(idx) + ".pdb"
        _cp_(source, target)

        ## Copy protein
        source = "tmpprot/" + str(pdb) + ".pdb"
        target = saveProtFinal + "top_" + str(idx) + ".pdb"
        _cp_(source, target)

        idx += 1

    ## Clean and save results
    _mkdir_(saveDir)
    _mv_("cluster.log", saveDir)
    clusterResult.to_csv(saveDir + "clusterResult.tsv", sep="	")
    dockResult.to_csv(saveDir + "dockResult.tsv", sep="	")

    if not flag_fast_grid:
        explicitResult = explicitResult.sort_values(by=[sort_energy], ignore_index=True)
        explicitResult.to_csv(saveDir + "explicitResult.tsv", sep="	")
    if flag_save_placement:
        _mv_(placementDir, saveDir)
    if flag_save_all:
        _mkdir_(saveDir + "dock_pose/")
        _mv_(saveLig + " " + saveProt, saveDir + "dock_pose/")
    if flag_save_cluster:
        _mkdir_(saveDir + "cluster/")
        _mv_(saveLigFinal, saveDir + "cluster/ligand/")
        _mv_(saveProtFinal, saveDir + "cluster/protein/")
    if flag_save_crossover and generation > 1:
        _mkdir_(saveDir + "crossover/")
        _mv_(crossoverLig, saveDir + "crossover/")
        _mv_(crossoverProt, saveDir + "crossover/")
    if flag_delete_grid:
        _rm_(nativeGridFile + " " + softGridFile + " " + hardGridFile)

    _rm_(
        saveLig
        + " "
        + saveProt
        + " "
        + saveLigFinal
        + " "
        + saveProtFinal
        + " "
        + crossoverLig
        + " "
        + crossoverProt
    )
    _rm_("cluster* [0-9]* tmp* *crd ligand_rotamer.pdb " + placementDir)

    return clusterResult, dockResult
