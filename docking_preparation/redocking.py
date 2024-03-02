#!/usr/bin/env python

import os
from vina import Vina
from pymol import cmd

def ligand_center(pdb_file): # 找到配体中心坐标
    cmd.reinitialize()
    cmd.load(pdb_file)
    center = cmd.centerofmass()
    xyz = [round(i,3) for i in center] # 坐标只保留3位小数
    return xyz

def redocking(pdb, receptor, ligand, pdb_ligand, outfile):
    v = Vina() #调用vina
    v.set_receptor(rigid_pdbqt_filename = receptor)
    v.set_ligand_from_file(ligand)
    v.compute_vina_maps(ligand_center(pdb_ligand), box_size=[30, 30, 30])
    v.dock(exhaustiveness=32, n_poses=1)
    v.write_poses(outfile, n_poses=1, overwrite= True)
    docking_score = v.energies(n_poses=1)[0][0]
    with open(outfile, "a") as f:
        f.write(pdb + "\t"+str(docking_score) + "\n")

def pdbqt2pdb(pdbqt_file, pdb_file): #用openbabel将对接结果pdbqt文件转为pdb文件
    os.system(f"obabel -ipdbqt {pdbqt_file} -opdb -O {pdb_file}")

def cal_rmsd(pdb_pose, redocking_pose, outfile):
    os.system(f"python isoRMSD.py -r {pdb_pose} -p {redocking_pose} -o {outfile}")


if __name__ == "__main__":
    pdblist = ["8D6D","8D6F","8D6C","8D6E","5VD0","5VD1","5VCW","5VCX","5VCZ","5VD3","8WJY","5VCV","5VCY"]
    for pdb in pdblist:
        pdbqtfile = os.path.join("pdb", pdb, pdb.lower() + "_redocking.pdbqt") 
        pdbfile = os.path.join("pdb", pdb, pdb.lower() + "_redocking.pdb")
        pdbqt2pdb(pdbqtfile, pdbfile)  