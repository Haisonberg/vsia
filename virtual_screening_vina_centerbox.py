#!/usr/bin/env python

import os
from vina import Vina

dirs = os.listdir("/home/qihaixiang/pkmyt1_VS/benchmark/aNd_prepared")
for file in dirs:
    v = Vina()
    v.set_receptor(rigid_pdbqt_filename="/home/qihaixiang/pkmyt1_VS/8D6EB_H.pdbqt")
    input_molecule_file = os.path.join("/home/qihaixiang/pkmyt1_VS/benchmark/aNd_prepared",file)
    v.set_ligand_from_file(input_molecule_file)

    v.compute_vina_maps(center=[5.584, -12.326, 10.915], box_size=[30, 30, 30])
    v.dock(exhaustiveness=24,n_poses=1)
    v.write_poses(pdbqt_filename="/home/qihaixiang/pkmyt1_VS/benchmark/results_vina_centerbox/"+str(file)[:-6]+"_out"+".pdbqt",n_poses=1)
    docking_score = v.energies(n_poses=1)[0][0]
   
    with open("/home/qihaixiang/pkmyt1_VS/benchmark/output_vina_centerbox.txt", "a") as f:
        f.write(str(file)[:-6]+"\t"+str(docking_score)+"\n")
