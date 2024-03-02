import os
from rdkit import Chem
from rdkit.Chem import AllChem

def pdb2sdf(pdb_file, out_file): #人工提取得到小分子pdb文件，openbabel转为sdf格式，错误需人工纠正
    os.system(f"obabel -ipdb {pdb_file} -osdf -O {out_file}")

def sdf2pdbqt(sdf_file, out_file): #meeko准备配体的pdbqt文件
    os.system(f"mk_prepare_ligand.py -i {sdf_file} -o {out_file}")

def add_fixed_atoms_helper(ff, mol): #索引重原子
  for i in range(0, mol.GetNumAtoms()):
    atom = mol.GetAtomWithIdx(i)
    if ((atom.GetAtomicNum() != 1) and (i < mol.GetNumAtoms())):
      ff.AddFixedPoint(i)

def H_minimize(sdf_file, out_file): #为配体加氢，并用MMFF94s力场做重原子约束的能量最小化
  mol = Chem.MolFromMolFile(sdf_file)
  mol_H = Chem.AddHs(mol, addCoords=True)
  # MMFF94s力场优化
  mp = AllChem.MMFFGetMoleculeProperties(mol_H, mmffVariant='MMFF94s')
  ff = AllChem.MMFFGetMoleculeForceField(mol_H, mp)
  add_fixed_atoms_helper(ff, mol_H)
  ff.Minimize(maxIts=1000)
  #输出sdf文件
  writer = Chem.SDWriter(out_file)
  writer.write(mol_H)

if __name__ == '__main__':
    pdblist = ["8D6D","8D6F","8D6C","8D6E","5VD0","5VD1","5VCW","5VCX","5VCZ","5VD3","8WJY","5VCV","5VCY"]
    for pdb in pdblist:
        sdffile = os.path.join("pdb", pdb, pdb.lower() + "_lig_H.sdf") #pdb文件绝对路径
        outfile = os.path.join("pdb", pdb, pdb.lower() + "_lig.pdbqt") #输出文件路径
        sdf2pdbqt(sdffile, outfile)


