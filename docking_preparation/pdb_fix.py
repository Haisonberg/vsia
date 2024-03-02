from pdbfixer import PDBFixer
from openmm.app import PDBFile
from pymol import cmd
import os
from pathlib import Path

def pdb_down(pdblist): #使用pymol批量下载pdb文件
    for pdb in pdblist:
        if not os.path.exists(pdb):
            os.mkdir(pdb)
        if not "%s.pdb" %pdb in os.listdir(pdb):
            print("Downloading %s" %pdb)
            cmd.fetch(pdb, type="pdb", path=pdb)

def pdb_clean(pdb_file, out_file=None, ext="_clean.pdb"):
    """Extract all ATOM and TER records in a PDB file and write them to a new file.

    Args:
        pdb_file (str): Path of the PDB file from which ATOM and TER records
            will be extracted
        out_file (str): Optional argument to specify a particular output filename.
            Defaults to <pdb_file>.clean.pdb.
        ext (str): File extension to use for output file. Defaults to ".clean.pdb"
    """
    # find all ATOM and TER lines
    with open(pdb_file, "r") as fid:
        good = [l for l in fid if l.startswith(("ATOM", "TER"))]

    # default output file to <pdb_file>_clean.pdb
    if out_file is None:
        out_file = os.path.splitext(pdb_file)[0] + ext

    # write the selected records to a new file
    with open(out_file, "w") as fid:
        fid.writelines(good)
    return Path(out_file)

def add_residues(pdb_file, out_file):
    '''补全pdb文件中缺失的残基、重原子
       蛋白质加氢 
    '''
    fixer = PDBFixer(filename = pdb_file)
    print("Finding missing residues...")
    fixer.findMissingResidues()
    chains = list(fixer.topology.chains())
    keys = fixer.missingResidues.keys()
    # 开头和末端缺失的氨基酸残基不用补全
    for key in list(keys):
        chain = chains[key[0]]
        if key[1] == 0 or key[1] == len(list(chain.residues())):
            print("ok")
            del fixer.missingResidues[key]
    print("Finding missing atoms...")
    fixer.findMissingAtoms()
    print("Adding missing atoms...")
    fixer.addMissingAtoms()
    print("Adding missing hydrogens...")
    fixer.addMissingHydrogens(7.4)
    with open(out_file, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)

def remove_chains(pdb_file, out_file):
    fixer = PDBFixer(filename = pdb_file)
    chains = list(fixer.topology.chains())
    best_chain = chains[0]
    # 只保留氨基酸残基相对最完整的一条链
    for chain in chains:
        if len(list(chain.residues())) > len(list(best_chain.residues())):
            best_chain = chain
    for chain in chains:
        if chains.index(chain) != chains.index(best_chain):
            fixer.removeChains(chainIndices=[chains.index(chain)])
    with open(out_file, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)    

def pdb_fix(pdb_file, out_file): #应该叫add_atoms
    '''只补全缺失的重原子
       不补全缺失的残基，不对蛋白加氢 
    '''
    fixer = PDBFixer(filename = pdb_file)
    fixer.missingResidues = {} # do not want any residues at all to be added
    # 补全pdb文件中缺失的重原子
    print("Finding missing atoms...")
    fixer.findMissingAtoms()
    print("Adding missing atoms...")
    fixer.addMissingAtoms()
    with open(out_file, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)

def addHs(pdb_file, out_file): #只对蛋白加氢
    fixer = PDBFixer(filename = pdb_file)
    print("Adding missing hydrogens...")
    fixer.addMissingHydrogens(7.4)
    with open(out_file, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)

def removeHs(pdb_file, out_file): #使用pymol删除蛋白原有的氢原子
    cmd.reinitialize()
    cmd.load(pdb_file)
    cmd.remove("hydrogens")
    cmd.save(out_file)

def pdb2pdbqt(pdb_file, out_file): #用AutoDockTools生成pdbqt文件
    os.system(f"pythonsh /home/qihaixiang/opt/MGLTools-1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py \
              -r {pdb_file} -o {out_file}")
    
if __name__ == '__main__':
    pdblist = ["3P1A"]
    for pdb in pdblist:
        pdbfile = os.path.join("pdb", pdb, pdb.lower() + "_clean_1chain.pdb") #pdb文件绝对路径
        outfile = os.path.join("pdb", pdb, pdb.lower() + "_clean_1chain_allfixed.pdb") #输出文件路径
        add_residues(pdbfile, outfile)


