#AF3 lrms calculation
import pymol
from pymol import cmd
import os
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS
import pandas as pd

#might have to make folder with all the pdbs before this
#could just load in pdb and sdf and merge into one object
def align(reference,struc,residues,offset,pdbid,save):
    resi_ids = residues.strip().split('|')
    resi_ids = resi_ids[:-1]
    ref_resis = ""
    struc_resis = ""
    
    for residue in resi_ids:
        struc_resi = int(residue)-offset
        ref_resis += f"{residue}+"
        struc_resis += f"{struc_resi}+"
    ref_resis = ref_resis[:-1]
    struc_resis = struc_resis[:-1]
    
    cmd.load(reference,"ref")
    cmd.remove("not alt ''+A")
    cmd.alter("all","alt=''")
    cmd.load(struc,"struc")
    cmd.select("ca_residues_ref", f"ref and resi {ref_resis} and name CA")
    #cmd.iterate("ca_residues_ref", "print(resn, resi, name)")
    cmd.select("ca_residues_struc", f"struc and resi {struc_resis} and name CA")
    #cmd.iterate("ca_residues_struc", "print(resn, resi, name)")
    cmd.align("ca_residues_struc","ca_residues_ref")

    cmd.save(f"{save}/{pdbid}_carb.pdb","struc and organic")
    cmd.reinitialize()
def get_mcs_rmsd(ref_file, mol_file):
    # Load molecules
    ref_mol = Chem.MolFromPDBFile(ref_file) #returns mol obj
    mol_mol = Chem.MolFromPDBFile(mol_file,sanitize=False)
    #print(ref_mol,mol_mol)
    #print(mol_file)
    mcs = rdFMCS.FindMCS([ref_mol, mol_mol])
    mcs_smarts = mcs.smartsString
    mcs_mol = Chem.MolFromSmarts(mcs_smarts)

    # Get the atom indices of the MCS in both molecules
    ref_match = ref_mol.GetSubstructMatch(mcs_mol)# returns a tuple of integers mol.GetAtoms()
    mol_match = mol_mol.GetSubstructMatch(mcs_mol)
    mmap = list(zip(mol_match, ref_match))
    


    # Calculate the RMSD for the MCS atoms
    return Chem.rdMolAlign.CalcRMS(mol_mol, ref_mol,map=[mmap])#, map=[list(mol_match),[list(ref_match)]]) 
    #takes symmetry into account! maybe try without atomIds


# To align stucs
ref_path = "/home/stakesh1/scr4_jgray21/stakesh1/ubcapin"
struc_path = "/home/stakesh1/scr4_jgray21/stakesh1/AF3_benchmark/AF3_predictions"
save_path = "/home/stakesh1/scr4_jgray21/stakesh1/benchmark_analysis/results/AF3_lrms.csv"
df = pd.DataFrame(columns=["PDB_ID", "LRMS"])

with open("/home/stakesh1/scr4_jgray21/stakesh1/benchmark_analysis/interacting_resi.txt",'r') as f:
    for line in f:
        
        ref, residues, offset = line.split(',')
        ref = ref[2:]
        ref_file = os.path.join(ref_path, ref)
        pdbid = ref.replace("refined_xtal_pdb/","")[:4]
        print(pdbid)
        for i in range(1,6):
            for j in range(5):

                struc1 = os.path.join(struc_path,f"{pdbid}/seed-{i}_sample-{j}/model.cif")
                carb_path = os.path.join(struc_path,f"{pdbid}/seed-{i}_sample-{j}")
                align(ref_file, struc1, residues,int(offset),pdbid,carb_path)
                #print(pdbid,carb_path)
                carb_file = f"{carb_path}/{pdbid}_carb.pdb"
                lrms = get_mcs_rmsd(f"/home/stakesh1/scr4_jgray21/stakesh1/ubcapin/xtal_carb/{pdbid}_xtalcarb.pdb", carb_file)
                df.loc[len(df)] = [f"{pdbid}-{i}-{j}", lrms]
                struc2 = f"/home/stakesh1/scr4_jgray21/stakesh1/ubcapin/xtal_carb/{pdbid}_xtalcarb_2.pdb"
                struc3 = f"/home/stakesh1/scr4_jgray21/stakesh1/ubcapin/xtal_carb/{pdbid}_xtalcarb_3.pdb"
                
                if os.path.isfile(struc2):
                    lrms2 = get_mcs_rmsd(struc2, carb_file)
                    df.loc[len(df)] = [f"{pdbid}-{i}-{j}_2", lrms2]
                if os.path.isfile(struc3):
                    lrms3 = get_mcs_rmsd(struc3, carb_file)
                    df.loc[len(df)] = [f"{pdbid}-{i}-{j}_3", lrms3]

df.to_csv(save_path, index=False)