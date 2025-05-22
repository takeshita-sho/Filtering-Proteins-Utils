import os
from glob import glob
import pymol
import pandas as pd

def contacts(file):
    """
    Returns an array of residue IDs that are within 10 angstroms of organic molecules.
    
    Args:
        file (str): Path to the PDB file
        
    Returns:
        list: List of residue IDs within 10Å of organic molecules
    """
    # Initialize PyMOL and load structure
    pymol.cmd.reinitialize()
    pymol.cmd.load(file, "structure")
    pymol.cmd.remove("hydrogens")
    # Select organic molecules (not protein/water/ions)
    pymol.cmd.select("organics", "organic")

    
    # Select residues within 10Å of organic molecules
    pymol.cmd.select("nearby", f"byres polymer within 5 of organics")
    
    # Get the residue IDs
    myspace = {'resi_ids': []}
    pymol.cmd.iterate("nearby and name CA", "resi_ids.append(resi)", space=myspace)
    
    # Clean up
    pymol.cmd.reinitialize()
    
    # Return unique, sorted residue IDs
    return [int(i) for i in myspace['resi_ids']]


ref_path = "/home/stakesh1/scr4_jgray21/stakesh1/AF3_benchmark/AF3_predictions"
df = pd.DataFrame(columns=["PDB_ID", "FNAT"])

with open("/home/stakesh1/scr4_jgray21/stakesh1/benchmark_analysis/fnat_resi.txt",'r') as f:
    for line in f:
        
        ref, residues, offset = line.split(',')
        #For diffdock holo offset is 0
        #offset=0
        offset = int(offset)
        resi_ids = residues.strip().split('|')
        resi_ids = resi_ids
        resi_ids = [int(resi)-offset for resi in resi_ids]
        ref = ref[2:]
        pdbid = ref.replace("refined_xtal_pdb/","")[:4]

        struc_path = f"{ref_path}/{pdbid}"
        for i in range(1,6):
            for j in range(0,5):
                filename = f"{struc_path}/seed-{i}_sample-{j}/model.cif"
                pred_resi = contacts(filename)
                fnat = len(set(pred_resi) & set(resi_ids)) / len(set(resi_ids))
                df.loc[len(df)] = [pdbid, fnat]
print("Done")
df.to_csv("/home/stakesh1/scr4_jgray21/stakesh1/benchmark_analysis/results/AF3_fnat.csv", index=False)
