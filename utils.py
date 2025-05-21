import numpy as np
import pandas as pd
import pyrosetta
import Bio
from Bio import BiopythonExperimentalWarning
from Bio.PDB import PDBParser
from Bio.PDB.PDBParser import PDBConstructionWarning
import warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    warnings.simplefilter('ignore', PDBConstructionWarning)
    
from pyrosetta import pose_from_pdb, init
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector
from pyrosetta.rosetta.core.pack.guidance_scoreterms.sap import calculate_sap
import os
import shutil
"""
Defines 2 hyperplanes and determines if the binder is between them
    resi1, resi2, resi3: Residues defining the first hyperplane
    resi4, resi5, resi6: Residues defining the second hyperplane
    Returns False if the binder is between the hyperplanes, True otherwise
"""
def tmcross(pdb,resi1,resi2,resi3,resi4,resi5,resi6,binder_chain,tm_chain,offset=0,):
#def tmcross(pdb,resi1=73,resi2=295,resi3=199,resi4=86,resi5=285,resi6=189,offset=52,chain_id='A'):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb)
    l = len(structure[0][binder_chain])
    
    e1 = structure[0][tm_chain][resi1+l-offset+1]["CA"].get_coord()
    e2 = structure[0][tm_chain][resi2+l-offset+1]["CA"].get_coord()
    e3 = structure[0][tm_chain][resi3+l-offset+1]["CA"].get_coord()
    i1 = structure[0][tm_chain][resi4+l-offset+1]["CA"].get_coord()
    i2 = structure[0][tm_chain][resi5+l-offset+1]["CA"].get_coord()
    i3 = structure[0][tm_chain][resi6+l-offset+1]["CA"].get_coord()
    
    enorm = np.cross(e2 - e1, e3 - e1)
    inorm = np.cross(i2 - i1, i3 - i1)
    for model in structure:
        for chain in model:
            if chain.id == binder_chain:
                for residue in chain:
                    for atom in residue:
                        i = atom.get_coord()
                        xe = i-e1
                        xi = i-i1
                        edot = np.dot(enorm,xe)
                        idot = np.dot(inorm,xi)
                        if (edot > 0 and idot < 0) or (edot < 0 and idot > 0):
                            
                            return False    
    return True

def sap(pdb):
    # Initialize PyRosetta
    init()

    # Create a pose from the PDB file
    pose = pose_from_pdb(pdb)

    # Define residue selectors for chain A
    chain_selector = ChainSelector("A")  # Select chain A
    score_sel = chain_selector
    sap_calculate_sel = chain_selector
    sasa_sel = chain_selector
   

    # Call calculate_sap with the required arguments
    sap_score = calculate_sap(
        pose,
        score_sel,
        sap_calculate_sel,
        sasa_sel
    )

    return sap_score
"""
read in csv and move files to destination
"""
def move_files(source, destination):
    df = pd.read_csv(source)
    for i in range(len(df)):
        file = df['file'][i]
        shutil.copy(file, destination)


#makes list of files in directory that do not appear in the file
def missing_files(dir, file, out):
    df = pd.read_csv(file, sep='\s+',usecols=range(1,11))
    f = open(out, 'w')
    
    for i in os.listdir(dir):
        hit = not df['description'].str.contains(i[:-4]).any()
        if hit:
            f.write(i[:-4]+'\n')
    f.close()

#calculate minimum distance between 2 chains
def calc_residue_dist(residue_one, residue_two) :
    """Returns the C-alpha distance between two residues"""
    diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
    return np.sqrt(np.sum(diff_vector * diff_vector))

def min_dist(pdb, dist=8.0) :
    structure = PDBParser(QUIET=True).get_structure('struc', pdb)
    model = structure[0]
    chain_one = model['A']
    chain_two = model['B']
    """Returns a matrix of C-alpha distances between two chains"""
    
    for row, residue_one in enumerate(chain_one) :
        for col, residue_two in enumerate(chain_two) :
            if calc_residue_dist(residue_one, residue_two) < dist:
                return True
            
    return False


def get_seq(pdb):
    # Just an example input pdb
    d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
            'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
            'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
            'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    # run parser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('struct', pdb)    
    chain = structure[0]['A']

    seq = []
    for residue in chain:
        seq.append(d3to1[residue.resname])
    return ''.join(seq)
