#!/usr/bin/env python3
# extract certein chain from pdb file

import re
import argparse
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import three_to_one
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq

_hydrogen = re.compile("[123 ]*H.*")
class select_chain(object):
    def __init__(self, chain_list, model_id=0): 
        """Initialize the class.""" 
        self.chain_list = chain_list 
        self.model_id = model_id
    
    def accept_model(self, model): 
        """Verify if model match the model identifier.""" 
        # model - only keep model 0 
        if model.get_id() == self.model_id: 
            return 1 
        return 0

    def accept_chain(self, chain): 
        """Verify if chain match chain identifier.""" 
        if chain.get_id() in self.chain_list: 
            return 1
        return 0
    
    def accept_residue(self, residue):
        hetatm_flag, resseq, icode = residue.get_id() 
        if hetatm_flag != " ": 
        # skip HETATMS 
            return 0 
        if icode != " ": 
            print("WARNING: Icode %s at position %s" % (icode, resseq))
        return 1

    def accept_atom(self, atom):
        name = atom.get_id() 
        if _hydrogen.match(name): 
            return 0 
        else: 
            return 1

def extract_sequence(in_pdb):
    import os.path
    parser = PDBParser(PERMISSIVE=1)
    struct = parser.get_structure("input", in_pdb)
    seq_list = []
    for chain in struct.get_chains():
        chain_id = f"{os.path.basename(in_pdb).split('.')[0]}_{chain.id}"
        chain_sequence = ""
        for residue in chain.get_residues():
            chain_sequence += three_to_one(residue.get_resname())
        record = SeqRecord(Seq(chain_sequence), id=chain_id, name=chain_id, description=chain_id)
        #print(record)
        seq_list.append(record)
    SeqIO.write(seq_list, in_pdb.split('.')[0]+".fasta", "fasta")
    return 0

def renumber(struct):
    #loop 1: renumber residues to negative number to avoid errors
    chain_id = ""
    for residue in struct.get_residues():
        chain = residue.get_parent()
        if chain_id != chain.get_id():
            chain_id = chain.get_id()
            residue_id = -1
        #print(chain.get_id(), residue_id)
        residue.id=(' ',residue_id,' ')
        residue_id -= 1
    #loop 2
    chain_id = ""
    for residue in struct.get_residues():
        chain = residue.get_parent()
        if chain_id != chain.get_id():
            chain_id = chain.get_id()
            residue_id = 1
        #print(chain.get_id(), residue_id)
        residue.id=(' ',residue_id,' ')
        residue_id += 1
    return struct

def extract(in_pdb, chain_ids, out_filename): 
    #Write out selected portion to filename. 
    if out_filename == None:
        out_filename = f"{in_pdb.split('.')[0]}_{chain_ids}.pdb"
    chain_list = list(chain_ids)
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    structure = parser.get_structure("input", in_pdb)
    sele = select_chain(chain_list)
    structure = renumber(structure)
    io = PDBIO()
    io.set_structure(structure)
    io.save(out_filename, sele)
    return out_filename

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="extract chains from pdb file and renumber it")
    parser.add_argument("input_pdb", help="input pbd file that need dice")
    parser.add_argument("chain", help="chain names that need to be extracted")
    parser.add_argument("output_pdb", nargs="?", default=None, help="output pdb file name, default = <input file name>_<chain>.pdb")
    args = parser.parse_args()
    output_file = extract(args.input_pdb, args.chain, args.output_pdb)
    extract_sequence(output_file)
    print(output_file)
