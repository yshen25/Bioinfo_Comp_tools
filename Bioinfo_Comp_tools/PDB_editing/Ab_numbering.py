# renumber pdb structures

from Bio import PDB
import os.path
import anarci
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import index_to_one, three_to_index

class select_chain(object):
    def __init__(self, chain_list=None, model_id=0): 
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
        if self.chain_list is None:
            return 1
        else:
            """the following code is used to verify if chain match chain identifier.""" 
            if chain.get_id() in self.chain_list: 
                return 1
            return 0
    
    def accept_residue(self, residue):
        hetatm_flag, resseq, icode = residue.get_id() 
        if hetatm_flag != " ": 
        # skip HETATMS 
            return 0
        if resseq < 0:
        # skip residues that are not recognized by ANARCI
            return 0
        if icode != " ": 
            print("WARNING: Icode %s at position %s" % (icode, resseq))
        return 1

    def accept_atom(self, atom):
        return 1
        # the following code removes hydrogens
        # name = atom.get_id() 
        # if _hydrogen.match(name): 
        #     return 0 
        # else: 
        #     return 1

class renumber_mAb:
    def __init__(self, scheme=None, CDR_indexes=None, flat=False, rename_chain=False):
        if not scheme:
            self.scheme = "imgt"
        else:
            self.scheme = scheme
        
        if not CDR_indexes:
            self.CDR_indexes = [27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117]
        else:
            self.CDR_indexes = CDR_indexe
        
        self.flat = flat
        self.rename_chain = rename_chain

    def anarci_parser(self, sequence):
        # scheme choose from [imgt, chothia, kabat, martin]
        scheme=self.scheme
        resnum_list = []
        # print(anarci.anarci([("chain", sequence)], scheme=scheme))
        anarci_result = anarci.anarci([("chain", sequence)], scheme=scheme)
        
        if not anarci_result[0][0]:
            return None, None
        else:
            chain_id = anarci_result[1][0][0]['chain_type']
            anarci_seq = anarci_result[0][0][0][0]
        for resi in anarci_seq:
            if resi[1] == "-":
                continue
            resnum_list.append(resi[0])
        return chain_id, resnum_list

    def gen_replace_dict(self, struct, chain_list=None) -> dict:
        
        # parser = PDBParser(PERMISSIVE=1)
        # struct = parser.get_structure("input", in_pdb)
        resnum_dict = {}
        for chain in struct.get_chains():
            if (chain.id in chain_list) or (chain_list is None):
            # chain_id = f"{os.path.basename(in_pdb).split('.')[0]}_{chain.id}"
                chain_sequence = ""
                for residue in chain.get_residues():
                    if residue._id[0] == ' ': # skip HETATM
                        chain_sequence += index_to_one(three_to_index(residue.get_resname()))
                chain_id, resnum_list = self.anarci_parser(chain_sequence)
                if not resnum_list:
                    continue
                resnum_dict[chain.id] = (chain_id, resnum_list)
        # SeqIO.write(seq_list, in_pdb.split('.')[0]+".fasta", "fasta")
        return resnum_dict

    def process(self, in_pdb, out_pdb=None, chain_list:list=None):
        #Write out selected portion to filename.
        if not out_pdb:
            out_pdb = f"{in_pdb.split('.')[0]}_{self.scheme}.pdb"
        # chain_list = list(chain_ids)
        parser = PDBParser(PERMISSIVE=1, QUIET=True)
        structure = parser.get_structure("input", in_pdb)
        # sele = select_chain(chain_list)
        
        new_resnum = self.gen_replace_dict(structure, chain_list)
        structure = self._renumber(structure, new_resnum)
        
        io = PDBIO()
        sele = select_chain()
        io.set_structure(structure)
        io.save(out_pdb, sele)
        return out_pdb  

    def _renumber(self, struct, new_resnum:dict):
        """
        flat: if True, numbering from 1 with no altloc, otherwise, number using specified scheme (e.g. IMGT)
        """
        # print(new_resnum)
        #loop 1: renumber residues to negative number to avoid errors
        chain_id = ""
        for residue in struct.get_residues():
            chain = residue.get_parent()
            if not chain.id in new_resnum.keys():
                # keep non-antibody chains untouched
                continue

            if chain_id != chain.get_id():
                chain_id = chain.get_id()
                residue_id = -1
            #print(chain.get_id(), residue_id)
            residue.id=(' ',residue_id,' ')
            residue_id -= 1
        
        #loop 2: start new numbering

        if self.flat:
            for chain in struct.get_chains():
                if chain.id in new_resnum.keys():
                    num_list = new_resnum[chain.id][1]
                    
                    for i, (resnum, residue) in enumerate(zip(num_list, chain.get_residues())):
                        residue.id = (' ', i+1, ' ') # this is for renumber using index
                        if resnum[0] in self.CDR_indexes:
                            residue.segid = "CDR"
                else:
                    continue

        else:
            for chain in struct.get_chains():
                if chain.id in new_resnum.keys():
                    num_list = new_resnum[chain.id][1]

                    for resnum, residue in zip(num_list, chain.get_residues()):
                        residue.id=(' ', resnum[0], resnum[1]) # this is for renumber using scheme
                else:
                    continue
        
        # loop 3: rename chain to H and L
        if self.rename_chain:
            for chain in struct.get_chains():
                if chain.id in new_resnum.keys():
                    chain.id = new_resnum[chain.id][0]
                else:
                    continue
        
        return struct

# low level parser
def block_by_chain_and_resnum(pdb_in, pdb_out, chain_resnum_dict):

    with open(pdb_in, "r") as h_in, open(pdb_out, "w") as h_out:
        for line in h_in.readlines():
            if line.startswith("ATOM"):
                chain_id = line[21]
                resnum = int(line[22:26])
                restype = line[55:57]
                # print(chain_id, resnum, restype)
                if chain_id in chain_resnum_dict.keys():
                    if not resnum in chain_resnum_dict[chain_id]:
                        line = line[:55] + "19" + line[57:]
                # print(line)
            else:
                pass
            h_out.write(line)

def block_by_seqid(pdb_in, pdb_out):

    with open(pdb_in, "r") as h_in, open(pdb_out, "w") as h_out:
        for line in h_in.readlines():
            if line.startswith("ATOM"):
                # chain_id = line[21]
                # resnum = int(line[22:26])
                restype = line[72:76].strip()
                # print(chain_id, resnum, restype)
                if restype != "CDR":
                    line = line[:55] + "19" + line[57:]
                # print(line)
            else:
                pass
            h_out.write(line)

if __name__ == "__main__":
    import argparse

    argparser = argparse.ArgumentParser(description="renumber mAb pdb")
    argparser.add_argument("in_pdb", help="input pdb")
    argparser.add_argument("--out_pdb", required=False, default=None, help="output pdb")
    argparser.add_argument("--scheme", required=False, default="imgt", help="imgt, chothia, kabat, martin")
    argparser.add_argument("--flat", action="store_true", help="flat numbering", default=False)
    argparser.add_argument("--rename_chain", action="store_true", help="rename chain to H and L", default=False)
    args = argparser.parse_args()
    # print(args)
    in_pdb = args.in_pdb
    out_pdb = args.out_pdb
    scheme = args.scheme
    flat = args.flat
    
    renumber_mAb(flat=flat, scheme=scheme).process(in_pdb, out_pdb)