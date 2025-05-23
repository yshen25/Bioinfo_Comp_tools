#!/usr/bin/env python3
#simple threading script replace applications in ROSETTA
#thread on one chain and keep other chains untouched
from Bio import AlignIO
from Bio.PDB.PDBParser import PDBParser
from os.path import splitext
from Bio.PDB.Polypeptide import one_to_index, index_to_three

def main(tempalte_pdb, chain_id, align_file, output_pdb=None):
    if output_pdb == None:
        output_pdb = splitext(tempalte_pdb)[0] + "_thread.pdb"
    pdb_parser = PDBParser(PERMISSIVE=1)
    template_struct = pdb_parser.get_structure("template", tempalte_pdb)
    alignment = AlignIO.read(align_file, "fasta")
    
    with open(output_pdb, "w") as target_struct:
        atom_id = 1 #serial number for atoms in pdb file
        for chain in template_struct.get_chains():
            resi_id = 1 #residue id for target sequence
            if chain.id == chain_id:
                resi_id2 = 1 #residue id for template sequence
                for i in range(len(alignment[0])):
                    # iterate through alignment, alignment[:,i][0] is template, [1] is target
                    # check consistency of pdb and alignment
                    # print(alignment[:,i][0], three_to_one(chain[resi_id2].resname))
                    if alignment[:,i][0] == alignment[:,i][1]:
                        #if one location is same in both target and template, copy all information from template_pdb
                        for atom in chain[resi_id2]:
                            new_line = f"ATOM  {str(atom_id).rjust(5)} {atom.fullname}{atom.altloc}{chain[resi_id2].resname} {chain.id}{str(resi_id).rjust(4)}{chain[resi_id2].id[2]}   {str(atom.coord[0]).rjust(8)}{str(atom.coord[1]).rjust(8)}{str(atom.coord[2]).rjust(8)}  {atom.occupancy:.2f} {str(atom.bfactor).rjust(5)}          {str(atom.element).rjust(2)} \n"
                            target_struct.write(new_line)
                            atom_id += 1
                        resi_id += 1
                        resi_id2 += 1
                    elif alignment[:,i][1] == "-":
                        #if there's a gap in target sequence, skip that location
                        print(i+1,alignment[:,i],"skip")
                        resi_id2 += 1
                        '''
                        for atom in chain[resi_id2]:
                            new_line = f"ATOM  {str(atom_id).rjust(5)} {atom.fullname}{atom.altloc}{chain[resi_id2].resname} {chain.id}{str(resi_id).rjust(4)}{chain[resi_id2].id[2]}   {str(atom.coord[0]).rjust(8)}{str(atom.coord[1]).rjust(8)}{str(atom.coord[2]).rjust(8)}  {atom.occupancy:.2f} {str(atom.bfactor).rjust(5)}          {str(atom.element).rjust(2)} \n"
                            target_struct.write(new_line)
                            atom_id += 1
                        resi_id2 += 1'''
                    elif alignment[:,i][0] == "-":
                        #if there's a gap in template sequence, build a loop
                        print(i+1,alignment[:,i],"loop")
                        resi_id += 1
                    else:
                        #if the location is difference in target and template, only keep the backbone atoms
                        print(i+1,alignment[:,i],"backbone")
                        for atom in sorted(chain[resi_id2])[0:4]:
                            new_line = f"ATOM  {str(atom_id).rjust(5)} {atom.fullname}{atom.altloc}{index_to_three(one_to_index(alignment[:,i][1]))} {chain.id}{str(resi_id).rjust(4)}{chain[resi_id2].id[2]}   {str(atom.coord[0]).rjust(8)}{str(atom.coord[1]).rjust(8)}{str(atom.coord[2]).rjust(8)}  {atom.occupancy:.2f} {str(atom.bfactor).rjust(5)}          {str(atom.element).rjust(2)} \n"
                            target_struct.write(new_line)
                            atom_id += 1
                        resi_id += 1
                        resi_id2 += 1
            else:
                for residue in chain.get_residues():
                    for atom in residue.get_atoms():
                        new_line = f"ATOM  {str(atom_id).rjust(5)} {atom.fullname}{atom.altloc}{residue.resname} {chain.id}{str(resi_id).rjust(4)}{residue.id[2]}   {str(atom.coord[0]).rjust(8)}{str(atom.coord[1]).rjust(8)}{str(atom.coord[2]).rjust(8)}  {atom.occupancy:.2f} {str(atom.bfactor).rjust(5)}          {str(atom.element).rjust(2)} \n"
                        target_struct.write(new_line)
                        atom_id += 1
                    resi_id += 1
            new_line = f"TER   {str(atom_id-1).rjust(5)}      {new_line[17:20]} {chain.id}{str(resi_id-1).rjust(4)}\n"
            target_struct.write(new_line)
        target_struct.write("END\n")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="thread protein sequence to template pdb")
    parser.add_argument("template_pdb", help="template structure used for scaffold")
    parser.add_argument("chain", help="single chain identifier")
    parser.add_argument("alignment_file", help="alignment in fasta format with target sequence on top")
    parser.add_argument("output_pdb", nargs="?", default=None, help="output pdb file name, default = <input_pdb>_thread.pdb")
    args = parser.parse_args()
    main(args.template_pdb, args.chain, args.alignment_file, args.output_pdb)
