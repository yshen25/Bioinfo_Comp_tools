#!/usr/bin/env python3
#used to trim pdb file
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import extract
def dice_pdb(filename, start, stop, outfile=None):
    if outfile == None:
        outfile = f"trim_{filename}"
    parser = PDBParser(PERMISSIVE=1)
    structure_id = "thread"
    structure = parser.get_structure(structure_id, filename)
    extract(structure, "A", start, stop, outfile)
    return outfile

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="dicing pdb between start and stop position")
    parser.add_argument("input_pdb", help="input pbd file that need dice")
    parser.add_argument("START", type=int, help="START position of peptide being extracted, in pdb numbering(start from 1)")
    parser.add_argument("STOP", type=int, help="STOP position of peptide being extracted, in pdb numbering(start from 1)")
    parser.add_argument("output_pdb", nargs="?", default=None, help="output pdb file name, default = trim_<input file name>")
    args = parser.parse_args()
    print(dice_pdb(args.input_pdb, args.START, args.STOP, args.output_pdb))