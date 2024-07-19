#!/usr/bin/env python3
# fix missing heavry atoms and hydrogen atoms for docking

from pdbfixer import PDBFixer
from openmm.app import PDBFile
from pathlib import Path

# def auto_naming(ori_filename, target_path, out_filename):
    
#     path = Path(target_path).joinpath(Path(ori_filename).name)
    
#     return str(path.with_stem(f"{path.stem}_noH")), str(path.with_stem(f"{path.stem}_H"))

def fix(input_pdb, output_pdb, no_H=False):

    fixer = PDBFixer(filename=input_pdb)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    if not no_H:
        fixer.addMissingHydrogens(7.0)
    
    PDBFile.writeFile(fixer.topology, fixer.positions, open(output_pdb, 'w'), keepIds=True)

    return

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="fix missing heavry atoms and hydrogen atoms for docking")
    parser.add_argument("input_pdb", help="input pbd file that need dice")
    parser.add_argument("output_path", nargs="?", default=None, help="path to output pdb file, default = input file path")
    parser.add_argument("output_filename", nargs="?", default=None, help="output pdb file name, default = <input file name>_fix.pdb")
    parser.add_argument("-NH", "--no_hydrogen", action='store_true', help="don't add hydrogen atoms")
    args = parser.parse_args()

    if args.output_path is None:
        output_path = Path(args.input_pdb).parent
    else:
        output_path = args.output_path

    if args.output_filename is None:
        output_filename = Path(args.input_pdb).stem + "_fix.pdb"
    else:
        output_filename = args.output_filename

    output_file = str(Path(output_path).joinpath(output_filename))
    
    fix(args.input_pdb, output_file, args.no_hydrogen)