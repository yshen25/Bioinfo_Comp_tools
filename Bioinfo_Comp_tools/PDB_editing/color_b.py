#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description: change b factor of PDB file to customized value, to show color range in visualization software
-----------------------------------------------------
Version: 1.1, read residue number and value from CSV file
Author: Yue (Shawn) Shen
Date: Mar 2022
"""
import pandas as pd

def csv_to_dict(InCSV):
    in_df = pd.read_csv(InCSV)
    return dict(in_df.values)

def resi2b(replace_dict, resiNum):
    if resiNum in replace_dict:
        return replace_dict[resiNum]
    else:
        return 0

def replace_b(replace_dict, PDB, outPDB):

    with open(PDB, "r") as in_f, open(outPDB, "w") as out_f:
        for line in in_f:
            if line.startswith("ATOM"):
                resi = int(line[23:26].strip())
                line = f"{line[0:60]}{str(round(resi2b(replace_dict, resi), 4)).rjust(6)}{line[66:]}"
            elif line.startswith("ANISOU"):
                continue
            out_f.write(line)

    return

def main(InCSV, PDB, outPDB=None):
    if outPDB is None:
        outPDB = PDB.split('.')[0] + "_colored.pdb"

    resi_value_dict = csv_to_dict(InCSV)

    replace_b(resi_value_dict, PDB, outPDB)

    return

if __name__ == "__main__":
    main("Book1.csv", "1i4f_Crown.pdb")