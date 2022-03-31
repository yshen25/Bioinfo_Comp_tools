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
# resi_to_b_dict = {1:0,2:0,3:0,4:0,5:0,6:0,7:0.00199601,8:0,9:0,10:0,11:0,12:0.00106821,13:0.00106821,14:0.00199626,15:0.00213641,16:0,17:0,18:0,19:0.00306465,20:0,21:0,22:0,23:0,24:0.00199626,25:0,26:0,27:0,28:0.00106814,29:0,30:0,31:0,32:0,33:0.00835438,34:0.00306542,35:0,36:0,37:0,38:0,39:0,40:0,41:0.00532326,42:0.00213655,43:0,44:0,45:0,46:0,47:0.00106828,48:0,49:0.00106828,50:0,51:0,52:0,53:0.00106828,54:0,55:0.0602486,56:0,57:0,58:0,59:0,60:0.00106821,61:0,62:0,63:0,64:0,65:0,66:0,67:0,68:0,69:0.00106814,70:0,71:0.02291717,72:0,73:0.00106814,74:0.00609982,75:0,76:0,77:0.00106814,78:0,79:0,80:0.0045274,81:0.00106814,82:0,83:0,84:0,85:0,86:0.00199626,87:0.00306484,88:0,89:0,90:0,91:0.00908489,92:0.00306465,93:0,94:0.00199663,95:0,96:0.00607517,97:0.00287086,98:0.00199588,99:0,100:0.00106814,101:0,102:0,103:0.00974945,104:0.00477946,105:0,106:0.00106807,107:0.00638924,108:0.01662186,109:0.00106801,110:0,111:0,112:0,113:0.00199588,114:0.00287104,115:0.00287122,116:0,117:0.00106814,118:0.00371194,119:0.00199601,120:0.00399249,121:0,122:0,123:0.00371218,124:0,125:0.00106807}
def csv_to_dict(InCSV):
    in_df = pd.read_csv()
    return

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
                line = f"{line[0:60]}{str(round(replace_dict[resi], 4)).rjust(6)}{line[66:]}"
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
    main()