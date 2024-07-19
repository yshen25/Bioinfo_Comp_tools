#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description: search PubChem ID against Entrez database for corresponding CAS ID
-----------------------------------------------------
Version: 1.0
Author: Yue (Shawn) Shen
Date: Apr 2020
"""

from Bio import Entrez
import time
import re
Entrez.email = "" # your email here to use Entrez

PubChem_file = "" # file of query PubChem IDs, one ID each line
ID_file = "" # file of hit CAS IDs

def id_to_cas():
    cas_pattern_top = re.compile(r"^CAS-\d+-\d+-\d$")
    cas_pattern_alt = re.compile(r"^\d+-\d+-\d$")
    with open(PubChem_file, "r") as list_file, open(ID_file, "w") as id_file:
        last_line = ""
        for line in list_file:
            cas = ""
            hit_name = ""
            hit_id = line.strip()
            if hit_id == "0":
                hit_name = "0"
                #print("0")
            elif hit_id != last_line:
                try:
                    handle2 = Entrez.esummary(db="pccompound", id=hit_id)
                    record2 = Entrez.read(handle2)
                    try:
                        hit_name = record2[0]['MeSHHeadingList'][0]
                    except:
                        try:
                            hit_name = record2[0]['SynonymList'][0]
                        except:
                            pass
                except:
                    try:
                        handle2 = Entrez.esummary(db="pcsubstance", id=hit_id)
                        record2 = Entrez.read(handle2)
                    except:
                        pass
                #print(f"{hit_id} - {hit_name}")
                for item in record2[0]['SynonymList']:
                    match_top = re.match(cas_pattern_top, item)
                    if match_top:
                        cas = match_top.group(0)
                        break
                    match_alt = re.match(cas_pattern_alt, item)
                    if match_alt:
                        cas = match_alt.group(0)
                        break
            else:
                pass
                #print("...")
            last_line = hit_id
            print(f"{hit_id}\t{cas}\t{hit_name}")
            id_file.write(f"{hit_id}\t{cas}\t{hit_name}\n")
        #sys.exit()
        time.sleep(1)
    handle2.close()
    return 0

id_to_cas()