#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description: re-order multiple sequence alignment according to phylogenetic tree, moving samilar sequences together
-----------------------------------------------------
Version: 1.1, Upgrade to Python3
Author: Yue (Shawn) Shen
Date: Mar 2022
"""

# libs
from pathlib import Path
from Bio import AlignIO, Phylo
from Bio.Align import MultipleSeqAlignment
import re, argparse

def get_option():
    #extract user input in command line
    parser = argparse.ArgumentParser(description = '[%(prog)s] a program that rearanging sequence aliment file according to phylogenetic tree', usage = '%(prog)s -a <alignment file> -t <tree file> -f <tree file format>\n [--help] for more information\n [--version] for version information', epilog = "Shawn, Oct 2017")
    parser.add_argument('-a',metavar='alignfile',help='alignment file name',required=True)
    parser.add_argument('-t',metavar='treefile',help='phylogenetic tree file name',required=True)
    parser.add_argument('-f',help='tree file format, select from the list',choices=['newick','nexus','nexml','phyloxml','cdao'],required=True)
    parser.add_argument('--version', action='version', version='--align2tree 0.2--')
    
    opt = parser.parse_args()
    return opt.a,opt.t,opt.f

def get_tree_label(TreeFile, TreeFormat, id_format=None):
    #read tree file and acquire the name of each leaf
    thandle = Phylo.read(TreeFile,TreeFormat)
    leaf_list = []
    if id_format:
        print(f"== regular expression to catch seg id: {id_format} ==")
        seq_id = re.compile(id_format)
    for leaf in thandle.get_terminals():
        if id_format:
            leaf_name = seq_id.match(leaf.name).group()
        else:
            leaf_name = leaf.name

        leaf_list.append(leaf_name)

    return leaf_list

def main(AlnFile, AlnFormat, TreeFile, TreeFormat, Output):
    # automatically assign output name if empty
    if Output is None:
        Output = Path(AlnFile).resolve().with_name(Path(AlnFile).resolve().stem + "_reindex" + Path(AlnFile).resolve().suffix)

    # read MSA file and build dictionary {id:record}
    alignments = AlignIO.read(AlnFile, AlnFormat)
    aln_dict = {}
    for aln in alignments:
        aln_dict[aln.id] = aln

    # the tree uses description rather than id as leaf name, use re to extract id name
    tree_label_format = r'\w*.\d{0,2}' # use None if no re needed

    # read tree file and get labels in order
    tree_labels = get_tree_label(TreeFile, TreeFormat, tree_label_format)
    
    # extract from dictionary using tree labels as keys
    reindex_aln = MultipleSeqAlignment([aln_dict[i] for i in tree_labels])
    print(reindex_aln)

    # write into new MSA file
    AlignIO.write(reindex_aln, Output, AlnFormat)
    
    print("========================\njob done")
    print(f"reindexed alignment file saved as: {Output}")

    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = 'reorder sequence aliment file according to phylogenetic tree, both files must contain same labels', epilog = "Shawn, Mar 2022")
    parser.add_argument("alignment", help="input alignment file that need reindex")
    parser.add_argument("aln_format", choices=['clustal','fasta','maf','mauve','nexus','phylip','stockholm'], help="alignment file format")
    parser.add_argument("tree_file", help="phylogenetic tree file used as template")
    parser.add_argument("tree_format", choices=['newick','nexus','nexml','phyloxml','cdao'], help="phylogenetic tree file format")
    parser.add_argument('--version', action='version', version='--align2tree 1.0--')
    parser.add_argument("out_aln", nargs="?", default=None, help="output alignment file name, default = <input file name>_reindex.aln")
    args = parser.parse_args()

    main(args.alignment, args.aln_format, args.tree_file, args.tree_format, args.out_aln)