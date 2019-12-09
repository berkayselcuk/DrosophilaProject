# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 14:41:27 2019

@author: Berkay Selcuk
"""
from ete3 import PhyloTree
sample_list=[r"D:\Users\suuser\Desktop\ADRB2_Trimmed.txt",
             r"D:\Users\suuser\Desktop\ACM2_Trimmed.txt",
             r"D:\Users\suuser\Desktop\AA2AR_Trimmed.txt",
             r"D:\Users\suuser\Desktop\OPRM_Trimmed.txt"]
file_name=["ADRB2","ACM2","AA2AR","OPRM"]
for sample in sample_list:
    idx=sample_list.index(sample)
    file=file_name[idx]
    t=open(sample,"r")
    line=t.readline()
    tree=PhyloTree(line)
    R = tree.get_midpoint_outgroup()
    tree.set_outgroup(R) #for rooting of the tree
    for node in tree.traverse():
        leaves=node.get_children()
        #print (leaves)
        if len(leaves)<2 or node.is_leaf():
            continue
        #leaf1=leaves[0].get_leaves()
        #print (leaf1)
        leaf2=leaves[1].get_leaves()
        for i in leaf2:
            leaf_node=str(i).split("|")[-1]
            print (leaf_node)
            if leaf_node=="9606":
                node.swap_children()
                break
#tree.link_to_alignment(alignment="_subtree_MA17112019.fasta", alg_format="fasta")
#print ("Tree is combined with the alignment")
#tree.show()
    tree.write(outfile=f"{file}_Ordered.txt")
    print (tree)
