# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 17:32:56 2019

@author: Berkay Selcuk
"""
from ete3 import PhyloTree
def etealign(tree,MA):
    t=tree
    treefix=open(t,"r")
    t=treefix.readline().replace("'", "")
    tree=PhyloTree(t)
    print (tree)
    tree.link_to_alignment(alignment=MA, alg_format="fasta") 
    tree.show()
    
    
gene_list=["ADRB2","ACM2","AA2AR","OPRM"]    
for gene in gene_list:
    etealign(gene+"_Ordered",gene+"_Ordered_MA08112019")
    