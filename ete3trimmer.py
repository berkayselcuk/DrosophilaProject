# -*- coding: utf-8 -*-
"""
Created on Sat Nov  9 23:58:01 2019

@author: Berkay Selcuk
"""
#USAGE:
#When you run this on a terminal the format is:
#ete3trimmer.py "tree_directory" "gene_of_interest "alignment directory"
#The species is set to HUMAN but you can change it by changing the given taxID and _HUMAN.
#import sys
#from ete3 import PhyloTree
#if __name__ =="__main__":
#        tree=sys.argv[1]
#        gene_of_interest=sys.argv[2]
#        align_doc=sys.argv[3]
    
gene_of_interest="OPRM"
tree=r"D:\Users\suuser\Desktop\RAxML_bestTree.P35372_tree_27112019"
#In this part we obtain the largest clade that contains our gene of interest.    
t=PhyloTree(tree)
R = t.get_midpoint_outgroup()
t.set_outgroup(R) #for rooting of the tree
for node in t.traverse():
    if node.is_leaf():
        continue
    name = node.name
    leaf_list=[]
    gene_check=0
    human_check=0
    ext=0
    for leaf in node.get_leaves():
        leaf_list.append(str(leaf)[3:])
        taxID=str(leaf).split("|")[-1]
        if taxID=="9606":
            human_check+=1
            if human_check>1:
                print ("More than one human")
                break
        gene=str(leaf).split("|")[-2]
        if gene=="{}_HUMAN".format(gene_of_interest):
            gene_check=1
        if gene_check==1 and human_check==1:
            print (name)
            print (len(node.get_leaves()))
            print (len(leaf_list))
            t.prune(leaf_list)
            print ("Sub tree is obtained")
            #print (t)
            ext=1
            break
    if ext==1:
        break


#To obtain all of the leaf nodes from the root
for node in t.traverse():
    total_leave_list=node.get_leaves()
    break

#This is the part where we trim
remove_list=[]
for node in t.traverse():
    prune_list=[]
    children=node.get_children()
    if len(children)==1 or node.is_leaf():
        continue
    leaf1=children[0].get_leaves()
    tax1=[]
    leaf2=children[1].get_leaves()
    tax2=[]
    for i in leaf1:
        taxID=str(i).split("|")[-1]
        if taxID=="\n--":
            continue
        if taxID not in tax1:
            tax1.append(taxID)
    for a in leaf2:
        taxID=str(a).split("|")[-1]
        if taxID=="\n--":
            continue
        if taxID not in tax2:
            tax2.append(taxID)
    go=0 #if there are same species in child nodes this becomes 1 and we continue
    tax_list=[]
    for tax in tax1:
        if tax in tax2:
            go=1
            tax_list.append(tax)
    if go==1:
        print ("Correct!")
        print (tax_list)
        ref_count=0
        tar_count=0
        ref_list=[]
        tar_list=[]
        for tax in tax_list:
            for leaf in leaf1:
                if str(leaf).split("|")[-1]==tax:
                    reference=leaf
                    ref_list.append(reference)
                    break
            for leaf in leaf2:
                if str(leaf).split("|")[-1]==tax:
                    target=leaf
                    tar_list.append(target)
                    break
            if reference.get_distance(node)>target.get_distance(node):
                tar_count+=1
            else:
                ref_count+=1
        print (ref_count)
        if ref_count>tar_count and 1.2*len(tax1)>=len(tax2):
            for lif in leaf2:
                print ("{} is the leaf.".format(lif))
                remove_list.append(lif)
        elif ref_count>tar_count:
            for i in tar_list:
                remove_list.append(i)
  
        elif tar_count>ref_count and 1.2*len(tax2)>=len(tax1):
            for lif in leaf1:
                print ("{} is the leaf.".format(lif))
                remove_list.append(lif)
        elif tar_count>ref_count:
            for i in ref_list:
                remove_list.append(i)

    #print (remove_list)
    for leaf in remove_list:
        try:
            total_leave_list.remove(leaf)
        except ValueError:
            continue

    t.prune(total_leave_list)
print ("Trimming is done!")
t.write(outfile="{}_Trimmed.txt".format(gene_of_interest))
print (t)
print ("Tree is saved to: DESKTOP!!!" )