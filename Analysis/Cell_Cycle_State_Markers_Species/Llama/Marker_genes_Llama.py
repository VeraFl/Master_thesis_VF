# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 11:04:12 2022

@author: veraf
"""

file1 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cell_Cycle_State_Markers_Species/Llama/G2M_phase_genes.txt", "r")
G2M_genes = file1.readlines()
file1.close()

file2 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cell_Cycle_State_Markers_Species/Llama/S_phase_genes.txt", "r")
S_genes = file2.readlines()
file2.close()

file3 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cell_Cycle_State_Markers_Species/Llama/genes_list_llama_camel229E_inf.txt", "r")
inf_genes = file3.readlines()
file2.close()

file4 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cell_Cycle_State_Markers_Species/Llama/genes_list_llama_camel229E_mock.txt", "r")
mock_genes = file4.readlines()
file2.close()


##############################################################################
##INFECTED
#compare the G2M genes with the infected genes:
n = 0
same_marker_G2M_inf = []
for gene in G2M_genes:
    for gene1 in inf_genes:
        if gene == gene1:
            if gene not in same_marker_G2M_inf:
                same_marker_G2M_inf.append(gene)
                n += 1
            
#compare the S genes with the infected genes:
n1 = 0
same_marker_S_inf = []
for gene in S_genes:
    for gene1 in inf_genes:
        if gene == gene1:
            if gene not in same_marker_S_inf:
                same_marker_S_inf.append(gene)
                n1 += 1                 

###############################################################################
##MOCK
#compare the G2M genes with the mock genes:
n = 0
same_marker_G2M_mock = []
for gene in G2M_genes:
    for gene1 in mock_genes:
        if gene == gene1:
            if gene not in same_marker_G2M_mock:
                same_marker_G2M_mock.append(gene)
                n += 1
            
#compare the S genes with the mock genes:
n1 = 0
same_marker_S_mock = []
for gene in S_genes:
    for gene1 in mock_genes:
        if gene == gene1:
            if gene not in same_marker_S_mock:
                same_marker_S_mock.append(gene)
                n1 += 1 

###############################################################################
#How many of the cell cycle state markers are in the infected and mock gene lists?

same_marker_G2M_inf_count = 0
for el in same_marker_G2M_inf:
    same_marker_G2M_inf_count += 1
    
same_marker_S_inf_count = 0
for el in same_marker_S_inf:
    same_marker_S_inf_count += 1

same_marker_G2M_mock_count = 0
for el in same_marker_G2M_mock:
    same_marker_G2M_mock_count += 1

same_marker_S_mock_count = 0
for el in same_marker_S_mock:
    same_marker_S_mock_count += 1

###############################################################################
# Which cell cycle genes are not in the datasets? - As a control

to_be_excluded_in_S_genes = []
for gen1 in S_genes:
    if gen1 not in same_marker_S_mock:
        to_be_excluded_in_S_genes.append(gen1)

to_be_excluded_in_G2M_genes = []
for gen2 in G2M_genes:
    if gen2 not in same_marker_G2M_mock:
        to_be_excluded_in_G2M_genes.append(gen2)


###############################################################################
# What are the new S and G2M gene lists? - Write the to txt file to use in R

f = open("new_G2M_phase_genes.txt", "w")
for el1 in same_marker_G2M_mock:
    f.write(el1)
f.close()

fi = open("new_S_phase_genes.txt", "w")
for el2 in same_marker_S_mock:
    fi.write(el2)
fi.close()
