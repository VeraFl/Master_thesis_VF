# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 11:04:12 2022

@author: veraf
"""

file1 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cell_Cycle_State_Markers_Species/G2M_phase_genes", "r")
Human_markers = file1.readlines()
file1.close()


file2 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cell_Cycle_State_Markers_Species/S_phase_genes", "r")
Camel_markers = file2.readlines()
file2.close()

file3 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cell_Cycle_State_Markers_Species/S_phase_genes", "r")
Camel_markers = file2.readlines()
file2.close()

#compare the two unique contig lists gtf and fasta
n = 0
same_marker = []
for gene in Human_markers:
    for gene1 in Camel_markers:
        if gene == gene1:
            if gene not in same_marker:
                same_marker.append(gene)
                n += 1
            
                    


same_count = 0
for el in same_marker:
    same_count += 1
    


f = open("Comparison_Dromedary_Human_Cluster_1.txt", "w")
for el in same_marker:
    f.write(el)
f.close()
