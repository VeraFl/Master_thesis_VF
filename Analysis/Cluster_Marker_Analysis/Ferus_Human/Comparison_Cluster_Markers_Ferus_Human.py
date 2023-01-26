# -*- coding: utf-8 -*-
"""
Created on Thu May 19 13:57:17 2022

@author: veraf
"""

#Import modules
import os
################################################################################################################################
# Import Files for the Alpaca Clusters
file1 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Ferus_Human/Single_cluster_data/Ferus_Cluster_0.txt", "r")
Ferus_markers_0 = file1.readlines()
file1.close()

file2 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Ferus_Human/Single_cluster_data/Ferus_Cluster_1.txt", "r")
Ferus_markers_1 = file2.readlines()
file2.close()

file3 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Ferus_Human/Single_cluster_data/Ferus_Cluster_2.txt", "r")
Ferus_markers_2 = file3.readlines()
file3.close()

file4 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Ferus_Human/Single_cluster_data/Ferus_Cluster_3.txt", "r")
Ferus_markers_3 = file4.readlines()
file4.close()

file5 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Ferus_Human/Single_cluster_data/Ferus_Cluster_4.txt", "r")
Ferus_markers_4 = file5.readlines()
file5.close()

file6 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Ferus_Human/Single_cluster_data/Ferus_Cluster_5.txt", "r")
Ferus_markers_5 = file6.readlines()
file6.close()

"""
file6 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Ferus_Human/Single_cluster_data/Ferus_Cluster_6.txt", "r")
Ferus_markers_6 = file6.readlines()
file6.close()

file6 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Ferus_Human/Single_cluster_data/Ferus_Cluster_7.txt", "r")
Ferus_markers_7 = file6.readlines()
file6.close()

file6 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Ferus_Human/Single_cluster_data/Ferus_Cluster_8.txt", "r")
Ferus_markers_8 = file6.readlines()
file6.close()
"""


# Import Files for the Human Clusters
file7 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Ferus_Human/Single_cluster_data/Human_Cluster_0.txt", "r")
Human_markers_0 = file7.readlines()
file7.close()

file8 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Ferus_Human/Single_cluster_data/Human_Cluster_1.txt", "r")
Human_markers_1 = file8.readlines()
file8.close()

file9 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Ferus_Human/Single_cluster_data/Human_Cluster_2.txt", "r")
Human_markers_2 = file9.readlines()
file9.close()

file10 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Ferus_Human/Single_cluster_data/Human_Cluster_3.txt", "r")
Human_markers_3 = file10.readlines()
file10.close()

file11 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Ferus_Human/Single_cluster_data/Human_Cluster_4.txt", "r")
Human_markers_4 = file11.readlines()
file11.close()

file12 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Ferus_Human/Single_cluster_data/Human_Cluster_5.txt", "r")
Human_markers_5 = file12.readlines()
file12.close()  

file12 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Ferus_Human/Single_cluster_data/Human_Cluster_6.txt", "r")
Human_markers_6 = file12.readlines()
file12.close() 

      
################################################################################################################################
#compare the marker list of the Ferus cluster X to all the human clusters
#Cluster 0
def compare_clusters(Ferus, Human):
    n = 0
    same_marker = []
    for gene in Ferus[1:]:
        for gene1 in Human[1:]:
            if gene == gene1:
                if gene not in same_marker:
                    same_marker.append(gene)
                    n += 1
    same_count = 0
    for el in same_marker:
        same_count += 1
    return(same_marker, same_count)


#Cluster 0 Ferus
Ferus_0_Human_0 = compare_clusters(Ferus_markers_0, Human_markers_0)
Ferus_0_Human_1 = compare_clusters(Ferus_markers_0, Human_markers_1)
Ferus_0_Human_2 = compare_clusters(Ferus_markers_0, Human_markers_2)
Ferus_0_Human_3 = compare_clusters(Ferus_markers_0, Human_markers_3)
Ferus_0_Human_4 = compare_clusters(Ferus_markers_0, Human_markers_4)
Ferus_0_Human_5 = compare_clusters(Ferus_markers_0, Human_markers_5)
Ferus_0_Human_6 = compare_clusters(Ferus_markers_0, Human_markers_6)


#Cluster 1 Ferus
Ferus_1_Human_0 = compare_clusters(Ferus_markers_1, Human_markers_0)
Ferus_1_Human_1 = compare_clusters(Ferus_markers_1, Human_markers_1)
Ferus_1_Human_2 = compare_clusters(Ferus_markers_1, Human_markers_2)
Ferus_1_Human_3 = compare_clusters(Ferus_markers_1, Human_markers_3)
Ferus_1_Human_4 = compare_clusters(Ferus_markers_1, Human_markers_4)
Ferus_1_Human_5 = compare_clusters(Ferus_markers_1, Human_markers_5)
Ferus_1_Human_6 = compare_clusters(Ferus_markers_1, Human_markers_6)


#Cluster 2 Ferus
Ferus_2_Human_0 = compare_clusters(Ferus_markers_2, Human_markers_0)
Ferus_2_Human_1 = compare_clusters(Ferus_markers_2, Human_markers_1)
Ferus_2_Human_2 = compare_clusters(Ferus_markers_2, Human_markers_2)
Ferus_2_Human_3 = compare_clusters(Ferus_markers_2, Human_markers_3)
Ferus_2_Human_4 = compare_clusters(Ferus_markers_2, Human_markers_4)
Ferus_2_Human_5 = compare_clusters(Ferus_markers_2, Human_markers_5)
Ferus_2_Human_6 = compare_clusters(Ferus_markers_2, Human_markers_6)


#Cluster 3 Ferus
Ferus_3_Human_0 = compare_clusters(Ferus_markers_3, Human_markers_0)
Ferus_3_Human_1 = compare_clusters(Ferus_markers_3, Human_markers_1)
Ferus_3_Human_2 = compare_clusters(Ferus_markers_3, Human_markers_2)
Ferus_3_Human_3 = compare_clusters(Ferus_markers_3, Human_markers_3)
Ferus_3_Human_4 = compare_clusters(Ferus_markers_3, Human_markers_4)
Ferus_3_Human_5 = compare_clusters(Ferus_markers_3, Human_markers_5)
Ferus_3_Human_6 = compare_clusters(Ferus_markers_3, Human_markers_6)


#Cluster 4 Ferus
Ferus_4_Human_0 = compare_clusters(Ferus_markers_4, Human_markers_0)
Ferus_4_Human_1 = compare_clusters(Ferus_markers_4, Human_markers_1)
Ferus_4_Human_2 = compare_clusters(Ferus_markers_4, Human_markers_2)
Ferus_4_Human_3 = compare_clusters(Ferus_markers_4, Human_markers_3)
Ferus_4_Human_4 = compare_clusters(Ferus_markers_4, Human_markers_4)
Ferus_4_Human_5 = compare_clusters(Ferus_markers_4, Human_markers_5)
Ferus_4_Human_6 = compare_clusters(Ferus_markers_4, Human_markers_6)


#Cluster 5 Ferus
Ferus_5_Human_0 = compare_clusters(Ferus_markers_5, Human_markers_0)
Ferus_5_Human_1 = compare_clusters(Ferus_markers_5, Human_markers_1)
Ferus_5_Human_2 = compare_clusters(Ferus_markers_5, Human_markers_2)
Ferus_5_Human_3 = compare_clusters(Ferus_markers_5, Human_markers_3)
Ferus_5_Human_4 = compare_clusters(Ferus_markers_5, Human_markers_4)
Ferus_5_Human_5 = compare_clusters(Ferus_markers_5, Human_markers_5)
Ferus_5_Human_6 = compare_clusters(Ferus_markers_5, Human_markers_6)


"""
#Cluster 6 Ferus
Ferus_6_Human_0 = compare_clusters(Ferus_markers_6, Human_markers_0)
Ferus_6_Human_1 = compare_clusters(Ferus_markers_6, Human_markers_1)
Ferus_6_Human_2 = compare_clusters(Ferus_markers_6, Human_markers_2)
Ferus_6_Human_3 = compare_clusters(Ferus_markers_6, Human_markers_3)
Ferus_6_Human_4 = compare_clusters(Ferus_markers_6, Human_markers_4)
Ferus_6_Human_5 = compare_clusters(Ferus_markers_6, Human_markers_5)
Ferus_6_Human_6 = compare_clusters(Ferus_markers_6, Human_markers_6)


#Cluster 7 Ferus
Ferus_7_Human_0 = compare_clusters(Ferus_markers_7, Human_markers_0)
Ferus_7_Human_1 = compare_clusters(Ferus_markers_7, Human_markers_1)
Ferus_7_Human_2 = compare_clusters(Ferus_markers_7, Human_markers_2)
Ferus_7_Human_3 = compare_clusters(Ferus_markers_7, Human_markers_3)
Ferus_7_Human_4 = compare_clusters(Ferus_markers_7, Human_markers_4)
Ferus_7_Human_5 = compare_clusters(Ferus_markers_7, Human_markers_5)
Ferus_7_Human_6 = compare_clusters(Ferus_markers_7, Human_markers_6)


#Cluster 8 Ferus
Ferus_8_Human_0 = compare_clusters(Ferus_markers_8, Human_markers_0)
Ferus_8_Human_1 = compare_clusters(Ferus_markers_8, Human_markers_1)
Ferus_8_Human_2 = compare_clusters(Ferus_markers_8, Human_markers_2)
Ferus_8_Human_3 = compare_clusters(Ferus_markers_8, Human_markers_3)
Ferus_8_Human_4 = compare_clusters(Ferus_markers_8, Human_markers_4)
Ferus_8_Human_5 = compare_clusters(Ferus_markers_8, Human_markers_5)
Ferus_8_Human_6 = compare_clusters(Ferus_markers_8, Human_markers_6)
"""




print(Ferus_0_Human_0[0])


os.chdir("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Ferus_Human/Ferus_Human_Comparison") 
###############################################################################
#Cluster 0 Ferus
f = open("Ferus_0_Human_0_markers.txt", "w")
for el in Ferus_0_Human_0[0]:
    f.write(el)
f.close()

f = open("Ferus_0_Human_1_markers.txt", "w")
for el in Ferus_0_Human_1[0]:
    f.write(el)
f.close()

f = open("Ferus_0_Human_2_markers.txt", "w")
for el in Ferus_0_Human_2[0]:
    f.write(el)
f.close()

f = open("Ferus_0_Human_3_markers.txt", "w")
for el in Ferus_0_Human_3[0]:
    f.write(el)
f.close()

f = open("Ferus_0_Human_4_markers.txt", "w")
for el in Ferus_0_Human_4[0]:
    f.write(el)
f.close()

f = open("Ferus_0_Human_5_markers.txt", "w")
for el in Ferus_0_Human_5[0]:
    f.write(el)
f.close()

f = open("Ferus_0_Human_6_markers.txt", "w")
for el in Ferus_0_Human_6[0]:
    f.write(el)
f.close()


###############################################################################
#Cluster 1 Ferus
f = open("Ferus_1_Human_0_markers.txt", "w")
for el in Ferus_1_Human_0[0]:
    f.write(el)
f.close()

f = open("Ferus_1_Human_1_markers.txt", "w")
for el in Ferus_1_Human_1[0]:
    f.write(el)
f.close()

f = open("Ferus_1_Human_2_markers.txt", "w")
for el in Ferus_1_Human_2[0]:
    f.write(el)
f.close()

f = open("Ferus_1_Human_3_markers.txt", "w")
for el in Ferus_1_Human_3[0]:
    f.write(el)
f.close()

f = open("Ferus_1_Human_4_markers.txt", "w")
for el in Ferus_1_Human_4[0]:
    f.write(el)
f.close()

f = open("Ferus_1_Human_5_markers.txt", "w")
for el in Ferus_1_Human_5[0]:
    f.write(el)
f.close()

f = open("Ferus_1_Human_6_markers.txt", "w")
for el in Ferus_1_Human_6[0]:
    f.write(el)
f.close()


###############################################################################
#Cluster 2 Ferus
f = open("Ferus_2_Human_0_markers.txt", "w")
for el in Ferus_2_Human_0[0]:
    f.write(el)
f.close()

f = open("Ferus_2_Human_1_markers.txt", "w")
for el in Ferus_2_Human_1[0]:
    f.write(el)
f.close()

f = open("Ferus_2_Human_2_markers.txt", "w")
for el in Ferus_2_Human_2[0]:
    f.write(el)
f.close()

f = open("Ferus_2_Human_3_markers.txt", "w")
for el in Ferus_2_Human_3[0]:
    f.write(el)
f.close()

f = open("Ferus_2_Human_4_markers.txt", "w")
for el in Ferus_2_Human_4[0]:
    f.write(el)
f.close()

f = open("Ferus_2_Human_5_markers.txt", "w")
for el in Ferus_2_Human_5[0]:
    f.write(el)
f.close()

f = open("Ferus_2_Human_6_markers.txt", "w")
for el in Ferus_2_Human_6[0]:
    f.write(el)
f.close()



###############################################################################
#Cluster 3 Ferus
f = open("Ferus_3_Human_0_markers.txt", "w")
for el in Ferus_3_Human_0[0]:
    f.write(el)
f.close()

f = open("Ferus_3_Human_1_markers.txt", "w")
for el in Ferus_3_Human_1[0]:
    f.write(el)
f.close()

f = open("Ferus_3_Human_2_markers.txt", "w")
for el in Ferus_3_Human_2[0]:
    f.write(el)
f.close()

f = open("Ferus_3_Human_3_markers.txt", "w")
for el in Ferus_3_Human_3[0]:
    f.write(el)
f.close()

f = open("Ferus_3_Human_4_markers.txt", "w")
for el in Ferus_3_Human_4[0]:
    f.write(el)
f.close()

f = open("Ferus_3_Human_5_markers.txt", "w")
for el in Ferus_3_Human_5[0]:
    f.write(el)
f.close()

f = open("Ferus_3_Human_6_markers.txt", "w")
for el in Ferus_3_Human_6[0]:
    f.write(el)
f.close()


###############################################################################
#Cluster 4 Ferus
f = open("Ferus_4_Human_0_markers.txt", "w")
for el in Ferus_4_Human_0[0]:
    f.write(el)
f.close()

f = open("Ferus_4_Human_1_markers.txt", "w")
for el in Ferus_4_Human_1[0]:
    f.write(el)
f.close()

f = open("Ferus_4_Human_2_markers.txt", "w")
for el in Ferus_4_Human_2[0]:
    f.write(el)
f.close()

f = open("Ferus_4_Human_3_markers.txt", "w")
for el in Ferus_4_Human_3[0]:
    f.write(el)
f.close()

f = open("Ferus_4_Human_4_markers.txt", "w")
for el in Ferus_4_Human_4[0]:
    f.write(el)
f.close()

f = open("Ferus_4_Human_5_markers.txt", "w")
for el in Ferus_4_Human_5[0]:
    f.write(el)
f.close()

f = open("Ferus_4_Human_6_markers.txt", "w")
for el in Ferus_4_Human_6[0]:
    f.write(el)
f.close()


###############################################################################
#Cluster 5 Ferus
f = open("Ferus_5_Human_0_markers.txt", "w")
for el in Ferus_5_Human_0[0]:
    f.write(el)
f.close()

f = open("Ferus_5_Human_1_markers.txt", "w")
for el in Ferus_5_Human_1[0]:
    f.write(el)
f.close()

f = open("Ferus_5_Human_2_markers.txt", "w")
for el in Ferus_5_Human_2[0]:
    f.write(el)
f.close()

f = open("Ferus_5_Human_3_markers.txt", "w")
for el in Ferus_5_Human_3[0]:
    f.write(el)
f.close()

f = open("Ferus_5_Human_4_markers.txt", "w")
for el in Ferus_5_Human_4[0]:
    f.write(el)
f.close()

f = open("Ferus_5_Human_5_markers.txt", "w")
for el in Ferus_5_Human_5[0]:
    f.write(el)
f.close()

f = open("Ferus_5_Human_6_markers.txt", "w")
for el in Ferus_5_Human_6[0]:
    f.write(el)
f.close()

"""
###############################################################################
#Cluster 6 Ferus
f = open("Ferus_6_Human_0_markers.txt", "w")
for el in Ferus_6_Human_0[0]:
    f.write(el)
f.close()

f = open("Ferus_6_Human_1_markers.txt", "w")
for el in Ferus_6_Human_1[0]:
    f.write(el)
f.close()

f = open("Ferus_6_Human_2_markers.txt", "w")
for el in Ferus_6_Human_2[0]:
    f.write(el)
f.close()

f = open("Ferus_6_Human_3_markers.txt", "w")
for el in Ferus_6_Human_3[0]:
    f.write(el)
f.close()

f = open("Ferus_6_Human_4_markers.txt", "w")
for el in Ferus_6_Human_4[0]:
    f.write(el)
f.close()

f = open("Ferus_6_Human_5_markers.txt", "w")
for el in Ferus_6_Human_5[0]:
    f.write(el)
f.close()

f = open("Ferus_6_Human_6_markers.txt", "w")
for el in Ferus_6_Human_6[0]:
    f.write(el)
f.close()

###############################################################################
#Cluster 7 Ferus
f = open("Ferus_7_Human_0_markers.txt", "w")
for el in Ferus_7_Human_0[0]:
    f.write(el)
f.close()

f = open("Ferus_7_Human_1_markers.txt", "w")
for el in Ferus_7_Human_1[0]:
    f.write(el)
f.close()

f = open("Ferus_7_Human_2_markers.txt", "w")
for el in Ferus_7_Human_2[0]:
    f.write(el)
f.close()

f = open("Ferus_7_Human_3_markers.txt", "w")
for el in Ferus_7_Human_3[0]:
    f.write(el)
f.close()

f = open("Ferus_7_Human_4_markers.txt", "w")
for el in Ferus_7_Human_4[0]:
    f.write(el)
f.close()

f = open("Ferus_7_Human_5_markers.txt", "w")
for el in Ferus_7_Human_5[0]:
    f.write(el)
f.close()

f = open("Ferus_7_Human_6_markers.txt", "w")
for el in Ferus_7_Human_6[0]:
    f.write(el)
f.close()

###############################################################################
#Cluster 8 Ferus
f = open("Ferus_8_Human_0_markers.txt", "w")
for el in Ferus_8_Human_0[0]:
    f.write(el)
f.close()

f = open("Ferus_8_Human_1_markers.txt", "w")
for el in Ferus_8_Human_1[0]:
    f.write(el)
f.close()

f = open("Ferus_8_Human_2_markers.txt", "w")
for el in Ferus_8_Human_2[0]:
    f.write(el)
f.close()

f = open("Ferus_8_Human_3_markers.txt", "w")
for el in Ferus_8_Human_3[0]:
    f.write(el)
f.close()

f = open("Ferus_8_Human_4_markers.txt", "w")
for el in Ferus_8_Human_4[0]:
    f.write(el)
f.close()

f = open("Ferus_8_Human_5_markers.txt", "w")
for el in Ferus_8_Human_5[0]:
    f.write(el)
f.close()

f = open("Ferus_8_Human_6_markers.txt", "w")
for el in Ferus_8_Human_6[0]:
    f.write(el)
f.close()

###############################################################################
"""