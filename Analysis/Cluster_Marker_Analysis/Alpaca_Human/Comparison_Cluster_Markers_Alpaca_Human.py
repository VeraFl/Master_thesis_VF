# -*- coding: utf-8 -*-
"""
Created on Thu May 19 13:57:17 2022

@author: veraf
"""

#Import modules
import os
################################################################################################################################
# Import Files for the Alpaca Clusters
file1 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Human/Single_cluster_data/Alpaca_Cluster_0.txt", "r")
Alpaca_markers_0 = file1.readlines()
file1.close()

file2 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Human/Single_cluster_data/Alpaca_Cluster_1.txt", "r")
Alpaca_markers_1 = file2.readlines()
file2.close()

file3 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Human/Single_cluster_data/Alpaca_Cluster_2.txt", "r")
Alpaca_markers_2 = file3.readlines()
file3.close()

file4 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Human/Single_cluster_data/Alpaca_Cluster_3.txt", "r")
Alpaca_markers_3 = file4.readlines()
file4.close()

file5 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Human/Single_cluster_data/Alpaca_Cluster_4.txt", "r")
Alpaca_markers_4 = file5.readlines()
file5.close()

file6 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Human/Single_cluster_data/Alpaca_Cluster_5.txt", "r")
Alpaca_markers_5 = file6.readlines()
file6.close()

"""
file6 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Human/Single_cluster_data/Alpaca_Cluster_6.txt", "r")
Alpaca_markers_6 = file6.readlines()
file6.close()

file6 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Human/Single_cluster_data/Alpaca_Cluster_7.txt", "r")
Alpaca_markers_7 = file6.readlines()
file6.close()

file6 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Human/Single_cluster_data/Alpaca_Cluster_8.txt", "r")
Alpaca_markers_8 = file6.readlines()
file6.close()

file6 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Human/Single_cluster_data/Alpaca_Cluster_9.txt", "r")
Alpaca_markers_9 = file6.readlines()
file6.close()

file6 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Human/Single_cluster_data/Alpaca_Cluster_10.txt", "r")
Alpaca_markers_10 = file6.readlines()
file6.close()
"""
# Import Files for the Human Clusters
file7 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Human/Single_cluster_data/Human_Cluster_0.txt", "r")
Human_markers_0 = file7.readlines()
file7.close()

file8 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Human/Single_cluster_data/Human_Cluster_1.txt", "r")
Human_markers_1 = file8.readlines()
file8.close()

file9 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Human/Single_cluster_data/Human_Cluster_2.txt", "r")
Human_markers_2 = file9.readlines()
file9.close()

file10 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Human/Single_cluster_data/Human_Cluster_3.txt", "r")
Human_markers_3 = file10.readlines()
file10.close()

file11 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Human/Single_cluster_data/Human_Cluster_4.txt", "r")
Human_markers_4 = file11.readlines()
file11.close()

file12 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Human/Single_cluster_data/Human_Cluster_5.txt", "r")
Human_markers_5 = file12.readlines()
file12.close()  

file12 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Human/Single_cluster_data/Human_Cluster_6.txt", "r")
Human_markers_6 = file12.readlines()
file12.close() 


      
################################################################################################################################
#compare the marker list of the alpaca cluster X to all the human clusters
#Cluster 0
def compare_clusters(Alpaca, Human):
    n = 0
    same_marker = []
    for gene in Alpaca[1:]:
        for gene1 in Human[1:]:
            if gene == gene1:
                if gene not in same_marker:
                    same_marker.append(gene)
                    n += 1
    same_count = 0
    for el in same_marker:
        same_count += 1
    return(same_marker, same_count)


#Cluster 0 Alpaca
Alpaca_0_Human_0 = compare_clusters(Alpaca_markers_0, Human_markers_0)
Alpaca_0_Human_1 = compare_clusters(Alpaca_markers_0, Human_markers_1)
Alpaca_0_Human_2 = compare_clusters(Alpaca_markers_0, Human_markers_2)
Alpaca_0_Human_3 = compare_clusters(Alpaca_markers_0, Human_markers_3)
Alpaca_0_Human_4 = compare_clusters(Alpaca_markers_0, Human_markers_4)
Alpaca_0_Human_5 = compare_clusters(Alpaca_markers_0, Human_markers_5)
Alpaca_0_Human_6 = compare_clusters(Alpaca_markers_0, Human_markers_6)


#Cluster 1 Alpaca
Alpaca_1_Human_0 = compare_clusters(Alpaca_markers_1, Human_markers_0)
Alpaca_1_Human_1 = compare_clusters(Alpaca_markers_1, Human_markers_1)
Alpaca_1_Human_2 = compare_clusters(Alpaca_markers_1, Human_markers_2)
Alpaca_1_Human_3 = compare_clusters(Alpaca_markers_1, Human_markers_3)
Alpaca_1_Human_4 = compare_clusters(Alpaca_markers_1, Human_markers_4)
Alpaca_1_Human_5 = compare_clusters(Alpaca_markers_1, Human_markers_5)
Alpaca_1_Human_6 = compare_clusters(Alpaca_markers_1, Human_markers_6)



#Cluster 2 Alpaca
Alpaca_2_Human_0 = compare_clusters(Alpaca_markers_2, Human_markers_0)
Alpaca_2_Human_1 = compare_clusters(Alpaca_markers_2, Human_markers_1)
Alpaca_2_Human_2 = compare_clusters(Alpaca_markers_2, Human_markers_2)
Alpaca_2_Human_3 = compare_clusters(Alpaca_markers_2, Human_markers_3)
Alpaca_2_Human_4 = compare_clusters(Alpaca_markers_2, Human_markers_4)
Alpaca_2_Human_5 = compare_clusters(Alpaca_markers_2, Human_markers_5)
Alpaca_2_Human_6 = compare_clusters(Alpaca_markers_2, Human_markers_6)



#Cluster 3 Alpaca
Alpaca_3_Human_0 = compare_clusters(Alpaca_markers_3, Human_markers_0)
Alpaca_3_Human_1 = compare_clusters(Alpaca_markers_3, Human_markers_1)
Alpaca_3_Human_2 = compare_clusters(Alpaca_markers_3, Human_markers_2)
Alpaca_3_Human_3 = compare_clusters(Alpaca_markers_3, Human_markers_3)
Alpaca_3_Human_4 = compare_clusters(Alpaca_markers_3, Human_markers_4)
Alpaca_3_Human_5 = compare_clusters(Alpaca_markers_3, Human_markers_5)
Alpaca_3_Human_6 = compare_clusters(Alpaca_markers_3, Human_markers_6)



#Cluster 4 Alpaca
Alpaca_4_Human_0 = compare_clusters(Alpaca_markers_4, Human_markers_0)
Alpaca_4_Human_1 = compare_clusters(Alpaca_markers_4, Human_markers_1)
Alpaca_4_Human_2 = compare_clusters(Alpaca_markers_4, Human_markers_2)
Alpaca_4_Human_3 = compare_clusters(Alpaca_markers_4, Human_markers_3)
Alpaca_4_Human_4 = compare_clusters(Alpaca_markers_4, Human_markers_4)
Alpaca_4_Human_5 = compare_clusters(Alpaca_markers_4, Human_markers_5)
Alpaca_4_Human_6 = compare_clusters(Alpaca_markers_4, Human_markers_6)



#Cluster 5 Alpaca
Alpaca_5_Human_0 = compare_clusters(Alpaca_markers_5, Human_markers_0)
Alpaca_5_Human_1 = compare_clusters(Alpaca_markers_5, Human_markers_1)
Alpaca_5_Human_2 = compare_clusters(Alpaca_markers_5, Human_markers_2)
Alpaca_5_Human_3 = compare_clusters(Alpaca_markers_5, Human_markers_3)
Alpaca_5_Human_4 = compare_clusters(Alpaca_markers_5, Human_markers_4)
Alpaca_5_Human_5 = compare_clusters(Alpaca_markers_5, Human_markers_5)
Alpaca_5_Human_6 = compare_clusters(Alpaca_markers_5, Human_markers_6)


"""
#Cluster 6 Alpaca
Alpaca_6_Human_0 = compare_clusters(Alpaca_markers_6, Human_markers_0)
Alpaca_6_Human_1 = compare_clusters(Alpaca_markers_6, Human_markers_1)
Alpaca_6_Human_2 = compare_clusters(Alpaca_markers_6, Human_markers_2)
Alpaca_6_Human_3 = compare_clusters(Alpaca_markers_6, Human_markers_3)
Alpaca_6_Human_4 = compare_clusters(Alpaca_markers_6, Human_markers_4)
Alpaca_6_Human_5 = compare_clusters(Alpaca_markers_6, Human_markers_5)
Alpaca_6_Human_6 = compare_clusters(Alpaca_markers_6, Human_markers_6)



#Cluster 7 Alpaca
Alpaca_7_Human_0 = compare_clusters(Alpaca_markers_7, Human_markers_0)
Alpaca_7_Human_1 = compare_clusters(Alpaca_markers_7, Human_markers_1)
Alpaca_7_Human_2 = compare_clusters(Alpaca_markers_7, Human_markers_2)
Alpaca_7_Human_3 = compare_clusters(Alpaca_markers_7, Human_markers_3)
Alpaca_7_Human_4 = compare_clusters(Alpaca_markers_7, Human_markers_4)
Alpaca_7_Human_5 = compare_clusters(Alpaca_markers_7, Human_markers_5)
Alpaca_7_Human_6 = compare_clusters(Alpaca_markers_7, Human_markers_6)


#Cluster 8 Alpaca
Alpaca_8_Human_0 = compare_clusters(Alpaca_markers_8, Human_markers_0)
Alpaca_8_Human_1 = compare_clusters(Alpaca_markers_8, Human_markers_1)
Alpaca_8_Human_2 = compare_clusters(Alpaca_markers_8, Human_markers_2)
Alpaca_8_Human_3 = compare_clusters(Alpaca_markers_8, Human_markers_3)
Alpaca_8_Human_4 = compare_clusters(Alpaca_markers_8, Human_markers_4)
Alpaca_8_Human_5 = compare_clusters(Alpaca_markers_8, Human_markers_5)
Alpaca_8_Human_6 = compare_clusters(Alpaca_markers_8, Human_markers_6)



#Cluster 9 Alpaca
Alpaca_9_Human_0 = compare_clusters(Alpaca_markers_9, Human_markers_0)
Alpaca_9_Human_1 = compare_clusters(Alpaca_markers_9, Human_markers_1)
Alpaca_9_Human_2 = compare_clusters(Alpaca_markers_9, Human_markers_2)
Alpaca_9_Human_3 = compare_clusters(Alpaca_markers_9, Human_markers_3)
Alpaca_9_Human_4 = compare_clusters(Alpaca_markers_9, Human_markers_4)
Alpaca_9_Human_5 = compare_clusters(Alpaca_markers_9, Human_markers_5)
Alpaca_9_Human_6 = compare_clusters(Alpaca_markers_9, Human_markers_6)


#Cluster 10 Alpaca
Alpaca_10_Human_0 = compare_clusters(Alpaca_markers_10, Human_markers_0)
Alpaca_10_Human_1 = compare_clusters(Alpaca_markers_10, Human_markers_1)
Alpaca_10_Human_2 = compare_clusters(Alpaca_markers_10, Human_markers_2)
Alpaca_10_Human_3 = compare_clusters(Alpaca_markers_10, Human_markers_3)
Alpaca_10_Human_4 = compare_clusters(Alpaca_markers_10, Human_markers_4)
Alpaca_10_Human_5 = compare_clusters(Alpaca_markers_10, Human_markers_5)
Alpaca_10_Human_6 = compare_clusters(Alpaca_markers_10, Human_markers_6)
"""
print(Alpaca_0_Human_0[0])


os.chdir("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Human/Alpaca_Human_Comparison") 
###############################################################################
#Cluster 0 Alpaca
f = open("Alpaca_0_Human_0_markers.txt", "w")
for el in Alpaca_0_Human_0[0]:
    f.write(el)
f.close()

f = open("Alpaca_0_Human_1_markers.txt", "w")
for el in Alpaca_0_Human_1[0]:
    f.write(el)
f.close()

f = open("Alpaca_0_Human_2_markers.txt", "w")
for el in Alpaca_0_Human_2[0]:
    f.write(el)
f.close()

f = open("Alpaca_0_Human_3_markers.txt", "w")
for el in Alpaca_0_Human_3[0]:
    f.write(el)
f.close()

f = open("Alpaca_0_Human_4_markers.txt", "w")
for el in Alpaca_0_Human_4[0]:
    f.write(el)
f.close()

f = open("Alpaca_0_Human_5_markers.txt", "w")
for el in Alpaca_0_Human_5[0]:
    f.write(el)
f.close()

f = open("Alpaca_0_Human_6_markers.txt", "w")
for el in Alpaca_0_Human_6[0]:
    f.write(el)
f.close()



###############################################################################
#Cluster 1 Alpaca
f = open("Alpaca_1_Human_0_markers.txt", "w")
for el in Alpaca_1_Human_0[0]:
    f.write(el)
f.close()

f = open("Alpaca_1_Human_1_markers.txt", "w")
for el in Alpaca_1_Human_1[0]:
    f.write(el)
f.close()

f = open("Alpaca_1_Human_2_markers.txt", "w")
for el in Alpaca_1_Human_2[0]:
    f.write(el)
f.close()

f = open("Alpaca_1_Human_3_markers.txt", "w")
for el in Alpaca_1_Human_3[0]:
    f.write(el)
f.close()

f = open("Alpaca_1_Human_4_markers.txt", "w")
for el in Alpaca_1_Human_4[0]:
    f.write(el)
f.close()

f = open("Alpaca_1_Human_5_markers.txt", "w")
for el in Alpaca_1_Human_5[0]:
    f.write(el)
f.close()

f = open("Alpaca_1_Human_6_markers.txt", "w")
for el in Alpaca_1_Human_6[0]:
    f.write(el)
f.close()



###############################################################################
#Cluster 2 Alpaca
f = open("Alpaca_2_Human_0_markers.txt", "w")
for el in Alpaca_2_Human_0[0]:
    f.write(el)
f.close()

f = open("Alpaca_2_Human_1_markers.txt", "w")
for el in Alpaca_2_Human_1[0]:
    f.write(el)
f.close()

f = open("Alpaca_2_Human_2_markers.txt", "w")
for el in Alpaca_2_Human_2[0]:
    f.write(el)
f.close()

f = open("Alpaca_2_Human_3_markers.txt", "w")
for el in Alpaca_2_Human_3[0]:
    f.write(el)
f.close()

f = open("Alpaca_2_Human_4_markers.txt", "w")
for el in Alpaca_2_Human_4[0]:
    f.write(el)
f.close()

f = open("Alpaca_2_Human_5_markers.txt", "w")
for el in Alpaca_2_Human_5[0]:
    f.write(el)
f.close()

f = open("Alpaca_2_Human_6_markers.txt", "w")
for el in Alpaca_2_Human_6[0]:
    f.write(el)
f.close()



###############################################################################
#Cluster 3 Alpaca
f = open("Alpaca_3_Human_0_markers.txt", "w")
for el in Alpaca_3_Human_0[0]:
    f.write(el)
f.close()

f = open("Alpaca_3_Human_1_markers.txt", "w")
for el in Alpaca_3_Human_1[0]:
    f.write(el)
f.close()

f = open("Alpaca_3_Human_2_markers.txt", "w")
for el in Alpaca_3_Human_2[0]:
    f.write(el)
f.close()

f = open("Alpaca_3_Human_3_markers.txt", "w")
for el in Alpaca_3_Human_3[0]:
    f.write(el)
f.close()

f = open("Alpaca_3_Human_4_markers.txt", "w")
for el in Alpaca_3_Human_4[0]:
    f.write(el)
f.close()

f = open("Alpaca_3_Human_5_markers.txt", "w")
for el in Alpaca_3_Human_5[0]:
    f.write(el)
f.close()

f = open("Alpaca_3_Human_6_markers.txt", "w")
for el in Alpaca_3_Human_6[0]:
    f.write(el)
f.close()



###############################################################################
#Cluster 4 Alpaca
f = open("Alpaca_4_Human_0_markers.txt", "w")
for el in Alpaca_4_Human_0[0]:
    f.write(el)
f.close()

f = open("Alpaca_4_Human_1_markers.txt", "w")
for el in Alpaca_4_Human_1[0]:
    f.write(el)
f.close()

f = open("Alpaca_4_Human_2_markers.txt", "w")
for el in Alpaca_4_Human_2[0]:
    f.write(el)
f.close()

f = open("Alpaca_4_Human_3_markers.txt", "w")
for el in Alpaca_4_Human_3[0]:
    f.write(el)
f.close()

f = open("Alpaca_4_Human_4_markers.txt", "w")
for el in Alpaca_4_Human_4[0]:
    f.write(el)
f.close()

f = open("Alpaca_4_Human_5_markers.txt", "w")
for el in Alpaca_4_Human_5[0]:
    f.write(el)
f.close()

f = open("Alpaca_4_Human_6_markers.txt", "w")
for el in Alpaca_4_Human_6[0]:
    f.write(el)
f.close()



###############################################################################
#Cluster 5 Alpaca
f = open("Alpaca_5_Human_0_markers.txt", "w")
for el in Alpaca_5_Human_0[0]:
    f.write(el)
f.close()

f = open("Alpaca_5_Human_1_markers.txt", "w")
for el in Alpaca_5_Human_1[0]:
    f.write(el)
f.close()

f = open("Alpaca_5_Human_2_markers.txt", "w")
for el in Alpaca_5_Human_2[0]:
    f.write(el)
f.close()

f = open("Alpaca_5_Human_3_markers.txt", "w")
for el in Alpaca_5_Human_3[0]:
    f.write(el)
f.close()

f = open("Alpaca_5_Human_4_markers.txt", "w")
for el in Alpaca_5_Human_4[0]:
    f.write(el)
f.close()

f = open("Alpaca_5_Human_5_markers.txt", "w")
for el in Alpaca_5_Human_5[0]:
    f.write(el)
f.close()

f = open("Alpaca_5_Human_6_markers.txt", "w")
for el in Alpaca_5_Human_6[0]:
    f.write(el)
f.close()


"""
###############################################################################
#Cluster 6 Alpaca
f = open("Alpaca_6_Human_0_markers.txt", "w")
for el in Alpaca_6_Human_0[0]:
    f.write(el)
f.close()

f = open("Alpaca_6_Human_1_markers.txt", "w")
for el in Alpaca_6_Human_1[0]:
    f.write(el)
f.close()

f = open("Alpaca_6_Human_2_markers.txt", "w")
for el in Alpaca_6_Human_2[0]:
    f.write(el)
f.close()

f = open("Alpaca_6_Human_3_markers.txt", "w")
for el in Alpaca_6_Human_3[0]:
    f.write(el)
f.close()

f = open("Alpaca_6_Human_4_markers.txt", "w")
for el in Alpaca_6_Human_4[0]:
    f.write(el)
f.close()

f = open("Alpaca_6_Human_5_markers.txt", "w")
for el in Alpaca_6_Human_5[0]:
    f.write(el)
f.close()

f = open("Alpaca_6_Human_6_markers.txt", "w")
for el in Alpaca_6_Human_6[0]:
    f.write(el)
f.close()



###############################################################################
#Cluster 7 Alpaca
f = open("Alpaca_7_Human_0_markers.txt", "w")
for el in Alpaca_7_Human_0[0]:
    f.write(el)
f.close()

f = open("Alpaca_7_Human_1_markers.txt", "w")
for el in Alpaca_7_Human_1[0]:
    f.write(el)
f.close()

f = open("Alpaca_7_Human_2_markers.txt", "w")
for el in Alpaca_7_Human_2[0]:
    f.write(el)
f.close()

f = open("Alpaca_7_Human_3_markers.txt", "w")
for el in Alpaca_7_Human_3[0]:
    f.write(el)
f.close()

f = open("Alpaca_7_Human_4_markers.txt", "w")
for el in Alpaca_7_Human_4[0]:
    f.write(el)
f.close()

f = open("Alpaca_7_Human_5_markers.txt", "w")
for el in Alpaca_7_Human_5[0]:
    f.write(el)
f.close()

f = open("Alpaca_7_Human_6_markers.txt", "w")
for el in Alpaca_7_Human_6[0]:
    f.write(el)
f.close()




###############################################################################
#Cluster 8 Alpaca
f = open("Alpaca_8_Human_0_markers.txt", "w")
for el in Alpaca_8_Human_0[0]:
    f.write(el)
f.close()

f = open("Alpaca_8_Human_1_markers.txt", "w")
for el in Alpaca_8_Human_1[0]:
    f.write(el)
f.close()

f = open("Alpaca_8_Human_2_markers.txt", "w")
for el in Alpaca_8_Human_2[0]:
    f.write(el)
f.close()

f = open("Alpaca_8_Human_3_markers.txt", "w")
for el in Alpaca_8_Human_3[0]:
    f.write(el)
f.close()

f = open("Alpaca_8_Human_4_markers.txt", "w")
for el in Alpaca_8_Human_4[0]:
    f.write(el)
f.close()

f = open("Alpaca_8_Human_5_markers.txt", "w")
for el in Alpaca_8_Human_5[0]:
    f.write(el)
f.close()

f = open("Alpaca_8_Human_6_markers.txt", "w")
for el in Alpaca_8_Human_6[0]:
    f.write(el)
f.close()


###############################################################################

#Cluster 9 Alpaca
f = open("Alpaca_9_Human_0_markers.txt", "w")
for el in Alpaca_9_Human_0[0]:
    f.write(el)
f.close()

f = open("Alpaca_9_Human_1_markers.txt", "w")
for el in Alpaca_9_Human_1[0]:
    f.write(el)
f.close()

f = open("Alpaca_9_Human_2_markers.txt", "w")
for el in Alpaca_9_Human_2[0]:
    f.write(el)
f.close()

f = open("Alpaca_9_Human_3_markers.txt", "w")
for el in Alpaca_9_Human_3[0]:
    f.write(el)
f.close()

f = open("Alpaca_9_Human_4_markers.txt", "w")
for el in Alpaca_9_Human_4[0]:
    f.write(el)
f.close()

f = open("Alpaca_9_Human_5_markers.txt", "w")
for el in Alpaca_9_Human_5[0]:
    f.write(el)
f.close()

f = open("Alpaca_9_Human_6_markers.txt", "w")
for el in Alpaca_9_Human_6[0]:
    f.write(el)
f.close()

###############################################################################

#Cluster 10 Alpaca
f = open("Alpaca_10_Human_0_markers.txt", "w")
for el in Alpaca_10_Human_0[0]:
    f.write(el)
f.close()

f = open("Alpaca_10_Human_1_markers.txt", "w")
for el in Alpaca_10_Human_1[0]:
    f.write(el)
f.close()

f = open("Alpaca_10_Human_2_markers.txt", "w")
for el in Alpaca_10_Human_2[0]:
    f.write(el)
f.close()

f = open("Alpaca_10_Human_3_markers.txt", "w")
for el in Alpaca_10_Human_3[0]:
    f.write(el)
f.close()

f = open("Alpaca_10_Human_4_markers.txt", "w")
for el in Alpaca_10_Human_4[0]:
    f.write(el)
f.close()

f = open("Alpaca_10_Human_5_markers.txt", "w")
for el in Alpaca_10_Human_5[0]:
    f.write(el)
f.close()

f = open("Alpaca_10_Human_6_markers.txt", "w")
for el in Alpaca_10_Human_6[0]:
    f.write(el)
f.close()
"""