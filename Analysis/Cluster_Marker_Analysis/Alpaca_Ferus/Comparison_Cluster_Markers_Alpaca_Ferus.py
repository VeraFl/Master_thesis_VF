# -*- coding: utf-8 -*-
"""
Created on Thu May 19 13:57:17 2022

@author: veraf
"""

#Import modules
import os
################################################################################################################################
# Import Files for the Alpaca Clusters
file1 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Ferus/Single_cluster_data/Alpaca_Cluster_0.txt", "r")
Alpaca_markers_0 = file1.readlines()
file1.close()

file2 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Ferus/Single_cluster_data/Alpaca_Cluster_1.txt", "r")
Alpaca_markers_1 = file2.readlines()
file2.close()

file3 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Ferus/Single_cluster_data/Alpaca_Cluster_2.txt", "r")
Alpaca_markers_2 = file3.readlines()
file3.close()

file4 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Ferus/Single_cluster_data/Alpaca_Cluster_3.txt", "r")
Alpaca_markers_3 = file4.readlines()
file4.close()

file5 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Ferus/Single_cluster_data/Alpaca_Cluster_4.txt", "r")
Alpaca_markers_4 = file5.readlines()
file5.close()

file6 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Ferus/Single_cluster_data/Alpaca_Cluster_5.txt", "r")
Alpaca_markers_5 = file6.readlines()
file6.close()



# Import Files for the Ferus Clusters
file7 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Ferus/Single_cluster_data/Ferus_Cluster_0.txt", "r")
Ferus_markers_0 = file7.readlines()
file7.close()

file8 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Ferus/Single_cluster_data/Ferus_Cluster_1.txt", "r")
Ferus_markers_1 = file8.readlines()
file8.close()

file9 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Ferus/Single_cluster_data/Ferus_Cluster_2.txt", "r")
Ferus_markers_2 = file9.readlines()
file9.close()

file10 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Ferus/Single_cluster_data/Ferus_Cluster_3.txt", "r")
Ferus_markers_3 = file10.readlines()
file10.close()

file11 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Ferus/Single_cluster_data/Ferus_Cluster_4.txt", "r")
Ferus_markers_4 = file11.readlines()
file11.close()

file12 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Ferus/Single_cluster_data/Ferus_Cluster_5.txt", "r")
Ferus_markers_5 = file12.readlines()
file12.close()  

"""
file12 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Ferus/Single_cluster_data/Ferus_Cluster_6.txt", "r")
Ferus_markers_6 = file12.readlines()
file12.close() 

file11 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Ferus/Single_cluster_data/Ferus_Cluster_7.txt", "r")
Ferus_markers_7 = file11.readlines()
file11.close()

file12 = open("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Ferus/Single_cluster_data/Ferus_Cluster_8.txt", "r")
Ferus_markers_8 = file12.readlines()
file12.close()  
"""
      
################################################################################################################################
#compare the marker list of the alpaca cluster X to all the Ferus clusters
#Cluster 0
def compare_clusters(Alpaca, Ferus):
    n = 0
    same_marker = []
    for gene in Alpaca[1:]:
        for gene1 in Ferus[1:]:
            if gene == gene1:
                if gene not in same_marker:
                    same_marker.append(gene)
                    n += 1
    same_count = 0
    for el in same_marker:
        same_count += 1
    return(same_marker, same_count)


#Cluster 0 Alpaca
Alpaca_0_Ferus_0 = compare_clusters(Alpaca_markers_0, Ferus_markers_0)
Alpaca_0_Ferus_1 = compare_clusters(Alpaca_markers_0, Ferus_markers_1)
Alpaca_0_Ferus_2 = compare_clusters(Alpaca_markers_0, Ferus_markers_2)
Alpaca_0_Ferus_3 = compare_clusters(Alpaca_markers_0, Ferus_markers_3)
Alpaca_0_Ferus_4 = compare_clusters(Alpaca_markers_0, Ferus_markers_4)
Alpaca_0_Ferus_5 = compare_clusters(Alpaca_markers_0, Ferus_markers_5)
#Alpaca_0_Ferus_6 = compare_clusters(Alpaca_markers_0, Ferus_markers_6)
#Alpaca_0_Ferus_7 = compare_clusters(Alpaca_markers_0, Ferus_markers_7)
#Alpaca_0_Ferus_8 = compare_clusters(Alpaca_markers_0, Ferus_markers_8)


#Cluster 1 Alpaca
Alpaca_1_Ferus_0 = compare_clusters(Alpaca_markers_1, Ferus_markers_0)
Alpaca_1_Ferus_1 = compare_clusters(Alpaca_markers_1, Ferus_markers_1)
Alpaca_1_Ferus_2 = compare_clusters(Alpaca_markers_1, Ferus_markers_2)
Alpaca_1_Ferus_3 = compare_clusters(Alpaca_markers_1, Ferus_markers_3)
Alpaca_1_Ferus_4 = compare_clusters(Alpaca_markers_1, Ferus_markers_4)
Alpaca_1_Ferus_5 = compare_clusters(Alpaca_markers_1, Ferus_markers_5)
#Alpaca_1_Ferus_6 = compare_clusters(Alpaca_markers_1, Ferus_markers_6)
#Alpaca_1_Ferus_7 = compare_clusters(Alpaca_markers_1, Ferus_markers_7)
#Alpaca_1_Ferus_8 = compare_clusters(Alpaca_markers_1, Ferus_markers_8)



#Cluster 2 Alpaca
Alpaca_2_Ferus_0 = compare_clusters(Alpaca_markers_2, Ferus_markers_0)
Alpaca_2_Ferus_1 = compare_clusters(Alpaca_markers_2, Ferus_markers_1)
Alpaca_2_Ferus_2 = compare_clusters(Alpaca_markers_2, Ferus_markers_2)
Alpaca_2_Ferus_3 = compare_clusters(Alpaca_markers_2, Ferus_markers_3)
Alpaca_2_Ferus_4 = compare_clusters(Alpaca_markers_2, Ferus_markers_4)
Alpaca_2_Ferus_5 = compare_clusters(Alpaca_markers_2, Ferus_markers_5)
#Alpaca_2_Ferus_6 = compare_clusters(Alpaca_markers_2, Ferus_markers_6)
#Alpaca_2_Ferus_7 = compare_clusters(Alpaca_markers_2, Ferus_markers_7)
#Alpaca_2_Ferus_8 = compare_clusters(Alpaca_markers_2, Ferus_markers_8)



#Cluster 3 Alpaca
Alpaca_3_Ferus_0 = compare_clusters(Alpaca_markers_3, Ferus_markers_0)
Alpaca_3_Ferus_1 = compare_clusters(Alpaca_markers_3, Ferus_markers_1)
Alpaca_3_Ferus_2 = compare_clusters(Alpaca_markers_3, Ferus_markers_2)
Alpaca_3_Ferus_3 = compare_clusters(Alpaca_markers_3, Ferus_markers_3)
Alpaca_3_Ferus_4 = compare_clusters(Alpaca_markers_3, Ferus_markers_4)
Alpaca_3_Ferus_5 = compare_clusters(Alpaca_markers_3, Ferus_markers_5)
#Alpaca_3_Ferus_6 = compare_clusters(Alpaca_markers_3, Ferus_markers_6)
#Alpaca_3_Ferus_7 = compare_clusters(Alpaca_markers_3, Ferus_markers_7)
#Alpaca_3_Ferus_8 = compare_clusters(Alpaca_markers_3, Ferus_markers_8)



#Cluster 4 Alpaca
Alpaca_4_Ferus_0 = compare_clusters(Alpaca_markers_4, Ferus_markers_0)
Alpaca_4_Ferus_1 = compare_clusters(Alpaca_markers_4, Ferus_markers_1)
Alpaca_4_Ferus_2 = compare_clusters(Alpaca_markers_4, Ferus_markers_2)
Alpaca_4_Ferus_3 = compare_clusters(Alpaca_markers_4, Ferus_markers_3)
Alpaca_4_Ferus_4 = compare_clusters(Alpaca_markers_4, Ferus_markers_4)
Alpaca_4_Ferus_5 = compare_clusters(Alpaca_markers_4, Ferus_markers_5)
#Alpaca_4_Ferus_6 = compare_clusters(Alpaca_markers_4, Ferus_markers_6)
#Alpaca_4_Ferus_7 = compare_clusters(Alpaca_markers_4, Ferus_markers_7)
#Alpaca_4_Ferus_8 = compare_clusters(Alpaca_markers_4, Ferus_markers_8)



#Cluster 5 Alpaca
Alpaca_5_Ferus_0 = compare_clusters(Alpaca_markers_5, Ferus_markers_0)
Alpaca_5_Ferus_1 = compare_clusters(Alpaca_markers_5, Ferus_markers_1)
Alpaca_5_Ferus_2 = compare_clusters(Alpaca_markers_5, Ferus_markers_2)
Alpaca_5_Ferus_3 = compare_clusters(Alpaca_markers_5, Ferus_markers_3)
Alpaca_5_Ferus_4 = compare_clusters(Alpaca_markers_5, Ferus_markers_4)
Alpaca_5_Ferus_5 = compare_clusters(Alpaca_markers_5, Ferus_markers_5)
#Alpaca_5_Ferus_6 = compare_clusters(Alpaca_markers_5, Ferus_markers_6)
#Alpaca_5_Ferus_7 = compare_clusters(Alpaca_markers_5, Ferus_markers_7)
#Alpaca_5_Ferus_8 = compare_clusters(Alpaca_markers_5, Ferus_markers_8)





print(Alpaca_0_Ferus_0[0])


os.chdir("C:/Users/veraf/OneDrive - Universitaet Bern/Single_Cell_CoV/Analysis/Cluster_Markers/Alpaca_Ferus/Alpaca_Ferus_Comparison") 
###############################################################################
#Cluster 0 Alpaca
f = open("Alpaca_0_Ferus_0_markers.txt", "w")
for el in Alpaca_0_Ferus_0[0]:
    f.write(el)
f.close()

f = open("Alpaca_0_Ferus_1_markers.txt", "w")
for el in Alpaca_0_Ferus_1[0]:
    f.write(el)
f.close()

f = open("Alpaca_0_Ferus_2_markers.txt", "w")
for el in Alpaca_0_Ferus_2[0]:
    f.write(el)
f.close()

f = open("Alpaca_0_Ferus_3_markers.txt", "w")
for el in Alpaca_0_Ferus_3[0]:
    f.write(el)
f.close()

f = open("Alpaca_0_Ferus_4_markers.txt", "w")
for el in Alpaca_0_Ferus_4[0]:
    f.write(el)
f.close()

f = open("Alpaca_0_Ferus_5_markers.txt", "w")
for el in Alpaca_0_Ferus_5[0]:
    f.write(el)
f.close()

"""
f = open("Alpaca_0_Ferus_6_markers.txt", "w")
for el in Alpaca_0_Ferus_6[0]:
    f.write(el)
f.close()

f = open("Alpaca_0_Ferus_7_markers.txt", "w")
for el in Alpaca_0_Ferus_7[0]:
    f.write(el)
f.close()

f = open("Alpaca_0_Ferus_8_markers.txt", "w")
for el in Alpaca_0_Ferus_8[0]:
    f.write(el)
f.close()

"""

###############################################################################
#Cluster 1 Alpaca
f = open("Alpaca_1_Ferus_0_markers.txt", "w")
for el in Alpaca_1_Ferus_0[0]:
    f.write(el)
f.close()

f = open("Alpaca_1_Ferus_1_markers.txt", "w")
for el in Alpaca_1_Ferus_1[0]:
    f.write(el)
f.close()

f = open("Alpaca_1_Ferus_2_markers.txt", "w")
for el in Alpaca_1_Ferus_2[0]:
    f.write(el)
f.close()

f = open("Alpaca_1_Ferus_3_markers.txt", "w")
for el in Alpaca_1_Ferus_3[0]:
    f.write(el)
f.close()

f = open("Alpaca_1_Ferus_4_markers.txt", "w")
for el in Alpaca_1_Ferus_4[0]:
    f.write(el)
f.close()

f = open("Alpaca_1_Ferus_5_markers.txt", "w")
for el in Alpaca_1_Ferus_5[0]:
    f.write(el)
f.close()

"""
f = open("Alpaca_1_Ferus_6_markers.txt", "w")
for el in Alpaca_1_Ferus_6[0]:
    f.write(el)
f.close()

f = open("Alpaca_1_Ferus_7_markers.txt", "w")
for el in Alpaca_1_Ferus_7[0]:
    f.write(el)
f.close()

f = open("Alpaca_1_Ferus_8_markers.txt", "w")
for el in Alpaca_1_Ferus_8[0]:
    f.write(el)
f.close()
"""


###############################################################################
#Cluster 2 Alpaca
f = open("Alpaca_2_Ferus_0_markers.txt", "w")
for el in Alpaca_2_Ferus_0[0]:
    f.write(el)
f.close()

f = open("Alpaca_2_Ferus_1_markers.txt", "w")
for el in Alpaca_2_Ferus_1[0]:
    f.write(el)
f.close()

f = open("Alpaca_2_Ferus_2_markers.txt", "w")
for el in Alpaca_2_Ferus_2[0]:
    f.write(el)
f.close()

f = open("Alpaca_2_Ferus_3_markers.txt", "w")
for el in Alpaca_2_Ferus_3[0]:
    f.write(el)
f.close()

f = open("Alpaca_2_Ferus_4_markers.txt", "w")
for el in Alpaca_2_Ferus_4[0]:
    f.write(el)
f.close()

f = open("Alpaca_2_Ferus_5_markers.txt", "w")
for el in Alpaca_2_Ferus_5[0]:
    f.write(el)
f.close()

"""
f = open("Alpaca_2_Ferus_6_markers.txt", "w")
for el in Alpaca_2_Ferus_6[0]:
    f.write(el)
f.close()

f = open("Alpaca_2_Ferus_7_markers.txt", "w")
for el in Alpaca_2_Ferus_7[0]:
    f.write(el)
f.close()

f = open("Alpaca_2_Ferus_8_markers.txt", "w")
for el in Alpaca_2_Ferus_8[0]:
    f.write(el)
f.close()
"""



###############################################################################
#Cluster 3 Alpaca
f = open("Alpaca_3_Ferus_0_markers.txt", "w")
for el in Alpaca_3_Ferus_0[0]:
    f.write(el)
f.close()

f = open("Alpaca_3_Ferus_1_markers.txt", "w")
for el in Alpaca_3_Ferus_1[0]:
    f.write(el)
f.close()

f = open("Alpaca_3_Ferus_2_markers.txt", "w")
for el in Alpaca_3_Ferus_2[0]:
    f.write(el)
f.close()

f = open("Alpaca_3_Ferus_3_markers.txt", "w")
for el in Alpaca_3_Ferus_3[0]:
    f.write(el)
f.close()

f = open("Alpaca_3_Ferus_4_markers.txt", "w")
for el in Alpaca_3_Ferus_4[0]:
    f.write(el)
f.close()

f = open("Alpaca_3_Ferus_5_markers.txt", "w")
for el in Alpaca_3_Ferus_5[0]:
    f.write(el)
f.close()


"""
f = open("Alpaca_3_Ferus_6_markers.txt", "w")
for el in Alpaca_3_Ferus_6[0]:
    f.write(el)
f.close()

f = open("Alpaca_3_Ferus_7_markers.txt", "w")
for el in Alpaca_3_Ferus_7[0]:
    f.write(el)
f.close()

f = open("Alpaca_3_Ferus_8_markers.txt", "w")
for el in Alpaca_3_Ferus_8[0]:
    f.write(el)
f.close()
"""

###############################################################################
#Cluster 4 Alpaca
f = open("Alpaca_4_Ferus_0_markers.txt", "w")
for el in Alpaca_4_Ferus_0[0]:
    f.write(el)
f.close()

f = open("Alpaca_4_Ferus_1_markers.txt", "w")
for el in Alpaca_4_Ferus_1[0]:
    f.write(el)
f.close()

f = open("Alpaca_4_Ferus_2_markers.txt", "w")
for el in Alpaca_4_Ferus_2[0]:
    f.write(el)
f.close()

f = open("Alpaca_4_Ferus_3_markers.txt", "w")
for el in Alpaca_4_Ferus_3[0]:
    f.write(el)
f.close()

f = open("Alpaca_4_Ferus_4_markers.txt", "w")
for el in Alpaca_4_Ferus_4[0]:
    f.write(el)
f.close()

f = open("Alpaca_4_Ferus_5_markers.txt", "w")
for el in Alpaca_4_Ferus_5[0]:
    f.write(el)
f.close()

"""
f = open("Alpaca_4_Ferus_6_markers.txt", "w")
for el in Alpaca_4_Ferus_6[0]:
    f.write(el)
f.close()

f = open("Alpaca_4_Ferus_7_markers.txt", "w")
for el in Alpaca_4_Ferus_7[0]:
    f.write(el)
f.close()

f = open("Alpaca_4_Ferus_8_markers.txt", "w")
for el in Alpaca_4_Ferus_8[0]:
    f.write(el)
f.close()

"""

###############################################################################
#Cluster 5 Alpaca
f = open("Alpaca_5_Ferus_0_markers.txt", "w")
for el in Alpaca_5_Ferus_0[0]:
    f.write(el)
f.close()

f = open("Alpaca_5_Ferus_1_markers.txt", "w")
for el in Alpaca_5_Ferus_1[0]:
    f.write(el)
f.close()

f = open("Alpaca_5_Ferus_2_markers.txt", "w")
for el in Alpaca_5_Ferus_2[0]:
    f.write(el)
f.close()

f = open("Alpaca_5_Ferus_3_markers.txt", "w")
for el in Alpaca_5_Ferus_3[0]:
    f.write(el)
f.close()

f = open("Alpaca_5_Ferus_4_markers.txt", "w")
for el in Alpaca_5_Ferus_4[0]:
    f.write(el)
f.close()

f = open("Alpaca_5_Ferus_5_markers.txt", "w")
for el in Alpaca_5_Ferus_5[0]:
    f.write(el)
f.close()

"""
f = open("Alpaca_5_Ferus_6_markers.txt", "w")
for el in Alpaca_5_Ferus_6[0]:
    f.write(el)
f.close()

f = open("Alpaca_5_Ferus_7_markers.txt", "w")
for el in Alpaca_5_Ferus_7[0]:
    f.write(el)
f.close()

f = open("Alpaca_5_Ferus_8_markers.txt", "w")
for el in Alpaca_5_Ferus_8[0]:
    f.write(el)
f.close()
"""

