# 2020-05-14

# Kiyoto Aramis Tanemura

# After assigning conformers to communities, choose centroid based on various criteria

import numpy as np
import pandas as pd

def centroid_weighted_degree(fileList, communityAssignment, affinityMatrix):
    # Provided with list of conformers assigned to communities, choose representative centroid by conformers of maximum in community weighted degree
    # inputs
    # fileList: (list) names of xyz files for each conformer
    # communityAssignment: (list) list of community assignment correspoinding to the index of fileList
    # affinityMatrix: {np.array) affinity matrix for each pairwise conformer similarity
    numFiles = len(fileList)
    communityList = list(set(communityAssignment))
    centralNodes = []
    comm_size = []
    for C in communityList:                                                                                                                                                                           
        C_members = [x for x in range(numFiles) if communityAssignment[x] == C]
        C_member_files = [fileList[x] for x in C_members]
        comm_size.append(len(C_members))
        community_subgraph = affinityMatrix[C_members, :][:, C_members]                                                                                                                               
        weightedDegrees = np.sum(community_subgraph, axis = 1).ravel().tolist()                                                                                                                       
        max_wdegree = np.max(weightedDegrees)                                                                                                                                                         
        max_index = weightedDegrees.index(max_wdegree) 
        centralNodes.append(C_member_files[max_index])

    # Sort centers by size of clusters in descending order
    centralDf = pd.DataFrame({'size': comm_size}, index = centralNodes)
    centralDf.sort_values(by = 'size', ascending = False, inplace = True)

    return list(centralDf.index)

def centroid_energy(fileList, communityAssignment, Epath, E_label = 'energy'):
    # Provided with list of conformers assigned to communities, choose representative centroid by lowest computed ANI energy
    # inputs                                                                                                                                                                                              
    # fileList: (list) names of xyz files for each conformer                                                                                                                                              
    # communityAssignment: (list) list of community assignment correspoinding to the index of fileList                                                                                                    
    # Epath: (string) path of energy values saved in csv
    # E_label: (string) column name for the energy values. 
    E_data = pd.read_csv(Epath, index_col = 0)
    comm_df = pd.DataFrame({'cluster': communityAssignment}, index = fileList)
    comm_df = comm_df.join(E_data, how = 'inner', sort = True) # inner join to avoid missing energy values
    communityList = list(set(communityAssignment))
    centralNodes = []
    for C in communityList:
        sub_df = comm_df[comm_df['cluster'] == C]
        min_val = sub_df[E_label].min()
        centralNodes.append(sub_df.loc[sub_df[E_label] == min_val].index[0])

    # Sort centers by energy in ascending order
    centralDf = comm_df.loc[centralNodes]
    centralDf.sort_values(by = E_label, inplace = True)

    return list(centralDf.index)

def centroid_medoid(fileList, communityAssignment, rmsdMatrix):
    # Choose centroids for each cluster by their medoid. 
    # inputs
    # fileList: (list) names of xyz or pdb files for each conformer
    # communityAssignment: (list) list of community assignment correspoinding to the index of fileList
    # rmsdMatrix: (numpy array) Matrix containing pairwise atomic RMSD between all conformers

    numFiles = len(fileList)
    communityList = list(set(communityAssignment))
    centralNodes = []
    comm_size = []
    for C in communityList:
        C_members = [x for x in range(numFiles) if communityAssignment[x] == C]
        C_member_files = [fileList[x] for x in C_members]
        comm_size.append(len(C_members))
        community_subgraph = rmsdMatrix[C_members, :][:, C_members]
        dist_sum = np.sum(community_subgraph, axis = 1)
        min_index = np.argmin(dist_sum)
        centralNodes.append(C_member_files[min_index])

    # Sort centers by size of clusters in descending order
    centralDf = pd.DataFrame({'size': comm_size}, index = centralNodes)
    centralDf.sort_values(by = 'size', ascending = False, inplace = True)

    return list(centralDf.index)


# added 2020-09-08 by KAT. My suspicion is that the selection by weighted degree is unsuitable for MD frames. 
# Will select centroids by the lowest in-cluster eccentricity. Eccentricity can be calculatd using Dijkstra's shortest path algorithm.

def diekstra(filtered_rmsd_matrix_community, i):
    '''Use Dijkstra's algorithm to find shortest path from index i to all other nodes'''
    # initialize lists
    visited = []
    unvisited = [x for x in range(filtered_rmsd_matrix_community.shape[0])]
    record = [np.inf for x in range(filtered_rmsd_matrix_community.shape[0])]#{x: np.inf for x in range(graph.shape[0])}
    record[i] = 0
    lastNode = [-1 for x in record]
    # repeat until all nodes have been visited
    while len(unvisited) > 0:
        visit_index = unvisited[np.argmin([record[x] for x in unvisited])]
        unvisited_neighbors = [x for x in unvisited if filtered_rmsd_matrix_community[visit_index, x] > 0]
        # Calculate distance to unvisited neighbor. If value is shorter than recorded, update distance.
        updateDist = filtered_rmsd_matrix_community[visit_index, :] + record[visit_index]
        for j in unvisited_neighbors:
            record[j] = np.min([updateDist[j], record[j]])
            if updateDist[j] < record[j]:
                lastNode[j] = visit_index
        # update visited/unvisited node list
        unvisited.remove(visit_index)
        visited.append(visit_index)
    return record, lastNode

def centroid_eccentricity(fileList, communityAssignment, filtered_rmsd_matrix):
    # Provided with list of conformers assigned to communities, choose representative centroid by conformers of minimum in community eccentricity                                                      
    # inputs                                                                                                                                                                                               
    # fileList: (list) names of xyz files for each conformer                                                                                                                                               
    # communityAssignment: (list) list of community assignment correspoinding to the index of fileList                                                                                                     
    # filtered_rmsd_matrix: {np.array) RMSD matrix between conformers, except assigning distances above threshold to zero   
    numFiles = len(fileList)
    communityList = list(set(communityAssignment))
    centralNodes = []
    comm_size = []
    for C in communityList:
        C_members = [x for x in range(numFiles) if communityAssignment[x] == C]
        C_member_files = [fileList[x] for x in C_members]
        comm_size.append(len(C_members))
        community_subgraph = filtered_rmsd_matrix[C_members, :][:, C_members]
        community_eccentricities = []
        for i in range(len(C_member_files)):
            community_eccentricities.append(np.max(diekstra(community_subgraph, i)[0]))
        min_eccentricity_index = np.argmin(community_eccentricities)
        centralNodes.append(C_member_files[min_eccentricity_index])

    # Sort centers by size of clusters in descending order                                                                                                                                                
    centralDf = pd.DataFrame({'size': comm_size}, index = centralNodes)
    centralDf.sort_values(by = 'size', ascending = False, inplace = True)

    return list(centralDf.index)

def centroid_betweenness(fileList, communityAssignment, filtered_rmsd_matrix):
    # Provided with list of conformers assigned to communities, choose representative centroid by conformers of maximum in community betweenness 
    # inputs
    # fileList: (list) names of xyz files for each conformer
    # communityAssignment: (list) list of community assignment correspoinding to the index of fileList
    # filtered_rmsd_matrix: {np.array) RMSD matrix between conformers, except assigning distances above threshold to zero
    numFiles = len(fileList)
    communityList = list(set(communityAssignment))
    centralNodes = []
    comm_size = []
    for C in communityList:
        C_members = [x for x in range(numFiles) if communityAssignment[x] == C]
        C_member_files = [fileList[x] for x in C_members]
        comm_size.append(len(C_members))
        community_subgraph = filtered_rmsd_matrix[C_members, :][:, C_members]
        community_betweenness = np.zeros(len(C_members))
        for i in range(len(C_member_files)):
            record, lastnode = diekstra(community_subgraph, i)
            for j in range(len(C_member_files)):
                previous_node = lastnode[j]
                while previous_node != -1:
                    community_betweenness[previous_node] += 1
                    previous_node = lastnode[previous_node]
        max_betweenness_index = np.argmax(community_betweenness)
        centralNodes.append(C_member_files[max_betweenness_index])

    # Sort centers by size of clusters in descending order         
    centralDf = pd.DataFrame({'size': comm_size}, index = centralNodes)
    centralDf.sort_values(by = 'size', ascending = False, inplace = True)

    return list(centralDf.index)
