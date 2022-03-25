# 2020-04-17

# Kiyoto Aramis Tanemura

# Code for Louvain algorithm was getting lengthy. I dedicate its own module

import numpy as np

def getModularity(affinityMatrix, communityAssignments, resolution = 1.0):
    # Provided an affinity matrix (np.array) and nodes specifying their assigned communities (list),
    # return the modularity of the whole graph
    arrayDim = affinityMatrix.shape[0]
    communities = list(set(communityAssignments))
    Q = 0
    affinityMatrix[range(arrayDim), range(arrayDim)] += affinityMatrix[range(arrayDim), range(arrayDim)]
    two_m = np.sum(affinityMatrix)
    for C in communities:
        communityIndices = [i for i in range(len(communityAssignments)) if communityAssignments[i] == C]
        sigma_in = np.sum(affinityMatrix[communityIndices, :][:, communityIndices])
        sigma_tot = np.sum(affinityMatrix[communityIndices, :])
        Q += sigma_in / two_m - resolution * (sigma_tot / two_m) ** 2
    return Q

def LouvainPhase1(affinityMatrix, communityAssignment, Q_threshold, max_iter, resolution):
    # inputs: affinityMatrix as numpy array. Clear diagonals prior to entering and filter edges with weights below threshold
    # communityAssignments: list of length of nodes. Values are community IDs
    # Q_threshold: minimum global modularity to terminate phase 1
    # max_iter: number of iteration to force termination of phase 1
    # output: updated community Assignment list

    num_nodes = affinityMatrix.shape[0]

    affinityMatrix[range(num_nodes), range(num_nodes)] += affinityMatrix[range(num_nodes), range(num_nodes)]

    two_m = np.sum(affinityMatrix)

    modularity_current = getModularity(affinityMatrix, communityAssignment, resolution)
    changeModularity = 1

    adjacencyMatrix = affinityMatrix > 0
    iterations = 0

    while changeModularity > Q_threshold and iterations < max_iter:
        modularity_prev = modularity_current
        for i in range(num_nodes): # for each node
            neighborIndices = [j for j in range(num_nodes) if adjacencyMatrix[i,j] == 1]
            communities = set([communityAssignment[k] for k in neighborIndices])
            communities.discard(communityAssignment[i])
            communities = list(communities)
            modularityList = []
            # ki: sum of edges incident to node i 
            ki = np.sum(affinityMatrix[i, :])
            for C in communities:
                C_members = [x for x in range(num_nodes) if communityAssignment[x] == C]
                # Sum weights of all edges incident to nodes in C
                sigma_tot = np.sum(affinityMatrix[C_members, :])
                # ki_in: sum of weight of edges incident to node i in C
                ki_in = np.sum(affinityMatrix[i, :][C_members])
                deltaQ = ki_in/two_m - resolution * sigma_tot * ki / (two_m ** 2 / 2) # simplified formula. Derivation from (https://hal.archives-ouvertes.fr/hal-01231784/document)
                modularityList.append(deltaQ)
            modularityList.append(0) # append zero to avoid error by empty list
            maxQgain = np.max(modularityList)
            if maxQgain > 0:
                communityToJoin = communities[modularityList.index(maxQgain)]
                communityAssignment[i] = communityToJoin
        modularity_current = getModularity(affinityMatrix, communityAssignment, resolution)
        changeModularity = modularity_current - modularity_prev
        iterations += 1

    return communityAssignment

def LouvainPhase2(affinityMatrix, communityAssignment):
    # Merge each communities into supernodes
    # input: affinityMatrix or condensed graph
    # output: condensed graph and list of communities corresponding to the axis of graph
    num_nodes = affinityMatrix.shape[0]
    communities = list(set(communityAssignment))
    num_communities = len(communities)
    phase2graph = np.zeros([num_communities, num_communities])
    for i in range(num_communities):
        i_members = [x for x in range(num_nodes) if communityAssignment[x] == communities[i]]
        for j in range(num_communities):
            j_members = [x for x in range(num_nodes) if communityAssignment[x] == communities[j]]
            phase2graph[i, j] = np.sum(affinityMatrix[i_members, :][:, j_members]) * (1/2) ** (i == j)
            # if same communities, edges are double counted. If i == j, divide by two

    return phase2graph, communities

def Louvain(affinityMatrix, Q_threshold = 0.001, max_iter = 50, resolution = 1.0):
    # perform the two phases of Louvain community iteratively
    graph = affinityMatrix
    comm = list(range(affinityMatrix.shape[0]))
    # communityAssignmentRecord and commRefList are lists of lists, with same lengths.
    # communityAssignmentRecord stores assingment after phase 1
    # commRefList keeps the individual communities before assignment
    communityAssignmentRecord = []
    commRefList = []
    changeModularity = 1
    iterations = 0
    while changeModularity > Q_threshold and iterations < max_iter:
        graph, comm = LouvainPhase2(graph, comm)
        commRefList.append(list(comm))
        modularity_past = getModularity(graph, comm, resolution) # Note: we want to compare Q before and after reassignments in phase 1
        comm = LouvainPhase1(graph, comm, Q_threshold, max_iter, resolution)
        communityAssignmentRecord.append(list(comm))
        modularity_curr = getModularity(graph, comm, resolution)
        changeModularity = modularity_curr - modularity_past
        iterations += 1
#        print('changeModularity', changeModularity)

    for i in range(len(communityAssignmentRecord) - 2, 0, -1):
        changeDict = {}
        for j in range(len(communityAssignmentRecord[i])):
            newComm = communityAssignmentRecord[i][j]
            oldComm = commRefList[i][j]
            if newComm != oldComm:
                changeDict[oldComm] = newComm # Change key to value
        
        for key in changeDict:
            oldIndices = [x for x in range(len(communityAssignmentRecord[i - 1])) if communityAssignmentRecord[i - 1][x] == key]
            for j in oldIndices:
                communityAssignmentRecord[i - 1][j] = changeDict[key]
    return communityAssignmentRecord[0]
