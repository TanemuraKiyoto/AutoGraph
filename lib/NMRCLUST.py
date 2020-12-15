# 2020-04-29

# Kiyoto Aramis Tanemura

# NMRCLUST algorithm implemented in python using numpy. Original algorithm DOI: 10.1093/protein/9.11.1063
# To use in your code, first compute rmsdMatrix
# Then 'communityAssignment = NRMCLUST(rmsdMatrix)'
# The file list corresponding to the axis of the rmsdMatrix will be assigned to clusters specified in communityAssignment

import numpy as np

def averageLinkage(rmsdMatrix, clusterXmembers, clusterYmembers):
    subgraph = rmsdMatrix[clusterXmembers,:][:,clusterYmembers]
    return np.mean(subgraph)

def spread(rmsdMatrix, members):
    subgraph = rmsdMatrix[members,:][:,members]
    N = len(members)
    offDiagonalSum = np.sum(subgraph) / 2 # note: diagonals are zero
    return offDiagonalSum / (N * (N-1) / 2)

def memberIndices(communityID, communityAssignment):
    indices = [x for x in range(len(communityAssignment)) if communityAssignment[x] == communityID]
    return indices

def averageSpread(rmsdMatrix, communityAssignment):
    communities = list(set(communityAssignment))
    spreads = [spread(rmsdMatrix, memberIndices(C, communityAssignment)) for C in communities]
    return np.mean(spreads)

def mergeClusters(rmsdMatrix, communityAssignment, aveLinkArray = None):
    communities = list(set(communityAssignment))
    numCommunities = len(communities)
    if type(aveLinkArray) == type(None):
        aveLinkArray = np.zeros([numCommunities, numCommunities])
        for i in range(numCommunities - 1):
            for j in range(i + 1, numCommunities):
                i_members = memberIndices(communities[i], communityAssignment)
                j_members = memberIndices(communities[j], communityAssignment)
                aveLinkVal = averageLinkage(rmsdMatrix, i_members, j_members)
                aveLinkArray[i,j] = aveLinkVal
                aveLinkArray[j,i] = aveLinkVal

    # Populate diagonals (self similarity) with max average linkage value to remove it from consideration.
    aveLinkArray[range(numCommunities), range(numCommunities)] = np.max(aveLinkArray) + 0.01
    minVal = np.min(aveLinkArray)
    # Find value of minimum average linkage. Keep only the first index
    i, j = np.where(aveLinkArray == minVal)
    i = i[0]
    j = j[0]
    # Merge communities as recorded on communityAssignment
    C = communities[i]
    G = communities[j]
    i_members=memberIndices(communities[i], communityAssignment)
    j_members=memberIndices(communities[j], communityAssignment)
    for mem in j_members:
        communityAssignment[mem] = C

    for k in range(numCommunities):
        if communities[k] in [C, G]:
            continue
        k_members = memberIndices(communities[k], communityAssignment)
        aveLinkVal = averageLinkage(rmsdMatrix, i_members+j_members, k_members)
        aveLinkArray[i, k] = aveLinkVal
        aveLinkArray[k, i] = aveLinkVal

    aveLinkArray = np.delete(aveLinkArray, j, 0)
    aveLinkArray = np.delete(aveLinkArray, j, 1)

    return communityAssignment, aveLinkArray

def normalizeAvSpVal(AvSpVal, AvSpMax, AvSpMin, N):
    return (N - 1) / (AvSpMax - AvSpMin) * (AvSpVal - AvSpMin) + 1

def normalizeAvSp(AvSpList, N):
    AvSpMax = np.max(AvSpList)
    AvSpMin = np.min(AvSpList)
    if AvSpMax == AvSpMin:
        return [1 for x in AvSpList]
    return [normalizeAvSpVal(x, AvSpMax, AvSpMin, N) for x in AvSpList]

def NMRCLUST_(rmsdMatrix):
    N = rmsdMatrix.shape[0]
    AvSpList = []
    assignList = []
    aveLinkArray = None
    communityAssignment = list(range(N))
    singletonPresent = True # spread cannot be calculated if size of cluster is 1. Avoid spread calculation until each cluster has at least 2 members
    while singletonPresent:
        communityAssignment, aveLinkArray = mergeClusters(rmsdMatrix, communityAssignment, aveLinkArray)
        communities = list(set(communityAssignment))
        commcount = [len([x for x in communityAssignment if x == C]) for C in communities]
        if 1 not in commcount:
            singletonPresent = False
        elif len(communities) == 2:
            singletonPresent = False
    
    # Begin recording Average spread and cumminities assignment once singletons are absent
    # Continue recording until all are merged into one cluster
    while len(set(communityAssignment)) > 1:
        communityAssignment, aveLinkArray = mergeClusters(rmsdMatrix, communityAssignment, aveLinkArray)
        AvSpList.append(float(averageSpread(rmsdMatrix, communityAssignment)))
        assignList.append(list(communityAssignment))

    nClustList = [len(set(x)) for x in assignList] # number of clusters from when it was began recoded to one cluster, decreasing by one cluster at each step
    AvSpNormList = normalizeAvSp(AvSpList, N)
    penaltyVals = [AvSpNormList[x] + nClustList[x] for x in range(len(AvSpNormList))]
    minPenalty = np.min(penaltyVals)
    minPenaltyIndex = penaltyVals.index(minPenalty)

    return assignList[minPenaltyIndex]
