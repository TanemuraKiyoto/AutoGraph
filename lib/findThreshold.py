# 2020-04-17

# Kiyoto Aramis Tanemura

# Find threshold weight such that it is the maximum value while maintaining exactly one component graph
# input: affinity matrix as numpy array
# output: threshold value as float

import numpy as np
from .functions import rbfKernel

def findThreshold(affinityMatrix):
    numFiles = affinityMatrix.shape[0]
    # Get all values in affinityMatrix. Sort in descending order
    vals = np.array([])
    for i in range(numFiles - 1):
        vals = np.append(vals, affinityMatrix[i, i + 1:])

    vals = np.sort(vals, axis = None)[::-1]

    # Initialize index values
    upperIndex = 0
    lowerIndex = len(vals) - 1

    while upperIndex != lowerIndex - 1: # Until the indices are consecutive, iterate
        midIndex = int(np.mean([lowerIndex, upperIndex]))
        adjacencyMatrix = affinityMatrix > vals[midIndex]
        # Will tally nodes visited from first node in graph by BFS. If not all nodes were visited, there are more than one component.
        nodesVisited = [0] + [i for i in range(numFiles) if adjacencyMatrix[0, i] == 1]
#        nodesVisited = [i for i in range(numFiles) if adjacencyMatrix[0, i] == 1]
        for i in nodesVisited:
            newNodes = [j for j in range(numFiles) if adjacencyMatrix[i, j] == 1 and j not in nodesVisited]
            nodesVisited += newNodes
        while len(newNodes) > 0:
            newNodes1list = []
            for i in newNodes:
                newNodes1 = [j for j in range(numFiles) if adjacencyMatrix[i, j] == 1 and j not in nodesVisited]
                newNodes1list += newNodes1
                nodesVisited += newNodes1
        
            newNodes = newNodes1list

        # If every node was visited, then we have exactly one component. Otherwise, there are disconnected graphs
        if len(nodesVisited) < numFiles:
            upperIndex = midIndex
        else:
            lowerIndex = midIndex

    threshold = vals[lowerIndex]

    return threshold
