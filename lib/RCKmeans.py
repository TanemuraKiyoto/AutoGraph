# 2020-05-20

# Kiyoto Aramis Tanemura

# We consider the use of representative conformation K-means to benchmark against AutoGraph as a method which does not require specification of number of clusters or threshold.
# Original algorithm found at DOI: 10.1186/s13321-017-0208-0

import numpy as np
from random import sample
from math import factorial

def kmedoid(rmsdMatrix, k = 2):
    # Generic k-medoid function with rmsd matrix as input
    n = rmsdMatrix.shape[0]
    medoids = sample(range(n), k)
    prev_medoids = []
    classification = np.zeros(n, dtype = int).tolist()
    while medoids != prev_medoids:
        for i in range(n):
            prev_medoids = medoids
            min_index = np.argmin(rmsdMatrix[i,:][medoids])
            classification[i] = min_index
        for j in range(k):
            members = [x for x in range(n) if classification[x] == j]
            sub_rmsd = rmsdMatrix[members,:][:, members]
            center_index = np.argmin(np.sum(sub_rmsd, axis = 0))
            medoids[j] = center_index

    return classification, medoids

def comb(n, r):
    if n < r:
        return 1
    return factorial(n) / (factorial(r) * factorial(n-r))

def MSQb(rmsdMatrix, medoids):
    sub_rmsd = rmsdMatrix[medoids, :][:, medoids]
    return np.sum(sub_rmsd) / (2 * comb(len(medoids), 2))

def MSQw(rmsdMatrix,classification):
    tally = 0
    n = len(classification)
    for i in set(classification):
        members = [x for x in range(n) if classification[x] == i]
        sub_rmsd = rmsdMatrix[members, :][:, members]
        tally += np.sum(sub_rmsd) / (2 * comb(len(members), 2))
    return tally / len(set(classification))

def SMA(MSQb_list, W = 10):
    if len(MSQb_list) >= W:
        return np.mean(MSQb_list[-10:])
    return -1

def RCKmeans_(rmsdMatrix):
    m = rmsdMatrix.shape[0]
    K_MSQb = [0, 0]
    prevSMA = -1
    for k in range(2,m):
        MSQw_list = np.zeros(100, dtype = int).tolist()
        MSQb_list = np.zeros(100, dtype = int).tolist()
        for i in range(100):
            classification, medoids = kmedoid(rmsdMatrix, k)
            MSQw_list[i] = MSQw(rmsdMatrix, classification)
            MSQb_list[i] = MSQb(rmsdMatrix, medoids)
        K_MSQb.append(MSQb_list[np.argmin(MSQw_list)])
        currSMA = SMA(K_MSQb, 10)
        if currSMA < prevSMA:
            return kmedoid(rmsdMatrix, np.argmax(K_MSQb))
        prevSMA = currSMA
