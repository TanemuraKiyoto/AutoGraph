# 2020-04-15

# Kiyoto Aramis Tanemura

# supporting functions for clustering by RMSD.
# I followed a numpy implementation of Kabsch algorithm to superimpose two coordinates (https://en.wikipedia.org/wiki/Kabsch_algorithm)

import numpy as np
from .functions import fileToArray

def getCentroid(coordMatrix):
    # return xyz coordinates of mean coordinate of provided coordinates
    # input: xyz coordinates of molecule, atoms in rows and x, y, z for columns
    return np.mean(coordMatrix, axis = 0)

def getDistance(coord0, coord1):
    # calculate Cartesian distance between two points in R3
    # inputs: two xyz coordinates as lists
    dif = coord0 - coord1
    sqdif = dif**2
    sumSqDif = np.sum(sqdif)
    return np.sqrt(sumSqDif)

def getCovarianceMatrix(coord0, coord1):
    return coord0.transpose() @ coord1

def getRotationMatrix(covMatrix):
    v, s, wt = np.linalg.svd(covMatrix)
    is_reflection = (np.linalg.det(v) * np.linalg.det(wt)) < 0
    if is_reflection:
        v[:, -1] = -v[:, -1]
    R = v @ wt
    return R

def superimpose(coord0, coord1, return_rot_matrix = False):
    coord0 -= getCentroid(coord0)
    coord1 -= getCentroid(coord1)
    H = getCovarianceMatrix(coord1, coord0)
    R = getRotationMatrix(H)
    if return_rot_matrix:
        return coord0, coord1 @ R, R
    return coord0, coord1 @ R

def rmsd(coord0, coord1):
    # input coordinate matrix of two structures
    coord0, coord1 = superimpose(coord0, coord1)
    num_atoms = coord0.shape[0]
    distances = []
    for i in range(num_atoms):
        distances.append(getDistance(coord0[i], coord1[i]))
    distSq = [x ** 2 for x in distances]
    return np.sqrt(np.mean(distSq))

def rmsdMatrix(fileList, inpath, hetatm = True):
    # Return numpy array of atomic RMSD values between conformers of each file.
    # inputs
    # fileList: (list of string) list of file names considered
    # inpath: (string) path to the input files
    # hetatm: (Boolean) if input files are pdb format, whether or not to read HETATM coordinates along with ATOM
    numFiles = len(fileList)
    coords = [fileToArray(inpath + theFile, hetatm) for theFile in fileList]
    rmsdMatrix = np.zeros([numFiles, numFiles])
    for i in range(numFiles):
        for j in range(i + 1, numFiles):
            rmsdVal = rmsd(coords[i], coords[j])
            rmsdMatrix[i, j] = rmsdVal
            rmsdMatrix[j, i] = rmsdVal

    return rmsdMatrix

def assign_remaining_files(master_fileList, fileList, centroids, communityAssignment, inpath):
    # Function for clustering with subset. Once clustering is finished on subset of files,                                                                                                                 
    # then assign remaining files to clusters by their proximity to centroids.                                                                                                                             
    # inputs                                                                                                                                                                                               
    # master_fileList: (list of string) list of all file names                                                                                                                                             
    # fileList: (list of string) list of file names sampled from master_fileList                                                                                                                           
    # centroids: (list of string) list of file names selected as community centroids                                                                                                                       
    # communityAssignment: (list of int) list of assinged community corresponding to fileList by index                                                                                                     
    # output: master community assignment. Similar to communityAssignment, but corresponding to index of master_fileList                                                                                   

    centroid_coords = [fileToArray(inpath + theFile) for theFile in centroids]
    num_centroids = len(centroids)
    centroid_communities = [communityAssignment[fileList.index(x)] for x in centroids]

    # remove already sampled files from yet to be assigned files.                                                                                                                                         
    remainder_indices = [x for x in range(len(master_fileList)) if master_fileList[x] not in fileList]
    remainder_coords = [fileToArray(inpath + master_fileList[x]) for x in remainder_indices]

    master_communityAssignment = [-1 for x in range(len(master_fileList))]

    # Fill master_communityAssignment first with already assigned conformers.                                                                                                                             
    for i in range(len(fileList)):
        sample_index = master_fileList.index(fileList[i])
        master_communityAssignment[sample_index] = communityAssignment[i]

    # Now assign the remainders.                                                                                                                                                                          
    for i in range(len(remainder_indices)):
        # get rmsd values of each remainder against all centroids                                                                                                                                         
        rmsdVals = [rmsd(remainder_coords[i], centroid_coords[j]) for j in range(num_centroids)]
        # find centroid corresponding to lowest RMSD                                                                                                                                                      
        minVal = np.min(rmsdVals)
        minIndex = rmsdVals.index(minVal)
        master_communityAssignment[remainder_indices[i]] = centroid_communities[minIndex]

    return master_communityAssignment
