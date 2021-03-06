# 2020-04-15

# Kiyoto Aramis Tanemura

# supporting functions for clustering by RMSD.
# I followed a numpy implementation of Kabsch algorithm to superimpose two coordinates (https://en.wikipedia.org/wiki/Kabsch_algorithm)

import numpy as np
import pandas as pd
import os

def xyzToArray(xyzFilePath, keepH = False):
    # Read xyz coordinates into numpy array. Remove hydrogen atoms
    with open(xyzFilePath, 'r') as f:
        lines = f.readlines()
    lines = [x.split() for x in lines]
    lines = [x for x in lines if len(x) == 4]
    if keepH:
        coords = [x[1:4] for x in lines]
    else:
        coords = [x[1:4] for x in lines if x[0] != 'H']
    for i in range(len(coords)):
        coords[i] = [float(x) for x in coords[i]]
    return np.array(coords)

def pdbToArray(pdbfilePath, hetatm = True, keepH = False):
    # Read pdb coordinates into numpy array. Remove hydrogen atoms
    # hetatm: (Boolean) if True, read coordinates for ATOM and HETATM. If False, read only coordinates ATOM
    with open(pdbfilePath, 'r') as f:
        lines = f.readlines()
#### Modified 2021-03-05. Bug reading PDB files. updated so that compatible with PDB files formatted as according to Chimera (https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html)
    lines = [x for x in lines if len(x) > 4]
    if hetatm:
        lines = [x for x in lines if x[:4] in ['ATOM', 'HETA']]
    else:
        lines = [x for x in lines if x[:4] == 'ATOM']
    if not keepH:
        lines = [x for x in lines if x[13] != 'H']
    coords = [[x[30:38], x[38:46], x[46:54]] for x in lines]
####
#    lines = [x.split() for x in lines]
#    if hetatm:
#        lines = [x for x in lines if x[0] in ['ATOM', 'HETATM']]
#    else:
#        lines =[x for x in lines if x[0] == 'ATOM']
#    if keepH:
#        lines = [x for x in lines]
#    else:
#        lines = [x for x in lines if x[2][0] != 'H']
#    coords = [x[6:9] for x in lines]
    for i in range(len(coords)):
        coords[i] = [float(x) for x in coords[i]]
    return np.array(coords)

def molToArray(molfilePath, keepH = False):
    with open(molfilePath, 'r') as f:
        lines = f.readlines()
    lines = [x.strip() for x in lines[4:]]
    lines = [x.split() for x in lines]
    lines = [x for x in lines if len(x)>10]
    if keepH:
        coords = [x[:3] for x in lines]
    else:
        coords = [x[:3] for x in lines if x[3] != 'H']
    for i in range(len(coords)):
        coords[i] = [float(x) for x in coords[i]]
    return np.array(coords)

def fileToArray(filepath, hetatm = True, keepH=False):
    # Read xyz or pdb files into numpy arrays
    file_extension = filepath.split('.')[-1]
    if file_extension == 'xyz':
        return xyzToArray(filepath, keepH)
    elif file_extension == 'pdb':
        return pdbToArray(filepath, hetatm, keepH)
    elif file_extension == 'mol':
        return molToArray(filepath, keepH)
    else:
        print('File extension not recognized when reading. Only files with extensions "xyz", "pdb", or "mol" are read and considered.')
        print('Program terminated. Either remove exception or modify file path and be sure they correspond to standard xyz or pdb format')
        quit()

def rbfKernel(r, epsilon = 1.0):
    return np.exp(-(epsilon*r)**2)

def clusterStats(RMSDmatrix, communityAssignment, centers, fileList, cluster_names):
    # Report on the speads in each community. Use RMSD matrix as a complete, weighted graph. 
    # inputs
    # RMSDmatrix: (numpy.array) matrix of RMSD values between all conformers
    # communityAssignment: (list of int) list specifying the assigned community to each conformer
    # centers: (list of string) file names of centers
    # fileList: (list of string) file names of all conformers
    # refer to the definitions of radius and diameter in graphs (https://mathworld.wolfram.com/GraphEccentricity.html)
    # output: (pandas.DataFrame) DF containing radius, diameter, mean RMSD of each cluster.
    
    # get communities in the order of centroids
    centerIndices = [fileList.index(x) for x in centers]
    communities = [communityAssignment[x] for x in centerIndices]
    
    sizes = []
    diameters = []
    meanRMSDs = []

    for C in communities:
        C_members = [x for x in range(RMSDmatrix.shape[0]) if communityAssignment[x] == C]
        num_members = len(C_members)
        community_matrix = RMSDmatrix[C_members,:][:,C_members]
        diameter = np.max(community_matrix)
        size = len(C_members)
        meanRMSD = np.sum(community_matrix) / (size * (size - 1)) # 2 * (n + 1) * (n / 2) simplified

        sizes.append(size)
        diameters.append(diameter)
        meanRMSDs.append(meanRMSD)

    diameter = np.max(RMSDmatrix)
    size = len(communityAssignment)
    meanRMSD = np.sum(RMSDmatrix) / (size * (size - 1))

    sizes.append(size)
    diameters.append(diameter)
    meanRMSDs.append(meanRMSD)

    outDf = pd.DataFrame({'size': sizes, 'diameter': diameters, 'mean_RMSD': meanRMSDs}, index = cluster_names + ['global'])

    center_matrix = RMSDmatrix[centerIndices,:][:,centerIndices]
    diameter = np.max(center_matrix)
    meanRMSD = np.mean(center_matrix)
    size = len(center_matrix)
    meanRMSD = np.sum(center_matrix) / (size * (size - 1))

    centerRow = pd.DataFrame({'size': [size], 'diameter': [diameter], 'mean_RMSD': [meanRMSD]}, index = ['centers'])
    outDf = outDf.append(centerRow)

    return outDf

def arrayToXyz(coordList, fileList, inpath, outpath, clus, keepH=False):
    # make directory populated with individual xyz structures, after processing.
    # inputs: coordList: list of numpy arrays with xyz coordinates
    # fileList: list of the names of files
    # inpath: path of a reference xyz file. needed only for list of atoms corresponding to each coordinate
    # oupath: destination of xyz file
    # clus: (string) name of the cluster. If name is 'centers', then append 'center_' to the file name.

    # First, get elements in the order corresponding to the coordinates
    with open(inpath, 'r') as f:
        elements = f.readlines()
    
    elements = [x.split() for x in elements]
    if keepH:
        elements = [x[0] for x in elements if len(x) == 4]
    else:
        elements = [x[0] for x in elements if len(x) == 4 and x[0] != 'H']

    for i in range(len(fileList)):
        coord = coordList[i]
        commentLine = fileList[i].split('.')[0][4:]
        if clus == 'centers':
            fileList[i] = 'center_' + fileList[i]
        content = coord.tolist()
#        print(len(content), len(elements))
        for j in range(len(content)):
            content[j] = [str(x) for x in content[j]]
            content[j] = elements[j] + '\t' + '\t'.join(content[j])
        content = '\n'.join(content)
        with open(outpath + fileList[i], 'w') as g:
            g.write(str(len(elements)) + '\n')
            g.write(commentLine + '\n')
            g.write(content)

def arrayToXyzFromMol(coordList, fileList, referenceFilePath, outpath):
    with open(referenceFilePath, 'r') as f:
        lines = f.readlines()
    lines = [x.strip() for x in lines[4:]]
    lines = [x.split() for x in lines]
    lines = [x for x in lines if len(x)>10]
    elements = [x[3] for x in lines if  x[3] != 'H']
    for i in range(len(fileList)):
        coord = coordList[i]
        content = coord.tolist()
        for j in range(len(content)):
            content[j] = [str(x) for x in content[j]]
            content[j] = elements[j] + '\t' + '\t'.join(content[j])
        content = '\n'.join(content)
        with open(outpath + fileList[i][:-3]+'xyz', 'w') as g:
            g.write(content)

def save_outputs(inpath, outpath, fileList, communityAssignments, centralNodes, cluster_names):
    centerIndices = [fileList.index(x) for x in centralNodes]
    community_list = [communityAssignments[x] for x in centerIndices]
    for i in range(len(community_list)):
        clusterpath = outpath + cluster_names[i]
        os.system('mkdir ' + clusterpath)
        C_members = [x for x in range(len(fileList)) if communityAssignments[x] == community_list[i]]
        for mem in C_members:
            os.system('cp ' + inpath + fileList[mem] + ' ' + clusterpath)
        
    clusterpath = outpath + 'centers/'
    os.system('mkdir ' + clusterpath)
    for centroid in centralNodes:
        os.system('cp ' + inpath + centroid + ' ' + clusterpath)

def cluster_summary(outpath, fileList, communityAssignments, centralNodes, cluster_names):
    # Save csv file in which the classification of each conformer is specified.
    centerIndices = [fileList.index(x) for x in centralNodes]
    community_list = [communityAssignments[x] for x in centerIndices]

    outDf = pd.DataFrame()
    for i in range(len(community_list)):
        C_indices = [x for x in range(len(fileList)) if communityAssignments[x] == community_list[i]]
        C_members = [fileList[x] for x in C_indices]
        Cdf = pd.DataFrame(cluster_names[i], index = C_members, columns = ['cluster'])
        outDf = outDf.append(Cdf, sort = False)
        
    outDf['center'] = 0
    outDf.loc[centralNodes, 'center'] = 1

    with open(outpath + 'cluster_summary.csv', 'w') as f:
        outDf.to_csv(f)
