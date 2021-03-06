# 2020-05-25

# Kiyoto Aramis Tanemura

# Perfrom clustering by DynamicTreeCut on Ward dendrogram.
# If used in a publication, please reference the original publication of the algorithms:
# Ward's method: Ward, J. H., Jr. J. Am. Stat. Assoc. 1963, 58, 236.
# Dynamic Tree Cut algorithm: Langfelder, P.; Zhang, B.; Horvath, S. Bioinformatics 2008, 24,719.
# Ward + DynamicTreeCut applied to find representative conformers: Kim, H.; Jang, C.; Yadav, D. K.; Kim, M. H. J. Cheminform. 2017. 9, 1. DOI: 10.1186/s13321-017-0208-0

# Import the DynamicTreeCut class and supply the relevant parameters at declaration.
# Deploy the algorithm using the run(...) method.

# e.g.
# >>> from AutoGraph import DynamicTreeCut
# >>> dtc = DynamicTreeCut()
# >>> dtc.run('data/conformers/', 'results/test_DTC_output/')

import numpy as np
from pandas import DataFrame, read_csv
import os
if __name__ == '__main__':
    from lib.functions import *
    from lib.rmsd import rmsdMatrix
    from lib.ward_dynamicTreeCut import ward_dynamicTreeCut
    from lib.select_centroids import centroid_medoid
else:
    from .lib.functions import *
    from .lib.rmsd import rmsdMatrix
    from .lib.ward_dynamicTreeCut import ward_dynamicTreeCut
    from .lib.select_centroids import centroid_medoid

class DynamicTreeCut:
    '''The comparison of automated clustering algorithms for resampling representative conformer ensembles with RMSD matrix. DOI: 10.1186/s13321-017-0208-0'''

    def __init__(self, tau = 5, copy_conformers = True, silence = False, hetatm = True):
        '''tau: (int) threshold forward run length to consider a breakpoint significant. Refer to original publication for clarification.\ncopy_conformers: (Boolean) make subdirectories for clusters populated by their conformers. Turn to False if copying all structures compromises device memory\nsilence: (Boolean) if True, forgo all print statements. Consider if calling AutoGraph in a loop to avoid crowding the screen. \nhetatm: (Boolean) if the input files are PDB, specify whether or not to read the coordinates of HETATM'''
        self.tau = tau
        self.copy_conformers = copy_conformers
        self.silence = silence
        self.hetatm = hetatm

    def formatPaths(self, inpath, outpath):
        '''Append "/" to end of path if absent'''
        if inpath[-1] != '/':
            inpath = inpath + '/'
        if outpath[-1] != '/':
            outpath = outpath + '/'
        return inpath, outpath

    def set_up(self, inpath, outpath):
        '''Set up the output directory and read input data'''
        inpath, outpath = self.formatPaths(inpath, outpath)
        self.fileList = [x for x in os.listdir(inpath) if x[-3:] in ['xyz', 'pdb', 'mol']]
        self.numFiles = len(self.fileList)
        if not os.path.exists(outpath):
            os.system('mkdir ' + outpath)
        self.rmsdMatrix = None
        if 'rmsdMatrix.csv' in os.listdir(outpath):
            self.rmsdMatrix = read_csv(outpath + 'rmsdMatrix.csv', index_col = 0)
            self.rmsdMatrix = self.rmsdMatrix.reindex(index = self.fileList, columns = self.fileList)
            self.rmsdMatrix = self.rmsdMatrix.to_numpy()

    def getRmsdMatrix(self, inpath, outpath):
        '''If RMSD matrix has not been read, compute and save it'''
        inpath, outpath = self.formatPaths(inpath, outpath)
        if type(self.rmsdMatrix) == type(None):
            self.rmsdMatrix = rmsdMatrix(self.fileList, inpath, self.hetatm)
            rmsdDf = DataFrame(self.rmsdMatrix, index = self.fileList, columns = self.fileList)
            with open(outpath + 'rmsdMatrix.csv', 'w') as g:
                rmsdDf.to_csv(g)
            del rmsdDf

    def report_save(self, inpath, outpath, communityAssignment):
        '''Report clustering results and save files'''
        inpath, outpath = self.formatPaths(inpath, outpath)
        cluster_names = ['cluster' + str(x) for x in range(len(self.centroids))]
        stats = clusterStats(self.rmsdMatrix, communityAssignment, self.centroids, self.fileList, cluster_names)
        with open(outpath + 'communityStats.csv', 'w') as h:
            stats.to_csv(h)
        if not self.silence:
            print('Cluster statistics: ')
            print(stats)

        cluster_summary(outpath, self.fileList, communityAssignment, self.centroids, cluster_names)

        if self.copy_conformers:
            save_outputs(inpath, outpath, self.fileList, communityAssignment, self.centroids, cluster_names)

        if not self.silence:
            print('Number of structures: ', len(self.fileList))
            print('Number of clusters: ', len(self.centroids))
            print('Job complete')

    def run(self, inpath, outpath):
        '''Perform RCKmeans conformational clustering'''
        inpath, outpath = self.formatPaths(inpath, outpath)
        self.set_up(inpath, outpath)
        self.getRmsdMatrix(inpath, outpath)
        communityAssignment = ward_dynamicTreeCut(self.rmsdMatrix)
        self.centroids = centroid_medoid(self.fileList, communityAssignment, self.rmsdMatrix)
        self.report_save(inpath, outpath, communityAssignment)

# Interactive program to perform Representative Conformer K-means on files of conformations.
# If used in a publication, please reference the original article of the algorithm: Kim, H.; Jang, C.; Yadav, D. K.; Kim, M. H. J. Cheminform. 2017. 9, 1. DOI: 10.1186/s13321-017-0208-0

if __name__ == '__main__':

    print('Welcome to our implementation of the Ward/DynamicTreeCut conformational clustering program.\nProvide the inputs to run DynamicTreeCut on your data.\nTo begin, enter "continue" to the prompt.')
    continueToProgram = input()
    if continueToProgram != 'continue':
        print('End program.')
        quit()
    print('What value to use as the threshold forward run length tau? (Default 5)')
    tau = int(input())
    print('Do you wish to have files copied in subdirectories corresponding to each cluster? (y/n)')
    copy_them = input()
    if copy_them.lower().strip() in ['y', 'yes', 'yeah', 'yup', 'sure', 'okay']:
        copy_conformers = True
    elif copy_them.lower().strip() in ['n', 'no', 'nope', 'nah', 'no thanks']:
        copy_conformers = False
    else:
        print('Response not recognized. Conformer files will be copied into respective subdirectories.')
        copy_conformers = True
    dtc = DynamicTreeCut(tau, copy_conformers)
    print('Enter input path to the directory containing XYZ/PDB/MOL files to cluster. (Where your files are located on your computer.)')
    inpath = input()
    print('Enter output path to the directory to save results.')
    outpath = input()
    print('DynamicTreeCut protocol is now running on your data.')
    dtc.run(inpath, outpath)
