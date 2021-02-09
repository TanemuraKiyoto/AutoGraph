# 2020-05-20

# Kiyoto Aramis Tanemura

# Our implementation of NMRCLUST conformational clustering algorithm.                        
# If used in a publication, please reference the original article of the algorithm: Kelley, L.; Gardner, S.; Sutcliffe, M. Prot. Eng. 1996, 9, 1063. DOI: 10.1093/protein/9.11.1063
# Import the NMRCLUST class and supply the relevant parameters at declaration.
# Deploy the algorithm using the run(...) method.

# e.g.                                                         
# >>> from AutoGraph import NMRCLUST                                                 
# >>> nmr = NMRCLUST()                    
# >>> nmr.run('data/conformers/', 'results/test_NC_output/')                                  

import numpy as np
from pandas import DataFrame, read_csv
import os
if __name__ == '__main__':
    from lib.functions import *
    from lib.rmsd import rmsd, rmsdMatrix
    from lib.NMRCLUST import NMRCLUST_
    from lib.select_centroids import centroid_medoid
else:
    from .lib.functions import *
    from .lib.rmsd import rmsd, rmsdMatrix
    from .lib.NMRCLUST import NMRCLUST_
    from .lib.select_centroids import centroid_medoid

class NMRCLUST:
    '''An automated approach for clustering an ensemble of NMR- derived protein structures into conformationally related subfamilies - DOI: 10.1093/protein/9.11.1063'''

    def __init__(self, copy_conformers = True, silence = False, hetatm = True):
        '''copy_conformers: (Boolean) make subdirectories for clusters populated by their conformers. Turn to False if copying all structures compromises device memory\nsilence: (Boolean) if True, forgo all print statements. Consider if calling AutoGraph in a loop to avoid crowding the screen\nhetatm: (Boolean) if the input files are PDB, specify whether or not to read the coordinates of HETATM'''
        self.copy_conformers = copy_conformers
        self.silence = silence
        self.hetatm = hetatm

    def formatPaths(self, inpath, outpath):
        '''Append "/" to end of path if not present'''
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
        '''Report clustering results and save outputs'''
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
            print('Number of structures: ', self.numFiles)
            print('Number of clusters: ', len(self.centroids))
            print('Job complete')

    def run(self, inpath, outpath):
        '''Perform NMRCLUST conformational clustering'''
        inpath, outpath = self.formatPaths(inpath, outpath)
        self.set_up(inpath, outpath)
        self.getRmsdMatrix(inpath, outpath)
        communityAssignment = NMRCLUST_(self.rmsdMatrix)
        self.centroids = centroid_medoid(self.fileList, communityAssignment, self.rmsdMatrix)
        self.report_save(inpath, outpath, communityAssignment)

# Interactive program to perform NMRCLUST on files of conformations without drafting a script.
# If used in a publication, please reference the original publication of the algorthm: Kelley, L.; Gardner, S.; Sutcliffe, M. Prot. Eng. 1996, 9, 1063. DOI: 10.1093/protein/9.11.1063

if __name__ == '__main__':
    print('Welcome to our implementation of the NMRCLUST conformational clustering program.\nProvide the inputs to run NMRCLUST on your data.\nTo begin, enter "continue" to the prompt.')
    continueToProgram = input()
    if continueToProgram != 'continue':
        print('End program.')
        quit()
    print('Do you wish to have files copied in subdirectories corresponding to each cluster? (y/n)')
    copy_them = input()
    if copy_them.lower().strip() in ['y', 'yes', 'yeah', 'yup', 'sure', 'okay']:
        copy_conformers = True
    elif copy_them.lower().strip() in ['n', 'no', 'nope', 'nah', 'no thanks']:
        copy_conformers = False
    else:
        print('Response not recognized. Conformer files will be copied into respective subdirectories.')
        copy_conformers = True
    nmr = NMRCLUST(copy_conformers)
    print('Enter input path to the XYZ or PDB files to cluster. (Where your files are located on your computer.)')
    inpath = input()
    print('Enter output path to the directory to save results.')
    outpath = input()
    print('NMRCLUST protocol is now running on your data.')
    nmr.run(inpath, outpath)
