# 2020-12-15

# Kiyoto Aramis Tanemura

# Cluster metabolite conformers by the AutoGraph protocol.

import numpy as np
from pandas import DataFrame, read_csv
import os, sys
from random import sample

if True:
if __name__ == '__main__':
    from lib.functions import *
    from lib.rmsd import rmsd, assign_remaining_files, rmsdMatrix
    from lib.LouvainClustering import Louvain
    from lib.findThreshold import findThreshold
    from lib.select_centroids import centroid_energy, centroid_weighted_degree, centroid_eccentricity, centroid_betweenness
else:
    from .lib.functions import *
    from .lib.rmsd import rmsd, assign_remaining_files, rmsdMatrix
    from .lib.LouvainClustering import Louvain
    from .lib.findThreshold import findThreshold
    from .lib.select_centroids import centroid_energy, centroid_weighted_degree, centroid_eccentricity, centroid_betweenness

class AutoGraph:
    '''AutoGraph: Autonomous graph based conformational clustering.\n'''
    def __init__(self, randomize = False, subset = 0, copy_conformers = True, silence = False, hetatm = True, centroid_selection='betweenness'):
        '''randomize: (Boolean) shuffle file order\nsubset: (int) perform protocol only on a randomly sampled subset\ncopy_conformers: {Boolean) save clustered conformer files together in a dedicated directory\nsilence: (Boolean) ignore print statements\nhetatom: (Boolean) if PDB is read, consider coordinates for HETATM as well as ATOM'''
        self.randomize = randomize
        self.subset = subset
        self.copy_conformers = copy_conformers
        self.silence = silence
        self.hetatm = hetatm
        self.threshold = 0.0
        self.centroid_selection = centroid_selection

    def formatPaths(self, inpath, outpath):
        '''Append "/" to path endings if not present.'''
        if inpath[-1] != '/':
            inpath = inpath + '/'
        if outpath[-1] != '/':
            outpath = outpath + '/'
        return inpath, outpath

    def set_up(self, inpath, outpath):
        '''Set up the output directory and read input data'''
        inpath, outpath = self.formatPaths(inpath, outpath)
        self.master_fileList = [x for x in os.listdir(inpath) if x[-3:] in ['xyz', 'pdb', 'mol']]
        self.sample_fileList = self.master_fileList
        self.numFiles = len(self.sample_fileList)
        if not os.path.exists(outpath):
            os.system('mkdir ' + outpath)
        if self.subset > 0:
            self.sample_fileList = sample(self.master_fileList, self.subset)
            self.sample_fileList = sample(self.master_fileList, self.subset)
            self.numFiles = self.subset
        elif self.randomize:
            self.sample_fileList = sample(self.master_fileList, self.numFiles)
        self.rmsdMatrix = None
        self.affinityMatrix = None
        if 'rmsdMatrix.csv' in os.listdir(outpath) and self.subset == 0:
            self.rmsdMatrix = read_csv(outpath + 'rmsdMatrix.csv', index_col = 0)
            self.rmsdMatrix = self.rmsdMatrix.reindex(index = self.sample_fileList, columns = self.sample_fileList)
            self.rmsdMatrix = self.rmsdMatrix.to_numpy()

        if 'affinityMatrix.csv' in os.listdir(outpath) and self.subset == 0:
            self.affinityMatrix = read_csv(outpath + 'affinityMatrix.csv', index_col = 0)
            self.affinityMatrix = self.affinityMatrix.reindex(index = self.sample_fileList, columns = self.sample_fileList)
            self.affinityMatrix = self.affinityMatrix.to_numpy()

    def getRmsdMatrix(self, inpath, outpath):
        '''If RMSD matrix has not been read, compute and save it'''
        inpath, outpath = self.formatPaths(inpath, outpath)
        if type(self.rmsdMatrix) == type(None):
            self.rmsdMatrix = rmsdMatrix(self.sample_fileList, inpath, self.hetatm)
            rmsdDf = DataFrame(self.rmsdMatrix, index = self.sample_fileList, columns = self.sample_fileList)
            with open(outpath + 'rmsdMatrix.csv', 'w') as g:
                rmsdDf.to_csv(g)
            del rmsdDf

    def getAffinityMatrix(self, inpath, outpath):
        '''If affinity matrix has not been read, compute and save it'''
        inpath, outpath = self.formatPaths(inpath, outpath)
        if type(self.affinityMatrix) == type(None):
            self.affinityMatrix = rbfKernel(self.rmsdMatrix)
            self.affinityMatrix[range(self.numFiles), range(self.numFiles)] = 0
            affinityDf = DataFrame(self.affinityMatrix, index = self.sample_fileList, columns = self.sample_fileList)
            with open(outpath + 'affinityMatrix.csv', 'w') as h:
                affinityDf.to_csv(h)
            del affinityDf

    def getThreshold(self):
        '''Compute the maximum threshold to have exactly one component.'''
        if self.threshold == 0:# If threshold has not yet been computed
            self.threshold = findThreshold(self.affinityMatrix)
            if not self.silence:
                print('threshold: %f' %self.threshold)
                print('threshold RMSD: %f' %( np.sqrt(-np.log(self.threshold)) ))

    def getFilteredAffinityMatrix(self, filter_threshold, outpath):
        '''Filter edges with weights below threshold'''
        adjacencyMatrix = self.affinityMatrix > filter_threshold
        filteredAffinityMatrix = self.affinityMatrix * adjacencyMatrix
        filteredAffinityDf = DataFrame(filteredAffinityMatrix, index = self.sample_fileList, columns = self.sample_fileList)
        with open(outpath + 'filteredAffinityMatrix.csv', 'w') as j:
            filteredAffinityDf.to_csv(j)
        del filteredAffinityDf
        return filteredAffinityMatrix

    def getCentroids(self, communityAssignment, Epath='', E_label='energy', filteredAffinityMatrix=None):
        '''Return file names of conformers designated as centroids. If energy is provided, find lowest energy conformers in each cluster. Otherwise choose by maximum in-cluster weighted degree'''
        if Epath == '':
            if self.centroid_selection == 'degree':
                self.centroids = centroid_weighted_degree(self.sample_fileList, communityAssignment, filteredAffinityMatrix)
            elif self.centroid_selection == 'eccentricity':
                self.centroids = centroid_eccentricity(self.sample_fileList, communityAssignment, self.rmsdMatrix * self.rmsdMatrix < np.sqrt(-np.log(self.threshold))) # add filtered rmsd matrix
            elif self.centroid_selection == 'betweenness':
                self.centroids = centroid_betweenness(self.sample_fileList, communityAssignment, self.rmsdMatrix * self.rmsdMatrix < np.sqrt(-np.log(self.threshold))) # add filtered rmsd matrix
            else:
                print('centroid criterion not recognized. Use keywords "degree", "eccentricity", or "betweenness" for centroid_selection or provide an energy output to base the selection')
        else:
            self.centroids = centroid_energy(self.sample_fileList, communityAssignment, Epath, E_label)

    def report_save(self, inpath, outpath, communityAssignment):
        '''Report clustering results in save files.'''
        inpath, outpath = self.formatPaths(inpath, outpath)
        cluster_names = ['cluster' + str(x) for x in range(len(self.centroids))]
        if self.subset > 0:
            sample_communityAssignment = list(communityAssignment)
            communityAssignment = assign_remaining_files(self.master_fileList, self.sample_fileList, self.centroids, sample_communityAssignment, inpath)
        else:
            # Note: cluster statistics requires the RMSD matrix. If subset is used, then RMSD matrix is for the subset, not all files.
            # Cluster stats require the full RMSD matrix, so will be computed only for clustering without subsetting.
            stats = clusterStats(self.rmsdMatrix, communityAssignment, self.centroids, self.sample_fileList, cluster_names)
            with open(outpath + 'communityStats.csv', 'w') as h:
                stats.to_csv(h)
            if not self.silence:
                print('Cluster statistics: ')
                print(stats)

        cluster_summary(outpath, self.master_fileList, communityAssignment, self.centroids, cluster_names)

        if self.copy_conformers:
            save_outputs(inpath, outpath, self.master_fileList, communityAssignment, self.centroids, cluster_names)

        if not self.silence:
            print('Number of structures: ', len(self.master_fileList))
            print('Number of clusters: ', len(self.centroids))
            print('Job complete')

    def run(self, inpath, outpath, Epath = '', E_label = 'energy'):
        '''Perform AutoGraph conformational clustering'''
        inpath, outpath = self.formatPaths(inpath, outpath)
        self.set_up(inpath, outpath)
        self.getRmsdMatrix(inpath, outpath)
        self.getAffinityMatrix(inpath, outpath)
        self.getThreshold()
        filtAff = self.getFilteredAffinityMatrix(self.threshold, outpath)
        communityAssignment = Louvain(filtAff, Q_threshold = 0.0, max_iter = 50, resolution = 1.0)
        self.getCentroids(communityAssignment, Epath, E_label, filtAff)
        self.report_save(inpath, outpath, communityAssignment)

# If not imported, but rather executed as the main source code, begin interactive program.

if __name__ == '__main__':

    print('Welcome to the AutoGraph conformational clustering program.\nProvide the inputs to run AutoGraph on your data.\nTo begin, enter "continue" to the prompt.')
    continueToProgram = input()
    if continueToProgram != 'continue':
        print('End program.')
        quit()
    print('Do you wish to use default parameters and proceed to AutoGraph? (y/n) If so, you will directly proceed to entering AutoGraph inputs. If not, you will first be prompted to enter parameters.')
    defaults = input()
    if defaults.lower().strip() in ['y', 'yes', 'yeah', 'yup', 'sure', 'okay']:
        randomize = False
        subset = 0
        copy_conformers = True
    else:
        print('Enter AutoGraph parameters')
        print('1. Do you wish to randomize file order? (y/n)')
        rand_order = input()
        if rand_order.lower().strip() in ['y', 'yes', 'yeah', 'yup', 'sure', 'okay']:
            randomize = True
        elif rand_order.lower().strip() in ['n', 'no', 'nope', 'nah', 'no thanks']:
            randomize = False
        else:
            print('Response not recognized. File order not randomized.')
            randomize = False
        print('2. Do you wish to run AutoGraph on a randomized subset of your data? (y/n)')
        only_subset = input()
        if only_subset.lower().strip() in ['y', 'yes', 'yeah', 'yup', 'sure', 'okay']:
            print('Specify number of files to subset (e.g. enter "100" for hundred files)')
            subset = int(input())
        elif only_subset.lower().strip() not in ['n', 'no', 'nope', 'nah', 'no thanks']:
            print('Response not recognized. AutoGraph will be performed on all files.')
            subset = 0
        else:
            subset = 0
        print('3. Do you wish to have files copied in subdirectories corresponding to each cluster? (y/n)')
        copy_them = input()
        if copy_them.lower().strip() in ['y', 'yes', 'yeah', 'yup', 'sure', 'okay']:
            copy_conformers = True
        elif copy_them.lower().strip() in ['n', 'no', 'nope', 'nah', 'no thanks']:
            copy_conformers = False
        else:
            print('Response not recognized. Conformer files will be copied into respective subdirectories.')
            copy_conformers = True
    ag = AutoGraph(randomize, subset, copy_conformers)
    print('Enter input path to the XYZ or PDB files to cluster. (Where your files are located on your computer.)')
    inpath = input()
    print('Enter output path to the directory to save results.')
    outpath = input()
    print('To choose representative conformers per cluster, AutoGraph chooses the lowest energy structure. Do you have these energy values stored in the format specified in README.md? (y/n)')
    provide_E = input()
    if provide_E.lower().strip() in ['y', 'yes', 'yeah', 'yup', 'sure', 'okay']:
        print('Enter path to the energy csv file')
        Epath = input()
        print('Enter label of the column containing energy values.')
        E_label = input()
        if E_label == '':
            E_label = 'energy'
    elif provide_E.lower().strip() not in ['n', 'no', 'nope', 'nah', 'no thanks']:
        print('Response not recognized. Centroids will be selected by in-cluster weighted degree, not energy.')
        Epath = ''
        E_label = 'energy'
    else:
        print('Centroids will be selected by in-cluster weighted degree, not energy.')
        Epath = ''
        E_label = 'energy'
    print('AutoGraph protocol is now running on your data.')
    ag.run(inpath, outpath, Epath, E_label)
