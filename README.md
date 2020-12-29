# Welcome to AutoGraph version 1.0

Perform autonomous graph based clustering on metabolite conformers to reduce redundancy from configuration space. AutoGraph prioritizes autonomation, simplicity, and practicality.

## Brief explanation of AutoGraph:
Atomic RMSD is computed between all conformers (or a randomized sample) in a directory. RMSD values are transformed using a Gaussian kernel function to build an affinity matrix between conformers. Edges with low weights are removed by applying the maximum threshold to yield a graph that has exactly one component. Clusters are identified using the Louvain algorithm. Centroids are selected from each cluster as the lowest energy conformer.

## Quick start guide:
1. Install any missing dependencies (Pandas, Numpy, Scipy)
2. Copy this package to a convenient location (e.g. in src directory of your working directory). Specify paths if located elsewhere
3. Consolidate conformer files (XYZ, PDB, or MOL) in one directory
4. If energy values have been computed for each conformers, save it in a csv file. The indices located at the left-most column must be file names (e.g. 'opt-ani_geom_123.xyz'). The column name should be 'energy' or specified when you call the AutoGraph function. If not computed, the cluster centers are chosen using a graph based metric
5. Choose whether you want to run AutoGraph interactively through a program or in your own python script. The program will prompt you for the inputs and perform AutoGraph on your data. Choose the program if you have only a few systems to consider or for a demo. Choose the script if you need to automate the protocol over many systems or you are recording the metrics over many clustering protocols.\
----interactive program route----
6. Run AutoGraph.py located in the package ("python AutoGraph/AutoGraph.py")
7. Answer the series of prompts and the AutoGraph protocol will be deployed on your data. The message 'Job complete' indicates completion of job. Check results saved in the output path\
-----automated script route----
6. Import the function in your python script >>> from AutoGraph import AutoGraph
7. Specify the parameters when constructing the AutoGraph instance by supplying the relevant arguments (see below) e.g. >>> ag = AutoGraph(randomize = True)
8. Run the clustering by executing the run(...) method. e.g. >>> ag.run('data/conformers/', 'results/testAG', Epath = 'data/energy_summary.csv', E_label = 'HF_3-21')
9. The message 'Job complete' indicates completion of job. Check results saved in the path to the output, specified in step 6.

## the AutoGraph class:
AutoGraph(randomize = False, subset = 0, copy_conformers = True, silence = False, hetatm = True)\
--------parameters-------\
randomize: (Boolean) randomize file order. The Louvain clustering algorithm is sensitive to order\
subset: (int) RMSD calculation scales O(n(n-1)/2). If too many conformers present, take a subset of files to perform clustering, then assign remainder by lowest RMSD to centroids\
copy_conformers: (Boolean) make subdirectories for clusters populated by their conformers. Turn to False if copying all structures compromises device memory\
silence: (Boolean) if True, forgo all print statements\
hetatm: (Boolean) if True, read ATOM and HETATM for PDB files. If False, only read ATOM

## the AutoGraph.run(...) method
AutoGraph.run(inpath, outpath, Epath = '', E_label = 'energy')\
----------inputs---------\
inpath: (string) path to the input xyz files (e.g. 'data/xyzfiles/')\
outpath: (string) path to the output. Directory does not need to already exist\
-----optional inputs-----\
Epath: (string) path to energy data csv file for choosing centroids. If not provided, choose centroid by max in-community weighted degree\
E_label: (string) column name of energy values for the file located at Epath.

Once you run the method, it will save relevant files to the path specified by outpath

## Description of output files and directories:
### files
- affinityMatrix.csv: The symmetric affinity matrix made by applying a Gaussian kernel to the RMSD matrix. Indices and columns correspond to the file names.
- filteredAffinityMatrix.csv: The affinity matrix after applying the adaptive threshold to remove low weight edges.
- communityStats.csv: Descriptive statistics applied to each cluster as well as the whole data ('global') and collection of centroids ('centers')
- cluster_summary.csv: Clustering output by each file. The cluster column specifies the cluster to which the file is assigned. The centroid column specifies whether the file is a centroid (1) or not (0)
- rmsdMatrix.csv: The symmetric RMSD matrix made by calculating the atomic RMSD between all conformers considered. Indices and columns correspond	 to the file names.
### directories
- clusterX: contains all files assigned to clusterX
- centers: Contains all files chosen as centroids. 

## Possible questions/problems:
- Do I need to know Python to use AutoGraph?
  No, you can run the program interactively, with the program walking you through the necessary inputs. Refer to the Quick Start Guide 
- Why can't I cluster the conformers saved as PDB files?
  The variability in PDB file format limits the deployment of AutoGraph on all PDB files. Please consider drafting a script to convert your PDB file to XYZ files to perform AutoGraph. Run AutoGraph with copy_conformers = False. Then use the cluster_summary.csv saved in the output directory to inform your downstream use of clustered conformers
- Why do I get different clustering results depending on the iteration?
  The Louvain algorithm used for clustering is deterministic, yet the result depends on file order (https://arxiv.org/pdf/0803.0476.pdf). A change in file order may result in a change in clustering result.
- I want to take the average over the AutoGraph output iterated over several randomization
  Call the AutoGraph function inside a loop. Multiple calls, however, may duplicate files in the outputs or overwrite summary files. Therefore, specify copy_conformers=False in the function, and rename the output files (e.g. "cluster_summary.csv" -> "iter3_cluster_summary.csv") between calls. Also, consider specifying silence=True in the argument to avoid crowding the screen
- The clustering terminates with an error when finding the threshold weight.
  If the values for weights in the affinity matrix are too small, even the complete graph is registered as empty.
- Can I run AutoGraph on protein conformations?
  Yes. You may want to extract the backbone atoms only and save files in a dedicated directory as PDB file (e.g. >>> grep 'CA' original1.pdb > backbone/extracted1.pdb ). Be sure all atoms are consistent between all files. Then run AutoGraph on the files in the directory containing only backbone atoms.
- Can I run AutoGraph on protein-protein complexes?
  With some work, yes. AutoGraph computes an atomic RMSD matrix based on superimposing all atoms, which is not typically used for PPI complexes. If the RMSD matrix is provided as a csv in the outpath, AutoGraph will read this matrix and perform clustering based on it. Run AutoGraph on a small subset of your data to see the formatting. Put your desired similarity value (ligand RMSD, interface RMSD, etc) in this format, save it as 'rmsdMatrix.csv' in the outpath, then run AutoGraph. AutoGraph will not include PPI based metrics in the future due to the variablibity in how to read the files, which compromises the integrity of its results.

## Notes:
For comparison purpose, we included our implementation of the following conformational clustering algorithms. If these algorithms were used in your publication, please reference the original publications of the algorithms:
- NMRCLUST: Kelley, L.; Gardner, S.; Sutcliffe, M. Prot. Eng. 1996, 9, 1063. DOI: 10.1093/protein/9.11.1063
- Representative Conformer Kmeans: Kim, H.; Jang, C.; Yadav, D. K.; Kim, M. H. J. Cheminform. 2017. 9, 1. DOI: 10.1186/s13321-017-0208-0
- Ward's method: Ward, J. H., Jr. J. Am. Stat. Assoc. 1963, 58, 236.
- Dynamic Tree Cut algorithm: Langfelder, P.; Zhang, B.; Horvath, S. Bioinformatics 2008, 24,719.
- Ward + DynamicTreeCut applied to find representative conformers: Kim, H.; Jang, C.; Yadav, D. K.; Kim, M. H. J. Cheminform. 2017. 9, 1. DOI: 10.1186/s13321-017-0208-0

If you use AutoGraph in your publication, please cite as follows:
Tanemura, Kiyoto; Das, Susanta; M. Merz Jr., Kenneth (2020): AutoGraph: Autonomous Graph Based Clustering of Small-Molecule Conformations. ChemRxiv. Preprint. https://doi.org/10.26434/chemrxiv.13491543.v1 

Kiyoto Aramis Tanemura
Kenneth Merz Research Group
Department of Chemistry
Michigan State University
tanemur1@msu.edu

2020-05-17
Last updated 2020-12-15