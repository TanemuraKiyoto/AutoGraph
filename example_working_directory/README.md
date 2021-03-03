2021-03-02

Kiyoto Aramis Tanemura

The example working directory provides example inputs/outputs for the AutoGraph protocol. The data directory contains a very small subset of conformers of O-succinyl-L-homoserine. Navigate to this directory and run AutoGraph either 1) interactively, or 2) by scripts.

1) Interactively

>>> python ../AutoGraph.py
Welcome to the AutoGraph conformational clustering program.
Provide the inputs to run AutoGraph on your data.
To begin, enter "continue" to the prompt.
>>> continue
Do you wish to use default parameters and proceed to AutoGraph? (y/n) If so, you will directly proceed to entering AutoGraph inputs. If not, you will first be prompted to enter parameters.
>>> y
Enter input path to the directory containing XYZ or PDB files to cluster. (Where your files are located on your computer.)
>>> data/conformers_oslh
Enter output path to the directory to save results.
>>> results/oslh
To choose representative conformers per cluster, AutoGraph chooses the lowest energy structure. Do you have these energy values stored in the format specified in README.md? (y/n)
>>> y
Enter path to the energy csv file
>>> data/energy.csv
Enter label of the column containing energy values.
>>> ANI
AutoGraph protocol is now running on your data.
threshold: 0.108667
threshold RMSD: 1.489788
Cluster statistics: 
          size  diameter  mean_RMSD
cluster0     3  1.559443   1.164474
cluster1     3  1.785357   1.410702
cluster2     1  0.000000        NaN
cluster3     3  0.789062   0.596409
global      10  2.239881   1.549548
centers      4  2.239881   1.774936
Number of structures:  10
Number of clusters:  4
Job complete

2) By scripts

>>> python src/cluster_oslh.py

The output should be similar to the code provided in the interactive mode.