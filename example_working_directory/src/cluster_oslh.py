# 2021-03-02

# Kiyoto Aramis Tanemura

# Short example script to cluster O-succinyl-L-homoserine conformers by the AutoGraph protocol.

# import the AutoGraph package
import sys
sys.path.append('/mnt/home/tanemur1/2020-05-28')
from AutoGraph import AutoGraph

# specify paths
inpath = 'data/conformers_oslh'
outpath = 'results/oslh'
energyFile = 'data/energy.csv'

# assing an AutoGraph instance with relevant arguments
ag = AutoGraph(copy_conformers = False)
# run the protocol. Provide the in/out paths, as well as the conformers' energies
ag.run(inpath, outpath, Epath = energyFile, E_label = 'ANI')

# This will save the files to the output path
