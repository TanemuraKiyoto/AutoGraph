# 2021-03-03

# Kiyoto Aramis Tanemura

# Show how to run benchmark clustering algorithms. Please cite the original algorithms if used in your research.

import sys
sys.path.append('/mnt/home/tanemur1/2020-05-28')
from AutoGraph import DynamicTreeCut
from AutoGraph import NMRCLUST
from AutoGraph import RCKmeans

inpath = 'data/conformers_oslh'
outpath = 'results/benchmarks/dynamicTreeCut_oslh'

# Run the Dynamic Tree Cut algorithm on the OSLH conformers

dtc = DynamicTreeCut()
dtc.run(inpath, outpath)

# Run the NMRCLUST algorithm on the OSLH conformers

outpath = 'results/benchmarks/nmrclust_oslh'

nc = NMRCLUST()
nc.run(inpath, outpath)

# Run the Representative Confomer K-Means algorithm on the OSLH conformers
# The run method raises an error. This is because the window for calculating simple moving average is set to ten,
# which is the current length of the data subset. If the dataset is of sufficient samples, the syntax below can be used.

outpath = 'results/benchmarks/rckmeans_oslh'

rck = RCKmeans()
#rck.run(inpath, outpath)
