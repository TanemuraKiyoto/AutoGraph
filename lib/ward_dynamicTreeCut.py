# 2020-05-22

# Kiyoto Aramis Tanemura

# The Ward algorithm for conformational clustering is used in building Markov State Model (DOI: 10.1021/acs.jctc.6b01238). To circumvent threshold selection, we will apply the dynamic tree cut method, used in conjunction to the Ward dendogram for automated conformational clustering (DOI: 10.1186/s13321-017-0208-0).

import numpy as np
from pandas import read_csv
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import squareform

def get_ward_dendrogram(rmsd_mat):
    '''Use Scipy functions to obtain dendrogram using Ward method
    Returns the dendrogram (refer to Scipy linkage output format) as a ClusterNode object'''
    dend = linkage(squareform(rmsd_mat),method = 'ward', optimal_ordering = True)
    return dend

def goLeftmost(node, path, path_record):
    '''provided a node, travel left until a leaf is reached'''
    curr_node = node
    if path in path_record:
        return path
    while not curr_node.is_leaf():
        curr_node = curr_node.get_left()
        path += ('l')
    return path

def travelDown(ref_node, rel_path):
    '''reach a node below, provided a starting node an path'''
    curr_node = ref_node
    for i in rel_path:
        if i == 'l':
            curr_node = curr_node.get_left()
        elif i == 'r':
            curr_node = curr_node.get_right()
    return curr_node
    
def getHeights(root):
    '''return heights of nonleaf nodes and their corresponding paths'''
    heights = []
    path_recorded = []
    path = '.'
    curr_node = root
    if root.is_leaf():
        return [0], [path]
    while 'l' in list(path) or not curr_node.is_leaf():
        path = goLeftmost(curr_node, path, path_recorded)
        path = path[:-1]
        curr_node = travelDown(root, path)
        if path not in path_recorded:
            heights.append(curr_node.dist)
            path_recorded.append(path)
            path += 'r'
            curr_node = travelDown(root, path)
    return heights, path_recorded

def treeCutCore(H, I, tau = 5):
    '''Determine significant clusters provided one calibration value'''
    H_hat = H - I
    trans_indices = [x for x in range(len(H)-1) if H_hat[x] > 0 and H_hat[x+1] < 0]
    breakpoints = []
    for i in trans_indices:
        back_index = 1
        while H_hat[i-back_index] > 0 and back_index <= i:
            back_index += 1
        breakpoints.append(i - back_index + 1)
    # Find significant breakpoints
    significant_breakpoints = [breakpoints[x] for x in range(len(breakpoints)) if trans_indices[x] - breakpoints[x] > tau]
    return significant_breakpoints

def adaptiveTreecutCore(H, tau = 5):
    '''Perform treeCutCore at mean height. If no significant breakpoints are detected,
    continue the operation below and above the mean.'''
    if len(H) == 0:
        return []
    lm = np.mean(H)
    lu = np.mean([lm, np.max(H)])
    ld = np.mean([lm, np.min(H)])
    bps = treeCutCore(H, lm, tau)
    if len(bps) == 0:
        bps = treeCutCore(H, ld, tau)
        if len(bps) == 0:
            bps = treeCutCore(H, lu, tau)
    return bps

def getClusterNodeIndices(comprehensive_path_list, cluster_substring):
    # Given a root string, return all permutations as indices
    return [x for x in range(len(comprehensive_path_list)) if comprehensive_path_list[x].startswith(cluster_substring)]

# After looking at the Java implementation, I suspect the breakpoints, corresponding to the indices of heights (distance value of nonleaf nodes),
# also correspond to the indices of leaves. Try clustering only by collecting all breakpoints, then using the indices to subset leaves.
def dynamicTreeCut(tree, n, tau = 5):
    allHeights, allPaths = getHeights(tree)
    allHeights = [0] + allHeights + [0] # Sandwitch the heights with zero so that the ends can be included or left out from major clusters.
    breakpoints = [0,-1]
    updateList = [-1]
    while len(updateList) > 0:
        updateList = []
        for i in range(len(breakpoints) - 1):
            Hi = allHeights[breakpoints[i]:breakpoints[i+1]]
            cutpoints = adaptiveTreecutCore(Hi, tau)
            updateList += [x for x in cutpoints if x not in breakpoints]
        breakpoints += updateList
    return breakpoints

def report_assingments(breakpoints, tree):
    '''To standardize outputs with the other clustering algorithms, this function takes the list of ClusterNodes produced in dynamicTreeCut and returns a list of ints specifying the cluster assignment'''
    leaves = tree.pre_order(lambda x: x.id)
    n = len(leaves)
    comm_assing = np.zeros(n, int)
    for i in range(len(breakpoints) - 1):
        members = leaves[breakpoints[i]:breakpoints[i+1]]
        comm_assing[members] = i
    return comm_assing.tolist()

def ward_dynamicTreeCut(rmsd_mat, tau = 5):
    n = rmsd_mat.shape[0]
    dend = get_ward_dendrogram(rmsd_mat)
    tree = to_tree(dend)
    breakpoints = dynamicTreeCut(tree, n, tau)
    comm_assing = report_assingments(breakpoints, tree)
    return comm_assing
