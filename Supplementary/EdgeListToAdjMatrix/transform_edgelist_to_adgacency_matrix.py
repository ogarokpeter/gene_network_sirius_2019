import numpy as np


filename = "/home/user/Sirius/output1w.sif"


def edgelist_to_adjacency_matrix(filename):
    g1 = []
    g2 = []
    prob = []
    with open(filename, 'r') as f:
        for line in f:
            try:
                ls = line.split()
                g1.append(ls[0])
                g2.append(ls[2])
                prob.append(float(ls[1]))
            except:
                pass
    n_genes = 0
    genes = dict()
    for i in range(len(g1)):
        a = 0
        if g1[i] not in genes.keys():
            genes[g1[i]] = n_genes
            a = n_genes
            n_genes += 1
        else:
            a = genes[g1[i]]
        b = 0
        if g2[i] not in genes.keys():
            genes[g2[i]] = n_genes
            b = n_genes
            n_genes += 1
        else:
            b = genes[g2[i]]
        g1[i] = a
        g2[i] = b

    adjacency_matrix = np.zeros(shape=(n_genes, n_genes))
    for i in range(len(g1)):
        adjacency_matrix[g1[i]][g2[i]] = prob[i]
        adjacency_matrix[g2[i]][g1[i]] = prob[i]
    return adjacency_matrix

print(edgelist_to_adjacency_matrix(filename))
