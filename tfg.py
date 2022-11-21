import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

data = pd.read_csv("TFG-Ohmnet_tissue-function-gene.tsv")
data = list(np.array(data))

for i in range(len(data)):
    data[i] = data[i][0].split("\t")[1:]

V1 = [data[i][0] for i in range(len(data))]
V2 = [data[i][1] for i in range(len(data))]
V3 = []

for i in V1:
    if i not in V3:
        V3.append(i)

for i in V2:
    if i not in V3:
        V3.append(i)

G = nx.Graph()
G.add_nodes_from(V3)
G.add_edges_from(data)
L = nx.normalized_laplacian_matrix(G)
e = np.linalg.eigvals(L.toarray())
nx.draw_networkx(G, node_size=10, with_labels=False)
plt.savefig("TFG_Graph.png")
plt.show()

plt.hist(e, bins = 100)
plt.title("Distribution of Eigenvalues (Tissue-specific protein-function associations)")
plt.xlabel("Eigenvalue")
plt.savefig("TFG.png")
plt.show()