import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import random
import itertools
from tqdm import tqdm

if __name__ == "__main__":
    con_name = "Starlink"
    np.random.seed(40)
    G1 = nx.grid_2d_graph(72, 22, periodic=True)
    G2 = nx.grid_2d_graph(36, 20, periodic=True)
    # nx.draw(G1)
    # plt.show()
    G = nx.union(G1, G2, rename=("G1-", "G2-"))
    rho_list = np.linspace(0.0, 2.0, 51)
    node2_list = [f"G2-{node_i}" for node_i in G2.nodes]
    node1_list = [f"G1-{node_i}" for node_i in G1.nodes]
    random.shuffle(node2_list)
    random.shuffle(node1_list)
    k = 0
    delta_rho = (rho_list[-1] - rho_list[0]) / (len(rho_list) - 1)
    diameter_list = [nx.diameter(G1)]
    aver_path_length_list = [nx.average_shortest_path_length(G1)]
    rho_edge_list = []
    for rho in tqdm(rho_list[1:]):
        for u in node2_list:
            if (rho - delta_rho < np.random.rand() < rho) and (rho <= 1):
                v = node1_list[k]
                k += 1
                G.add_edge(v, u)
            # elif (rho > 1) and (rho - 1 - delta_rho < np.random.rand() < rho - 1):
            #     v = node1_list[k]
            #     k += 1
            #     G.add_edge(v, u)
            if k == len(node1_list) - 1:
                break
        Total_ILL = nx.number_of_edges(G)
        rho_edge_list.append(k / Total_ILL)
        diameter_list.append(nx.diameter(G))
        aver_path_length_list.append(nx.average_shortest_path_length(G))
    np.save(f"diameter_network/{con_name}_aver_length.npy", np.array(aver_path_length_list))
    np.save(f"diameter_network/{con_name}_diameter.npy", np.array(diameter_list))
