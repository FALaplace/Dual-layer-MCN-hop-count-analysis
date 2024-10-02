import numpy as np
import pickle
from Network_Class import Network
import re
import networkx as nx

R_e = 6371
w_e = 2 * np.pi / 86400
d2r = lambda x: x * np.pi / 180
r2d = lambda x: x * 180 / np.pi
lla2xyz = lambda a, b, r: np.array([[r * np.cos(a) * np.cos(b)], [r * np.cos(a) * np.sin(b)], [r * np.sin(a)]])

if __name__ == "__main__":
    con_name = "Starlink"
    con_file = f"constellation/{con_name}/LLA&LLR.npy"
    sat_dict_file = f"constellation/{con_name}/satdict.pkl"
    con_info_file = f"constellation/{con_name}/con_info.txt"
    alldata = np.load(con_file)
    f = open(sat_dict_file, 'rb')
    sat_dict = pickle.load(f)
    sat_list = list(sat_dict.keys())
    f.close()
    N1 = Network(sat_list, alldata, con_info_file)
    rho_list = np.linspace(0.01, 1, 100)
    diam_G = np.zeros_like(rho_list)
    for r in range(rho_list.size):
        rho = np.round(rho_list[r], 2)
        for t in np.random.randint(0, 240, size=4):
            N1.GraphShift(t, rho=rho)
            if not nx.is_connected(N1.G):
                d1 = 100
            else:
                d1 = nx.average_shortest_path_length(N1.G)
            diam_G[r] += d1 / 4
        pass
    filename = f"diameter_network/{con_name}.npy"
    np.save(filename, diam_G)
