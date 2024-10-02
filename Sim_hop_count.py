"""
Create Simulation results at directory Hop_dis_twolayers
Change user1 latitude in [50, 40, 30, 20, 10, 0]
Change rho in [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
Global Hop distribution with size [39, 119]
"""

import numpy as np
import pickle
from Network_Class import Network
import os
import re

CON_NAME_CONFIG = ["Kuiper", "OneWeb", "Starlink", "TeleSat", "Synthetic3", "Starlink_2"]
con_name = "Synthetic3"
con_file = f"constellation/{con_name}/LLA&LLR.npy"
sat_dict_file = f"constellation/{con_name}/satdict.pkl"
con_info_file = f"constellation/{con_name}/con_info.txt"
alldata = np.load(con_file)
f = open(sat_dict_file, 'rb')
sat_dict = pickle.load(f)
sat_list = list(sat_dict.keys())
f.close()
f = open(con_info_file, 'r')
[n1, m1] = f.readline().replace('\n', '').split(',')
[n2, m2] = f.readline().replace('\n', '').split(',')
f.close()
N1 = Network(sat_list, alldata, con_info_file)
dir_name = "Hop_dis_twolayers"
sub_dir_name = f"({n1},{m1})&({n2},{m2})_global_hop_distribution"
check_dir = os.path.join(dir_name, sub_dir_name)
if not os.path.isdir(check_dir):
    os.mkdir(check_dir)
    print(f"{dir_name}/{sub_dir_name} has been created!")

if __name__ == "__main__":
    # N1.GraphInitialization(rho=0.0)
    lat1_list = [50]
    rho_list = [0.0, 0.02, 0.04, 0.06, 0.08] + [i / 10 for i in range(1, 11)]
    for lat1 in lat1_list:
        for rho in rho_list:
            h_dis_location = f"{dir_name}/{sub_dir_name}/"
            h_dis_name = f'lat1={lat1}_rho={rho}.npy'
            print(f"creating {h_dis_name} under {h_dis_location[:-1]}")
            h_aver_dis = np.zeros([39, 119])
            for tk in range(0, 60, 4):
                N1.GraphShift(tk, rho=rho)
                h_dis = N1.hop_count_cal_single_pos(tk, lat1)
                h_aver_dis = h_aver_dis + h_dis / len(range(0, 60, 4))
            np.save(h_dis_location + h_dis_name, h_aver_dis)
