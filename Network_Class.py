"""
定义卫星网络的类
"""
import numpy as np
import networkx as nx
import re
from tqdm import tqdm, trange
from numpy.linalg import norm

R_e = 6371
w_e = 2 * np.pi / 86400
d2r = lambda x: x * np.pi / 180
r2d = lambda x: x * 180 / np.pi
lla2xyz = lambda a, b, r: np.array([[r * np.cos(a) * np.cos(b)], [r * np.cos(a) * np.sin(b)], [r * np.sin(a)]])


class Network:
    def __init__(self, Nodelist, NodeData, Con_info):
        super(Network, self).__init__()
        self.T = NodeData.shape[1]
        self.t = 0
        self.G = nx.Graph()  # 空图
        self.NodeNum: int = len(Nodelist)
        self.Nodelist: list = Nodelist
        self.G.add_nodes_from(Nodelist)
        self.__alldata = NodeData[:, :, 1:]
        self.EdgeNum = 0
        self.layout: dict = {}
        self.N_list: list = []
        self.M_list: list = []
        self.Node_dict: dict = {self.Nodelist[k]: k for k in range(self.NodeNum)}
        self.IOL_list: list = []
        self.ILL_list: list = []
        self.All_Edge_list: list = []
        self.degree_distribution: list
        for k in range(0, self.NodeNum):
            now_node = self.Nodelist[k]
            now_layer_num = list(map(int, re.findall(r"-?\d+\.?\d*", now_node)))[0]
            now_plane_num = list(map(int, re.findall(r"-?\d+\.?\d*", now_node)))[1]
            now_sat_num = list(map(int, re.findall(r"-?\d+\.?\d*", now_node)))[2]
            next_k = (k + 1) % self.NodeNum
            next_sat = self.Nodelist[next_k]
            next_layer_num = list(map(int, re.findall(r"-?\d+\.?\d*", next_sat)))[0]
            if next_layer_num != now_layer_num:
                self.N_list.append(now_plane_num + 1)
                self.M_list.append(now_sat_num + 1)
            elif (now_node == self.Nodelist[-1]) and (not self.N_list):
                self.N_list.append(now_plane_num + 1)
                self.M_list.append(now_sat_num + 1)
        self.LayerNum = len(self.N_list)
        f = open(Con_info, 'r')
        self.con_info = f.readlines()
        f.close()
        self.__alt_list = list(map(float, self.con_info[2].split(',')))
        self.__inc_list = list(map(float, self.con_info[3].split(',')))
        self.__F_list = list(map(float, self.con_info[4].split(',')))
        self.__semi_cone_angle = float(self.con_info[5])
        pass

    @staticmethod
    def __CalDistance_LLR(info1, info2):
        V1 = np.array(info1[0:3])
        V2 = np.array(info2[0:3])
        distance = np.linalg.norm(V1 - V2)
        return distance

    def __CreateISLs(self, t):
        for k in range(0, self.NodeNum):
            now_node = self.Nodelist[k]
            now_layer_num = list(map(int, re.findall(r"-?\d+\.?\d*", now_node)))[0]
            now_plane_num = list(map(int, re.findall(r"-?\d+\.?\d*", now_node)))[1]
            now_sat_num = list(map(int, re.findall(r"-?\d+\.?\d*", now_node)))[2]
            N = self.N_list[now_layer_num - 1]
            M = self.M_list[now_layer_num - 1]
            now_node_xyz = self.__alldata[k, t, 0:6]
            Dis_List = []
            for ki in range(M):
                node_k = 'Sat_' + str(now_layer_num) + '_' + str((now_plane_num - 1) % N) + '_' + str(ki)
                k_index = self.Node_dict[node_k]
                node_k_data = self.__alldata[k_index, t, 0:6]
                dis_k = self.__CalDistance_LLR(now_node_xyz, node_k_data)
                Dis_List.append(dis_k)
            opt_index = np.argmin(Dis_List)
            left_node = 'Sat_' + str(now_layer_num) + '_' + str((now_plane_num - 1) % N) + '_' + str(opt_index)
            if nx.degree(self.G, left_node) >= 4:
                opt_index = sorted(enumerate(Dis_List), key=lambda x: x[1])[1][0]
                left_node = 'Sat_' + str(now_layer_num) + '_' + str((now_plane_num - 1) % N) + '_' + str(opt_index)
            backward_node = 'Sat_' + str(now_layer_num) + '_' + str(now_plane_num) + '_' + str((now_sat_num + 1) % M)
            b_index = self.Node_dict[backward_node]
            b_xyz = self.__alldata[b_index, t, 0:6]
            b_dis = self.__CalDistance_LLR(now_node_xyz, b_xyz)
            l_dis = min(Dis_List)
            self.G.add_edge(now_node, backward_node, weight=b_dis)
            self.All_Edge_list.append((now_node, backward_node, b_dis))
            if l_dis < 5000:
                self.G.add_edge(now_node, left_node, weight=l_dis)
                self.IOL_list.append((now_node, left_node, l_dis))
                self.All_Edge_list.append((now_node, left_node, l_dis))

    def __UpdateIOLs(self, t):
        assert bool(len(self.IOL_list)), "please initialization network before shift network"
        old_new_dict = {}
        for i in range(len(self.IOL_list)):
            u, v, w_old = self.IOL_list[i]
            uk = self.Node_dict[u]
            vk = self.Node_dict[v]
            llr_u = self.__alldata[uk, t]
            llr_v = self.__alldata[vk, t]
            w = self.__CalDistance_LLR(llr_u, llr_v)
            self.IOL_list[i] = (u, v, w)
            old_new_dict[(u, v, w_old)] = (u, v, w)
            pass
        self.All_Edge_list = [old_new_dict[i] if i in old_new_dict.keys() else i for i in self.All_Edge_list]
        self.G.add_weighted_edges_from(self.IOL_list)

    def __CreateILLs(self, t, max_link_num, rho):
        np.random.seed(45)
        assert self.LayerNum > 1, "this constellation cannot create inter layers link, " \
                                  "please check the layer number if bigger than 1. "
        M_list = [sum(self.M_list[:-1]), self.M_list[-1]]
        N_list = [sum(self.N_list[:-1]), self.N_list[-1]]
        wait_list = self.Nodelist[0: M_list[0] * N_list[0]]
        for s in range(M_list[0] * N_list[0], self.NodeNum):
            if np.random.rand() < rho:
                Node_Name = self.Nodelist[s]
                Pos_sat = self.__alldata[s, t][0:3]  # x,y,z
                # Vel_sat = self.__alldata[s, t][3:6]  # vx, vy, vz
                # Mon_sat = np.cross(Pos_sat, Vel_sat)  # (x, y, z)×(vx,vy,vz)
                # R_sat = np.linalg.norm(Pos_sat)  # 当前卫星轨道半径
                # V_sat = np.linalg.norm(Vel_sat)  # 当前卫星速度大小
                # H_sat = np.linalg.norm(Mon_sat)  # 当前卫星的角动量大小
                # d1 = 2 * np.pi * (R_e + 550) / self.N_list[-1]  # 轨道间距
                # d2 = 2 * np.pi * (R_e + 550) / self.M_list[-1]  # 轨道内卫星间距
                spare_dis_list = []
                spare_sat_dict = {}
                for tar_sat in wait_list:
                    q = self.Node_dict[tar_sat]
                    dis_Vec = np.array(self.__alldata[q, t][0:3]) - np.array(Pos_sat)  # 相对位置矢量
                    dis = np.linalg.norm(dis_Vec)
                    # ang1 = r2d(np.arccos(abs(np.dot(dis_Vec, Pos_sat)) / (R_sat * dis)))  # <dis_Vec, R_sat>
                    # ang2 = r2d(np.arccos(abs(np.dot(dis_Vec, Vel_sat)) / (V_sat * dis)))  # <dis_Vec, V_sat>
                    # ang3 = r2d(np.arccos(abs(np.dot(dis_Vec, Mon_sat)) / (H_sat * dis)))  # <dis_Vec, H_sat>
                    # flag1 = (ang2 < 10) or (ang3 < 10) or (ang1 < 56)  # 角度约束
                    flag2 = dis < 3000  # range constraint
                    flag3 = nx.degree(self.G, tar_sat) < 5 and nx.degree(self.G, Node_Name) < 5  # degree constraint
                    flag4 = self.__alldata[q, t, 9] * self.__alldata[s, t, 9] > 0  # same direction constraint
                    if flag2 and flag3 and flag4:
                        spare_dis_list.append(dis)
                        spare_sat_dict[dis] = tar_sat
                if bool(len(spare_dis_list)):
                    spare_dis_list = sorted(spare_sat_dict)
                    counter = 0
                    for dis_i in spare_dis_list:
                        sat_i = spare_sat_dict[dis_i]
                        self.ILL_list.append((Node_Name, sat_i, dis_i))
                        self.All_Edge_list.append((Node_Name, sat_i, dis_i))
                        wait_list.remove(sat_i)
                        counter += 1
                        if counter >= max_link_num:
                            break
        self.G.add_weighted_edges_from(self.ILL_list)

    def GraphInitialization(self, ml=1, rho=1.0):
        self.t = 0
        self.G.remove_edges_from(self.All_Edge_list)
        self.ILL_list = []
        self.IOL_list = []
        self.All_Edge_list = []
        self.__CreateISLs(0)
        try:
            self.__CreateILLs(0, ml, rho)
        except AssertionError as e:
            if self.LayerNum == 1:
                print("Warning: Single Layer Constellation Dose Not Need to Create Inter Layer Links!")
            else:
                print(e)

    def GraphShift(self, t, ml=1, rho=1.0):
        # self.__UpdateIOLs(t)
        # self.G.remove_edges_from(self.ILL_list)
        # for li in self.ILL_list:
        #     self.All_Edge_list.remove(li)
        # self.ILL_list.clear()
        # try:
        #     self.__CreateILLs(t, ml, rho)
        # except AssertionError as e:
        #     if self.LayerNum == 1:
        #         print("Warning: Single Layer Constellation Dose Not Need to Create Inter Layer Links!")
        #     else:
        #         print(e)
        self.t = t

    def hop_count_cal(self, t, pos1, pos2):
        assert t == self.t, "the Graph is at time{0}, not time{1}, " \
                            "please shift Graph to time{2} first.".format(self.t, t, t)
        acc_sat_dict1 = self.get_access_sat(t, pos1[0], pos1[1])
        acc_sat_dict2 = self.get_access_sat(t, pos2[0], pos2[1])
        hop_count_list = []
        if not bool(acc_sat_dict1.keys()):
            print("Warning: cannot find access satellite for location({0},{1})".format(pos1[0], pos1[1]))
            return np.nan
        if not bool(acc_sat_dict2.keys()):
            print("Warning: cannot find access satellite for location({0},{1})".format(pos2[0], pos2[1]))
            return np.nan
        for dis1 in acc_sat_dict1.keys():
            sat1 = acc_sat_dict1[dis1]
            for dis2 in acc_sat_dict2.keys():
                sat2 = acc_sat_dict2[dis2]
                if nx.has_path(self.G, sat1, sat2):
                    PATH = nx.shortest_path(self.G, sat1, sat2)
                    # PATH_LEN = nx.shortest_path_length(self.G, sat1, sat2, weight="weight")
                    PATH_LEN = len(PATH) - 1
                    hop_count_list.append(PATH_LEN)
        if len(hop_count_list) <= 3:
            hop_count = np.average(hop_count_list)
        else:
            # hop_count_list = sorted(hop_count_list)[0:3]
            hop_count_list = sorted(hop_count_list)[0:4]
            hop_count = np.average(hop_count_list)
        return hop_count

    def get_access_sat(self, t, lat, lon, R=R_e):
        assert t == self.t, "the Graph is at time {0}, not time {1}, " \
                            "please shift Graph to time {2} first.".format(self.t, t, t)
        lon = d2r(lon) + w_e * t * 60
        lat = d2r(lat)
        vec = lla2xyz(lat, lon, R).reshape([3, ])
        acc_sat_dict = {"A1": [], "A2": [], "D1": [], "D2": []}
        asc_sat_dis1, dec_sat_dis1, asc_sat_dis2, dec_sat_dis2 = 2703, 2703, 2703, 2703
        for k in range(self.NodeNum):
            sat_xyz = self.__alldata[k, t, :3]
            sat_lat_rate = self.__alldata[k, t, 9]
            dis_vec = sat_xyz - vec
            dis = np.round(norm(dis_vec), 2)
            cosang1 = abs(np.dot(dis_vec, sat_xyz)) / (norm(sat_xyz) * norm(dis_vec))
            cosang1 = np.clip(cosang1, 0, 1)
            ang1 = np.arccos(cosang1)
            ang1 = r2d(ang1)
            flag1 = ang1 < self.__semi_cone_angle * 0.95
            flag2 = dis < np.sqrt(norm(sat_xyz) ** 2 - R ** 2)
            if flag1 and flag2:
                if k < self.N_list[0] * self.M_list[0]:
                    # on the first layer
                    if dis < asc_sat_dis1 and sat_lat_rate > 0:
                        # asc_sat_dis1 = dis
                        acc_sat_dict["A1"].append(self.Nodelist[k])
                    elif dis < dec_sat_dis1 and sat_lat_rate < 0:
                        # dec_sat_dis1 = dis
                        acc_sat_dict["D1"].append(self.Nodelist[k])
                else:
                    if dis < asc_sat_dis2 and sat_lat_rate > 0:
                        # asc_sat_dis2 = dis
                        acc_sat_dict["A2"].append(self.Nodelist[k])
                    elif dis < dec_sat_dis2 and sat_lat_rate < 0:
                        # dec_sat_dis2 = dis
                        acc_sat_dict["D2"].append(self.Nodelist[k])
                    pass
                pass
            pass
        return acc_sat_dict

    def hop_count_cal_single_pos(self, t, lat):
        """
        Coverage Area:
            latitude range: [57, -57]
            longitude range: [-180, 170]
        Implementation method:
            latitude: np.linspace(57, -57, 39)
            longitude: np.linspace(-177, 177, 119)
        Return:
            Global Hop Count Distribution H: np.array, size=[39, 119]
        """
        global_lat_list = np.round(np.linspace(-57, 57, 39))
        global_lon_list = np.round(np.linspace(-177, 177, 119))
        global_hop_distribution = np.zeros([39, 119])
        pos1 = [lat, 0.0]
        i = 0
        for lat2 in tqdm(global_lat_list):
            j = 0
            for lon2 in global_lon_list:
                pos2 = [lat2, lon2]
                h = self.hop_count_cal(t, pos1, pos2)
                global_hop_distribution[i, j] = h
                j += 1
            i += 1
        return global_hop_distribution

    @property
    def alt_list(self):
        return self.__alt_list

    @property
    def inc_list(self):
        return self.__inc_list

    @property
    def F_list(self):
        return self.__F_list

    @property
    def semi_cone_angle(self):
        return self.__semi_cone_angle

    def Latency(self, t, pos1, pos2):
        acc_sat1 = self.get_access_sat(t, pos1[0], pos1[1])
        acc_sat2 = self.get_access_sat(t, pos2[0], pos2[1])
        latency = 999
        for dis1 in acc_sat1.keys():
            sat1 = acc_sat1[dis1]
            for dis2 in acc_sat2.keys():
                sat2 = acc_sat2[dis2]
                if nx.has_path(self.G, sat1, sat2):
                    PATH_LEN = nx.dijkstra_path_length(self.G, sat1, sat2)
                    PATH = nx.dijkstra_path(self.G, sat1, sat2)
                    now_latency = PATH_LEN / 300 + len(PATH) - 1
                    latency = min(now_latency, latency)
                    pass
                pass
            pass
        return np.round(latency, 2)


def CSV_to_NPY(csv):
    t = csv.shape[0]
    csv_keys = csv.keys()[2:]
    sat_num = len(csv_keys)
    sat_lla = np.zeros([t, sat_num, 6])
    for ti in trange(t):
        for ki in range(sat_num):
            now_sat = csv_keys[ki]
            sat_info = csv.iloc[ti][now_sat][1:-1]
            info = list(map(float, sat_info.split(',')))
            sat_lla[ti, ki] = info
    return sat_lla


class SingleLayerNetwork:
    def __init__(self, Nodelist, NodeData):
        super(SingleLayerNetwork, self).__init__()
        self.T = NodeData.shape[1]
        self.t = 0
        self.G = nx.Graph()
        self.NodeNum: int = len(Nodelist)
        self.Nodelist: list = Nodelist
        self.G.add_nodes_from(Nodelist)
        self.__alldata = NodeData[:, :, 1:]
        self.EdgeNum = 0
        self.layout: dict = {}
        self.Node_dict: dict = {self.Nodelist[k]: k for k in range(self.NodeNum)}
        self.degree_distribution: list
        last_node = self.Nodelist[-1]
        self.N = list(map(int, re.findall(r"-?\d+\.?\d*", last_node)))[1] + 1
        self.M = list(map(int, re.findall(r"-?\d+\.?\d*", last_node)))[2] + 1
        self.All_Edge_list: list = []
        self.IOL_list: list = []
        self.RISL_list: list = []
        self.Layer1_n_list: list = []
        self.Layer2_n_list: list = []
        self.scheme: int = 0

    def GraphInitialization(self, rho=0.0, scheme=2):
        self.t = 0
        self.G.remove_edges_from(self.All_Edge_list)
        self.RISL_list = []
        self.IOL_list = []
        self.All_Edge_list = []
        self.scheme = 0
        self.Layer1_n_list = []
        self.Layer2_n_list = []
        self.__CreateISLs(0, scheme)
        self.__CreateRandomISLs(0, rho)
        pass

    def GraphShift(self, t, rho):
        self.__UpdateIOLs(t)
        self.G.remove_edges_from(self.RISL_list)
        for li in self.RISL_list:
            self.All_Edge_list.remove(li)
        self.RISL_list.clear()
        self.__CreateRandomISLs(t, rho)
        self.t = t
        pass

    @staticmethod
    def __CalDistance_LLR(info1, info2):
        V1 = np.array(info1[0:3])
        V2 = np.array(info2[0:3])
        distance = np.linalg.norm(V1 - V2)
        return distance

    def __CreateISLs(self, t, gm):
        assert self.N % gm == 0, "this scheme cannot divide the Starlink to two layers constellation."
        self.scheme = gm
        gn = self.N // gm
        layer_1_n_list = [gm * i for i in range(gn)]
        layer_2_n_list = list(set(range(self.N)) ^ set(layer_1_n_list))
        self.Layer1_n_list = layer_1_n_list
        self.Layer2_n_list = layer_2_n_list
        # create layer1 ISLs
        for turn1 in range(gn):
            n1 = layer_1_n_list[turn1]
            left_n1 = layer_1_n_list[turn1 - 1]
            for m1 in range(self.M):
                now_sat = "Sat_1" + '_' + str(n1) + '_' + str(m1)
                now_k = self.Node_dict[now_sat]
                now_xyz = self.__alldata[now_k, t, 0:6]
                Dis_List_left = []
                for ki in range(self.M):
                    node_ki = 'Sat_1' + '_' + str(left_n1) + '_' + str(ki)
                    ki_index = self.Node_dict[node_ki]
                    node_ki_xyz = self.__alldata[ki_index, t, 0:6]
                    dis_k = self.__CalDistance_LLR(now_xyz, node_ki_xyz)
                    if nx.degree(self.G, node_ki) < 4:
                        Dis_List_left.append(dis_k)
                    else:
                        Dis_List_left.append(9999)
                left_m1 = np.argmin(Dis_List_left)  # the optimized m
                left_sat = 'Sat_1' + '_' + str(left_n1) + '_' + str(left_m1)  # left neighbor name
                left_dis = min(Dis_List_left)  # distance between now sat and the left sat
                backward_sat = 'Sat_1' + '_' + str(n1) + '_' + str((m1 + 1) % self.M)  # backward neighbor name
                backward_k = self.Node_dict[backward_sat]  # the index of the backward neighbor
                backward_xyz = self.__alldata[backward_k, t, 0:6]  # related information of the backward neighbor
                backward_dis = self.__CalDistance_LLR(now_xyz, backward_xyz)
                self.G.add_edge(now_sat, backward_sat, weight=backward_dis)
                self.All_Edge_list.append((now_sat, backward_sat, backward_dis))
                self.G.add_edge(now_sat, left_sat, weight=left_dis)
                self.IOL_list.append((now_sat, left_sat, left_dis))
                self.All_Edge_list.append((now_sat, left_sat, left_dis))
        # create layer2 ISLs
        for turn2 in range(self.N - gn):
            n2 = layer_2_n_list[turn2]
            left_n2 = layer_2_n_list[turn2 - 1]
            for m2 in range(self.M):
                now_sat = "Sat_1" + '_' + str(n2) + '_' + str(m2)
                now_k = self.Node_dict[now_sat]
                now_xyz = self.__alldata[now_k, t, 0:6]
                Dis_List_left = []
                for ki in range(self.M):
                    node_ki = 'Sat_1' + '_' + str(left_n2) + '_' + str(ki)
                    ki_index = self.Node_dict[node_ki]
                    node_ki_xyz = self.__alldata[ki_index, t, 0:6]
                    dis_k = self.__CalDistance_LLR(now_xyz, node_ki_xyz)
                    if nx.degree(self.G, node_ki) < 4:
                        Dis_List_left.append(dis_k)
                    else:
                        Dis_List_left.append(9999)
                    pass
                left_m2 = np.argmin(Dis_List_left)  # the optimized m
                left_sat = 'Sat_1' + '_' + str(left_n2) + '_' + str(left_m2)  # left neighbor name
                left_dis = min(Dis_List_left)  # distance between now sat and the left sat
                backward_sat = 'Sat_1' + '_' + str(n2) + '_' + str((m2 + 1) % self.M)  # backward neighbor name
                backward_k = self.Node_dict[backward_sat]  # the index of the backward neighbor
                backward_xyz = self.__alldata[backward_k, t, 0:6]  # related information of the backward neighbor
                backward_dis = self.__CalDistance_LLR(now_xyz, backward_xyz)
                self.G.add_edge(now_sat, backward_sat, weight=backward_dis)
                self.All_Edge_list.append((now_sat, backward_sat, backward_dis))
                self.G.add_edge(now_sat, left_sat, weight=left_dis)
                self.IOL_list.append((now_sat, left_sat, left_dis))
                self.All_Edge_list.append((now_sat, left_sat, left_dis))
        return self.IOL_list

    def __UpdateIOLs(self, t):
        assert bool(len(self.IOL_list)), "please initialization network before shift network"
        old_new_dict = {}
        for i in range(len(self.IOL_list)):
            u, v, w_old = self.IOL_list[i]
            uk = self.Node_dict[u]
            vk = self.Node_dict[v]
            llr_u = self.__alldata[uk, t]
            llr_v = self.__alldata[vk, t]
            w = self.__CalDistance_LLR(llr_u, llr_v)
            self.IOL_list[i] = (u, v, w)
            old_new_dict[(u, v, w_old)] = (u, v, w)
            pass
        self.All_Edge_list = [old_new_dict[i] if i in old_new_dict.keys() else i for i in self.All_Edge_list]
        self.G.add_weighted_edges_from(self.IOL_list)
        return self.IOL_list

    def __CreateRandomISLs(self, t, rho):
        np.random.seed(45)
        Graph = self.G.copy()
        for n in self.Layer1_n_list:
            for m in range(self.M):
                if np.random.random() < rho:
                    now_sat = "Sat_1" + "_" + str(n) + "_" + str(m)
                    now_sat_k = self.Node_dict[now_sat]
                    now_sat_xyz = self.__alldata[now_sat_k, t, 0:6]
                    # delete origin left and right neighbors
                    now_neighbors = list(Graph.neighbors(now_sat))
                    right_sat_origin = "Sat_1" + "_" + str((n + self.scheme) % self.N)  # old left neighbors
                    left_sat_origin = "Sat_1" + "_" + str((n - self.scheme) % self.N)  # old right neighbors
                    for neighbor in now_neighbors:
                        if "Sat_1" + "_" + str((n - self.scheme) % self.N) in neighbor:
                            left_sat_origin = neighbor
                            break
                    for neighbor in now_neighbors:
                        if "Sat_1" + "_" + str((n + self.scheme) % self.N) in neighbor:
                            right_sat_origin = neighbor
                            break
                    if self.G.has_edge(now_sat, left_sat_origin):
                        self.G.remove_edge(now_sat, left_sat_origin)
                    if self.G.has_edge(now_sat, right_sat_origin):
                        self.G.remove_edge(now_sat, right_sat_origin)
                    # find new left and right neighbors
                    Dis_List_left = []
                    for ki in range(self.M):
                        node_k = 'Sat_1' + '_' + str((n - 1) % self.N) + '_' + str(ki)
                        k_index = self.Node_dict[node_k]
                        node_k_data = self.__alldata[k_index, t, 0:6]
                        dis_k = self.__CalDistance_LLR(now_sat_xyz, node_k_data)
                        Dis_List_left.append(dis_k)
                        pass
                    opt_index = np.argmin(Dis_List_left)
                    left_sat = "Sat_1" + "_" + str((n - 1) % self.N) + "_" + str(opt_index)  # new left neighbors
                    right_sat = "Sat_1" + "_" + str((n + 1) % self.N)  # new right neighbors
                    for neighbor in list(Graph.neighbors(left_sat)):
                        if "Sat_1" + "_" + str((n + 1) % self.N) in neighbor:
                            right_sat = neighbor
                            break
                    if right_sat == 'Sat_1_1':
                        print("failed!")
                        pass
                    if self.G.has_edge(left_sat, right_sat):
                        self.G.remove_edge(left_sat, right_sat)
                    # create new links
                    left_dis = min(Dis_List_left)
                    right_k = self.Node_dict[right_sat]
                    right_xyz = self.__alldata[right_k, t, 0:6]
                    right_dis = self.__CalDistance_LLR(now_sat_xyz, right_xyz)
                    self.G.add_edge(now_sat, left_sat, weight=left_dis)
                    self.G.add_edge(now_sat, right_sat, weight=right_dis)
                    self.All_Edge_list.append((now_sat, left_sat, left_dis))
                    self.All_Edge_list.append((now_sat, right_sat, right_dis))
                    self.RISL_list.append((now_sat, left_sat, left_dis))
                    self.RISL_list.append((now_sat, right_sat, right_dis))
                    pass
                pass
            pass
        return self.RISL_list

    def hop_count_cal(self, t, pos1, pos2):
        assert t == self.t, "the Graph is at time{0}, not time{1}, " \
                            "please shift Graph to time{2} first.".format(self.t, t, t)
        acc_sat_dict1 = self.get_access_sat(t, pos1[0], pos1[1])
        acc_sat_dict2 = self.get_access_sat(t, pos2[0], pos2[1])
        hop_count_list = []
        if not bool(acc_sat_dict1.keys()):
            print("Warning: cannot find access satellite for location({0},{1})".format(pos1[0], pos1[1]))
            return np.nan
        if not bool(acc_sat_dict2.keys()):
            print("Warning: cannot find access satellite for location({0},{1})".format(pos2[0], pos2[1]))
            return np.nan
        for dis1 in acc_sat_dict1.keys():
            sat1 = acc_sat_dict1[dis1]
            for dis2 in acc_sat_dict2.keys():
                sat2 = acc_sat_dict2[dis2]
                if nx.has_path(self.G, sat1, sat2):
                    PATH = nx.shortest_path(self.G, sat1, sat2, weight='weight')
                    # PATH_LEN = nx.shortest_path_length(self.G, sat1, sat2, weight="weight")
                    PATH_LEN = len(PATH) - 1
                    hop_count_list.append(PATH_LEN)
        if not hop_count_list:
            avr_hop_count = np.nan
        elif len(hop_count_list) <= 5:
            avr_hop_count = np.average(hop_count_list)
        else:
            hop_count_list = sorted(hop_count_list)[0:5]
            avr_hop_count = np.average(hop_count_list)
        return avr_hop_count

    def get_access_sat(self, t, lat, lon, R=R_e):
        assert t == self.t, "the Graph is at time{0}, not time{1}, " \
                            "please shift Graph to time{2} first.".format(self.t, t, t)
        lon = d2r(lon) + w_e * t * 60
        lat = d2r(lat)
        vec = lla2xyz(lat, lon, R).reshape([3, ])
        acc_sat_dict = {}
        for k in range(self.NodeNum):
            sat_xyz = self.__alldata[k, t, 0:3]
            dis_vec = sat_xyz - vec
            dis = np.round(norm(dis_vec), 2)
            ang1 = np.arccos(abs(np.dot(dis_vec, sat_xyz)) / (norm(sat_xyz) * dis))
            ang1 = r2d(ang1)
            flag1 = ang1 < 56
            flag2 = dis < np.sqrt(norm(sat_xyz) ** 2 - R ** 2)
            if flag1 and flag2:
                acc_sat_dict[dis] = self.Nodelist[k]
                pass
            pass
        return acc_sat_dict

    def hop_count_cal_single_pos(self, t, lat):
        """
        Coverage Area:
            latitude range: [53, -53]
            longitude range: [-180, 170]
        Implementation method:
            latitude: np.linspace(53.2, -53.2, 12)
            longitude: np.linspace(-175, 175, 36) - 5.0
        Return:
            Hop Count Distribution H: np.array, size=[12, 36]
        """
        hop_count_distribution = np.zeros([12, 36])
        # hop_count_distribution = np.zeros([5, 36])
        pos1 = [lat, 0.0]
        i = 0
        tbar = tqdm(np.linspace(53, -53, 12))
        tbar.set_description("Time={0}s".format(self.t * 60))
        for lat2 in tbar:
            j = 0
            for lon2 in np.linspace(-175, 175, 36):
                pos2 = [lat2, lon2 - 5.0]
                h = self.hop_count_cal(t, pos1, pos2)
                hop_count_distribution[i, j] = h
                j += 1
            i += 1
        return hop_count_distribution
