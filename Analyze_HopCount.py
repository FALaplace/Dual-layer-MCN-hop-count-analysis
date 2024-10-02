import numpy as np
from numpy import arcsin, arctan, sin, cos, tan, arccos, sign, exp
from math import pi, sqrt, factorial
import re
from visualization_analyze import draw_error_lat1, draw_error_rho
import os

Re = 6371.0
d2r = lambda x: x * np.pi / 180
r2d = lambda x: x * 180 / np.pi
distance = lambda lat1, lat2, lon1, lon2: Re * arccos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2))
pp = lambda sigma, n: exp(-sigma) * sigma ** n / factorial(n)
di = d2r(3)


def Hop_Count_Analyze_Origin(fai1, fai2, lam, m, n, i, f):
    """
    i:轨道倾角(deg)
    f:相位因子
    """
    fai1 = np.random.choice([0.002, -0.002]) if fai1 == 0 else fai1
    # f = 0
    fai2 = np.clip(fai2, -d2r(i), d2r(i))
    delta_f = f / (m * n) * 2 * pi
    u1 = arcsin(sin(fai1) / sin(d2r(i)))
    u2 = arcsin(sin(fai2) / sin(d2r(i)))
    l1 = arctan(np.cos(d2r(i)) * np.tan(u1))
    l2 = arctan(np.cos(d2r(i)) * np.tan(u2))
    H_list = []
    for mode in ['A2A', 'A2D', 'D2A', 'D2D']:
        if mode == 'A2A':
            u11 = u1
            u12 = u2
            l11 = l1
            l12 = l2
        elif mode == 'A2D':
            u11 = u1
            u12 = np.sign(fai2) * pi - u2
            l11 = l1
            l12 = pi + arctan(np.cos(d2r(i)) * np.tan(u12))
        elif mode == 'D2A':
            u12 = u2
            u11 = np.sign(fai1) * pi - u1
            l12 = l2
            l11 = pi + arctan(np.cos(d2r(i)) * np.tan(u11))
        else:
            u11 = np.sign(fai1) * pi - u1
            u12 = np.sign(fai2) * pi - u2
            l11 = pi + arctan(np.cos(d2r(i)) * np.tan(u11))
            l12 = pi + arctan(np.cos(d2r(i)) * np.tan(u12))
        delta_L0 = lam + l11 - l12
        if abs(delta_L0) > pi:
            delta_L0 = np.mod(delta_L0 + pi, 2 * pi) - pi
        H1 = delta_L0 * n / (2 * pi)
        delta_U0 = u12 - u11 - H1 * delta_f
        if abs(delta_U0) > pi:
            delta_U0 = np.mod(delta_U0 + pi, 2 * pi) - pi
        H2 = delta_U0 * m / (2 * pi)
        H_list.append(round(abs(H1) + abs(H2)))
    return np.min(np.round(H_list)), H_list


def Coverage_Cal(N, M, i, lat, Ro):
    f = lambda x: M * N / (2 * pi ** 2) * (sqrt(2) * cos(x) / sqrt(cos(2 * x) - cos(2 * (i + di))))
    Sat_Num = N * M
    dlat = Ro / Re
    dlon = Ro / Re / cos(lat)
    if abs(lat) < i:
        lat1 = min(lat + dlat, i)
        lat2 = max(lat - dlat, -i)
        N = (2 * dlon) * abs(lat2 - lat1) * (f(lat1) + 4 * f((lat2 + lat1) / 2) + f(lat2)) / 6 * pi / 4
        return N, 0
    elif abs(lat) < i + dlat:
        lat2 = dlat - (abs(lat) - i)
        u1 = pi / 2
        u2 = arcsin(sin(d2r(i - lat2)) / sin(d2r(i)))
        b = d2r(abs(lat) - i) * Re
        R_ = sqrt(Ro ** 2 - b ** 2)
        dlon = r2d(R_ / Re) / cos(d2r(i))
        N = abs(u1 - u2) / (2 * pi) * Sat_Num * 2 * dlon / 360 * pi / 4
        return N, 1
    else:
        return 0, 1


def Probability_Model(H_list: list, p1_list, p2_list, h2, h3):
    x_list = []
    p_list = []
    x_list.append(min(H_list))
    p_list.append(p1_list[0] * p2_list[0])
    x_list.append(min(H_list[0::2] + [h3]))
    p_list.append(p1_list[0] * p2_list[1])
    x_list.append(min(H_list[1::2] + [h3]))
    p_list.append(p1_list[0] * p2_list[2])
    x_list.append(min(min(h2), h3))
    p_list.append(p1_list[0] * p2_list[3])
    x_list.append(min(H_list[0:2] + [h3]))
    p_list.append(p1_list[1] * p2_list[0])
    x_list.append(min(min(h2), h3, H_list[0]))
    p_list.append(p1_list[1] * p2_list[1])
    x_list.append(min(H_list[1], min(h2), h3))
    p_list.append(p1_list[1] * p2_list[2])
    x_list.append(min(h3, h2[1]))
    p_list.append(p1_list[1] * p2_list[3])
    x_list.append(min(H_list[2:] + [h3]))
    p_list.append(p1_list[2] * p2_list[0])
    x_list.append(min(H_list[2], min(h2), h3))
    p_list.append(p1_list[2] * p2_list[1])
    x_list.append(min(H_list[3], min(h2), h3))
    p_list.append(p1_list[2] * p2_list[2])
    x_list.append(min(h2[3], h3))
    p_list.append(p1_list[2] * p2_list[3])
    x_list.append(min(min(h2), h3))
    p_list.append(p1_list[3] * p2_list[0])
    x_list.append(min(h2[0], h3))
    p_list.append(p1_list[3] * p2_list[1])
    x_list.append(min(h2[1], h3))
    p_list.append(p1_list[3] * p2_list[2])
    x_list.append(h3)
    p_list.append(p1_list[3] * p2_list[3])
    return np.dot(x_list, p_list)


def Cal_H12(h1: int, h2: int, r):
    if r == 0:
        return h2
    if h1 == 0:
        return h2 / 2 + 1
    h12_list = [h2 * (1 - k / h1) + k + 1 if h2 * (1 - k / h1) + k + 1 < h2 else h2 for k in range(0, h1 + 1)]
    h12_list = sorted(h12_list)
    p12_list = [r * (1 - r) ** k for k in range(h1)]
    p12_list.append(1 - np.sum(p12_list))
    h12 = np.dot(h12_list, p12_list)
    return h12


def Cal_H22(h1: int, h2: int, r):
    if r == 0:
        return h2
    if h1 < h2 - 2:
        kmax = int(h2 * (1 - 2 / (h2 - h1)))
        h22_list = [h1 * (1 - k / h2) + k + 2 for k in range(kmax + 1)]
        h22_list.append(h2)
        p22_list = [r ** 2 * (k + 1) * (1 - r) ** k for k in range(kmax + 1)]
        p22_list.append(1 - sum(p22_list))
        h22 = np.dot(h22_list, p22_list)
        return h22
    else:
        return h2


def Two_Layer_Hop_Count_Analyze(r, x1, x2_list, y_list, con_name):
    """
    :param r: rho, interlayer links
    :param x1: user1's latitude
    :param x2_list: user2's latitude range
    :param y_list: longitude difference range
    :param con_name:tuple : walker constellation name, the name must be
                            ["Kuiper", "OneWeb", "TeleSat", "Starlink", "Starlink_1", "Starlink_2"]
    :return: analyze results of global hop count distribution
    """
    x1 = x1 if x1 != 0 else np.random.choice([-1, 1]) * 0.002
    con_info_path = f"constellation/{con_name}/con_info.txt"
    f = open(con_info_path, 'r')
    con_info = f.readlines()
    [N2, M2] = list(map(int, con_info[0].split(",")))
    [N1, M1] = list(map(int, con_info[1].split(",")))
    [alt2, alt1] = list(map(float, con_info[2].split(",")))
    [i2, i1] = list(map(float, con_info[3].split(",")))
    [F2, F1] = list(map(int, con_info[4].split(",")))
    semi_cone_angle = float(con_info[5])
    rho1 = r * N1 * M1 / (N2 * M2)
    z = np.zeros([x2_list.size, y_list.size], dtype=int)
    z1 = np.zeros_like(z)
    R1 = tan(d2r(semi_cone_angle)) * alt1
    R2 = tan(d2r(semi_cone_angle)) * alt2
    int1, fg1 = Coverage_Cal(N1, M1, d2r(i1), d2r(x1), R1)
    p1_list = [1 - pp(int1, 1) - pp(int1, 0), pp(int1, 1) / 2, pp(int1, 1) / 2, pp(int1, 0)]
    p1 = min(int1 / 2, 1)
    for i in range(x2_list.size):
        now_lat = d2r(x2_list[i]) if x2_list[i] != 0 else np.random.choice([-1, 1]) * 0.002
        int2, fg2 = Coverage_Cal(N1, M1, d2r(i1), d2r(x2_list[i]), R1)
        p2_list = [1 - pp(int2, 1) - pp(int2, 0), pp(int2, 1) / 2, pp(int2, 1) / 2, pp(int2, 0)]
        p2 = min(int2 / 2, 1)
        for j in range(y_list.size):
            now_lon = d2r(y_list[j])
            h1, h1_list = Hop_Count_Analyze_Origin(d2r(x1), now_lat, now_lon, M1, N1, i1, F1)
            h2, h2_list = Hop_Count_Analyze_Origin(d2r(x1), now_lat, now_lon, M2, N2, i2, F2)
            h1_int = int(h1) if r >= 0.1 else int(np.average((sorted(h1_list)[1::2]) - h1) * (0.1 - r) ** 2 * 100 + h1)
            h2_int = int(h2) if r >= 0.1 else int(np.average((sorted(h2_list)[1::2]) - h2) * (0.1 - r) ** 2 * 100 + h2)
            H12_A1 = Cal_H12(min(h1_list[:2]), min(h2_list[:2]), r)
            H12_D1 = Cal_H12(min(h1_list[2:]), min(h2_list[2:]), r)
            H12_A2 = Cal_H12(min(h1_list[::2]), min(h2_list[::2]), r)
            H12_D2 = Cal_H12(min(h1_list[1::2]), min(h2_list[1::2]), r)
            H12 = [H12_A1, H12_D1, H12_A2, H12_D2]
            H22 = Cal_H22(h1_int, h2_int, rho1)
            H4 = Probability_Model(h1_list, p1_list, p2_list, H12, H22)
            z[i, j] = np.round(H4 * (1.05 + 15 * (0.1 - r) ** 2)) if r <= 0.1 else np.round(H4 * 1.05)
            z1[i, j] = h1
            pass
        pass
    return z, np.average(z1)


def Error_Analyze(z_a, z_r):
    # z_r = z_r[8:-8]
    # z_a = z_a[8:-8]
    # z_r = np.round(z_r)
    e = np.abs(z_r - z_a)
    # e = z_r - z_a
    e = np.array(e, dtype=int)
    # er = e / z_r
    er_sum = np.sum(e) / np.sum(z_r)
    return np.round(er_sum, 4), np.round(np.average(z_r), 4), np.round(np.average(z_a), 4)


if __name__ == "__main__":
    ConName = "Starlink"
    ConInfo = f"constellation/{ConName}/con_info.txt"
    f1 = open(ConInfo, 'r')
    [N2, M2] = list(map(int, f1.readline().split(",")))
    [N1, M1] = list(map(int, f1.readline().split(",")))
    f1.close()
    lat1 = 20
    rho = 0.2
    dir_name = "Hop_dis_twolayers"
    sub_dir_name = f"({N2},{M2})&({N1},{M1})_global_hop_distribution"
    lat2_list = np.round(np.linspace(-57, 57, 39))
    lon_list = np.round(np.linspace(-177, 177, 119))
    rho_list = [0.0, 0.02, 0.04, 0.06, 0.08] + [i / 10 for i in range(1, 11)]
    lat1_list = [0, 10, 20, 30, 40, 50]
    error = []
    ave_hop = []
    ave_hop2 = []
    ave_hop_w1 = []
    filename1 = "filename"
    try:
        # for lat1 in lat1_list:
        for rho in rho_list:
            filename1 = f"lat1={lat1}_rho={rho}.npy"
            Z_real1 = np.load(f"{dir_name}/{sub_dir_name}/{filename1}")
            Z_anal1, avhw1 = Two_Layer_Hop_Count_Analyze(rho, lat1, lat2_list, lon_list, ConName)
            ela1, avh1, avh2 = Error_Analyze(Z_anal1, Z_real1)
            # print(ela1, avh1, avh2)
            error.append(ela1)
            ave_hop.append(avh1)
            ave_hop2.append(avh2)
            ave_hop_w1.append(0.95 * avhw1 + avhw1 * d2r(lat1) ** 2 * 0.15)
    except FileNotFoundError:
        print(f"No Such File {filename1}! The data has not been proposed!")
    else:
        # draw_error_lat1(lat1_list, error, ave_hop, ave_hop2, "Starlink", y4=ave_hop_w1)
        draw_error_rho(rho_list, error, ave_hop, ave_hop2, ConName)
    finally:
        print(np.round(np.average(error), 4), np.round(np.average(ave_hop), 4))
        print(np.round(np.average(error) * np.average(ave_hop), 4))
        pass
