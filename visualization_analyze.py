import numpy as np
import matplotlib.pyplot as plt

import seaborn as sns
import re
import os


def draw_error_lat1(x1, y1, y2, y3, con_name, y4=None):
    fig1 = plt.figure(1, figsize=(10, 6))
    axes = fig1.add_subplot(111)
    axes.tick_params(axis="both", labelsize=14)
    axes.plot(x1, y2, marker="s", c="dimgrey", label="Average Hop(Simulation Result)", markersize=13, lw=2.5, alpha=1.0)
    axes.plot(x1, y3, marker="v", c="red", label="Average Hop(Analysis Result)", markersize=13, lw=2.5, alpha=0.6)
    # axes.plot(x1, y4, marker="v", c=sns.xkcd_rgb["hazel"], label=r"W$_1$ Layer Average Hop(Analysis Result)")
    axes.set_yticks([i * 2 for i in range(6)])
    axes.set_ylim([0.0, 10])
    axes.set_xlabel(r"$\varphi_1$", fontsize=15)
    axes.set_ylabel(r"Hop Count", fontsize=15)
    axes.grid(True)
    axes2 = axes.twinx()
    axes2.tick_params(axis="both", labelsize=14)
    axes2.plot(x1, y1, marker="x", c='navy', label="Relative Error", markersize=13, lw=2.5)
    axes2.set_yticks([i * 0.05 for i in range(6)])
    axes2.set_xticks(x1)
    axes2.set_xlim([x1[0], x1[-1]])
    axes2.set_ylabel("Relative Error", fontsize=15)
    fig1.legend(bbox_to_anchor=(0.9, 0.6), fontsize=13)
    # plt.title(f"{con_name} Hop-Count and Error Distribution under User1's Latitude")
    plt.savefig("fig4.1.png")
    # plt.show()
    # plt.close("all")


def draw_error_rho(x1, y1, y2, y3, con_name):
    fig1 = plt.figure(1, figsize=(10, 6))
    axes = fig1.add_subplot(111)
    axes.tick_params(axis="both", labelsize=14.5)
    xmin = 0.0
    xmax = 1.0
    major_ticks = np.linspace(xmin, xmax, 11)
    minor_ticks = x1
    axes.set_xticks(major_ticks)
    axes.set_xticks(minor_ticks, minor=True)
    axes.set_yticks([i * 2 for i in range(11)])
    axes.set_yticks(list(range(20)), minor=True)
    axes.plot(x1, y2, marker="s", c="dimgrey", label="Average Hop(Simulation Result)", markersize=10, lw=2.5, alpha=1.0)
    axes.plot(x1, y3, marker="v", c="red", label="Average Hop(Analysis Result)", markersize=10, lw=2.5, alpha=0.6)
    # axes.plot(x1, y2, marker="", c="grey", label="Average Hop(Simulation Result)")
    # axes.plot(x1, y3, marker="v", c=sns.xkcd_rgb["hazel"], label="Average Hop(Analysis Result)")
    axes.grid(which="major", alpha=0.8)
    axes.grid(which="minor", alpha=0.5, linestyle="-")
    axes.set_ylim([0.0, 15])
    axes.set_xlabel(r"$\rho_1$", fontsize=18)
    axes.set_ylabel(r"Hop Count", fontsize=18)
    axes2 = axes.twinx()
    axes2.tick_params(axis="both", labelsize=14.5)
    axes2.plot(x1, y1, marker="x", c='navy', label="Relative Error", markersize=10, lw=2.5)
    axes2.set_yticks([i * 0.05 for i in range(10)])
    axes.set_xlim([x1[0], x1[-1]])
    axes2.set_ylim([0.0, 0.375])
    axes2.set_ylabel("Relative Error", fontsize=18)
    # plt.title(fr"{con_name} Hop-Count and Error under $\rho$")
    fig1.legend(bbox_to_anchor=(0.9, 0.88), fontsize=14)
    plt.savefig("fig4.2.png")
    # plt.show()
    plt.close("all")


def draw_aver_hop_rho(x, y1, y2=None, y3=None, y4=None, y5=None, con_fig=None):
    import matplotlib as mpl
    mpl.rcParams["xtick.labelsize"] = 13
    mpl.rcParams["ytick.labelsize"] = 13
    plt.figure(1, figsize=(10, 6))
    plt.plot(x, y1, label=con_fig[0], marker="o", markersize=9, c="red", lw=3, ls=(1, (5, 1)))
    plt.plot(x, y2, label=con_fig[1], marker="s", markersize=9, c="green", alpha=0.6, lw=3, ls=(1, (3, 1, 1, 1, 1, 1)))
    plt.plot(x, y3, label=con_fig[2], marker="^", markersize=9, c="orange", alpha=0.6, lw=3, ls=(1, (3, 1, 1, 1)))
    plt.plot(x, y4, label=con_fig[3], marker="v", markersize=9, c="navy", alpha=0.6, lw=3, ls="-.")
    plt.plot(x, y5, label=con_fig[4], marker="d", markersize=9, c="k", alpha=0.6, lw=2)
    plt.yticks([5 + i for i in range(12)])
    plt.xticks([i / 10 for i in range(11)])
    plt.xlim([0, 1])
    plt.legend(fontsize=13)
    plt.grid(True, "major", alpha=0.5)
    plt.xlabel(r"$\rho_1$", fontsize=15)
    plt.ylabel(r"Average $H$", fontsize=15)
    # plt.show()
    plt.savefig("fig4.3.png")


def draw_aver_rho_dis(con_name, lat1, rho):
    save_path = os.path.join("Hop_dis_image", con_name)
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    con_info = f"constellation/{con_name}/con_info.txt"
    f = open(con_info, "r")
    [N1, M1] = list(map(int, f.readline().split(',')))
    [N2, M2] = list(map(int, f.readline().split(',')))
    f.close()
    data_file = f"Hop_dis_twolayers/({N1},{M1})&({N2},{M2})_global_hop_distribution/lat1={lat1}_rho={rho}.npy"
    save_path = f"Hop_dis_image/{con_name}/lat1={lat1}_rho={rho}.png"
    z = np.load(data_file)
    plt.rc("xtick", labelsize=15)
    plt.rc("ytick", labelsize=15)
    fig1 = plt.figure(1, figsize=(20, 11))
    ax = fig1.add_subplot(111)
    img = plt.imread("Earth Surface.jpg")
    # fg, ax = plt.subplots()
    ax.imshow(img, extent=[0, 118, -10, 50])
    sns.set(font_scale=1.5)
    sns.heatmap(z, cmap="plasma", vmax=20, vmin=0, ax=ax, linewidths=0.00, square=True, center=None,
                cbar_kws={"label": "Hop", "ticks": [0, 5, 10, 15, 20], "shrink": 0.6},
                robust=True,
                alpha=0.85)
    plt.ylim([-10, 50])
    lat = [-90 + 10 * i for i in range(19)]
    lon = [180 - 20 * i for i in range(19)]
    shift_lat = lambda xi: 13 * (xi + 90) / 40 - 9.75
    shift_lon = lambda yi: z.shape[1] / 360 * (yi + 180)
    lat_origin = [shift_lat(latk) for latk in lat]
    lon_origin = [shift_lon(lonk) for lonk in lon]
    plt.yticks(lat_origin, lat)
    plt.xticks(lon_origin, lon, rotation=0)
    plt.xlabel(r"$\lambda$(°)", fontsize=20)
    plt.ylabel(r"$\varphi_2$(°)", fontsize=20)
    # plt.savefig(save_path)
    plt.show()
    plt.close("all")


def draw_aver_path_len(x1, y1, y2, y3):
    fig, axes = plt.subplots()
    xmin = 0.0
    xmax = 1.0
    major_ticks_top = np.linspace(xmin, xmax, 11)
    axes.set_xticks(major_ticks_top)
    axes.set_yticks(np.linspace(0, 80, 5))
    axes.set_yticks(np.linspace(0, 80, 9))
    axes.grid(which="major", alpha=0.8)
    axes.grid(which="minor", alpha=0.4)
    plt.plot(x1, y1, marker="o", label="(72,22)&(15,20)", markersize=2)
    plt.plot(x1, y2, marker="o", label="(72,22)&(25,20)", markersize=2)
    plt.plot(x1, y3, marker="o", label="(72,22)&(36,20)", markersize=2)
    plt.legend()
    plt.xlim([xmin - 0.05, xmax + 0.05])
    plt.ylim([5, 80])
    plt.xlabel(r"$\rho$")
    plt.ylabel(r"diameter")
    plt.title("Two-Layer Constellation Network Diameter")
    plt.savefig("Images/network_diameter.png")
    print("pass")
    # plt.show()


def draw_average_hop_optimized(x, y1, y2, y3, y4):
    fig, axes = plt.subplots()
    xmin = x[0]
    xmax = x[-1]
    major_ticks_top = np.linspace(xmin, xmax, 11)
    axes.set_xticks(major_ticks_top)
    axes.set_yticks([5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, y1[0]])
    axes.grid(which="major", alpha=0.6)
    axes.grid(which="minor", alpha=0.3)
    plt.plot(x, y1, linestyle="--", label="Traditional ISLs", linewidth=1, c='r')
    plt.plot(x, y2, marker="^", label="(36,22)&(36,22),s=2", c='b')
    plt.plot(x, y3, marker="^", label="(24,22)&(48,22),s=3", c='r')
    plt.plot(x, y4, marker="^", label="(18,22)&(54,22),s=4", c='g')
    plt.legend(loc="center right")
    plt.xlim([xmin, xmax])
    plt.ylim([5, 10])
    plt.xlabel(r"$\rho$")
    plt.ylabel(r"average $H$")
    plt.title("Starlink-L1")
    plt.text(np.argmin(y2) / 10 - 0.015, np.round(min(y2), 2) + 0.1, np.round(min(y2), 2))
    plt.text(np.argmin(y3) / 10 - 0.015, np.round(min(y3), 2) + 0.1, np.round(min(y3), 2))
    plt.text(np.argmin(y4) / 10 - 0.015, np.round(min(y4), 2) - 0.18, np.round(min(y4), 2))
    plt.savefig("Images/average_hop_2.png")
    # plt.show()


shift_lat = lambda xi: 13 * (xi + 90) / 38 - 11.29
shift_lon = lambda yi: 119 / 360 * (yi + 180)


def draw_rho_global(z, path: str):
    lat1 = re.findall(r"-?\d+\.?\d*", path.split("/")[-1])[-1]
    city_dict = {"0": "Singapore", "10": None, "20": None, "30": "Shanghai", "40": "Beijing", "50": "London"}
    con_dict = {}
    con_params = os.listdir(path.split("/")[0])
    con_names = ["OneWeb", "TeleSat", "(72,22)&(15,20)", "(72,22)&(25,20)", "Starlink"]
    for con_param, con_name in zip(con_params, con_names):
        con_dict[con_param] = con_name
    background = f"EarthSurface_{city_dict[lat1]}.png" if city_dict[lat1] is not None else "EarthSurface_London.png"
    plt.rc("xtick", labelsize=40)
    plt.rc("ytick", labelsize=40)
    fig1 = plt.figure(1, figsize=(25, 10))
    ax1 = fig1.add_axes([0.06, 0.00, 1.0, 1])
    img = plt.imread(background)
    ax1.imshow(img, extent=[0, 119, -11.29, 50.29])
    draw_subpoint_track(con_dict[path.split("/")[1]], ax1)
    sns.set(font_scale=3.0)
    sns.heatmap(z, cmap="plasma", vmax=50, vmin=0, ax=ax1, linewidths=0.00, square=True, center=None,
                cbar_kws={"label": "Reduction(%)", "ticks": [0, 10, 20, 30, 40, 50], "shrink": 0.8, "pad": 0.01},
                robust=True,
                alpha=0.8)

    # set ticks of x-axis and y-axis
    lat = [57, 0, -57]
    ax1.plot(shift_lon(0), shift_lat(float(lat1)), marker="s", markersize=15, c="navy")
    ax1.text(shift_lon(2), shift_lat(float(lat1) - 1.5), "User1", fontsize=35)
    lat_origin = [shift_lat(latk) for latk in lat]
    lat_tick = [f"{latk}°N" if latk > 0 else (f"{abs(latk)}°S" if latk < 0 else "0") for latk in lat]
    plt.yticks(lat_origin, lat_tick)
    plt.xticks([])
    plt.xlabel(r"$\lambda$", fontsize=40)
    plt.ylabel(r"$\varphi_2$", fontsize=40, labelpad=-40)
    plt.ylim([shift_lat(70), shift_lat(-70)])
    plt.gca().invert_yaxis()
    plt.savefig(path)
    # plt.show()
    plt.close("all")


def draw_rho_lat(x: np.ndarray, y: np.ndarray, path):
    lst = re.findall(r"-?\d+\.?\d*", path.split("/")[-1])
    lst = list(map(int, lst))
    if lst[-1] < 0:
        lst[-1] = f"{abs(lst[-1])}°S"
    else:
        lst[-1] = f"{abs(lst[-1])}°N"
    x_new = x[int((x.size - 1) / 2):]
    y_new = np.zeros([y.shape[0], x_new.size])
    for k in range(x_new.size):
        y_new[:, x_new.size - k - 1] = (y[:, k] + y[:, -k - 1]) / 2
    y0, y4 = y_new[0, :], y_new[4, :]
    plt.figure(1, figsize=(30, 16))
    plt.rc("axes", labelsize=20)
    plt.rc("xtick", labelsize=20)
    plt.rc("ytick", labelsize=20)
    dis_count = (y0 - y4) / y0 * 100
    plt.plot(x_new, y_new[0, :], label=r"$\rho=0.0$", linewidth=2.5)
    plt.plot(x_new, y_new[1, :], label=r"$\rho=0.1$", linewidth=2.5)
    plt.plot(x_new, y_new[2, :], label=r"$\rho=0.2$", linewidth=2.5)
    plt.plot(x_new, y_new[3, :], label=r"$\rho=0.3$", linewidth=2.5)
    plt.plot(x_new, y_new[4, :], label=r"$\rho=0.4$", linewidth=2.5)
    plt.plot(x_new, y_new[9, :], label=r"$\rho=0.9$", linewidth=2.5)
    plt.legend(loc="upper right")
    plt.xticks([10 * i for i in range(19)])
    plt.yticks([2 * i for i in range(11)])
    plt.xlim([-1.4, 180])
    plt.ylim([0, 20])
    plt.xlabel(r"$\Delta\lambda$(°)")
    plt.ylabel("H")
    ax2 = plt.twinx()
    plt.bar(x=x_new, height=dis_count, label='Reduction', color='Coral', alpha=0.8, width=2.8)
    plt.yticks([10 * i for i in range(11)])
    plt.ylim([0, 100])
    plt.ylabel(r"Reduction(%)")
    plt.grid(True)
    plt.title(path.split("/")[1] + rf" $\varphi_1$={lst[1]}°N $\varphi_2$={lst[-1]}", fontsize=25)
    plt.savefig(path)
    # plt.show()
    plt.close("all")
    pass


def draw_subpoint_track(con_name, ax: plt.Axes, legend=True):
    sat_data = np.load(f"constellation/{con_name}/LLA&LLR.npy")
    with open(f"constellation/{con_name}/con_info.txt", "r") as f:
        [N, M] = f.readline().replace('\n', '').split(",")
    LLA_sat1: np.ndarray = sat_data[0, :, [7, 8]]
    LLA_sat2: np.ndarray = sat_data[int(N) * int(M), :, [7, 8]]
    Lon_sat1, Lon_sat2 = shift_lon(LLA_sat1[1, :]), shift_lon(LLA_sat2[1, :])
    Lat_sat1, Lat_sat2 = shift_lat(LLA_sat1[0, :]), shift_lat(LLA_sat2[0, :])
    for i in range(1, len(Lon_sat1)):
        if abs(LLA_sat1[1, i] - LLA_sat1[1, i - 1]) > 180:
            Lon_sat1[i - 1] = None
        if abs(LLA_sat2[1, i] - LLA_sat2[1, i - 1]) > 180:
            Lon_sat2[i - 1] = None
    if legend:
        ax.plot(Lon_sat2, Lat_sat2, linestyle="-.", alpha=0.25, linewidth=4.5, label=r"ground track of $W_1$")
        ax.plot(Lon_sat1, Lat_sat1, linestyle=":", alpha=0.6, linewidth=4.5, label=r"ground track of $W_2$", c='coral')
    else:
        ax.plot(Lon_sat2, Lat_sat2, linestyle="-.", alpha=0.25, linewidth=4.5)
        ax.plot(Lon_sat1, Lat_sat1, linestyle=":", alpha=0.6, linewidth=4.5, c='coral')
    if legend:
        ax.legend(loc='upper right', fontsize=35)
    return ax
