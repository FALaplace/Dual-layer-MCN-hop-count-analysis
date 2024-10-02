import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from visualization_analyze import draw_subpoint_track

shift_lat = lambda xi: 13 * (xi + 90) / 38 - 11.29
shift_lon = lambda yi: 119 / 360 * (yi + 180)


def draw_hop(con_name, lat1, rho):
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
    # set background image
    city_dict = {0: "Singapore", 10: None, 20: None, 30: 'Shanghai', 40: "Beijing", 50: "London"}
    background = f"EarthSurface_{city_dict[lat1]}.png" if city_dict[lat1] is not None else "EarthSurface_London.png"
    plt.rc("xtick", labelsize=40)
    plt.rc("ytick", labelsize=40)
    fig1 = plt.figure(1, figsize=(25, 10))
    ax = fig1.add_axes([0.06, 0.00, 1.0, 1])
    img = plt.imread(background)
    ax.imshow(img, extent=[0, 119, -11.29, 50.29])
    draw_subpoint_track(con_name, ax)
    sns.set(font_scale=3.5)
    sns.heatmap(z, cmap="YlGnBu", vmax=20, vmin=0, ax=ax, linewidths=0.00, square=True, center=None,
                cbar_kws={"label": "Hop", "ticks": [0, 5, 10, 15, 20], "shrink": 0.8, "pad": 0.01},
                robust=True,
                alpha=0.85)
    lat = [57, 0, -57]
    ax.plot(shift_lon(0), shift_lat(float(lat1)), marker="s", markersize=15, c="navy")
    ax.text(shift_lon(2), shift_lat(float(lat1) - 1.5), "User1", fontsize=35)
    lat_origin = [shift_lat(latk) for latk in lat]
    lat_tick = [f"{latk}°N" if latk > 0 else (f"{abs(latk)}°S" if latk < 0 else "0") for latk in lat]
    plt.yticks(lat_origin, lat_tick)
    plt.xticks([])
    plt.xlabel(r"$\lambda$", fontsize=40)
    plt.ylabel(r"$\varphi_2$", fontsize=40, labelpad=-40)
    plt.ylim([shift_lat(-70), shift_lat(70)])
    plt.savefig(save_path)
    plt.close("all")


if __name__ == "__main__":
    for lat1 in [0, 10, 20, 30, 40, 50]:
        draw_hop("Starlink", lat1, 0.0)
        draw_hop("Starlink", lat1, 0.4)
