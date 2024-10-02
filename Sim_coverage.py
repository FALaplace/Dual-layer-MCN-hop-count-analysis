import matplotlib.pyplot as plt
import numpy as np
import pickle
from Network_Class import Network
import os
import re
from tqdm import trange
import json

CON_NAME_CONFIG = ["TeleSat", "OneWeb", "Kuiper", "Starlink", "Synthetic1", "Synthetic2", "Synthetic3"]
con_name = "TeleSat"
con_file = f"constellation/{con_name}/LLA&LLR.npy"
sat_dict_file = f"constellation/{con_name}/satdict.pkl"
con_info_file = f"constellation/{con_name}/con_info.txt"
alldata = np.load(con_file)
f = open(sat_dict_file, 'rb')
sat_dict = pickle.load(f)
sat_list = list(sat_dict.keys())
f.close()
N1 = Network(sat_list, alldata, con_info_file)


def res():
    W1 = np.load(f"access_sat_record/{con_name}/W1.npy")
    W2 = np.load(f"access_sat_record/{con_name}/W2.npy")
    print(np.max(W1), np.min(W1), np.average(W1))
    print(np.max(W2), np.min(W2), np.average(W2))


def sim():
    l1_re, l2_re = [], []
    N1.GraphInitialization(rho=0.0)
    for ti in trange(240):
        N1.GraphShift(ti)
        ass_info = N1.get_access_sat(ti, 0, 0)
        l1_re.append(len(ass_info["A1"] + ass_info["D1"]))
        l2_re.append(len(ass_info["A2"] + ass_info["D2"]))
    np.save(f"access_sat_record/{con_name}/W1.npy", l1_re)
    np.save(f"access_sat_record/{con_name}/W2.npy", l2_re)


def draw():
    lw = 5.0
    sr = 400
    xticklabel = [r"40$\times$33" + "\n TeleSat", r"27$\times$13" + "\n TeleSat", r"32$\times$72" + "\n OneWeb",
                  r"18$\times$36" + "\n OneWeb", r"34$\times$34" + "\n Kuiper", r"28$\times$28" + "\n Kuiper",
                  r"72$\times$22" + "\n Starlink", r"36$\times$20" + "\n Starlink", r"15$\times$20", r"25$\times$20",
                  r"5$\times$100"]
    with open('access_sat_record/intensity.json', 'r') as file:
        json_data = file.read()
    intensity = json.loads(json_data)
    print(intensity)
    fig: plt.Figure = plt.figure(1, figsize=(22, 11))
    ax: plt.Axes = fig.add_axes([0.08, 0.08, 0.84, 0.84])
    xi = 1
    for cons in CON_NAME_CONFIG:
        i1, i2 = intensity[cons]
        W1, W2 = np.load(f"access_sat_record/{cons}/W1.npy"), np.load(f"access_sat_record/{cons}/W2.npy")
        max1, max2 = np.max(W1), np.max(W2)
        min1, min2 = np.min(W1), np.min(W2)
        av1, av2 = np.average(W1), np.average(W2)
        if xi == 1:
            ax.scatter([xi, xi + 1], [max1, max2], marker="^", c="brown", s=sr,
                       label="maximum observation satellite number")
            ax.scatter([xi, xi + 1], [av1, av2], marker="s", c="brown", s=sr,
                       label="average observation satellite number")
            ax.scatter([xi, xi + 1], [i1, i2], marker="*", c="blue", s=sr,
                       label=r"intensity parameter $\delta$")
            ax.scatter([xi, xi + 1], [min1, min2], marker="v", c="brown", s=sr,
                       label="minimum observation satellite number")
            ax.plot([xi, xi], [min1, max1], linestyle="--", c="c", linewidth=lw)
            ax.plot([xi + 1, xi + 1], [min2, max2], linestyle="--", c="c", linewidth=lw)
            xi += 2
        else:
            if cons in ["Synthetic1", "Synthetic2", "Synthetic3"]:
                ax.plot([xi, xi], [min2, max2], linestyle="--", c="c", linewidth=lw)
                ax.scatter([xi], [max2], marker="^", c="brown", s=sr)
                ax.scatter([xi], [av2], marker="s", c="brown", s=sr)
                ax.scatter([xi], [i2], marker="*", c="blue", s=sr)
                ax.scatter([xi], [min2], marker="v", c="brown", s=sr)
                xi += 1
            else:
                ax.plot([xi, xi], [min1, max1], linestyle="--", c="c", linewidth=lw)
                ax.plot([xi + 1, xi + 1], [min2, max2], linestyle="--", c="c", linewidth=lw)
                ax.scatter([xi, xi + 1], [max1, max2], marker="^", c="brown", s=sr)
                ax.scatter([xi, xi + 1], [av1, av2], marker="s", c="brown", s=sr)
                ax.scatter([xi, xi + 1], [i1, i2], marker="*", c="blue", s=sr)
                ax.scatter([xi, xi + 1], [min1, min2], marker="v", c="brown", s=sr)
                xi += 2

        pass
    ax.legend(fontsize=25)
    ax.set_xticks(list(range(1, len(xticklabel) + 1)))
    ax.set_yticks([0, 2, 4, 6, 8, 10, 12, 14, 16])
    ax.set_xticklabels(xticklabel, fontsize=25)
    ax.set_yticklabels([0, 2, 4, 6, 8, 10, 12, 14, 16], fontsize=25)
    ax.set_ylabel(r"Observation satellite number at latitude 0$^{\circ}$", fontsize=25)
    plt.show()
    # plt.savefig("fig3.10.jpg")
    pass


def draw_asssat_num():
    synthetic3_w2: np.ndarray = np.load("access_sat_record/Synthetic3/W2.npy")
    synthetic2_w2: np.ndarray = np.load("access_sat_record/Synthetic2/W2.npy")
    print(synthetic2_w2.shape, synthetic3_w2.shape)
    t_seq = list(range(synthetic3_w2.size))
    fig: plt.Figure = plt.figure(1, figsize=(16, 8))
    ax: plt.Axes = fig.add_axes([0.08, 0.1, 0.84, 0.8])
    ax.hlines(y=1.606, xmin=0, xmax=240, linestyles="--", colors="k", label=r"intensity parameter $\delta$")
    ax.plot(t_seq, synthetic2_w2[:240], marker="o", c="red", label=r"$W_1$ of Synthetic2: 25$\times$20")
    ax.plot(t_seq, synthetic3_w2, marker="o", c="blue", label=r"$W_1$ of Synthetic3: 5$\times$100")
    ax.tick_params(labelsize=15)
    ax.set_ylabel(r"Observation satellite number at latitude 0$^{\circ}$", fontsize=18)
    ax.set_xlabel("time steps", fontsize=18)
    ax.legend(fontsize=18)
    ax.set_xlim(0, 240)
    plt.show()
    pass


if __name__ == "__main__":
    # sim()
    # res()
    # draw_asssat_num()
    draw()
    pass
