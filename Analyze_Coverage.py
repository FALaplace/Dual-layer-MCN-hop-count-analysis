import numpy as np
import matplotlib.pyplot as plt


def draw_LLA(filename: str):
    all_data: np.ndarray = np.load(filename)
    LLA_data = all_data[1584:, 0, [7, 8]]
    background = "EarthSurface_London.png"
    fig1: plt.Figure = plt.figure(figsize=(12, 6))
    ax: plt.Axes = fig1.add_axes([0, 0, 1, 1])
    img = plt.imread(background)
    ax.imshow(img, extent=[-180, 180, -90, 90])
    ax.scatter(LLA_data[:, 1], LLA_data[:, 0], color="maroon", s=20)
    ax.set_xticks([])
    ax.set_yticks([])
    plt.show()


if __name__ == "__main__":
    f = "constellation/Synthetic2/LLA&LLR.npy"
    draw_LLA(f)
