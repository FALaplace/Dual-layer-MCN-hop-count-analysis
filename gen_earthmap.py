import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

fig1 = plt.figure(1, figsize=(16, 8))
ax1 = fig1.add_axes([0, 0, 1, 1])
ax1.set_axis_off()
city = "London"
longitude = 0

if __name__ == "__main__":
    # 创建Basemap对象
    map1 = Basemap(projection='cyl', resolution='c', lon_0=longitude, ax=ax1)

    # 绘制地图
    map1.drawcoastlines(linewidth=0.8)
    # map1.drawcountries()
    map1.drawmeridians(range(-180, 181, 20), labels=[1, 0, 0, 1], linewidth=0.3)
    map1.drawparallels(range(-90, 91, 20), labels=[1, 0, 0, 1], linewidth=0.3)
    for i in range(-180, 181, 20):
        if i < longitude - 180:
            i_label = i + 360 - 5
            label = f"{abs(i)}°E" if i > 0 else (f"{abs(i)}°W" if i < 0 else "0")
        elif i > longitude + 180:
            i_label = i - 360 - 5
            label = f"{abs(i)}°E" if i > 0 else (f"{abs(i)}°W" if i < 0 else "0")
        else:
            i_label = i - 5
            label = f"{abs(i)}°E" if i > 0 else (f"{abs(i)}°W" if i < 0 else "0")
        if i == 0 or abs(i) == 180:
            i_label = i if i == 0 else i - 5
            label = f"{abs(i)}"
        plt.text(i_label, -65, label, fontsize=16)
    # plt.show()
    plt.savefig(f"EarthSurface_{city}.png")
