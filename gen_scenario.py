"""
create double layers constellation
"""
import pickle
from comtypes.gen import STKObjects
from comtypes.gen import STKUtil
from comtypes.client import CreateObject
import time
from tqdm import trange
import numpy as np
import re
import os

# 卫星walker星座的参数需要用列表(list)数据类型
# 卫星名称命名为——Sat_layernum_obritnum_satnum
#         例子：Sat_2_0_21表示第二轨道层第0个轨道面第21颗卫星
#         取值范围:
#           layernum: 1,2,3...
#           obritnum: 0,1,2...(N_i - 1)
#           satnum: 0,1,2...(M_i - 1)
constellation_name = "Synthetic3"
start_time = time.time()
N_list = [72, 5]
M_list = [22, 100]
h_list = [550, 580]
F_list = [12, 1]
a_list = [53.2, 60]
use_engine = True
read_scenario = False
scenario_name = "MCN"
LLA_LLR_DATA_FILEPATH = f"constellation/{constellation_name}/LLA&LLR.npy"
SAT_INFO_FILEPATH = f"constellation/{constellation_name}/satdict.pkl"
CON_INFO_FILEPATH = f"constellation/{constellation_name}/con_info.txt"


def CreateSingleLayerWalkerConstellation(N, M, h, F, a, name):
    if not re.findall(r'\d', name):
        name = name + '_1'
    constellation = scen1.Children.New(STKObjects.eConstellation, name)
    constellation2 = constellation.QueryInterface(STKObjects.IAgConstellation)
    # Insert the constellation of Satellites
    tbar = trange(N)
    tbar.set_description(f'Create {name} Constellation')
    for n in tbar:
        AoP = np.mod(np.random.normal(scale=120 / N), 360)
        for m in range(M):
            satellite = scen1.Children.New(STKObjects.eSatellite, f"{name}_{n}_{m}")
            satellite2 = satellite.QueryInterface(STKObjects.IAgSatellite)
            satellite2.SetPropagatorType(STKObjects.ePropagatorTwoBody)
            # Set initial state
            twoBodyPropagator = satellite2.Propagator.QueryInterface(STKObjects.IAgVePropagatorTwoBody)
            keplarian = twoBodyPropagator.InitialState.Representation.ConvertTo(
                STKUtil.eOrbitStateClassical).QueryInterface(STKObjects.IAgOrbitStateClassical)
            keplarian.SizeShapeType = STKObjects.eSizeShapeAltitude  # 定义远地点高度
            keplarian.SizeShape.QueryInterface(STKObjects.IAgClassicalSizeShapeAltitude).ApogeeAltitude = h  # 远地点高度
            keplarian.SizeShape.QueryInterface(STKObjects.IAgClassicalSizeShapeAltitude).PerigeeAltitude = h  # 近地点高度
            keplarian.Orientation.Inclination = a  # 轨道倾角
            keplarian.Orientation.ArgOfPerigee = AoP  # 近地点角距
            keplarian.Orientation.AscNodeType = STKObjects.eAscNodeRAAN
            RAAN = 360 * n / N  # 升交点赤经
            keplarian.Orientation.AscNode.QueryInterface(STKObjects.IAgOrientationAscNodeRAAN).Value = RAAN
            keplarian.LocationType = STKObjects.eLocationTrueAnomaly
            trueAnomaly = np.mod(360 * m / M + 360 * n * F / (M * N), 360)  # 真近角
            keplarian.Location.QueryInterface(STKObjects.IAgClassicalLocationTrueAnomaly).Value = trueAnomaly
            # Propagate
            satellite2.Propagator.QueryInterface(STKObjects.IAgVePropagatorTwoBody).InitialState.Representation.Assign(
                keplarian)
            satellite2.Propagator.QueryInterface(STKObjects.IAgVePropagatorTwoBody).Propagate()
            # Add to constellation object
            constellation2.Objects.AddObject(satellite)  # 向星座加入该卫星
            pass
        pass
    return constellation


def CreateMultilayerWalkerConstellation(N, M, h, F, a, name):
    layernum = len(N)
    for l in range(layernum):
        CreateSingleLayerWalkerConstellation(N[l], M[l], h[l], F[l], a[l], f"{name}_{l + 1}")
    pass


# 导出地心惯性系下卫星坐标和速度
def get_all_data(satlist):
    """
    Create Numpy Array to store the satellites' information
    ALL_DATA.shape = (satellite_number, time, 13)
    ALL_DATA include data: [t ,x, y, z, vx, vy, vz, lat, lon, alt, lat_rate, lon_rate, alt_rate] total 13 parameters
    SAT_DICT = {satellite_name: satellite_ID}
    satellite_ID_list = [0, 1, 2, ..., satellite_number-1]
    """
    elems1 = ["x", "y", "z", "Velocity x", "Velocity y", "Velocity z"]
    elems2 = ['Lat', 'Lon', 'Alt', 'Lat Rate', 'Lon Rate', 'Alt Rate']
    timedp1 = satlist[0].DataProviders.GetDataPrvTimeVarFromPath("LLA State//Fixed")
    time_info = timedp1.ExecElements(scen2.StartTime, scen2.StopTime, 60, ["Time"]).DataSets
    time_list = list(time_info.GetDataSetByName('Time').GetValues())
    tbar = trange(satlist.Count)
    tbar.set_description('Get satellite data')
    ALL_DATA = np.zeros([satlist.Count, len(time_list), 13])
    ALL_DATA[:, :, 0] = np.round(time_list)
    SAT_DICT = {}
    for i in tbar:
        sat_name = satlist[i].InstanceName
        satdp1 = satlist[i].DataProviders.GetDataPrvTimeVarFromPath("Points(Inertial)/Center")
        satdp2 = satlist[i].DataProviders.GetDataPrvTimeVarFromPath("LLA State/Fixed")
        xyz_info = satdp1.ExecElements(scen2.StartTime, scen2.StopTime, 60, elems1).DataSets
        LLA_info = satdp2.ExecElements(scen2.StartTime, scen2.StopTime, 60, elems2).DataSets
        for item in range(6):
            ALL_DATA[i, :, item + 1] = np.round(xyz_info.GetDataSetByName(elems1[item]).GetValues(), 4)
        for item in range(6):
            ALL_DATA[i, :, item + 7] = np.round(LLA_info.GetDataSetByName(elems2[item]).GetValues(), 4)
        SAT_DICT[sat_name] = i
    return ALL_DATA, SAT_DICT


if __name__ == '__main__':
    if use_engine:
        print("Launching STK Engine...")
        stkxApp = CreateObject("STKX11.Application")
        # Disable graphics. The NoGraphics property must be set to true before the root object is created.
        stkxApp.NoGraphics = True
        # Create root object
        root = CreateObject('AgStkObjects11.AgStkObjectRoot')
    else:
        print("Launching STK GUI...")
        uiApp = CreateObject("STK11.Application")
        uiApp.Visible = True
        uiApp.UserControl = True
        # Get root object
        root = uiApp.Personality2
    # 设置时间日期格式
    # 创建新场景
    print("正在创建场景........")
    if not read_scenario:
        root.NewScenario(scenario_name)
        print("场景[{}]已经创建成功".format(scenario_name))
    scen1 = root.CurrentScenario
    scen2 = scen1.QueryInterface(STKObjects.IAgScenario)
    scen2.StartTime = '26 Jan 2023 00:00:00.000'
    scen2.StopTime = '26 Jan 2023 04:00:00.000'
    root.UnitPreferences.Item('DateFormat').SetCurrentUnit('EpSec')
    # todo!!!! 星座名称不能含有数字!
    CreateMultilayerWalkerConstellation(N_list, M_list, h_list, F_list, a_list, 'Sat')
    sat_list = root.CurrentScenario.Children.GetElements(STKObjects.eSatellite)
    all_data, sat_dict = get_all_data(sat_list)
    np.save(LLA_LLR_DATA_FILEPATH, all_data)
    f = open(SAT_INFO_FILEPATH, 'wb')
    pickle.dump(sat_dict, f)
    f.close()
    f = open(CON_INFO_FILEPATH, "w")
    f.write(f"{N_list[0]},{M_list[0]}\n")
    f.write(f"{N_list[1]},{M_list[1]}\n")
    f.write(f"{h_list[0]},{h_list[1]}\n")
    f.write(f"{a_list[0]},{a_list[1]}\n")
    f.write(f"{F_list[0]},{F_list[1]}\n")
    f.close()
    print('finished')
    # get_LLA_data(sat_list, LLA_data_file)
