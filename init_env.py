"""
    Python 首次连接STK（此脚本只运行一次）
    初始化STK-Python连接环境
    此文件运行后，会在目录[Anaconda3]\Lib\site-packages\comtypes\gen下创建以下模块：
        - AgUiApplicationLib.py
        - AgUiCoreLib.py
        - AgSTKGraphicsLib.py
        - AgSTKVgtLib.py
        - STKObjects.py
        - STKUtil.py
        - AgStkGatorLib
"""
import comtypes
from comtypes.client import CreateObject

# Get reference to running STK instance
uiApplication = CreateObject("STK11.Application")
uiApplication.Visible = True
uiApplication.UserControl = True
print('app的类型为：', type(uiApplication))
# Get our IAgStkObjectRoot interface
# 获取Object Model的根对象：IAgStkObjectRoot
# 此接口为Object Model中的最顶层接口，由此接口可创建场景、地面站、卫星等
root = uiApplication.Personality2
print('root的类型为：', type(root))
"""
Note: When 'root=uiApplication.Personality2' is executed, 
the comtypes library automatically creates a gen folder that 
contains STKUtil and STK Objects. After running this at 
least once on your computer, the following two lines should 
be moved before the 'uiApplication=CreateObject("STK11.Application")' 
line for improved performance.  
"""
#  创建Astrogator相关的模块：AgStkGatorLib
comtypes.client.GetModule((comtypes.GUID("{090D317C-31A7-4AF7-89CD-25FE18F4017C}"), 1, 0))

print('python 首次连接STK完成！')
print('STK Object Model API 的python模块已在comtypes\gen目录下创建！')
print('请关闭已打开的STK!')
