import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
from scipy.spatial import Voronoi, voronoi_plot_2d
from math import sqrt, atan, acos

# 定义需要写入的文本文件的读写通道
out = open('./voronoi_out/Voronoi_vertices.out', 'w')
out1 = open('./voronoi_out/plot_1.out', 'w')
out2 = open('./voronoi_out/final_plot.p', 'w')
out3 = open('./voronoi_out/original_points.out', 'w')
out4 = open('./voronoi_out/cell_1.out', 'w')
out5 = open('./voronoi_out/grain_25.inp', 'w')

# 定义初始变量
npoin = 25      # 25个晶粒
xmax, ymax = 32.0, 32.0 # 生成多晶的区域为一个32x32的正方形
x0, y0 = 0.0, 0.0     # 区域起点为(0.0, 0.0)
extra = 6.0
# 定义两个数组存储voronoi中心点的位置
x, y = np.zeros(9*npoin), np.zeros(9*npoin)
xtemp, ytemp = np.zeros(npoin), np.zeros(npoin)
x[0:npoin] = np.random.random(npoin)*xmax
y[0:npoin] = np.random.random(npoin)*ymax
# 通过平移来达成周期性边界条件
for ipoin in range(npoin):
  jpoin = npoin + ipoin
  x[jpoin] = x[ipoin]
  y[jpoin] = y[ipoin] - ymax
for ipoin in range(npoin):
  jpoin = 2*npoin + ipoin
  x[jpoin] = x[ipoin] + xmax
  y[jpoin] = y[ipoin] - ymax
for ipoin in range(npoin):
  jpoin = 3*npoin + ipoin
  x[jpoin] = x[ipoin] + xmax
  y[jpoin] = y[ipoin]
for ipoin in range(npoin):
  jpoin = 4*npoin + ipoin
  x[jpoin] = x[ipoin] + xmax
  y[jpoin] = y[ipoin] + ymax
for ipoin in range(npoin):
  jpoin = 5*npoin + ipoin
  x[jpoin] = x[ipoin]
  y[jpoin] = y[ipoin] + ymax
for ipoin in range(npoin):
  jpoin = 6*npoin + ipoin
  x[jpoin] = x[ipoin] - xmax
  y[jpoin] = y[ipoin] + ymax
for ipoin in range(npoin):
  jpoin = 7*npoin + ipoin
  x[jpoin] = x[ipoin] - xmax
  y[jpoin] = y[ipoin]
for ipoin in range(npoin):
  jpoin = 8*npoin + ipoin
  x[jpoin] = x[ipoin] - xmax
  y[jpoin] = y[ipoin] - ymax
# 将平移之后的点填入points中
points = np.zeros((9*npoin, 2))
points[:, 0] = x
points[:, 1] = y
# 构建voronoi对象
vor = Voronoi(points)
# 绘制voronoi图像
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
fig = voronoi_plot_2d(vor=vor, ax=ax, show_vertices=False)
ax.set_xlim(-30, 60)
ax.set_ylim(-30, 60)
ax.add_patch(Rectangle((x0, y0), width=xmax, height=ymax, 
                       linewidth = 3, edgecolor='r', facecolor='none'))
ax.axvline(x=x0, linewidth = 2, linestyle='--')
ax.axvline(x=xmax, linewidth = 2, linestyle='--')
ax.axhline(y=y0, linewidth = 2, linestyle='--')
ax.axhline(y=ymax, linewidth = 2, linestyle='--')
plt.savefig("./voronoi_out/Voronoi_figure.png", dpi=600)

# 将模拟盒子的信息写入
out4.write(f"{x0:14.6f} {y0:14.6f}\n")
out4.write(f"{xmax:14.6f} {y0:14.6f}\n")
out4.write(f"{xmax:14.6f} {ymax:14.6f}\n")
out4.write(f"{x0:14.6f} {ymax:14.6f}\n")
out4.write(f"{x0:14.6f} {y0:14.6f}\n")

# 获取voronoi多边形的信息
# c: 一个二维数组, 存储了voronoi所有多边形的顶点的xy坐标
# f: 一个二维的列表, 表明voronoi多边形由哪些顶点构成, 比如第一个元素为[0,3,6,8]
#    表明第一个voronoi多边形由数组c里面的第0,3,6,8个点构成
# prs: 一个整数的列表, 构建起了voronoi多边形与其对应顶点之间的对应关系, 具体而言
#      f[prs[1]], 表示第2个voronoi多边形的顶点索引
# pts: 一个包含了所有用来构建voronoi多边形的中心点坐标的列表
c, f, prs, pts = vor.vertices, vor.regions, vor.point_region, vor.points
nvelem = len(f) # 获取voronoi多边形的数量
ncount = 0
lnods = []
polygons = []
# 遍历所有中心点, 将voronoi中心点落在有限区域内的voronoi多边形的顶点全部加入到lnods中
for i,pt in enumerate(pts):
  if pt[0] >= x0-extra and pt[0] <= x0 + xmax+extra and pt[1] >= y0-extra and pt[1] <= y0 + ymax+extra:
    region_id = prs[i]
    temp_polygon = []
    temp_nods = []
    for index in f[region_id]:
      temp_polygon.append([c[index, 0], c[index, 1]])
      temp_nods.append(index)
    polygons.append(temp_polygon)
    lnods.append(temp_nods)
    ncount += 1
# 将结果进行可视化
nnode_ = 0
# 找到voronoi多边形单元中最多有几个顶点
for polygon in polygons:
  if len(polygon) > nnode_:
    nnode_ = len(polygon)
  ax.add_patch(Polygon(polygon))
plt.savefig("./voronoi_out/Voronoi_reference.png", dpi=600)

# 将voronoi单元对应哪个晶粒进行筛选
twopi = 8.0*atan(1.0)
epsilon = 1.0e-4
nelem = len(lnods)  # 获取有几个符合条件的voronoi单元
igrain = np.zeros(nelem, dtype=np.int32)
for isector in range(9):
  for ipoin in range(npoin):
    jpoin = isector*npoin + ipoin
    for ielem in range(nelem):
      theta = 0.0
      nnode = len(lnods[ielem]) # 获取第ielem个voronoi多边形里面有几个顶点
      # 遍历所有顶点
      for inode in range(nnode):
        kk = lnods[ielem][inode]
        xv1, yv1 = c[kk, 0], c[kk, 1]
        jnode = inode + 1
        if inode == nnode - 1:
          jnode = 0
        jj = lnods[ielem][jnode]
        xv2, yv2 = c[jj, 0], c[jj, 1] 
        p2x, p2y = xv1 - x[jpoin], yv1 - y[jpoin]
        p1x, p1y = xv2 - x[jpoin], yv2 - y[jpoin]
        x1, x2 = sqrt(p1x*p1x + p1y*p1y), sqrt(p2x*p2x + p2y*p2y)
        if x1*x2 <= epsilon:
          theta = twopi
        else:
          tx1 = (p1x*p2x + p1y*p2y) / (x1*x2)
          if tx1 >= 1.0:
            tx1 = 0.99999999
          theta += acos(tx1)
      if abs(theta - twopi) < epsilon:
        igrain[ielem] = ipoin
# 将结果输入到inp文件中以用于后续的仿真
nn1 = len(c)  # 获取voronoi多边形的所有顶点
out5.write(f'{nn1} {nnode_} {ncount} {npoin}\n')
# 将voronoi多边形所有顶点的坐标写入文件out5中
for i in range(nn1):
  out5.write(f'{i+1} {c[i, 0]} {c[i, 1]}\n')
for i, polygon in enumerate(lnods):
  out5.write(f'{i+1} ')
  for j in range(nnode_):
    if j <= len(lnods[i])-1:
      out5.write(f'{lnods[i][j]+1} ')
    else:
      out5.write('0 ')
  out5.write(f' {igrain[i]+1}\n')
# 关闭所有文件
out.close()
out1.close()
out2.close()
out3.close()
out4.close()
out5.close()
# print(c[307])
# print(c[308])
# print(c[309])
# print(c[310])