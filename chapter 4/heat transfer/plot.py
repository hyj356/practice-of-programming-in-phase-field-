import numpy as np
import matplotlib.pyplot as plt

initial_temperature = np.loadtxt('./initial_temp.txt')
timestep_200 = np.loadtxt('./timestep_200.txt')
timestep_400 = np.loadtxt('./timestep_400.txt')
timestep_600 = np.loadtxt('./timestep_600.txt')
# 创建包含4个子图的图形对象，2行2列布局
fig, axs = plt.subplots(2, 2)

# 在第一个子图中绘制初始温度
axs[0, 0].plot(initial_temperature)
axs[0, 0].set_title('initial_temperature')

# 在第二个子图中绘制经过200步迭代之后的温度
axs[0, 1].plot(timestep_200)
axs[0, 1].set_title('timestep_200')

# 在第三个子图中绘制经过400步迭代之后的温度
axs[1, 0].plot(timestep_400)
axs[1, 0].set_title('timestep_400')

# 在第四个子图中绘制经过600步迭代之后的温度
axs[1, 1].plot(timestep_600)
axs[1, 1].set_title('timestep_600')

# 调整子图之间的间距
plt.tight_layout()
plt.savefig("temperature.png", dpi=600)
# 显示图形
plt.show()