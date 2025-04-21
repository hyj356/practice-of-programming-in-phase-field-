import numpy as np
import matplotlib.pyplot as plt

initial = np.loadtxt('./temp.txt')
fig, axs = plt.subplots(1, 1)
axs.imshow(initial, vmin=0.0, vmax=1.0)
# 显示图形
plt.show()