import numpy as np
import matplotlib.pyplot as plt

initial = np.loadtxt('./data/con_initial.txt')
eta1 = np.loadtxt('./data/eta_grain1.txt')
eta2 = np.loadtxt('./data/eta_grain2.txt')
fig, axs = plt.subplots(1, 3)
axs[0].imshow(initial, vmin=0.0, vmax=1.0)
axs[1].imshow(eta1, vmin=0.0, vmax=1.0)
axs[2].imshow(eta2, vmin=0.0, vmax=1.0)
# 显示图形
plt.show()