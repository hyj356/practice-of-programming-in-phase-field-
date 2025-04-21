import numpy as np
import matplotlib.pyplot as plt

temp = np.loadtxt("./final_temp.txt")
temp0 = np.loadtxt("./initial_temp.txt")
fig, ax = plt.subplots(1, 2)
ax[0].plot(temp0)
ax[0].set_ylim(0, 1)
ax[1].plot(temp)
ax[1].set_ylim(0, 1)
plt.show()