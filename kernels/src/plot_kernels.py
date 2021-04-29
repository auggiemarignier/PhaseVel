import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import sys
import os

directory = sys.argv[1]
kernels_dir = os.path.join(directory, "kernelsasc")
kernel_files = os.listdir(kernels_dir)

modes = ["S"]
ns = [0]
ls = [32, 43, 62, 76, 97, 132, 167, 202, 225, 254, 291, 320, 342, 382, 415]

mpl.rc("lines", linewidth=1)
mpl.rc("axes", grid=True)

i = 0
fig, axs = plt.subplots(1, 4, sharey=True, figsize=(10, 5))
cmap = cm.plasma(np.linspace(0, 1, len(ls)))
for file in kernel_files:
    mode_type = file.split(".")[0]
    n_order = int(file.split(".")[1])
    l_order = int(file.split(".")[2])
    if mode_type not in modes or n_order not in ns or l_order not in ls:
        continue
    rad, kkappa, kmu, kalpha, kbeta = np.loadtxt(
        os.path.join(kernels_dir, file), usecols=np.arange(5), unpack=True, skiprows=1
    )
    dpth = (6371000 - rad) / 1000
    axs[0].plot(kkappa, dpth, c=cmap[i], label=f"$\ell = ${l_order}")
    axs[1].plot(kmu, dpth, c=cmap[i], label=f"$\ell = ${l_order}")
    axs[2].plot(kalpha, dpth, c=cmap[i], label=f"$\ell = ${l_order}")
    axs[3].plot(kbeta, dpth, c=cmap[i], label=f"$\ell = ${l_order}")
    i += 1

ax = axs[0]
ax.set_ylim([1000, 0])
ax.set_ylabel("Depth (km)")
ax.set_xlabel("$K_{\kappa}$")

axs[1].set_xlabel("$K_{\mu}$")
axs[2].set_xlabel("$K_{\\alpha}$")
axs[3].set_xlabel("$K_{\\beta}$")

plt.legend(bbox_to_anchor=(1.05, 1))
plt.tight_layout()
plt.show()
