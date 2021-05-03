import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import sys
import os

plt.style.use("dark_background")

directory = sys.argv[1]
kernels_dir = os.path.join(directory, "kernelsasc")
kernel_files = os.listdir(kernels_dir)

modes = ["S"]
ns = [0]
ls = [32, 43, 62, 76, 97, 132, 167, 202, 225, 254, 291, 320, 342, 382, 415]
kernels = ["beta"]

mpl.rc("lines", linewidth=1)
mpl.rc("axes", grid=True)

i = 0
fig, axs = plt.subplots(1, len(kernels), sharey=True, subplot_kw={"box_aspect": 2.3})
if not isinstance(axs, np.ndarray):
    axs = np.array([axs])
cmap = cm.plasma(np.linspace(0, 1, len(ls)))
for file in sorted(kernel_files):
    mode_type = file.split(".")[0]
    n_order = int(file.split(".")[1])
    l_order = int(file.split(".")[2])
    if mode_type not in modes or n_order not in ns or l_order not in ls:
        continue
    rad, kkappa, kmu, kalpha, kbeta = np.loadtxt(
        os.path.join(kernels_dir, file), unpack=True, skiprows=1
    )
    kernels_dict = {"kappa": kkappa, "mu": kmu, "alpha": kalpha, "beta": kbeta}
    dpth = (6371000 - rad) / 1000
    for kernel, ax in zip(kernels, axs):
        ax.plot(kernels_dict[kernel], dpth, c=cmap[i], label=f"$\ell = ${l_order}")
    i += 1

axs[0].set_ylim([1000, 0])
axs[0].set_ylabel("Depth (km)")

kernel_labels = {
    "kappa": "$K_{\kappa}$",
    "mu": "$K_{\mu}$",
    "alpha": "$K_{\\alpha}$",
    "beta": "$K_{\\beta}$",
}
for kernel, ax in zip(kernels, axs):
    ax.set_xlabel(kernel_labels[kernel])

plt.legend(bbox_to_anchor=(1.05, 1))
plt.tight_layout()
plt.show()
