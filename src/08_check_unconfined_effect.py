# -*- coding: utf-8 -*-
"""Checking the effect of unconfined conditions."""
import os
import numpy as np
from matplotlib import pyplot as plt
# plotting style
plt.style.use("ggplot")
# increase fontsize of plots, prevent type 3 fonts in pdf output
plt.rcParams.update({"font.size": 16, "pdf.fonttype": 42, "ps.fonttype": 42})

# file extension of the saved plots
file_exts = [".pdf", ".png"]

# paths
here = os.path.abspath(os.path.dirname(__file__))
result = os.path.normpath(os.path.join(here, "..", "results", "03_unconfined"))
data = os.path.normpath(os.path.join(here, "..", "data"))
files = ["horkheim_p05p44_analysis.txt", "lauswiesen_B3B4_analysis.txt"]
os.makedirs(result, exist_ok=True)

for file in files:
    # loading
    ana = np.loadtxt(os.path.join(data, file))
    dd_time = ana[:, 0][~np.isnan(ana[:, 0])]
    dd_obs = ana[:, 1][~np.isnan(ana[:, 1])]
    dd_dx_time = ana[:, 2][~np.isnan(ana[:, 2])]
    dd_dx_obs = ana[:, 3][~np.isnan(ana[:, 3])]
    th_dx_time = ana[:, 4][~np.isnan(ana[:, 4])]
    th_dx_obs = ana[:, 5][~np.isnan(ana[:, 5])]
    # plotting
    fig, ax = plt.subplots(dpi=75, figsize=[7.5, 5.5])
    ax.scatter(dd_time, dd_obs, marker=".", color="C1", label="drawdown")
    ax.plot(dd_dx_time, dd_dx_obs, label="time derivative")
    ax.plot(th_dx_time, th_dx_obs, color="k", label="theis time derivative")
    ax.set_xscale("symlog", linthresh=1)
    ax.set_yscale("symlog", linthresh=1e-4)
    ax.set_xlim([1, 1e5])
    ax.set_ylim([1e-4 if file.startswith("horkheim") else 1e-3, 1e-1])
    ax.set_xlabel("$t$ in [s]", fontsize=16)
    ax.set_ylabel("$h$ and $dh/dx$ in [m]", fontsize=16)
    ax.legend(loc="upper left")
    fig.tight_layout()
    for file_ext in file_exts:
        fig.savefig(
            os.path.join(result, os.path.splitext(file)[0] + file_ext),
            bbox_inches="tight",
            dpi=300,
        )
    fig.show()
