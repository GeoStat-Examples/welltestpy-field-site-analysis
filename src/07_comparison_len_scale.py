# -*- coding: utf-8 -*-
"""Plot length scale comparison."""

import os
import numpy as np
from matplotlib import pyplot as plt
from anaflow import theis, ext_theis_2d

# plotting style
plt.style.use("ggplot")
# increase fontsize of plots, prevent type 3 fonts in pdf output
plt.rcParams.update({"font.size": 16, "pdf.fonttype": 42, "ps.fonttype": 42})

# file extension of the saved plots
file_ext = ".pdf"

# paths
here = os.path.abspath(os.path.dirname(__file__))
results = os.path.normpath(os.path.join(here, "..", "results", "02_len_scale"))
os.makedirs(results, exist_ok=True)


def dashes(i=1, max_d=12, space=1):
    """Dashes for matplotlib."""
    return i * [space, space] + [max_d - 2 * i * space, space]


time_labels = ["10 s", "10 min", "10 h"]
time = [10, 600, 36000]  # 10s, 10min, 10h
rad = np.geomspace(0.02, 4)  # radius from the pumping well in [0, 4]
var = 0.5  # variance of log-transmissivity
len_scales = [5, 10, 20, 40]  # correlation length of log-transmissivity
TG = 1e-4  # the geometric mean transmissivity
TH = TG * np.exp(-var / 2.0)  # the harmonic mean transmissivity
S = 1e-4  # storativity
rate = -1e-4  # pumping rate

tg_str = "Theis($T_G$)"
th_str = "Theis($T_H$)"
ex_str = "$h^{{\\mathrm{{eff}}}}(\\ell={})$"

head_TG = theis(time, rad, S, TG, rate)
head_TH = theis(time, rad, S, TH, rate)
head_ef = []
for i, len_scale in enumerate(len_scales):
    head_ef.append(ext_theis_2d(time, rad, S, TG, var, len_scale, rate))
time_ticks = []

fig, ax = plt.subplots(figsize=[7.5, 6])

for i, step in enumerate(time):
    col = "C" + str(1)
    label_TG = tg_str if i == 0 else None
    label_TH = th_str if i == 0 else None
    for j, len_scale in enumerate(len_scales):
        label_ef = ex_str.format(len_scale) if i == 0 else None
        ax.plot(
            rad,
            head_ef[j][i],
            label=label_ef,
            color="k",
            dashes=dashes(j + 1),
            linewidth=2,
            alpha=0.5,
        )
    ax.plot(
        rad, head_TG[i], label=label_TG, color=col, linestyle="--", linewidth=3
    )
    ax.plot(rad, head_TH[i], label=label_TH, color=col, linewidth=3)
    text_v = (head_TG[i][-1] + head_TH[i][-1]) / 2
    ax.annotate(
        time_labels[i],
        xy=(rad[-1], text_v),
        xytext=(rad[-1] + 0.1, text_v),
        arrowprops=dict(arrowstyle="-", color="black", linewidth=2),
        verticalalignment="center",
    )
ax.set_xlabel("$r$ in [m]")
ax.set_ylabel("$h$ in [m]")
ax.set_xlim([0, rad[-1]])
ax.set_ylim([-2, 0])
ax.locator_params(axis="y", nbins=5)
ax.legend(loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.45))
fig.tight_layout()
fig.subplots_adjust(bottom=0.3)
fig.savefig(
    os.path.join(results, "len_scale_comparison" + file_ext),
    bbox_inches="tight",
    dpi=150,
)
fig.show()
