# -*- coding: utf-8 -*-
"""Post processing sensitivities."""
import os
import numpy as np
from scipy.stats import gmean as gm
import matplotlib.pyplot as plt

# plotting style
plt.style.use("ggplot")
# increase fontsize of plots, prevent type 3 fonts in pdf output
plt.rcParams.update({"font.size": 16, "pdf.fonttype": 42, "ps.fonttype": 42})

# file extension of the saved plots
file_ext = ".pdf"

# paths
here = os.path.abspath(os.path.dirname(__file__))
results = os.path.normpath(os.path.join(here, "..", "results", "01_estimate"))

# well names at each site
wells = {"hh": ["p05", "p40", "p42", "p44"], "lw": ["B2", "B3", "B4", "B5"]}


def plot(site, root=None):
    """Plot sensitivities."""
    root = results if root is None else root
    site = site.lower()
    ptests = wells[site]
    estim = os.path.join(root, "all_" + site, "sens_estimate.txt")
    testdict = {}
    print(estim)
    # load data from file
    testdict["all"] = np.loadtxt(estim)

    for p in ptests:
        estim = os.path.join(root, p, "sens_estimate.txt")
        print(estim)
        # load data from file
        testdict[p] = np.loadtxt(estim)

    print(testdict)
    lin = np.ones(2)

    # keys = ["all", "mean"]+ptests#+["mean"]#, "geo-mean"]
    keys = ["all"] + ptests

    varnames = [r"$T^G$", r"$\sigma^2$", r"$\ell$", r"$S$"]

    mean = np.zeros(len(varnames))
    gmean = np.zeros(len(varnames))
    temp = np.zeros((len(ptests), len(varnames)))
    for i, res in enumerate(ptests):
        temp[i, :] = testdict[res]

    for i in range(len(mean)):
        mean[i] = np.mean(temp[:, i])
        gmean[i] = gm(temp[:, i])

    testdict["mean"] = mean
    testdict["geo-mean"] = gmean

    fig = plt.figure(dpi=75, figsize=[7.5, 4])
    for j, var in enumerate(varnames):

        ax = fig.add_subplot(1, len(varnames), j + 1)
        for i, res in enumerate(keys):
            if i < 1:
                ax.plot(
                    testdict[res][j] * lin,
                    label=res,
                    linewidth=3,
                    color="k",
                    alpha=0.7,
                    dashes=max(i - 1, 0) * (1, 1) + np.sign(i) * (3, 1),
                )
            else:
                ax.bar(
                    0.125 + 0.25 * (i - 1),
                    testdict[res][j],
                    0.25,
                    label=res,
                    color="C" + str(i + 1),
                    alpha=0.8,
                    linewidth=4 - np.sign(i),
                )
        ax.set_xlabel(var, fontsize=18)
        ax.set_ylim([-0.05, 1.05])
        ax.locator_params(axis="y", nbins=6)
        if j == 0:
            ax.set_ylabel("total sensitivity", fontsize=16)
        else:
            ax.set_yticklabels([])
        ax.set_xticks([])

    legend = ax.get_legend_handles_labels()
    fig.legend(*legend, loc="lower center", ncol=6, bbox_to_anchor=(0.5, 0.05))
    fig.tight_layout()
    fig.subplots_adjust(wspace=0.1, bottom=0.3)
    fig.show()
    fig.savefig(
        os.path.join(root, site.upper() + "_sensitivity" + file_ext),
        bbox_inches="tight",
    )


plot("HH")
plot("LW")
