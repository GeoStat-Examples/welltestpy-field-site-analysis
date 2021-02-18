# -*- coding: utf-8 -*-
"""Post-processing results."""
import os
import numpy as np
from scipy.stats import gmean as gm
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

# ploting style
plt.style.use("ggplot")
# increase fontsize of plots, prevent type 3 fonts in pdf output
plt.rcParams.update({"font.size": 16, "pdf.fonttype": 42, "ps.fonttype": 42})

# file extension of the saved plots
file_ext = ".png"

# paths
here = os.path.abspath(os.path.dirname(__file__))
results = os.path.normpath(os.path.join(here, "..", "results", "01_estimate"))

# well names at each site
wells = {"hh": ["p05", "p40", "p42", "p44"], "lw": ["B2", "B3", "B4", "B5"]}


def plot(site, root=None):
    """Plot estimation results."""
    root = results if root is None else root
    site = site.lower()
    ptests = wells[site]
    estim = os.path.join(root, "all_" + site, "estimate.txt")
    testdict = {}
    print(estim)
    # load data from file
    testdict["all"] = np.loadtxt(estim)

    for p in ptests:
        estim = os.path.join(root, p, "estimate.txt")
        print(estim)
        # load data from file
        testdict[p] = np.loadtxt(estim)
        # if var is 0 --> len_scale is 0
        if testdict[p][1] < 0.02:
            print(testdict[p][1])
            testdict[p][2] = 0

    lin = np.ones(2)
    keys = ["all", "mean"] + ptests  # +["mean"]#, "geo-mean"]
    varnames = [
        r"$T^G_{all}$",
        r"$\sigma^2_{all}$",
        r"$\ell_{all}$",
        r"$S_{all}$",
    ]
    varunits = [
        r" $\frac{m^2}{s}$",
        r"",
        r" $m$",
        r"",
    ]

    for res in testdict:
        testdict[res][0] = np.exp(testdict[res][0])
        testdict[res][3] = np.exp(testdict[res][3])

    # rescale to common y-range
    max_y = 2.0
    for res in ptests:
        for i in range(4):
            testdict[res][i] /= testdict["all"][i]
            max_y = max(max_y, testdict[res][i])
    allsave = testdict["all"]
    testdict["all"] = np.ones(4)
    # allsave[2] *= np.sqrt(np.pi) * 0.5  # rescale len_scale to integral scale

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

    # decimal precision for ticks
    tick_prec = [3, 1, 1, 3]
    labels = ["{:03." + str(prec) + "f}" for prec in tick_prec]

    fig = plt.figure(dpi=75, figsize=[7.5, 4])
    for j, var in enumerate(varnames):
        ax = fig.add_subplot(1, len(varnames), j + 1)
        for i, res in enumerate(keys):
            if i < 2:
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
                    0.125 + 0.25 * (i - 2),
                    testdict[res][j],
                    0.25,
                    label=res,
                    color="C" + str(i),
                    alpha=0.8,
                    linewidth=4 - np.sign(i),
                )
        ax.set_ylim([-0.1, max_y + 0.1])
        ax.set_xlabel(
            var + " = " + labels[j].format(allsave[j]) + varunits[j],
            fontsize=16,
        )
        ax.set_xticks([])
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
        if j == 0:
            ax.set_ylabel("multiples of $all$-result", fontsize=16)
        else:
            ax.set_yticklabels([])

    legend = ax.get_legend_handles_labels()
    fig.legend(
        *legend,
        loc="lower center",
        ncol=6,
        bbox_to_anchor=(0.5, 0),
        handlelength=1
    )
    fig.tight_layout()
    fig.subplots_adjust(wspace=0.1, bottom=0.3)
    fig.show()
    fig.savefig(
        os.path.join(root, site.upper() + "_results" + file_ext),
        bbox_inches="tight",
    )


plot("HH")
plot("LW")
