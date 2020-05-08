# -*- coding: utf-8 -*-
"""Post-processing results."""
import os
import glob
import numpy as np
from scipy.stats import gmean as gm
import matplotlib.pyplot as plt

# ploting style
plt.style.use("ggplot")
plt.rcParams.update({"font.size": 12})

# file extension of the saved plots
file_ext = ".png"

# paths
here = os.path.abspath(os.path.dirname(__file__))
results = os.path.normpath(os.path.join(here, "..", "results", "01_estimate"))


def plot(site, root=None):
    """Plot estimation results."""
    if root is None:
        root = results
    if site == "HH":
        ptests = ["p05", "p40", "p42", "p44"]
        estfile = glob.glob(
            os.path.join(root, "all_hh", "*[0-9]_estimate.txt")
        )
    else:
        ptests = ["B2", "B3", "B4", "B5"]
        estfile = glob.glob(
            os.path.join(root, "all_lw", "*[0-9]_estimate.txt")
        )
    testdict = {}
    # get the latest estimation-file
    estfile.sort()
    assert estfile, "No results found"
    estfile = estfile[-1]
    print(estfile)
    # load data from file
    testdict["all"] = np.loadtxt(estfile)

    for p in ptests:
        # get the latest estimation-file
        estfile = glob.glob(os.path.join(root, p, "*[0-9]_estimate.txt"))
        # sort by date
        estfile.sort()
        # take the newest one
        estfile = estfile[-1]
        print(estfile)
        # load data from file
        testdict[p] = np.loadtxt(estfile)
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
        r" $\left[\frac{m^2}{s}\right]$",
        r" $\left[-\right]$",
        r" $\left[m\right]$",
        r" $\left[-\right]$",
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

    fig = plt.figure(dpi=75, figsize=[12, 4])
    for j, var in enumerate(varnames):
        if j == 0:
            ax1 = fig.add_subplot(1, len(varnames), j + 1)
            for i, res in enumerate(keys):
                if i < 2:
                    ax1.plot(
                        testdict[res][j] * lin,
                        label=res,
                        linewidth=3,
                        color="k",
                        alpha=0.7,
                        dashes=max(i - 1, 0) * (1, 1) + np.sign(i) * (3, 1),
                    )
                else:
                    ax1.bar(
                        0.125 + 0.25 * (i - 2),
                        testdict[res][j],
                        0.25,
                        label=res,
                        color="C" + str(i),
                        alpha=0.8,
                        linewidth=4 - np.sign(i),
                    )
            ax1.set_xlabel(
                var + " = {:05.3f}".format(allsave[j]) + varunits[j],
                fontsize=16,
            )
            ax1.set_ylabel("deviation", fontsize=16)
            ax1.set_ylim([-0.1, max_y + 0.1])
            ticks = ax1.get_yticks() * allsave[j]
            ax1.set_yticklabels(np.round(ticks, 3))
        else:
            ax = fig.add_subplot(1, len(varnames), j + 1)  # , sharey=ax1)
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
            ax.set_xlabel(
                var + " = {:05.3f}".format(allsave[j]) + varunits[j],
                fontsize=16,
            )
            ax.set_ylim([-0.1, max_y + 0.1])
            ticks = ax.get_yticks() * allsave[j]
            ax.set_yticklabels(np.round(ticks, 3))

        plt.xticks([], [])
        if j == 3:
            ax.legend(loc="upper left", handlelength=3, bbox_to_anchor=(1, 1))
    fig.subplots_adjust(wspace=0.35, left=0.08)
    plt.tight_layout()
    plt.show()
    fig.savefig(
        os.path.join(root, site + "_results" + file_ext), bbox_inches="tight"
    )


plot("HH")
plot("LW")
