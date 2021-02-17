# -*- coding: utf-8 -*-
"""Post processing sensitivities."""
import os
import glob
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


def plot(site, root=None):
    """Plot sensitivities."""
    if root is None:
        root = results
    if site == "HH":
        ptests = ["p05", "p40", "p42", "p44"]
        estfile = glob.glob(os.path.join(root, "all_hh", "*FAST_estimate.txt"))
    else:
        ptests = ["B2", "B3", "B4", "B5"]
        estfile = glob.glob(os.path.join(root, "all_lw", "*FAST_estimate.txt"))

    testdict = {}
    # get the latest estimation-file
    estfile.sort()
    assert estfile, "No sensitivity results found"
    estfile = estfile[-1]
    print(estfile)
    # load data from file
    testdict["all"] = np.loadtxt(estfile)

    for p in ptests:
        # get the latest estimation-file
        estfile = glob.glob(os.path.join(root, p, "*FAST_estimate.txt"))
        # sort by date
        estfile.sort()
        # take the newest one
        estfile = estfile[-1]
        print(estfile)
        # load data from file
        testdict[p] = np.loadtxt(estfile)

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

    fig = plt.figure(dpi=75, figsize=[10, 4])
    for j, var in enumerate(varnames):

        ax = fig.add_subplot(1, len(varnames), j + 1)  # , sharey=ax1)
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

        plt.xticks([], [])
        if j == 3:
            ax.legend(loc="upper left", handlelength=3, bbox_to_anchor=(1, 1))

    fig.subplots_adjust(wspace=0.5)
    fig.show()
    fig.savefig(
        os.path.join(root, site + "_sensitivity" + file_ext),
        bbox_inches="tight",
    )


plot("HH")
plot("LW")
