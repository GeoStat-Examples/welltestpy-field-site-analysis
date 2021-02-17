# -*- coding: utf-8 -*-
"""Plot HH+LW."""

import os
import welltestpy as wtp
import matplotlib.pyplot as plt

# increase fontsize of plots, prevent type 3 fonts in pdf output
plt.rcParams.update({"font.size": 16, "pdf.fonttype": 42, "ps.fonttype": 42})

# file extension of the saved plots
file_ext = ".pdf"

# paths
here = os.path.abspath(os.path.dirname(__file__))
results = os.path.normpath(os.path.join(here, "..", "results", "00_overview"))
data = os.path.normpath(os.path.join(here, "..", "data"))
os.makedirs(results, exist_ok=True)

# list of campaigns
campaigns = [
    os.path.join(data, "lauswiesen.cmp"),
    os.path.join(data, "horkheim.cmp"),
]

abbrev = {"lauswiesen": "LW", "horkheim": "HH"}

for cmp_file in campaigns:
    cmp = wtp.load_campaign(cmp_file)
    fig1 = cmp.plot(title=False)
    fig2 = cmp.plot_wells().get_figure()  # here the axes are returned
    fig1.savefig(os.path.join(results, abbrev[cmp.name] + "_tests" + file_ext))
    fig2.savefig(os.path.join(results, abbrev[cmp.name] + "_wells" + file_ext))
