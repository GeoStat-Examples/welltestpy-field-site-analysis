# -*- coding: utf-8 -*-
"""Estimate HH+LW: single test and combined test estimation."""

import os
from mpi4py import MPI
import welltestpy as wtp
import matplotlib.pyplot as plt

# increase fontsize of plots, prevent type 3 fonts in pdf output
plt.rcParams.update({"font.size": 16, "pdf.fonttype": 42, "ps.fonttype": 42})

# rank is the actual core-number, size is total number of cores
rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size()

# saving estimation results in results
here = os.path.abspath(os.path.dirname(__file__))
results = os.path.normpath(os.path.join(here, "..", "results", "02_estimate"))
data = os.path.normpath(os.path.join(here, "..", "data"))
os.makedirs(results, exist_ok=True)

# setups for single test and combined estimation
setups = [
    [os.path.join(data, "lauswiesen.cmp"), "all_lw"],
    [os.path.join(data, "lauswiesen.cmp"), "B2"],
    [os.path.join(data, "lauswiesen.cmp"), "B3"],
    [os.path.join(data, "lauswiesen.cmp"), "B4"],
    [os.path.join(data, "lauswiesen.cmp"), "B5"],
    [os.path.join(data, "horkheim.cmp"), "all_hh"],
    [os.path.join(data, "horkheim.cmp"), "p05"],
    [os.path.join(data, "horkheim.cmp"), "p40"],
    [os.path.join(data, "horkheim.cmp"), "p42"],
    [os.path.join(data, "horkheim.cmp"), "p44"],
]

for i, setup in enumerate(setups):
    if i % size != rank:
        continue
    campaign_file, well = setup
    print(well + " on core {}/{} started".format(rank + 1, size))
    # load the lauswiesen pumpingtest
    cmp = wtp.load_campaign(campaign_file)
    # get the included tests
    is_all = well.startswith("all")
    testinclude = {well: cmp.tests[well].wells} if not is_all else None
    # set up the estimation
    estimation = wtp.estimate.Theis(
        "est", cmp, testinclude=testinclude, generate=True
    )
    # run the estimation
    estimation.run(
        rep=2500,
        folder=os.path.join(results, well),
        dbname="db",
        traceplotname="paratrace.pdf",
        fittingplotname="fit.pdf",
        interactplotname="parainteract.pdf",
        estname="estimate.txt",
        # run=False,
    )
    # estimate the sensitivites
    estimation.sensitivity(
        folder=os.path.join(results, well),
        dbname="sens_db",
        plotname="sens.pdf",
        traceplotname="sens_trace.pdf",
        sensname="sens_estimate.txt",
    )
