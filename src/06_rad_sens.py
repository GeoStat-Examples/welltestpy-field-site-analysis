# -*- coding: utf-8 -*-
"""Radial depending sensitivities."""
import os
import glob
import numpy as np
from scipy.interpolate import UnivariateSpline as uvs
import matplotlib.pyplot as plt
from mpi4py import MPI
import welltestpy as wtp
import anaflow as ana

# plotting style
plt.style.use("ggplot")
plt.rcParams.update({"font.size": 16})

# file extension of the saved plots
file_ext = ".pdf"

# rank is the actual core-number, size is total number of cores
rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size()


# paths
here = os.path.abspath(os.path.dirname(__file__))
results = os.path.normpath(
    os.path.join(here, "..", "results", "01_estimate", "rad_sens")
)


def dashes(i=1, max_d=12, space=1):
    """Dashes for matplotlib."""
    return i * [space, space] + [max_d - 2 * i * space, space]


def lin(xp, fp):
    """
    Linear interpolation of given values as a callable function.

    Parameters
    ----------
    xp : list of float
        The x values of the data points.
    fp : list of float
        The function values of the data points.

    Returns
    -------
    callable
        The linear interpolation.
    """

    def inter(x):
        return np.interp(x, xp, fp)

    return inter


def harm_mu(mu, var):
    """
    Recalculate mu to get the exponent for the harmonic mean in log-norm.

    Parameters
    ----------
    mu : float
        Mean of the log-normal distribution.
    var : TYPE
        Variance of the log-normal distribution.

    Returns
    -------
    float
        Recalculated mean.

    """
    var = 0 if var is None else var
    return mu - 0.5 * var


def est_sens(
    rad,
    mu=-9,
    lnS=-7,
    fix=True,
    dummy=False,
    fix_var=None,
    harm=False,
    var_mean=5,
):
    """
    Estimate radial depending sensitivities.

    Parameters
    ----------
    rad : float
        Radial distance to calculate the sensitivity at.
    mu : float, optional
        Mean of the log-normal distribution. The default is -9.
    lnS : float, optional
        log-storativity. The default is -7.
    fix : bool, optional
        Whether to fix mu and lnS (only var and len_scale estimated).
        The default is True.
    dummy : bool, optional
        Whether to use a dummy paramter. The default is False.
    fix_var : float, optional
        Whether to fix the variance. The default is None.
    harm : bool, optional
        Whether to use the harmonic mean for the reference type-curve.
        Harmonic mean is calculated from fix_var, if given, or var_mean.
        The default is False.
    """
    var_flag = fix_var is not None
    root = results
    root += "_all" if not fix else ""
    root += "_dummy" if dummy else ""
    root += "_fix_var" if var_flag else ""
    root += "_harm" if harm else ""
    # recalculate mu to get the harm-mean if wanted
    harm_var = var_mean if (harm and not var_flag) else fix_var
    mu = harm_mu(mu, harm_var)
    val_fix = {"mu": mu, "lnS": lnS} if fix else {}
    if var_flag:
        val_fix["var"] = fix_var
    rad = float(rad)
    prate = -1
    time = np.geomspace(10, 7200, 10)
    drawdown = ana.theis(time, rad, np.exp(lnS), np.exp(mu), prate)
    campaign = wtp.data.Campaign(name="sens-campaign")
    campaign.add_well(name="well_0", radius=0, coordinates=(0.0, 0.0))
    campaign.add_well(name="well_1", radius=0, coordinates=(rad, 0.0))
    pumptest = wtp.data.PumpingTest("well_0", "well_0", prate)
    pumptest.add_transient_obs("well_1", time, drawdown)
    campaign.addtests(pumptest)
    estimation = wtp.estimate.ExtTheis2D(
        "est", campaign, val_fix=val_fix, generate=True
    )
    if dummy:
        estimation.gen_setup(dummy=dummy)
    estimation.sensitivity(folder=os.path.join(root, "rad_" + str(rad)))


def post_all_sens(
    save=True,
    smooth=None,
    fix=True,
    dummy=False,
    fix_var=None,
    harm=False,
    plt_dummy=False,
    typ="ST",
):
    """
    Post-process sensitivity analysis.

    Parameters
    ----------
    save : bool, optional
        Whether to save the plot. The default is True.
    smooth : bool, optional
        Whether to smooth the result. The default is None.
    fix : bool, optional
        Whether mu and ln(S) were fixed. The default is True.
    dummy : bool, optional
        Whether a dummy paramter was used. The default is False.
    fix_var : float, optional
        The used fixed variance if any. The default is None.
    harm : bool, optional
        Whether the Theis(T_harm) solution was used as reference.
        The default is False.
    plt_dummy : bool, optional
        Whether to plot the dummy result. The default is False.
    typ : str, optional
        The type of the FAST result. Either "ST" for total sensitivity,
        or "S1" for first order sensitivity. The default is "ST".
    """
    var_flag = fix_var is not None
    root = results
    root += "_all" if not fix else ""
    root += "_dummy" if dummy else ""
    root += "_fix_var" if var_flag else ""
    root += "_harm" if harm else ""
    radii_dir = glob.glob(os.path.join(root, "rad*"))
    radii_dir.sort()  # sorting by radii
    radii = []
    rad_Si = []
    sig2_sen = []
    corr_sen = []
    mu_sen = []
    lnS_sen = []
    dum_sen = []
    for rad_dir in radii_dir:
        radii.append(float(os.path.basename(rad_dir)[4:]))
        ST_files = glob.glob(os.path.join(rad_dir, "*_FAST_estimate.txt"))
        ST_files.sort()  # use the latest estimation
        S1_files = glob.glob(os.path.join(rad_dir, "*_FAST_estimate_S1.txt"))
        S1_files.sort()  # use the latest estimation
        ST = np.loadtxt(ST_files[-1])
        S1 = np.loadtxt(S1_files[-1])
        Si = {"S1": S1, "ST": ST}
        rad_Si.append(Si)
        if fix:
            if not var_flag:
                sig2_sen.append(Si[typ][0])
                corr_sen.append(Si[typ][1])
            else:
                sig2_sen.append(0)
                corr_sen.append(Si[typ][0])
            mu_sen.append(0)
            lnS_sen.append(0)
        else:
            if not var_flag:
                mu_sen.append(Si[typ][0])
                sig2_sen.append(Si[typ][1])
                corr_sen.append(Si[typ][2])
                lnS_sen.append(Si[typ][3])
            else:
                mu_sen.append(Si[typ][0])
                sig2_sen.append(0)
                corr_sen.append(Si[typ][1])
                lnS_sen.append(Si[typ][2])
        if dummy:
            dum_sen.append(Si[typ][-1])
        else:
            dum_sen.append(0)

    order = np.argsort(radii)
    radii = np.array(radii)[order]
    sig2_sen = np.array(sig2_sen)[order]
    corr_sen = np.array(corr_sen)[order]
    mu_sen = np.array(mu_sen)[order]
    lnS_sen = np.array(lnS_sen)[order]
    dum_sen = np.array(dum_sen)[order]
    # smoove the lines
    if smooth is not None:
        spmu = uvs(radii, mu_sen, s=smooth)
        spsig = uvs(radii, sig2_sen, s=smooth)
        spcorr = uvs(radii, corr_sen, s=smooth)
        splnS = uvs(radii, lnS_sen, s=smooth)  # s=0.1 0.01
        spdum = uvs(radii, dum_sen, s=smooth)
    else:
        spmu = lin(radii, mu_sen)
        spsig = lin(radii, sig2_sen)
        spcorr = lin(radii, corr_sen)
        splnS = lin(radii, lnS_sen)
        spdum = lin(radii, dum_sen)
    radii2 = np.linspace(radii[0], radii[-1])
    fig, ax0 = plt.subplots(figsize=[10, 4])
    dash_cnt = 1
    if not fix:
        ax0.plot(
            radii2,
            splnS(radii2),
            linewidth=3,
            dashes=dashes(dash_cnt),
            label=r"$\ln(S)$",
        )
        dash_cnt += 1
    if not fix:
        ax0.plot(
            radii2,
            spmu(radii2),
            linewidth=3,
            dashes=dashes(dash_cnt),
            label=r"$\mu$",
        )
        dash_cnt += 1
    if not var_flag:
        ax0.plot(
            radii2,
            spsig(radii2),
            linewidth=3,
            dashes=dashes(dash_cnt),
            label=r"$\sigma^{2}$",
        )
        dash_cnt += 1
    ax0.plot(
        radii2,
        spcorr(radii2),
        linewidth=3,
        dashes=dashes(dash_cnt),
        label=r"$\ell$",
    )
    dash_cnt += 1
    if dummy and plt_dummy:
        ax0.plot(radii2, spdum(radii2), linewidth=3, label="dummy")
    # ax0.legend()
    ax0.legend(loc="upper left", bbox_to_anchor=(1, 1))

    ax0.set_ylim([-0.1, 1.1])
    ax0.set_xlabel(r"Radius of observation $r$ in $[m]$")
    # if typ == "ST":
    #     ax0.set_ylabel(r"FAST total-sensitivity")
    # else:
    #     ax0.set_ylabel(r"FAST first-order-sensitivity")
    fig.tight_layout()
    fig.show()
    if save:
        file_name = os.path.join(root, "Rad_plot_{}".format(typ))
        file_name += "_smooth" if smooth is not None else "_raw"
        file_name += "_dummy" if dummy and plt_dummy else ""
        file_name += file_ext
        fig.savefig(file_name, bbox_inches="tight", dpi=300)


if __name__ == "__main__":
    # size = -1  # hack for serial execution of estimation
    # setup
    fix = False  # whether to fix mu and ln(S) in the estimation
    dummy = False  # whether to use a dummy parameter
    plt_dummy = False  # whether to plot the dummy result
    fix_var = None  # use a fixed value for variance
    harm = False  # use Theis(T_harm) as reference using fix_var or var_mean
    smooth = 1.0  # 0.01 0.1 1.0 None
    if size != 1:  # estimation in parallel (use 04_est_sens.sh)
        radii = np.linspace(0.0, 20.0, 41)
        radii[0] = 0.1
        for i, rad in enumerate(radii):
            if i % size != rank:
                continue
            print("{} on core {}/{} started".format(rad, rank + 1, abs(size)))
            est_sens(rad, fix=fix, dummy=dummy, fix_var=fix_var, harm=harm)
    else:  # plotting in serial (run this file)
        # post_all_sens(
        #     fix=fix,
        #     dummy=dummy,
        #     fix_var=fix_var,
        #     harm=harm,
        #     plt_dummy=plt_dummy,
        # )
        post_all_sens(
            smooth=smooth,
            fix=fix,
            dummy=dummy,
            fix_var=fix_var,
            harm=harm,
            plt_dummy=plt_dummy,
        )
        # post_all_sens(
        #     typ="S1",
        #     fix=fix,
        #     dummy=dummy,
        #     fix_var=fix_var,
        #     harm=harm,
        #     plt_dummy=plt_dummy,
        # )
        # post_all_sens(
        #     typ="S1",
        #     smooth=smooth,
        #     fix=fix,
        #     dummy=dummy,
        #     fix_var=fix_var,
        #     harm=harm,
        #     plt_dummy=plt_dummy,
        # )
