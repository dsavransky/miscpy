# %pylab --no-import-all
import pandas
import os

import astropy.units as u
import astropy.constants as const
import re

import matplotlib
import numpy as np
import matplotlib.pyplot as plt

matplotlib.rcParams.update({"font.size": 16})

from EXOSIMS.PlanetPhysicalModel.ForecasterMod import ForecasterMod

fmod = ForecasterMod()
RfromM = fmod.calc_radius_from_mass
MfromR = fmod.calc_mass_from_radius

# get data
from EXOSIMS.util.getExoplanetArchive import getExoplanetArchivePSCP

data = getExoplanetArchivePSCP()

# this turncates to only those methods with > 30 detections (skip if you want all of them)
methods, methods_inds, methods_counts = np.unique(
    data["discoverymethod"].values, return_index=True, return_counts=True
)
methods = methods[methods_counts > 30]
methods_counts = methods_counts[methods_counts > 30]
methodorder = np.argsort(methods_counts)[::-1]

inds = data["discoverymethod"] == methods[0]
for method in methods:
    inds = inds | (data["discoverymethod"] == method)
data = data.loc[inds].reset_index(drop=True)

# imaging discoveries mean some photometry available
data.loc[data["discoverymethod"] == "Imaging", "pl_nespec"] = 1

#### drop everything after 2010
data2 = data.loc[data["disc_year"] < 2011].reset_index(drop=True)


#### spectra refs
datestr = "2021.04.21"
basepath = os.path.join(os.environ["HOME"], "Downloads")

tspec = pandas.read_csv(
    os.path.join(basepath, "transitspec_" + datestr + os.extsep + "csv"), comment="#"
)
espec = pandas.read_csv(
    os.path.join(basepath, "emissionspec_" + datestr + os.extsep + "csv"), comment="#"
)

pnames = np.hstack([tspec["plntname"].values, espec["plntname"].values])
pubs = np.hstack([tspec["plntranreflink"].values, espec["plntreflink"].values])

pubdates = np.array([re.search("(20\d\d)", feh).groups(0)[0] for feh in pubs])
pubdates = pubdates.astype(float)

upnames = np.unique(pnames)
firstpubdate = np.zeros(upnames.shape)
for j, n in enumerate(upnames):
    firstpubdate[j] = pubdates[pnames == n].min()

nopre2010dat = upnames[firstpubdate > 2010]

for n in nopre2010dat:
    data2.loc[data2["pl_name"] == n, "pl_nespec"] = 0
    data2.loc[data2["pl_name"] == n, "pl_ntranspec"] = 0

######### Plot setup
syms = "os^pvD<"
cmap = [
    [0, 0, 1.0000],
    [1.0000, 0, 0],
    [0, 0.4000, 0],
    [0, 0.4000, 0.4000],
    [0.7500, 0.500, 0],
    [0.2500, 0.2500, 0.2500],
    [0.7500, 0, 0.7500],
]

cmap = ["blue", "silver", "orange", "darkmagenta"]
# solar system planets
GMsun = 1.32712440018e20  # %m^3/s^2
G = 6.67428e-11  # m^3/kg/s^2
Msun = GMsun / G
MEarth = Msun / 328900.56 / (1 + 1 / 81.30059)
Ms = np.hstack(
    (
        Msun / np.array([6023600.0, 408523.71]),
        MEarth,
        Msun / np.array([3098708.0, 1047.3486, 3497.898, 22902.98, 19412.24, 1.35e8]),
    )
)
Ms = Ms / Ms[4]
Mse = Ms / Ms[2]
Rs = np.array(
    [2440.0, 6052.0, 6378.0, 3397.0, 71492.0, 60268.0, 25559.0, 24766.0, 24766.0]
)
Rs = Rs / Rs[2]
smas = np.array([0.3871, 0.7233, 1, 1.524, 5.203, 9.539, 19.19, 30.06, 39.48])
planetnames = [
    "Mercury",
    "Venus",
    "Earth",
    "Mars",
    "Jupiter",
    "Saturn",
    "Uranus",
    "Neptune",
    "Pluto",
]
has = ["left", "left", "left", "left", "left", "left", "left", "left", "left"]
offs = [
    (4, -3),
    (-5, -16),
    (4, -8),
    (4, -4),
    (6, -4),
    (5, -2),
    (0, -15),
    (-5, 4),
    (0, 0),
]

p2sma = lambda mu, T: ((mu * T**2 / (4 * np.pi**2)) ** (1 / 3.0)).to("AU")


############# First plot (just plot everything)
def plot1(data, earthmassonly=False, title=None):
    fig, ax = plt.subplots(figsize=(10, 5))
    fig.subplots_adjust(bottom=0.15, top=0.95, left=0.125, right=0.9)
    for m, s, c, o in zip(methods, syms, cmap, methodorder):
        inds = data["discoverymethod"] == m
        mj = data.loc[inds, "pl_bmassj"]
        sma = data.loc[inds, "pl_orbsmax"]

        # fill in missing masses from radii
        radj = data.iloc[mj.index[mj.isna()]]["pl_radj"]
        mfills = MfromR(radj.values * u.R_jup).to(u.M_jup).value
        mj.loc[mj.isna()] = mfills

        # fill in missing smas from period
        orbper = data.iloc[sma.index[sma.isna()]]["pl_orbper"]
        stmass = data.iloc[sma.index[sma.isna()]]["st_mass"]

        GMs = const.G * (stmass.values * u.solMass)  # units of solar mass
        T = orbper.values * u.day
        smafill = p2sma(GMs, T)
        sma.loc[sma.isna()] = smafill.value

        if earthmassonly:
            mj /= Ms[2]
        ax.scatter(
            sma,
            mj,
            marker=s,
            s=60,
            zorder=o,
            facecolors=c,
            edgecolors="k",
            alpha=0.75,
            label=m,
        )

    if earthmassonly:
        tmp = Mse
    else:
        tmp = Ms
    ax.scatter(
        smas,
        tmp,
        marker="o",
        s=60,
        facecolors="yellow",
        edgecolors="k",
        alpha=1,
        zorder=methodorder.max(),
    )
    for a, m, n, ha, off in zip(smas, tmp, planetnames, has, offs):
        ax.annotate(n, (a, m), ha=ha, xytext=off, textcoords="offset points")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim([1e-2, 1e3])
    ax.set_xlabel("Semi-Major Axis (AU)")

    if earthmassonly:
        ax.set_ylim(np.array([1e-4, 40]) / Ms[2])
        ax.set_ylabel("(Minimum) Mass (Earth Masses)")
    else:
        ax.set_ylim([1e-4, 40])
        ax.set_ylabel("(Minimum) Mass (M$_J$)")
    ax.legend(loc="lower right", scatterpoints=1, fancybox=True, prop={"size": 14})

    if not earthmassonly:
        ax2 = ax.twinx()
        ax2.set_yscale("log")
        ax2.set_ylim(np.array(ax.get_ylim()) / Ms[2])
        ax2.set_ylabel("M$_\oplus$")
        plt.subplots_adjust(right=0.88)
    else:
        plt.subplots_adjust(right=0.95)

    if title:
        plt.title(title)
        plt.subplots_adjust(top=0.93)


plot1(data, earthmassonly=True, title="2020: {} planets".format(len(data)))
plot1(data2, earthmassonly=True, title="2010: {} planets".format(len(data2)))
############# Second plot (highlingt things with spectral info)


def plot2(data, earthmassonly=False, title=None):
    fig, ax = plt.subplots(figsize=(10, 5))
    fig.subplots_adjust(bottom=0.15, top=0.95, left=0.125, right=0.9)
    for m, s, c, o in zip(methods, syms, cmap, methodorder):
        inds = data["discoverymethod"] == m
        mj = data.loc[inds, "pl_bmassj"]
        sma = data.loc[inds, "pl_orbsmax"]
        nspecs = data.loc[inds, "pl_nespec"] + data.loc[inds, "pl_ntranspec"]

        # fill in missing masses from radii
        radj = data.iloc[mj.index[mj.isna()]]["pl_radj"]
        mfills = MfromR(radj.values * u.R_jup).to(u.M_jup).value
        mj.loc[mj.isna()] = mfills

        # fill in missing smas from period
        orbper = data.iloc[sma.index[sma.isna()]]["pl_orbper"]
        stmass = data.iloc[sma.index[sma.isna()]]["st_mass"]

        GMs = const.G * (stmass.values * u.solMass)  # units of solar mass
        T = orbper.values * u.day
        smafill = p2sma(GMs, T)
        sma.loc[sma.isna()] = smafill.value

        if earthmassonly:
            mj /= Ms[2]

        if m == "Imaging":
            ax.scatter(
                sma,
                mj,
                marker=s,
                s=60,
                zorder=o,
                facecolors=c,
                edgecolors="k",
                alpha=0.75,
                label=m,
            )
        else:
            ax.scatter(
                sma.loc[nspecs > 0],
                mj.loc[nspecs > 0],
                marker=s,
                s=60,
                zorder=o,
                facecolors=c,
                edgecolors="k",
                alpha=0.95,
                label=m,
            )
            ax.scatter(
                sma.loc[nspecs == 0],
                mj.loc[nspecs == 0],
                marker=s,
                s=60,
                zorder=-1,
                facecolors=c,
                edgecolors="k",
                alpha=0.1,
                label=None,
            )

    if earthmassonly:
        tmp = Mse
    else:
        tmp = Ms

    ax.scatter(
        smas,
        tmp,
        marker="o",
        s=60,
        facecolors="yellow",
        edgecolors="k",
        alpha=1,
        zorder=methodorder.max(),
    )

    for a, m, n, ha, off in zip(smas, tmp, planetnames, has, offs):
        ax.annotate(n, (a, m), ha=ha, xytext=off, textcoords="offset points")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim([1e-2, 1e3])
    ax.set_xlabel("Semi-Major Axis (AU)")

    if earthmassonly:
        ax.set_ylim(np.array([1e-4, 40]) / Ms[2])
        ax.set_ylabel("(Minimum) Mass (Earth Masses)")
    else:
        ax.set_ylim([1e-4, 40])
        ax.set_ylabel("(Minimum) Mass (M$_J$)")

    ax.legend(loc="lower right", scatterpoints=1, fancybox=True, prop={"size": 14})

    if not earthmassonly:
        ax2 = ax.twinx()
        ax2.set_yscale("log")
        ax2.set_ylim(np.array(ax.get_ylim()) / Ms[2])
        ax2.set_ylabel("M$_\oplus$")
        plt.subplots_adjust(right=0.88)
    else:
        plt.subplots_adjust(right=0.95)

    if title:
        plt.title(title)
        plt.subplots_adjust(top=0.93)


plot2(data, earthmassonly=True, title="2020")
plot2(
    data,
    earthmassonly=True,
    title="2020: {} planets with spectroscopic measurements".format(
        np.where(data["pl_nespec"] + data["pl_ntranspec"] > 0)[0].size
    ),
)
plot2(
    data2,
    earthmassonly=True,
    title="2010: {} planets with spectroscopic measurements".format(
        np.where(data2["pl_nespec"] + data2["pl_ntranspec"] > 0)[0].size
    ),
)
