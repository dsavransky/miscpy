import numpy as np
import matplotlib
import __main__ as main

if hasattr(main, "__file__"):
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation
import mpl_toolkits.axes_grid1.axes_size as Size
from mpl_toolkits.axes_grid1 import Divider
from mpl_toolkits.basemap import Basemap
import time

from EXOSIMS.util.getExoplanetArchive import getExoplanetArchivePSCP

data = getExoplanetArchivePSCP()

# sort by discovery year
data = data.sort_values(by=["disc_year"]).reset_index(drop=True)

# rename methods with fewer than 30 detections
methods, methods_inds, methods_counts = np.unique(
    data["discoverymethod"].values, return_index=True, return_counts=True
)

othermethods = methods[methods_counts < 50]
for m in othermethods:
    data.loc[data["discoverymethod"] == m, "discoverymethod"] = "Other"

methods = methods[methods_counts >= 50]
methods_counts = methods_counts[methods_counts > 50]
methodorder = np.argsort(methods_counts)[::-1]


ra = data["ra"].values
dec = data["dec"].values
# mthds = np.unique(data['discoverymethod'])
# want the order to be: RV, Transit, Imag, mulens, other
# mthds = np.roll(mthds,2)
mthds = np.hstack([methods[methodorder], "Other"])

# syms = 'o^Dsv';
# cmap = [[1.0000,         0,         0],
#        [0.7500,         0,    0.7500],
#        [0,         0,    1.0000],
#        [0.7500,    0.500,         0],
#        [0.2500,    0.2500,   0.2500],
#        [0,    0.4000,    0],
#        [0,    0.4000,         0.400]]


def setup_fig(syms, cmap, mthds):
    # setup axes
    fignum = 1313
    xsz = 10
    ysz = 6
    ssz = 0.25
    fontsize = 16

    plt.close(fignum)
    fig = plt.figure(fignum, figsize=(xsz, ysz))
    fig.clf()
    rect = (0, 0, 1, 1)  # rect = L,B,W,H
    main_ax = fig.add_axes(rect, label="main")
    lg_ax = fig.add_axes(rect, label="legend")

    # calculate dimentions for all axes
    shm = Size.Fixed(ssz)
    ys = np.array((0.15, 0.85)) * (ysz - ssz * 3)
    horiz = [
        shm,
        Size.Scaled(xsz - ssz * 2),
        shm,
    ]
    vert = [shm, Size.Fixed(ys[0]), shm, Size.Scaled(ys[1]), shm]
    divider = Divider(fig, rect, horiz, vert, aspect=False)

    # distribute axes
    main_ax.set_axes_locator(divider.new_locator(nx=1, ny=3))
    lg_ax.set_axes_locator(divider.new_locator(nx=1, ny=1))

    mp = Basemap(projection="hammer", lat_0=0, lon_0=0, ax=main_ax)
    mp.drawmapboundary(fill_color="white")
    parallels = np.arange(-75, 76, 15)
    mp.drawparallels(parallels, labels=[1, 0, 0, 0], fontsize=fontsize)
    meridians = np.arange(-150, 151, 30)
    mp.drawmeridians(meridians)
    vnudge = 1
    hnudge = 0
    # Run through (lon, lat) pairs, with lat=0 in each pair.
    lons = meridians.copy()
    lats = len(lons) * [0.0]
    for lon, lat in zip(lons, lats):
        x, y = mp(lon + hnudge, lat + vnudge)
        main_ax.text(
            x,
            y,
            u"%0.0f\u00b0" % lon,
            fontsize=fontsize,
            verticalalignment="bottom",
            horizontalalignment="center",
        )

    # set up legend
    lg_ax.get_xaxis().set_visible(False)
    lg_ax.get_yaxis().set_visible(False)
    lg_ax.axis((0, 1, 0, 1))

    linds = ([0, 1], [2, 3], [4])
    lhts = ([0.75, 0.25], [0.75, 0.25], [0.5])
    for j, hts in enumerate(lhts):
        for k, ht in enumerate(hts):
            lg_ax.plot(
                0.1 + 0.8 / 2.5 * j,
                ht,
                syms[j * 2 + k],
                color=cmap[j * 2 + k],
                markersize=12,
            )
            lg_ax.text(
                0.15 + 0.8 / 2.5 * j,
                ht,
                mthds[j * 2 + k],
                fontsize=fontsize,
                horizontalalignment="left",
                verticalalignment="center",
            )

    # and create the year/numstars text
    yrtxt = main_ax.text(
        0,
        0,
        "1900",
        fontsize=fontsize,
        horizontalalignment="left",
        verticalalignment="bottom",
    )
    nstartxt = main_ax.text(
        main_ax.get_xlim()[1],
        0,
        "0",
        fontsize=fontsize,
        horizontalalignment="right",
        verticalalignment="bottom",
    )

    return fig, mp, yrtxt, nstartxt


syms = "os^pvD<"
cmap = ["darkgreen", "red", "silver", "blue", "gray"]

fig, mp, yrtxt, nstartxt = setup_fig(syms, cmap, mthds)

runtime = 45  # sec
fps = 30
nframes = runtime * fps

# generate discovery date array
[years, inds] = np.unique(data["disc_year"].values, return_index=True)
disc_dates = np.zeros(len(dec))
for j in np.arange(0, len(inds) - 1):
    n = inds[j + 1] - inds[j]
    disc_dates[inds[j] : inds[j + 1]] = years[j] + np.linspace(0, 1 - 1.0 / n, n)

# and up to today
currtoy = float(time.strftime("%j")) / 365.25
disc_dates[inds[j + 1] :] = years[j + 1] + np.arange(
    0, currtoy, currtoy / disc_dates[inds[j + 1] :].size
)

#if (int(time.strftime('%Y')) == years[j+1]):
#    #and up to today
#    disc_dates[inds[j+1]:] = years[j+1] + np.arange(0,currtoy,currtoy/disc_dates[inds[j+1]:].size)
#else:
#    n = len(disc_dates) - inds[j+1]
#    disc_dates[inds[j+1]:] = years[j+1] + np.linspace(0,1-1./n,n)


# generate all sizes
hostnames = data["hostname"].values
starnames = np.unique(hostnames)
ssizes = np.zeros(len(starnames)) + 50.0
msizes = np.zeros(len(dec))
for j in range(len(msizes)):
    ii = np.where(starnames == hostnames[j])[0]
    msizes[j] = ssizes[ii]
    # ssizes[ii] *= 2
    ssizes[ii] *= 2 - np.log2(ssizes[ii] / 50.0) * 0.1

xs, ys = mp(ra, dec)
cmap = np.array(cmap)
mthdind = np.array([np.where(mthds == m)[0][0] for m in data["discoverymethod"].values])
syms = np.array(list(syms))

# draw all stars up to frame j
def drawStars(todate, fromdate=0):
    inds = np.where((disc_dates <= todate) & (disc_dates > fromdate))[0]
    ps = ()
    for j in inds:
        p = mp.scatter(
            xs[j],
            ys[j],
            s=msizes[j],
            marker=syms[mthdind[j]],
            linewidths=1,
            edgecolors=cmap[mthdind[j]],
            facecolors="none",
        )
        ps += (p,)
    yrtxt.set_text(u"%0.1f" % todate)
    nstartxt.set_text(u"%0.0f" % len(np.where((disc_dates <= todate))[0]))
    return ps


def drawStarFrames(j):
    print(j)
    if j == 0:
        ps = drawStars(1995, 0)
    else:
        fac = (float(time.strftime("%Y")) + currtoy - 1995) / (nframes - 1)
        ps = drawStars(fac * j + 1995, fac * (j - 1) + 1995)
    return ps


def frameAudio():
    from miscpy.utils.audio import audio
    from scipy.ndimage import gaussian_filter1d
    from EXOSIMS.PlanetPhysicalModel.ForecasterMod import ForecasterMod
    import astropy.units as u

    fm = ForecasterMod()
    a = audio()
    fac = (float(time.strftime("%Y")) + currtoy - 1995) / (nframes - 1)

    """
    #version 1
    freqs = np.array([a.noteFreq('C'),a.noteFreq('E'),a.noteFreq('G'),
                     a.noteFreq('C',octave=1),a.noteFreq('E',octave=1)])

    channels = np.zeros((nframes,len(mthds)))
    for j in range(1,nframes):
        todate = fac*j + 1995
        fromdate = fac*(j-1)+1995
        inds = np.where((disc_dates <= todate) & (disc_dates > fromdate))[0]
        if len(inds) > 0:
            channels[j] = np.histogram(mthdind[inds],range(len(mthds)+1))[0]

    channels = channels.transpose()

    framesamps = np.ones(int(1/fps * a.samplingFrequency))
    out = np.zeros(runtime*a.samplingFrequency)
    for j,c in enumerate(channels):
        o = []
        while len(c) > 0:
            #look for noise
            ind = np.argmax(c > 0)
            if ind == 0:
                ind = len(c)
            o.append(np.zeros(int(np.round(ind/fps * a.samplingFrequency))))
            c = c[ind:]
            if len(c) == 0:
                 continue
            #look for noise
            ind = np.argmax(c == 0)
            if ind == 0:
                ind = len(c)
            tmp = [framesamps*np.sqrt(t) for t in c[:ind]]
            o.append(np.hstack(tmp)*a.tone(freqs[j],duration=ind/fps,e=0.5))
            c = c[ind:]
        out += np.hstack(o)
        out = gaussian_filter1d(out,10)
    a.writeWav(out,'history_animation.wav')
    """

    # version 2
    me = data["pl_bmasse"]
    rade = data.iloc[me.index[me.isna()]]["pl_rade"]
    mfills = fm.calc_mass_from_radius(rade.values * u.Rearth).value
    me.loc[me.isna()] = mfills
    me.loc[me.isna()] = (1 * u.M_jupiter / u.M_earth).decompose().value
    me = me.values
    octaves = 4 - np.digitize(me, np.logspace(np.log10(np.min(me)), 4, 8))

    silence = np.zeros(int(1 / fps * a.samplingFrequency))
    notes = np.array(["C", "E", "G", "B", "D"])
    flats = np.array([False, False, False, True, False])

    out = [silence]
    for j in range(1, nframes):
        todate = fac * j + 1995
        fromdate = fac * (j - 1) + 1995
        inds = np.where((disc_dates <= todate) & (disc_dates > fromdate))[0]
        if len(inds) == 0:
            out.append(silence)
        else:
            x = silence.copy()
            for ii in inds:
                x += a.tone(
                    a.noteFreq(
                        notes[mthdind[ii]], flat=flats[mthdind[ii]], octave=octaves[ii]
                    ),
                    duration=1 / fps,
                )
            out.append(x / np.sqrt(len(inds)))

    out = np.hstack(out)
    out = gaussian_filter1d(out, 10)
    a.writeWav(out, "history_animation2.wav")


if __name__ == "__main__":
    anim = animation.FuncAnimation(
        fig, drawStarFrames, frames=nframes, interval=1.0 / fps * 1000, blit=False
    )
    Writer = animation.writers["ffmpeg"]
    writer = Writer(fps=30)
    anim.save("history_animation.mp4", writer=writer)


"""
to combine audio and video:

ffmpeg -i history_animation.mp4 -i history_animation.wav -c:v copy -map 0:v:0 -map 1:a:0 -c:a aac -b:a 192k output.mp4
"""
