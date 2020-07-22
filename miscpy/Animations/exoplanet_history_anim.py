import requests
import pandas
import numpy as np
import matplotlib
import __main__ as main
if hasattr(main, '__file__'):
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import animation
import mpl_toolkits.axes_grid1.axes_size as Size
from mpl_toolkits.axes_grid1 import Divider
from mpl_toolkits.basemap import Basemap
import time
try:
    from StringIO import StringIO
except ImportError:
    from io import BytesIO as StringIO
import glob

#from astropy.io.votable import parse
#votable = parse('planets.votable')
#table = votable.get_first_table()
#data = table.array


files = glob.glob('data*.pkl')
if files:
    files = np.sort(np.array(files))[-1]
    data = pandas.read_pickle(files)
    print("Loaded data from %s"%files)
else:
    query = """https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=exoplanets&select=*&format=csv"""
    r = requests.get(query)
    data = pandas.read_csv(StringIO(r.content))
    data.to_pickle('data'+time.strftime('%Y%m%d')+'.pkl')


#sort by discovery year
data = data.sort_values(by=['pl_disc']).reset_index(drop=True)

#rename methods with fewer than 30 detections
mthds = np.unique(data['pl_discmethod'])
for m in mthds:
    n = data['pl_discmethod'] == m
    if len(np.where(n)[0]) < 30:
        data.loc[n,'pl_discmethod']  = 'Other'

ra = data['ra'].values
dec = data['dec'].values
mthds = np.unique(data['pl_discmethod'])
#want the order to be: RV, Transit, Imag, mulens, other
mthds = np.roll(mthds,2) 

syms = 'o^Dsv';
cmap = [[1.0000,         0,         0],
        [0.7500,         0,    0.7500],
        [0,         0,    1.0000],
        [0.7500,    0.500,         0],
        [0.2500,    0.2500,   0.2500],        
        [0,    0.4000,    0],
        [0,    0.4000,         0.400]]


def setup_fig(syms,cmap):
    #setup axes
    fignum = 1313
    xsz = 10; ysz = 6; ssz = 0.25
    fontsize = 16

    plt.close(fignum)
    fig = plt.figure(fignum,figsize=(xsz,ysz))
    fig.clf()
    rect = (0,0,1,1) # rect = L,B,W,H
    main_ax = fig.add_axes(rect,label="main")
    lg_ax = fig.add_axes(rect,label="legend")

    #calculate dimentions for all axes
    shm = Size.Fixed(ssz)
    ys = np.array((0.15,0.85))*(ysz - ssz*3)
    horiz = [shm, Size.Scaled(xsz - ssz*2), shm,]
    vert = [shm, Size.Fixed(ys[0]), shm, Size.Scaled(ys[1]), shm]
    divider = Divider(fig, rect, horiz, vert, aspect=False)

    #distribute axes
    main_ax.set_axes_locator(divider.new_locator(nx=1, ny=3))
    lg_ax.set_axes_locator(divider.new_locator(nx=1, ny=1))

    mp = Basemap(projection='hammer',lat_0=0,lon_0=0,ax=main_ax)
    mp.drawmapboundary(fill_color='white')
    parallels = np.arange(-75,76,15)
    mp.drawparallels(parallels, labels=[1,0,0,0],fontsize=fontsize)
    meridians = np.arange(-150,151,30)
    mp.drawmeridians(meridians)
    vnudge=1; hnudge = 0
    # Run through (lon, lat) pairs, with lat=0 in each pair.
    lons = meridians.copy()
    lats = len(lons)*[0.]
    for lon,lat in zip(lons, lats):
        x, y = mp(lon+hnudge, lat+vnudge)
        main_ax.text(x, y, u"%0.0f\u00b0" % lon,
                 fontsize=fontsize,
                 verticalalignment='bottom',
                 horizontalalignment='center')


    #set up legend
    lg_ax.get_xaxis().set_visible(False)
    lg_ax.get_yaxis().set_visible(False)
    lg_ax.axis((0,1,0,1))
            

    linds = ([0,1],[2,3],[4])
    lhts = ([0.75,0.25],[0.75,0.25],[0.5])
    for j,hts in enumerate(lhts):
       for k,ht in enumerate(hts):
           lg_ax.plot(0.1+0.8/2.5*j,ht,syms[j*2+k],color = cmap[j*2+k],markersize=12)
           lg_ax.text(0.15+0.8/2.5*j,ht,mthds[j*2+k],fontsize=fontsize,
                      horizontalalignment='left',verticalalignment='center')

           
    #and create the year/numstars text
    yrtxt = main_ax.text(0,0,'1900',fontsize=fontsize,
                      horizontalalignment='left',verticalalignment='bottom')
    nstartxt = main_ax.text(main_ax.get_xlim()[1],0,'0',fontsize=fontsize,
                      horizontalalignment='right',verticalalignment='bottom')


    return fig,mp,yrtxt,nstartxt

fig,mp,yrtxt,nstartxt = setup_fig(syms,cmap)

runtime = 30 #sec
fps = 30
nframes = runtime*fps

#generate discovery date array
[years,inds] = np.unique(data['pl_disc'].values,return_index=True)
disc_dates = np.zeros(len(dec))
for j in np.arange(0,len(inds)-1):
    n = inds[j+1] - inds[j]
    disc_dates[inds[j]:inds[j+1]] = years[j] + np.linspace(0,1-1./n,n)

#and up to today
currtoy = float(time.strftime('%j'))/365.25
disc_dates[inds[j+1]:] = years[j+1] + np.arange(0,currtoy,currtoy/disc_dates[inds[j+1]:].size)

#generate all sizes
hostnames = data['pl_hostname'].values
starnames = np.unique(hostnames)
ssizes = np.zeros(len(starnames))+50.
msizes = np.zeros(len(dec))
for j in range(len(msizes)):
    ii = np.where(starnames == hostnames[j])[0]
    msizes[j] = ssizes[ii]
    ssizes[ii] *= 2

xs,ys = mp(ra,dec)
cmap = np.array(cmap)
mthdind = np.array([np.where(mthds == m)[0][0] for m in data['pl_discmethod'].values])
syms = np.array(list(syms))

#draw all stars up to frame j
def drawStars(todate,fromdate=0):
    inds = np.where((disc_dates <= todate) & (disc_dates > fromdate))[0]
    ps = ()
    for j in inds:
        p = mp.scatter(xs[j],ys[j], s=msizes[j], marker=syms[mthdind[j]], linewidths=1, edgecolors= cmap[mthdind[j]],facecolors = 'none')
        ps += p,
    yrtxt.set_text(u"%0.1f" %todate)
    nstartxt.set_text(u"%0.0f" %len(np.where((disc_dates <= todate))[0]))
    return ps
    

def drawStarFrames(j):
    print(j)
    if j == 0:
        ps = drawStars(1995,0)
    else:
        fac = (float(time.strftime('%Y'))+currtoy - 1995)/(nframes-1)
        ps = drawStars(fac*j + 1995,fac*(j-1)+1995)
    return ps


if __name__ == '__main__':
    anim = animation.FuncAnimation(fig, drawStarFrames, frames=nframes, interval=1./fps*1000, blit=False)
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=30)
    anim.save('history_animation.mp4', writer=writer)

    
