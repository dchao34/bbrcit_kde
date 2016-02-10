import numpy as np
import matplotlib.pyplot as plt

def plot_data2d(data2d_fname, ax=None, xlim=None, ylim=None,
                title=None, axis_fontsize=20, tick_labelsize=16):

    # read in the specialized data file
    data = np.genfromtxt(data2d_fname)

    # plot the data
    if ax is None: ax = plt.gca()
    ax.plot(data[:,0],data[:,1],',')

    if xlim:
        ax.set_xlim(xlim)

    if ylim:
        ax.set_ylim(ylim)

    # customize axis labels
    ax.set_xlabel(r'$X_1$', fontsize=axis_fontsize)
    ax.set_ylabel(r'$X_2$', fontsize=axis_fontsize)
    ax.tick_params(axis='both', which='major', labelsize=tick_labelsize)
    if title:
        ax.set_title(title, fontsize=axis_fontsize)

    return


def plot_contour2d(kde2d_fname, ax=None, colorbar=True,
                   vmin=None, vmax=None,
                   nlevels=10, nlevels_f=50,
                   linewidths=0.4,
                   title=None, axis_fontsize=20, tick_labelsize=16):

    # read in the specialized file generated from kde_scan
    X, Y = None, None
    with open(kde2d_fname, 'r') as f:
        x = np.array(map(float, f.next().strip().split()))
        y = np.array(map(float, f.next().strip().split()))
        X, Y = np.meshgrid(x, y)
    Z = np.genfromtxt(kde2d_fname,skip_header=2)

    # plot the data
    if ax is None: ax = plt.gca()
    csf = ax.contourf(X, Y, Z, nlevels_f,
                      vmin=vmin,vmax=vmax,
                      cmap=plt.cm.Blues)
    cs = ax.contour(X, Y, Z, nlevels,linewidths=linewidths, colors='k', linestyles='--')


    # custimize contour labels
    if colorbar:
        cbar = plt.gcf().colorbar(csf)
        cbar.ax.tick_params(labelsize=axis_fontsize)
    ax.clabel(cs, inline=1, fmt='%1.2f', fontsize=axis_fontsize);

    # customize axis labels
    ax.set_xlabel(r'$X_1$', fontsize=axis_fontsize)
    ax.set_ylabel(r'$X_2$', fontsize=axis_fontsize)
    ax.tick_params(axis='both', which='major', labelsize=tick_labelsize)
    if title:
        ax.set_title(title, fontsize=axis_fontsize)

    return
