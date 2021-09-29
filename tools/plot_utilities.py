import sys
sys.path.insert(0,"\\".join(sys.path[0].split("\\")[:-1]))
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import ticker
import numpy as np
from scipy.interpolate import griddata
plt.rcParams.update({'font.size': 13, 'figure.figsize': (10.5, 7), 'savefig.format': 'pdf'})

def plot_single(x_data, y_data, x_label, y_label, fname, xlog=False, ylog=False, scatter=False, equal_ax=False):
    """
    Simple plot
    """
    fig, ax = plt.subplots()
    # Plot
    if scatter:
        ax.scatter(x_data, y_data)
    else:
        ax.plot(x_data, y_data)
    # Set labels
    ax.set_xlabel(x_label), ax.set_ylabel(y_label)
    # Save space
    fig.tight_layout()
    ax.grid()
    # Set log axis if requested
    if xlog:
        plt.xscale("log")
    if ylog:
        plt.yscale("log")
    if equal_ax:
        ax.set_aspect("equal")
    # Save the plot in the figure folder, as a pdf
    plt.savefig("figures/%s.pdf" % fname)
    plt.close()

def plot_multiple(x_datas, y_datas, x_label, y_label, fname, legends=None, colors=None, xlog=False,\
     ylog=False, ylim=None, xlim=None, legend_loc=0, lstyle="solid"):
    """
    Plot multiple lines
    """
    fig, ax = plt.subplots()
    # Plot
    for i in range(len(x_datas)):
        # If specified, use the appropriate color
        if colors is not None and len(colors) > i:
            ax.plot(x_datas[i], y_datas[i], color=colors[i])
        else:
            ax.plot(x_datas[i], y_datas[i], linestyle=lstyle)
    # Set labels
    ax.set_xlabel(x_label), ax.set_ylabel(y_label)
    # Set legend
    if legends is not None:
        ax.legend(legends, loc=legend_loc)
    # Save space
    fig.tight_layout()
    ax.grid()
    # Set log axis if requested
    if xlog:
        plt.xscale("log")
    if ylog:
        plt.yscale("log")
    # Set the axis limits if requested
    if ylim is not None:
        ax.set_ylim(ylim)
    if xlim is not None:
        ax.set_xlim(xlim)
    # Save the plot in the figure folder, as a pdf
    plt.savefig("figures/%s.pdf" % fname)
    plt.close()


def plot_dual(x_data, y_data_1, y_data_2, x_label, y_label_1, y_label_2, fname, diff_x=False):
    """
    Plot with two y axis. Inputs should be self-explanatory.
    """
    # Check if the x_data is the same for both y or not
    if diff_x:
        x_data_1, x_data_2 = x_data[0], x_data[1]
    else:
        x_data_1, x_data_2 = x_data, x_data
    fig, ax1 = plt.subplots()
    # Set the color of the left y axis
    color = "tab:red"
    ax1.set_xlabel(x_label)
    ax1.set_ylabel(y_label_1, color=color)
    # Plot the data for the left axis
    ax1.plot(x_data_1, y_data_1, color=color)
    ax1.tick_params(axis="y", labelcolor=color)
    # Add the grid based on the left axis (the horizontal lines are based on the left axis)
    ax1.grid()
    # Create a second y axis that shares the same x axis
    ax2 = ax1.twinx()
    # Set the color of the right y axis
    color = "tab:blue"
    ax2.set_ylabel(y_label_2, color=color)
    # Plot the data on the right axis
    ax2.plot(x_data_2, y_data_2, color=color)
    ax2.tick_params(axis="y", labelcolor=color)
    # Save space
    fig.tight_layout()
    # Save the plot in the figure folder, as a pdf
    plt.savefig("figures/%s.pdf" % fname)
    plt.close()

def plot_4d(x, y, z, h, labels, fname):
    fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
    grid_x, grid_y = np.mgrid[min(x):max(x):200j, min(y):max(y):200j]
    grid_z = griddata((x, y), z, (grid_x, grid_y), method='nearest')
    grid_c = griddata((x, y), h, (grid_x, grid_y), method='nearest')
    scamap = plt.cm.ScalarMappable(cmap='rainbow')
    fcolors = scamap.to_rgba(grid_c)
    ax.plot_surface(grid_x, grid_y, grid_z, cmap='rainbow', facecolors=fcolors)
    ax.view_init(30, 120)
    clb = fig.colorbar(scamap)
    ax.set_xlabel(labels[0], fontsize=12)
    ax.set_ylabel(labels[1], fontsize=12)
    ax.set_zlabel(labels[2], fontsize=12)
    clb.ax.set_title(labels[3], fontsize=12)
    plt.tight_layout()
    plt.savefig("figures/%s.pdf" % fname)
    plt.close()