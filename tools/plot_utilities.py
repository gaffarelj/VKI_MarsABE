import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
plt.rcParams.update({'font.size': 13, 'figure.figsize': (10.5, 7), 'savefig.format': 'pdf'})
import sys
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "VKI_MarsABE":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))


def pareto_optimums(points):
    points = np.asarray(points).T
    pareto_optimal = np.ones(points.shape[0], dtype=bool)
    for i, c in enumerate(points):
        if pareto_optimal[i]:
            pareto_optimal[pareto_optimal] = np.any(points[pareto_optimal]<=c, axis=1)
    return pareto_optimal

def comp_pareto(X, Y, front_sign=[1, 1]):
    sl = sorted([[X[i]*front_sign[0], Y[i]*front_sign[1]] for i in range(len(X))])
    pf = [sl[0]]
    for xy in sl[1:]:
        if xy[1] <= pf[-1][1]:
            pf.append(xy)
    x, y = [c[0]*front_sign[0] for c in pf], [c[1]*front_sign[1] for c in pf]
    return x, y

def plot_single(x_data, y_data, x_label, y_label, fname, xlog=False, ylog=False, scatter=False, \
    equal_ax=False, add_front=False, front_sign=[1, 1], z_data=None, z_label="", marker="o", cmap="rainbow", \
    cticks=None, clabels=None, fig=None, ax=None, NB=None, title=None, markersize=None, annot=None, \
    xline=None, yline=None, size=None):
    """
    Simple plot
    """
    # Setup a new figure
    if None in [fig, ax]:
        if size is not None:
            fig, ax = plt.subplots(figsize=size)
        else:
            fig, ax = plt.subplots()
    sys.path = [p for p in sys.path if p != ""]
    # Plot
    if scatter:
        if z_data is not None:
            if markersize is None:
                plt.scatter(x_data, y_data, c=z_data, cmap=cmap, marker=marker)
            else:
                plt.scatter(x_data, y_data, c=z_data, cmap=cmap, marker=marker, s=markersize)
            if cticks is not None:
                if NB is None:
                    cbar = plt.colorbar(label=z_label, ticks=cticks)
                else:
                    sm = plt.cm.ScalarMappable(cmap=cmap, norm=NB[0])
                    cbar = plt.colorbar(sm, label=z_label, ticks=cticks, boundaries=NB[1])
                if clabels is not None:
                    cbar.ax.set_yticklabels(clabels)
            else:
                plt.colorbar(label=z_label)
        else:
            if markersize is None:
                ax.scatter(x_data, y_data, marker=marker)
            else:
                ax.scatter(x_data, y_data, marker=marker, s=markersize)
    else:
        ax.plot(x_data, y_data)
    if annot is not None:
        for i, text in enumerate(annot):
            ax.annotate(text, (x_data[i], y_data[i]), textcoords="offset points", xytext=(0,10), ha='center')
    # Add a pareto front (front_sign is used to specify whether it's best to be high or low
    # If [1, -1]: best is to be high for first objective, low for second)
    if add_front:
        x, y = comp_pareto(x_data, y_data, front_sign)
        ax.step(x, y, where='post', color=(0.35, 0.7, 0.5))
    # Add lines
    if xline is not None:
        plt.axvline(xline, color="black")
    if yline is not None:
        plt.axhline(yline, color="black")
    # Set labels
    ax.set_xlabel(x_label), ax.set_ylabel(y_label)
    # Add title
    if title is not None:
        ax.set_title(title)
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
    # Show the plot
    if fname == "SHOW":
        plt.show()
    # Return the figure and axis
    elif fname == "RETURN":
        return fig, ax
    # Save the plot in the figure folder, as a pdf
    else:
        plt.savefig(sys.path[0]+"/figures/%s.pdf" % fname)
    plt.close()

def plot_multiple(x_datas, y_datas, x_label, y_label, fname, legends=None, colors=None, xlog=False,\
     ylog=False, ylim=None, xlim=None, legend_loc=0, lstyle="solid", fig=None, ax=None, title=None):
    """
    Plot multiple lines
    """
    # Setup a new figure
    if None in [fig, ax]:
        fig, ax = plt.subplots()
    sys.path = [p for p in sys.path if p != ""]
    if type(lstyle) != list:
        lstyle = [lstyle] * len(x_datas)
    # Plot
    for i in range(len(x_datas)):
        # If specified, use the appropriate color
        if colors is not None and len(colors) > i:
            ax.plot(x_datas[i], y_datas[i], color=colors[i], linestyle=lstyle[i])
        else:
            ax.plot(x_datas[i], y_datas[i], linestyle=lstyle[i])
    # Set labels
    ax.set_xlabel(x_label), ax.set_ylabel(y_label)
    # Set legend
    if legends is not None:
        ax.legend(legends, loc=legend_loc)
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
    # Add title
    if title is not None:
        ax.set_title(title)
    # Save space
    fig.tight_layout()
    # Show the plot
    if fname == "SHOW":
        plt.show()
    # Return the figure and axis
    elif fname == "RETURN":
        return fig, ax
    # Save the plot in the figure folder, as a pdf
    else:
        plt.savefig(sys.path[0]+"/figures/%s.pdf" % fname)
    plt.close()


def plot_dual(x_data, y_data_1, y_data_2, x_label, y_label_1, y_label_2, fname, diff_x=False, fig=None, ax1=None, title=None):
    """
    Plot with two y axis. Inputs should be self-explanatory.
    """
    # Setup a new figure
    if None in [fig, ax1]:
        fig, ax1 = plt.subplots()
    # Check if the x_data is the same for both y or not
    if diff_x:
        x_data_1, x_data_2 = x_data[0], x_data[1]
    else:
        x_data_1, x_data_2 = x_data, x_data
    sys.path = [p for p in sys.path if p != ""]
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
    # Add title
    if title is not None:
        ax2.set_title(title)
    # Save space
    fig.tight_layout()
    # Show the plot
    if fname == "SHOW":
        plt.show()
    # Return the figure and axis
    elif fname == "RETURN":
        return fig, ax1, ax2
    # Save the plot in the figure folder, as a pdf
    else:
        plt.savefig(sys.path[0]+"/figures/%s.pdf" % fname)
    plt.close()

def plot_4d(x, y, z, h, labels, fname, fig=None, ax=None):
    # Setup a new figure
    if None in [fig, ax]:
        fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
    sys.path = [p for p in sys.path if p != ""]
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
    # Show the plot
    if fname == "SHOW":
        plt.show()
    # Return the figure and axis
    elif fname == "RETURN":
        return fig, ax
    # Save the plot in the figure folder, as a pdf
    else:
        plt.savefig(sys.path[0]+"/figures/%s.pdf" % fname)
    plt.close()