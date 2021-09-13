import sys
sys.path.insert(0,"\\".join(sys.path[0].split("\\")[:-1]))
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 13, 'figure.figsize': (10.5, 7), 'savefig.format': 'pdf'})

def plot_single(x_data, y_data, x_label, y_label, fname, xlog=False, ylog=False):
    """
    Simple plot
    """
    fig, ax = plt.subplots()
    # Plot
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
    # Save the plot in the figure folder, as a pdf
    plt.savefig("figures/%s.pdf" % fname)

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


def plot_dual(x_data, y_data_1, y_data_2, x_label, y_label_1, y_label_2, fname):
    """
    Plot with two y axis. Inputs should be self-explanatory.
    """
    fig, ax1 = plt.subplots()
    # Set the color of the left y axis
    color = "tab:red"
    ax1.set_xlabel(x_label)
    ax1.set_ylabel(y_label_1, color=color)
    # Plot the data for the left axis
    ax1.plot(x_data, y_data_1, color=color)
    ax1.tick_params(axis="y", labelcolor=color)
    # Add the grid based on the left axis (the horizontal lines are based on the left axis)
    ax1.grid()
    # Create a second y axis that shares the same x axis
    ax2 = ax1.twinx()
    # Set the color of the right y axis
    color = "tab:blue"
    ax2.set_ylabel(y_label_2, color=color)
    # Plot the data on the right axis
    ax2.plot(x_data, y_data_2, color=color)
    ax2.tick_params(axis="y", labelcolor=color)
    # Save space
    fig.tight_layout()
    # Save the plot in the figure folder, as a pdf
    plt.savefig("figures/%s.pdf" % fname)

