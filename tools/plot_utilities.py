import sys
sys.path.insert(0,"\\".join(sys.path[0].split("\\")[:-1]))
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 13, 'figure.figsize': (10.5, 7), 'savefig.format': 'pdf'})

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