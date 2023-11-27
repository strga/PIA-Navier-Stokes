import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def plot_contour(data, title, filename):
    # Determine the number of data points along one dimension
    N = data.shape[0]

    # Create linearly spaced arrays for x and y coordinates in the domain
    x_dom = np.linspace(0, 1, N)
    y_dom = 1 - np.linspace(0, 1, N)

    # Generate a meshgrid to create a 2D grid from the x and y arrays
    X, Y = np.meshgrid(x_dom, y_dom)

    # Create a new figure for plotting
    plt.figure()

    # Create a filled contour plot with the data
    plt.contourf(X, Y, data, 21)

    # Add a colorbar to the plot for reference
    plt.colorbar()

    # Set the colormap to 'jet'
    plt.set_cmap('jet')

    # Add a title and labels to the plot
    plt.title(title)
    plt.xlabel('x')
    plt.ylabel('y')

    # Save the plot to a file
    plt.savefig(filename)

    # Display the plot
    plt.show()

# Read u velocity data from CSV and convert to numpy array
u_final = pd.read_csv('u_data.csv', header=None).to_numpy()

# Read v velocity data from CSV and convert to numpy array
v_final = pd.read_csv('v_data.csv', header=None).to_numpy()

# Generate and save contour plot for u velocity data
plot_contour(u_final, 'Contour plot of u in the domain', 'u_contour_plot.png')

# Generate and save contour plot for v velocity data
plot_contour(v_final, 'Contour plot of v in the domain', 'v_contour_plot.png')
