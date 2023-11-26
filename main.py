import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def plot_contour(data, title, filename):
    N = data.shape[0]
    x_dom = np.linspace(0, 1, N)
    y_dom = 1 - np.linspace(0, 1, N)
    X, Y = np.meshgrid(x_dom, y_dom)

    plt.figure()
    plt.contourf(X, Y, data, 21)
    plt.colorbar()
    plt.set_cmap('jet')
    plt.title(title)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig(filename)
    plt.show()

# Read data from CSV
u_final = pd.read_csv('u_data.csv', header=None).to_numpy()
v_final = pd.read_csv('v_data.csv', header=None).to_numpy()

# Generate plots
plot_contour(u_final, 'Contour plot of u in the domain', 'u_contour_plot.png')
plot_contour(v_final, 'Contour plot of v in the domain', 'v_contour_plot.png')
