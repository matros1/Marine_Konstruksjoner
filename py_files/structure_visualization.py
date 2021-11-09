# This visualization is based on the original Matlab code by Josef Kiendl, and has been modified to fit TMR4176 by Jon Arnt Kårstad
# NB! Denne filen krever at du har installert Python-pakkene: NumPy, SciPy og Matplotlib
# More detailed information regarding Python (matrices, visualizations etc.) may be found at 'https://www.ntnu.no/wiki/display/imtsoftware'

from importedLibraries import *

def setup_plots():
    '''
    Function that defines initial and deformed axis and figures.
    :return: Figures and axis
    '''
    fig_init, ax_init = plt.subplots()
    fig_def, ax_def = plt.subplots()
    ax_init.set_title('Initialramme')
    ax_def.set_title('Deformert ramme')
    ax_init.axes.set_aspect('equal')
    ax_def.axes.set_aspect('equal')
    return fig_init, ax_init, fig_def, ax_def

def plot_structure(ax, punkt, elem, numbers, index_start):
    '''
    A function that plots a somewhat realistic image of the initial noedes and beams before..
    :param ax:
    :param punkt:
    :param elem:
    :param numbers:
    :param index_start:
    :return: Nothing.
    '''
    # This is a translation of the original function written by Josef Kiendl in Matlab
    # It has been slightly modified in order to be used in TMR4176

    # This function plots the beam structure defined by nodes and elements
    # The bool (0 or 1) 'numbers' decides if node and element numbers are plotted or not

    # Change input to the correct format
    nodes = np.array(punkt[:, 0:2], copy = 1, dtype = int)
    el_nod = np.array(elem[:, 0:2], copy=1, dtype=int) + 1

    # elNodToNodesKonstant is a subtracting factor to go from el_nod[i,0] to a nodes element.
    # This constant is made by Hugo to make this script compatable with our nodes and beams arrays.
    elNodToNodesKonstant = 2

    # Start plotting part
    for iel in range(0, el_nod.shape[0]):
        # Plot element
        ax.plot([nodes[el_nod[iel, 0] - elNodToNodesKonstant, 0], nodes[el_nod[iel, 1] - elNodToNodesKonstant, 0]],
                [nodes[el_nod[iel, 0] - elNodToNodesKonstant, 1], nodes[el_nod[iel, 1] - elNodToNodesKonstant, 1]], '-k', linewidth = 2)

        if numbers == 1:
            # Plot element numbers. These are not plotted in the midpoint to
            # avoid number superposition when elements cross in the middle
            ax.text(nodes[el_nod[iel, 0] - elNodToNodesKonstant, 0] +
                    ( nodes[el_nod[iel, 1] - elNodToNodesKonstant, 0] - nodes[el_nod[iel, 0] - elNodToNodesKonstant, 0] ) / 2.5,
                    nodes[el_nod[iel, 0] - elNodToNodesKonstant, 1] +
                    ( nodes[el_nod[iel, 1] - elNodToNodesKonstant, 1] - nodes[el_nod[iel, 0] - elNodToNodesKonstant, 1] ) / 2.5,
                    str(iel + index_start), color = 'blue', fontsize = 16)

    if numbers == 1:
        # Plot node number
        for inod in range(0, nodes.shape[0]):
            ax.text(nodes[inod, 0], nodes[inod, 1], str(inod + index_start), color = 'red', fontsize = 16)


def plot_structure_def(ax, punkt, elem, numbers, index_start, r):
    '''
    Makes a plot of the deformed beams. Assumes nodes to be fixed.
    :param ax:
    :param punkt:
    :param elem:
    :param numbers:
    :param index_start:
    :param r:
    :return: A plot of the structures nodes and deformed beams.
    '''
    # This is a translation of the original function written by Josef Kiendl in Matlab
    # This function plots the deformed beam structure defined by nodes and elements
    # The bool (0 or 1) 'numbers' decides if node and element numbers are plotted or not

    # Change input to the correct format
    nodes = np.array(punkt[:, 0:2], copy = 1, dtype = int)
    el_nod = np.array(elem[:, 0:2], copy=1, dtype=int)
    nod_dof = np.arange(1, nodes.shape[0] + 1, 1, dtype=int)

    if numbers == 1:
        # Plot node number
        for inod in range(0, nodes.shape[0]):
            ax.text(nodes[inod, 0], nodes[inod, 1], str(inod + index_start), color = 'red', fontsize = 16)

    elNodToNodesKonstant = 1

    for iel in range(0, el_nod.shape[0]):
        delta_x = nodes[el_nod[iel, 1] - elNodToNodesKonstant, 0] - nodes[el_nod[iel, 0] - elNodToNodesKonstant, 0]
        delta_z = nodes[el_nod[iel, 1] - elNodToNodesKonstant, 1] - nodes[el_nod[iel, 0] - elNodToNodesKonstant, 1]
        L = np.sqrt(delta_x ** 2 +  delta_z ** 2)
        if delta_z >= 0:
            psi = np.arccos(delta_x / L)
        else:
            psi = -np.arccos(delta_x / L)

        phi = np.zeros((2, 1))
        for inod in range(0, 2):
            if nod_dof[el_nod[iel, inod] - elNodToNodesKonstant] > 0:
                phi[inod] = r[(el_nod[iel, inod] - elNodToNodesKonstant)*3 + 2]
        x = np.array([0, L])
        z = np.array([0, 0])
        xx = np.arange(0, 1.01, 0.01)*L
        cs = CubicSpline(x, z, bc_type = ((1, -phi[0, 0]), (1, -phi[1, 0])))
        zz = cs(xx)

        # Rotate
        xxzz = np.array([[np.cos(psi), -np.sin(psi)], [np.sin(psi), np.cos(psi)]]) @ np.vstack([xx, zz])

        # Displace
        xx2 = xxzz[0, :] + nodes[el_nod[iel, 0] - elNodToNodesKonstant, 0]
        zz2 = xxzz[1, :] + nodes[el_nod[iel, 0] - elNodToNodesKonstant, 1]
        ax.plot(xx2, zz2, '-k', linewidth = 2)

        if numbers == 1:
            # Plot element numbers. These are not plotted in the midpoint to
            # avoid number superposition when elements cross in the middle
            ax.text(xx2[round(xx2.size / 2.5)], zz2[round(xx2.size / 2.5)], str(iel + index_start), color = 'blue', fontsize = 16)


