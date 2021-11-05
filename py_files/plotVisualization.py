from importedLibraries import *
from structure_visualization import *

'''
This subfile controlls plotting. It is dependent on the given code from data studass, "structure_visualization.py" 
'''

def nodeArrayDisplacedPosition(nodeArray, r):
    for i,node in enumerate(nodeArray):
        node[0] += r[3*i]
        node[1] += r[3*i + 1]
    return nodeArray

def plot(NODE, BEAM, r):
    '''
    Makes a plot using given structure visualization plots.
    :param NODE: Array of nodes
    :param BEAM: Array of beams
    :param MATERIAL:
    :param NODELOAD:
    :param BEAMLOAD:
    :return: A plot
    '''
    # Konstant that says whats the nr of the first beam and node is.
    indexStart = 1

    fig_init, ax_init, fig_def, ax_def = setup_plots()

    plot_structure(ax_init, NODE, BEAM, 0, indexStart)

    plot_structure_def(ax_def, nodeArrayDisplacedPosition(NODE, r), BEAM, 0, indexStart, r)

    plt.show()

