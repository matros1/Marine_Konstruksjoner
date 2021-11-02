# Imports sub python files
from functions import *
from plotVisualization import *
from readTextToArray import *
from structure_visualization import *


# This file reads in all subfiles.
# This is the command room, and here we controll the program.


def main():
    # Reads some files and makes np arrays
    nodeArray, beamArray, materialArray, nodeloadArray, beamloadArray, pipeLibrary, IPELibrary = readINputFile("InputDataPortalFrame.txt")

    # Makes lists of node objects and beam objects from nodes and beams np arrays
    nodesObjectList, beamsObjectList = makeListOfNodeAndBeamClasses(nodeArray,beamArray, materialArray, nodeloadArray, beamloadArray,pipeLibrary, IPELibrary)

    # Makes the Resulting Load Vector
    R = makeResultingLoadVector(nodesObjectList)

    # Global stiffness matrix with only zeros
    K = makeGlobalStiffnessMatrix(beamsObjectList, nodesObjectList)

    # Calculate the displacement vector
    r = np.linalg.solve(K,R)
    print(r)

    #This does not work yet :(
    #calculateBeamReactionForces(beamsObjectList, r)

    # Plots. Here we use the imported library structure visualization to visualize our frame.
    plot(nodeArray,beamArray, r)

    return 0


main()
