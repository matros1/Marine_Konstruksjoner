# Imports sub python files
from functions import *
from plotVisualization import *
from readTextToArray import *
from structure_visualization import *


# This file reads in all subfiles.
# This is the command room, and here we controll the program.


def main():
    # Reads some files and makes np arrays
    # Reads some files and makes np arrays
    nodeArray, beamArray, materialArray, nodeloadArray, beamdistributedloadArray, beamPointLoadArray, pipeLibrary, IPELibrary = readAll()


    # Makes lists of node objects and beam objects from nodes and beams np arrays
    nodesObjectList, beamsObjectList = makeListOfNodeAndBeamClasses(nodeArray, beamArray, materialArray, nodeloadArray,
                                                                    beamdistributedloadArray, beamPointLoadArray,
                                                                    pipeLibrary, IPELibrary)



    # Makes the Resulting Load Vector
    R = makeResultingLoadVector(nodesObjectList)

    # Global stiffness matrix
    K = makeGlobalStiffnessMatrix(beamsObjectList, nodesObjectList)

    # Calculate the displacement vector
    r = np.linalg.solve(K, R)
    print(r[8*3+2], r[9*3])

    # TODO: All node are assumed to be fixed, which is not physical. Especially node 11.

    # Plots. Here we use the imported library structure visualization to visualize our frame.
    # the code runs, but deformations look too small or non existent
    plot(nodeArray, beamArray, r)

    return 0


main()
