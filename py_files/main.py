# Imports sub python files
from functions import *
from plotVisualization import *
from readTextToArray import *
from structure_visualization import *


# This file reads in all subfiles.
# This is the command room, and here we controll the program.


def main():
    # Reads some files and makes np arrays
    nodeArray, beamArray, materialArray, nodeloadArray, beamloadArray, pipeLibrary, IPELibrary = readInputFile(
        "InputDataJacket.txt")

    # Inizializes beams and nodes objects and appends to list.
    nodesObjectList, beamsObjectList = initializeNodesAndBeamsList(nodeArray, beamArray)

    # Uses spesifications from text files to create correct geometry for each beam.
    beamsObjectList = makeBeamsGeometry(beamArray, beamsObjectList, pipeLibrary, IPELibrary)

    # Gives correct E modulus to each beam.
    beamsObjectList = giveEmodulToBeams(beamsObjectList, materialArray, beamArray)

    # Enumerates beams
    beamsObjectList = giveNumberToObjects(beamsObjectList)

    # Enumerates nodes
    nodesObjectList = giveNumberToObjects(nodesObjectList)

    # Creates the local stiffness matrix of each beam.
    beamsObjectList = giveLocalStiffnessMatrixToBeamsLocalOrientation(beamsObjectList)

    # Orients the local stiffness matrix to global coordinates
    beamsObjectList = giveLocalStiffnessMatrixToBeamsGlobalOrientation(beamsObjectList)

    # Scales the wave forces according to beam diameter
    beamloadArray = dimentionizeLoadsOnBeams(beamArray, beamloadArray, IPELibrary, pipeLibrary)

    # Connects the distributed loads to the beam objects,
    # and calculates Fixed Clamping Moment (FastInnspenningsmomenter) for each beam affected by the distributed loads
    beamsObjectList = connectDistributedNormalLoadsToBeamsAndCalculateFIM(beamsObjectList, beamloadArray, IPELibrary,
                                                        pipeLibrary)

    # Connects nodeloads to the nodes
    nodesObjectList = connectNodeLoadsToNodes(nodesObjectList, nodeloadArray)

    # Makes the Resulting Load Vector
    R = makeResultingLoadVector(nodesObjectList)

    # Global stiffness matrix
    K = makeGlobalStiffnessMatrix(beamsObjectList, nodesObjectList)

    # Calculate the displacement vector
    r = np.linalg.solve(K, R)

    # This works for FixedBeam, and partly works for PortalFrame
    beamsObjectList = calculateBeamReactionForces(beamsObjectList, r)

    printBeam(beamsObjectList)


    # Plots. Here we use the imported library structure visualization to visualize our frame.
    # The plot only shows rotations, which is scaled by a factor of 30
    plot(nodeArray,beamArray, r * 30)

    return False


main()
