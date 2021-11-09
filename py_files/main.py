'''
Developed at the Norwegian University of Science and Technology, Department of Marine Technology
08.11.2021 as part of TMR4167 Marin teknikk - Konstruksjoner
'''

# Imports sub python files
from functions import *
from plotVisualization import *
from readTextToArray import *
from structure_visualization import *


# This is the command room, here we control the program.

def main():
    # Reads from inputfile and sorts the data in appropriate numpy arrays.
    # Make sure the file path is correct for your enviroment's working-directory.
    nodeArray, beamArray, materialArray, nodeloadArray, beamloadArray, pipeLibrary, IPELibrary, referenceDiameter = readInputFile(
        "InputDataJacket.txt")

    # Initialize beams and nodes objects and appends to list.
    nodesObjectList, beamsObjectList = initializeNodesAndBeamsList(nodeArray, beamArray)

    # Uses specifications from input files to create geometry for each beam.
    beamsObjectList = makeBeamsGeometry(beamArray, beamsObjectList, pipeLibrary, IPELibrary)

    # Gives correct Young's modulus to each beam.
    beamsObjectList = giveYoungsModulusToBeams(beamsObjectList, materialArray, beamArray)

    # Enumerates each beam.
    beamsObjectList = giveNumberToObjects(beamsObjectList)

    # Enumerates each node.
    nodesObjectList = giveNumberToObjects(nodesObjectList)

    # Creates the local stiffness matrix for each beam.
    beamsObjectList = giveLocalStiffnessMatrixToBeamsLocalOrientation(beamsObjectList)

    # Orients the local stiffness matrix to global coordinates for each beam.
    beamsObjectList = giveLocalStiffnessMatrixToBeamsGlobalOrientation(beamsObjectList)

    # Connects distributed loads to corresponding beams.
    beamsObjectList = addDistributedLoadsToBeam(beamsObjectList, beamloadArray)

    # Scale distributed loads with respect to reference diameter.
    beamsObjectList = scaleDistributedBeamLoads(beamsObjectList, referenceDiameter)

    # Calculates fixed support moment and forces for each beam.
    beamsObjectList = calculateFixedSupportMomentAndForces(beamsObjectList)

    # Connects node loads to corresponding nodes.
    nodesObjectList = connectNodeLoadsToNodes(nodesObjectList, nodeloadArray)

    # Makes the resulting load vector.
    resultingLoadVector = makeResultingLoadVector(nodesObjectList)

    # Makes the global stiffness matrix.
    globalStiffnessMatrix = makeGlobalStiffnessMatrix(beamsObjectList, nodesObjectList)

    # Uses numpy.linalg.solve() to solve the system relation for the global displacement vector.
    globalDisplacementVector = np.linalg.solve(globalStiffnessMatrix, resultingLoadVector)

    # Calculates reaction forces for each beam.
    beamsObjectList = calculateBeamReactionForces(beamsObjectList, globalDisplacementVector)

    # Calculates max moment and bending tension for each beam.
    beamsObjectList = calculateMaxMomentAndBendingTension(beamsObjectList)

    # Prints specifications for each beam to terminal for debugging.
    printBeamSpecsToTerminal(beamsObjectList)

    # Makes and plots moment diagrams.
    localMaxMomentForDistributedLoadBeams = plotMomentDiagram(beamsObjectList)
    print(localMaxMomentForDistributedLoadBeams)

    # The handed out code structure_visualization.py is altered to match our code and used to plot our jacket
    # construction. Displacements er scaled by 20.
    plot(nodeArray, beamArray, globalDisplacementVector * 20)

    # Export data to txt-files
    outputDataToFile(beamsObjectList)

    return 0


main()
