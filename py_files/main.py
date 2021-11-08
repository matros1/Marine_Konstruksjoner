# Imports sub python files
from functions import *
from plotVisualization import *
from readTextToArray import *
from structure_visualization import *

# This is the command room, here we controll the program.

def main():
    # Reads some files and makes np arrays
    nodeArray, beamArray, materialArray, nodeloadArray, beamloadArray, pipeLibrary, IPELibrary, loadScale= readInputFile("InputDataJacket.txt")
    #nodeArray, beamArray, materialArray, nodeloadArray, beamloadArray, pipeLibrary, IPELibrary, loadScale = readInputFile("InputDataFixedBeam.txt")
    #nodeArray, beamArray, materialArray, nodeloadArray, beamloadArray, pipeLibrary, IPELibrary, loadScale = readInputFile("InputDataPortalFrame.txt")

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

    # Connects the distributed loads to the beam objects,
    # and calculates Fixed Clamping Moment (FastInnspenningsmomenter) for each beam affected by the distributed loads
    beamsObjectList = connectAndScaleDistributedLoadsAndCalculateFixedSupport(beamsObjectList, beamloadArray, loadScale)

    # Connects nodeloads to the nodes
    nodesObjectList = connectNodeLoadsToNodes(nodesObjectList, nodeloadArray)

    # Makes the Resulting Load Vector
    R = makeResultingLoadVector(nodesObjectList)

    # Global stiffness matrix with only zeros
    K = makeGlobalStiffnessMatrix(beamsObjectList, nodesObjectList)

    # Calculate the displacement vector
    r = np.linalg.solve(K,R)

    #Calculate reaction forces
    beamsObjectList = calculateBeamReactionForces(beamsObjectList, r)

    #Calculate stress-levels and 
    beamsObjectList = calculateMaxMomentAndTension(beamsObjectList)

    # Print data to terminal 
    printBeams(beamsObjectList)

    # Plots. Here we use the imported library structure visualization to visualize our frame.
    # The plot only shows rotations, which is scaled by a factor of 20
    plot(nodeArray,beamArray, r * 20)

    # Export data to txt-files
    outputMomentsToFile('moments', beamsObjectList)
    outputSheerToFile('sheer', beamsObjectList)
    outputStressToFile('stress',beamsObjectList)
    outputNormalForceToFile('normalForce', beamsObjectList)

    return 0

main()
