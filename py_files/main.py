#Imports sub python files
from structure_visualization import *
from readTextToArray import *
from functions import *
from plotVisualization import *

#This file reads in all subfiles.
#This is the command room, and here we controll the program.
#


def main():
    #Reads some files and makes np arrays
    nodeArray, beamArray, materialArray, nodeloadArray, beamloadArray = readAll()

    #Makes lists of node objects and beam objects from nodes and beams np arrays
    nodesObjectList, beamsObjectList = makeListOfNodeAndBeamClasses(nodeArray, beamArray)

    #Connects the distributed loads to the beam objects,
    #and calculates Fixed Clamping Moment (FastInnspenningsmomenter) for each beam affected by the distributed loads
    connectDistributedNormalLoadsToBeamsAndCalculateFIM(beamsObjectList, beamloadArray)

    #Connects nodeloads to the nodes
    connectNodeLoadsToNodes(nodesObjectList, nodeloadArray)

    #Makes the Resulting Load Vector
    makeResultingLoadVector(nodesObjectList)

    #SystemStivhetsMatrise med bare 0-er
    SSM = getGlobalStiffnessMatrixOfZeros(len(nodesObjectList)*3)

    #Plots. Here we use the imported library structure visualization to visualize our frame.
    plot(nodeArray, beamArray, materialArray, nodeloadArray, beamloadArray)
    #TODO. Make a visualization of deformations.

    return 0

main()