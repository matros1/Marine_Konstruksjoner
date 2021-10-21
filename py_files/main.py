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


    #SystemStivhetsMatrise med bare 0-er
    SSM = getGlobalStiffnessMatrixOfZeros(len(nodesObjectList)*3)

    # TODO - ElementStivhetsMatrise
    ESM = 0

    # TODO -Fastinnspenningsmomenter
    # As all our load is located at nodes, this vector will be 0 (?).
    FIM = 0
    # TODO: createElementStiffnessMatrix()

    #Plots. Here we use the imported library structure visualization to visualize our frame.
    #TODO. Make a visualization of deformations.
    print(nodeArray,beamArray)

    plot(nodeArray, beamArray, materialArray, nodeloadArray, beamloadArray)

    #dette er en forandring

    return 0

main()