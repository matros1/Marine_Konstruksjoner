from classes import *
from importedLibraries import *


# Funksjoner
# TODO: create the element matrices. For this we actually have to understand some
# course theory, not just python coding ;)

# Testing


def initializeNodesAndBeamsList(nodeArray, beamArray):
    nodesObjectList = []
    beamsObjectList = []

    for i in range(len(nodeArray)):
        nodesObjectList.append(Node(nodeArray[i][0], nodeArray[i][1], 0, 0, 0))

    for j in range(len(beamArray)):
        N1 = nodesObjectList[int(beamArray[j][0]) - 1]
        N2 = nodesObjectList[int(beamArray[j][1]) - 1]
        beamsObjectList.append(Beam(N1, N2))

    return nodesObjectList, beamsObjectList


def makeBeamsGeometry(beamList, beamsObjectList, pipeLibrary, IPELibrary):
    i = 0
    TempBeamsList = []
    for beam in beamsObjectList:
        if beamList[i][3] == 1:  # If beam geometry is IPE
            k = beamList[i][4] - 1
            beam.makeIPE(IPELibrary[k][0], IPELibrary[k][1], IPELibrary[k][2],
                         IPELibrary[k][3], IPELibrary[k][4], IPELibrary[k][5])
            # hight, w top, w bot, w mid, t top, t bot
        elif beamList[i][3] == 2:  # If beam geometry is pipe
            k = beamList[i][4] - 1
            beam.makePipe(pipeLibrary[k][0], pipeLibrary[k][1])
            # arguments: radius, ratio air

        i += 1
        TempBeamsList.append(beam)

    return TempBeamsList


def giveEmodulToBeams(beamsObjectList, materialArray, beamData):
    i = 0
    TempBeamsList = []
    for beam in beamsObjectList:
        k = beamData[i][2] - 1
        beam.makeStiffness(materialArray[k][0])
        TempBeamsList.append(beam)
        i += 1

    return TempBeamsList


def giveNumberToObjects(objectList):
    for i in range(len(objectList)):
        objectList[i].giveNumber(i + 1)
    return objectList


def giveLocalStiffnessMatrixToBeamsLocalOrientation(beamsObjectList):
    for i in range(len(beamsObjectList)):
        beamsObjectList[i].makeLocalStiffnessMatrix()
    return beamsObjectList


def giveLocalStiffnessMatrixToBeamsInGlobalCoordinates(beamsObjectList):
    for i in range(len(beamsObjectList)):
        beamsObjectList[i].makeTransformedStiffnessMatrix()
    return beamsObjectList


def connectDistributedNormalLoadsToBeamsAndCalculateFIM(beamsObjectList, beamloadArray):
    '''
    Connects the distributed loads to the beam objects,
        and calculates Fixed Clamping Moment (FastInnspenningsmomenter)
        for each beam affected by the distributed loads
    :param beamsObjectList: list of all Beam-objects
    :param beamloadArray: np array of all distributed loads
    :return: uses member funtion Beam.addDistributedLoad(beamload),
        where beamload is a list with data for the distributed load
        also uses member function Beam.calculateFIM()
    '''
    for i in range(len(beamloadArray)):
        for j in range(len(beamsObjectList)):
            if (beamsObjectList[j].number == beamloadArray[i][0]):
                beamsObjectList[j].addDistributedNormalLoad(beamloadArray[i])
                beamsObjectList[j].calculateFIM()
                # beamsObjectList[j].calculateFISheer(), see classes.py
    return beamsObjectList


def connectNodeLoadsToNodes(nodesObjectList, nodeloadArray):
    '''
    Conncts the nodeloads to the node objects
    param nodesObjectList: list of all Node objects
    param nodeloadArray: np array of all node loads
    '''
    for i in range(len(nodeloadArray)):
        for j in range(len(nodesObjectList)):
            if (nodesObjectList[j].number == nodeloadArray[i][0]):
                nodesObjectList[j].addNodeLoad(nodeloadArray[i])
    return nodesObjectList


def makeListOfNodeAndBeamClasses(nodeArray, beamArray, materialArray, nodeloadArray, beamloadArray, pipeLibrary,
                                 IPELibrary):
    '''
    Takes one np array of beams and one of nodes a turns them into a list of node and beam objects.
    :param NODE: np array of all nodes
    :param BEAM: np array of all beams
    :return: list of node and beam classes
    '''
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
    beamsObjectList = giveLocalStiffnessMatrixToBeamsInGlobalCoordinates(beamsObjectList)

    # Connects the distributed loads to the beam objects,
    # and calculates Fixed Clamping Moment (FastInnspenningsmomenter) for each beam affected by the distributed loads
    connectDistributedNormalLoadsToBeamsAndCalculateFIM(beamsObjectList, beamloadArray)

    # Connects nodeloads to the nodes
    connectNodeLoadsToNodes(nodesObjectList, nodeloadArray)

    return nodesObjectList, beamsObjectList


def getMomentFromTriangleLoad(q1, q2, L):
    '''
    Returns a array of the moment distribution along a beam from a triangle pluss rectangle load.
    This can also ofcorse return moment at L/2 as required from the project description.
    :param q1: Initial load
    :param q2: Final load (larger than or equal to q1)
    :param L: length of beam
    :return: a moment distribution list.
    '''

    N = 100
    x = np.linspace(0, L, N + 1)
    dx = L / N
    momentList = [0] * (N + 1)
    q_const = np.abs(min(q1, q2))

    for i in range(N + 1):
        # Constant load
        momentList[i] += q_const * (i * dx) * (L - i * dx) / 2
        # Triangle load
        if q1 >= q2:
            momentList[i] += (q1 - q2) * (i * dx) / (6 * L) * (2 * L * L - 3 * L * (i * dx) + (i * dx) ** 2)
        else:
            momentList[i] += (q2 - q1) * (L - i * dx) / (6 * L) * (2 * L * L - 3 * L * (L - i * dx) + (L - i * dx) ** 2)

    return np.array(momentList), x


def makeResultingLoadVector(nodesObjectList):
    '''
    Makes list of Fixed Clamping Forces and Vectors
    K*r = R
    K: Global Stiffness Matrix
    r: Deformation vector
    R: Loadvector
    '''
    R = []
    for i in range(len(nodesObjectList)):
        R.append(-nodesObjectList[i].Fx)
        R.append(-nodesObjectList[i].Fz)
        R.append(-nodesObjectList[i].M)
    return np.array(R)


def makeGlobalStiffnessMatrix(beamsObjectList, nodesObjectList):
    n = len(nodesObjectList)
    M = np.zeros((3*n, 3*n))

    for beam in beamsObjectList:
        n1 = beam.node1.number - 1
        n2 = beam.node2.number - 1
        K = beam.transformedStiffnessMatrix
        for i in range(3):
            for j in range(3):
                M[n1 * 3 + i][n1 * 3 + j] += K[i][j]
                M[n2 * 3 + i][n2 * 3 + j] += K[3 + i][3 + j]
                M[n1 * 3 + i][n2 * 3 + j] += K[i][3 + j]
                M[n2 * 3 + i][n1 * 3 + j] += K[3 + i][j]

    print(M)
    return M
