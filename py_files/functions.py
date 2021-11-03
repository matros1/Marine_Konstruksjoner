from classes import *
from importedLibraries import *


# Funksjoner

def initializeNodesAndBeamsList(nodeArray, beamArray):
    nodesObjectList = []
    beamsObjectList = []

    for i in range(len(nodeArray)):
        nodesObjectList.append(
            Node(nodeArray[i][0], nodeArray[i][1], nodeArray[i][2], nodeArray[i][3], nodeArray[i][4]))

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
    return beamsObjectList


def connectPointLoadsToNodesAndCalculateFIM(beamsObjectList, beamPointloadArray, nodesObjectList):
    for beam in beamsObjectList:
        for i in range(len(beamPointloadArray)):
            if beam.number == beamPointloadArray[i][0]:
                P = beamPointloadArray[i][1]
                dL = beamPointloadArray[i][2]
                theta = beam.orientation
                n1 = beam.node1.number
                n2 = beam.node2.number
                P1 = P * (1 - dL)
                P2 = P * dL
                # Fast clamping moment from 2 times statically indetermined beam.
                M1 = -P * dL * (1 - dL) ** 2 / beam.length ** 2
                M2 = P * dL ** 2 * (1 - dL) / beam.length ** 2
                load1 = [0, P1 * np.cos(theta), P1 * np.sin(theta), M1]
                load2 = [0, P2 * np.cos(theta), P2 * np.sin(theta), M2]
                nodesObjectList[n1 - 1].addNodeLoad(load1)
                nodesObjectList[n2 - 1].addNodeLoad(load2)
    return nodesObjectList


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


def makeListOfNodeAndBeamClasses(nodeArray, beamArray, materialArray, nodeloadArray, beamDistributedloadArray,
                                 beamPointloadArray, pipeLibrary, IPELibrary):
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
    connectDistributedNormalLoadsToBeamsAndCalculateFIM(beamsObjectList, beamDistributedloadArray)

    # Connects nodeloads to the nodes
    nodesObjectList = connectNodeLoadsToNodes(nodesObjectList, nodeloadArray)

    # Calculates FIM in node1 and node2 of a beam given by beamPointloadArray and appends these to corresponding nodes.
    nodesObjectList = connectPointLoadsToNodesAndCalculateFIM(beamsObjectList, beamPointloadArray, nodesObjectList)

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
    r: Displacement vector
    R: Loadvector
    '''
    R = []
    for i in range(len(nodesObjectList)):
        R.append(-nodesObjectList[i].Fx)  # u
        R.append(-nodesObjectList[i].Fz)  # w
        R.append(-nodesObjectList[i].M)  # ø
    return np.array(R)


def addNode13Load(Fx, Fz, M, nodesObjectList, beamsObjectList):
    nodesObjectList[8].Fx += Fx / 2
    nodesObjectList[9].Fx += Fx / 2
    nodesObjectList[8].Fz += Fz / 2
    nodesObjectList[9].Fz += Fz / 2
    nodesObjectList[8].M += Fz * beamsObjectList[8].length / 2
    nodesObjectList[9].M -= Fz * beamsObjectList[8].length / 2


def makeGlobalStiffnessMatrix(beamsObjectList, nodesObjectList):
    M = np.zeros((3 * len(nodesObjectList), 3 * len(nodesObjectList)))

    # Make stiffnessmatrix
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

    # Account for fixing point conditions
    for n, node in enumerate(nodesObjectList):
        for k in range(3):
            if k == 0:
                displacement = node.u
            elif k == 1:
                displacement = node.w
            elif k == 2:
                displacement = node.ø

            if displacement == 1:  # Do nothing
                pass
            elif displacement == 2:  # Multiply diagonal with 10^6
                M[n * 3 + k][n * 3 + k] = M[n * 3 + k][n * 3 + k] * 10 ** 6
            elif displacement == 0:  # Make diagonal 1, rest of row/collumn 0
                for l in range(len(nodesObjectList)):
                    M[n * 3 + k][l] = 0
                    M[l][n * 3 + k] = 0
                M[n + k][n + k] = 1

    return M
