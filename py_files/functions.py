'''
Developed by Matias Rosenlund, Christian Lindahl Elseth og Hugo Furnes
Norwegian University of Science and Technology, Department of Marine Technology
08.11.2021 as part of TMR4167 Marin teknikk - Konstruksjoner
'''


from classes import *
from importedLibraries import *


# This file includes relevant functions.

def initializeNodesAndBeamsList(nodeArray, beamArray):
    '''
    Makes and initializes the beams and nodes object lists
    param nodeArray: matrix containing node-data from input file
    param beamArray: matrix containing beam-data from input file
    return: Object lists for nodes, beams
    '''
    nodesObjectList = []
    beamsObjectList = []

    for i in range(len(nodeArray)):
        nodesObjectList.append(
            Node(nodeArray[i][0], nodeArray[i][1], nodeArray[i][2], nodeArray[i][3], nodeArray[i][4]))

    for j in range(len(beamArray)):
        node1 = nodesObjectList[int(beamArray[j][0]) - 1]
        node2 = nodesObjectList[int(beamArray[j][1]) - 1]
        beamsObjectList.append(Beam(node1, node2))

    return nodesObjectList, beamsObjectList


def makeBeamsGeometry(beamArray, beamsObjectList, pipeLibrary, IPELibrary):
    '''
    Decides from the input file which geometry to add to each beam object
    param beamArray: matrix containing beam-data from input file
    param beamsObjectList: list containing all beam-objects
    param pipeLibrary: matrix containing pipe-data from input file
    param IPELibrary: matrix containing IPE-data from input file
    return: Beams object list
    '''
    i = 0
    for beam in beamsObjectList:
        if beamArray[i][3] == 1:  # If beam geometry is IPE
            k = int(beamArray[i][4] - 1)
            beam.makeIPE(IPELibrary[k][0], IPELibrary[k][1], IPELibrary[k][2],
                         IPELibrary[k][3])
            # hight, w top, w bot, w mid, t top, t bot
        elif beamArray[i][3] == 2:  # If beam geometry is pipe
            k = int(beamArray[i][4] - 1)
            beam.makePipe(pipeLibrary[k][0], pipeLibrary[k][1])
            # arguments: radius, ratio air

        i += 1

    return beamsObjectList


def giveYoungsModulusToBeams(beamsObjectList, materialArray, beamArray):
    '''
    This function appends Young's modulus to elements
    :param beamsObjectList: List of beam objects
    :param materialArray: Library of materials
    :param beamArray: Array of beam data from input file
    :return: Beams object list
    '''
    i = 0
    for beam in beamsObjectList:
        materialType = int(beamArray[i][2] - 1)
        # Matrial array includes aluminium and steel
        beam.makeStiffness(materialArray[materialType])
        i += 1

    return beamsObjectList


def giveNumberToObjects(objectList):
    '''
    Enumerates objects in a list
    :param objectList: A list of objects, can be both nodes and beams
    :return: A list of objects
    '''
    for i in range(len(objectList)):
        objectList[i].giveNumber(i + 1)
    return objectList


def giveLocalStiffnessMatrixToBeamsLocalOrientation(beamsObjectList):
    '''
    Uses the Beam member function makeLocalStiffnessMatrix() to make a local stiffness
    matrix to each beam in the list of beam objects
    :param beamsObjectList: A list of beam objects
    :return: A list of beam objects
    '''
    for beam in beamsObjectList:
        beam.makeLocalStiffnessMatrix()
    return beamsObjectList


def giveLocalStiffnessMatrixToBeamsGlobalOrientation(beamsObjectList):
    '''
    Transforms the element stiffness matrix.
    :param beamsObjectList: A list of beam objects
    :return: A list of beam objects
    '''
    for i in range(len(beamsObjectList)):
        beamsObjectList[i].transformElementStiffnessMatrix()
    return beamsObjectList


def addDistributedLoadsToBeam(beamsObjectList, beamloadArray):
    '''
    Adds distributed loads from input file to beam objects
    :param beamsObjectList: A list of beam objects
    :param beamloadArray: An array of beam loads
    :return: A list of beam objects
    '''
    for i in range(len(beamloadArray)):
        for beam in beamsObjectList:
            if (beam.number == beamloadArray[i][0]):
                beam.addDistributedNormalLoad(beamloadArray[i])
    return beamsObjectList


def scaleDistributedBeamLoads(beamsObjectList, referenceDiameter):
    '''
    Scales distributed loads on beam objects according to their diameter with respect
    to the reference diameter.
    :param beamsObjectList: A list of beam objects
    :param referenceDiameter: A reference diameter
    :return: A list of beam objects
    '''
    if referenceDiameter == -1:
        pass
    else:
        for beam in beamsObjectList:
            beam.scaleDistributedLoad(referenceDiameter)
        return beamsObjectList


def calculateFixedSupportMomentAndForces(beamsObjectList):
    '''
    Calculates the each beam's fixed support moment and forces at both ends.
    :param beamsObjectList: A list of beam objects
    :return: A list of beam objects
    '''
    for beam in beamsObjectList:
        beam.calculateFixedSupport()
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
    Makes a list of fixed support forces and moments that satisfy
    the system relation, K*r = R.
    :param nodesObjectList: A list of nodes
    :return: An 1x3n array
    '''
    R = []
    for i in range(len(nodesObjectList)):
        R.append(-nodesObjectList[i].Fx)  # u
        R.append(-nodesObjectList[i].Fz)  # w
        R.append(-nodesObjectList[i].M)  # ø
    return np.array(R)


def makeGlobalStiffnessMatrix(beamsObjectList, nodesObjectList):
    '''
    Makes the global stiffness matrix by adding each beams transformed stiffness matrix in position
    given by the nodes each beam is connected to.
    We then account for fixed support conditions.
    param beamsObjectList: List of beam objects
    param nodesObjectList: List of node objects
    return: Global stiffness matrix, 3n x 3n
    '''
    M = np.zeros((3 * len(nodesObjectList), 3 * len(nodesObjectList)))

    # Make stiffness matrix
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

    # Account for fixed support conditions
    for n, node in enumerate(nodesObjectList):
        for k in range(3):
            if k == 0:
                displacement = node.supportX
            elif k == 1:
                displacement = node.supportZ
            elif k == 2:
                displacement = node.supportTheta

            if displacement == 1:  # Do nothing
                pass
            elif displacement == 2:  # Multiply diagonal with 10^6
                M[n * 3 + k][n * 3 + k] = M[n * 3 + k][n * 3 + k] * 10 ** 6
            elif displacement == 0:  # Make diagonal 1, rest of row/collumn 0
                for l in range(len(nodesObjectList)):
                    M[n * 3 + k][l] = 0
                    M[l][n * 3 + k] = 0
                M[n * 3 + k][n * 3 + k] = 1

    return M


def calculateBeamReactionForces(beamsObjectList, r):
    '''
    Calculates beams reaction force and adds these as member-variables to each beam
    param beamsObjectList: List of beam objects
    param r: Array containing resulting loads
    return: A list of beam objects
    '''
    for beam in beamsObjectList:
        localDisplacements = np.zeros(6)
        n1 = beam.node1.number - 1
        n2 = beam.node2.number - 1
        for i in range(3):
            localDisplacements[i] = r[3 * n1 + i]
            localDisplacements[3 + i] = r[3 * n2 + i]
        beam.localDisplacements = np.matmul(beam.T_transponent, localDisplacements)
        beam.reactionForces = np.matmul(beam.localStiffnessMatrix, beam.localDisplacements)
        try:
            m1 = (1 / 20) * beam.q1 * (beam.length) ** 2 + (1 / 30) * beam.q2 * (beam.length) ** 2
            m2 = -(1 / 30) * beam.q1 * (beam.length) ** 2 - (1 / 20) * beam.q2 * (beam.length) ** 2

            v2 = (m1 + m2 - (beam.q1 * beam.length ** 2) / 6 - (beam.q2 * beam.length ** 2) / 3) / beam.length
            v1 = -(beam.length / 2) * (beam.q1 + beam.q2) - v2

            beam.reactionForces += np.array([0, v1, m1, 0, v2, m2], dtype=float)
        except AttributeError:
            pass
    return beamsObjectList


def printBeamSpecsToTerminal(beamsObjectList, s=0, n=999999999):
    '''
    Prints the reaction forces for each beam. This is mostly used for debugging.
    :param beamsObjectList:
    :param s: A list of beam objects
    :param n: if n is given, the function will print forces of the n first beams, else, print for all beams
    :return: prints to console
    '''

    if n > len(beamsObjectList) - s:
        n = len(beamsObjectList) - s
    for i in range(n):
        beamsObjectList[i + s].printBeam()
        beamsObjectList[i + s].printBeamMoments()
        beamsObjectList[i + s].printSecurityFactor()


def calculateMaxMomentAndTension(beamsObjectList):
    '''
    Calculates max moment and bending tension.
    :param beamsObjectList: A list of beam objects
    :return: Nothing
    '''
    for i in range(len(beamsObjectList)):
        beamsObjectList[i].calculateMaxMoment()
        beamsObjectList[i].calculateMaxBendingTension()
    return beamsObjectList

def outputMomentsToFile(filename, beamsObjectList):
    f = open(filename+'.txt', 'w')
    for beam in beamsObjectList:
        f.write(f'Beam: {beam.number},\t M1: {round(beam.reactionForces[2])}[Nm], \t M2: {round(beam.reactionForces[5])}[Nm]\n')
    f.close

def outputSheerToFile(filename, beamsObjectList):
    f = open(filename+'.txt', 'w')
    for beam in beamsObjectList:
        f.write(f'Beam: {beam.number},\t V1: {round(beam.reactionForces[1])}[N], \t V2: {round(beam.reactionForces[4])}[N]\n')
    f.close

def outputStressToFile(filename, beamsObjectList):
    f = open(filename+'.txt', 'w')
    for beam in beamsObjectList:
        f.write(f'Beam: {beam.number},\t sigma_x: {round(beam.sigmax/10**6,2)}[MPa]\n')
    f.close

def outputNormalForceToFile(filename, beamsObjectList):
    f = open(filename+'.txt', 'w')
    for beam in beamsObjectList:
        f.write(f'Beam: {beam.number},\t N: {round(beam.reactionForces[0])}[N]\n')
    f.close

def outputDataToFile(beamsObjectList):
    outputMomentsToFile('moments', beamsObjectList)
    outputSheerToFile('sheer', beamsObjectList)
    outputStressToFile('stress',beamsObjectList)
    outputNormalForceToFile('normalForce', beamsObjectList)


#Saved for later


def getMomentDiagram(L, M1, V1, q1, q2, N):
    X = np.linspace(0, L, N + 1)
    dx = L / N
    momentList = []

    # for i in range(N + 1):
    #     m1 = q1 * (i * dx) / (6 * L) * (2 * L * L - 3 * L * (i * dx) + (i * dx) ** 2)
    #     m2 = q2 * (L - i * dx) / (6 * L) * (2 * L * L - 3 * L * (L - i * dx) + (L - i * dx) ** 2)
    #     momentList.append(-M1 + m1 + m2)

    # for j in range(N + 1):
    #     x = j * dx
    #     qs = q1 * (1 - x / L) + q1 * x / L
    #
    #     if (q2 > q1):
    #         ms = - M1 - V1*x + (-q1 * x ** 2) / 2 + (-q2 + q1) * (x ** 2) * 2/ (3 * L)
    #     else:
    #         ms = - M1 - V1*x + (-q2 * x ** 2) / 2 + (-q1+q2)*(1-x/L)*(x ** 2)/3 - (q1-q2)*x**2/6
    #     momentList.append(ms)

    for j in range(N+1):
        x = j * dx
        qs = q1 * (1 - x / L) + q1 * x / L

        if (q2 > q1):
            ms = - M1 - V1 * x + (-q1 * x ** 2) / 2 + (-qs + q1) * (x ** 2) / (6)
        else:
            ms = - M1 - V1 * x + (-qs * x ** 2) / 2 + (-q1 + qs) * (x ** 2) / 3
        momentList.append(ms)

    return np.array(X), np.array(momentList)


def printMomentDiagram(beamsObjectList):
    '''
    Prints a moment diagram to screen for each beam exposed to a distributed load.
    :param beamsObjectList: A list of beam objects
    :return: nothing
    '''
    N = 100
    Micro = 10 ** -6

    plotList = []
    tempBeamList = []
    for beam in beamsObjectList:
        try:

            plotList.append(getMomentDiagram(beam.length, beam.reactionForces[2], beam.reactionForces[1], beam.q1, beam.q2, N))
            tempBeamList.append(beam)
            print(beam.q1, beam.q2, beam.reactionForces[2], beam.reactionForces[5], beam.reactionForces[1])
        except AttributeError:
            pass

    for i in range(len(plotList)):
        plt.figure(i)
        plt.plot(plotList[i][0], plotList[i][1] * Micro, c='r', label='M(x)')
        plt.xlabel('x [ m ]')
        plt.ylabel('Moment [ MNm ]')
        plt.title('Bending moment diagram for element ' + str(tempBeamList[i].number))
        plt.grid()

        plt.plot()
        if (plotList[i][1][0] < 0):
            maxValue = np.max(plotList[i][1])
            iMax = np.argmax(plotList[i][1])
        else:
            maxValue = np.min(plotList[i][1])
            iMax = np.argmin(plotList[i][1])
        xMax = iMax * tempBeamList[i].length / N
        maxValueMNm = maxValue * Micro
        startValueMNm = plotList[i][1][0] * Micro
        endValueMNm = plotList[i][1][-1] * Micro
        plt.plot(xMax, maxValueMNm, 'o', c='black',
                 label='Extream value = ' + str(round(maxValueMNm, 2)) + ' MNm at ' + str(round(xMax, 2)) + ' m')
        plt.plot(plotList[i][0][0], startValueMNm, 'o', c='blue',
                 label='Start value = ' + str(round(startValueMNm, 2)) + ' MNm at ' + str(
                     round(plotList[i][0][0], 2)) + ' m')
        plt.plot(plotList[i][0][-1], endValueMNm, 'o', c='blue',
                 label='End value = ' + str(round(endValueMNm, 2)) + ' MNm at ' + str(
                     round(plotList[i][0][-1], 2)) + ' m')
        plt.legend()

    plt.show()
