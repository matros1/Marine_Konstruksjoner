from classes import *
from importedLibraries import *


# Funksjoner

def initializeNodesAndBeamsList(nodeArray, beamArray):
    '''
    Makes and initializes the object lists, for both beams and nodes
    param nodeArray: matrix containing node-data from inputfile
    param beamArray: matrix containing beam-data from inputfile
    return: Objectlists for beams and nodes
    '''
    nodesObjectList = []
    beamsObjectList = []

    for i in range(len(nodeArray)):
        nodesObjectList.append(Node(nodeArray[i][0], nodeArray[i][1], nodeArray[i][2], nodeArray[i][3], nodeArray[i][4]))

    for j in range(len(beamArray)):
        N1 = nodesObjectList[int(beamArray[j][0]) - 1]
        N2 = nodesObjectList[int(beamArray[j][1]) - 1]
        beamsObjectList.append(Beam(N1, N2))

    return nodesObjectList, beamsObjectList


def makeBeamsGeometry(beamList, beamsObjectList, pipeLibrary, IPELibrary):
    '''
    Decides from the input file which geometry to add to each beam object
    param beamList: matrix containing beam-data from input file
    param beamsObjectList: list containing all beam-objects
    param pipeLibrary: matrix containing pipe-data from input file
    param IPELibrary: matrix containing IPE-data from input file
    return: new beamsObjectList with added geometry added to the beams
    '''
    i = 0
    TempBeamsList = []
    for beam in beamsObjectList:
        if beamList[i][3] == 1:  # If beam geometry is IPE
            k = int(beamList[i][4] - 1)
            beam.makeIPE(IPELibrary[k][0], IPELibrary[k][1], IPELibrary[k][2],
                         IPELibrary[k][3])
            # hight, w top, w bot, w mid, t top, t bot
        elif beamList[i][3] == 2:  # If beam geometry is pipe
            k = int(beamList[i][4] - 1)
            beam.makePipe(pipeLibrary[k][0], pipeLibrary[k][1])
            # arguments: radius, ratio air

        i += 1
        TempBeamsList.append(beam)

    return TempBeamsList


def giveEmodulToBeams(beamsObjectList, materialArray, beamData):
    i = 0
    TempBeamsList = []
    for beam in beamsObjectList:
        k = int(beamData[i][2] - 1)
        beam.makeStiffness(materialArray[k])
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


def giveLocalStiffnessMatrixToBeamsGlobalOrientation(beamsObjectList):
    for i in range(len(beamsObjectList)):
        beamsObjectList[i].makeTransformedStiffnessMatrix()
    return beamsObjectList


def connectAndScaleDistributedLoadsAndCalculateFixedSupport(beamsObjectList, beamloadArray, base = 1):
    '''
    Connects the distributed loads to the beam objects,
        and calculates Fixed Clamping Moment (FastInnspenningsmomenter)
        for each beam affected by the distributed loads
    :param beamsObjectList: list of all Beam-objects
    :param beamloadArray: np array of all distributed loads
    :return: uses member funtion Beam.addDistributedLoad(beamload),
        where beamload is a list with data for the distributed load
        also uses member function Beam.FixedSupport()
    '''
    for i in range(len(beamloadArray)):
        for j in range(len(beamsObjectList)):
            if (beamsObjectList[j].number == beamloadArray[i][0]):
                beamsObjectList[j].addDistributedNormalLoad(beamloadArray[i])
                beamsObjectList[j].scaleDistributedLoad(base)
                beamsObjectList[j].calculateFixedSupport()
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
    Makes list of Fixed Clamping Forces and Vectors
    K*r = R
    K: Global Stiffness Matrix
    r: Displacement vector
    R: Loadvector
    '''
    R = []
    for i in range(len(nodesObjectList)):
        R.append(-nodesObjectList[i].Fx)    # u
        R.append(-nodesObjectList[i].Fz)    # w
        R.append(-nodesObjectList[i].M)     # ø
    return np.array(R)


def makeGlobalStiffnessMatrix(beamsObjectList, nodesObjectList):
    '''
    Makes the global stiffnessmatrix by adding the transformed stiffnessmatrix for
    each beam into the global, in position given by the nodes each beam is connected to.
    We then account for fixing point coditions (opplagerbetingelser)
    param beamsObjectList: List of all beam-objects
    param nodesObjectList: List of all node-objects
    return: Global stiffness-matrix
    '''
    M = np.zeros((3*len(nodesObjectList), 3*len(nodesObjectList)))

    #Make stiffnessmatrix
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

    #Account for fixing point conditions
    for n,node in enumerate(nodesObjectList):
        for k in range(3):
            if k == 0:
                displacement = node.u
            elif k == 1:
                displacement = node.w
            elif k == 2:
                displacement = node.ø

            if displacement == 1: # Do nothing
                pass
            elif displacement == 2: # Multiply diagonal with 10^6
                M[n * 3 + k][n *3 + k] = M[n * 3 + k][n * 3 + k]*10**6
            elif displacement == 0: # Make diagonal 1, rest of row/collumn 0
                for l in range(len(nodesObjectList)):
                    M[n * 3 + k][l] = 0
                    M[l][n * 3 + k] = 0
                M[n * 3 + k][n * 3 + k] = 1
    #print(M)
    return M

def calculateBeamReactionForces(beamsObjectList, r):
    '''
    NB!! This works for FixedBeam, and partly works for PortalFrame
    Calculates the beams reaction forces and moments
    param beamsObjectList: List of all beam-objects
    param r: vector containing all displacements
    return: adds the reactionforces as member-variables to each beam and returns beamsObjectList
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
            m1 = (1/20)*beam.q1*(beam.length)**2 + (1/30)*beam.q2*(beam.length)**2
            m2 = -(1/30)*beam.q1*(beam.length)**2 - (1/20)*beam.q2*(beam.length)**2

            v2 = (m1 + m2 - (beam.q1*beam.length**2)/6 - (beam.q2*beam.length**2)/3)/beam.length
            v1 = -(beam.length/2)*(beam.q1 + beam.q2) - v2

            beam.reactionForces += np.array([0,v1,m1,0,v2,m2], dtype=float)
        except AttributeError:
            pass
    return beamsObjectList


def printBeams(beamsObjectList, s = 0, n = 999999999):
    '''
    prints the Reactionforces for each beam
    param beamsObjectList: list of all the beam-objects
    param n: if n is given, the function will print forces of the n first beams, else, print for all beams
    return: prints to console
    '''
    if n > len(beamsObjectList)-s:
        n = len(beamsObjectList)-s
    for i in range(n):
        beamsObjectList[i+s].printBeam()
        beamsObjectList[i+s].printBeamMoments()
        beamsObjectList[i+s].printSecurityFactor()

def calculateMaxMomentAndTension(beamsObjectList):
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
