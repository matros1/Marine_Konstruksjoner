'''
Developed at the Norwegian University of Science and Technology, Department of Marine Technology
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
        IPEorPipe = beamArray[i][3]
        if IPEorPipe == 1:
            # If beam geometry is IPE
            geo = int(beamArray[i][4] - 1)
            # IPE library includes hight, w top, w bot, w mid, t top, t bot
            beam.makeIPE(IPELibrary[geo][0], IPELibrary[geo][1], IPELibrary[geo][2],
                         IPELibrary[geo][3])

        elif IPEorPipe == 2:
            # Else if beam geometry is pipe
            geo = int(beamArray[i][4] - 1)
            # Pipe library includes radius, ratio air
            beam.makePipe(pipeLibrary[geo][0], pipeLibrary[geo][1])

        i += 1

    return beamsObjectList


def giveYoungsModulusToBeams(beamsObjectList, materialArray, beamArray):
    '''
    This function appends Youngs modulus to elements
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
    Calculates the each beams fixed support moment and forces at both ends.
    :param beamsObjectList: A list of beam objects
    :return: A list of beam objects
    '''
    for beam in beamsObjectList:
        beam.calculateFixedSupport()
    return beamsObjectList


def connectNodeLoadsToNodes(nodesObjectList, nodeloadArray):
    '''
    Conncts the nodeloads to the node objects
    :param nodesObjectList: a list of node objects
    :param nodeloadArray: an array node data from input file
    :return: a list of node objects
    '''

    for i in range(len(nodeloadArray)):
        for j in range(len(nodesObjectList)):
            if (nodesObjectList[j].number == nodeloadArray[i][0]):
                nodesObjectList[j].addNodeLoad(nodeloadArray[i])
    return nodesObjectList


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
    :param beamsObjectList: A list of beam objects
    :param nodesObjectList: A list of node objects
    :return: A 3n x 3n array
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
                M[3 * n + k][3 * n + k] = 1

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


def appendMomentDiagramToBeams(beamsObjectList):
    '''
    Appends moment diagram and related values as member variabels to beam.
    :param beamsObjectList: A list of beam objects
    :return: A list of beam objects
    '''
    for beam in beamsObjectList:
        if beam.hasDistributedLoad:
            beam.appendMomentDiagramToBeam()
    return beamsObjectList


def printBeamSpecsToTerminal(beamsObjectList, s=0, n=999):
    '''
    Prints the reaction forces for beam s to s + n. This is mostly used for debugging.
    :param beamsObjectList: A list of beam objects
    :param s: First beam number
    :param n: First beam + n is last beam to be printed
    :return: prints to console
    '''

    if n > len(beamsObjectList) - s:
        n = len(beamsObjectList) - s
    for i in range(n):
        beamsObjectList[i + s].printBeam()
        beamsObjectList[i + s].printBeamMoments()
        beamsObjectList[i + s].printSecurityFactor()


def calculateMaxMomentAndBendingTension(beamsObjectList):
    '''
    Calculates max moment and bending tension by using beam class member functions.
    :param beamsObjectList: A list of beam objects
    :return: A list of beam objects
    '''
    for beam in beamsObjectList:
        beam.calculateMaxTension()
    return beamsObjectList


def plotMomentDiagram(beamsObjectList):
    Micro = 10 ** -6
    i = 0
    for beam in beamsObjectList:
        if (beam.hasDistributedLoad):
            plt.figure(i)
            plt.plot(np.array(beam.equidistributedX), np.array(beam.momentDiagram) * Micro, c='r', label='M(x)')
            plt.xlabel('x [ m ]')
            plt.ylabel('Moment [ MNm ]')
            plt.title('Bending moment diagram for element ' + str(beam.number))
            plt.grid()
            plt.plot()

            # Extract extreme values and makes the graphs more readable.
            xMax = beam.localMaxMomentX
            maxValueMNm = beam.localMaxMoment * Micro
            startValueMNm = beam.momentDiagram[0] * Micro
            endValueMNm = beam.momentDiagram[-1] * Micro
            plt.plot(xMax, maxValueMNm, 'o', c='black',
                     label='Extream value = ' + str(round(maxValueMNm, 2)) + ' MNm at ' + str(round(xMax, 2)) + ' m')
            plt.plot(beam.equidistributedX[0], startValueMNm, 'o', c='blue',
                     label='Start value = ' + str(round(startValueMNm, 2)) + ' MNm at ' + str(
                         round(beam.equidistributedX[0], 2)) + ' m')
            plt.plot(beam.equidistributedX[-1], endValueMNm, 'o', c='blue',
                     label='End value = ' + str(round(endValueMNm, 2)) + ' MNm at ' + str(
                         round(beam.equidistributedX[-1], 2)) + ' m')
            plt.legend()

            i += 1

    plt.show()


def outputResultsToFile(filename, beamsObjectList, nodesObjectList, r):
    f = open(filename + '.txt', 'w')
    f.write('Node displacements\nNode\tu [mm]\tw [mm]\ttheta [rad]\n')
    for i in range(len(nodesObjectList)):
        f.write(f' {i+1}\t\t{round(r[i*3]*10**3)}\t\t{round(r[i*3+1]*10**3)}\t\t{round(r[i*3+2],4)}\n')
    f.write('\nMoments\nBeam\tM1 [MNm]\tM2 [MNm]\tM_field_max\n')
    for beam in beamsObjectList:
        if(beam.hasDistributedLoad):
            f.write(f' {beam.number}\t\t{round(beam.reactionForces[2]/10**6,2)} \t\t{round(beam.reactionForces[5]/10**6,2)}\t\t{round(beam.localMaxMoment/10**6,2)}\n')
        else:
            f.write(f' {beam.number}\t\t{round(beam.reactionForces[2]/10**6,2)} \t\t{round(beam.reactionForces[5]/10**6,2)}\n')
    f.write('\nShear\nBeam\tQ1 [MN]\tQ2 [MN]\n')
    for beam in beamsObjectList:
        f.write(f' {beam.number}\t\t{round(beam.reactionForces[1]/10**6,2)} \t\t{round(beam.reactionForces[4]/10**6,2)}\n')
    f.write('\nAxial Force\nBeam\tN [MN]\n')
    for beam in beamsObjectList:
        f.write(f' {beam.number}\t\t{round(beam.reactionForces[3]/10**6,2)}\n')
    f.write('\nStress\nBeam\tSigma_x [MPa]\t%fy\n')
    for beam in beamsObjectList:
        f.write(f' {beam.number}\t\t\t{round(beam.sigmax / 10 ** 6,1)}\t\t\t{int(round(beam.securityFactor * 100))}%\n')
    