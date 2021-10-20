from importedLibraries import *
from classes import *

#Funksjoner
#TODO: create the element matrices. For this we actually have to understand some
#course theory, not just python coding ;)


def getGlobalStiffnessMatrixOfZeros(n):
    '''
    Makes an n by n array of zeros.
    :return: n by n zero matrix
    '''
    SSM = []
    for i in range(n):
        SSM.append(np.zeros(n))
    return np.array(SSM)

def createElementStiffnessMatrixOfZeroes():
    '''
    Makes and returns an element stiffness matrix.
    :return:
    '''
    K = []
    for i in range(6):
        K.append(np.zeros(6))
    K = np.array(K)
    return K

def getMomentFromTriangleLoad(q1,q2,L):
    '''
    Returns a array of the moment distribution along a beam from a triangle pluss rectangle load.
    This can also ofcorse return moment at L/2 as required from the project description.
    :param q1: Initial load
    :param q2: Final load (larger than or equal to q1)
    :param L: length of beam
    :return: a moment distribution list.
    '''

    N=100
    x = np.linspace(0,L,N+1)
    dx = L/N
    momentList = [0]*(N+1)
    q_const = np.abs(min(q1,q2))

    for i in range(N+1):
        #Constant load
        momentList[i] += q_const*(i*dx)*(L-i*dx)/2
        #Triangle load
        if q1 >= q2:
            momentList[i] += (q1 - q2) * (i * dx) / (6 * L) * (2 * L * L - 3 * L * (i * dx) + (i * dx) ** 2)
        else:
            momentList[i] += (q2 - q1) * (L - i * dx) / (6 * L) * (2 * L * L - 3 * L * (L - i * dx) + (L- i * dx) ** 2)

    return np.array(momentList), x

def makeListOfNodeAndBeamClasses(NODE, BEAM):
    '''
    Takes one np array of beams and one of nodes a turns them into a list of node and beam objects.
    :param NODE: np array of all nodes
    :param BEAM: np array of all beams
    :return: list of node and beam classes
    '''

    nodesObjectList = []
    beamsObjectList = []

    for i in range(len(NODE)):
        nodesObjectList.append(Node(NODE[i][0],NODE[i][1],NODE[i][2],NODE[i][3],NODE[i][4]))

    for j in range(len(BEAM)):
        N1 = nodesObjectList[int(BEAM[j][0])-1]
        N2 = nodesObjectList[int(BEAM[j][1])-1]
        beamsObjectList.append(Beam(N1,N2))

    return nodesObjectList, beamsObjectList
