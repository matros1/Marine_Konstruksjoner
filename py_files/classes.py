'''
Developed by Matias Rosenlund, Christian Lindahl Elseth og Hugo Furnes
Norwegian University of Science and Technology, Department of Marine Technology
08.11.2021 as part of TMR4167 Marin teknikk - Konstruksjoner
'''

from functions import *
from importedLibraries import *


# This file handles the class decleration for Nodes and Beams

class Beam:
    '''
    This is a beam class. It includes multiple in house functions. Some is called when initializing and some
    can be called if you want to add info to the beam, such as geometry or number.
    '''

    def __init__(self, NODE1, NODE2):
        '''
        Initializes a Beam object. Calculates orientation and length in global coordinates.
        :param NODE1: start node
        :param NODE2: end node
        '''
        self.node1 = NODE1
        self.x1 = NODE1.x
        self.z1 = NODE1.z
        self.node2 = NODE2
        self.x2 = NODE2.x
        self.z2 = NODE2.z
        self.orientation = self.getGlobalOrientation()
        self.length = self.getLength()

    def getGlobalOrientation(self):
        '''
        Calculates the objects orientation and makes a variable called orientation.
        :return:
        '''

        dz = self.z2 - self.z1
        dx = self.x2 - self.x1

        if dz != 0:
            if dx != 0:
                theta_ref_global = -np.arctan2(dz, dx)
            elif self.z1 > self.z2:
                theta_ref_global = np.pi / 2
            else:
                theta_ref_global = -np.pi / 2
        else:
            if self.x1 < self.x2:
                theta_ref_global = 0
            else:
                theta_ref_global = np.pi
        return theta_ref_global

    def getLength(self):
        '''
        Calculates the objects length. This function is used when initializing the beam.
        :return: length of beam.
        '''
        dx = self.x1 - self.x2
        dz = self.z1 - self.z2
        return np.sqrt(dx * dx + dz * dz)

    def makePipe(self, r, ratioAir):
        '''
        Turns the beam object into a pipe, and calculates corresponding area and moment of inertia.
        :param r: radius
        :param ratioAir: fraction of radius that is air (hollow pipe)
        :return: nothing
        '''
        I = np.pi / 4 * (r ** 4 - ratioAir ** 4)
        A = np.pi * (r ** 2 - ratioAir ** 2)
        self.area = A
        self.momentOfInertiaStrong = I
        self.momentOfInertiaWeak = I
        self.diameter = 2 * r
        self.Zc = r

    def makeIPE(self,H, b,t_mid,t_f):
        '''
        Turns the beam object into a IPE profle,
        and calculates corresponding area and moment of inertia.
        :param H: height
        :param w_top: width top
        :param w_bot: width bottom
        :param w_mid: width middle
        :param t_top: thickness top
        :param t_bot: thickness bottom
        :return: nothing
        '''
        A_f = b * t_f
        Amid = (H - 2 * t_f) * t_mid

        I_f= b * t_f**3 / 12
        Imid = (H- 2 * t_f)**3 * t_mid / 12

        area = 2 * A_f + Amid
        momInertiaStrong = 2 * I_f + Imid + 2 * A_f * (H/2 - t_f/2)**2
        momInertiaWeak = 2 * b**3 * t_f + t_mid**3 * (H - 2 * t_f)

        self.area = area
        self.momentOfInertiaStrong = momInertiaStrong
        self.momentOfInertiaWeak = momInertiaWeak
        self.diameter = b
        self.Zc = H / 2

    def giveNumber(self, beamNumber):
        '''
        Enumerates the beam.
        :param beamNumber: beam number
        :return: nothing
        '''
        self.number = beamNumber

    def makeStiffness(self, materialArray):
        '''
        Appends Young's modulus, stiffness and yield stress to beam.
        :param E: Young's modulus
        :return: nothing
        '''
        self.E = materialArray[0]
        self.sigmafy = materialArray[1]
        self.stiffnessStrongAxis = self.momentOfInertiaStrong * self.E
        self.stiffnessWeakAxis = self.momentOfInertiaWeak * self.E

    def makeLocalStiffnessMatrix(self):
        '''
        Makes an element stiffness matrix, which requires predefined area, moment of inertia and E modulus.
        :return: nothing
        '''
        localStiffnessMatrix = np.zeros((6, 6))

        # Uses local variables to make variables that simplify later calculations
        EA_L = self.area * self.E / self.length
        EI_L = self.momentOfInertiaStrong * self.E / self.length
        EI_L2 = self.momentOfInertiaStrong * self.E / (self.length ** 2)
        EI_L3 = self.momentOfInertiaStrong * self.E / (self.length ** 3)

        # All E*A/L dependent matrix elements
        localStiffnessMatrix[0][0] = EA_L
        localStiffnessMatrix[3][3] = EA_L
        localStiffnessMatrix[0][3] = -1 * EA_L
        localStiffnessMatrix[3][0] = -1 * EA_L

        # All E*I/L dependent matrix elements
        localStiffnessMatrix[2][2] = 4 * EI_L
        localStiffnessMatrix[5][5] = 4 * EI_L
        localStiffnessMatrix[5][2] = 2 * EI_L
        localStiffnessMatrix[2][5] = 2 * EI_L

        # All E*A/L^2 dependent matrix elements
        localStiffnessMatrix[4][5] = 6 * EI_L2
        localStiffnessMatrix[5][4] = 6 * EI_L2
        localStiffnessMatrix[1][5] = -6 * EI_L2
        localStiffnessMatrix[5][1] = -6 * EI_L2
        localStiffnessMatrix[4][2] = 6 * EI_L2
        localStiffnessMatrix[2][4] = 6 * EI_L2
        localStiffnessMatrix[1][2] = -6 * EI_L2
        localStiffnessMatrix[2][1] = -6 * EI_L2

        # All E*A/L^3 dependent matrix elements
        localStiffnessMatrix[1][1] = 12 * EI_L3
        localStiffnessMatrix[4][4] = 12 * EI_L3
        localStiffnessMatrix[1][4] = -1 * 12 * EI_L3
        localStiffnessMatrix[4][1] = -1 * 12 * EI_L3

        self.localStiffnessMatrix = localStiffnessMatrix

    def transformElementStiffnessMatrix(self):
        '''
        Transforms the beams local stiffnessmatrix to the global orientation, 
        making it ready to be put in the global stiffness matrix.
        :return: applies the stiffness matrix in the global orientation to the beams
        '''
        a = self.orientation

        self.T = np.array([
            [np.cos(a), np.sin(a), 0, 0, 0, 0],
            [-np.sin(a), np.cos(a), 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, np.cos(a), np.sin(a), 0],
            [0, 0, 0, -np.sin(a), np.cos(a), 0],
            [0, 0, 0, 0, 0, 1]])

        self.T_transponent = np.array([
            [np.cos(a), -np.sin(a), 0, 0, 0, 0],
            [np.sin(a), np.cos(a), 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, np.cos(a), -np.sin(a), 0],
            [0, 0, 0, np.sin(a), np.cos(a), 0],
            [0, 0, 0, 0, 0, 1]])
        self.transformedStiffnessMatrix = list(
            np.matmul(self.T, np.matmul(self.localStiffnessMatrix, self.T_transponent)))

    def addDistributedNormalLoad(self, load):
        '''
        Adds the distributed load from imput file to the beam. Adds q1 and q2 to the beam object,
        where q1 is the normal load working on NODE1, and q2 for NODE2.
        :param load: list with load data, in global orientation
        :return:
        '''
        if (np.sin(self.orientation) == 0):
            self.q1 = load[2]
            self.q2 = load[4]
        elif (np.cos(self.orientation) == 0):
            self.q1 = load[1]
            self.q2 = load[3]
        else:
            self.q1 = load[1] / np.sin(self.orientation) + load[2] / np.cos(self.orientation)
            self.q2 = load[3] / np.sin(self.orientation) + load[4] / np.cos(self.orientation)

    def calculateFixedSupport(self):
        '''
        Calculates Fixed Support for beams affected by distributed loads
        Adds Fixed Support to the nodes affected
        Based on table found in "TMR4167 Marin teknikk 2 – Konstruksjoner - Del 1", page 281
        :return: nothing
        '''

        # If q1 or q2 is not defined an Attribute error will occur
        try:
            m1 = (1 / 20) * self.q1 * (self.length) ** 2 + (1 / 30) * self.q2 * (self.length) ** 2
            m2 = -(1 / 30) * self.q1 * (self.length) ** 2 - (1 / 20) * self.q2 * (self.length) ** 2

            self.node1.M += m1
            self.node2.M += m2

            v2 = (m1 + m2 - (self.q1 * self.length ** 2) / 6 - (self.q2 * self.length ** 2) / 3) / self.length
            v1 = -(self.length / 2) * (self.q1 + self.q2) - v2

            self.node1.Fx += v1 * np.sin(self.orientation)
            self.node1.Fz += v1 * np.cos(self.orientation)
            self.node2.Fx += v2 * np.sin(self.orientation)
            self.node2.Fz += v2 * np.cos(self.orientation)
        except AttributeError:
            pass

    def printBeam(self):
        '''
        Prints beam data in a more readable way
        :return: nothing
        '''

        string = f'Beam {self.number} from node {self.node1.number} to {self.node2.number}, Ø: {round(self.orientation * 180 / np.pi, 2)}, L: {round(self.length, 2)}\n'
        for i in range(2):
            string += f'u{i + 1}: {round(self.localDisplacements[i * 3], 4)} \tN{i + 1}: {round(self.reactionForces[i * 3], 2)}\n'
            string += f'w{i + 1}: {round(self.localDisplacements[i * 3 + 1], 4)} \tV{i + 1}: {round(self.reactionForces[i * 3 + 1], 2)}\n'
            string += f'ø{i + 1}: {round(self.localDisplacements[i * 3 + 2], 4)} \tM{i + 1}: {round(self.reactionForces[i * 3 + 2], 2)}\n'
        print(string)

    def scaleDistributedLoad(self, referenceDiameter):
        '''
        Scales forces linearly with respect to reference diameter
        :param referenceDiameter: a float thar reperesent the reference diameter
        :return: nothing
        '''
        try:
            self.q1 *= self.diameter / referenceDiameter
            self.q2 *= self.diameter / referenceDiameter
        except AttributeError:
            pass

    def calculateMaxMoment(self):
        '''
            Denne er feil.
        :return:
        '''
        try:
            # Max moment where sheer = 0, and where the resultant force of the distributed load attacks
            self.L_R = ((self.q1 / 2 + self.q2) * 2 * self.length) / ((self.q1 + self.q2) * 3)
            q_R = (self.q1 - self.q1 * self.L_R / self.length) + self.q2 * self.L_R / self.length
            if abs(self.q1) < abs(q_R):
                self.M_Max = -self.reactionForces[2] - self.reactionForces[1] * self.L_R - (
                            self.q1 * self.L_R ** 2) / 2 - (q_R - self.q1) * (self.L_R ** 2) / 6
            else:
                self.M_Max = -self.reactionForces[2] - self.reactionForces[1] * self.L_R - (q_R * self.L_R ** 2) / 2 - (
                            self.q1 - q_R) * (self.L_R ** 2) / 3

        except AttributeError:
            pass

    def printBeamMoments(self):
        '''
        Prints moment values to terminal.
        :return: nothing
        '''
        try:
            string = f'M1: {round(self.reactionForces[2])}, M_max: {round(self.M_Max)} at {round(self.L_R, 3)}, M2: {round(self.reactionForces[5])}'
        except AttributeError:
            string = f'M1: {round(self.reactionForces[2])}, No ditributed load, M2: {round(self.reactionForces[5])}'
        print(string)

    def calculateMaxBendingTension(self):
        '''
        Calculates maximun bending tension and security factor.
        :return: nothing
        '''
        try:
            self.sigmax = max([abs(self.reactionForces[2]), abs(self.M_R),
                               abs(self.reactionForces[5])]) * self.Zc / self.momentOfInertiaStrong
        except AttributeError:
            self.sigmax = max(
                [abs(self.reactionForces[2]), abs(self.reactionForces[5])]) * self.Zc / self.momentOfInertiaStrong
        self.securityFactor = self.sigmax / self.sigmafy

    def printSecurityFactor(self):
        '''
        Prints the security factor to terminal.
        :return: nothing
        '''
        print(f'Sigma_x: {round(self.sigmax / 10 ** 6)} [MPa], {round(self.securityFactor * 100)}% of fy\n\n')

    def getMomentDiagram(self):
        '''
        Denne er feil.
        :return:
        '''
        try:
            L = self.length
            N = 20
            x = np.linspace(0, L, N + 1)
            dx = L / N
            momentList = [0] * (N + 1)

            for i in range(N + 1):
                M1 = self.q1 * (i * dx) / (6 * L) * (2 * L * L - 3 * L * (i * dx) + (i * dx) ** 2)
                M2 = self.q2 * (L - i * dx) / (6 * L) * (2 * L * L - 3 * L * (L - i * dx) + (L - i * dx) ** 2)
                momentList[i] += self.reactionForce[2] + M1 + M2

            return momentList, np.array(x)
        except AttributeError:
            pass


class Node:
    '''
    This class is used to create beams and store
    data about displacements and fixed supported forces and moments.
    We use these to create the R-vector.
    '''

    def __init__(self, x, z, supportX, supportZ, supportTheta):
        '''
        Initializes a Node object.
        :param x: x coordinate
        :param z: z coordinate
        :param supportX: The node's support condition in X
        :param supportZ: The node's support condition in Z
        :param supportTheta: The node's rotational support condition
        '''
        # Coordinates
        self.x = x
        self.z = z
        # Displacement and rotation support conditions
        self.supportX = supportX
        self.supportZ = supportZ
        self.supportTheta = supportTheta
        # Fixed support forces and moment
        self.Fx = 0
        self.Fz = 0
        self.M = 0

    def giveNumber(self, nodeNumber):
        '''
        Enumerates the node. Not sure if we will use this or not.
        :param nodeNumber: number of node
        :return: applies number to node.
        '''
        self.number = nodeNumber

    def addNodeLoad(self, nodeLoad):
        '''
        Adds the nodeloads directly to the nodes
        param nodeLoad: vector containing the nodeloads
        -= because we want the reactionforces from the loads
        '''
        self.Fx -= nodeLoad[1]
        self.Fz -= nodeLoad[2]
        self.M -= nodeLoad[3]
