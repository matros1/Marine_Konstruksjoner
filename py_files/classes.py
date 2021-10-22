from importedLibraries import *
from functions import *

'''
This py file creates classes that is handy to use when dealing with larger programs
It gets rid of all the mess when dealing with multiple lists,
One for all beam lengths, coords, orientations, geometry, area, moment of inertia, etc
Now all this intormation can be stored in one list of Beam objects. Much better.
'''

class Beam:
    '''
    This is a beam class. It includes multiple in house functions. Some is called when initializing and some
    can be called if you want to add info to the beam, such as geometry or number.
    '''
    def __init__(self, N1, N2):
        self.x1 = N1.x
        self.x2 = N2.x
        self.z1 = N1.z
        self.z2 = N2.z
        self.orientation = self.get_global_oriantation()
        self.length = self.get_length()
        # The following parameters might be useful later on.
        # self.load
        # self.N1
        # self.V1
        # self.M1
        # self.N2
        # self.V2
        # self.M2

    def makePipe(self, r, ratioAir):
        '''
        Turns the beam object into a pipe, and calculates corresponding area and moment of inertia.
        :param r: radius
        :param ratioAir: how much of the radius is air (thin walled pipe)
        :return: applies area and moment of inertia to the beam
        '''
        I = np.pi / 4 * (r * (1 - ratioAir)) ** 2
        A = I * 4
        self.area = A
        self.momentOfInertiaStrong = I
        self.momentOfInertiaWeak = I

    def makeIPE(self,H,w_top,w_bot,w_mid,t_top,t_bot):
        '''
        Turns the beam object into a IPE profle,
        and calculates corresponding area and moment of inertia.
        :param H: height
        :param w_top: width top
        :param w_bot: width bottom
        :param w_mid: width middle
        :param t_top: thickness top
        :param t_bot: thickness bottom
        :return: applies area and moment of inertia to the beam.
        '''
        Atop = w_top*t_top
        Abot = w_bot*t_bot
        Amid = (H-t_top-t_bot)*w_mid

        Itop = w_top * t_top**3 / 12
        Ibot = w_bot * t_bot ** 3 / 12
        Imid = (H-t_top-t_bot)**3 * w_mid / 12

        area = Atop + Abot + Amid
        momInertiaStrong = Itop + Ibot + Imid + Atop*(H/2 - t_top/2)**2 + Abot*(H/2 - t_bot/2)**2
        momInertiaWeak = (w_top**3*t_top + w_mid**3 * (H-t_top-t_bot) + w_bot**3 * t_bot)/12

        self.area = area
        self.momentOfInertiaStrong = momInertiaStrong
        self.momentOfInertiaWeak = momInertiaWeak

    def get_global_oriantation(self):
        '''
        Calculates the objects orientation and makes a variable called orientation.
        :return:
        '''

        dz = self.z2 - self.z1
        dx = self.x2 - self.x1

        if dz != 0:
            if dx != 0:
                theta_ref_global = -np.arctan(dz/dx)
            elif self.z1 > self.z2:
                theta_ref_global = np.pi / 2
            else:
                theta_ref_global = -np.pi / 2
        else:
            theta_ref_global = 0
        return theta_ref_global

    def get_length(self):
        '''
        Calculates the objects length. This function is used when initializing the beam.
        :return: length of beam.
        '''
        dx = self.x1 - self.x2
        dz = self.z1 - self.z2
        return np.sqrt(dx*dx + dz*dz)

    def giveNumber(self, beamNumber):
        '''
        Enumerates the beam. Not sure if we will use this or not.
        :param beamNumber: beam number
        :return: applies number to beam.
        '''
        self.number = beamNumber

    def makeTransformedStiffnessMatrix(self):
        '''
        Transforms the beams local stiffnessmatrix to the global orientation, 
        making it ready to be put in the global system-stiffnessmatrix.
        :return: applies the stiffnessmatrix in the global orientation to the beam

        NB!! Denne må testes
        '''
        a = self.orientation
        T = np.array([
            [np.cos(a), -np.sin(a), 0,  0,          0,          0],
            [np.sin(a), np.cos(a),  0,  0,          0,          0],
            [0,         0,          1,  0,          0,          0],
            [0,         0,          0,  np.cos(a),  -np.sin(a), 0],
            [0,         0,          0,  np.sin(a),  np.cos(a),  0],
            [0,         0,          0,  0,          0,          1]])

        T_transponent = np.array([
            [(10*np.cos(a)+5*np.cos(3*a)+np.cos(5*a))/(8*np.cos(2*a)+2*np.cos(4*a)+6),  (np.sin(a)+np.sin(3*a))/(2*np.cos(2*a)+2),  0,  0,                                                              0,                              0],
            [(-np.sin(a)-np.sin(3*a))/(2*np.cos(2*a)+2),                                np.cos(a),                                  0,  0,                                                              0,                              0],
            [0,                                                                         0,                                          1,  0,                                                              0,                              0],
            [0,                                                                         0,                                          0,  (4*np.cos(2*a)+np.cos(4*a)+3)/(6*np.cos(a)+2*np.cos(3*a)),      np.sin(2*a)/(2*np.cos(a)),      0],
            [0,                                                                         0,                                          0,  -np.sin(2*a)/(2*np.cos(a)),                                     (np.cos(2*a)+1)/(2*np.cos(a)),  0],
            [0,                                                                         0,                                          0,  0,                                                              0,                              1]])

        self.transformedStiffnessMatrix = np.matmul(T, np.matmul(self.localStiffnessMatrix, T_transponent))

    def addDistributedNormalLoad(self, load):
        '''
        Adds the distributed load to the beam
        param: load: list with load data, in global orientation
        return: adds q1 and q2 to the beam object, where q1 is the normal
            load working on N1, and q2 for N2
            Note:   positive load is difined as upwards(positive z) in the beams local orientation
        '''
        if(self.orientation == 0):
            self.q1 = load[2]
            self.q2 = load[4]
        elif(np.cos(self.orientation) == 0):
            self.q1 = load[1]
            self.q2 = load[3]
        else:
            self.q1 = load[1]/np.sin(self.orientation) + load[2]/np.cos(self.orientation)
            self.q2 = load[3]/np.sin(self.orientation) + load[4]/np.cos(self.orientation)

    def calculateFIM(self):
        '''
        Calculates FIM for beams affected by distributed loads
        Based on table found in "TMR4167 Marin teknikk 2 – Konstruksjoner - Del 1", page 281
        '''
        if(self.q1 == self.q2):
            self.m1 = (1/12)*self.q1*(self.length)**2
            self.m2 = -(1/12)*self.q1*(self.length)**2
        elif(abs(self.q1) < abs(self.q2)):
            self.m1 = (1/12)*self.q1*(self.length)**2 + (1/30)*self.q2*(self.length)**2
            self.m2 = -(1/12)*self.q1*(self.length)**2 - (1/20)*self.q2*(self.length)**2 
        else:
            self.m1 = (1/12)*self.q1*(self.length)**2 + (1/20)*self.q2*(self.length)**2
            self.m2 = -(1/12)*self.q1*(self.length)**2 - (1/30)*self.q2*(self.length)**2 
            
            
    #TODO


class Node:
    '''
    This class is mostly used to create beams, as beams require 2 nodes to be created.
    It can also be used to give displacements. But so far I dont completely understand
    how to make the element matrises. When these materices are created we can use node objects to
    add displacements in the global coordinate system.
    '''
    def __init__(self, x, z, u, w, ø):
        self.x = x
        self.z = z
        self.u = u
        self.w = w
        self.ø = ø

    def giveNumber(self,nodeNumber):
        '''
        Enumerates the node. Not sure if we will use this or not.
        :param nodeNumber: number of node
        :return: applies number to node.
        '''
        self.number = nodeNumber

    def giveDisplacement(self, u, w, theta):
        '''
        Not sure exactly how to use this yet. This gives a virtual displacement to the node
        in the global coordinate system.
        :param u: x displacement
        :param w: z displacement
        :param theta: rotational displacement (clockwise positive)
        :return: applies displacements to node.
        '''
        self.u = u
        self.w = w
        self.theta = theta

