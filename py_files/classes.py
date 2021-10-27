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
    def __init__(self, NODE1, NODE2):
        self.node1 = NODE1
        self.x1 = NODE1.x
        self.z1 = NODE1.z
        self.node2 = NODE2
        self.x2 = NODE2.x
        self.z2 = NODE2.z
        self.orientation = self.get_global_oriantation()
        self.length = self.get_length()

    def makePipe(self, r, ratioAir):
        '''
        Turns the beam object into a pipe, and calculates corresponding area and moment of inertia.
        :param r: radius
        :param ratioAir: how much of the radius is air (thin walled pipe)
        :return: applies area and moment of inertia to the beam
        '''
        I = np.pi / 4 * (r**4 - ratioAir**4)
        A = np.pi * (r**2 - ratioAir**2)
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

        zc = (Atop*(H - t_top/2) + Amid*((H-t_top-t_bot)/2 + t_bot) + Abot*t_bot/2)/(Atop + Amid + Abot)

        Itop = w_top * t_top**3 / 12
        Ibot = w_bot * t_bot ** 3 / 12
        Imid = (H-t_top-t_bot)**3 * w_mid / 12

        area = Atop + Abot + Amid
        momInertiaStrong = Itop + Ibot + Imid + Atop*(H - t_top/2 - zc)**2 + Amid*((H-t_top-t_bot)/2 + t_bot - zc)**2 + Abot*(t_bot/2 - zc)**2
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
            load working on NODE1, and q2 for NODE2
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
        Adds FIM to the nodes affected
        Based on table found in "TMR4167 Marin teknikk 2 – Konstruksjoner - Del 1", page 281
        '''
        self.node1.M += (1/20)*self.q1*(self.length)**2 + (1/30)*self.q2*(self.length)**2
        self.node2.M += -(1/30)*self.q1*(self.length)**2 - (1/20)*self.q2*(self.length)**2

    #TODO
    """ def calculateFISheer(self):
        '''
        Calculates Fixed Clamping Sheerforces from distributed loads
        '''
        self.node1.Fx += np.sin(self.orientation)*(something)
        self.node1.Fy += np.cos(self.orientation)*(something)

        self.node2.Fx += np.sin(self.orientation)*(something)
        self.node2.Fy += np.cos(self.orientation)*(something) 
    """

class Node:
    '''
    This class is mostly used to create beams, but also to store
    data about displacements and fixed clamping forces and moments
    We can use this to create the R-vector
    '''
    def __init__(self, x, z, u, w, ø):
        #Coordinates
        self.x = x
        self.z = z
        #Displacements and rotations
        self.u = u
        self.w = w
        self.ø = ø
        #Fixed Clamping Forces and Moment
        self.Fx = 0
        self.Fz = 0
        self.M = 0

    def giveNumber(self,nodeNumber):
        '''
        Enumerates the node. Not sure if we will use this or not.
        :param nodeNumber: number of node
        :return: applies number to node.
        '''
        self.number = nodeNumber

    def addNodeLoad(self, nodeLoad):
        self.Fx += nodeLoad[1]
        self.Fz += nodeLoad[2]
        self.M += nodeLoad[3]

