from importedLibraries import *

def readCSV_float(filename):
    '''
    Reads from file and makes an float array of list elements.
    :param filename: .txt file
    :return: np.array() of flots
    '''
    f = open(filename,'r')
    NODE = f.readlines()
    for i,linje in enumerate(NODE):
        NODE[i] = linje.split(',')
    f.close()
    return np.array(NODE, float)

def readCSV_int(filename):
    '''
    Reads from file and makes an int array of list elements.
    :param filename: .txt file
    :return: np.array() of flots
    '''
    f = open(filename,'r')
    NODE = f.readlines()
    for i,linje in enumerate(NODE):
        NODE[i] = linje.split(',')
    f.close()
    return np.array(NODE, int)


def readAll():
    '''
    Reads nodedata, beamdata, materialdata, nodeloaddata and beamloaddata to np array.
    :return: np array of nodedata, beamdata, materialdata, nodeloaddata and beamloaddata.
    '''
    NODE = readCSV_float("NodeData.txt")
    # x, z, u, w, theta
    # På u,w,theta er 1="Denne er ukjent, må tas med i matrisen", 0="Denne er 0"

    BEAM = readCSV_int("BeamData.txt")
    # N1, N2, M, G
    # N1-Node lengst til venstre
    # N2-Node lengst til høyre
    #  M-Materiale, 1=stål, 2=aluminium
    #  Geometri. 1 = pipe, 2 = IPE
    #  Geometri utvalgt fra biblioteket

    MATERIAL = readCSV_float("MaterialData.txt")
    # E-modul, flytspenning, tverrkontraksjonstallet
    # 1 er stål og 2 er aluminium

    NODELOAD = readCSV_int("NodeLoadData.txt")
    # Nodenr, Fx, Fy, Fz, Mx, My, Mz

    BEAMLOAD = readCSV_int("BeamLoadData.txt")


    PipeData = readCSV_float("PipeData")
    # radius, ratio air

    IPEData = readCSV_float("IPEData")
    #total height, width top, bottom, mid, thickness topp, bot

    return NODE,BEAM,MATERIAL,NODELOAD,BEAMLOAD, PipeData, IPEData
