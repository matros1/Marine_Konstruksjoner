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
    NODE = readCSV_float("txt_files/NodeData.txt")
    # x, z, u, w, theta
    # På u,w,theta er 1 = "Denne er ukjent, må tas med i matrisen", 0 = "Denne skal være 0"

    BEAM = readCSV_int("txt_files/BeamData.txt")
    # N1, N2, M, G
    # N1-Node lengst til venstre
    # N2-Node lengst til høyre
    #  M-Materiale, 1=stål, 2=aluminium
    #  G-Tverrsnittsgeometri

    MATERIAL = readCSV_float("txt_files/MaterialData.txt")
    # E-modul, f_y, tverrkontraksjonstallet
    # 1 er stål og 2 er aluminium

    NODELOAD = readCSV_int("txt_files/NodeLoadData.txt")
    # Nodenr, Fx, Fy, Fz, Mx, My, Mz

    BEAMLOAD = readCSV_int("txt_files/BeamLoadData.txt")

    return NODE,BEAM,MATERIAL,NODELOAD,BEAMLOAD
