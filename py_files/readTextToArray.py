from importedLibraries import *


def readCSV_float(filename):
    '''
    Reads from file and makes an float array of list elements.
    :param filename: .txt file
    :return: np.array() of flots
    '''
    f = open(filename, 'r')
    NODE = f.readlines()
    for i, linje in enumerate(NODE):
        NODE[i] = linje.split(',')
    f.close()
    return np.array(NODE, float)


def readCSV_int(filename):
    '''
    Reads from file and makes an int array of list elements.
    :param filename: .txt file
    :return: np.array() of flots
    '''
    f = open(filename, 'r')
    NODE = f.readlines()
    for i, linje in enumerate(NODE):
        NODE[i] = linje.split(',')
    f.close()
    return np.array(NODE, int)


def readAll():
    '''
    Reads nodedata, beamdata, materialdata, nodeloaddata and beamloaddata to np array.
    :return: np array of nodedata, beamdata, materialdata, nodeloaddata and beamloaddata.
    '''
    nodeData = readCSV_float("NodeData.txt")
    # x, z, u, w, theta
    # På u,w,theta er 1="Denne er ukjent, må tas med i matrisen", 0="Denne er 0"

    beamData = readCSV_int("BeamData.txt")
    # N1, N2, M, G
    # N1-Node lengst til venstre
    # N2-Node lengst til høyre
    #  M-Materiale, 1=stål, 2=aluminium
    #  Geometri. 1 = pipe, 2 = IPE
    #  Geometri utvalgt fra biblioteket

    materialData = readCSV_float("MaterialData.txt")
    # E-modul, flytspenning, tverrkontraksjonstallet
    # 1 er stål og 2 er aluminium

    nodeLoad = readCSV_int("NodeLoadData.txt")
    # Nodenr, Fx, Fy, Fz, Mx, My, Mz

    beamDistributedLoadData = readCSV_int("BeamLoadData.txt")

    PipeData = readCSV_float("PipeData")
    # radius, ratio air

    IPEData = readCSV_float("IPEData")
    # total height, width top, bottom, mid, thickness topp, bot

    return nodeData, beamData, materialData, nodeLoad, beamDistributedLoadData, PipeData, IPEData


def readInputFile(filepath):
    '''
    Reads input file
    param filepath: filepath
    return: Lists containing apropriate data
    '''
    f = open(filepath, 'r')
    lineList = f.readlines()
    NODE, BEAM, MATERIAL, NODELOAD, BEAMLOAD, PIPE, IPE, POINTLOAD = [], [], [], [], [], [], [], []
    for i, line in enumerate(lineList):
        line = line.split(',')
        if line[0] == "NODE":
            NODE.append(line[1:])
        elif line[0] == "BEAM":
            BEAM.append(line[1:])
        elif line[0] == "MATERIAL":
            MATERIAL.append(line[1:])
        elif line[0] == "NODELOAD":
            NODELOAD.append(line[1:])
        elif line[0] == "BEAMLOAD":
            BEAMLOAD.append(line[1:])
        elif line[0] == "PIPE":
            PIPE.append(line[1:])
        elif line[0] == "IPE":
            IPE.append(line[1:])
        elif line[0] == "POINTLOAD":
            POINTLOAD.append(line[1:])
        elif line[0] == "---\n":
            pass
        else:
            print(f"Error reading line {i + 1} in input file!")
    f.close()
    return np.array(NODE, dtype=float), np.array(BEAM, dtype=float), np.array(MATERIAL, dtype=float), np.array(NODELOAD,dtype=int), \
           np.array(BEAMLOAD, dtype=int), np.array(PIPE, dtype=float), np.array(IPE, dtype=float)
