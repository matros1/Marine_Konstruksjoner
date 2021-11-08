'''
Developed by Matias Rosenlund, Christian Lindahl Elseth og Hugo Furnes
Norwegian University of Science and Technology, Department of Marine Technology
08.11.2021 as part of TMR4167 Marin teknikk - Konstruksjoner
'''



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
    :return: np.array() of int's
    '''
    f = open(filename, 'r')
    NODE = f.readlines()
    for i, linje in enumerate(NODE):
        NODE[i] = linje.split(',')
    f.close()
    return np.array(NODE, int)


def readInputFile(filepath):
    '''
    Reads input file and sorts info into appropriate array
    :param filepath: Lists containing appropriate data
    :return: Arrays of input file data
    '''
    f = open(filepath, 'r')
    lineList = f.readlines()
    NODE, BEAM, MATERIAL, NODELOAD, BEAMLOAD, PIPE, IPE = [], [], [], [], [], [], []
    # Set default reference diameter
    REFERENCEDIAMETER = 1
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
        elif line[0] == "REFERENCEDIAMETER":
            REFERENCEDIAMETER = float(line[1].strip())
        elif line[0] == "---\n":
            pass
        else:
            print(f"Error reading line {i + 1} in input file!")
    f.close()
    return np.array(NODE, dtype=float), np.array(BEAM, dtype=float), np.array(MATERIAL, dtype=float), \
           np.array(NODELOAD, dtype=int), np.array(BEAMLOAD, dtype=int), np.array(PIPE, dtype=float), \
           np.array(IPE, dtype=float), REFERENCEDIAMETER
