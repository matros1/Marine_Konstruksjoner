'''
Developed at the Norwegian University of Science and Technology, Department of Marine Technology
08.11.2021 as part of TMR4167 Marin teknikk - Konstruksjoner
'''

from importedLibraries import *

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
    REFERENCEDIAMETER = -1
    for i, line in enumerate(lineList):
        line = line.split(',')
        if line[0] == "NODE":
            NODE.append(line[1:])
            # x, z, u, w, theta
            # u,w,theta: 
            # 0 = "Not held in place, free to rotate/move"
            # 1 = "Unknown, needs to be calculated"
            # 2 = "Held in place, multiplies stiffness by 10^6"
        elif line[0] == "BEAM":
            BEAM.append(line[1:])
            # N1, N2, M, G, n
            # N1-Node 1
            # N2-Node 2
            # M-Material, 1=steel, 2=aluminium
            # Geometry. 2 = pipe, 1 = IPE
            # n geometry from pipe or IPE library
        elif line[0] == "MATERIAL":
            MATERIAL.append(line[1:])
            # E-modulus, yield stress
        elif line[0] == "NODELOAD":
            NODELOAD.append(line[1:])
            # Nodenr, Fx, Fz, My
        elif line[0] == "BEAMLOAD":
            BEAMLOAD.append(line[1:])
            # Beamnr, Fx1, Fz1, Fx2, Fz2
        elif line[0] == "PIPE":
            PIPE.append(line[1:])
            # outer radius, inner radius
        elif line[0] == "IPE":
            IPE.append(line[1:])
            # h, b, t_w, t_f
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
