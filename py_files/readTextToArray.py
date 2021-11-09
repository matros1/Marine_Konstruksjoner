'''
Developed by Matias Rosenlund, Christian Lindahl Elseth og Hugo Furnes
Norwegian University of Science and Technology, Department of Marine Technology
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
            # På u,w,theta er 
            # 0 = "Denne er ikke fastholdt"
            # 1 = "Denne er ukjent, må tas med i matrisen"
            # 2 = "Denne er fastholdt"
        elif line[0] == "BEAM":
            BEAM.append(line[1:])
            # N1, N2, M, G, n
            # N1-Node lengst til venstre
            # N2-Node lengst til høyre
            # M-Materiale, 1=stål, 2=aluminium
            # Geometri. 2 = pipe, 1 = IPE
            # n utvalgt fra biblioteket
        elif line[0] == "MATERIAL":
            MATERIAL.append(line[1:])
            # E-modul, flytspenning
            # 1 er stål og 2 er aluminium
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
