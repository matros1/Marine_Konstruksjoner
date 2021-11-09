

def secondDegreePolynomial(c0, c1, c2, x):
    return c0 + c1 * x + c2 * x ** 2

def thiredDegreePolynomial(c0, c1, c2, c3, x):
    return c0 + c1 * x + c2 * x ** 2 + c3 * x ** 3


def printSubplotOfMomentDiagram(beamsObjectList):
    plotList = []
    beamList = []
    for beam in beamsObjectList:
        try:
            X = [0, beam.L_R, beam.length]
            Y = [-beam.reactionForces[2], beam.M_Max, beam.reactionForces[5]]

            polyfit = np.polyfit(X, Y, 3)
            if (beam.M_Max>0):
                ffit = np.clip(np.polynomial.polynomial.polyval(X,Y), a_min=beam.reactionForces[5], a_max=beam.M_Max)
            else:
                ffit = np.clip(np.polynomial.polynomial.polyval(X, Y), a_min=beam.M_Max, a_max=beam.reactionForces[5])
            plotList.append(polyfit)
            beamList.append(beam)
        except AttributeError:
            pass

    # print(plotList)
    fig, axs = plt.subplots(1,len(plotList))
    j = 0
    N = 30

    for beam in beamList:
        x = np.linspace(0, beam.length, N)
        print(plotList[j])
        axs[j].plot(x, 10**(-6)*thiredDegreePolynomial(plotList[j][3],plotList[j][2], plotList[j][1], plotList[j][0], x))
        axs[j].set_title('Beam number ' + str(beam.number))
        axs[j].set_xlabel('x [ m ]')
        axs[j].plot(beam.L_R,beam.M_Max * 10 ** (-6), 'o', c='r', label='Max moment')

        j += 1
    axs[0].set_ylabel('Moment [ MNm ]')
    plt.show()