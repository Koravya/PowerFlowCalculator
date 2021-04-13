from sympy import I, symbols, conjugate, im, re, diff, sin, cos, Matrix
from sympy.abc import x, y, z
import math

sym = {
    'A' : 0,
    'B' : 1,
    'C' : 2,
    'D' : 3,
    'E' : 4,
    'F' : 5,
    'G' : 6,
    'H' : 7,
    'I' : 8,
    'J' : 9,
    'K' : 10,
    'L' : 11,
    'M' : 12,
    'N' : 13,
    'O' : 14,
    'P' : 15,
    'Q' : 16,
    'R' : 17,
    'S' : 18,
    'T' : 19,
    'U' : 20,
    'V' : 21,
    'W' : 22,
    'X' : 23,
    'Y' : 24,
    'Z' : 25
}

def arrayaddr(name):
    temp = name.split('_')
    return [sym[x] for x in temp[1]]

def magnitude(complexR):
    return math.sqrt(math.pow(re(complexR),2)+math.pow(im(complexR),2))

def angle(complexR):
    return math.atan(im(complexR)/re(complexR))

def toRect(magnitude, angle):
    return magnitude*(cos(angle) + I*sin(angle))

def vars(key, count, start=0):
    return [symbols(key+list(sym.keys())[x]) for x in range(start, len(Vbus))]

def powerFlow(Ybus, Vbus, k, reactive=False, varV=None, varY=None, varS=None, varT=None):
    busCount = len(Vbus)

    varV = varV if varV is not None else vars('v', busCount)
    varY = varY if varY is not None else vars('y', busCount)
    varS = varS if varS is not None else vars('s', busCount)
    varT = varT if varT is not None else vars('t', busCount)

    p = 0

    for i in range(busCount):
        mag = varV[k]*varV[i]*varY[i]
        ang = varT[i]+varS[i]-varS[k]
        p += -1*mag*sin(ang) if reactive else mag*cos(ang)

    return p, [diff(p,x) for x in varS[1:].union(varV[1:])]

def getMismatch(schedule, busCount):
    varSchedP = vars('sp', busCount, 1)
    varSchedQ = vars('sq', busCount, 1)

    return [re(schedule[x+1]) - varSchedP[x] for x in range(len(varSchedP))] + [im(schedule[x+1]) - varSchedQ[x] for x in range(len(varSchedQ))], varSchedP, varSchedQ

Ybus = [
        [3-I*9, -2+I*6, -1+I*3, 0],
        [-2+I*6, 3.67-I*11, -0.67+I*2, -1+I*3],
        [-1+I*3, -0.67+I*2, 3.67-I*11, -2+I*6],
        [0, -1+I*3, -2+I*6, 3-I*9]
]

Vbus = [1.04, 1, 1, 1]

schedule = [None, 0.5-I*0.2, -1.0+I*0.5, 0.3-I*0.1]

runs = 1

print('<Given Data>','\n')
print('Ybus matrix: ', Ybus)
print('Vbus Matrix (With inital Guesses): ', Vbus)
print('Scheduled Power: ', schedule)

realEQ = list()
reactEQ = list()

jp = list()
jq = list()

syms = lambda equation, keyword: [x for x in equation.free_symbols if keyword in str(x)]

varV = None
varY = None
varS = None
varT = None

for i in range(1,len(Vbus)):
    p, dp = powerFlow(Ybus, Vbus, i, False, varV=varV, varY=varY, varS=varS, varT=varT)

    varV = syms(realEQ[0], 'v') if varV is None else varV
    varY = syms(realEQ[0], 'y') if varY is None else varY
    varS = syms(realEQ[0], 's') if varS is None else varS
    varT = syms(realEQ[0], 't') if varT is None else varT

    q, dq = powerFlow(Ybus, Vbus, i, True, varV=varV, varY=varY, varS=varS, varT=varT)

    realEQ.append(p[0])
    reactEQ.append(q[0])

    jp.append(p[1]+p[2])
    jq.append(q[1]+q[2])

print('\n', '<General Data>','\n')
print('Real Power: ', realEQ)
print('Reactive Power: ', reactEQ)
print('Jacobian Matrix: ', jp+jq)

for l in range(runs):
    jpCalc = jp
    jqCalc = jq
    realEqCalc = realEQ
    reactEQCalc = reactEQ

    for k in range(0, len(jpCalc)):
        svars = [(x, magnitude(Ybus[sym[str(x)[-1]]][k+1])) for x in varY]
        svars += [(x, angle(Ybus[sym[str(x)[-1]]][k+1])) for x in varT]
        svars += [(x, magnitude(Vbus[sym[str(x)[-1]]])) for x in varV]
        svars += [(x, angle(Vbus[sym[str(x)[-1]]])) for x in varS]

        realEqCalc[k] = realEqCalc[k].subs(svars)
        reactEQCalc[k] = reactEQCalc[k].subs(svars)

        for b in range(len(jpCalc[k])):
            jpCalc[k][b] = jpCalc[k][b].subs(svars)
            jqCalc[k][b] = jqCalc[k][b].subs(svars)

    Jacobian = jpCalc + jqCalc

    print('\n', '<Data for run: {}>'.format(l), '\n')
    print('Real Power: ', realEqCalc)
    print('Reactive Power: ', reactEQCalc)
    print('Jacobian Matrix: ', Jacobian)

    Jacobian = Matrix(Jacobian).inv()

    Mismatch, RealCalculated, ReactiveCalculated = getMismatch(schedule, len(Vbus))

    val = [(x, realEqCalc[sym[str(x)[-1]]-1]) for x in RealCalculated]
    val += [(x, reactEQCalc[sym[str(x)[-1]]-1]) for x in ReactiveCalculated]

    correct = Jacobian*Matrix([x.subs(val) for x in Mismatch])
    correct = list(correct.col(0))
    half = int(len(correct)/2)

    corMag = correct[0:half]
    corAngle = correct[half:]

    Vbus = [Vbus[0]] + [Vbus[x+1]+toRect(corMag[x], corAngle[x]*(180/math.pi)) for x in range(len(corMag))]

    print('Voltage Correction: ', correct)
    print('Bus Voltage: ', ['{}<{}'.format(magnitude(Vbus[x]), angle(Vbus[x])) for x in range(len(Vbus))])
