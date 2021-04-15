from sympy import symbols, diff, sin, cos, Matrix
from sympy.abc import x, y, z
import math, cmath

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

def vars(key, count, start=0):
    return [symbols(key+list(sym.keys())[x]) for x in range(start, len(Vbus))]

def powerFlow(Ybus, Vbus, k, reactive=False):
    busCount = len(Vbus)

    v = vars('v', busCount)
    y = vars('y', busCount)
    s = vars('s', busCount)
    t = vars('t', busCount)

    p = 0

    for i in range(busCount):
        mag = v[k]*v[i]*y[i]
        ang = t[i]+s[i]-s[k]
        p += mag*sin(ang) if reactive else mag*cos(ang)

    p = -1*p if reactive else p

    return p, [diff(p,x) for x in s[1:]+v[1:]]

def getMismatch(schedule, busCount):
    varSchedP = vars('sp', busCount, 1)
    varSchedQ = vars('sq', busCount, 1)

    return [schedule[x+1].real - varSchedP[x] for x in range(len(varSchedP))] + [schedule[x+1].imag - varSchedQ[x] for x in range(len(varSchedQ))], varSchedP, varSchedQ

Ybus = [
        [complex('3-9j'), complex('-2+6j'), complex('-1+3j'), complex('0')],
        [complex('-2+6j'), complex('3.67-11j'), complex('-0.67+2j'), complex('-1+3j')],
        [complex('-1+3j'), complex('-0.67+2j'), complex('3.67-11j'), complex('-2+6j')],
        [complex(0), complex('-1+3j'), complex('-2+6j'), complex('3-9j')]
]

Vbus = [complex(1.04,0), complex(1,0), complex(1,0), complex(1,0)]

schedule = [None, complex('0.5-0.2j'), complex('-1.0+0.5j'), complex('0.3-0.1j')]

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

for i in range(1,len(Vbus)):
    p, dp = powerFlow(Ybus, Vbus, i)
    q, dq = powerFlow(Ybus, Vbus, i, True)

    realEQ.append(p)
    reactEQ.append(q)

    jp.append(dp)
    jq.append(dq)

varV = [x for x in syms(realEQ[0], 'v')]
varY = [x for x in syms(realEQ[0], 'y')]
varS = [x for x in syms(realEQ[0], 's')]
varT = [x for x in syms(realEQ[0], 't')]

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
        svars = [(x, abs(Ybus[sym[str(x)[-1]]][k+1])) for x in varY]
        svars += [(x, cmath.phase(Ybus[sym[str(x)[-1]]][k+1])) for x in varT]
        svars += [(x, abs(Vbus[sym[str(x)[-1]]])) for x in varV]
        svars += [(x, cmath.phase(Vbus[sym[str(x)[-1]]])) for x in varS]

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

    Jacobian = Matrix(Jacobian)
    print(Jacobian)

    Jacobian = Jacobian.inv()

    Mismatch, RealCalculated, ReactiveCalculated = getMismatch(schedule, len(Vbus))

    val = [(x, realEqCalc[sym[str(x)[-1]]-1]) for x in RealCalculated]
    val += [(x, reactEQCalc[sym[str(x)[-1]]-1]) for x in ReactiveCalculated]

    print(Matrix([[x.subs(val) for x in Mismatch]]))

    correct = Jacobian.T*Matrix([x.subs(val) for x in Mismatch])
    correct = list(correct.col(0))
    half = int(len(correct)/2)

    corMag = correct[0:half]
    corAngle = correct[half:]

    Vbus = [Vbus[0]] + [Vbus[x+1]+cmath.rect(corMag[x], corAngle[x]) for x in range(len(corMag))]

    print('Voltage Correction: ', correct)
    print('Bus Voltage: ', ['{}<{}'.format(abs(Vbus[x]), cmath.phase(Vbus[x])) for x in range(len(Vbus))])