import argparse, yaml, math

def rect(magnitude, angle):
    return complex(magnitude*math.cos(angle),magnitude*math.sin(angle))

def getVal(entry, key):
    raw = str(entry[key] if key in list(entry.keys()) else entry[list(entry.keys())[1]]).replace(' ', '').split('<')
    val = rect(float(raw[0]), float(raw[1])) if len(raw) > 1 else complex(raw[0].replace('i' 'j').replace('I','j'))
    return val if key in list(entry.keys()) else 1/val

def getMatrix(system, key='admittance', precision=3):
    size = len(system['busses'])
    matrix = [[0 for x in range(size)] for y in range(size)]

    for bus in system['busses']:
        k = int(bus['bus'])-1
        for j in range(size):
            if j == k:
                val = 0
                for x in bus['loads']:
                    val += getVal(x,key)
                for x in bus['connections']:
                    val -= getVal(x,key)
                matrix[k][k] = complex(round(val.real,precision), round(val.imag,precision)) if type(val) is complex else round(val,precision)
            else:
                val = 0
                for con in bus['connections']:
                    if con['bus']-1 == j:
                        val = -1*getVal(con, key)
                        break
                matrix[k][j] = complex(round(val.real,precision), round(val.imag,precision)) if type(val) is complex else round(val,precision)
    
    return matrix

parser  = argparse.ArgumentParser(description='Power System Matrix Generator')
parser.add_argument('file', help='System Definition file in yaml format')
parser.add_argument('--Ybus', action='store_const', const=True, default=False, help='Solve for Ybus matrix')
parser.add_argument('--Zbus', action='store_const', const=True, default=False, help='Solve for Zbus Matrix')
parser.add_argument('-p', '--Precision', default=3, type=int, help='Round result to given decimal place')

args = vars(parser.parse_args())

if not args['Zbus'] and not args['Ybus']:
    args['Zbus'] = True
    args['Ybus'] = True

with open(args['file']) as fs:
    sysDef = yaml.load(fs, Loader=yaml.FullLoader)


if args['Ybus']:
    ybus = getMatrix(sysDef, 'admittance', args['Precision'])
    for row in ybus:
        print(row)

if args['Zbus']:
    zbus = getMatrix(sysDef, 'impedance', args['Precision'])
    for row in zbus:
        print(row)