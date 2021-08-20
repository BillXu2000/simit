import sys
with open(sys.argv[1] + ".1.node", "r") as fi:
    a = fi.readlines()[1: -1]

for i in a:
    i = i.split(" ")
    k = []
    for j in i:
        if j.strip() != '': k.append(j.strip())
    print('', end = "v ")
    [print(i, end = " ") for i in k[1:]]
    print('')

with open(sys.argv[1] + ".1.face", "r") as fi:
    a = fi.readlines()[1: -1]

for i in a:
    i = i.split(" ")
    k = []
    for j in i:
        if j.strip() != '': k.append(j.strip())
    print('', end = "f ")
    [print(int(i) + 1, end = " ") for i in k[1:-1]]
    print('')
