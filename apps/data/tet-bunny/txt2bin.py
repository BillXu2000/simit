import numpy as np, sys

if sys.argv[1][-3:] == 'ele':
    with open(sys.argv[1], 'r') as fi:
        lines = fi.readlines()[1:-1]
        nodes = np.zeros((len(lines), 4), np.int32)
        for i, line in enumerate(lines):
            splits = line.split(' ')
            num = []
            for j in splits:
                if len(j) > 0: num.append(int(j))
            for j, k in enumerate(num[1:]):
                nodes[i, j] = k
elif sys.argv[1][-4:] == 'node':
    with open(sys.argv[1], 'r') as fi:
        lines = fi.readlines()[1:-1]
        nodes = np.zeros((len(lines), 3), np.float32)
        for i, line in enumerate(lines):
            splits = line.split(' ')
            num = []
            for j in splits:
                if len(j) > 0: num.append(float(j))
            for j, k in enumerate(num[1:]):
                nodes[i, j] = k

np.save(sys.argv[1], nodes)
