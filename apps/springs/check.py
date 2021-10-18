import numpy as np, re
for i in range(1, 5):
    with open(f'{i}.obj', 'r') as fi:
        lines = fi.readlines()
    ans = []
    for line in lines:
        num = re.findall(r'\S+', line)
        if num[0] == 'v':
            for j in num[1:]:
                ans.append(float(j))
    arr = np.array(ans)
    print(arr.shape)
    print(i, arr.mean(), (arr**2).mean())


