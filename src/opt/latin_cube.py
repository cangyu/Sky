import numpy as np


def latin_sample(entry_list):
    """
    拉丁超立方抽样
    """

    n = len(entry_list)
    if n == 0:
        raise ValueError("Invalid input.")

    '''Defensive check'''
    m = len(entry_list[0])
    for k, sp in enumerate(entry_list):
        if len(sp) != m:
            raise AssertionError("Inconsistent length of entry {}.".format(k))

    rp = np.empty((n, m), int)
    for i in range(n):
        rp[i] = np.random.permutation(m)

    ret = []
    for i in range(m):
        cur_comb = []
        for j in range(n):
            cur_comb.append(entry_list[j][rp[j][i]])
        ret.append(cur_comb)

    return ret
