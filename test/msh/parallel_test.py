import os
import multiprocessing
import time
import numpy as np

test_list = []


def f(idx, blk, dim):
    u, v, w = dim
    print("PID:", os.getpid())
    print("Building block {} ...".format(idx))
    print("Current block dimension: {} x {} x {}".format(u, v, w))
    a = np.arange(2000 ** 2).reshape((2000, 2000))
    b = np.linalg.inv(a)
    test_list.append(b)  # resources for process are independent!!!
    print("Block {} calculation done.".format(idx))
    return b


print('Main:', os.getpid())
for k in range(20):
    core_num = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=core_num)
    for i in range(13):
        pool.apply_async(f, (i, None, (i, i + 1, i + 2)))
    pool.close()
    pool.join()
    print(len(test_list))
    print("Another stuff after pool.")
