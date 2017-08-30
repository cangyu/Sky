import os
import multiprocessing
import time


def f(idx, blk, dim):
    u, v, w = dim
    print("PID:", os.getpid())
    print("Building block {} ...".format(idx))
    print("Current block dimension: {} x {} x {}".format(u, v, w))
    time.sleep(10)
    print("Block {} calculation done.".format(idx))


if __name__ == '__main__':
    print('Main:', os.getpid())
    core_num = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=core_num)
    for i in range(13):
        pool.apply_async(f, (i, None, (i, i + 1, i + 2)))
    pool.close()
    pool.join()

    print("Another stuff after pool.")
