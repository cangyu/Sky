import os
import multiprocessing
import time


def f(idx, blk, dim):
    u, v, w = dim
    print("PID:", os.getpid())
    print("Current blk dimension: {} x {} x {}".format(u, v, w))
    print("Building...")
    time.sleep(20)
    print("Block {} calculation done.".format(idx))


if __name__ == '__main__':
    for k in range(10):
        print("Sample ", k)
        print('Main:', os.getpid())
        cores = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes=cores)
        for i in range(3):
            pool.apply_async(f, (i, None, (i, i + 1, i + 2)))
        pool.close()
        pool.join()

        print("Another stuff after pool.")
