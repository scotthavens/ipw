

import multiprocessing as mp

def cube(x):
    return x**3

pool = mp.Pool(processes=4)
results = [pool.apply_async(cube, args=(x,)) for x in range(1,70000)]
print('done')
# print(results)


