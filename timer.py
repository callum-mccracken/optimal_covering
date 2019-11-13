"""
Contains a simple timer function, to be called as a decorator.

Usage:

# import
from timer import timeit

# define function with timeit used as decorator
@timeit
def some_func(x):
    ...

# call the function
some_func(5)

# prints function name and execution time
>>> "some_func"  12.09 ms
"""
from time import time


def timeit(method):
    def timed(*args, **kw):
        ts = time()
        result = method(*args, **kw)
        te = time()
        if 'log_time' in kw:
            name = kw.get('log_name', method.__name__.upper())
            kw['log_time'][name] = int((te - ts) * 1000)
        else:
            print('%r  %2.2f ms' % (method.__name__, (te - ts) * 1000))
        return result
    return timed
