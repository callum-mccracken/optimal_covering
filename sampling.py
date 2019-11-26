"""a few functions to do with random sampling"""

import numpy as np
import random


def random_sample(iterable, k=1):
    """
    Returns k randomly selected values from iterable.

    The reason for using this rather than random.sample is that
    random.sample requires a sequence, not an iterable,
    so we'd have to store the iterable in a list or something which could get
    enormous.
    """
    result = [None] * k
    for i, item in enumerate(iterable):
        if i < k:
            result[i] = item
        else:
            j = int(random.random() * (i+1))
            if j < k:
                result[j] = item
    random.shuffle(result)
    return result


def sampled(a, b, N):
    """
    Make N random combinations of len(a) values taken from a and b,
    (assuming len(a)=len(b)),
    with a random number taken from each of a and b, without replacement.

    Note a and b may not have all unique elements from one another, but we
    will make sure the output has only unique elements.
    """
    assert len(a) == len(b)
    n = len(a)

    combos = np.empty((N, n), dtype=object)

    # to hold "taken" for each sample,
    # and to make sure we don't have two identical ones
    taken_arr = np.empty((N,n,2), dtype=object)

    # now loop through and and make N random combinations
    for i in range(N):
        # the values we've selected for this sample
        values = np.empty(n, dtype=object)
        # to ensure we don't repeat elements in future samples, keep track
        taken = np.array([[False, False]]*n)
        keep_picking = True
        while keep_picking:
            keep_picking = True
            # pick n different elements
            for j in range(n):
                # select value
                a_or_b = np.random.randint(0, 2)  # 1 or 0
                take_from = a if a_or_b else b
                index = np.random.randint(0, n)
                while taken[index][a_or_b]==True or (take_from[index] in values):
                    a_or_b = np.random.randint(0, 2)
                    take_from = a if a_or_b else b
                    index = np.random.randint(0, n)
                    if np.all(taken):
                        raise ValueError("No more points to take!")
                # "take" value
                values[j] = take_from[index]
                taken[index][a_or_b] = True
            # put taken value in taken_arr
            if not any([(taken == t).all() for t in taken_arr]):
                taken_arr[i, :, :] = taken
                keep_picking = False
        combos[i] = values
    return combos


if __name__ == "__main__":
    a = [1,2,3,4]
    b = [5,6,7,8]
    c = sampled(a, b, 10)
    print(c)

    print(random_sample(range(1000), 10))