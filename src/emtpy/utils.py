__author__ = 'pard'


def heaviside(x):
    if x == 0.0:
        return 0.5

    return 0.0 if x < 0.0 else 1.0

def deltaij(i,j):
    return 1.0 if i==j else 0.0


def product(a_list):
    prod = 1
    for item in a_list:
        prod *= item
    return prod


def ravel_index(x, dims):
    """
        From http://stackoverflow.com/questions/5777058/mapping-a-3-dimensional-coordinate-to-a-1-dimensional-index
        Used to perform the inverse of numpy.unravel_index() in mapping an index tuple onto a single flat index.
    """
    i = 0
    for dim, j in zip(dims, x):
        i *= dim
        i += j
    return i