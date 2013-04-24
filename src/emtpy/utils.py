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