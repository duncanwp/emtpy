__author__ = 'pard'
from nose.tools import istest, eq_


@istest
def test_mod():
    """
        To compare with fortran mod
    """
    eq_(8 % 5, 3)
    eq_(-8 % 5, 2) # Fortran gives -3
    eq_(8 % -5 , -2) # Fortran gives 3
    eq_(-8 % -5, -3)
    # So if my array is 10 long, if I want the previous element to 0 I want to end up with an index of 9:
    eq_((0-1) % 10, 9)
    # and if I want the previous element to 8 I want to end up with an index of 7:
    eq_((8-1) % 10, 7)
    # if I want the previous element to 9 I want to end up with an index of 8:
    eq_((9-1) % 10, 8)
    # if I want the next element after 9 I want to end up with an index of 0:
    eq_((9+1) % 10, 0)
    # if I want the element two back from 0 I want to end up with an index of 8:
    eq_((0-2) % 10, 8)

@istest
def test_mod_indexing():
    """
        To show that actually I don't need mod at all, as numpy arrays allow negative indexing, you can't index out of
        positive bounds but you can easily get around this
    """
    import numpy as np
    test_array = np.arange(10)
    eq_(test_array[0], 0)
    eq_(test_array[9], 9)
    # So if my array is 10 long, if I want the previous element to 0 I want to end up with an index of 9:
    eq_(test_array[0-1], 9)
    # and if I want the previous element to 8 I want to end up with an index of 7:
    eq_(test_array[8-1], 7)
    # if I want the previous element to 9 I want to end up with an index of 8:
    eq_(test_array[9-1], 8)
    # if I want the next element after 9 I want to end up with an index of 0:
    eq_(test_array[9+1-len(test_array)], 0)
    # if I want the element two back from 0 I want to end up with an index of 8:
    eq_(test_array[0-2], 8)
