from numpy.core.numeric import array_equal

__author__ = 'duncan'
from nose.tools import istest, eq_


@istest
def test_constructing_ysm():
    import numpy as np
    from emtpy.yale_sparse import YaleSparse
    dense = np.array([[3.0,0.0,1.0,0.0,0.0],
                      [0.0,4.0,0.0,0.0,0.0],
                      [0.0,7.0,5.0,9.0,0.0],
                      [0.0,0.0,0.0,0.0,2.0],
                      [0.0,0.0,0.0,6.0,5.0]])
    sp = YaleSparse.from_dense(dense,15)
    ref_ija = np.array([7,8,8,10,11,12,3,2,4,5,4])
    ref_sa = np.array([3.0,4.0,5.0,0.0,5.0,0.0,1.0,7.0,9.0,2.0,6.0])
    assert array_equal(ref_ija, sp.ija)
    assert array_equal(ref_sa, sp.sa)
