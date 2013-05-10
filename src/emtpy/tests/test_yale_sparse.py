__author__ = 'duncan'
from nose.tools import istest, eq_
from numpy.testing.utils import assert_array_equal
import numpy as np
from emtpy.yale_sparse import YaleSparse, YaleSparseSymmetric


class TestYaleSparse(object):

    def __init__(self):
        self.dense = np.array([[3.0, 0.0, 1.0, 0.0, 0.0],
                               [0.0, 4.0, 0.0, 0.0, 0.0],
                               [0.0, 7.0, 5.0, 9.0, 0.0],
                               [0.0, 0.0, 0.0, 0.0, 2.0],
                               [0.0, 0.0, 0.0, 6.0, 5.0]])
    @istest
    def test_constructing_ysm(self):
        sp = YaleSparse.from_dense(self.dense, 11)
        ref_ija = np.array([6, 7, 7, 9, 10, 11, 2, 1, 3, 4, 3])
        ref_sa = np.array([3.0, 4.0, 5.0, 0.0, 5.0, 0.0, 1.0, 7.0, 9.0, 2.0, 6.0])
        assert_array_equal(ref_ija, sp.ija)
        assert_array_equal(ref_sa, sp.sa)

    @istest
    def test_deconstructing_ysm(self):
        sp = YaleSparse.from_dense(self.dense, 11)
        new_dense = sp.to_dense((5, 5))
        assert_array_equal(new_dense, self.dense)


    @istest
    def test_ysm_right_multiplication(self):
        test_vector = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        sp = YaleSparse.from_dense(self.dense, 11)
        assert_array_equal(sp.av(test_vector), np.dot(self.dense, test_vector))

class TestYaleSparseSymmetric(object):

    def __init__(self):
        self.dense = np.array([[3.0, 0.0, 1.0, 0.0, 0.0],
                               [0.0, 4.0, 7.0, 0.0, 0.0],
                               [1.0, 7.0, 5.0, 9.0, 0.0],
                               [0.0, 0.0, 9.0, 0.0, 6.0],
                               [0.0, 0.0, 0.0, 6.0, 5.0]])
    @istest
    def test_constructing_ysm(self):
        sp = YaleSparseSymmetric.from_dense(self.dense, 10)
        ref_ija = np.array([6, 7, 8, 9, 10, 10, 2, 2, 3, 4])
        ref_sa = np.array([3.0, 4.0, 5.0, 0.0, 5.0, 0.0, 1.0, 7.0, 9.0, 6.0])
        assert_array_equal(ref_ija, sp.ija)
        assert_array_equal(ref_sa, sp.sa)

    @istest
    def test_deconstructing_ysm(self):
        sp = YaleSparseSymmetric.from_dense(self.dense, 10)
        new_dense = sp.to_dense((5, 5))
        assert_array_equal(new_dense, self.dense)


    @istest
    def test_ysm_right_multiplication(self):
        test_vector = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        sp = YaleSparseSymmetric.from_dense(self.dense, 10)
        assert_array_equal(sp.av(test_vector), np.dot(self.dense, test_vector))

