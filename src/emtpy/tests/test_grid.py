__author__ = 'duncan'
from nose.tools import istest, eq_


class TestPhysicalGrid(object):

    def __init__(self):
        from emtpy.grid import PhysicalGrid
        import numpy as np
        values = np.arange(9*10*11)
        self.test_grid = PhysicalGrid((9,10,11), (5,5,5), values.reshape((9,10,11)))

    @istest
    def test_getn(self):
        eq_(self.test_grid.getn(0,0,0), 0)
        eq_(self.test_grid.getn(8,9,10), self.test_grid.no_elements-1)
        eq_(self.test_grid.getn((8 + 1) % 9, 9, 10), 100)