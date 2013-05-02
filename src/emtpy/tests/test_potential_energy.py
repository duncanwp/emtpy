__author__ = 'pard'
from nose.tools import istest, eq_
from emtpy.potential_energy import OneDWell


class TestOneDWell(object):

    def __init__(self):
        self.pot_energy = OneDWell(2.0, 1.0, (10, 10, 10), (10.0, 0.0, 0.0))

    @istest
    def test_shape(self):
        eq_(self.pot_energy[(0,0,0)],0.0)

        eq_(self.pot_energy[(0,0,0)],0.0)