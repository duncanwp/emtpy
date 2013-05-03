__author__ = 'pard'
from nose.tools import istest, eq_
from emtpy.potential_energy import OneDWell, Harmonic


class TestOneDWell(object):

    def __init__(self):
        self.pot_energy = OneDWell(2.0, 1.0, (10,), (10.0,))

    @istest
    def test_shape(self):
        eq_(self.pot_energy.values[(0,)],0.0)

        eq_(self.pot_energy.values[(5,)],1.0)



class TestHarmonic(object):

    def __init__(self):
        self.pot_energy = Harmonic(1.0, (10,10,10), (10.0,10.0,10.0))

    @istest
    def test_shape(self):
        from emtpy.constants import eV, me
        const = ((0.5*me*1E-18)/eV)

        eq_(self.pot_energy.values[(5,5,5)],0.0)
        eq_(self.pot_energy.values[(0,0,0)],25.0*const)
        eq_(self.pot_energy.values[(10,10,10)],25.0*const)
