__author__ = 'pard'
from nose.tools import istest, eq_
from emtpy.potential_energy import Harmonic, OneDWell
from emtpy.material_distribution import MaterialDistribution

class MockMaterialDistribution(MaterialDistribution):

    def __init__(self, eff_mass=None):
        from emtpy.constants import me
        self.eff_mass = (1.0*me, 1.0*me) if eff_mass is None else eff_mass
        self.increments = (1, 1, 1)
        self.size = (1, 1, 1)

    def inv_mass(self, idx, dim):
        return 1.0/self.eff_mass[dim]


class EngineTests(object):

    def __init__(self, engine):
        self.TestEngine = engine

    def setup_harmonic(self):
        self.potential = Harmonic(1.0, (100, 100, 100), (10.0, 10.0, 10.0))
        self.mat_dist = MockMaterialDistribution()
        self.engine = self.TestEngine(self.mat_dist, self.potential)

    def setup_one_d_well(self):
        self.potential = OneDWell(2.0, -1.0, (100,), (10.0,))
        self.mat_dist = MockMaterialDistribution()
        self.engine = self.TestEngine(self.mat_dist, self.potential)

    @istest
    def test_harmonic_oscilator_energies(self):
        from emtpy.constants import hbarev
        self.setup_harmonic()
        vals, vectors = self.engine.solve()
        eq_(vals[0], hbarev*1)
        eq_(vals[1], hbarev*2)

    @istest
    def test_one_d_well_energies(self):
        from emtpy.constants import hbarev
        import matplotlib.pyplot as plt
        self.setup_one_d_well()
        vals, vectors = self.engine.solve()
        plt.plot(self.potential.values)
        plt.plot(vectors[:, 0]+vals[0])
        plt.show()
        eq_(vals[0], hbarev*1)
        eq_(vals[1], hbarev*2)


class ArpackTests(EngineTests):

    def __init__(self):
        from emtpy.engine import SparseSolver
        super(ArpackTests, self).__init__(SparseSolver)

    @istest
    def test_one_d_well_energies(self):
        super(ArpackTests, self).test_one_d_well_energies()