__author__ = 'pard'
from nose.tools import istest, eq_
from emtpy.potential_energy import Harmonic, OneDWell
from emtpy.material_distribution import MaterialDistribution

class MockMaterialDistribution(MaterialDistribution):

    def __init__(self, eff_mass=None):
        from emtpy.constants import me
        self.eff_mass = (0.1, 0.1) if eff_mass is None else eff_mass
        self.increments = (1, 1, 10)
        self.size = (1, 1, 10)
        self.shape = (1, 1, 100)

    def inv_mass(self, idx, dim):
        return 1.0/self.eff_mass[dim]

    def inv_mass_xy(self, idx):
        return self.inv_mass(idx, 0)

    def inv_mass_z(self, idx):
        return self.inv_mass(idx, 1)


class EngineTests(object):

    def __init__(self, engine):
        self.TestEngine = engine

    def setup_harmonic(self):
        self.potential = Harmonic(1.0, (100, 100, 100), (10.0, 10.0, 10.0))
        self.mat_dist = MockMaterialDistribution()
        self.engine = self.TestEngine(self.mat_dist, self.potential)

    def setup_one_d_well(self):
        self.potential = OneDWell(4.0, 1.0, (100,), (10.0,))
        self.mat_dist = MockMaterialDistribution()
        self.engine = self.TestEngine(self.mat_dist, self.potential)

    def setup_three_d_one_d_well(self):
        self.potential = OneDWell(2.0, 1.0, (1, 1, 100), (1.0, 1.0, 10.0))
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
        self.setup_three_d_one_d_well()
        vals, vectors = self.engine.solve()
        #pot = self.potential.values[0,0,:]
        plt.plot(self.potential.values[0,0,:])
        plt.plot(vectors[:, 0]+vals[0])
        plt.plot(vectors[:, 1]+vals[1])
        plt.show()
        eq_(vals[0], hbarev*1)
        eq_(vals[1], hbarev*2)



class ArpackTests(EngineTests):

    def __init__(self):
        from emtpy.engine import ARPACKSolver
        super(ArpackTests, self).__init__(ARPACKSolver)

    @istest
    def test_one_d_well_energies(self):
        super(ArpackTests, self).test_one_d_well_energies()


class ArpackOriginalTests(EngineTests):

    def __init__(self):
        from emtpy.engine import ARPACKOriginal
        super(ArpackOriginalTests, self).__init__(ARPACKOriginal)

    @istest
    def test_one_d_well_energies(self):
        super(ArpackOriginalTests, self).test_one_d_well_energies()