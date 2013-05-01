__author__ = 'pard'
from nose.tools import istest, eq_
from emtpy.potential_energy import APotentialEnergy
from emtpy.material_distribution import MaterialDistribution


class OneDWell(APotentialEnergy):

    def __init__(self, width, depth, shape, size):
        from numpy import fromfunction

        # Create a new grid for the potential
        super(Harmonic, self).__init__(shape, size)

        # Define r_0, the center of the well, to be the middle of the grid
        r_0 = map(lambda x: x / 2.0, size[0])

        def _single_well(i, j, k):
            """
                Function to evaluate the quantum well potential on any index in a 3-D array
            """
            from emtpy.utils import heaviside
            x = self.coord((i, j, k))
            return 0.0 - depth*heaviside(abs(x-r_0)-width/2.0)

        self.values = fromfunction(_single_well, shape)


class Harmonic(APotentialEnergy):

    def __init__(self, omega, shape, size):
        from emtpy.constants import eV, me
        from numpy import fromfunction

        # Create a new grid for the potential
        super(Harmonic, self).__init__(shape, size)

        # Define the normalization constant
        const = ((0.5*me*1E-18)/eV)

        # Define r_0, the center of the parabola, to be the middle of the grid
        r_0 = map(lambda x: x / 2.0, size)

        def _harmonic_oscillator(i, j, k):
            """
                Function to evaluate the harmonic potential on any index in a 3-D array
            """
            x = self.coord((i, j, k))
            pos = map(lambda x, r: (x-r)**2, x, r_0)
            return const*(omega**2)*sum(pos)

        self.values = fromfunction(_harmonic_oscillator, shape)


class MockMaterialDistribution(MaterialDistribution):

    def __init__(self, eff_mass=(1.0, 1.0)):
        self.eff_mass = eff_mass
        self.increments = (1, 1, 1)
        self.size = (1, 1, 1)

    def inv_mass_xy(self, idx):
        return 1.0/self.eff_mass[0]

    def inv_mass_z(self, idx):
        return 1.0/self.eff_mass[1]


class EngineTests(object):

    def __init__(self, engine):
        self.TestEngine = engine

    def setup_harmonic(self):
        self.potential = Harmonic(1.0, (100, 100, 100), (10.0, 10.0, 10.0))
        self.mat_dist = MockMaterialDistribution()
        self.engine = self.TestEngine(self.mat_dist, self.potential)

    def setup_one_d_well(self):
        self.potential = OneDWell(2.0, 1.0, (100, 1, 1), (10.0, 0.0, 0.0))
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
        vals, vectors = self.engine.solve()
        eq_(vals[0], hbarev*1)
        eq_(vals[1], hbarev*2)


class ArpackTests(EngineTests):

    def __init__(self):
        from emtpy.engine import SparseSolver
        super(ArpackTests, self).__init__(SparseSolver)
