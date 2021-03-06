__author__ = 'pard'
from emtpy.tests.semi_analytic_solution import finite_square_box_energies, infinite_square_box_energy
from nose.tools import istest, nottest, eq_
from emtpy.utils import assert_almost_equal
from emtpy.potential_energy import Harmonic, OneDWell, ThreeDOneDWell


class MockMaterialDistribution(object):

    def __init__(self,grid,  eff_mass=None):
        self.__eff_mass = (1.0, 1.0) if eff_mass is None else eff_mass
        self.chi = grid
        # self.increments = (1, 1, 10)
        # self.size = (1, 1, 10)
        # self.shape = (1, 1, 100)

    def inv_mass(self, idx, dim):
        return 1.0/self.__eff_mass[dim]

    def inv_mass_xy(self, idx):
        return self.inv_mass(idx, 0)

    def inv_mass_z(self, idx):
        return self.inv_mass(idx, 1)


class EngineTests(object):

    def __init__(self, engine):
        self.TestEngine = engine
        self.well_depth = 2.0
        self.well_width = 1.0

    def setup_3d_harmonic(self):
        from emtpy.grid import PhysicalGrid
        self.potential = Harmonic(1.0, (10, 10, 10), (3.0, 3.0, 3.0))
        self.mat_dist = MockMaterialDistribution(PhysicalGrid((10, 10, 10), (3.0, 3.0, 3.0)))
        self.engine = self.TestEngine(self.mat_dist, self.potential)

    def setup_one_d_well(self):
        from emtpy.grid import PhysicalGrid
        self.potential = OneDWell(2.0, 1.0, (100,), (10.0,))
        self.mat_dist = MockMaterialDistribution(PhysicalGrid((100,), (10.0,)))
        self.engine = self.TestEngine(self.mat_dist, self.potential)

    def setup_three_d_one_d_well(self, middle=5.0):
        from emtpy.grid import PhysicalGrid
        self.potential = ThreeDOneDWell(self.well_width, self.well_depth, (1, 1, 500), (1.0, 1.0, 10.0), middle)
        self.mat_dist = MockMaterialDistribution(PhysicalGrid((1, 1, 500), (1.0, 1.0, 10.0)))
        self.engine = self.TestEngine(self.mat_dist, self.potential)

    @istest
    def test_three_d_one_d_well_energies(self):
        comp_energies = finite_square_box_energies(self.well_width, self.well_depth)
        self.setup_three_d_one_d_well()
        vals, vectors = self.engine.solve()

        for ref_en, en in zip(comp_energies, vals):
            assert_almost_equal(abs(en), ref_en)

    @istest
    def test_periodic_bcs(self):
        """
            Test the energies are right when the well stradles the system boundary
        @return:
        """
        import matplotlib.pyplot as plt
        from numpy import dot
        from itertools import combinations

        self.setup_three_d_one_d_well(0.0)
        plt.plot(self.potential.values[0,0,:])

        plt.show()

        comp_energies = finite_square_box_energies(self.well_width, self.well_depth)

        vals, vectors = self.engine.solve()

        for ref_en, en in zip(comp_energies, vals):
            assert_almost_equal(abs(en), ref_en)


        vecs = [vectors[:, i] for i in range(len(vals))]
        vector_pairs = combinations(vecs, 2)

        for pair in vector_pairs:
            assert_almost_equal(dot(*pair), 0.0)


    @istest
    def test_three_d_one_d_well_wavefuncs(self):
        from numpy import dot
        from itertools import combinations
        self.setup_three_d_one_d_well()
        vals, vectors = self.engine.solve()

        vecs = [vectors[:, i] for i in range(len(vals))]
        vector_pairs = combinations(vecs, 2)

        for pair in vector_pairs:
            assert_almost_equal(dot(*pair), 0.0)

    @nottest
    def plot_three_d_one_d_well_wavefuncs(self):
        import matplotlib.pyplot as plt
        self.setup_three_d_one_d_well()
        vals, vectors = self.engine.solve()

        plt.plot(self.potential.values[0,0,:])
        plt.plot(vectors[:, 0]+vals[0])
        plt.plot(vectors[:, 1]+vals[1])

        plt.show()

    @istest
    def test_harmonic_oscilator_energies(self):
        from emtpy.constants import hbarev
        self.setup_3d_harmonic()
        vals, vectors = self.engine.solve()
        eq_(vals[0], hbarev*1)
        eq_(vals[1], hbarev*2)

    @istest
    def test_one_d_well_energies(self):
        import matplotlib.pyplot as plt
        comp_energies = finite_square_box_energies(2.0E-9, 1000.0)
        self.setup_one_d_well()
        vals, vectors = self.engine.solve()
        #pot = self.potential.values[0,0,:]
        assert_almost_equal(1000.0-abs(vals[0]), infinite_square_box_energy(1, 2.0E-9), 1)
        assert_almost_equal(1000.0-abs(vals[1]), infinite_square_box_energy(2, 2.0E-9), 1)
        plt.plot(self.potential.values[0,0,:])
        plt.plot(vectors[:, 0]+vals[0])
        plt.plot(vectors[:, 1]+vals[1])

        plt.show()



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
        super(ArpackOriginalTests, self).test_periodic_bcs()


    @istest
    def test_harmonic_oscilator_energies(self):
        from emtpy.constants import hbarev
        self.setup_3d_harmonic()
        vals, vectors = self.engine.solve()
        eq_(vals[0], hbarev*1)
        eq_(vals[1], hbarev*2)



    # @istest
    # def test_one_d_well_energies(self):
    #     super(ArpackOriginalTests, self).test_one_d_well_energies()