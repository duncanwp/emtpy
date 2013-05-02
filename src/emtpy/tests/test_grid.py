__author__ = 'duncan'
import unittest


class TestPhysicalGrid(unittest.TestCase):

    def setUp(self):
        from emtpy.grid import PhysicalGrid
        import numpy as np
        values = np.arange(9*10*11)
        # The 3d grid deliberately has integer physical sizes to test the handling of it
        self.test_3d_grid = PhysicalGrid((9, 10, 11), (5, 5, 5), values.reshape((9, 10, 11)))
        self.test_2d_grid = PhysicalGrid((9, 10), (5.0, 5.0), np.arange(90).reshape((9, 10)))
        self.test_1d_grid = PhysicalGrid((9,), (5.0,), np.arange(9))

    def test_getn_1d(self):
        self.assertEqual(self.test_1d_grid.getn((0,)), 0)
        self.assertEqual(self.test_1d_grid.getn((8,)), self.test_1d_grid.no_elements-1)
        self.assertEqual(self.test_1d_grid.getn(((8 + 1) % 9,)), 0)

    def test_getn_2d(self):
        self.assertEqual(self.test_2d_grid.getn((0, 0)), 0)
        self.assertEqual(self.test_2d_grid.getn((8, 9)), self.test_2d_grid.no_elements-1)
        self.assertEqual(self.test_2d_grid.getn(((8 + 1) % 9, 9)), 9)

    def test_getn_3d(self):
        self.assertEqual(self.test_3d_grid.getn((0, 0 ,0)), 0)
        self.assertEqual(self.test_3d_grid.getn((8, 9, 10)), self.test_3d_grid.no_elements-1)
        self.assertEqual(self.test_3d_grid.getn(((8 + 1) % 9, 9, 10)), 109)

    def test_get_ijk_1d(self):
        self.assertEqual(self.test_1d_grid.getijk(0), (0,))
        self.assertEqual(self.test_1d_grid.getijk(self.test_1d_grid.no_elements-1), (8,))
        self.assertEqual(self.test_1d_grid.getijk(0), ((8 + 1) % 9,))

    def test_get_ijk_2d(self):
        self.assertEqual(self.test_2d_grid.getijk(0), (0, 0))
        self.assertEqual(self.test_2d_grid.getijk(self.test_2d_grid.no_elements-1), (8, 9))
        self.assertEqual(self.test_2d_grid.getijk(9), ((8 + 1) % 9, 9))

    def test_get_ijk_3d(self):
        self.assertEqual(self.test_3d_grid.getijk(0), (0, 0, 0))
        self.assertEqual(self.test_3d_grid.getijk(self.test_3d_grid.no_elements-1), (8, 9, 10))
        self.assertEqual(self.test_3d_grid.getijk(109), ((8 + 1) % 9, 9, 10))

    def test_find_index_1d(self):
        self.assertTupleEqual(self.test_1d_grid.find_index((2.5,)), (5,))
        self.assertTupleEqual(self.test_1d_grid.find_index((0.0,)), (0,))
        self.assertTupleEqual(self.test_1d_grid.find_index((5.0,)), (9,))

    def test_find_index_3d(self):
        self.assertTupleEqual(self.test_3d_grid.find_index((2.5, 2.5, 2.5)), (5, 5, 6))
        self.assertTupleEqual(self.test_3d_grid.find_index((0.0, 0.0, 0.0)), (0, 0, 0))
        self.assertTupleEqual(self.test_3d_grid.find_index((5.0, 5.0, 5.0)), (9, 10, 11))