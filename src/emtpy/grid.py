__author__ = 'pard'
from numpy import array


class GridError(RuntimeError):
    pass


class PhysicalGrid(object):

    def __init__(self, shape, physical_size, values=None):
        from utils import product
        #super(PhysicalGrid, self).__init__(shape=shape)
        if len(shape) != len(physical_size):
            raise GridError()
        self.shape = shape
        self.size = map(float, physical_size)
        self.volume = product(physical_size)
        self.increments = map(lambda x, y: x/y, self.size, shape)
        self.values = values
        self.no_elements = product(shape)

    def compatible_grid(self, other):
        """
            Check the compatibility of two grids
        """
        return self.shape == other.shape and self.size == other.size

    def coord_array(self, dim):
        """
            Return an array of physical coordinates along the given dimension
        """
        from numpy import arange
        return arange(self.shape[dim])*self.increments[dim]

    def find_index(self, coord):
        """
            Given a physical coordinate return the corresponding index
        """
        return tuple(map(lambda x, y: int(round(x/y)), coord, self.increments))

    def coord(self, idx):
        """
            Given an index return the corresponding physical coordinate
        """
        from numpy import array
        return array(map(lambda x, y: x*y, self.increments, idx))

    def normalization(self):
        from math import sqrt
        return 1.0/sqrt(self.volume)

    def rect_integral(self):
        return self.values.sum()/(self.normalization()**2)

    def adim(self, dir):
        if (self.direction == 3):
           adim = dir
        elif (self.direction == 1):
           if (dir == 1): adim = 3
           if (dir == 2): adim = 2
           if (dir == 3): adim = 1
        elif (self.direction == 2):
           if (dir == 1): adim = 1
           if (dir == 2): adim = 3
           if (dir == 3): adim = 2
        else:
           raise GridError('ERROR in well direction')
        return adim

    def getijk(self, n):
        """
            Given a number n representing the flattened grid, return the original indices
        """
        from numpy import unravel_index
        return unravel_index(n, self.shape)

    def getn(self, idx):
        """
            Given an index return a single number n representing the corresponding index on the flattened grid
        """
        from utils import ravel_index
        return ravel_index(idx, self.shape)

    def fourier_coord(self, idx):
        """
            Return the fourier index of idx. That is; the index of the corresponding point in the fourier transformed
            grid.

            From the Fortran code:
                n(1) = (mod((l + (sz(1)/2)),sz(1)) - sz(1)/2)
                n(2) = (mod((m + (sz(2)/2)),sz(2)) - sz(2)/2)
                n(3) = (mod((p + (sz(3)/2)),sz(3)) - sz(3)/2)
                k(1) = (2.0d0*pi*n(1))/(sz(1)*grid(ag,1)*1E-9)
                k(2) = (2.0d0*pi*n(2))/(sz(2)*grid(ag,2)*1E-9)
                k(3) = (2.0d0*pi*n(3))/(sz(3)*grid(ag,3)*1E-9)
        """
        from math import pi
        n = map(lambda x, y: ((x + (y/2.0)) % y) - y/2.0, idx, self.shape)
        k = array(map(lambda x, y: 2.0*pi*x/(y*1E-9), n, self.size))
        return k

