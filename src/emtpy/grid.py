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
        self.size = physical_size
        self.volume = product(physical_size)
        self.increments = map(lambda x, y: x/y, physical_size, shape)
        self.values = values
        self.no_elements = product(shape)

    def compatible_grid(self, other):
        return self.shape == other.shape and self.size == other.size

    def coord_array(self, dim):
        from numpy import array
        return array(map(lambda x, y: x*y, self.increments[dim], range(self.size)))

    def coord(self, idx):
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
        from math import floor
        sz = self.shape
        k = floor(float((n-1)/(sz[1]*sz[2])))
        j = floor(float((n-1-(k*sz[1]*sz[2]))/sz[1]))
        i = floor(float(n-1-(k*sz[1]*sz[2])-(j*sz[1])))
        return i, j, k

    def getn(self, i, j, k):
        return j*self.shape[1] + i + k*self.shape[1]*self.shape[2] + 1
