__author__ = 'pard'


class GridError(RuntimeError):
    pass

class Grid(object):

    def grid(self, no, dir):
        if (dir < 3):
           grid = self.aBN/float(no(dir))
        elif (dir == 3):
           grid = self.cBN/float(no(dir))
        else:
           raise GridError('Invalid dimenion')
        return grid


    def normalization(self, no):
        from math import sqrt
        if (self.nodimensions == 1):
           normalization = 1.0/sqrt(self.grid(no,self.direction))
        else:
           normalization = 1.0/sqrt(self.grid(no,1)*self.grid(no,2)*self.grid(no,3))
        return normalization


    def integral(self, array, ag):
        return sum(array)/(self.normalization(ag)**2)


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

    def getijk(self, n, sz):
        from math import floor
        k = floor(float((n-1)/(sz(1)*sz(2))))
        j = floor(float((n-1-(k*sz(1)*sz(2)))/sz(1)))
        i = floor(float(n-1-(k*sz(1)*sz(2))-(j*sz(1))))
        return (i,j,k)

    def getn(self, n,i,j,k,sz):
        return j*sz(1) + i + k*sz(1)*sz(2) + 1

    def nvalue(self, i, j, k, sz):
        return j*sz(1) + i + k*sz(1)*sz(2) + 1
