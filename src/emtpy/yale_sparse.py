__author__ = 'duncan'

class YaleSparse(object):

    def __init__(self, ija, sa):
        self.ija = ija
        self.sa = sa

    @classmethod
    def from_dense(cls, dense_matrix, nmax, thresh=0.0):
        """
        Converts a square matrix a(1:n,1:n) with physical dimension np into row-indexed sparse
        storage mode. Only elements of a with magnitude ≥thresh are retained. Output is in
        two linear arrays with physical dimension nmax (an input parameter): sa(1:) contains
        array values, indexed by ija(1:). The logical sizes of sa and ija on output are both
        ija(ija(1)-1)-1 (see text).
        """
        import numpy as np
        sa = np.zeros(nmax, dtype=np.float64)
        ija = np.zeros(nmax, dtype=np.int8)

        n = dense_matrix.shape[0]
        # Store diagonal elements.

        sa[:n] = dense_matrix.diagonal()
        ija[1]=n+2
        # Index to 1st row off-diagonal element, if any.
        k=n+1
        # Loop over rows.
        for i in range(1,n):
            for j in range(1,n):
                # Loop over columns.
                if(abs(dense_matrix[i,j]).ge.thresh):
                    if(i != j):
                        # Store off-diagonal elements and their columns.
                        k=k+1
                        if(k > nmax): raise ValueError("nmax too small in sprsin")
                        sa[k]=dense_matrix[i,j]
                        ija[k]=j
            ija[i+1]=k+1
            # As each row is completed, store index to next.

        return cls(ija, sa)

    def av(self, v):
        import numpy as np
        n = len(v)
        y = np.zeros(n, dtype=np.float64)

        for i in range(n):
           y[i] = self.sa[i]*v[i]
           for j in range (self.ija[i], (self.ija[i+1]-2)):
              y[i] += self.sa[j]*v[self.ija[j]]

        for i in range(n):
           for j in range (self.ija[i], (self.ija[i+1]-2)):
              y[self.ija[j]] = y[self.ija[j]] + self.sa[j]*v[i] #Use for symetric storage

        return y
