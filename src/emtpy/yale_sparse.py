__author__ = 'duncan'


class YaleSparse(object):

    def __init__(self, ija, sa):
        self.ija = ija
        self.sa = sa

    @classmethod
    def from_dense(cls, dense_matrix, nmax, thresh=0.0):
        """
        Converts a square matrix a(1:n,1:n) with physical dimension np into row-indexed sparse
        storage mode. Only elements of a with magnitude greater than thresh are retained. Output is in
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
        ija[0] = n + 1

        # Index to 1st row off-diagonal element, if any.
        k = n

        # Loop over rows.
        for i in range(n):
            for j in range(n):
                # Loop over columns.
                if(abs(dense_matrix[i, j]) > thresh):
                    if(i != j):
                        # Store off-diagonal elements and their columns.
                        k = k + 1
                        if(k > nmax): raise ValueError("nmax too small in sprsin")
                        sa[k] = dense_matrix[i, j]
                        ija[k] = j
            ija[i+1] = k + 1
            # As each row is completed, store index to next.

        return cls(ija, sa)

    def to_dense(self, shape):
        import numpy as np
        dense = np.zeros(shape)
        #diag = np.diag_indices_from(dense)

        n = shape[0]

        for i in range(n):
            dense[i, i] = self.sa[i]
            for k in range(self.ija[i], (self.ija[i+1])):
                dense[i, self.ija[k]] = self.sa[k]

        return dense

    def av(self, v):
        import numpy as np
        n = len(v)
        y = np.zeros(n, dtype=np.float64)

        for i in range(n):
            y[i] = self.sa[i] * v[i]
            for k in range(self.ija[i], (self.ija[i+1])):
                y[i] += self.sa[k] * v[self.ija[k]]

        # for i in range(n):
        #     for j in range (self.ija[i], (self.ija[i+1]-2)):
        #         y[self.ija[j]] = y[self.ija[j]] + self.sa[j]*v[i] #Use for symetric storage

        return y


class YaleSparseSymmetric(object):

    def __init__(self, ija, sa):
        self.ija = ija
        self.sa = sa

    @classmethod
    def from_dense(cls, dense_matrix, nmax, thresh=0.0):
        """
        Converts a square matrix a(1:n,1:n) with physical dimension np into row-indexed sparse
        storage mode. Only elements of a with magnitude greater than thresh are retained. Output is in
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
        ija[0] = n + 1

        # Index to 1st row off-diagonal element, if any.
        k = n

        # Loop over rows.
        for i in range(n):
            for j in range(i+1,n):
                # Loop over columns.
                if(abs(dense_matrix[i, j]) > thresh):
                    # Store off-diagonal elements and their columns.
                    k = k + 1
                    if(k > nmax): raise ValueError("nmax too small in sprsin")
                    sa[k] = dense_matrix[i, j]
                    ija[k] = j
            ija[i+1] = k + 1
            # As each row is completed, store index to next.

        return cls(ija, sa)

    def to_dense(self, shape):
        import numpy as np
        dense = np.zeros(shape)
        #diag = np.diag_indices_from(dense)

        n = shape[0]

        for i in range(n):
            dense[i, i] = self.sa[i]
            for k in range(self.ija[i], (self.ija[i+1])):
                j = self.ija[k]
                dense[i, j] = self.sa[k]
                dense[j, i] = dense[i, j]

        return dense

    def av(self, v):
        import numpy as np
        n = len(v)
        y = np.zeros(n, dtype=np.float64)

        if self.ija[0] != n + 1:
            raise ValueError("Mismatched vector and matrix")

        for i in range(n):
            y[i] = self.sa[i] * v[i]
            for k in range(self.ija[i], (self.ija[i+1])):
                y[i] += self.sa[k] * v[self.ija[k]]

        for i in range(n):
            for k in range(self.ija[i], (self.ija[i+1])):
                y[self.ija[k]] += self.sa[k]*v[i] #Use for symetric storage

        return y
