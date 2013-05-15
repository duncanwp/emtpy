__author__ = 'duncan'

from nose.tools import istest, eq_
from semi_analytic_solution import AnalyticSolution


class TestAnalyticSolution(object):

    def __init__(self):
        from emtpy.constants import eV
        self.AS = AnalyticSolution(2.0E-9, 10.0*eV)

    @istest
    def plot_symmetric_function(self):
        import matplotlib.pyplot as pyp
        import numpy as np
        import math
        x = np.linspace(0.0, math.sqrt(self.AS.u_0_sq), 100)

        pyp.plot(x, self.AS.symmetric_solution_vec(x), label='diff')
        pyp.plot(x, self.AS.a(x), label='a')
        pyp.plot(x, self.AS.b_symm(x), label='b')
        pyp.ylim(0.0, math.sqrt(self.AS.u_0_sq))
        pyp.xlim(0.0, math.sqrt(self.AS.u_0_sq))
        pyp.legend()
        pyp.show()

    @istest
    def plot_anti_symmetric_function(self):
        import matplotlib.pyplot as pyp
        import numpy as np
        import math

        x = np.linspace(0, math.sqrt(self.AS.u_0_sq), 100)

        pyp.plot(x, self.AS.anti_symmetric_solution_vec(x), label='diff')
        pyp.plot(x, self.AS.a(x), label='a')
        pyp.plot(x, self.AS.b_anti_symm(x), label='b')
        pyp.ylim(0.0, math.sqrt(self.AS.u_0_sq))
        pyp.xlim(0.0, math.sqrt(self.AS.u_0_sq))
        pyp.legend()
        pyp.show()

    @istest
    def plot_both(self):
        import matplotlib.pyplot as pyp
        import numpy as np
        import math

        x = np.linspace(0, math.sqrt(self.AS.u_0_sq), 100)

        pyp.plot(x, self.AS.a(x), label='a')
        pyp.plot(x, self.AS.b_anti_symm(x), label='b_anti')
        pyp.plot(x, self.AS.b_symm(x), label='b_symm')
        pyp.ylim(0.0, math.sqrt(self.AS.u_0_sq))
        pyp.xlim(0.0, math.sqrt(self.AS.u_0_sq))
        pyp.legend()
        pyp.show()

    @istest
    def test_symmetric_function(self):
        import matplotlib.pyplot as pyp
        import numpy as np
        import math

        eq_(self.AS.finite_square_box_energies(), 0.5)

    @istest
    def test_anti_symmetric_function(self):
        import matplotlib.pyplot as pyp
        import numpy as np
        import math
