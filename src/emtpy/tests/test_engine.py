__author__ = 'pard'
from nose.tools import istest, eq_
from emtpy.potential_energy import APotentialEnergy


class harmonic(APotentialEnergy):

    def __init__(self, omega, shape, size):
      #   implicit none
      #   integer, intent(in) :: sz(3), ag(3)
      #   complex,intent(inout),dimension(0:sz(1)-1,0:sz(2)-1,0:sz(3)-1)::V
      #   real, intent(in) :: omega
      #   real, intent(in), dimension(3) :: r
      #   integer :: l,m,p
      #   real :: const = ((0.5d0*me*1E-18)/eV)
      #   real, dimension(3) :: x
      #   do l = 0, sz(1)-1
      #      do m = 0, sz(2)-1
      #         do p = 0, sz(3)-1
      #            x(1) = (l*grid(ag,1))+Rsize(1)
      #            x(2) = (m*grid(ag,1))+Rsize(2)
      #            x(3) = (p*grid(ag,3))+Rsize(3)
      #            V(l,m,p) = const*(omega**2)*sum((x-r)**2)
      #         end do
      #      end do
      #   end do
      # end subroutine harmonic
        from emtpy.constants import eV, me
        from numpy import fromfunction

        # Create a new grid for the potential
        super(harmonic, self).__init__(shape, size)

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
            # pos = []
            # pos[0] = (i*shape[0]/size[0] - r_0[0])**2
            # pos[1] = (j*shape[1]/size[1] - r_0[1])**2
            # pos[2] = (k*shape[2]/size[2] - r_0[2])**2
            return const*(omega**2)*sum(pos)

        self.values = fromfunction(_harmonic_oscillator, shape)


def test_carrier():
    pass

class EngineTests(object):

    def __init__(self):
        self.potential = harmonic(1.0, (100, 100, 100), (10.0, 10.0, 10.0))
        self.carrier = test_carrier()

    @istest
    def test_harmonic_oscilator_energies(self):
        self.engine.solve(self.potential, self.carrier)


class ArpackTests(EngineTests):

    def __init__(self):
        from emtpy.engine import Arpack
        super(ArpackTests, self).__init__()
        self.engine = Arpack()