__author__ = 'duncan'
from grid import PhysicalGrid


class APotentialEnergy(PhysicalGrid):
    pass


class ThreeDOneDWell(APotentialEnergy):

    def __init__(self, width, depth, shape, size, middle=None):
        import numpy as np
        from utils import heaviside

        # Create a new grid for the potential
        super(ThreeDOneDWell, self).__init__(shape, size)

        if middle is None:
            # Define r_0, the center of the well, to be the middle of the grid
            r_0 = size[2] / 2.0
        else:
            # The middle is defined
            r_0 = middle

        one_d_vector = depth - depth*heaviside((r_0 + (width/2.0) - self.coord_array(2)) % self.size[2]) - \
                       depth*heaviside((r_0 - (width/2.0) - self.coord_array(2)) % self.size[2])
        # one_d_vector = depth - depth*heaviside(((width/2.0) - abs(self.coord_array(2)-r_0) + self.size[2]/2.0) % self.size[2])
        self.values = np.array([[one_d_vector]])
        #self.values = np.tile(one_d_vector,(1, shape[1], shape[2]))


class OneDWell(APotentialEnergy):

    def __init__(self, width, depth, shape, size):
        import numpy as np
        from utils import heaviside

        # Create a new grid for the potential
        super(OneDWell, self).__init__(shape, size)

        # Define r_0, the center of the well, to be the middle of the grid
        r_0 = size[0] / 2.0

        one_d_vector = depth - depth*heaviside(width/2.0-abs(self.coord_array(0)-r_0))
        # self.values = np.array([[one_d_vector]])
        self.values = one_d_vector
        #self.values = np.tile(one_d_vector,(1, shape[1], shape[2]))


class Harmonic(APotentialEnergy):

    def __init__(self, omega, shape, size):
        from constants import eV, me
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


class Offset(APotentialEnergy):

    def __init__(self, mat_dist):
        import numpy as np
        from material import Electron

        X = mat_dist.chi

        V = np.zeros(X.shape)
        if isinstance(mat_dist.carrier, Electron):
           off = mat_dist.cboffset
           zero = BE
           sign = 1.0
        else:
           off = mat_dist.vboffset
           zero = 0.0
           sign = -1.0

        V = sign*off*self.Eg(X) + zero

    def Eg(self, chi):
        return (chi*(AE-BE)) - ((b*chi)*(1.0-chi))


#   real function Eg(chi)
#     real, intent(in) :: chi
#     Eg = (chi*(AE-BE)) - ((b*chi)*(1.0d0-chi))
#   end function Eg
#

#   subroutine offset(X,Con,para,sz,ag)
#     implicit none
#     integer, intent(in) :: sz(3), ag(3)
#     complex,intent(inout),dimension(0:sz(1)-1,0:sz(2)-1,0:sz(3)-1)::X
#     complex,intent(out),dimension(0:sz(1)-1,0:sz(2)-1,0:sz(3)-1)::Con
#     character(2), intent(in) :: para
#     real :: off,zero,sign
#     integer :: l,m,p
#     if (para == '_e') then
#        off = cboffset
#        zero = BE
#        sign = 1.0d0
#     else
#        off = vboffset
#        zero = 0.0d0
#        sign = -1.0d0
#     end if
#     do p = 0, sz(3)-1
#        do m = 0, sz(2)-1
#           do l = 0, sz(1)-1
#              Con(l,m,p) = sign*(off*Eg(real(X(l,m,p)))) + zero
#           end do
#        end do
#     end do
#   end subroutine offset
#

class Deformation(APotentialEnergy):

    def __init__(self, strain, mat_param):
        import numpy as np

        # The choice between different deformation potential parameters based on carrier type should definitely be
        #   done by the material or material_distribution objects. This then reduces this class to one line:
        V = deformation = (strain[1,1,:,:,:]+strain[2,2,:,:,:])*(defBNa2+BND2) + strain[3,3,:,:,:]*(defBNa1+BND1)

        # This obviously won't work for one-d arrays so I wonder if we could make a strain class which takes care of the
        #  extra indices for us. I.e. we could ask for the trace or just certain diagonal components


#   subroutine defpot(U,X,deformation,para,sz,ag)
#     implicit none
#     integer, intent(in) :: sz(3), ag(3)
#     complex,intent(in)::U(3,3,0:sz(1)-1,0:sz(2)-1,0:sz(3)-1)
#     complex,intent(in):: X(0:sz(1)-1,0:sz(2)-1,0:sz(3)-1)
#     complex,intent(out):: deformation(0:sz(1)-1,0:sz(2)-1,0:sz(3)-1)
#     character(2), intent(in) :: para
#     !Assumes homogeneous deformation potential paramaters (GaN)
#     if (para =='_e') then
#        deformation = (U(1,1,:,:,:)+U(2,2,:,:,:))*(defBNa2+BND2) &
#             &+ U(3,3,:,:,:)*(defBNa1+BND1)
#     else if (para =='_h') then
#        deformation = (U(1,1,:,:,:)+U(2,2,:,:,:))*(BND2+BND4) &
#             &+ U(3,3,:,:,:)*(BND1+BND3)
#        !!Assumes DOS hole sees A/B hole deformation potential!!
#     else if (para =='sh') then
#        deformation = (U(1,1,:,:,:)+U(2,2,:,:,:))*(BND2) + U(3,3,:,:,:)*(BND1)
#     else if (para == 'hh' .or. para == 'lh') then
#        !For Heavy and Light holes
#        deformation = (U(1,1,:,:,:)+U(2,2,:,:,:))*(BND2+BND4) &
#             &+ U(3,3,:,:,:)*(BND1+BND3)
#     else
#        deformation = (U(1,1,:,:,:)+U(2,2,:,:,:))*(defBNa2-BND4) &
#             &+ U(3,3,:,:,:)*(defBNa1-BND3)
#     end if
#   end subroutine defpot
#


class SpontaneousPolarization(APotentialEnergy):

    def __init__(self, mat_dist):
        import numpy as np

        X = mat_dist.chi

        V = np.zeros(X.shape)
        for idx, chi in np.ndenumerate(X):
            k = X.fourier_coord(idx)
            V[idx] = self.__fphisp(k, chi)

    def __fphisp(self, kv, chi):
        ksq = sum(kv**2)
        if ksq == 0.0:
            const = 0.0
        else:
            const = (complex(0.0,1.0)*kv(3))/(epsilon0*epsilonr*ksq)
        return const*(chi*(-deltaP))

#   complex function fphisp(kv,chi)
#     real, intent(in), dimension(3) :: kv
#     complex, intent(in) :: chi!,chi2
#     real :: ksq
#     complex :: const
#     ksq = sum(kv**2)
#     if (ksq .eq. 0.0) then
#        const = 0.0d0
#     else
#        const = (cmplx(0.0d0,1.0d0)*kv(3))/(epsilon0*epsilonr*ksq)
#     end if
# !    fphisp = const*((chi*(bPsp-deltaP))-(chi2*bPsp))
#     fphisp = const*(chi*(-deltaP))
#   end function fphisp

#   subroutine spontpot(X,V,sz,ag)
#     implicit none
#     integer, intent(in) :: sz(3), ag(3)
#     complex,intent(in),dimension(0:sz(1)-1,0:sz(2)-1,0:sz(3)-1):: X
#     complex,intent(inout),dimension(0:sz(1)-1,0:sz(2)-1,0:sz(3)-1)::V
#     integer :: l,m,p
#     integer, dimension(3) :: n
#     real, dimension(3) :: k
#     real :: ksq
#     complex :: chi,chi2
#     do p = 0, sz(3)-1
#        do m = 0, sz(2)-1
#           do l = 0, sz(1)-1
#              n(1) = (mod((l + (sz(1)/2)),sz(1)) - sz(1)/2)
#              n(2) = (mod((m + (sz(2)/2)),sz(2)) - sz(2)/2)
#              n(3) = (mod((p + (sz(3)/2)),sz(3)) - sz(3)/2)
#              k(1) = (2.0d0*pi*n(1))/(sz(1)*grid(ag,1)*1E-9)
#              k(2) = (2.0d0*pi*n(2))/(sz(2)*grid(ag,2)*1E-9)
#              k(3) = (2.0d0*pi*n(3))/(sz(3)*grid(ag,3)*1E-9)
#              !Scale coordinates, in nm
#              chi = X(l,m,p)
# !             chi2 = X2(l,m,p)
#              V(l,m,p) = fphisp(k,chi)
#              !Calculates piezoelectric potential and spontaneous polarization
#           end do
#        end do
#     end do
#   end subroutine spontpot
#


class PiezoelectricPolarization(APotentialEnergy):

    def __init__(self, mat_dist):
        import numpy as np

        X = mat_dist.chi

        V = np.zeros(X.shape)
        for idx, chi in np.ndenumerate(X):
            k = X.fourier_coord(idx)
            V[idx] = self.__fphi(k, chi)

    def __fphi(self, kv, chi):
        ksq = sum(kv**2)
        if ksq == 0.0:
            return 0.0
        prefactor = -(chi*kv(3)*eps(1,1)*complex(0.0,1.0))/(epsilonr*epsilon0*(ksq**2))
        fphi = prefactor*(ksq*((2.0*e31) + e33) + gamma*((e31 + 2.0*e15)*((kv(1)**2) + (kv(2)**2)) + (e33*(kv(3)**2))))
        return fphi

#   complex function fphi(kv,chi)
#     real, intent(in), dimension(3) :: kv
#     complex, intent(in) :: chi
#     real :: ksq
#     complex :: prefactor
#     ksq = sum(kv**2)
#     if (ksq .eq. 0.0) then
#        fphi = 0.0d0
#     else
#        prefactor = -(chi*kv(3)*eps(1,1)*cmplx(0.0d0,1.0d0))&
#             &/(epsilonr*epsilon0*(ksq**2))
#        fphi = prefactor*(ksq*((2.0d0*e31) + e33) + gamma*((e31 + &
#             &2.0d0*e15)*((kv(1)**2) + (kv(2)**2)) + (e33*(kv(3)**2))))
#     end if
#   end function fphi

#   subroutine piezopot(X,V,sz,ag)
#     implicit none
#     integer, intent(in) :: sz(3), ag(3)
#     complex,intent(inout),dimension(0:sz(1)-1,0:sz(2)-1,0:sz(3)-1)::V
#     complex,intent(in),dimension(0:sz(1)-1,0:sz(2)-1,0:sz(3)-1):: X
#     integer :: l,m,p
#     integer, dimension(3) :: n
#     real, dimension(3) :: k
#     real :: ksq
#     complex :: chi
#     do p = 0, sz(3)-1
#        do m = 0, sz(2)-1
#           do l = 0, sz(1)-1
#              n(1) = (mod((l + (sz(1)/2)),sz(1)) - sz(1)/2)
#              n(2) = (mod((m + (sz(2)/2)),sz(2)) - sz(2)/2)
#              n(3) = (mod((p + (sz(3)/2)),sz(3)) - sz(3)/2)
#              k(1) = (2.0d0*pi*n(1))/(sz(1)*grid(ag,1)*1E-9)
#              k(2) = (2.0d0*pi*n(2))/(sz(2)*grid(ag,2)*1E-9)
#              k(3) = (2.0d0*pi*n(3))/(sz(3)*grid(ag,3)*1E-9)
#              !Scale coordinates, in nm
#              chi = X(l,m,p)
#              V(l,m,p) = fphi(k,chi)
#              !Calculates piezoelectric potential and spontaneous polarization
#           end do
#        end do
#     end do
#   end subroutine piezopot
#
