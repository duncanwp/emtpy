__author__ = 'duncan'
from grid import PhysicalGrid


class APotentialEnergy(PhysicalGrid):
    pass


class OneDWell(APotentialEnergy):

    def __init__(self, width, depth, shape, size):
        import numpy as np
        from utils import heaviside

        # Create a new grid for the potential
        super(OneDWell, self).__init__(shape, size)

        # Define r_0, the center of the well, to be the middle of the grid
        r_0 = size[0] / 2.0

        one_d_vector = 0.0 - depth*heaviside(width/2.0-abs(self.coord_array(0)-r_0))
        self.values = np.tile(one_d_vector,(1, shape[1], shape[2]))


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



#
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
#
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
#
# !!$  complex function Afracp(x,R)
# !!$    real, dimension(3,10000),intent(in) :: R
# !!$    complex :: const, sumation
# !!$    real,intent(in),dimension(3) :: x
# !!$    real,dimension(3) :: s
# !!$    real :: disp
# !!$    integer :: n,m
# !!$!    const = ncat*((2.0d0*pi))/(sqrt((2.0d0*pi)**3)*product(sigma))
# !!$    const = ncat/(sqrt((2.0d0*pi)**3)*product(sigma))
# !!$    sumation = 0.0d0
# !!$    do n = 1, Ano
# !!$       do m = 1, 3
# !!$          if (x(m)-R(m,n) .gt. ((Rsize(m+3)-Rsize(m))/2.0)) then
# !!$             disp = (x(m)-R(m,n)) - (Rsize(m+3)-Rsize(m))
# !!$          else if (x(m)-R(m,n) .le. -((Rsize(m+3)-Rsize(m))/2.0)) then
# !!$             disp = (x(m)-R(m,n)) + (Rsize(m+3)-Rsize(m))
# !!$          else
# !!$             disp = x(m)-R(m,n)
# !!$          end if
# !!$          s(m) = (disp**2)/(sigma(m)**2)
# !!$       end do
# !!$       sumation = sumation + exp(-0.5d0*sum(s))
# !!$    end do
# !!$    Afracp = sumation*const!!!*2.5d0
# !!$    !!!! Factor of 2.5 compensates for the 40% efficiency of detection !!!!!
# !!$  end function Afracp
#
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
