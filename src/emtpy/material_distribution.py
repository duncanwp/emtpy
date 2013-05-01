__author__ = 'duncan'
from grid import PhysicalGrid


class MaterialDistribution(object):

    def __init__(self, ternary_material, distribution):
        super(MaterialDistribution, self).__init__()
        self.chi = distribution
        self.materialA = ternary_material.materialA
        self.materialB = ternary_material.materialB

    def inv_mass_xy(self, idx):
        return (self.chi(idx)/self.materialA.effective_mass[0]) + ((1.0-self.chi(idx))/self.materialB.effective_mass[0])

    def inv_mass_z(self, idx):
        return (self.chi(idx)/self.materialA.effective_mass[1]) + ((1.0-self.chi(idx))/self.materialB.effective_mass[1])


class Distribution(PhysicalGrid):
    pass
#   subroutine conc(R,X,Ano,sz,ag)
#     implicit none
#     real, dimension(3,maxpos),intent(in) :: R
#     integer, intent(in) :: sz(3),ag(3),Ano
#     complex,intent(out),dimension(0:sz(1)-1,0:sz(2)-1,0:sz(3)-1)::X
#     complex,dimension(0:sz(1)-1,0:sz(2)-1,0:sz(3)-1)::temp
#     integer :: l,m,p,nonn,nolattp(3),n,q
#     real, dimension(3) :: y,s,disp,i
#     real :: const
#     const = ncat/(sqrt((2.0d0*pi)**3)*product(sigma))
#     nonn = 3
#     nolattp = nonn*ag
#     temp = 0.0d0
#     do n = 1, Ano
#        do p = -nolattp(3), nolattp(3)
#           do m = -nolattp(2), nolattp(2)
#              do l = -nolattp(1), nolattp(1)
#                 disp(1) = real(l)*grid(ag,1)
#                 disp(2) = real(m)*grid(ag,2)
#                 disp(3) = real(p)*grid(ag,3)
#                 do q = 1, 3
#                    y(q) = disp(q) + R(q,n)
#                    i(q) = nint(y(q)/grid(ag,q))
#                    s(q) = (disp(q)**2)/(sigma(q)**2)
#                    if (i(q) .lt. 0) then
#                       i(q) = i(q) + sz(q)
#                    else if (i(q) .ge. sz(q)) then
#                       i(q) = i(q) - sz(q)
#                    end if
#                 end do
#                 temp(i(1),i(2),i(3))=temp(i(1),i(2),i(3))+exp(-0.5d0*sum(s))
#              end do
#           end do
#        end do
#     end do
#     X = temp*const
#     return
#   end subroutine conc
#
