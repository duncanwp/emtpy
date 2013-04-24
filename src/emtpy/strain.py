__author__ = 'duncan'

#   subroutine strain(X,TrU,sz,ag)
#     implicit none
#     integer, intent(in) :: sz(3), ag(3)
#     complex,intent(out),dimension(0:sz(1)-1,0:sz(2)-1,0:sz(3)-1)::TrU
#     complex,intent(in),dimension(0:sz(1)-1,0:sz(2)-1,0:sz(3)-1):: X
#     integer :: l,m,p
#     integer, dimension(3) :: n
#     real, dimension(3) :: k
#     real :: ksq
#     do p = 0, sz(3)-1
#        do m = 0, sz(2)-1
#           do l = 0, sz(1)-1
#              n(1) = (mod((l + (sz(1)/2)),sz(1)) - sz(1)/2)
#              n(2) = (mod((m + (sz(2)/2)),sz(2)) - sz(2)/2)
#              n(3) = (mod((p + (sz(3)/2)),sz(3)) - sz(3)/2)
#              k(1) = (2.0d0*pi*n(1))/(sz(1)*grid(ag,1)*1E-9)
#              k(2) = (2.0d0*pi*n(2))/(sz(2)*grid(ag,2)*1E-9)
#              k(3) = (2.0d0*pi*n(3))/(sz(3)*grid(ag,3)*1E-9)
#              ksq = sum(k**2)
#              if (ksq .eq. 0.0) then
#                 TrU(l,m,p) = eps(1,1)*X(l,m,p)
#              else
#                 TrU(l,m,p)=eps(1,1)*X(l,m,p)*(1.0d0+(gamma*((k(3)*k(3))/ksq)))
#              end if
# !             TrU(l,m,p)=eps*(1.0+gamma*((k(3)*k(3))/sum(k**2)))*X(l,m,p)
# !             TrU(l,m,p)=eps*X(l,m,p)*(3.0d0+gamma)
#           end do
#        end do
#     end do
#   end subroutine strain
#


#   subroutine straintensor(X,U,sz,ag)
#     implicit none
#     integer, intent(in) :: sz(3), ag(3)
#     complex,intent(out)::U(3,3,0:sz(1)-1,0:sz(2)-1,0:sz(3)-1)
#     complex,intent(in):: X(0:sz(1)-1,0:sz(2)-1,0:sz(3)-1)
#     integer :: l,m,p,i,j
#     integer, dimension(3) :: n
#     real, dimension(3) :: k
#     real :: ksq
#     do p = 0, sz(3)-1
#        do m = 0, sz(2)-1
#           do l = 0, sz(1)-1
#              do j = 1, 3
#                 do i = 1, 3
#                    n(1) = (mod((l + (sz(1)/2)),sz(1)) - sz(1)/2)
#                    n(2) = (mod((m + (sz(2)/2)),sz(2)) - sz(2)/2)
#                    n(3) = (mod((p + (sz(3)/2)),sz(3)) - sz(3)/2)
#                    k(1) = (2.0d0*pi*n(1))/(sz(1)*grid(ag,1)*1E-9)
#                    k(2) = (2.0d0*pi*n(2))/(sz(2)*grid(ag,2)*1E-9)
#                    k(3) = (2.0d0*pi*n(3))/(sz(3)*grid(ag,3)*1E-9)
#                    ksq = sum(k**2)
#                    if (ksq .eq. 0.0) then
#                       if (i==j) then
#                          if (i==3 .and. j==3) then
#                             U(i,j,l,m,p) = eps(i,j)*(deltaij(i,j)+gamma)*X(l,m,p)
#                          else
#                             U(i,j,l,m,p) = eps(i,j)*X(l,m,p)
#                          end if
#                       else
#                          U(i,j,l,m,p)=0.0d0
#                       end if
# !        U(i,j,l,m,p)=eps*X(l,m,p)*(deltaij(i,j)*(1.0d0+(gamma/3.0d0)))
#                     else
#         U(i,j,l,m,p)=eps(i,j)*X(l,m,p)*(deltaij(i,j)+(gamma*((k(i)*k(j))/ksq)))
#                    end if
#                 end do
#              end do
#           end do
#        end do
#     end do
#   end subroutine straintensor
#
