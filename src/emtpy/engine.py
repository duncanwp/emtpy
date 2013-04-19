__author__ = 'pard'

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
#   subroutine Apos(AllR,R,rrun,Ano)
#     implicit none
#     real, dimension(3,maxpos),intent(in) :: AllR
#     real, dimension(3,maxpos),intent(out) :: R
#     real, dimension(6), intent(in) :: rrun
#     integer, intent(out) :: Ano
#     integer :: n,ifail, p
#     real :: t
#     real, dimension(3) :: x
#     p = 0
#     do n = 1, totA
#        if (AllR(1,n) .ge. rrun(1) .and. AllR(1,n) .lt. rrun(4)) then
#           if (AllR(2,n) .ge. rrun(2) .and. AllR(2,n) .lt. rrun(5)) then
#              if (AllR(3,n) .ge. rrun(3) .and. AllR(3,n) .lt. rrun(6)) then
#                    p = p + 1
#                    R(:,p) = AllR(:,n)-rrun(1:3)
#             end if
#           end if
#        end if
#     end do
#     Ano = p
#   end subroutine Apos
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
#   subroutine outputArr(Arr,unit,dimn,slice,sz,ag,val)
#     implicit none
#     integer, intent(in) :: sz(3), ag(3)
#     complex,intent(in),dimension(0:sz(1)-1,0:sz(2)-1,0:sz(3)-1) :: Arr
#     integer, intent(in) :: unit,dimn
#     integer, intent(in) :: slice
#     real, intent(in), optional :: val
#     !    dim = 1 > xslice, 2 > yslice, 3 > zslice
#     real,dimension(:,:), allocatable :: surfacepos
#     integer :: l,m,p, j, surfacearea, count
#     real :: x1,x2,x3
# 10  format(4G24.15E3)
# 11  format(3G24.15E3)
#     surfacearea = 2*((sz(1)*sz(2))+sz(1)*sz(3)+sz(2)*sz(3))
#     allocate(surfacepos(3,surfacearea))
#
#     if (present(val) .and. dimn .ne. 4) &
#          &write(6,*) 'ERROR in calling array output'
#     if (dimn == 1) write(unit,*) &
#          &'#  y/nm           z/nm           Value          Average in x'
#     if (dimn == 2) write(unit,*) &
#          &'#  x/nm           z/nm           Value          Average in y'
#     if (dimn == 3) write(unit,*) &
#          &'#  x/nm           y/nm           Value          Average in z'
#     if (dimn == 4) write(unit,*) '#  x/nm           y/nm           z/nm'
#     if (dimn == 1) then
#        l = slice
#        do p = 0, sz(3)-1
#           do m = 0, sz(2)-1
#              x2 = m*grid(ag,1)
#              x3 = p*grid(ag,3)
#              write(unit,10) x2, x3, real(Arr(l,m,p)), &
#                   &sum(real(Arr(:,m,p)))/sz(1)
#           end do
#        end do
#     else if (dimn == 2) then
#        m = slice
#        do p = 0, sz(3)-1
#           do l = 0, sz(1)-1
#              x1 = l*grid(ag,1)
#              x3 = p*grid(ag,3)
#              write(unit,10) x1, x3, real(Arr(l,m,p)), &
#                   &sum(real(Arr(l,:,p)))/sz(2)
#           end do
#        end do
#     else if (dimn == 3) then
#        p = slice
#        do m = 0, sz(2)-1
#           do l = 0, sz(1)-1
#              x1 = l*grid(ag,1)
#              x2 = m*grid(ag,1)
#              write(unit,10) x1, x2, real(Arr(l,m,p)), &
#                   &sum(real(Arr(l,m,:)))/sz(3)
#           end do
#        end do
#     else
#        call isosurface(real(Arr),val,surfacepos,surfacearea,count,sz,ag)
# !       write(6,*) 'No of isosurface points for psi^2 =',val,':',count
# !       write(6,*) 'Ratio of Total SA to psi SA:',real(count)/&
# !            &real(surfacearea)
#        do j = 1, count
#           write(unit,11) surfacepos(:,j)
#        end do
#
#     end if
#   end subroutine outputArr
#
#   real function invmassxy(chi,para)
#     implicit none
#     real, intent(in) :: chi
#     character(2), intent(in) :: para
#     if (para =='_e') then
#        invmassxy = (chi/mANxy) + ((1.0d0-chi)/mBNxy)
#     else if (para == 'hh') then
#        invmassxy = (chi/mhhANxy) + ((1.0d0-chi)/mhhBNxy)
#     else if (para == 'lh') then
#        invmassxy = (chi/mlhANxy) + ((1.0d0-chi)/mlhBNxy)
#     else if (para == 'sh') then
#        invmassxy = (chi/mshANxy) + ((1.0d0-chi)/mshBNxy)
#     else if (para == '_h') then
#        invmassxy = (chi/mhANDOS) + ((1.0d0-chi)/mhBNDOS)
#     else
#        write(6,*) "Missing carrier parameter"
#        call abort
#     end if
#   end function invmassxy
#
#   real function invmassz(chi,para)
#     implicit none
#     real, intent(in) :: chi
#     character(2), intent(in) :: para
#     if (para =='_e') then
#        invmassz = (chi/mANz) + ((1.0d0-chi)/mBNz)
#     else if (para == 'hh') then
#        invmassz = (chi/mhhANz) + ((1.0d0-chi)/mhhBNz)
#     else if (para == 'lh') then
#        invmassz = (chi/mlhANz) + ((1.0d0-chi)/mlhBNz)
#     else if (para == 'sh') then
#        invmassz = (chi/mshANz) + ((1.0d0-chi)/mshBNz)
#     else if (para == '_h') then
#        invmassz = (chi/mhANDOS) + ((1.0d0-chi)/mhBNDOS)
#     else
#        write(6,*) "Missing carrier parameter"
#        call abort
#     end if
#   end function invmassz
#
#   real function an3(i,j,k,X,para,sz,ag)
#     implicit none
#     integer, intent(in) :: i, j, k
#     integer, intent(in) :: sz(3),ag(3)
#     complex, intent(in), dimension(0:sz(1)-1,0:sz(2)-1,0:sz(3)-1) :: X
#     character(2), intent(in) :: para
#     real :: chia,chib
#     chia = real(X(i,j,k))
#     chib = real(X(mod(i-1+sz(1),sz(1)),j,k))
#     an3 = (-units/(4.0d0*(grid(ag,1)**2)))*(invmassxy(chia,para)+invmassxy(chib,para))
#   end function an3
#
#   real function bn3(i,j,k,X,para,sz,ag)
#     implicit none
#     integer, intent(in) :: i, j, k
#     integer, intent(in) :: sz(3),ag(3)
#     complex, intent(in), dimension(0:sz(1)-1,0:sz(2)-1,0:sz(3)-1) :: X
#     character(2), intent(in) :: para
#     real :: chia,chib
#     chia = real(X(i,j,k))
#     chib = real(X(i,mod(j-1+sz(2),sz(2)),k))
#     bn3 = (-units/(4.0d0*(grid(ag,1)**2)))*(invmassxy(chia,para)+invmassxy(chib,para))
#   end function bn3
#
#   real function cn3(i,j,k,X,para,sz,ag)
#     implicit none
#     integer, intent(in) :: i, j, k
#     integer, intent(in) :: sz(3),ag(3)
#     complex, intent(in), dimension(0:sz(1)-1,0:sz(2)-1,0:sz(3)-1) :: X
#     character(2), intent(in) :: para
#     real :: chia,chib
#     chia = real(X(i,j,k))
#     chib = real(X(i,j,mod(k-1+sz(3),sz(3))))
#     cn3 = (-units/(4.0d0*(grid(ag,3)**2)))*(invmassz(chia,para) + invmassz(chib,para))
#   end function cn3
#
#   real function dn4(i,j,k,X,V,para,sz,ag)
#     implicit none
#     integer, intent(in) :: i, j, k
#     integer, intent(in) :: sz(3),ag(3)
#     complex,intent(in),dimension(0:sz(1)-1,0:sz(2)-1,0:sz(3)-1)::X,V
#     character(2), intent(in) :: para
#     real :: chi,chia,chib,alphazero,betazero,gammazero
#     chi = real(X(i,j,k))
#     chia = real(X(mod(i+1,sz(1)),j,k))
#     chib = real(X(mod(i-1+sz(1),sz(1)),j,k))
#     alphazero = 2.0d0*invmassxy(chi,para) + invmassxy(chia,para) + &
#          &invmassxy(chib,para)
#
#     chia = real(X(i,mod(j+1,sz(2)),k))
#     chib = real(X(i,mod(j-1+sz(2),sz(2)),k))
#     betazero = 2.0d0*invmassxy(chi,para) + invmassxy(chia,para) + &
#          &invmassxy(chib,para)
#
#     chia = real(X(i,j,mod(k+1,sz(3))))
#     chib = real(X(i,j,mod(k-1+sz(3),sz(3))))
#     gammazero = 2.0d0*invmassz(chi,para) + invmassz(chia,para) + &
#          &invmassz(chib,para)
#
#     dn4 = (units*0.25d0*((alphazero+betazero)/(grid(ag,1)**2) + &
#          &(gammazero)/(grid(ag,3)**2))) + real(V(i,j,k))
#   end function dn4
#
#   subroutine wavefncmplx(nev,ncv,tol,maxitr,evec,eval,confined,para,n,ag,calc)
#     implicit none
#     character(2),intent(in) :: para
#     integer, intent(in) :: nev, maxitr,ncv,n,ag(3),calc
#     real, intent(in) :: tol
#     complex, dimension(n,nev), intent(out) :: evec
#     real, dimension(nev,2), intent(out) :: eval
#     integer, intent(out) :: confined
#     real, allocatable, dimension(:,:) :: rd
#     complex, allocatable,dimension(:,:) :: v
#     complex, allocatable, dimension(:) :: workev
#     complex, allocatable, dimension(:) :: workl,workd,resid,ax,d
#     real, allocatable, dimension(:) :: rwork
#     logical :: slct(ncv)
#     integer :: iparam(11), ipntr(14),selv(2)
#     character :: bmat*1, which*2
#     integer :: ido, lworkl, info, ierr
#     integer :: i,j, ishfts, mode1, nconv,funit
#     integer :: logfil, ndigit, mgetv0, mnaupd, mnaup2, mnaitr
#     integer :: mneigt, mnapps, mngets, mneupd
#     integer :: paranum
#     common /debug/ logfil, ndigit, mgetv0, mnaupd, &
#          & mnaup2, mnaitr, mneigt, mnapps, mngets, mneupd
#
#     logical :: rvec
#     complex :: shift
#
#     !   BLAS & LAPACK routines used
#     Real dznrm2, dlapy2
#     external dznrm2, dlapy2, zaxpy
#     intrinsic abs
#
#     selv = (/ 1, 3/)
#
#     allocate(v(n,ncv))
#     allocate(d(ncv))
#     allocate(rd(ncv,3))
#     allocate(workl(3*ncv*ncv+5*ncv))
#     allocate(workev(ncv*3))
#     allocate(workd(3*n))
#     allocate(resid(n))
#     allocate(rwork(ncv))
#     allocate(ax(n))
#
#     resid = 0.0d0
#     info = 0
#
#     confined = 0
#
#     iparam = 0
#     ipntr = 0
#     v = 0.0d0
#     workl = 0.0d0
#     ax = 0.0d0
#     nconv = 0
#     d = 0.0d0
#
#     if (para == 'hh') paranum = 1
#     if (para == 'lh') paranum = 2
#     if (para == 'sh') paranum = 3
#     if (para == '_h') paranum = 0
#     if (para == '_e') paranum = -1
#
#     ndigit = -3 !specifies the number of decimal digits
#     logfil = 17  !Unit for output
#     mngets = 0
#     mnaitr = 0  !Information about each restart
#     mnapps = 0  !Information about deflation...
#     mnaupd = 1  !Information about Ritz values
#     mnaup2 = 0  !Further information about Ritz values
#     mneigt = 0  !Schir matrix and Hessenberg matrix
#     mneupd = 1  !Print final set of converged Ritz values
#     !Sets the level of output, see debug.doc for details
#
#     bmat  = 'I'
#     !computation mode, I is standard eigenvalue problem
#     which = 'SR'
#     !Chooses the Smallest Real eigenvalues
#     !Chooses which eigenvalues to compute
#     !other options SM, LA, SA, LI, SI
#     shift = 0.0d0
#
#     lworkl = 3*ncv*ncv+5*ncv
#     ido = 0
#     ishfts = 1
#     mode1 = 1
#     iparam(1) = ishfts
#     iparam(3) = maxitr
#     iparam(7) = mode1
#
#     if (para == '_e') then
#        write(17,*) '---------Electron wave function---------'
#     else
#        write(17,*) '---------',para,' wave function---------'
#     end if
#
# 10  continue
#     call znaupd ( ido, bmat, n, which, nev, tol, resid, &
#          & ncv, v, n, iparam, ipntr, workd, workl, &
#          & lworkl,rwork, info )
#     if (ido .eq. -1 .or. ido .eq. 1) then
#        call avc (n, workd(ipntr(1)), workd(ipntr(2)))
#        go to 10
#     end if
#     if ( info .lt. 0 ) then
#        print *, ' '
#        print *, ' Error with _saupd, info = ', info
#        print *, ' Check documentation in _saupd '
#        print *, ' '
#     else
#        rvec = .true.
#        call zneupd ( rvec, 'A', slct, d, v, n, shift, workev, &
#             &         bmat, n, which, nev, tol, resid, ncv, v, n, &
#             &         iparam, ipntr, workd, workl, lworkl, rwork, ierr )
#        if ( ierr .ne. 0) then
#           print *, ' '
#           print *, ' Error with _seupd, info = ', ierr
#           print *, ' Check the documentation of _seupd. '
#           print *, ' '
#        else
#           nconv =  iparam(5)
#           do 20 j=1, nconv
#              call avc(n, v(1,j), ax)
#              call zaxpy(n, -d(j), v(1,j), 1, ax, 1)
#              rd(j,1) = real(d(j))
#              rd(j,2) = aimag(d(j))
#              rd(j,3) = dznrm2(n,ax,1)
#              rd(j,3) = rd(j,3) / dlapy2(rd(j,1),rd(j,2))
# 20           continue
#              call dmout(17, nconv, 3, rd, ncv, -6,&
#                   &            'Ritz values and relative residuals')
#        end if
#        if ( info .eq. 1) then
#           print *, ' '
#           print *, ' Maximum number of iterations reached.'
#           print *, ' '
#        else if ( info .eq. 3) then
#           print *, ' '
#           print *, ' No shifts could be applied during implicit',&
#                &                ' Arnoldi update, try increasing NCV.'
#           print *, ' '
#        end if
#
#        write(17,*)
#        write(17,*)  ' '
#        write(17,*)   ' ARPACK Routine '
#        write(17,*)   ' ============== '
#        write(17,*) ' '
#        write(17,*) ' Size of the matrix is ', n
#        write(17,*) ' The number of Ritz values requested is ', nev
#        write(17,*) ' The number of Arnoldi vectors generated (NCV) is ', ncv
#        write(17,*) ' What portion of the spectrum: ', which
#        write(17,*) ' The number of converged Ritz values is ',nconv
#        write(17,*) ' The number of Implicit Arnoldi update iterations &
#             &taken is ', iparam(3)
#        write(17,*) ' The number of OP*x is ', iparam(9)
#        write(17,*) ' The convergence criterion is ', tol
#        write(17,*) ' '
#     end if
#
#     write(6,*) 'Converged eigenvalues:'
#     write(6,*) iparam(5),nconv
#     write(6,*) '(',nev,')'
#
#     write(6,*) slct
#
#     do i = 1, nconv
#        write(6,*) 'Max values of eigenvector',i
#        write(6,*) maxval(real(v(:,i))),maxval(aimag(v(:,i)))
#     end do
#
#     if (para .ne. '_e') then
#        rd(:,1) = -rd(:,1)
#        do i = 1, iparam(5)
#           if (rd(i,1) .gt. minVcon(calc,paranum+1)) then
#              confined = confined + 1
#           end if
#        end do
#        funit = 47
#        write(6,*) 'Confined holes: ',confined
#        write(17,*) 'Confined holes: ',confined
#     else
#        do i = 1, iparam(5)
#           if (rd(i,1) .lt. maxCcon) then
#              confined = confined + 1
#           end if
#        end do
#        funit = 46
#        write(17,*) 'Confined electrons: ',confined
#        write(6,*) 'Confined electrons: ',confined
#     end if
#
#     eval = rd(:nev,selv)
#     evec = v(:,:nev)*normalization(ag)
#
#   end subroutine wavefncmplx
#
#   subroutine wavefn(nev,ncv,tol,maxitr,evec,eval,confined,para,n,ag,calc)
#     implicit none
#     character(2),intent(in) :: para
#     integer, intent(in) :: nev, maxitr,ncv,n,ag(3),calc
#     real, intent(in) :: tol
#     complex, dimension(n,nev), intent(out) :: evec
#     real, dimension(nev,2), intent(out) :: eval
#     integer, intent(out) :: confined
#     real, allocatable,dimension(:,:) :: v, d
#     real, allocatable, dimension(:) :: workl,workd,resid,ax
#     logical :: select(ncv)
#     integer :: iparam(11), ipntr(11)
#     character :: bmat*1, which*2
#     integer :: ido, lworkl, info, ierr
#     integer :: i,j, ishfts, mode1, nconv,paranum
#     integer :: logfil, ndigit, mgetv0, msaupd, msaup2, msaitr
#     integer :: mseigt, msapps, msgets, mseupd
#     integer,dimension(2) :: selv
#
#     common /debug/ logfil, ndigit, mgetv0, msaupd, &
#          & msaup2, msaitr, mseigt, msapps, msgets, mseupd
#     logical :: rvec
#     Real :: shift, four
#
#     !   BLAS & LAPACK routines used
#     Real dnrm2
#     external dnrm2, daxpy
#     intrinsic abs
#
#     selv = (/ 1, 2 /)
#
#     allocate(v(n,ncv))
#     allocate(d(ncv,2))
#     allocate(workl(ncv*(ncv+8)))
#     allocate(workd(3*n))
#     allocate(resid(n))
#     allocate(ax(n))
#
#     if (para == 'hh') paranum = 1
#     if (para == 'lh') paranum = 2
#     if (para == 'sh') paranum = 3
#     if (para == '_h') paranum = 0
#     if (para == '_e') paranum = -1
#
#     resid = 0.0d0
#     info = 0
#
#     confined = 0
#
#     iparam = 0
#     ipntr = 0
#     v = 0.0d0
#     workl = 0.0d0
#     ax = 0.0d0
#     nconv = 0
#     d = 0.0d0
#
#     ndigit = -3 !specifies the number of decimal digits
#     logfil = 17  !Unit for output
#     msgets = 0
#     msaitr = 0  !Information about each restart
#     msapps = 0  !Information about deflation...
#     msaupd = 1  !Information about Ritz values
#     msaup2 = 0  !Further information about Ritz values
#     mseigt = 0  !Schir matrix and Hessenberg matrix
#     mseupd = 1  !Print final set of converged Ritz values
#     !Sets the level of output, see debug.doc for details
#
#     bmat  = 'I'
#     !computation mode, I is standard eigenvalue problem
#     which = 'SA'
#     !Chooses which eigenvalues to compute
#     !other options SM, LA, SA, LI, SI
#     shift = 0.0d0
#
#     lworkl = ncv*(ncv+8)
#     ido = 0
#     ishfts = 1
#     mode1 = 1
#     iparam(1) = ishfts
#     iparam(3) = maxitr
#     iparam(7) = mode1
#
#     if (para == '_e') then
#        write(17,*) '---------Electron wave function---------'
#     else
#        write(17,*) '---------',para,' wave function---------'
#     end if
#
# 10  continue
#     call dsaupd ( ido, bmat, n, which, nev, tol, resid, &
#          & ncv, v, n, iparam, ipntr, workd, workl, &
#          & lworkl, info )
#     if (ido .eq. -1 .or. ido .eq. 1) then
#        call av (n, workd(ipntr(1)), workd(ipntr(2)))
#        go to 10
#     end if
#     if ( info .lt. 0 ) then
#        print *, ' '
#        print *, ' Error with _saupd, info = ', info
#        print *, ' Check documentation in _saupd '
#        print *, ' '
#     else
#        rvec = .true.
#        call dseupd ( rvec, 'All', select, d, v, n, shift, &
#             &         bmat, n, which, nev, tol, resid, ncv, v, n, &
#             &         iparam, ipntr, workd, workl, lworkl, ierr )
#        if ( ierr .ne. 0) then
#           print *, ' '
#           print *, ' Error with _seupd, info = ', ierr
#           print *, ' Check the documentation of _seupd. '
#           print *, ' '
#        else
#           nconv =  iparam(5)
#           do 20 j=1, nconv
#              call av(n, v(1,j), ax)
#              call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
#              d(j,2) = dnrm2(n, ax, 1)
#              d(j,2) = d(j,2) / abs(d(j,1))
# 20           continue
#              call dmout(17, nconv, 2, d, ncv, -6,&
#                   &            'Ritz values and relative residuals')
#        end if
#        if ( info .eq. 1) then
#           print *, ' '
#           print *, ' Maximum number of iterations reached.'
#           print *, ' '
#        else if ( info .eq. 3) then
#           print *, ' '
#           print *, ' No shifts could be applied during implicit',&
#                &                ' Arnoldi update, try increasing NCV.'
#           print *, ' '
#        end if
#
#        write(17,*)
#        write(17,*)  ' '
#        write(17,*)   ' ARPACK Routine '
#        write(17,*)   ' ============== '
#        write(17,*) ' '
#        write(17,*) ' Size of the matrix is ', n
#        write(17,*) ' The number of Ritz values requested is ', nev
#        write(17,*) ' The number of Arnoldi vectors generated (NCV) is ', ncv
#        write(17,*) ' What portion of the spectrum: ', which
#        write(17,*) ' The number of converged Ritz values is ',nconv
#        write(17,*) ' The number of Implicit Arnoldi update iterations &
#             &taken is ', iparam(3)
#        write(17,*) ' The number of OP*x is ', iparam(9)
#        write(17,*) ' The convergence criterion is ', tol
#        write(17,*) ' '
#     end if
#
#     if (para .ne. '_e') then
#        d(:,1) = -d(:,1)
#        do i = 1, iparam(5)
#           if (d(i,1) .gt. minVcon(calc,paranum+1)) then
#              confined = confined + 1
#           end if
#        end do
#        write(6,*) 'Confined holes: ',confined
#        write(17,*) 'Confined holes: ',confined
#     else
#        do i = 1, iparam(5)
#           if (d(i,1) .lt. maxCcon) then
#              confined = confined + 1
#           end if
#        end do
#        write(17,*) 'Confined electrons: ',confined
#        write(6,*) 'Confined electrons: ',confined
#     end if
#
#     eval = d(:nev,selv)
#     evec = cmplx(v(:,:nev),0.0d0)*normalization(ag)
#   end subroutine wavefn


 #  subroutine sprsymqckprdc(V,X,para,sz,ag)
 #    implicit none
 #    integer, intent(in) :: sz(3),ag(3)
 #    complex,intent(in),dimension(0:sz(1)-1,0:sz(2)-1,0:sz(3)-1)::V,X
 #    character(2), intent(in) :: para
 #    integer :: i, j, k, l, m, n, p, q,nm, elements
 #    real :: H,sum
 #
 #    k = product(sz) + 1
 #    ija(1) = product(sz) + 2
 #
 #    do q = 0, sz(3)-1
 #       do m = 0, sz(2)-1
 #          do l = 0, sz(1)-1
 #             call getn(j,l,m,q,sz)
 #             n = 0
 #             sa(j) = dn4(l,m,q,X,V,para,sz,ag)
 #             do p = 1, 3
 #                H = 0.0d0
 #                if (p == 1) then
 #                   if (l == sz(1)-1) then
 #                      H = an3(0,m,q,X,para,sz,ag)
 #                      call getn(n,0,m,q,sz)
 #                   else
 #                      H = an3(l+1,m,q,X,para,sz,ag)
 #                      call getn(n,l+1,m,q,sz)
 #                   end if
 #                else if (p == 2) then
 #                   if (m == sz(2)-1) then
 #                      H = bn3(l,0,q,X,para,sz,ag)
 #                      call getn(n,l,0,q,sz)
 #                   else
 #                      H = bn3(l,m+1,q,X,para,sz,ag)
 #                      call getn(n,l,m+1,q,sz)
 #                   end if
 #                else if (p == 3) then
 #                   if (q == sz(3)-1) then
 #                      H = cn3(l,m,0,X,para,sz,ag)
 #                      call getn(n,l,m,0,sz)
 #                   else
 #                      H = cn3(l,m,q+1,X,para,sz,ag)
 #                      call getn(n,l,m,q+1,sz)
 #                   end if
 #                end if
 #                if (H .ne. 0.0) then
 #                   k = k + 1
 #                   sa(k) = H
 #                   ija(k) = n
 #                end if
 #             end do
 #             ija(j+1) = k + 1
 #          end do
 #       end do
 #    end do
 #    elements = k
 #
 #    write(17,*) "Elements used: ", elements
 #    write(17,*) "Space allocated: ", (product(sz)*4)+1
 #    write(17,*) "Space required without symmetric storage: ", (product(sz)*7)+1
 #    write(17,*) "Space required with standard array storage: ", &
 #         &real(product(sz))**2
 #
 #  end subroutine sprsymqckprdc
 #
 #  subroutine av(n,x,y)
 #    implicit none
 #    integer,intent(in) :: n
 #    real, dimension(n), intent(in) :: x
 #    real, dimension(n), intent(out) :: y
 #    integer :: i,j
 #
 #    if (ija(1).ne.(n+2)) then
 #       write(6,*) "Matrix and vector size mismatch"
 #       call abort
 #    end if
 #
 #    do i = 1, n
 #       y(i) = sa(i)*x(i)
 #       do j = ija(i), (ija(i+1)-1)
 #          y(i) = y(i) + sa(j)*x(ija(j))
 #       end do
 #    end do
 #
 #    do i = 1, n
 #       do j = ija(i), (ija(i+1)-1)
 #          y(ija(j)) = y(ija(j)) + sa(j)*x(i)!Use for symetric storage
 #       end do
 #    end do
 #
 #  end subroutine av
 #
 # subroutine avc(n,x,y)
 #    implicit none
 #    integer,intent(in) :: n
 #    complex, dimension(n), intent(in) :: x
 #    complex, dimension(n), intent(out) :: y
 #    integer :: i,j
 #
 #    if (ija(1).ne.(n+2)) then
 #       write(6,*) "Matrix and vector size mismatch"
 #       call abort
 #    end if
 #
 #    do i = 1, n
 #       y(i) = cmplx(sa(i),0.0)*x(i)
 #       do j = ija(i), (ija(i+1)-1)
 #          y(i) = y(i) + cmplx(sa(j),0.0)*x(ija(j))
 #       end do
 #    end do
 #
 #    do i = 1, n
 #       do j = ija(i), (ija(i+1)-1)
 #          y(ija(j)) = y(ija(j)) + cmplx(sa(j),0.0)*x(i)
 #          !Use for symetric storage
 #       end do
 #    end do
 #  end subroutine avc