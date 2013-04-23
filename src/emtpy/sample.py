__author__ = 'pard'
from emtpy.grid import Grid


class potential_energy(Grid):
    pass

def generate(sample, radius, clumploc, mlheight, genPos):
    """
        This is the routine which generates a given potential
    """
#     implicit none
#     integer, intent(in) :: sample,mlheight
# !    real, dimension(3),intent(in) :: maxsize
# !    real, intent(in) :: upperinterface,lowerinterface,rate,
#     real, intent(in) :: radius,clumploc
#     real, dimension(3,maxpos),optional,intent(out) :: genPos
# !   complex,optional,intent(out) :: chi(0:ksize(1)-1,0:ksize(2)-1,0:ksize(3)-1)
#     real, allocatable, dimension(:) :: disp,T
#     real,dimension(3) :: zerovec
#     real :: width, Es,p,leak,FWHM, temp
#     integer :: i,j,k, ul, ll,m,q,varsp(3)
#     real, allocatable,dimension(:,:,:) :: Xf, Xphi
#     integer, allocatable, dimension(:) :: bottom, stepprof,Xprof
#     integer, allocatable, dimension(:,:) :: top
#     character(40) :: filename
#     write(filename,'(A17,A11)') rundirroot,'profile.dat'
#     open(unit=25,file=filename,status='Unknown')
#     write(filename,'(A17,A7)') rundirroot,'pos.dat'
#     open(unit=36,file=filename,status='Unknown')
#     call randinit
#     !Open files and initialise random number generator
# !!$
# !!$    if (present(chi) .and. sample .ne. 8) then
# !!$       write(6,*) 'WARNING - Invalid argument in call to routine: generate'
# !!$       call abort
# !!$    end if
#
#     temp = 0.0d0
#
#     leak = background
#     ll = nint(qwa/latt(adim(3)))
#     ul = nint(qwb/latt(adim(3)))
#     spacings = nint(Rsize/latt)
#     do i = 1, 3
#        varsp(i) = spacings(adim(i))
#     end do
#     write(6,*) 'Dimensions of sample (in units of atomic sites): '
#     write(6,*) spacings
#     write(6,*) 'Total number of atomic sites: ', product(spacings)
#     allocate(Xf(0:spacings(1)-1,0:spacings(2)-1,0:spacings(3)-1))
#     allocate(Xphi(0:spacings(1)-1,0:spacings(2)-1,0:spacings(3)-1))
#     allocate(bottom(0:spacings(adim(1))-1))
#     allocate(Xprof(0:spacings(direction)-1))
#     allocate(T(0:spacings(direction)-1))
#     allocate(top(0:spacings(adim(1))-1,0:spacings(adim(2))-1))
#     Xphi = 0.0d0
#     Xprof = 0
#     Xf = 0
#     genPos = 0
#     totA = 0
#     top = ul
#     bottom = ll
#     if (sample == 4) read(38,*) temp
#     T = temp
# !    T(ll:ul) = 710
#     T = T + 273.0
# !    FWHM = 2.860d0
#     FWHM = qwb-qwa!upperinterface-lowerinterface
#
#     if (radius .gt. 0.0d0) call clumps(top,radius,clumploc,mlheight,varsp(1:2))
#
#     if (sample == 2 .or. sample == 12) then
#        call readprofile(Xphi,FWHM,nom)
#        XF=Xphi
#     else if (sample == 4) then
#        call flowrate(nom,leak,Xphi,bottom,top,varsp(1:2))
# !       call mydiffusion(Xf,Xphi,T)
#        call diffusionnew(Xf,Xphi,T)
# !       Xf=Xphi
#     else if (sample == 6) then
#        call flowrate(leak,nom,Xphi,bottom,top,varsp(1:2))
#        Xf = Xphi
#        !eg for AlxGa1-xN where the Al is in the buffer
# !!$    else if (sample == 8) then
# !!$!       call steps(bstep,tstep,bottom,top)
# !!$!      call clumps(top,radius,clumploc,1)
# !!$!       call voids(top,radius,clumploc,-1)
# !!$       call flowrate(rate,leak,Xphi,bottom,top)
# !!$       Xf=Xphi
# !!$       chi=Xphi
#     else if (sample == 0) then
#        call lumpywell(leak,nom,Xphi,bottom,top,radius)
#        Xf=Xphi
#     else
#        call flowrate(nom,leak,Xphi,bottom,top,varsp(1:2))
#        Xf = Xphi
#     end if
#
#     if (sample .ne. 8) then
#        do k = 0, spacings(3)-1
#           do j = 0, spacings(2)-1
#              do i = 0, spacings(1)-1
#                 call random_number(p)
#                 if (p.le.Xf(i,j,k)) then
#                    totA = totA + 1
#                    genPos(1,totA) = i*latt(1)!ap(1)
#                    genPos(2,totA) = j*latt(2)!ap(2)
#                    genPos(3,totA) = k*latt(3)!ap(3)
#                    if (direction == 1) Xprof(i) = Xprof(i) + 1
#                    if (direction == 2) Xprof(j) = Xprof(j) + 1
#                    if (direction == 3) Xprof(k) = Xprof(k) + 1
#                    write(36,'(3G24.15E3)') genPos(:,totA)
#                 end if
#              end do
#           end do
#        end do
#        !use 3-d dist as a prob to place atoms randomly
#     end if
#
#     do i = 0, spacings(direction)-1
#           write(25,12) i*latt(direction), real(Xprof(i))&
#                &/real((spacings(adim(1)))*(spacings(adim(2)))),&
#                &Xf(0,0,i),Xphi(0,0,i)
#     end do
#
#     deallocate(Xf)
#     deallocate(Xphi)
#     deallocate(bottom)
#     deallocate(Xprof)
#     deallocate(T)
#     deallocate(top)
#
#     close(25)
#     close(36)
# 12  format(4G24.15E3)
#   end subroutine generate
    pass


def lumpy_well():
    pass


def diffpara(profparam, T):
    from math import exp
    from emtpy.constants import kb, eV
    Es = profparam*eV
    return exp(Es/(kb*T))

def flowrate(rate,leak,Xphi,bottom,top,sp):
#     implicit none
#     integer, intent(in) :: sp(2)
#     real, intent(in) :: rate,leak
#     real,intent(out)::Xphi(0:spacings(1)-1,0:spacings(2)-1,0:spacings(3)-1)
# !    real,intent(out)::Xphi(:,:,:)
# !    integer, dimension(:),intent(in) :: bottom
# !    integer, dimension(:,:),intent(in) :: top
#     integer, dimension(0:sp(1)-1),intent(in) :: bottom
#     integer, dimension(0:sp(1)-1,0:sp(2)-1),intent(in) :: top
#     integer :: i,j
#     Xphi = leak
#     if (direction .eq. 1) then
#        do j = 0, spacings(3)-1
#           do i = 0, spacings(2)-1
#              Xphi(bottom(j):top(j,i),i,j) = rate
#           end do
#        end do
#     else if (direction .eq. 2) then
#        do j = 0, spacings(3)-1
#           do i = 0, spacings(1)-1
#              Xphi(i,bottom(i):top(i,j),j) = rate
#           end do
#        end do
#     else if (direction .eq. 3) then
#        do j = 0, spacings(2)-1
#           do i = 0, spacings(1)-1
#              Xphi(i,j,bottom(i):top(i,j)) = rate
#           end do
#        end do
#     end if
#   end subroutine flowrate
    pass

def diffusionnew(X,Xphi,T):
  #   implicit none
  #   real,intent(out)::X(0:spacings(1)-1,0:spacings(2)-1,0:spacings(3)-1)
  #   real,intent(in)::Xphi(0:spacings(1)-1,0:spacings(2)-1,0:spacings(3)-1)
  #   real,dimension(0:spacings(1)-1,0:spacings(2)-1,0:spacings(3)-1)::Xs,Xb
  #   real, dimension(0:spacings(3)-1),intent(in) :: T
  #   integer :: i,j,k
  #   real :: a,b,aminus
  #   Xs(0,0,0) = Xphi(0,0,0)
  #   Xb = 0.0d0
  #   write(6,*) diffpara(T(0))
  #   do k = 1, spacings(3)-1
  #      do j = 0, spacings(2)-1
  #         do i = 0, spacings(1)-1
  #            a = diffpara(T(k))
  #            aminus=a-1.0d0
  #            b = Xs(i,j,k-1) + XPhi(i,j,k)
  #            Xs(i,j,k) = (1.0d0+a-b+(a*b)-&
  #                 &sqrt((4.0d0*aminus*b)+(1.0d0+a+b-(a*b))**2))/(2.0d0*aminus)
  #            Xb(i,j,k) = -(1.0d0+a+b-(a*b)-&
  #                 &sqrt((4.0d0*aminus*b)+(1.0d0+a+b-(a*b))**2))/(2.0d0*aminus)
  #         end do
  #      end do
  #   end do
  #   X = Xb
  # end subroutine diffusionnew
    pass

def steps(bstep,tstep,bottom,top):
  #   implicit none
  #   integer, intent(in) ::  bstep, tstep
  #   integer, dimension(0:spacings(1)-1),intent(inout) :: bottom
  #   integer, dimension(0:spacings(1)-1,0:spacings(2)-1),intent(inout) :: top
  #   integer :: i
  #   do i = 0, spacings(1)-1
  #      if (mod(i+nint(70/latt(adim(1))),bstep) == 1 .and. i .ne. 1) then
  #         bottom(i:spacings(1)-1) = bottom(i:spacings(1)-1) + 1
  #      end if
  #      if (mod(i+nint(70/latt(adim(1))),tstep) == 1 .and. i .ne. 1) then
  #         top(i:spacings(1)-1,0:spacings(2)-1) = &
  #              &top(i:spacings(1)-1,0:spacings(2)-1)+1
  #      end if
  #   end do
  # end subroutine steps
    pass

def clumps(top,radius,clumploc,height,sp):
  #   implicit none
  #   integer, intent(in) :: height, sp(2)
  #   integer, dimension(0:sp(1)-1,0:sp(2)-1),intent(inout) :: top
  #   real, intent(in) :: radius,clumploc
  #   !This is the height of the WWF in ML, it can be negative for 'voids'
  #   integer, dimension(2) :: clumpcenter
  #   integer :: l,m
  #   real :: rs
  #   do l = -(nint(radius/latt(adim(1)))),(nint(radius/latt(adim(1))))
  #      do m = -nint(radius/latt(adim(2))), nint(radius/latt(adim(2)))
  #         rs = sqrt((real(l)*latt(adim(1)))**2 + (real(m)*latt(adim(2)))**2)
  #         if (rs .le. radius) then
  #            clumpcenter(1) = l+(nint(clumploc/latt(adim(1))))
  #            clumpcenter(2) = m+(nint(clumploc/latt(adim(2))))
  #            top(clumpcenter(1),clumpcenter(2))= &
  #                 &top(clumpcenter(1),clumpcenter(2)) + height
  #         end if
  #      end do
  #   end do
  # end subroutine clumps
    pass

def sqclumps(top,radius,clumploc,height,sp):
  #   implicit none
  #   integer, intent(in) :: height, sp(2)
  #   integer, dimension(0:sp(1)-1,0:sp(2)-1),intent(inout) :: top
  #   real, intent(in) :: radius,clumploc
  #   !This is the height of the WWF in ML, it can be negative for 'voids'
  #   integer, dimension(2) :: clumpcenter
  #   integer :: l,m
  #   do l = -(nint(radius/latt(adim(2)))),(nint(radius/latt(adim(2))))
  #      do m = -nint(radius/latt(adim(1))), nint(radius/latt(adim(1)))
  #         clumpcenter(1) = l+(nint(clumploc/latt(adim(2))))
  #         clumpcenter(2) = m+(nint(clumploc/latt(adim(1))))
  #         top(clumpcenter(1),clumpcenter(2))= &
  #              &top(clumpcenter(1),clumpcenter(2)) + height
  #      end do
  #   end do
  # end subroutine sqclumps
    pass

def samplepos(Rlist,arsize,sampler):
#     implicit none
#     !!!TAKE CARE THINKING ABOUT WETHER TO
#     !!!USE A OR A0 WHEN USING REAL DATA
#     real, dimension(3,maxpos),intent(out) :: Rlist
#     integer,dimension(3), intent(in) :: arsize
#     integer, intent(in) :: sampler
#     integer :: n,i
# !    integer, dimension(4) :: ifail
#     integer :: ifail
#     integer, dimension(arsize(3)) :: profhist
#     real :: t
#     real, dimension(3) :: r,offset
#     character(28) :: filename
#     write(filename,'(A17,A11)') rundirroot,'profile.dat'
#     open(unit=25,file=filename,status='Unknown')
#     profhist = 0
#     n = 0
#     offset = 0.0
#     if (sampler == 1) then
#        offset(1:2) = 5.0
#        offset(3) = -18.0
#     else! if (sampler == 2) then
#        offset(3) = 7.0
#     end if
#     do
#        r(1) = datumpos(ifail)
#        r(2) = datumpos(ifail)
#        r(3) = datumpos(ifail)
#        t = datumpos(ifail)
#        !Read positions from file
#        if (ifail == 1) exit
#        r = r + offset
#        if (t .ge. 110.0 .and. t .le. 120.0) then
#           !If the A atom lies in the region we want it is added to R
#           n = n + 1
# !!$          R(1,n) = a0*anint(x(1)/a0)
# !!$          R(2,n) = a0*anint(x(2)/a0)
# !!$          R(3,n) = c0*anint(x(3)/c0)
#           Rlist(1,n) = aBN*anint(r(1)/aBN)
#           Rlist(2,n) = aBN*anint(r(2)/aBN)
#           Rlist(3,n) = cBN*anint(r(3)/cBN)
#           !Position vectors placed on the lattice
#           profhist(nint(r(3)/cBN)) = profhist(nint(r(3)/cBN)) + 1
#        else if (t .ge. 55.0 .and. t .le. 60.0) then
#           !If the A atom lies in the region we want it is added to R
#           n = n + 1
# !!$!          R(1,n) = a0*anint(x(1)/a0)
# !!$!          R(2,n) = a0*anint(x(2)/a0)
# !!$!          R(3,n) = c0*anint(x(3)/c0)
#           Rlist(1,n) = aBN*anint(r(1)/aBN)
#           Rlist(2,n) = aBN*anint(r(2)/aBN)
#           Rlist(3,n) = cBN*anint(r(3)/cBN)
#           !Position vectors placed on the lattice
#           profhist(nint(r(3)/cBN)) = profhist(nint(r(3)/cBN)) + 1
#        end if
#     end do
#     do i = 1, arsize(3)
#        write(25,*) real(i)*cBN, real(profhist(i))/real(arsize(1)*arsize(2))
#     end do
#     totA = n
#     close(25)
#     return
#   end subroutine samplepos
    pass

def sqwella(R,za,zb,Ano,sz,ag):
  #   implicit none
  #   integer, intent(in) :: sz(3), ag(3)
  #   real, dimension(3,20000),intent(out) :: R
  #   real, intent(in) :: za,zb
  #   integer, intent(out) :: Ano
  #   integer :: n,l,m,p
  #   real,dimension(3) :: x
  #   n = 0
  #   do l = 0, sz(1)-1,ag(1)
  #      do m = 0, sz(2)-1,ag(2)
  #         do p = 0, sz(3)-1,ag(3)
  #            x(1) = l*grid(ag,1)
  #            x(2) = m*grid(ag,1)
  #            x(3) = p*grid(ag,3)
  #            if ((x(3) .gt. za) .and. (x(3) .lt. zb)) then
  #               n = n + 1
  #               R(1,n) = grid(ag,1)*anint(x(1)/grid(ag,1))
  #               R(2,n) = grid(ag,1)*anint(x(2)/grid(ag,1))
  #               R(3,n) = grid(ag,3)*anint(x(3)/grid(ag,3))
  #            end if
  #         end do
  #      end do
  #   end do
  #   Ano = n
  # end subroutine sqwella
    pass

def sqwellb(A,za,zb,valwell,valbarr,sz,ag):
  #   implicit none
  #   integer, intent(in) :: sz(3), ag(3)
  #   complex,intent(inout),dimension(0:sz(1)-1,0:sz(2)-1,0:sz(3)-1)::A
  #   real, intent(in) :: za,zb,valwell,valbarr
  #   integer :: l,m,p
  #   real,dimension(3) :: x
  #   do l = 0, sz(1)-1
  #      do m = 0, sz(2)-1
  #         do p = 0, sz(3)-1
  #            x(1) = l*grid(ag,1)
  #            x(2) = m*grid(ag,1)
  #            x(3) = p*grid(ag,3)
  #            if ((x(direction) .ge. za) .and. (x(direction) .le. zb)) then
  #               A(l,m,p) = valwell
  #            else
  #               A(l,m,p) = valbarr
  #            end if
  #         end do
  #      end do
  #   end do
  # end subroutine sqwellb
    pass

def sqwell(A,wella,wellb,valwell,valbarr,sz,ag):
#     implicit none
#     integer, intent(in) :: sz(3), ag(3)
#     complex,intent(inout),dimension(0:sz(1)-1,0:sz(2)-1,0:sz(3)-1)::A
#     real, intent(in) :: wella,wellb,valwell,valbarr
#     integer :: l,m,p
#     real,dimension(3) :: x
#     do l = 0, sz(1)-1
#        do m = 0, sz(2)-1
#           do p = 0, sz(3)-1
#              x(1) = l*grid(ag,1)
#              x(2) = m*grid(ag,1)
#              x(3) = p*grid(ag,3)
# !             if (direction == 1) then
#              if ((x(direction) .ge. wella) .and. &
#                   &(x(direction) .le. wellb)) then
#                 A(l,m,p) = valwell
#              else
#                 A(l,m,p) = valbarr
#              end if
#           end do
#        end do
#     end do
#   end subroutine sqwell
    pass

def dblwell(A,widthwell,widthsep,valsep,valbarr,sz,ag):
#     implicit none
#     integer, intent(in) :: sz(3), ag(3)
#     complex,intent(inout),dimension(0:sz(1)-1,0:sz(2)-1,0:sz(3)-1)::A
#     real, intent(in) :: widthwell,widthsep,valsep,valbarr
#     integer :: l,m,p
#     integer,dimension(3) :: n
#     real :: y
#     do l = 0, sz(1) - 1
#        do m = 0, sz(2)-1
#           do p = 0, sz(3)-1
#              n(1) = l
#              n(2) = m
#              n(3) = p
# !             x(1) = (l-sz(1)/2)*grid(ag,1)
# !             x(2) = (m-sz(2)/2)*grid(ag,1)
# !             x(3) = (p-sz(3)/2)*grid(ag,3)
#              y = (n(direction)-sz(direction)/2)*grid(ag,direction)
#              if (abs(y) .gt. (widthsep/2.0)) then
#                 if (abs(y) .le. ((widthsep/2.0)+widthwell+&
#                      &(heaviside(y)*1.0))) then
#                    A(l,m,p) = 0.0d0 !valwell
#                 else
#                    A(l,m,p) = valbarr
#                 end if
#              else
#                 A(l,m,p) = valsep
#              end if
#           end do
#        end do
#     end do
#   end subroutine dblwell
    pass

def dblwellsym(A,widthwell,widthsep,valsep,valbarr,sz,ag):
  #   implicit none
  #   integer, intent(in) :: sz(3), ag(3)
  #   complex,intent(inout),dimension(0:sz(1)-1,0:sz(2)-1,0:sz(3)-1)::A
  #   real, intent(in) :: widthwell,widthsep,valsep,valbarr
  #   integer :: l,m,p
  #   real,dimension(3) :: x
  #   do l = 0, sz(1) - 1
  #      do m = 0, sz(2)-1
  #         do p = 0, sz(3)-1
  #            x(1) = (l-sz(1)/2)*grid(ag,1)
  #            x(2) = (m-sz(2)/2)*grid(ag,1)
  #            x(3) = (p-sz(3)/2)*grid(ag,3)
  #            if (abs(x(direction)) .gt. (widthsep/2.0)) then
  #               if (abs(x(direction)) .le. ((widthsep/2.0)+widthwell)) then
  #                  A(l,m,p) = 0.0d0 !valwell
  #               else
  #                  A(l,m,p) = valbarr
  #               end if
  #            else
  #               A(l,m,p) = valsep
  #            end if
  #         end do
  #      end do
  #   end do
  # end subroutine dblwellsym
    pass

def mqwell(A,widthwell,widthbarr,edge,valwell,valbarr,sz,ag):
#     implicit none
#     integer, intent(in) :: sz(3), ag(3)
#     complex,intent(inout),dimension(0:sz(1)-1,0:sz(2)-1,0:sz(3)-1)::A
#     real, intent(in) :: widthwell,widthbarr,edge,valwell,valbarr
#     integer :: i,j,k, wwell,wbarr
#     real,dimension(3) :: x
#     real :: disp,wells
#     wwell = nint(widthwell/grid(ag,direction))
# !    wbarr = nint(widthbarr/latt(direction))
# !    edgeint = nint(edge/latt(direction))
#     wells = 0.0d0
# !    write(6,*) edge, Rsize(direction+3)-edge, valwell,valbarr
#     write(6,*) widthwell,widthbarr,edge,valwell,valbarr
#     write(6,*) 'well width: ',wwell, widthwell,latt(direction)
#     do k = 0, sz(3)-1
#        do j = 0, sz(2)-1
#           do i = 0, sz(1)-1
#              x(1) = i*grid(ag,1)
#              x(2) = j*grid(ag,1)
#              x(3) = k*grid(ag,3)
#              if ((x(direction) .le. ((nint(Rsize(direction+3)/latt(direction))&
#                   &*latt(direction))-edge)) .and. &
#                   &(x(direction) .ge. edge)) then
#                 disp = x(direction)-edge-&
#                      &aint(nint(10000.0*wells)/10000.0)*(widthbarr+widthwell)
# !                write(6,*) x(direction),disp, wells,aint(wells)
#                 if (disp .ge. 0.0 .and. disp .le. widthwell) then
#                    A(i,j,k) = valwell
#                    wells = wells + (grid(ag,direction)/widthwell)
#                 else
#                    A(i,j,k) = valbarr!*0.6
#                 end if
#              else
#                 A(i,j,k) = valbarr
#              end if
#           end do
#        end do
#     end do
#   end subroutine mqwell
    pass

