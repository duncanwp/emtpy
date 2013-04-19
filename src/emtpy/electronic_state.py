__author__ = 'pard'

class state(object):
    pass

#
#   real function occprob(n,t,ntot,eigenval,leigenvec,reigenvec,neigenvals)
#     implicit none
#     integer, intent(in) :: n,ntot,neigenvals
#     real, intent(in) :: t,leigenvec(ntot,neigenvals),reigenvec(ntot,neigenvals)
#     complex, intent(in) :: eigenval(neigenvals)
#     integer :: i,m
# !    real :: sm
#     occprob = 0.0d0
#     do i = 1, neigenvals
# !       sm = 0.0d0
#        do m = 1, ntot
# !          sm = sm + (leigenvec(m,i))/real(ntot)
#           occprob = occprob + exp(-eigenval(i)*t)*&
#                &(reigenvec(n,i))*(leigenvec(m,i))/real(ntot)
#        end do
# !       occprob = occprob + exp(-eigenval(i)*t)*(reigenvec(n,i))*sm
#     end do
#   end function occprob
#
#   subroutine eigenout(eigenval,eigenvec,local,conf,para,X,sz,ag,off)
#     implicit none
#     integer, intent(in) :: sz(3), ag(3),off(3),conf
#     real, intent(in),dimension(conf,2) :: eigenval
#     complex, intent(inout), dimension(product(sz),conf) :: eigenvec
#     real,dimension(0:sz(1)-1,0:sz(2)-1,0:sz(3)-1) :: temp
#     complex,dimension(0:sz(1)-1,0:sz(2)-1,0:sz(3)-1) :: tempc
#     character(2),intent(in) :: para
#     real, dimension(product(sz)),intent(in) :: X
#     real, dimension(26) :: eigen
#     character(40) :: filename
#     real,dimension(4) :: intunder,volunder,chiunder,values
#     real,intent(out), dimension(conf,6) :: local
#     logical,dimension(product(sz)) :: mask
#     integer :: i,slice,n,m,paranum,unit
#     real :: psisqm,four,fg(3),ssize(3),psq(3),stddevp(3),econst
#     complex :: ksq(3),p(3),k(3),stddevk(3)
#
#     econst = hbar**2/(2.0d0*me*eV)
#
#     do i = 1, 3
#        fg(i) = grid(ag,i)
#     end do
#
#     ssize = fg*off
#
#     if (direction .eq.3) slice = 1
#     if (direction .eq. 1) slice = 2
#     if (direction .eq. 2) slice = 3
#
#     if (para == 'hh') paranum = 1
#     if (para == 'lh') paranum = 2
#     if (para == 'sh') paranum = 3
#     if (para == '_h') paranum = 0
#     if (para == '_e') paranum = -1
#
#     if (para == '_e') then
#        unit = 23
#     else
#        unit = 24
#     end if
#
#     do m = 1, conf
#
#        psisqm = maxval((abs(eigenvec(:,m)))**2)
#
#        values(2) = exp(-1.5)*psisqm
#        values(1) = 0.5*psisqm
#        values(3) = exp(-3.0)*psisqm
#        values(4) = 0.0d0
#        !FOR 3-D Guassians these give 1st & 2nd standard deviations (& FWHM)
#        !THIS DOESN'T SEEM TO BE TRUE, it gives something like them, but be
#        !carefull when saying anything like this as a multivariate normal
#        !distribution has a covariance MATRIX, as opposed to 3 values for
#        !sigma.Just quote the actual probability under each iso-surface.
#
#        do i = 1, 3
#           mask = (abs(eigenvec(:,m)))**2 .ge. values(i)
#           !conjg(eigenvec(:,m))*eigenvec(:,m)=(abs(eigenvec(:,m))**2)
#           !Except the abs way of doing it explicitly gives a real number
#           intunder(i) = sum((abs(eigenvec(:,m)))**2,mask)*product(fg)
#           volunder(i) = real(count(mask))*product(fg)
#           chiunder(i) = sum(X,mask)/real(count(mask))
#        end do
#
#        if (m .lt. 6 .or. nodimensions .eq. 1) then
#           write(filename,'(A24,A6,I2.2,A2,A4)') rundir, 'Eigenv',m,para,'.dat'
#           !Can also use: filename = 'Eigenv' // m // para // '.dat'
#           open(15,file=filename,status='Unknown',buffered='yes')
#
#           write(filename,'(A24,A7,I2.2,A2,A4)') rundir,'Eigenvw',m,para,'.dat'
#           open(16,file=filename,status='Unknown',buffered='yes')
#
# 8          write(filename,'(A24,A8,I2.2,A2,A4)') rundir,'Eigenvc1',m,para,'.dat'
#           open(39,file=filename,status='Unknown',buffered='yes')
#
#           write(filename,'(A24,A8,I2.2,A2,A4)') rundir,'Eigenvc2',m,para,'.dat'
#           open(50,file=filename,status='Unknown',buffered='yes')
#
#           write(filename,'(A24,A8,I2.2,A2,A4)') rundir,'Eigenvc3',m,para,'.dat'
#           open(51,file=filename,status='Unknown',buffered='yes')
#
#           call onedtothreed(tempc,eigenvec(:,m),sz)
#           !temp = aimag(tempc)
#           temp = real(tempc)
#
#           call outputwave(temp,15,slice,sz(slice)/2,eigenval(m,1),sz,ag)
#           call outputwave(temp,39,4,0,0.0,sz,ag,values(1))
#           call outputwave(temp,50,4,0,0.0,sz,ag,values(2))
#           call outputwave(temp,51,4,0,0.0,sz,ag,values(3))
#           call outputwave(temp,16,direction,sz(direction)/2,&
#                &eigenval(m,1),sz,ag)
#
#           close(15)
#           close(16)
#           close(39)
#           close(50)
#           close(51)
#        end if
#
#        four=sum(real(eigenvec(:,m))**4)/(normalization(ag)**2)
#
#        call expectedpos(eigenvec(:,m),local(m,1:3),local(m,4:6),sz,ag)
# !       call expectedk(eigenvec(:,m),eigenvec(:,m),k,sz,ag)
#        call expectedksq(eigenvec(:,m),eigenvec(:,m),k,ksq,stddevk,sz,ag)
# !       call pmelement(eigenvec(:,m),eigenvec(:,m),p,sz,ag)
#        call psqmelement(eigenvec(:,m),eigenvec(:,m),p,psq,stddevp,sz,ag)
#        eigen(1) = run
#        eigen(2) = m
#        eigen(3) = eigenval(m,1)
#        eigen(4) = eigenval(m,1) - (maxenergy/2.0d0)
#        eigen(5) = eigenval(m,2)
#        eigen(6:8) = local(m,1:3) + ssize
#        eigen(9:11) = local(m,4:6)
#        eigen(12) = abs(k(adim(3)))
#        eigen(13) = abs(p(adim(3)))
#        eigen(14) = (hbar)*eigen(12)/eigen(13)
#        eigen(15) = sqrt(local(m,adim(1)+3)*local(m,adim(2)+3))
#        eigen(16) = 1.0/(four)**(1.0/3.0)
# !       eigen(17:19) = intunder(1:3)
#        eigen(17) = ksq(adim(1)) + ksq(adim(2))
#        eigen(18) = econst*(ksq(1) + ksq(2))*invmassxy(nom,para)
#        eigen(19) = econst*ksq(3)*invmassz(nom,para)
#        !Assuming effective mass for nominal QW concentration
#        eigen(20:22) = volunder(1:3)
#        eigen(23:25) = chiunder(1:3)
#        eigen(26) = paranum
#        write(unit,50) (eigen(i),i=1,26)
#     end do
# 50  format(26G24.15E3)
#   end subroutine eigenout
#
#   subroutine melements(eeval,eevec,heval,hevec,en,hn,funit,sz,ag,hnums)
#     implicit none
#     integer, intent(in) :: sz(3), ag(3)
#     integer, intent(in) :: en, hn,funit
#     complex, intent(in), dimension(product(sz),en) :: eevec
#     real, intent(in), dimension(en,2) :: eeval
#     complex, intent(in), dimension(product(sz),hn) :: hevec
#     real, intent(in), dimension(hn,2) :: heval
#     integer,dimension(3),optional, intent(in) :: hnums
#     real, dimension(en,hn,8) :: matrix
#     complex, dimension(3) :: p
# !    complex :: ksq(3),p(3),k(3),stddevk(3)
#     real :: energy,sqoverl!,psq(3),stddevp(3)
#     integer :: n, m, paranum, j
#
#     do n = 1, en
#        do m = 1, hn
#
#           energy = eeval(n,1)-heval(m,1)
#
#           if (energy .gt. 0.0d0) then
#           !This stops the ouput of 'self overlaps'
#
#              if (present(hnums)) then
#                 if ((m .ge. hnums(1)) .and. m .lt. hnums(2)) then
#                    paranum  = 1
#                 else if ((m .ge. hnums(2)) .and. m .lt. hnums(3)) then
#                    paranum = 2
#                 else if (m .ge. hnums(3)) then
#                    paranum = 3
#                 else
#                    write(6,*) 'Error in choosing hole band'
#                    call abort
#                 end if
#              else
#                 paranum = 0
#              end if
#
#              call pmelement(eevec(:,n),hevec(:,m),p,sz,ag)
# !!$             call psqmelement(eevec(:,n),hevec(:,m),p,psq,stddevp,sz,ag)
# !!$             call expectedk(eevec(:,n),hevec(:,m),k,sz,ag)
# !!$             call expectedksq(eevec(:,n),hevec(:,m),k,ksq,stddevk,sz,ag)
#              !This is k_he = <e|k|h>
#              !p is the momentum matrix element, <e|p|h>
#              sqoverl = (abs((dot_product(eevec(:,n),hevec(:,m))/&
#                   &(normalization(ag)**2))))**2
#              !This is the squared overlap integeral, |Sz*z d^3r|^2
#
#              matrix(n,m,1) = energy
#              matrix(n,m,2) = sqrt(eeval(n,2)**2 + heval(m,2)**2)
#              ! Adding the error of Eh and Ee
#              matrix(n,m,3) = wavelength(energy)
#              matrix(n,m,4) = dot_product(p,p)
#              !The square of the momentum matrix element, |<f|p|i>|^2
#              matrix(n,m,5) = sqoverl
#              !Square of the overlap integral
#              matrix(n,m,6) = direlement(eeval(n,1),eevec(:,n),&
#                   &heval(m,1),hevec(:,m),1,sz,ag)
#              matrix(n,m,7) = direlement(eeval(n,1),eevec(:,n),&
#                   &heval(m,1),hevec(:,m),2,sz,ag)
#              !These calculate <e|M|h> = M_he
#              !Fourier overlap elements for a q given by the delta
#              !function for LA and TA phonons
#              matrix(n,m,8) = paranum
#              !Number relating to band of the final state
#              write(funit,31) n,m,(matrix(n,m,j),j = 1, 8)
#           end if
#        end do
#     end do
# 31  format(10G24.15E3)
#   end subroutine melements
#
#   subroutine melementsx(eeval,eevecin,esz,eag,en,heval,hevecin,hsz,hag,hn,&
#        &fun,offs)
#     implicit none
#     integer, intent(in) :: esz(3), eag(3), hsz(3), hag(3),offs(3)
#     integer, intent(in) :: en, hn,fun
#     complex, intent(in), dimension(product(esz),en) :: eevecin
#     real, intent(in), dimension(en,2) :: eeval
#     complex, intent(in), dimension(product(hsz),hn) :: hevecin
#     real, intent(in), dimension(hn,2) :: heval
#     complex, dimension(product(hsz)*holecalcs) :: eevec,hevec
#     real, dimension(en,hn,8) :: matrix
#     complex, dimension(3) :: p!, k, ksq, stddevk
#     real :: energy,sqoverl!,psq(3),stddevp(3)
#     integer :: n, m, paranum, j,nwsz(3)
#
#     nwsz(3) = hsz(3)
#     nwsz(1:2) = hsz(1:2)*ehsratio
#
#     do n = 1, en
#        call threedinterpolate(eevecin(:,n),esz,eag,eevec,nwsz,hag)
#        do m = 1, hn
#
#           energy = eeval(n,1)-heval(m,1)
#
#           call lineuparr(hevecin(:,m),hsz,hag,hevec,nwsz,offs)
#
# !!$          if (present(hnums)) then
# !!$             if ((m .ge. hnums(1)) .and. m .lt. hnums(2)) then
# !!$                paranum  = 1
# !!$             else if ((m .ge. hnums(2)) .and. m .lt. hnums(3)) then
# !!$                paranum = 2
# !!$             else if (m .ge. hnums(3)) then
# !!$                paranum = 3
# !!$             else
# !!$                write(6,*) 'Error in choosing hole band'
# !!$                call abort
# !!$             end if
# !!$          else
#              paranum = 0
# !!$          end if
#
# !          call pmelement(eevec,hevec,mome,nwsz,hag)
#           call pmelement(eevec,hevec,p,nwsz,hag)
# !          call expectedksq(eevec,hevec,k,ksq,stddevk,nwsz,hag)
#           !p is the momentum matrix element, <h|p|e> = P_eh = P_cv
#
#           !Using p=hk/2pi instead, as the FFT will be more accurate than
#           !the finite difference method used to find p
#           !Not any more as the FFT uses too much memory...
#           sqoverl = (abs((dot_product(eevec,hevec)/&
#                &(normalization(hag)**2))))**2
#           !This is the squared overlap integeral, |Sz*z d^3r|^2
#
#           matrix(n,m,1) = energy
#           matrix(n,m,2) = sqrt(eeval(n,2)**2 + heval(m,2)**2)
#           ! Adding the error of Eh and Ee
#           matrix(n,m,3) = wavelength(energy)
#           matrix(n,m,4) = sqoverl
#           !Square of the overlap integral
#           matrix(n,m,5) = (abs(p(1)))**2
#           matrix(n,m,6) = (abs(p(2)))**2
#           matrix(n,m,7) = (abs(p(3)))**2
#           !These are, |<f|p(n)|i>|^2 for n=1,2,3
#           matrix(n,m,8) = paranum
#           !Number relating to band of the final state
#           write(fun,31) n,m,(matrix(n,m,j),j = 1, 8)
#        end do
#     end do
# 31  format(10G24.15E3)
#   end subroutine melementsx
#
#   subroutine threedinterpolate(vec,vsz,vevals,vecout,vszn,vevalsn)
#     implicit none
#     integer, intent(in),dimension(3) :: vsz,vevals,vszn,vevalsn
#     complex, intent(in) :: vec(product(vsz))
#     complex, intent(out) :: vecout(product(vszn))
#     complex :: tempin(0:vsz(1)-1,0:vsz(2)-1,0:vsz(3)-1)
#     complex :: tempout(0:vszn(1)-1,0:vszn(2)-1,0:vszn(3)-1)
#     real :: A, B, C, D, E, F, G, H
#     real,dimension(3) :: t, x, xia, xib,grdspc
#     integer :: xi(3), i, j, k,xp(3), n
#
#     call onedtothreed(tempin,vec,vsz)
#
#     do i = 1, 3
#        grdspc(i) = grid(vevals,i)
#     end do
#
#     do k = 0, vszn(3)-1
#        do j = 0, vszn(2)-1
#           do i = 0, vszn(1)-1
#              x(1) = real(i)*grid(vevalsn,1)
#              x(2) = real(j)*grid(vevalsn,2)
#              x(3) = real(k)*grid(vevalsn,3)
#              xi = int(x/grdspc)
#              xia = xi*grdspc
#              xib = (xi+1)*grdspc
#              do n = 1, 3
#                 if ((xi(n)+1) .ge. vsz(n)) then
#                    xp(n) = 0
#                 else
#                    xp(n) = xi(n)+1
#                 end if
#              end do
#              t = (x-xia)/(xib-xia)
#              A = tempin(xi(1),xi(2),xi(3))
#              B = tempin(xi(1),xi(2),xp(3))
#              C = tempin(xi(1),xp(2),xp(3))
#              D = tempin(xi(1),xp(2),xi(3))
#              E = tempin(xp(1),xp(2),xi(3))
#              F = tempin(xp(1),xp(2),xp(3))
#              G = tempin(xp(1),xi(2),xp(3))
#              H = tempin(xp(1),xi(2),xi(3))
#              tempout(i,j,k) = &
#                   &(((1.0d0-t(1))*A + t(1)*H)*(1.0d0-t(2)) + &
#                   &t(2)*((1.0d0-t(1))*B + t(1)*G))*(1.0d0-t(3)) +&
#                   &(((1.0d0-t(1))*D + t(1)*E)*(1.0d0-t(2)) + &
#                   &t(2)*((1.0d0-t(1))*C + t(1)*F))*t(3)
#              !Linear interpolation onto finer matrix grid
#           end do
#        end do
#     end do
#
#     call threedtooned(tempout,vecout,vszn)
#
#   end subroutine threedinterpolate
#
#   subroutine lineuparr(vec,vsz,vevals,vecout,vszn,off)
#     implicit none
#     integer, intent(in),dimension(3) :: vsz,vevals,off,vszn
#     complex, intent(in) :: vec(product(vsz))
#     complex, intent(out),dimension(product(vszn)) :: vecout
#     complex :: tempin(0:vsz(1)-1,0:vsz(2)-1,0:vsz(3)-1)
#     complex :: tempout(0:vszn(1)-1,0:vszn(2)-1,0:vszn(3)-1)
#     real :: A, B, C, D, E, F, G, H
#     real :: t(3), x(3), xia(3),xib(3)
#     integer :: xi(3)
#
#     call onedtothreed(tempin,vec,vsz)
#
#     tempout = 0.0d0
#     tempout(off(1):off(1)+vsz(1)-1,&
#          &off(2):off(2)+vsz(2)-1,&
#          &off(3):off(3)+vsz(3)-1) = tempin
#
#     call threedtooned(tempout,vecout,vszn)
#
#   end subroutine lineuparr
#
# !!$  subroutine melementsx(eeval,eevecin,heval,hevecin,en,hn,funit,&
# !!$       &esz,eag,hsz,hag,calc,hnums)
# !!$    implicit none
# !!$    integer, intent(in) :: esz(3), eag(3), hsz(3), hag(3)
# !!$    integer, intent(in) :: en, hn,funit
# !!$    complex, intent(in), dimension(product(esz),en) :: eevecin
# !!$    real, intent(in), dimension(en,2) :: eeval
# !!$    complex, intent(in), dimension(product(hsz),hn) :: hevecin
# !!$    real, intent(in), dimension(hn,2) :: heval
# !!$    integer,dimension(3),optional, intent(in) :: hnums
# !!$    complex, allocatable :: eevec, hevec
# !!$    complex::eevec(product(hsz)*holecalcs,en),hevec(product(hsz)*holecalcs,hn)
# !!$    real, dimension(en,hn,8) :: matrix
# !!$    complex, dimension(3) :: mome
# !!$    real :: energy,sqoverl
# !!$    integer :: n, m, paranum, j, nwsz, off(3),newksz(3)
# !!$
# !!$    newksz = hsz*holecalcs
# !!$    nwsz = product(newksz)
# !!$
# !!$    off(1)=nint((rsize(1)/real(ehsratio))*(calc-(holecalcoffy(calc)*ehsratio))&
# !!$         &/grid(hag,1))
# !!$    off(2)=nint((rsize(2)/real(ehsratio))*(holecalcoffy(calc)-1)/grid(hag,2))
# !!$    off(3) = 0
# !!$
# !!$    call lineuparr(eevecin,esz,eag,hevecin,hsz,hag,eevec,hevec,nwsz,newksz,off)
# !!$    !This routine enables overlaps to be calculated for arrays with different
# !!$    !sizes and spacings. On output both arrays are evaluated using the spacing
# !!$    !of the second array, and is the physical size of the first array.
# !!$    !This routine interpolates the coarser array using ...........?!
# !!$
# !!$    do n = 1, en
# !!$       do m = 1, hn
# !!$
# !!$          energy = abs(eeval(n,1)-heval(m,1))
# !!$
# !!$          if (energy .gt. 0.0d0) then
# !!$          !This stops the ouput of 'self overlaps'
# !!$
# !!$             if (present(hnums)) then
# !!$                if ((m .ge. hnums(1)) .and. m .lt. hnums(2)) then
# !!$                   paranum  = 1
# !!$                else if ((m .ge. hnums(2)) .and. m .lt. hnums(3)) then
# !!$                   paranum = 2
# !!$                else if (m .ge. hnums(3)) then
# !!$                   paranum = 3
# !!$                else
# !!$                   write(6,*) 'Error in choosing hole band'
# !!$                   call abort
# !!$                end if
# !!$             else
# !!$                paranum = 0
# !!$             end if
# !!$
# !!$             call pmelement(eevec(:,n),hevec(:,m),mome)
# !!$             !Mome is now the momentum matrix element, <e|p|h>
# !!$             sqoverl = (abs((dot_product(eevec(:,n),hevec(:,m))/&
# !!$                  &(normalization**2))))**2
# !!$             !This is the squared overlap integeral, |Sz*z d^3r|^2
# !!$
# !!$             matrix(n,m,1) = energy
# !!$             matrix(n,m,2) = sqrt(eeval(n,2)**2 + heval(m,2)**2)
# !!$             ! Adding the error of Eh and Ee
# !!$             matrix(n,m,3) = wavelength(energy)
# !!$             matrix(n,m,4) = sqrt(dot_product(mome,mome))
# !!$             !The magnitude of the momentum matrix element
# !!$             matrix(n,m,5) = sqoverl
# !!$             !Square of the overlap integral
# !!$             matrix(n,m,6) = direlement(eeval(n,1),eevec(:,n),&
# !!$                  &heval(m,1),hevec(:,m),1)
# !!$             matrix(n,m,7) = direlement(eeval(n,1),eevec(:,n),&
# !!$                  &heval(m,1),hevec(:,m),2)
# !!$             !Fourier overlap elements for a q given by the delta
# !!$             !function for LA and TA phonons
# !!$             matrix(n,m,8) = paranum
# !!$             !Number relating to band of the final state
# !!$             write(funit,31) n,m,(matrix(n,m,j),j = 1, 8)
# !!$          end if
# !!$       end do
# !!$    end do
# !!$31  format(10G24.15E3)
# !!$  end subroutine melementsx
#
#   real function findq(energy1,energy2,pol)
#     implicit none
#     real, intent(in) :: energy1, energy2
#     integer, intent(in) :: pol
#     findq = abs(energy2-energy1)/(hbarev*ssound(pol))
#     !This gives q in units of per meter
#     !Polarization = 1 is for LA
#     !             = 2 is for TA
#   end function findq
#
# !!$  real function fourierelementavq(e1,psi1,e2,psi2,polar,sz,ag)
# !!$    implicit none
# !!$    integer, intent(in) :: sz(3), ag(3)
# !!$    real, intent(in) :: e1, e2
# !!$    complex, intent(in) :: psi1(product(sz)), psi2(product(sz))
# !!$    integer, intent(in) :: polar
# !!$    real :: overlap(product(sz)), q,s(3)
# !!$    real :: int
# !!$    integer :: i,j,k
# !!$    int = 0.0d0
# !!$    q = findq(e1,e2,polar)
# !!$    !This is the magnitude of q, |q|
# !!$    overlap = conjg(psi1)*psi2
# !!$    do k = 0, sz(3)-1
# !!$       do j = 0, sz(2)-1
# !!$          do i = 0, sz(1)-1
# !!$             s(1) = (i*grid(ag,1))
# !!$             s(2) = (j*grid(ag,1))
# !!$             s(3) = (k*grid(ag,3))
# !!$             int = int + ((abs(overlap(nvalue(i,j,k,sz))))**2)*&
# !!$                  &exp(-cmplx(0.0d0,(q/sqrt(3.0d0))*sum(s)))
# !!$             !See notes regarding factor of 1/sqrt(3)
# !!$          end do
# !!$       end do
# !!$    end do
# !!$    fourierelementavq = int/(sqrt(2.0d0*pi)**3)
# !!$  end function fourierelementavq
#
#   real function direlement(e1,psi1in,e2,psi2in,polar,sz,ag)
#     !M_if = <f|M|i> = <1|M|2> = M_21
#     implicit none
#     integer, intent(in) :: sz(3), ag(3)
#     real, intent(in) :: e1, e2
#     complex, intent(in) :: psi1in(product(sz)), psi2in(product(sz))
#     complex :: psi1(product(sz)), psi2(product(sz))
#     integer, intent(in) :: polar
#     real :: magq,q(3),theta,phi,qgrid(3),sig,fg(3), norm
#     complex :: overlap(product(sz)),sumation
#     integer :: l,m,p,n(3),i,j
#     sumation = 0.0d0
#     magq = findq(e1,e2,polar)
#     !This is the magnitude of q, |q|
#     !Note, this gives q in units of m^-1
#
# !    write(6,*) 'qvalue:',e1,e2,magq
#
#     if (nodimensions == 1) then
#        qgrid = 1.0d0
#        fg = 1.0d0
#        fg(direction) = grid(ag,direction)*1E-9
#        qgrid(direction) = (2.0d0*pi)/(sz(direction)*fg(direction))
#        norm = 1E-9
#     else
#        do i = 1, 3
#           fg(i) = grid(ag,i)*1E-9
#           !Coordinates in m
#           qgrid(i) = (2.0d0*pi)/(sz(i)*fg(i))
#           !Coordinates in m-1
#           norm = 1E-27
#        end do
#     end if
#
#     psi1 = psi1in/sqrt(norm)
#     psi2 = psi2in/sqrt(norm)
#
#     overlap = conjg(psi1)*psi2
#
# !    write(6,*) (abs(sum(overlap)*product(fg)))**2
#
# !    sig = (2.0d0*pi)/(sz(3)*grid(ag,3)*1E-9)
#     sig = sqrt(sum(qgrid**2))
# !    write(6,*) 'std dev of "delta" function', sig
#     !This is a rough estimate for the thickness, the std dev is one fourier
#     !space lattice spacing.
#
#     call FFT(overlap,sz,1)
# !    call FFT(psi1,sz,1)
# !    call FFT(psi2,sz,1)
#     !Call FFT
#
# !    write(6,*) (abs(sum(overlap)*product(fg)/product(sz)))**2
# !    write(6,*) (abs(sum(conjg(psi1)*psi2)*product(fg)/product(sz)))**2
# !    write(6,*) dot_product(psi1,psi2)*product(fg)/product(sz)
#
#     do p = 0, sz(3)-1
#        do m = 0, sz(2)-1
#           do l = 0, sz(1)-1
#              n(1) = (mod((l + (sz(1)/2)),sz(1)) - sz(1)/2)
#              n(2) = (mod((m + (sz(2)/2)),sz(2)) - sz(2)/2)
#              n(3) = (mod((p + (sz(3)/2)),sz(3)) - sz(3)/2)
#              q = qgrid*n
# !             q(1) = (2.0d0*pi*n(1))/(sz(1)*grid(ag,1)*1E-9)
# !             q(2) = (2.0d0*pi*n(2))/(sz(2)*grid(ag,1)*1E-9)
# !             q(3) = (2.0d0*pi*n(3))/(sz(3)*grid(ag,3)*1E-9)
#              !Scale coordinates, in nm
# !             sumation = sumation + &
# !                  &conjg(psi1(nvalue(l,m,p,sz)))*psi2(nvalue(l,m,p,sz))!*&
#              sumation = sumation + (abs(overlap(nvalue(l,m,p,sz)))**2)*&
#                   &onedgauss(sqrt(dot_product(q,q)),sig,magq)
#           end do
#        end do
#     end do
#
# !    write(6,*) 'Phonon element: ',&!(abs(sumation))**2,&
# !!         &(abs(sumation*product(qgrid)))**2,&
# !         &(abs(sumation*product(fg)/product(sz)))**2
#     direlement = sumation*product(fg)/product(sz)
#
#   end function direlement
#
#   real function gldnprb(e1,e2,melement,T,band)
#     implicit none
#     real, intent(in) :: e1, e2,melement,T
#     integer, intent(in) :: band
#     real :: q, M, sum, delta,probp(2),phonon
#     integer :: polar,i
#     !e1,psi1 is the initial state, e2,psi2 is the final state
#     delta = e2-e1
#     if (delta .lt. 0.0d0) then
#        !i.e. the inital state is at a higher energy than the final state,
#        !so  it would need to emit/create a phonon
#        phonon = 1.0d0
#     else
#        !i.e. the inital state is at a lower energy than the final state,
#        !so  it would need to absorb/anhilate a phonon
#        phonon = 0.0d0
#     end if
#     do i = 1, 2
#        polar = i
#        q = findq(e1,e2,polar)
#        probp(polar) = (bosefuncq(q,T,polar)+phonon)*&
#             &acousticeav(q,polar,band)*melement/Omega(q,polar)
#     end do
#     gldnprb = sum(probp)*pi/(rho*((2.0d0*pi)**3))
#     !The 2/hbar in the golden rule cancels with the hbar/2 from melement
#   end function gldnprb
#
# !!$  real function acoustice(wavevec,pol)
# !!$    implicit none
# !!$    real, intent(in) :: wavevec(3)
# !!$    integer, intent(in) :: pol
# !!$    real :: magqsqd
# !!$    magqsqd = dot_product(wavevec,wavevec)
# !!$    acoustice = ((hydrodef**2)*magqsqd) !+ piezoel(wavevec,pol)**2
# !!$  end function acoustice
# !!$
# !!$  real function piezoel(q,polar)
# !!$    implicit none
# !!$    real, intent(in) :: q(3)
# !!$    integer, intent(in) :: polar
# !!$    real :: magqsqd
# !!$    magqsqd = dot_product(q,q)
# !!$    piezoel=0.0d0!?????
# !!$  end function piezoel
#
#   real function acousticeav(wavevec,pol,b)
#     implicit none
#     real, intent(in) :: wavevec
#     integer, intent(in) :: pol,b
#     real :: magqsqd
#     acousticeav = ((hydrodef(b)**2)*(wavevec**2)) + &
#          &piezoav(pol)*((eV*2.0d0*pi)/(epsilonr*epsilon0))**2
#   end function acousticeav
#
#   real function fermifunc(energy,ef,T)
#     implicit none
#     real, intent(in) :: energy,ef,T
#     !Energies in eV, T in K.
#     fermifunc = 1.0d0/(1.0d0+exp((energy-ef)/(KbeV*T)))
#     !Fermi-Dirac dist.
#   end function fermifunc
#
#   real function bosefuncq(q,T,pol)
#     implicit none
#     real, intent(inout) :: q
#     real,intent(in) :: T
#     integer, intent(inout) :: pol
#     !Energies in eV, T in K.
#     !E(omega) depends on the polarization, pol
#     bosefuncq = 1.0d0/(exp(Eomega(q,pol)/(KbeV*T))-1.0d0)
#     !Bose-Einstein dist.
#   end function bosefuncq
#
#   real function Eomega(q,pol)
#     implicit none
#     real, intent(in) :: q
#     integer, intent(in) :: pol
#     if (pol == 0) then
#        Eomega = omega0
#     else
#        Eomega = hbarev*ssound(pol)*q
#     end if
#     !See Ashcroft and mermin (22.65)
#     !Here i'm assuming the speed of sound, ssound is isotropic - which
#     !it pretty much is.
#     !Polarization = 0 is for LO
#     !             = 1 is for LA
#     !             = 2 is for TA
#     !Eomega is really E(omega) as it is in units of eV
#   end function Eomega
#
#   real function Omega(q,pol)
#     implicit none
#     real, intent(in) :: q
#     integer, intent(in) :: pol
#     if (pol == 0) then
#        Omega = omega0/hbarev
#     else
#        Omega = ssound(pol)*q
#     end if
#     !See Ashcroft and mermin (22.65)
#     !Here i'm assuming the speed of sound, ssound is isotropic - which
#     !it pretty much is.
#     !Polarization = 0 is for LO
#     !             = 1 is for LA
#     !             = 2 is for TA
#     !Eomega is really E(omega) as it is in units of eV
#   end function Omega
#
#   real function Tempf(i)
#     implicit none
#     integer, intent(in) :: i
#     Tempf = 300.0d0*real(i-1)/real(NoTemps-1)
#   end function Tempf
#
#   real function specE(i)
#     implicit none
#     integer, intent(in) :: i
#     specE = (real(i)*(maxenergy)/real(specsize))
#   end function specE
#
#   integer function speci(E)
#     implicit none
#     real, intent(in) :: E
#     speci = nint(E*real(specsize)/maxenergy)
#   end function speci
#
#   real function wavelength(E)
#     implicit none
#     real, intent(in) :: E
#     wavelength = h*cvacuum*1E9_8/(E*eV)
#   end function wavelength
#
#   real function onedgauss(x,sigma,mu)
#     implicit none
#     real, intent(in) :: x,sigma,mu
#     onedgauss=(exp(-0.5d0*(((x-mu)**2)/(sigma**2))))/(sigma*sqrt(2.0d0*pi))
#   end function onedgauss
#
#   subroutine localint(Arr,sz,ag)
#     implicit none
#     integer, intent(in) :: sz(3), ag(3)
#     real, intent(in),dimension(0:sz(1)-1,0:sz(2)-1,0:sz(3)-1):: Arr
#     integer, dimension(3) :: bigloc
#     real, dimension(3) :: biglocr
#     integer, dimension(6) :: intarea,arrsize
#     bigloc = (lbound(Arr**2) + maxloc(Arr**2) - 1)
#     biglocr(1) = bigloc(1)*grid(ag,1) + Rsize(1)
#     biglocr(2) = bigloc(2)*grid(ag,1) + Rsize(2)
#     biglocr(3) = bigloc(3)*grid(ag,3) + Rsize(3)
#     intarea(1:2) = bigloc(1:2)-nint(1.0/grid(ag,1))
#     intarea(3) = bigloc(3)-nint(1.0/grid(ag,3))
#     intarea(4:5) = bigloc(1:2)+nint(1.0/grid(ag,1))
#     intarea(6) = bigloc(3)+nint(1.0/grid(ag,3))
#     arrsize(1:3) = 0
#     arrsize(4:6) = sz-1
#     write(6,*) intarea
#     write(6,*) arrsize
#
# !       write(6,*) 'Length of eigenvector: ', sqrt(sum(temp**2))
#     write(6,*) 'Max value of psi: ', maxval(Arr)
#     write(6,*) 'Location of max value: ', biglocr
#     write(6,*) 'Max value of psi^2: ', maxval(Arr**2)
#     write(6,*) 'Prob of finding particle within 1nm of max psi: ',&
#          & perint(Arr**2,intarea,arrsize,ag)
#     return
#   end subroutine localint
#
#   real function perint(Arr,limits,asize,ag)
#     implicit none
#     integer, intent(in) :: ag(3)
#     integer,dimension(6), intent(in) :: limits, asize
#     real,intent(in)::Arr(asize(1):asize(4),asize(2):asize(5),asize(3):asize(6))
#     real,dimension(3) :: bsize
#     real :: sum
#     integer :: l,m,p,ll,mm,pp, n
#     sum = 0.0d0
#     do n = 1, 3
#        bsize(n) = limits(n+3)-limits(n)
#     end do
#
#     do l = 1, bsize(1)
#        if (l + limits(1) .gt. asize(4)) then
#           ll = l - bsize(1)
#        else
#           ll = l
#        end if
#        do m = 1, bsize(2)
#           if (m + limits(2) .gt. asize(5)) then
#              mm = m - bsize(2)
#           else
#              mm = m
#           end if
#           do p = 1, bsize(3)
#              if (p + limits(3) .gt. asize(6)) then
#                 pp = p - bsize(3)
#              else
#                 pp = p
#              end if
#              sum = sum + Arr(ll,mm,pp)
#           end do
#        end do
#     end do
#     perint = sum*grid(ag,3)*(grid(ag,1)**2)
#   end function perint
#
#   subroutine pmelement(wavefnb,wavefna,pm,sz,ag)
#     !P_if = <f|p|i> = P_ba = <a|p|b>
#     !This routine could be sped up by actually working out what i+1 etx is in
#     !terms of n, I don't think its very complicated, but it would need testing.
#     implicit none
#     integer, intent(in) :: sz(3), ag(3)
#     complex, intent(in),dimension(product(sz))::wavefna
#     complex, intent(in),dimension(product(sz))::wavefnb
#     complex, dimension(3,product(sz)):: diffwavefnb
#     complex, intent(out), dimension(3) :: pm
#     integer :: i,j,k,n,p1,p2,p3,m1,m2,m3
#     integer, dimension(3) :: plus, minus
#     real :: fg(3),norm
#
#     if (nodimensions == 1) then
#        fg = 1.0d0
#        fg(direction) = grid(ag,direction)*1E-9
#        norm = 1E-9
#     else
#        do i = 1, 3
#           fg(i) = grid(ag,i)*1E-9
#           !Coordinates in m
#        end do
#        norm = 1E-27
#     end if
#
#     do k = 0, sz(3)-1
#        do j = 0, sz(2)-1
#           do i = 0, sz(1)-1
#              plus(1) = mod(i+1,sz(1))
#              plus(2) = mod(j+1,sz(2))
#              plus(3) = mod(k+1,sz(3))
#              minus(1) = mod(i-1+sz(1),sz(1))
#              minus(2) = mod(j-1+sz(2),sz(2))
#              minus(3) = mod(k-1+sz(3),sz(3))
#              call getn(m1,minus(1),j,k,sz)
#              call getn(m2,i,minus(2),k,sz)
#              call getn(m3,i,j,minus(3),sz)
#              call getn(n,i,j,k,sz)
#              call getn(p1,plus(1),j,k,sz)
#              call getn(p2,i,plus(2),k,sz)
#              call getn(p3,i,j,plus(3),sz)
#              diffwavefnb(1,n)=(wavefnb(p1)-wavefnb(m1))
#              diffwavefnb(2,n)=(wavefnb(p2)-wavefnb(m2))
#              diffwavefnb(3,n)=(wavefnb(p3)-wavefnb(m3))
#           end do
#        end do
#     end do
#     do n = 1, 3
#        pm(n) = cmplx(0.0d0,-1.0d0)*hbar*&
#             &dot_product(wavefna/sqrt(norm),&
#             &diffwavefnb(n,:)/(sqrt(norm)*2.0d0*fg(n)))*product(fg)
#     end do
# !    write(6,*) 'p: ', pm(adim(3))
# !    write(6,*) 'p^2: ', psqm(3)
# !    write(6,*) 'd(p): ',sqrt(psqm-pm**2)
#   end subroutine pmelement
#
#   subroutine psqmelement(wavefninb,wavefnina,pm,psqm,stddevp,sz,ag)
#     !P_if = <f|p|i> = P_ba = <a|p|b>
#     !This routine could be sped up by actually working out what i+1 etx is in
#     !terms of n, I don't think its very complicated, but it would need testing.
#     implicit none
#     integer, intent(in) :: sz(3), ag(3)
#     complex, intent(in),dimension(product(sz))::wavefnina
#     complex, intent(in),dimension(product(sz))::wavefninb
#     complex,dimension(product(sz))::wavefna, wavefnb
#     real, dimension(3), intent(out) :: psqm,stddevp
#     complex, dimension(3,product(sz)):: diffwavefnb,scnddiffwavefnb
#     complex, intent(out), dimension(3) :: pm
#     integer :: i,j,k,n,p1,p2,p3,m1,m2,m3
#     integer, dimension(3) :: plus, minus
#     real :: fg(3),norm
#
#     if (nodimensions == 1) then
#        fg = 1.0d0
#        fg(direction) = grid(ag,direction)*1E-9
#        norm = 1E-9
#     else
#        do i = 1, 3
#           fg(i) = grid(ag,i)*1E-9
#           !Coordinates in m
#        end do
#        norm = 1E-27
#     end if
#
#     wavefna=wavefnina/sqrt(norm)
#     wavefnb=wavefninb/sqrt(norm)
#
#     do k = 0, sz(3)-1
#        do j = 0, sz(2)-1
#           do i = 0, sz(1)-1
#              plus(1) = mod(i+1,sz(1))
#              plus(2) = mod(j+1,sz(2))
#              plus(3) = mod(k+1,sz(3))
#              minus(1) = mod(i-1+sz(1),sz(1))
#              minus(2) = mod(j-1+sz(2),sz(2))
#              minus(3) = mod(k-1+sz(3),sz(3))
#              call getn(n,i,j,k,sz)
#              call getn(p1,plus(1),j,k,sz)
#              call getn(p2,i,plus(2),k,sz)
#              call getn(p3,i,j,plus(3),sz)
#              call getn(m1,minus(1),j,k,sz)
#              call getn(m2,i,minus(2),k,sz)
#              call getn(m3,i,j,minus(3),sz)
#              scnddiffwavefnb(1,n)=wavefnb(p1)-(2.0d0*wavefnb(n))+wavefnb(m1)
#              scnddiffwavefnb(2,n)=wavefnb(p2)-(2.0d0*wavefnb(n))+wavefnb(m2)
#              scnddiffwavefnb(3,n)=wavefnb(p3)-(2.0d0*wavefnb(n))+wavefnb(m3)
#              diffwavefnb(1,n)=(wavefnb(p1)-wavefnb(m1))/(2.0d0*fg(1))
#              diffwavefnb(2,n)=(wavefnb(p2)-wavefnb(m2))/(2.0d0*fg(2))
#              diffwavefnb(3,n)=(wavefnb(p3)-wavefnb(m3))/(2.0d0*fg(3))
#           end do
#        end do
#     end do
#     do n = 1, 3
#        psqm(n) = -(hbar**2)*&
#             &dot_product(wavefna,scnddiffwavefnb(n,:)/(fg(n)**2))*product(fg)
#        pm(n) = cmplx(0.0d0,-1.0d0)*hbar*&
#             &dot_product(wavefna,diffwavefnb(n,:))*product(fg)
#     end do
#     stddevp = real(sqrt(psqm-pm**2))
# !    write(6,*) 'p: ', pm(adim(3))
# !    write(6,*) 'p^2: ', psqm(3)
# !    write(6,*) 'd(p): ',sqrt(psqm-pm**2)
#   end subroutine psqmelement
#
#   subroutine expectedk(psi2in,psi1in,expk,sz,ag)
#     !k_if = <f|k|i> = k_21 = <1|k|2>
#     implicit none
#     complex, intent(out) :: expk(3)
#     integer, intent(in) :: sz(3), ag(3)
#     complex, intent(in) :: psi1in(product(sz)), psi2in(product(sz))
#     complex :: psi1(product(sz)), psi2(product(sz))
#     real :: q(3),qgrid(3),fg(3),norm
#     complex :: sumation(3)
#     integer :: l,m,p,n(3),i
#     sumation = 0.0d0
#
#     if (nodimensions == 1) then
#        qgrid = 1.0d0
#        fg = 1.0d0
#        qgrid(direction) = (2.0d0*pi)/(sz(direction)*grid(ag,direction)*1E-9)
#        fg(direction) = grid(ag,direction)*1E-9
#        norm = 1E-9
#     else
#        do i = 1, 3
#           qgrid(i) = (2.0d0*pi)/(sz(i)*grid(ag,i)*1E-9)
#           !Coordinates in m-1
#           fg(i) = grid(ag,i)*1E-9
#           !Coordinates in m
#           norm = 1E-27
#        end do
#     end if
#
#     psi1=psi1in/sqrt(norm)
#     psi2=psi2in/sqrt(norm)
#
#     call FFT(psi1,sz,1)
#     call FFT(psi2,sz,1)
#     !Call FFT
#
#     do p = 0, sz(3)-1
#        do m = 0, sz(2)-1
#           do l = 0, sz(1)-1
#              n(1) = (mod((l + (sz(1)/2)),sz(1)) - sz(1)/2)
#              n(2) = (mod((m + (sz(2)/2)),sz(2)) - sz(2)/2)
#              n(3) = (mod((p + (sz(3)/2)),sz(3)) - sz(3)/2)
#              do i = 1, 3
#                 q(i) = (qgrid(i)*n(i))
#              end do
#              sumation = sumation + (conjg(psi1(nvalue(l,m,p,sz)))*q*&
#                   &psi2(nvalue(l,m,p,sz)))
#           end do
#        end do
#     end do
#
#     expk = sumation*product(fg)/product(sz)
# !    write(6,*) 'k: ', expk(adim(3))*hbar
#
#   end subroutine expectedk
#
#   subroutine expectedksq(psi2in,psi1in,expk,expksq,stddevk,sz,ag)
#     !k_if = <f|k|i> = k_21 = <1|k|2>
#     implicit none
#     complex, intent(out) :: expk(3),expksq(3),stddevk(3)
#     integer, intent(in) :: sz(3), ag(3)
#     complex, intent(in) :: psi1in(product(sz)), psi2in(product(sz))
#     complex :: phi1(product(sz)), phi2(product(sz))
#     real :: qsq(3),qgrid(3),fg(3),norm,q(3)
#     complex :: sumation(3),normal,sumationsq(3)
#     integer :: l,m,p,n(3),i
#     sumation = 0.0d0
#     sumationsq = 0.0d0
#
#     if (nodimensions == 1) then
#        qgrid = 1.0d0
#        fg = 1.0d0
#        qgrid(direction) = (2.0d0*pi)/(sz(direction)*grid(ag,direction)*1E-9)
#        fg(direction) = grid(ag,direction)*1E-9
#        norm = 1E-9
#     else
#        do i = 1, 3
#           qgrid(i) = (2.0d0*pi)/(sz(i)*grid(ag,i)*1E-9)
#           !Coordinates in m-1
#           fg(i) = grid(ag,i)*1E-9
#           !Coordinates in m
#        end do
#        norm = 1E-27
#     end if
#
#     phi1=psi1in/sqrt(norm)
#     phi2=psi2in/sqrt(norm)
#
#     call FFT(phi1,sz,1)
#     call FFT(phi2,sz,1)
#     !Call FFT
#
#     do p = 0, sz(3)-1
#        do m = 0, sz(2)-1
#           do l = 0, sz(1)-1
#              n(1) = (mod((l + (sz(1)/2)),sz(1)) - sz(1)/2)
#              n(2) = (mod((m + (sz(2)/2)),sz(2)) - sz(2)/2)
#              n(3) = (mod((p + (sz(3)/2)),sz(3)) - sz(3)/2)
#              do i = 1, 3
#                 qsq(i) = (qgrid(i)*n(i))**2
#                 q(i) = qgrid(i)*n(i)
#              end do
#              sumationsq = sumationsq + (conjg(phi1(nvalue(l,m,p,sz)))*qsq*&
#                   &phi2(nvalue(l,m,p,sz)))
#              sumation = sumation + (conjg(phi1(nvalue(l,m,p,sz)))*q*&
#                   &phi2(nvalue(l,m,p,sz)))
#           end do
#        end do
#     end do
#
#     expk = sumation*product(fg)/product(sz)
#     expksq = sumationsq*product(fg)/product(sz)
#     stddevk = sqrt(expksq-expk**2)
# !    write(6,*) 'k: ', expk(adim(3))*(hbar)
# !    write(6,*) 'k^2: ', expksq(3)*(hbar**2)
# !    write(6,*) 'Delta(k): ', sqrt(expksq-expk**2)
#   end subroutine expectedksq
#
#   subroutine expectedpos(Arr,r,stddevr,sz,ag)
#     implicit none
#     integer, intent(in) :: sz(3), ag(3)
#     complex, intent(in),dimension(product(sz)):: Arr
#     real, intent(out), dimension(3) :: r, stddevr
#     real, dimension(3) :: expr,s,sqvarexpr,z,exprsq
#     real,dimension(product(sz)):: Arrsq
#     integer, dimension(3) :: zeropos
#     integer, dimension(1) :: zerop
#     integer :: i,j,k,n,m
#     real :: fg(3),ssize(3)
#     do i = 1, 3
#        fg(i) = grid(ag,i)
#     end do
#     ssize = fg*sz
#     expr = 0.0d0
#     exprsq = 0.0d0
#     sqvarexpr = 0.0d0
#     Arrsq = (abs(Arr))**2
#     zerop = minloc(Arrsq)
#     call getijk(zerop(1),zeropos(1),zeropos(2),zeropos(3),sz)
# !    z(1) = real(zeropos(1))*a0
# !    z(2) = real(zeropos(2))*a0
# !    z(3) = real(zeropos(3))*c0
#     z = real(zeropos)*fg
#     write(17,*) 'position of minima: ', z
#     do k = 0, sz(3)-1
#        do j = 0, sz(2)-1
#           do i = 0, sz(1)-1
#              s(1) = (real(i)*fg(1)) - z(1)
#              s(2) = (real(j)*fg(2)) - z(2)
#              s(3) = (real(k)*fg(3)) - z(3)
#              do n = 1, 3
#                 if (s(n) .lt. 0.0) then
#                    s(n) = s(n) + ssize(n)
#                 end if
#              end do
#              expr = expr + s*Arrsq(nvalue(i,j,k,sz))
#              exprsq = exprsq + (s**2)*Arrsq(nvalue(i,j,k,sz))
#           end do
#        end do
#     end do
#     expr = expr*product(fg)
#     exprsq = exprsq*product(fg)
#     do m = 1, 3
#        if ((expr(m) + z(m)) .gt. ssize(m)) then
#           r(m) = expr(m) + z(m) - ssize(m)
#        else
#           r(m) =  expr(m) + z(m)
#        end if
#     end do
#     write(17,*) 'r: ',r
#     stddevr = sqrt(exprsq-expr**2)
#     return
#   end subroutine expectedpos
#
#   subroutine outputwave(Arr,unit,dimn,slice,off,sz,ag,val)
#     implicit none
#     integer, intent(in) :: sz(3), ag(3)
#     real, intent(in), dimension(0:sz(1)-1,0:sz(2)-1,0:sz(3)-1) :: Arr
#     integer, intent(in) :: unit,dimn,slice
#     real, intent(in) :: off
#     real,dimension(:,:), allocatable :: surfacepos
#     real, intent(in), optional :: val
#     !    dim = 1 > xslice, 2 > yslice, 3 > zslice
#     !val1>val2>val3
#     integer :: l,m,p,label,surfacearea,count,i,j
#     real :: x1,x2,x3,norm
# 20  format(9G24.15E3)
# 21  format(3G24.15E3)
#     norm = 1.0d0/product(Rsize)
# !    surfacearea = product(sz)
#     surfacearea = 2*((sz(1)*sz(2))+(sz(1)*sz(3))+(sz(2)*sz(3)))
#     allocate(surfacepos(3,surfacearea))
#
#     if (present(val) .and. (dimn .ne. 4)) write(6,*) &
#          &'ERROR(1) in calling array output'
#     if (dimn == 1) write(unit,*) '#    y/nm     z/nm        Psi...'
#     if (dimn == 2) write(unit,*) '#    x/nm     z/nm        Psi...'
#     if (dimn == 3) write(unit,*) '#    x/nm     y/nm        Psi...'
#     if (dimn == 4) write(unit,*) '#    x/nm     y/nm        z/nm'
#     if (dimn .lt. 4) then
#        if (dimn == 1) then
#           l = slice
#           do p = 0, sz(3)-1
#              do m = 0, sz(2)-1
#                 x2 = m*grid(ag,1)
#                 x3 = p*grid(ag,3)
#                 write(unit,20) x2, x3, Arr(l,m,p), sum(Arr(:,m,p))/sz(1), &
#                      & Arr(l,m,p)+off,(Arr(l,m,p))**2,(Arr(l,m,p))**2+off, &
#                      & (sum(Arr(:,m,p))/sz(1))+off,abs(Arr(l,m,p))
#              end do
#           end do
#        else if (dimn == 2) then
#           m = slice
#           do p = 0, sz(3)-1
#              do l = 0, sz(1)-1
#                 x1 = l*grid(ag,1)
#                 x3 = p*grid(ag,3)
#                 write(unit,20) x1, x3, Arr(l,m,p), sum(Arr(l,:,p))/sz(2), &
#                      & Arr(l,m,p)+off,(Arr(l,m,p))**2,(Arr(l,m,p))**2+off, &
#                      & (sum(Arr(l,:,p))/sz(2))+off,abs(Arr(l,m,p))
#              end do
#           end do
#        else if (dimn == 3) then
#           p = slice
#           do m = 0, sz(2)-1
#              do l = 0, sz(1)-1
#                 x1 = l*grid(ag,1)
#                 x2 = m*grid(ag,1)
#                 write(unit,20) x1, x2, Arr(l,m,p), sum(Arr(l,m,:))/sz(3), &
#                      & Arr(l,m,p)+off,(Arr(l,m,p))**2,(Arr(l,m,p))**2+off, &
#                      & (sum(Arr(l,m,:))/sz(3))+off,abs(Arr(l,m,p))
#              end do
#           end do
#        end if
#     else
#        call isosurface(Arr**2,val,surfacepos,surfacearea,count,sz,ag)
# !       write(6,*) 'No of isosurface points for psi^2 =',val,':',count
# !       write(6,*) 'Ratio of Total SA to psi SA:',real(count)/&
# !            &real(surfacearea)
#        do j = 1, count
#           write(unit,21) surfacepos(:,j)
#        end do
#
#     end if
#     return
#   end subroutine outputwave