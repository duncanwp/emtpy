__author__ = 'pard'

#
# program main
#   use routines
# !  use avroutines
#   use makesample
#
#   implicit none
#   real :: tol, t1, t2, totaltime
#   integer :: hours,minutes,seconds
#   real, dimension(3,maxpos) :: AllR
#   integer :: ifail, enev,maxitr,hnev, sample, number,i,j,k,ave
#   integer, dimension(3) :: eksz, hksz,eag,hag
#   integer, dimension(100) :: funits,runits
#   character(10) :: time
#   character(8) :: date
#   character(40) :: filename, name
#   logical :: opened
#   call mkl_set_num_threads(4)
#   AllR = 0.0d0
#   call sampleinit(tol,ifail,enev,hnev,maxitr,sample,AllR,funits,ave,&
#        &eksz,hksz,eag,hag)
#
#   if (ave .ne. 1) then
#      write(6,*) "Defining arrays and variables completed,",totA,&
#           &" In positions read."
#
#      call cpu_time(t1)
#
#      write(rundir,'(A17,A3,I3.3,A1)') rundirroot,'run',run,'/'
#
#      call fileopenr(runits)
#
#      call date_and_time(date,time)
#      write(17,*) 'Date: ',date,',time: ', time,', run: ',run
#
#      call maincalc(sample,AllR,tol,enev,hnev,maxitr,eksz,hksz,eag,hag)
#
#      call cpu_time(t2)
#      write(17,*) 'Run, ',run+1,' completed, time taken was: ',t2-t1, 'seconds.'
#      write(6,*) 'Run, ',run+1,' completed, time taken was: ',t2-t1, 'seconds.'
#      call closefiles(runits)
#   else
#      call reconstructsp(enev,hnev,sample)
#   end if
#
#   call closefiles(funits)
#
# contains
#
#   subroutine bandcalc(echi,hchi,vPot,cPot,tol,enev,hnev,&
#        &maxit,eksize,hksize,eag,hag,hnevb)
#     implicit none
#     integer, intent(in),dimension(3) :: eksize, hksize, eag, hag
#     integer, intent(in) :: hnevb
#     complex,intent(in),dimension(0:eksize(1)-1,0:eksize(2)-1,0:eksize(3)-1):: &
#          &echi, cPot
#     complex,intent(in),dimension(holecalcs,0:hksize(1)-1,0:hksize(2)-1,0:hksize(3)-1):: &
#          &hchi
#     complex,intent(in) :: &
#          &vPot(holecalcs,band+1,0:hksize(1)-1,0:hksize(2)-1,0:hksize(3)-1)
#     integer,intent(inout) :: enev,hnev,maxit
#     real, intent(inout) :: tol
# !    complex,allocatable,dimension(:,:) :: evecsmth, tothvecbig
#     real, dimension(enev,2) :: eeigenval
#     real, dimension(hnev,2) :: heigenval
#     real, dimension(hnevb,2) :: tothval
#     real, dimension(enev,6) :: elocal
#     real, dimension(hnev,6) :: hlocal
#     complex, dimension(product(hksize),hnevb) :: tothvec
#     complex, dimension(product(eksize),enev) :: eeigenvec
#     complex, dimension(product(hksize),hnev) :: heigenvec
#     real :: eX(product(eksize)),hX(holecalcs,product(hksize))
#     integer :: i, locale,enav,j,k,holeno,conf,hconfined!,newksz(3),nwsz
#     integer :: calc,funits(100),off(3),hnav
#     character(2) :: param
#
#     locale = 0
#     holeno = 0
#     off = 0
# !!$    newksz(3) = hksize(3)
# !!$    newksz(1:2) = hksize(1:2)*ehsratio
# !!$    nwsz = product(newksz)
#
#     if (enev .lt. 4) then
#        enav = 20
#     else
#        enav = enev*2
#     end if
#
#     if (hnev .lt. 4) then
#        hnav = 20
#     else
#        hnav = hnev*2
#     end if
#
#     do k = 0, eksize(3)-1
#        do j = 0, eksize(2)-1
#           do i = 0, eksize(1)-1
#              eX(nvalue(i,j,k,eksize)) = real(echi(i,j,k))
#           end do
#        end do
#     end do
#
#     do calc = 1, holecalcs
#        do k = 0, hksize(3)-1
#           do j = 0, hksize(2)-1
#              do i = 0, hksize(1)-1
#                 hX(calc,nvalue(i,j,k,hksize)) = real(hchi(calc,i,j,k))
#              end do
#           end do
#        end do
#     end do
#
#     if (enev .gt. 0) then
#        param = '_e'
#
#        allocate(ija((product(eksize)*4)+1))
#        allocate(sa((product(eksize)*4)+1))
#        ija = 0
#        sa = 0.0d0
#
#        call sprsymqckprdc(cPot,echi,param,eksize,eag)
#        write(17,*) "Electron Hamiltonian constructed"
#        write(6,*) "Electron Hamiltonian constructed"
# !    call wavefncmplx(enev,enav,tol,maxit,eeigenvec,eeigenval,conf,param,&
# !         &product(eksize),eag,1)
#        call wavefn(enev,enav,tol,maxit,eeigenvec,eeigenval,conf,param,&
#             &product(eksize),eag,1)
#        write(17,*) "Electron Wavefunctions calculated"
#        write(6,*) "Electron Wavefunctions calculated"
#        call eigenout(eeigenval(:conf,:),eeigenvec(:,:conf),elocal,conf,param,&
#             &eX,eksize,eag,off)
#        call melements(eeigenval(:conf,:),eeigenvec(:,:conf),&
#             &eeigenval(:conf,:),eeigenvec(:,:conf),conf,conf,46,eksize,eag)
#        write(17,*) "Relevant Electron data and matrix elements output"
#        write(6,*) "Relevant Electron data and matrix elements output"
#
#        deallocate(ija)
#        deallocate(sa)
#
#        do i = 1, conf
#           if (localised(elocal(i,4:6),1)) then
#              write(6,*) 'Electron state ', i,'is localised.'
#              write(6,*) elocal(i,4:6)
#              write(6,*) fraclocal*maxsigma(1,:)
#              locale = locale + 1
#           end if
#        end do
#        write(6,*) locale,'localised electrons'
#     end if
#
#     if (locale .gt. 0 .or. enev .eq. 0) then
#
#        do calc = 1, holecalcs
#           holeno = 0
# !!!!!!
# !!!!!!This line resets the total no. of holes every run, for outputting each
# !!!!!!hole calculation seperately
#           allocate(ija((product(hksize)*4)+1))
#           allocate(sa((product(hksize)*4)+1))
#           call fileopenh(funits,calc)
#           write(6,*) 'Calculating hole section: ', calc
#           off(1)=nint((rsize(1)/real(ehsratio))*&
#                &mod(calc-1,ehsratio)/grid(hag,1))
#           off(2)=nint((rsize(2)/real(ehsratio))*&
#                &(holecalcoffy(calc)-1)/grid(hag,2))
#           off(3) = 0
#           write(6,*) 'Offset: ', off
#           if (band .eq. 0) then
#              ija = 0
#              sa = 0.0d0
#              param = '_h'
#              call sprsymqckprdc(-vPot(calc,band+1,:,:,:),hchi(calc,:,:,:),&
#                   &param,hksize,hag)
#              write(17,*) "Hole Hamiltonian constructed"
#              write(6,*) "Hole Hamiltonian constructed"
# !             call wavefncmplx(hnev,2*hnev,tol,maxit,heigenvec,heigenval&
# !                  &,hconfined,param,product(hksize),hag,calc)
#              call wavefn(hnev,hnav,tol,maxit,heigenvec,heigenval&
#                   &,hconfined,param,product(hksize),hag,calc)
#              write(17,*) "Hole Wavefunctions calculated"
#              write(6,*) "Hole Wavefunctions calculated"
#              call eigenout(heigenval(:hconfined,:),heigenvec(:,:hconfined),&
#                   &hlocal,hconfined,param,hX(calc,:),hksize,hag,off)
#              tothvec(:,:hconfined) = heigenvec(:,:hconfined)
#              tothval(:hconfined,:) = heigenval(:hconfined,:)
#              holeno = hconfined
#           else
#              do i = 1, band
#                 ija = 0
#                 sa = 0.0d0
#                 if (i == 1) param = 'hh'
#                 if (i == 2) param = 'lh'
#                 if (i == 3) param = 'sh'
#                 call sprsymqckprdc(-vPot(calc,i+1,:,:,:),hchi(calc,:,:,:),&
#                      &param,hksize,hag)
#                 write(17,*) param, " Hamiltonian constructed"
#                 write(6,*) param, " Hamiltonian constructed"
# !                call wavefncmplx(hnev,2*hnev,tol,maxit,heigenvec,&
# !                     &heigenval,hconfined,param,product(hksize),hag,calc)
#                 call wavefn(hnev,hnav,tol,maxit,heigenvec,&
#                      &heigenval,hconfined,param,product(hksize),hag,calc)
#                 write(17,*) param, " Wavefunctions calculated"
#                 write(6,*) param, " Wavefunctions calculated"
#                 call eigenout(heigenval(:hconfined,:),heigenvec(:,:hconfined),&
#                      &hlocal,hconfined,param,hX(calc,:),hksize,hag,off)
#                 write(17,*) "Relevant Hole data output"
#                 write(6,*) "Relevant Hole data output"
#                 tothvec(:,holeno+1:holeno+hconfined) = heigenvec(:,:hconfined)
#                 tothval(holeno+1:holeno+hconfined,:) = heigenval(:hconfined,:)
#                 holeno = holeno + hconfined
#                 hconfined = 0
#              end do
#           end if
#           call melements(tothval(:holeno,:),tothvec(:,:holeno),&
#                &tothval(:holeno,:),tothvec(:,:holeno),holeno,holeno,47,&
#                &hksize,hag)
#           write(17,*) "Relevant Hole matrix elements output"
#           write(6,*) "Relevant Hole matrix elements output"
#
#           deallocate(ija)
#           deallocate(sa)
#
# !!$          allocate(evecsmth(product(hksize)*holecalcs,conf))
# !!$          allocate(tothvecbig(product(hksize)*holecalcs,holeno))
# !!$
# !!$          do i = 1, conf
# !!$          call threedinterpolate(eeigenvec(:,i),eksize,eag,evecsmth(:,i),&
# !!$                  &newksz,hag)
# !!$          end do
# !!$
# !!$          do i = 1, holeno
# !!$          call lineuparr(tothvec(:,i),hksize,hag,tothvecbig(:,i),newksz,off)
# !!$          end do
#
# !          call melements(eeigenval(:conf,:),evecsmth(:,:conf),&
# !               &tothval(:holeno,:),tothvecbig(:,:holeno),conf,holeno,&
# !               &20,newksz,hag)
#           call melementsx(eeigenval(:conf,:),eeigenvec(:,:conf),eksize,eag,&
#                &conf,tothval(:holeno,:),tothvec(:,:holeno),hksize,hag,holeno,&
#                &20,off)
#           write(17,*) "Relevant Electron/Hole matrix elements output"
#           write(6,*) "Relevant Electron/Hole matrix elements output"
# !!$          deallocate(evecsmth)
# !!$          deallocate(tothvecbig)
#           call closefiles(funits)
#        end do
#     else
#        do i = 1, conf
#           write(6,*) 'Electron state ', i,'is not localised.'
#           write(6,*) elocal(i,4:6)
#           write(6,*) fraclocal*maxsigma(1,:)
#        end do
#        write(6,*) 'NO LOCALISED ELECTRON STATE IN THIS SYSTEM'
#     end if
#
#   end subroutine bandcalc
#
#   subroutine maincalc(sample,AllR,tol,enev,hnev,maxit,eksize,hksize,eag,hag)
#     implicit none
#     integer,intent(inout) :: sample,enev,hnev,maxit
#     real, intent(inout) :: tol
#     real, intent(inout),dimension(3,maxpos) :: AllR
#     integer, intent(inout), dimension(3) :: eksize, hksize, eag, hag
#     complex :: vV(holecalcs,band+1,0:hksize(1)-1,0:hksize(2)-1,0:hksize(3)-1)
#     complex,dimension(0:eksize(1)-1,0:eksize(2)-1,0:eksize(3)-1):: eConfine,&
#          &epiezo, echi,cV, edef,eEffm,espont,wellprofe
#     complex :: hchi(holecalcs,0:hksize(1)-1,0:hksize(2)-1,0:hksize(3)-1)
#     complex,dimension(0:hksize(1)-1,0:hksize(2)-1,0:hksize(3)-1):: hConfine,&
#          &hpiezo, hdef,hEffm,hspont,wellprofh
#     complex :: eU(3,3,0:eksize(1)-1,0:eksize(2)-1,0:eksize(3)-1)
#     complex :: hU(3,3,0:hksize(1)-1,0:hksize(2)-1,0:hksize(3)-1)
#     logical :: wellmaskh(0:hksize(1)-1,0:hksize(2)-1,0:hksize(3)-1)
#     logical :: wellmaske(0:eksize(1)-1,0:eksize(2)-1,0:eksize(3)-1)
#     integer :: i,j,k,slice,ifail,Ano,Totn,eatoms,hnevb, ii, jj
#     integer,dimension(100) :: units
#     character(40) :: filename, name
#     character(2) :: param
#     real :: coff, voff,SelR(3,maxpos),hrsizerun(6)
#     real :: avgape, avdefe, avcone, avchie
#     real :: avgaph, avdefh, avconh, avchih
#
#     Totn = 0
#
#     call getchi(sample,AllR,echi,totA,eksize,eag)
#
#     write(6,*) 'maxval of chi:',maxval(real(echi))
#
#     call fileopen(units)
#
#     if (direction.eq.3) slice = 1
#     if (direction .eq. 1) slice = 2
#     if (direction .eq. 2) slice = 3
#
#     eConfine = 0.0d0
#     eEffM = 0.0d0
#     edef = 0.0d0
#     eU = 0.0d0
#     epiezo = 0.0d0
#     espont = 0.0d0
#     cV = 0.0d0
#     Vv = 0.0d0
#
#     call FFT3d(echi,eksize,1)
#     write(17,*) "Molar concentration calculated and fourier transformed"
#     write(6,*) "Molar concentration calculated and fourier transformed"
#
#     call straintensor(echi,eU,eksize,eag)
#
#     call spontpot(echi,espont,eksize,eag)
#     call piezopot(echi,epiezo,eksize,eag)
#     write(17,*) "Piezoelectric and polarization potentials calculated in &
#          &fourier space"
#     write(6,*) "Piezoelectric and polarization potentials calculated in &
#             &fourier space"
#
#     do j = 1, 3
#        do i= 1, 3
#           write(filename,'(A24,A6,2I1.1,A4)' ) rundir,'strain',i,j,'.dat'
#           open(27,file=filename,status='Unknown')
#           call FFT3d(eU(i,j,:,:,:),eksize,-1)
#           call outputArr(eU(i,j,:,:,:),27,slice,eksize(slice)/2,eksize,eag)
#           close(27)
#        end do
#     end do
#
#     call FFT3d(espont,eksize,-1)
#     call FFT3d(epiezo,eksize,-1)
#     write(17,*) "Piezoelectric and polarization potentials inverse fourier &
#          &transformed"
#     write(6,*) "Piezoelectric and polarization potentials inverse fourier &
#          &transformed"
#
#     call FFT3d(echi,eksize,-1)
#     write(17,*) "Molar concentration inverse fourier transformed"
#     write(6,*) "Molar concentration inverse fourier transformed"
#
#     call outputArr(espont,11,slice,eksize(slice)/2,eksize,eag)
#     call outputArr(epiezo,37,slice,eksize(slice)/2,eksize,eag)
#     write(17,*) 'Piezoelectric drop accross the well: ',&
#          &maxval(real(epiezo))-minval(real(epiezo))
#     write(17,*) 'Psp drop accross the well: ', maxval(real(espont)) - &
#          &minval(real(espont))
#     cV = espont + epiezo
#
#     write(17,*) 'Combined electric field drop accross the well: ',&
#          & maxval(real(cV))-minval(real(cV))
#     write(17,*) 'or, ', 10.0d0*(maxval(real(cV))-minval(real(cV)))/(qwb-qwa),&
#          &'MV/cm'
#     write(17,*) 'or, ', 1E3*(maxval(real(cV))-minval(real(cV)))/(qwb-qwa),&
#          &'MV/m'
#
#     call offset(echi,eConfine,'_e',eksize,eag)
#     call defpot(eU,echi,edef,'_e',eksize,eag)
#     cV = eConfine - cV + edef
#     maxCcon = maxval(real(eConfine+edef))
#     minCcon = minval(real(eConfine+edef))
# !    coff = maxCcon - minCcon
#
#     call sqwell(wellprofe,qwa,qwb,nom,background,eksize,eag)
#     wellmaske = real(wellprofe) .gt. 0.0
#
#     call sqwell(wellprofh,qwa,qwb,nom,background,hksize,hag)
#     wellmaskh = real(wellprofh) .gt. 0.0
#
#     avgape = abs(sum(eConfine,wellmaske)/real(count(wellmaske))-eConfine(0,0,0))
#     avdefe = abs(sum(edef,wellmaske)/real(count(wellmaske))-edef(0,0,0))
#     avcone = abs(sum(eConfine+edef,wellmaske)/real(count(wellmaske))-&
#          &eConfine(0,0,0))
#     avchie = sum(echi,wellmaske)/real(count(wellmaske))
#     write(17,*) "Av. chi in well: ", avchie
#     write(17,*) "average gap        average deformation     average confinement"
#     write(17,*) avgape,avdefe,avcone
#
#     coff = avcone
#
#     call outputArr(eConfine,21,slice,eksize(slice)/2,eksize,eag)
#     call outputArr(edef,28,slice,eksize(slice)/2,eksize,eag)
#     call outputArr(cV,10,direction,eksize(direction)/2,eksize,eag)
#     call outputArr(cV,12,slice,eksize(slice)/2,eksize,eag)
#     call outputArr(echi,8,slice,eksize(slice)/2,eksize,eag)
# !    call outputArr(vV(1,:,:,:),14,slice,ksize(slice)/2,eksize,eag)
#     eU(1,1,:,:,:) = eU(1,1,:,:,:)+eU(2,2,:,:,:)+eU(3,3,:,:,:)
#     call outputArr(eU(1,1,:,:,:),26,slice,eksize(slice)/2,eksize,eag)
#
#     write(17,*) "Total potential calculated and output"
#     write(6,*) "Total potential calculated and output"
#
#     write(6,*) 'Electron potential calculated'
#
#     do i = 1, holecalcs
#
#        hConfine = 0.0d0
#        hEffM = 0.0d0
#        hdef = 0.0d0
#        hU = 0.0d0
#        hpiezo = 0.0d0
#        espont = 0.0d0
#        SelR = 0.0d0
#
#        ii = int(real((i-1)/ehsratio)) + 1
#        jj = int(real((i-(ii-1)*ehsratio)))
#
#        hrsizerun(adim(1)) = (rsize(adim(1))/real(ehsratio))*(ii-1)
#        hrsizerun(adim(2)) = (rsize(adim(2))/real(ehsratio))*(jj-1)
#        hrsizerun(adim(3)) = 0.0d0
#        hrsizerun(adim(1)+3) = (rsize(adim(1))/real(ehsratio))*ii
#        hrsizerun(adim(2)+3) = (rsize(adim(2))/real(ehsratio))*jj
#        hrsizerun(adim(3)+3) = rsize(adim(3))
#        write(6,*) hrsizerun
#        call Apos(AllR,SelR,hrsizerun,Ano)
#        Totn = Totn + ano
#        write(6,*) 'Number of In atoms in hole sample',i,'is: ',Ano
# !     write(6,*) maxval(SelR(1,:)),maxval(SelR(2,:)),maxval(SelR(3,:))
# !     write(6,*) minval(SelR(1,:Ano)),minval(SelR(2,:Ano)),minval(SelR(3,:Ano))
#        call getchi(sample,SelR,hchi(i,:,:,:),Ano,hksize,hag)
# !       write(6,*) 'maxval of chi:',maxval(real(hchi(i,:,:,:)))
#
#        call FFT3d(hchi(i,:,:,:),hksize,1)
#        write(17,*) "Molar concentration calculated and fourier transformed"
#        write(6,*) "Molar concentration calculated and fourier transformed"
#
#        call straintensor(hchi(i,:,:,:),hU,hksize,hag)
#
#        call spontpot(hchi(i,:,:,:),hspont,hksize,hag)
#        call piezopot(hchi(i,:,:,:),hpiezo,hksize,hag)
#        write(17,*) "Piezoelectric and polarization potentials calculated in &
#             &fourier space"
#        write(6,*) "Piezoelectric and polarization potentials calculated in &
#             &fourier space"
#
#        do j = 1, 3
#           do k = 1, 3
#              write(filename,'(A24,A6,2I1.1,A4)' ) rundir,'strain',k,j,'.dat'
#              open(27,file=filename,status='Unknown')
#              call FFT3d(hU(k,j,:,:,:),hksize,-1)
#              call outputArr(hU(k,j,:,:,:),27,slice,hksize(slice)/2,hksize,hag)
#              close(27)
#           end do
#        end do
#
#        call FFT3d(hspont,hksize,-1)
#        call FFT3d(hpiezo,hksize,-1)
#        write(17,*) "Piezoelectric and polarization potentials inverse fourier &
#             &transformed"
#        write(6,*) "Piezoelectric and polarization potentials inverse fourier &
#             &transformed"
#
#        call FFT3d(hchi(i,:,:,:),hksize,-1)
#        write(17,*) "Molar concentration inverse fourier transformed"
#        write(6,*) "Molar concentration inverse fourier transformed"
#
#
#        write(17,*) 'Piezoelectric drop accross the well: ',&
#             &maxval(real(hpiezo))-minval(real(hpiezo))
#        write(17,*) 'Psp drop accross the well: ', &
#             &maxval(real(hspont)) - minval(real(hspont))
#        do j = 1, (band+1)
#           vV(i,j,:,:,:) = hspont + hpiezo
#        end do
#
#        call offset(hchi(i,:,:,:),hConfine,'_h',hksize,hag)
#        do j = 0, band
#           if (j == 0) param = '_h'
#           if (j == 1) param = 'hh'
#           if (j == 2) param = 'lh'
#           if (j == 3) param = 'sh'
#           call defpot(hU,hchi(i,:,:,:),hdef,param,hksize,hag)
#           vV(i,j+1,:,:,:) = hConfine - vV(i,j+1,:,:,:) - hdef
#           maxVcon(i,j+1) = maxval(real(hConfine-hdef))
#           minVcon(i,j+1) = minval(real(hConfine-hdef))
# !          voff = maxVcon(i,j+1) - minVcon(i,j+1)
#
#           avgaph = abs(sum(hConfine,wellmaskh)/real(count(wellmaskh))-hConfine(0,0,0))
#           avdefh = abs(sum(hdef,wellmaskh)/real(count(wellmaskh))-hdef(0,0,0))
#           avconh = abs(sum(hConfine-hdef,wellmaskh)/real(count(wellmaskh))-hConfine(0,0,0))
#           avchih = sum(hchi(i,:,:,:),wellmaskh)/real(count(wellmaskh))
#           write(17,*) "Av. chi in well: ", avchih
#           write(17,*) "average gap        average deformation     &
#                &average confinement"
#           write(17,*) avgaph,avdefh,avconh
#
#           voff = avconh
#
#           write(17,*) 'Offsets for',param
#           write(17,*) 'CB offset: ', coff
#           write(17,*) 'Or: %', 100.0d0*coff/(coff+voff)
#           write(17,*) 'VB offset: ', voff
#           write(17,*) 'Or: %', 100.0d0*voff/(coff+voff)
#        end do
#
#        if (i == 1) then
#           call outputArr(hspont,11,slice,hksize(slice)/2,hksize,hag)
#           call outputArr(hpiezo,37,slice,hksize(slice)/2,hksize,hag)
#           call outputArr(hConfine,22,slice,hksize(slice)/2,hksize,hag)
#           call outputArr(hdef,29,slice,hksize(slice)/2,hksize,hag)
#           call outputArr(hchi(i,:,:,:),38,slice,hksize(slice)/2,hksize,hag)
#           call outputArr(vV(i,1,:,:,:),14,slice,hksize(slice)/2,hksize,hag)
#        end if
#
#     end do
#
# !    write(6,*) minVcon(:,1)
#
#     write(6,*) 'Total number of In atoms in hole samples: ',Totn
#     write(6,*) 'THis should be: ', totA
#
#     write(17,*) "Max concentration: ", maxval(real(echi))
#     write(17,*) "Max value for conduction potential: ", &
#          &maxval(real(cV)), "eV"
#     write(17,*) "Max value for valence potential: ",&
#          &maxval(real(vV(:,1,:,:,:))) , "eV"
#     write(17,*) "Min value for conduction potential: ", &
#          &minval(real(cV)), "eV"
#     write(17,*) "Min value for valence potential: ", &
#          &minval(real(vV(:,1,:,:,:))), "eV"
#     write(17,*) "Energy gap: ",&
#          &minval(real(cV))-maxval(real(vV(:,1,:,:,:))) , "eV"
#
#     call closepos
#     call closefiles(units)
#
#     if (band .gt. 0) then
#        hnevb = hnev*band
#     else
#        hnevb = hnev
#     end if
#
#     call bandcalc(echi,hchi,vV,cV,tol,enev,hnev,maxitr,&
#          &eksize,hksize,eag,hag,hnevb)
#
#     return
#   end subroutine maincalc
#
#   subroutine getchi(sample,AllR,chi,noatoms,sz,ag)
#     implicit none
#     integer, intent(in) :: sample,sz(3),ag(3),noatoms
#     real,intent(inout), dimension(3,maxpos) :: AllR
#     complex,intent(inout)::chi(0:sz(1)-1,0:sz(2)-1,0:sz(3)-1)
#     real, dimension(3) :: rpos
#     real, dimension(5) :: nomarr = (/ 0.05, 0.12, 0.15, 0.19, 0.25 /)
#     real :: FWHM,edge,well,barrier
#     integer :: i,j,k,nQWs
#     integer,dimension(3) :: spacings
# !!$    if(sample .eq. 8) then
# !!$       latt(1:2) = a0
# !!$       latt(3) = c0
# !!$       call generate(8,real(Rsize(4:6)),qwb,qwa,nom,0.0,0.0,0,AllR,chi)
#      if(sample .eq. 7) then
#        if (samplever == 0) then
#           write(6,*) nomarr(run+1)
#           call sqwell(chi,qwa,qwb,nomarr(run+1),background,sz,ag)
#        else if (samplever == 4) then
#           !!!!
#           !!!Creates a 1-d well with a Cambridge profile
#           !!!!
#           FWHM = qwb-qwa!!!2.860d0
#           call profilechi(chi,FWHM,sz,ag)
#        else if (samplever == 9) then
#           !!!!
#           !!!Creates a 1-d well with ternary alloy in barriers
#           !!!!
#           well = qwa
#           barrier = qwb
#           call sqwell(chi,well,barrier,background,nom,sz,ag)
#        else
#           call sqwell(chi,qwa,qwb,nom,background,sz,ag)
#        end if
#     else if (sample .eq. 10) then
#        !!!!!!
#        !!! Creates 1-d MQWs, with the ternary alloy in the barrier
#        !!!!!!
#        well = nint(qwa/latt(direction))*latt(direction)
#        barrier = (nint(qwb/latt(direction))*latt(direction)) - &
#             &(run*latt(direction))
#        if (samplever == 2) then
#           call dblwell(chi,qwa,qwb-(run*0.25),background,nom,sz,ag)
#        else
#           call dblwellsym(chi,qwa,qwb-(run*0.25),background,nom,sz,ag)
#        end if
#     else if (nodimensions .eq. 1) then
#        call sqwell(chi,qwa,qwb,nom,background,sz,ag)
#     else
#        call conc(AllR,chi,noatoms,sz,ag)
#        write(6,*) 'Concentration profile calculated'
#     end if
#     return
#   end subroutine getchi