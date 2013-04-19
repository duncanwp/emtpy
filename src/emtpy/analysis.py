__author__ = 'pard'


#   subroutine spectra(eloc,hloc,width,sarr,evarr,hvarr,&
#        &earr,harr,rearr,frac,sample,calc)
#     implicit none
#     integer,intent(in) :: eloc,hloc,sample,calc
#     real, intent(in) :: width,frac
#     real,intent(in) :: sarr(eloc,hloc,8),earr(eloc,eloc,8),harr(hloc,hloc,8)
#     real,intent(in) :: evarr(eloc,26), hvarr(hloc,26)
#     real, dimension(8+(3*NoTemps),specsize),intent(inout) :: rearr
#     real, dimension(specsize) :: absarr,selarr
#     real, dimension(specsize,3) :: eDOS,hDOS !(step,normal and integrated)
#     real, allocatable,dimension(:,:) :: edist,hdist, PL
#     real, dimension(eloc,eloc,NoTemps) :: emaster
#     real, dimension(hloc,hloc,NoTemps) :: hmaster
#     integer :: n, m, nn, mm, i, j, k, nnn, mmm, status
#     integer :: rank, lwork
#     integer :: hstnum, estnum, totalt,en, hn
#     real :: preve, prevh, t, tstep, sm, norm
#     real :: eocc(eloc,eloc,NoTemps), hocc(hloc,hloc,NoTemps)
#     real :: hs(hloc), es(eloc),rcond
#     real :: epsi(eloc,eloc),hpsi(hloc,hloc)
#     real :: ewr(eloc),ewi(eloc),evl(eloc,eloc),evr(eloc,eloc)
#     real :: hwr(hloc),hwi(hloc),hvl(hloc,hloc),hvr(hloc,hloc)
#     real :: ewrt(NoTemps,eloc),hwrt(NoTemps,hloc)
#     real :: work(1,2*hloc**2)
#     character(40) :: filename
#     allocate(edist(specsize,NoTemps))
#     allocate(hdist(specsize,NoTemps))
#     allocate(PL(specsize,NoTemps))
#
#     totalt = 1000
#     tstep = 1.0
# !    write(6,*) 'Total time: ',tstep*real(totalt)
#
#     lwork = 2*hloc**2
#     absarr = 0.0d0; selarr = 0.0d0
#     eDOS = 0.0d0; hDOS = 0.0d0
#     edist = 0.0d0; hdist = 0.0d0
#     eocc = 0.0d0; hocc = 0.0d0
#     preve = 0.0d0; prevh = 0.0d0
#     PL = 0.0d0
#     status = 0
#     emaster = 0.0d0; hmaster = 0.0d0
#
#     do n = 1, eloc
#        do nn = 1, eloc
#           do j = 1, NoTemps
#              if (n .ne. nn) then
#                 emaster(n,nn,j) = gldnprb(evarr(n,3),evarr(nn,3),&
#                      &earr(nn,n,6),Tempf(j),1)&
#                      & + gldnprb(evarr(n,3),evarr(nn,3),&
#                      &earr(nn,n,7),Tempf(j),1)
#              else
#                 do nnn = 1, eloc
#                    if (n .ne. nnn) then
#                       emaster(n,nn,j) = emaster(n,nn,j) - &
#                            &(gldnprb(evarr(nnn,3),evarr(n,3),&
#                            &earr(n,nnn,6),Tempf(j),1)&
#                            & + gldnprb(evarr(nnn,3),evarr(n,3),&
#                            &earr(n,nnn,7),Tempf(j),1))
#                    end if
#                 end do
#              end if
#           end do
#        end do
#     end do
#
#     do m = 1, hloc
#        do mm = 1, hloc
#           do j = 1, NoTemps
#              if (m .ne. mm) then
#                 hmaster(m,mm,j) = gldnprb(hvarr(m,3),hvarr(mm,3),&
#                      &harr(mm,m,6),Tempf(j),2)&
#                      & + gldnprb(hvarr(m,3),hvarr(mm,3),&
#                      &harr(mm,m,7),Tempf(j),2)
#              else
#                 do mmm = 1, hloc
#                    if (m .ne. mmm) then
#                       hmaster(m,mm,j) = hmaster(m,mm,j) - &
#                            &(gldnprb(hvarr(mmm,3),hvarr(m,3),&
#                            &harr(m,mmm,6),Tempf(j),2)&
#                            & + gldnprb(hvarr(mmm,3),hvarr(m,3),&
#                            &harr(m,mmm,7),Tempf(j),2))
#                    end if
#                 end do
#              end if
#           end do
#        end do
#     end do
#
#     write(filename,'(A24,A5,I2.2,A4)' ) &
#          &rundir,'stepr',calc,'.dat'
#     open(unit=68,file=filename,status='Unknown')
#
#     write(filename,'(A24,A5,I2.2,A4)' ) &
#          &rundir,'sthpr',calc,'.dat'
#     open(unit=69,file=filename,status='Unknown')
#
# !!$    do j = 1, NoTemps
# !!$       hstnum = 0; estnum = 0
# !!$       work = 0.0d0
# !!$
# !!$       if ( j == 1) then
# !!$          do k = 1, eloc
# !!$             do i = 1, eloc
# !!$                write(68,30) i,k,emaster(i,k,j)
# !!$             end do
# !!$          end do
# !!$
# !!$          do k = 1, hloc
# !!$             do i = 1, hloc
# !!$                write(69,31) i,k,hmaster(i,k,j)
# !!$             end do
# !!$          end do
# !!$       end if
# !!$
# !!$       if (eloc .eq. 1) then
# !!$          eocc(1,1,j) = 1.0d0
# !!$       else
# !!$          call dgeev('V','V',eloc,emaster(:,:,j),eloc,ewr,ewi,evl,eloc,evr,&
# !!$               &eloc,work,lwork,status)
# !!$          if (status .ne. 0) then
# !!$             write(6,*) 'linear solver (electron) status for temp',j,':',status
# !!$             call abort
# !!$          end if
# !!$          ewr(:) = ewr(:)/sum(ewr)
# !!$          ewrt(j,:) = ewr(:)
# !!$          do i = 1, eloc
# !!$             if (abs(ewr(i)) .lt. 1E-8) then
# !!$                estnum = estnum + 1
# !!$                eocc(:,1,j) = eocc(:,1,j) + (evr(:,i))**2
# !!$             end if
# !!$          end do
# !!$          if (estnum .gt. 0) eocc(:,1,j) = eocc(:,1,j)/real(estnum)
# !!$       end if
# !!$       work = 0.0d0
# !!$
# !!$       if (hloc .eq. 1) then
# !!$          hocc(1,1,j) = 1.0d0
# !!$       else
# !!$          call dgeev('V','V',hloc,hmaster(:,:,j),hloc,hwr,hwi,hvl,hloc,hvr,&
# !!$               &hloc,work,lwork,status)
# !!$          if (status .ne. 0) then
# !!$             write(6,*) 'linear solver (hole) status for temp',j,':',status
# !!$             call abort
# !!$          end if
# !!$
# !!$          do i = 1, hloc
# !!$             do k = 1, hloc
# !!$                hvl(:,i) = hvl(:,i)/dot_product(hvl(:,i),hvr(:,i))
# !!$!             write(6,*) i,k,dot_product(hvl(:,i),hvr(:,k))!,&
# !!$!                  &dot_product(hvl(:,i),hvl(:,k))!,&
# !!$!             write(6,*) dot_product(evl(:,i)**2,evr(:,k)**2)!,&
# !!$!                  &dot_product(evr(:,i),evr(:,k))
# !!$             end do
# !!$          end do
# !!$
# !!$          sm = sum(hwr)
# !!$          hwr(:) = hwr(:)/sm
# !!$          hwrt(j,:) = hwr(:)
# !!$
# !!$          if (any(hwi .gt. 0.0d0)) write(6,*) 'COMPLEX NO! - ',cmplx(0.0d0,hwi)
# !!$          do i = 1, hloc
# !!$             if (abs(hwr(i)) .lt. 1E-8) then
# !!$                hstnum = hstnum + 1
# !!$                hocc(:,1,j) = hocc(:,1,j) + (hvr(:,i))**2
# !!$             end if
# !!$          end do
# !!$          if (hstnum .gt. 0) hocc(:,1,j) = hocc(:,1,j)/real(hstnum)
# !!$
# !!$       end if
# !!$
# !!$!       write(6,*) Tempf(j), estnum,hstnum
# !!$
# !!$       write(68,30) Tempf(j),(eocc(i,1,j) ,i=1,eloc)
# !!$       write(69,31) Tempf(j),(hocc(i,1,j),i=1,hloc)
# !!$
# !!$       if ( j == 1 ) then
# !!$          write(filename,'(A24,A5,I2.2,A1,I3.3,A4)' ) &
# !!$               &rundir,'eprob',calc,'T',nint(Tempf(j)),'.dat'
# !!$          open(unit=26,file=filename,status='Unknown')
# !!$          write(filename,'(A24,A5,I2.2,A1,I3.3,A4)' ) &
# !!$               &rundir,'hprob',calc,'T',nint(Tempf(j)),'.dat'
# !!$          open(unit=27,file=filename,status='Unknown')
# !!$
# !!$          do i = 1, totalt
# !!$             t = i*tstep
# !!$             write(26,30) t, (occprob(en,t,eloc,cmplx(ewr,ewi),&
# !!$                  &evl,evr,eloc),en=1,eloc)
# !!$             write(27,31) t, (occprob(hn,t,hloc,cmplx(hwr,hwi),&
# !!$                  &hvl,hvr,hloc),hn=1,hloc)
# !!$          end do
# !!$          close(26)
# !!$          close(27)
# !!$       end if
# !!$    end do
#
# !!$    do i = 1, eloc
# !!$       write(68,32) i,(ewrt(j,i),j=1,NoTemps)
# !!$    end do
# !!$
# !!$    do i = 1, hloc
# !!$       write(69,32) i,(hwrt(j,i),j=1,NoTemps)
# !!$    end do
#
#     close(68)
#     close(69)
#
#     do n = 1, eloc
#        do i = 1, specsize
#           eDOS(i,1) = eDOS(i,1) + heaviside(specE(i)-evarr(n,3))
#           eDOS(i,2) = eDOS(i,2) + onedgauss(specE(i),width,evarr(n,3))
#           eDOS(i,3) = eDOS(i,3) + onedgauss(specE(i),width,evarr(n,3))
#           do j = 1, NoTemps
#              edist(i,j) = edist(i,j) + &
#                   &eocc(n,1,j)*onedgauss(specE(i),width,evarr(n,3))
#           end do
#        end do
#           !Calculate the DOS and distributions
#     end do
#
#     do m = 1, hloc
#        do i = 1, specsize
#           hDOS(i,1) = hDOS(i,1) + heaviside(-(specE(i)-hvarr(m,3)))
#           hDOS(i,2) = hDOS(i,2) + onedgauss(specE(i),width,hvarr(m,3))
#           hDOS(i,3) = hDOS(i,3) + onedgauss(specE(i),width,hvarr(m,3))
#           do j = 1, NoTemps
#              hdist(i,j) = hdist(i,j) + &
#                   &hocc(m,1,j)*onedgauss(specE(i),width,hvarr(m,3))
#           end do
#        end do
#     end do
#     !Repeat for holes
#
#     do i = 2, specsize
#        preve = eDOS(i-1,3)
#        prevh = hDOS(specsize-i+2,3)
#        eDOS(i,3) = eDOS(i,3) + preve
#        hDOS(specsize+1-i,3) = hDOS(specsize+1-i,3) + prevh
#     end do
#
#     eDOS(:,2:3) = eDOS(:,2:3)/real(specsize)
#     hDOS(:,2:3) = hDOS(:,2:3)/real(specsize)
#     !Integrate and normalise DOS
#
#     do m = 1, hloc
#        do n = 1, eloc
#           do i = 1, specsize
#              absarr(i) = absarr(i) + sarr(n,m,5)*&
#                   &onedgauss(specE(i),width,sarr(n,m,1))
#              if (n.eq.m) selarr(i) = selarr(i) + sarr(n,m,5)*&
#                   &onedgauss(specE(i),width,sarr(n,m,1))
#              do j = 1, NoTemps
#                 PL(i,j) = PL(i,j) + sarr(n,m,5)*eocc(n,1,j)*hocc(m,1,j)*&
#                      &onedgauss(specE(i),width,sarr(n,m,1))
#              end do
#           end do
#        end do
#     end do
#     !Caluclate the PL spectra, and overlap spectra
#
# !    if (maxval(PL) .gt. 0.0d0) then
# !       write(6,*) 'Definately something in PL'
# !    end if
#
#     if (maxval(absarr) .gt. 0.0d0) then
#        norm = 1.0d0/maxval(absarr)
#     else
#        norm = 0.0d0
#     end if
#
#     rearr(1,:) = absarr*norm
#     rearr(2,:) = selarr
#     do i = 1, specsize
#        rearr(3:5,i) = eDOS(i,:)
#        rearr(6:8,i) = hDOS(i,:)
#        rearr(9:(8+NoTemps),i) = edist(i,:)
#        rearr((9+NoTemps):(8+(2*NoTemps)),i) = hdist(i,:)
#        rearr((9+(2*NoTemps)):(8+(3*NoTemps)),i) = PL(i,:)
#     end do
#
#     deallocate(edist)
#     deallocate(hdist)
#     deallocate(PL)
#     !Output these spectra
#     return
# 30  format(<(1+eloc)>G24.15E3)
# 31  format(<(1+hloc)>G24.15E3)
# 32  format(<(1+NoTemps)>G24.15E3)
# !31  format(<(1+hloc)>G24.15E3)
#   end subroutine spectra
