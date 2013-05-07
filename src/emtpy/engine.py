__author__ = 'pard'


class Engine(object):

    def __init__(self, mat_dist, potential_energy):
        # This check doesn't currently work with MockMaterialDistribution so I might need a better one
#        if not mat_dist.compatible_grid(potential_energy):
#            raise ValueError
        self.mat_dist = mat_dist
        self.pot_energy = potential_energy
        self.hamiltonian = self._create_hamiltonian()

    def _create_hamiltonian(self):
        pass

    def solve(self):
        pass


class ARPACK3DSolver(Engine):

    def _create_hamiltonian(self):
        import numpy as np
        from scipy.sparse.coo import coo_matrix

        i_index = np.zeros(self.pot_energy.no_elements*7)
        j_index = np.zeros(self.pot_energy.no_elements*7)
        data = np.zeros(self.pot_energy.no_elements*7)

        row = 0
        # Loop over the potential energy grid. Every element in this array corresponds to a row in the Hamiltonaian
        for idx, val in np.ndenumerate(self.pot_energy.values):
            non_zero_index = row*7

            # Calculate where the offset matrix elements are in the hamiltonian
            idx_plus_x = ((idx[0]+1) % self.pot_energy.size[0], idx[1], idx[2])
            idx_plus_y = (idx[0], (idx[1]+1) % self.pot_energy.size[1], idx[2])
            idx_plus_z = (idx[0], idx[1], (idx[2]+1) % self.pot_energy.size[2])
            n_x = self.pot_energy.getn(*idx_plus_x)
            n_y = self.pot_energy.getn(*idx_plus_y)
            n_z = self.pot_energy.getn(*idx_plus_z)

            # Calculate the diagonal
            data[non_zero_index] = self.dn4(idx)
            i_index[non_zero_index] = row
            j_index[non_zero_index] = row

            # Calculate the x offset diagonal
            data[non_zero_index+1] = self.an3(idx_plus_x)
            i_index[non_zero_index+1] = row
            j_index[non_zero_index+1] = n_x

            data[non_zero_index+4] = self.an3(idx_plus_x)
            i_index[non_zero_index+4] = n_x
            j_index[non_zero_index+4] = row

            # Calculate the y offset diagonal
            data[non_zero_index+2] = self.bn3(idx_plus_y)
            i_index[non_zero_index+2] = row
            j_index[non_zero_index+2] = n_y

            data[non_zero_index+5] = self.bn3(idx_plus_y)
            i_index[non_zero_index+5] = n_y
            j_index[non_zero_index+5] = row

            # Calculate the z offset diagonal
            data[non_zero_index+3] = self.cn3(idx_plus_z)
            i_index[non_zero_index+3] = row
            j_index[non_zero_index+3] = n_z

            data[non_zero_index+6] = self.cn3(idx_plus_z)
            i_index[non_zero_index+6] = n_z
            j_index[non_zero_index+6] = row

            row += 1

        # Create the sparse matrix
        coo = coo_matrix((data, (i_index, j_index)), shape=(self.pot_energy.no_elements, self.pot_energy.no_elements))

        # Convert to compressed row format
        return coo.tocsr()

    def an3(self, idx):
        from constants import units
        previous_x_index = (idx[0]-1, idx[1], idx[2])
        grid_spacing = self.mat_dist.increments[0]
        return (-units / (4.0*(grid_spacing**2)))*(self.mat_dist.inv_mass_xy(idx) + self.mat_dist.inv_mass_xy(previous_x_index))

    def bn3(self, idx):
        from constants import units
        previous_y_index = (idx[0], idx[1]-1, idx[2])
        grid_spacing = self.mat_dist.increments[1]
        return (-units / (4.0*(grid_spacing**2)))*(self.mat_dist.inv_mass_xy(idx) + self.mat_dist.inv_mass_xy(previous_y_index))

    def cn3(self, idx):
        from constants import units
        previous_z_index = (idx[0], idx[1], idx[2]-1)
        grid_spacing = self.mat_dist.increments[2]
        return (-units / (4.0*(grid_spacing**2)))*(self.mat_dist.inv_mass_z(idx) + self.mat_dist.inv_mass_z(previous_z_index))

    def dn4(self, idx):
        from constants import units
        previous_x_index = (idx[0]-1, idx[1], idx[2])
        next_x_index = (idx[0]+1-self.mat_dist.size[0], idx[1], idx[2])
        alphazero = (self.mat_dist.inv_mass_xy(previous_x_index) + 2.0*self.mat_dist.inv_mass_xy(idx) +
                     self.mat_dist.inv_mass_xy(next_x_index)) / self.mat_dist.increments[0]**2

        previous_y_index = (idx[0], idx[1]-1, idx[2])
        next_y_index = (idx[0], idx[1]+1-self.mat_dist.size[1], idx[2])
        betazero = (self.mat_dist.inv_mass_xy(previous_y_index) + 2.0*self.mat_dist.inv_mass_xy(idx) +
                    self.mat_dist.inv_mass_xy(next_y_index)) / self.mat_dist.increments[1]**2

        previous_z_index = (idx[0], idx[1], idx[2]-1)
        next_z_index = (idx[0], idx[1], idx[2]+1-self.mat_dist.size[2])
        gammazero = (self.mat_dist.inv_mass_z(previous_z_index) + 2.0*self.mat_dist.inv_mass_z(idx) +
                     self.mat_dist.inv_mass_z(next_z_index)) / self.mat_dist.increments[2]**2

        return units*0.25*(alphazero + betazero + gammazero) + self.pot_energy.values[idx]

    def solve(self):
        from scipy.sparse.linalg import eigsh
        return eigsh(self.hamiltonian)


class ARPACKSolver(Engine):

    def _create_hamiltonian(self):
        import numpy as np
        from scipy.sparse.coo import coo_matrix

        i_index = np.zeros(self.pot_energy.no_elements*4)
        j_index = np.zeros(self.pot_energy.no_elements*4)
        data = np.zeros(self.pot_energy.no_elements*4)

        off_diag_elements = [self.an3, self.bn3, self.cn3]

        row = 0
        # Loop over the potential energy grid. Every element in this array corresponds to a row in the Hamiltonaian
        for idx, val in np.ndenumerate(self.pot_energy.values):
            non_zero_index = row*4

            # Calculate where the offset matrix elements are in the hamiltonian
            idx_plus = []
            idx_minus = []
            central_diffs = []
            n = []
            for dim in range(len(idx)):
                idx_plus.append(list(idx))
                idx_plus[dim][dim] = (idx[dim] + 1) % self.pot_energy.size[dim]
                idx_minus.append(list(idx))
                idx_minus[dim][dim] = (idx[dim] - 1) % self.pot_energy.size[dim]

                n.append(self.pot_energy.getn(idx_plus[dim]))
            # idx_plus_x = ((idx[0]+1) % self.pot_energy.size[0], idx[1], idx[2])
            # idx_plus_y = (idx[0], (idx[1]+1) % self.pot_energy.size[1], idx[2])
            # idx_plus_z = (idx[0], idx[1], (idx[2]+1) % self.pot_energy.size[2])
            # n_x = self.pot_energy.getn(idx_plus_x)
            # n_y = self.pot_energy.getn(idx_plus_y)
            # n_z = self.pot_energy.getn(idx_plus_z)


            # Calculate the offset diagonals
            for dim, idx_p in enumerate(idx_plus):
                data[non_zero_index+dim] = off_diag_elements[dim](idx_p, idx_minus)
                central_diffs.append(self.central_dif(idx_minus, idx, idx_plus, dim))
                i_index[non_zero_index+dim] = row
                j_index[non_zero_index+dim] = n[dim]

            # Calculate the diagonal
            data[non_zero_index] = self.dn4(idx, central_diffs)
            i_index[non_zero_index] = row
            j_index[non_zero_index] = row


            # # Calculate the y offset diagonal
            # data[non_zero_index+2] = self.bn3(idx_plus_y)
            # i_index[non_zero_index+2] = row
            # j_index[non_zero_index+2] = n_y
            #
            # # Calculate the z offset diagonal
            # data[non_zero_index+3] = self.cn3(idx_plus_z)
            # i_index[non_zero_index+3] = row
            # j_index[non_zero_index+3] = n_z

            row += 1

        # Create the sparse matrix
        coo = coo_matrix((data, (i_index, j_index)), shape=(self.pot_energy.no_elements, self.pot_energy.no_elements))

        # Convert to compressed row format
        return coo.tocsr()

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

    def an3(self, idx, previous_x_index):
        from constants import units
        grid_spacing = self.mat_dist.increments[0]
        return (-units / (4.0*(grid_spacing**2)))*(self.mat_dist.inv_mass(idx, 0) + self.mat_dist.inv_mass(previous_x_index, 0))

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

    def bn3(self, idx, previous_y_index):
        from constants import units
        grid_spacing = self.mat_dist.increments[1]
        return (-units / (4.0*(grid_spacing**2)))*(self.mat_dist.inv_mass(idx, 0) + self.mat_dist.inv_mass(previous_y_index, 0))

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

    def cn3(self, idx, previous_z_index):
        from constants import units
        grid_spacing = self.mat_dist.increments[2]
        return (-units / (4.0*(grid_spacing**2)))*(self.mat_dist.inv_mass(idx, 1) + self.mat_dist.inv_mass(previous_z_index, 1))


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

    def central_dif(self, prev_idx, idx, next_idx, dim):
        return (self.mat_dist.inv_mass(prev_idx, dim) + 2.0*self.mat_dist.inv_mass(idx, dim) +
                self.mat_dist.inv_mass(next_idx, dim)) / self.mat_dist.increments[dim]**2

    def dn4(self, idx, central_diffs):
        from constants import units
        return units*0.25*(sum(central_diffs)) + self.pot_energy.values[idx]

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

    def solve(self):
        from scipy.sparse.linalg import eigsh
        return eigsh(self.hamiltonian, k=2, which='SM', ncv=100)


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