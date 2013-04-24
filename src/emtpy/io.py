__author__ = 'duncan'



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