program main
implicit none
integer :: i_hbond
real(8),allocatable,dimension(:)  :: sm_val
integer, allocatable, dimension(:)  :: sm_indx, sm_jndx
integer :: k,j

open(10, file='hbond.sm', form='unformatted')
do k=1,20
  read(10) i_hbond
  print *, i_hbond
  allocate(sm_val(i_hbond), sm_indx(i_hbond), sm_jndx(i_hbond))
  read(10) sm_val, sm_indx, sm_jndx
  !print *, sm_val, sm_indx, sm_jndx
  deallocate(sm_val, sm_indx, sm_jndx)
end do

close(10)
end program main
