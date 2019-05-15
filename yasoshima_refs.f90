      program main
      implicit none
      save

!input data : set_num; number of sumples
!           : set_line; file lines(weight_refs,engref) 

integer :: set_numref = SETNUM
print *, "===== refs ====="
print *, "weight start"
call ave_weight(set_numref)

contains
!=======================================================================
subroutine ave_weight(set_numref)
implicit none
save

integer :: k,l,i
integer,intent(in) :: set_numref
integer :: sample_num
integer :: access
integer :: set_line=5
integer, dimension(:),allocatable :: num
character(len=80) :: fname
real*8,dimension(:),allocatable :: ave_w,sum_w,ave_s
real*8,dimension(:,:),allocatable :: w,s

allocate(num(set_line),ave_w(set_line),sum_w(set_line),ave_s(set_line))
allocate(w(0:set_numref,set_line),s(0:set_numref,set_line))

 ave_w = 0
 sum_w = 0
 ave_s = 0
 sample_num = 0

 do l=0,set_numref-1

  write(fname,"('ERmod_',i4.4,'/refs/weight_soln')") l+1
  if (access( fname, " ") /= 0) cycle
    sample_num = sample_num + 1
    open(10,file=fname)
    do i=1,set_line
     read(10,1000) num(i),w(l,i),s(l,i)
     sum_w(i) = sum_w(i) + w(l,i)
     ave_s(i) = ave_s(i) + s(l,i) * w(l,i)
    end do
    close(10)
 end do

 open(34,file="aveERmod_01/refs/log.tt")
 write(34,*) "sample_num=", sample_num
 close(34)

 open(20,file="aveERmod_01/refs/weight_refs")
 do i=1,set_line
  ave_w(i) = sum_w(i)/dble(sample_num)
  ave_s(i) = ave_s(i)/sum_w(i)
  write(20,1000) num(i),ave_w(i),ave_s(i)
 end do
 close(20)
1000 format(I5,E20.8,E20.7)

print *, "eng start"
call ave_eng(set_numref,ave_w,w,sum_w)
print *, "all jobs end"

return
end subroutine ave_weight
!=======================================================================


!=======================================================================
subroutine ave_eng(set_numref,ave_w,w,sum_w)
implicit none
save

integer :: k,l,i
integer,intent(in) :: set_numref
integer :: sample_num
integer :: access
real*8,dimension(:),intent(in) :: ave_w,sum_w
real*8,dimension(0:,:),intent(in) :: w
integer :: set_line=1496
character(len=80) :: fname,fname2,fname3
character(len=25),dimension(:),allocatable :: line
real*8,dimension(:),allocatable :: ave_e
real*8,dimension(:,:),allocatable :: e_num

allocate(line(set_line),ave_e(set_line))
allocate(e_num(0:set_numref,set_line))

 sample_num = 0
 do k=1,5
  ave_e = 0.d0
  e_num = 0.0

  do l=0,set_numref-1
   write(fname,"('ERmod_',i4.4,'/refs/engsln.',i2.2)") l+1,k
   write(fname3,"('ERmod_',i4.4,'/refs/weight_soln')") l+1
   if (access( fname3, " ") /= 0) cycle
    sample_num = sample_num + 1
    open(10,file=fname)
    do i=1,set_line
     read(10,660) line(i), e_num(l,i)
     ave_e(i) = ave_e(i) + e_num(l,i)*w(l,k)
    end do
   close(10)
  end do

  write(fname2,"('aveERmod_01/refs/engref.'i2.2)") k
  open(20,file=fname2)
  do i=1,set_line
   ave_e(i) = ave_e(i)/sum_w(k)
   write(20,660) line(i), ave_e(i)
  end do
  close(20)
 end do
600 format(A20,E25.15e2)
660 format(A20,g25.15)

print *, "corref start"
call ave_cor(set_numref,ave_w,w,sum_w,set_line)

return
end subroutine ave_eng
!=======================================================================

!=======================================================================
subroutine ave_cor(set_numref,ave_w,w,sum_w,set_line)
implicit none
save

integer :: k,l,i,j
integer,intent(in) :: set_numref, set_line
integer :: sample_num
integer :: access
real*8,dimension(:),intent(in) :: ave_w,sum_w
real*8,dimension(0:,:),intent(in) :: w
real*8,dimension(:,:), allocatable :: cormat_temp,rdcor
character(len=80) :: fname,fname2,fname3

allocate(cormat_temp(set_line,set_line),rdcor(set_line,set_line))

 do k=1,5
  cormat_temp = 0
  rdcor = 0

  do l=0,set_numref-1
   write(fname,"('ERmod_',i4.4,'/refs/corsln.',i2.2)") l+1,k
   write(fname3,"('ERmod_',i4.4,'/refs/weight_soln')") l+1
   if (access( fname3, " ") /= 0) cycle
    sample_num = sample_num + 1
    open(10,file=fname,status='old',form="UNFORMATTED")
    read(10) cormat_temp

    do i=1,set_line
     do j=1,set_line
      rdcor(i,j) = rdcor(i,j) + cormat_temp(i,j)*w(l,k)
     end do
    end do
  end do

  write(fname2,"('aveERmod_01/refs/corref.'i2.2)") k
   if (access( fname, " ") /= 0) cycle
    sample_num = sample_num + 1
    open(20,file=fname2,form="UNFORMATTED")
    do i=1,set_line
     do j=1,set_line
      rdcor(i,j) = rdcor(i,j)/sum_w(k)
     end do
    end do
  write(20) rdcor
  close(20)
 end do

return
end subroutine ave_cor

!=======================================================================
end program main
