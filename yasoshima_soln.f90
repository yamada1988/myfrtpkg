      program main
      implicit none
      save

!input data : set_num; number of sumples
!           : set_line; file lines(weight_refs,engref) 

integer :: set_num=SETNUM
print *, "===== soln ====="
print *, "weight start"
call ave_weight_sln(set_num)

contains
!=======================================================================
subroutine ave_weight_sln(set_num)
implicit none
save

integer :: k,l,i
integer :: access
integer,intent(in) :: set_num
integer :: sample_num
integer :: set_line=5
integer, dimension(:),allocatable :: num
character(len=80) :: fname
real*8,dimension(:),allocatable :: ave_w,sum_w,ave_s
real*8,dimension(:,:),allocatable :: w,s

allocate(num(set_line),ave_w(set_line),sum_w(set_line),ave_s(set_line))
allocate(w(0:set_num,set_line),s(0:set_num,set_line))

 ave_w = 0
 sum_w = 0
 ave_s = 0
 sample_num = 0

 do l=0,set_num-1

  write(fname,"('ERmod_',i4.4,'/soln/weight_soln')") l+1
  if (access( fname, " ") /= 0) cycle
    sample_num = sample_num + 1
    open(10,file=fname)
    do i=1,set_line
     read(10,1200) num(i),w(l,i)
     sum_w(i) = sum_w(i) + w(l,i)
    end do
    close(10)
 end do

 open(34,file="aveERmod_01/soln/log.tt")
 write(34,*) "sample_num=", sample_num
 close(34)

 open(20,file="aveERmod_01/soln/weight_soln")
 do i=1,set_line
  ave_w(i) = sum_w(i)/dble(sample_num)
  ave_s(i) = ave_s(i)/sum_w(i)
  !print *, num(i), ave_w(i)
  write(20,1200) num(i),ave_w(i)
 end do
 close(20)
1200 format(I5,E20.8)

print *, "aveuv start"
call ave_uv_sln(set_num,ave_w,w,sum_w)

return
end subroutine ave_weight_sln
!=======================================================================

!=======================================================================
subroutine ave_uv_sln(set_num,ave_w,w,sum_w)
implicit none
save

integer :: k,l,i,j
integer,intent(in) :: set_num
integer :: sample_num
integer :: access
real*8,dimension(:),intent(in) :: ave_w,sum_w
real*8,dimension(0:,:),intent(in) :: w
integer :: set_line=5
character(len=80) :: fname,fname2
character(len=25),dimension(:),allocatable :: line
real*8,dimension(:),allocatable :: ave_e
real*8,dimension(:,:),allocatable :: e_num

allocate(ave_e(set_line))
allocate(e_num(0:set_num,set_line))

  ave_e = 0.d0
  e_num = 0.0

  do l=0,set_num-1
   write(fname,"('ERmod_',i4.4,'/soln/aveuv.tt')") l+1
   if (access( fname, " ") /= 0) cycle
     sample_num = sample_num + 1
     open(10,file=fname)
     do i=1,set_line
       read(10,*) j,e_num(l,i) 
       !print *, l+1, i, e_num(l,i)
       ave_e(i) = ave_e(i) + e_num(l,i)*w(l,i)
     end do
     close(10)
  end do
  
  write(fname2,"('aveERmod_01/soln/aveuv.tt')") 
  open(20,file=fname2)
  do i=1,set_line
   ave_e(i) = ave_e(i)/sum_w(i)
   print *, i, ave_e(i)
   write(20,*) i, ave_e(i)
  end do
  close(20)
!700 format(A12,f7.5)

print *, "eng start"
call ave_eng_sln(set_num,ave_w,w,sum_w)

return
end subroutine ave_uv_sln
!=======================================================================


!=======================================================================
subroutine ave_eng_sln(set_num,ave_w,w,sum_w)
implicit none
save

integer :: k,l,i
integer,intent(in) :: set_num
integer :: sample_num
integer :: access
real*8,dimension(:),intent(in) :: ave_w,sum_w
real*8,dimension(0:,:),intent(in) :: w
integer :: set_line=1496
character(len=80) :: fname,fname2
character(len=25),dimension(:),allocatable :: line
real*8,dimension(:),allocatable :: ave_e
real*8,dimension(:,:),allocatable :: e_num

allocate(line(set_line),ave_e(set_line))
allocate(e_num(0:set_num,set_line))

 do k=1,5
  ave_e = 0.d0
  e_num = 0.0

  do l=0,set_num-1
   write(fname,"('ERmod_',i4.4,'/soln/engsln.',i2.2)") l+1,k
   if (access( fname, " ") /= 0) cycle
     sample_num = sample_num + 1
     open(10,file=fname)
     do i=1,set_line
      read(10,660) line(i), e_num(l,i)
      ave_e(i) = ave_e(i) + e_num(l,i)*w(l,k)
     end do
     close(10)
  end do

  write(fname2,"('aveERmod_01/soln/engsln.'i2.2)") k
  open(20,file=fname2)
  do i=1,set_line
   ave_e(i) = ave_e(i)/sum_w(k)
   write(20,660) line(i), ave_e(i)
  end do
  close(20)
 end do
600 format(A20,E25.15e2)
660 format(A20,g25.15)

return
end subroutine ave_eng_sln
!=======================================================================



!=======================================================================
end program main
