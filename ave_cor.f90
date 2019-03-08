implicit none

integer :: k,l,i,j
integer :: lnum=32, set_line=1496
integer :: wset_line=5
real*8,dimension(:,:), allocatable :: cormat_temp,rdcor
character(len=80) :: fname, ofname, wname
integer, dimension(:),allocatable :: num
real*8,dimension(:),allocatable :: ave_w,sum_w,ave_s
real*8,dimension(:,:),allocatable :: w,s

allocate(cormat_temp(set_line,set_line),rdcor(set_line,set_line))
allocate(num(wset_line),ave_w(wset_line),sum_w(wset_line),ave_s(wset_line))
allocate(w(0:lnum,wset_line),s(0:lnum,wset_line))

 ave_w = 0
 sum_w = 0
 ave_s = 0

 do l=0,lnum-1

  write(wname,"('ERmod_',i4.4,'/refs/weight_soln')") l+1
  open(30,file=wname)
  do i=1,wset_line
   write(*,*) i
   read(30,1000) num(i),w(l,i),s(l,i)
   sum_w(i) = sum_w(i) + w(l,i)
   ave_s(i) = ave_s(i) + s(l,i) * w(l,i)
  end do
  close(30)
 end do

 open(40,file="aveERmod_01/refs/weight_refs")
 do i=1,5
  ave_w(i) = sum_w(i)/dble(lnum)
  ave_s(i) = ave_s(i)/sum_w(i)
  write(40,1000) num(i),ave_w(i),ave_s(i)
 end do
 close(40)

1000 format(I5,E20.8,E20.7)

do k=1,5
  cormat_temp = 0
  rdcor = 0
  do l = 0,lnum-1

   write(fname,"('ERmod_',i4.4,'/refs/corsln.',i2.2)") l+1, k
   open(10,file=fname,status='old',form="UNFORMATTED")
   read(10) cormat_temp
   close(10)

   do i=1,set_line
    do j=1,set_line
     rdcor(i,j) = rdcor(i,j) + cormat_temp(i,j)*w(l,k)
    end do
   end do
  end do

  write(ofname,"('aveERmod_01/refs/corref.'i2.2)") k
  open(20,file=ofname,form="UNFORMATTED")
  do i=1,set_line
   do j=1,set_line
    rdcor(i,j) = rdcor(i,j)/sum_w(k)
   end do
  end do
  write(20) rdcor
  close(20)
end do




stop
end
