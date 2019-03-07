implicit none

integer :: k,l,i,j
integer :: lnum=32, set_line=1496
real*8,dimension(:,:), allocatable :: cormat_temp,rdcor
character(len=80) :: fname, ofname 

allocate(cormat_temp(set_line,set_line),rdcor(set_line,set_line))

do k=1,5
  do l = 1,lnum
  cormat_temp = 0
  rdcor = 0

   write(fname,"('ERmod_',i4.4,'/refs/corsln.',i2.2)") l, k
   open(10,file=fname,status='old',form="UNFORMATTED")
   read(10) cormat_temp
   close(10)

   do i=1,set_line
    do j=1,set_line
     rdcor(i,j) = rdcor(i,j) + cormat_temp(i,j)
    end do
   end do
  end do

  write(ofname,"('aveERmod_01/refs/corref.'i2.2)") k
  open(20,file=ofname,form="UNFORMATTED")
  do i=1,set_line
   do j=1,set_line
    rdcor(i,j) = rdcor(i,j)/lnum
   end do
  end do
  write(20) rdcor
  close(20)
end do




stop
end
