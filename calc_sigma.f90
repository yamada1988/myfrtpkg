program calc_sigma
  implicit none
  integer i, i0, j, k
  double precision, allocatable,dimension(:,:) :: P
  double precision, allocatable,dimension(:) :: t
  integer n
  double precision, allocatable,dimension(:,:) :: G
  n = num
  allocate( P(3,n) )
  allocate( t(n) )
  allocate( G(3,n) )
  G = 0.0

  open(99, file='sigma.xvg' )
  do j = 1,26
  read (99,'()')
  end do

  do i = 1, n
    read (99,*) t(i), P(1,i), P(2,i), P(3,i)
  end do
  close(99)

  do i = 1, n
    do k = 1,3
      do i0 = 1, n-i
        G(k,i) = G(k,i) + P(k,i+i0)*P(k,i0)
      end do
      G(k,i) = G(k,i)/(n-i+1)
    end do
  end do


  open(98, file='Gxy.xvg')
  do i = 1, n
    write(98, "(f9.4, f20.11)") t(i), G(1,i)
  end do
  close(98)

  open(97, file='Gxz.xvg')
  do i = 1, n
    write(97, "(f9.4, f20.11)") t(i), G(2,i)
  end do
  close(97)

  open(96, file='Gyz.xvg')
  do i = 1, n
    write(96, "(f9.4, f20.11)") t(i), G(3,i)
  end do
  close(96)  

 
end
