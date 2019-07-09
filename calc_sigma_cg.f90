 program calc_sigma
  implicit none
  integer i, i0, j, k
  double precision, allocatable,dimension(:,:) :: P
  double precision, allocatable,dimension(:) :: t
  double precision, allocatable,dimension(:) :: t_cg, Geta
  integer n,m
  double precision, allocatable,dimension(:,:) :: G
  double precision, allocatable,dimension(:,:) :: P_cg
  integer blockn
  double precision :: V, kBT, alpha
  V = 4.725E-24 
  kBT = 4.114E-21
  alpha = V/kBT
  alpha = 1.0

  blockn = 50


  n = num
  m = int(n/blockn)
  allocate( P(3,n) )
  allocate( P_cg(3, m) )
  allocate( t(n) )
  allocate( G(3,0:n) )
  allocate( t_cg(m) )
  allocate( Geta(0:m) )

  G = 0.0
  P_cg = 0.0

  open(99, file='sigma.xvg' )
  do j = 1,26
  read (99,'()')
  end do

  do i = 1, n
    read (99,*) t(i), P(1,i), P(2,i), P(3,i)
  end do
  close(99)

  do i = 1, m
    do j=(i-1)*blockn+1, i*blockn
      t_cg(i) = t_cg(i) + t(j)
    end do
    t_cg(i) = t_cg(i)/blockn
  end do

  do i = 1, m
    do k = 1,3
      do j = (i-1)*blockn+1, i*blockn
        P_cg(k,i) = P_cg(k,i) + P(k,j)
      end do
      P_cg(k,i) = P_cg(k,i)/blockn
    end do
  end do

!$omp parallel
!$omp do
  do i = 0, m
    do k = 1,3
      do i0 = 1, m-i
        G(k,i) = G(k,i) + P_cg(k,i+i0)*P_cg(k,i0)
      end do
      G(k,i) = G(k,i)/(n-i+1)
    end do
  end do
!$omp end do
!$omp end parallel

  G = G * alpha

  open(98, file='Gxy_cg.xvg')
  do i = 1, m
    write(98, "(f9.4, f20.11)") t_cg(i), G(1,i)
  end do
  close(98)

  open(97, file='Gxz_cg.xvg')
  do i = 1, m
    write(97, "(f9.4, f20.11)") t_cg(i), G(2,i)
  end do
  close(97)

  open(96, file='Gyz_cg.xvg')
  do i = 1, m
    write(96, "(f9.4, f20.11)") t_cg(i), G(3,i)
  end do
  close(96)  

do i=0, m-1
  do k=1,3
    Geta(i) = Geta(i) + G(k,i)
  end do
  Geta(i) = Geta(i)/dble(3.0)
end do

  open(95, file='Geta_cg.xvg')
  do i = 0, m-1
    write(95, "(f9.4, f20.11)") t(i), Geta(i)
  end do 
  close(95)

end
