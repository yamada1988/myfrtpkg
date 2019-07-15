program main
  implicit none
  include 'fftw3.f'
  integer, parameter :: ix = 32
  real(8), dimension(ix) :: x
  real(8), dimension(ix) :: qq1,qq2
  complex(8), dimension(ix/2+1) :: fqq1,fqq2

  real(8) :: plan
  real(8) :: pi = 3.14159265359
  integer :: i
  real(8) :: xmin,xmax,dx

  xmin = 0
  xmax = 2.*pi
  dx = (xmax-xmin)/dble(ix)
  x(1) = 0.5d0*dx
  do i = 2,ix
     x(i) = x(i-1) + dx
  enddo

  do i = 1,ix
     qq1(i) = 1.d0
     qq2(i) = cos(2.d0*x(i))
  enddo

  call dfftw_plan_dft_r2c_1d(plan,ix,qq1,fqq1,FFTW_ESTIMATE)
  call dfftw_execute(plan)
  call dfftw_destroy_plan(plan)

  call dfftw_plan_dft_r2c_1d(plan,ix,qq2,fqq2,FFTW_ESTIMATE)
  call dfftw_execute(plan)
  call dfftw_destroy_plan(plan)

  do i = 1,ix/2+1
  write(*,*) fqq1(i)/dble(ix/2),fqq2(i)/dble(ix/2)
  enddo

  stop
end program main
