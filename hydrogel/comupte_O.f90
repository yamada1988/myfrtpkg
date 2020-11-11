!  XDR Fortran Interface XTC Example Program with Wrapper
!  2014 (c) James W. Barnett <jbarnet4@tulane.edu>
!  https://github.com/wesbarnett/
!

program example
    !$ use omp_lib
    use xdr, only: xtcfile
    implicit none
    type(xtcfile) :: xtc
    real(kind=8), parameter :: pi = 3.141592653589793238462643383
    real(kind=8), allocatable :: pos_PVA(:, :), pos_O(:, :), dist(:)
    real(kind=8) :: box_dim, xr(3), dxr
    integer :: i, j, it, jt, nt
    character(50) :: fname, dummy
    integer :: npos, npos_SOL, npos_O
    integer, parameter :: npos_PVA = 351000 
    integer, parameter :: Nmax = Nt
    real(kind=16), parameter :: dt = 1.0E0

    !call xtc % init("/data/yamada/Data/2020/ohtaka/PVA100/80wt%_N500_43nm/md01_01_skip.xtc")
    call xtc % init("../../../md01_01_skip.xtc")
    npos = xtc % NATOMS   !number of water molecules (NATOMS is obtained after calling init)
    npos_SOL = npos - npos_PVA
    npos_O = npos_SOL / 3
    allocate(pos_PVA(3, npos_PVA))
    allocate(pos_O(3, npos_O))
    allocate(dist(npos_O))

    do it=1, Nmax+100
      !write(*,*) it
      call xtc % read
      if ( it < Nmax ) cycle
      !if ( mod(it, 10) /= 0) cycle
      box_dim = xtc % box(1,1)

      do i = 1, npos_PVA
        pos_PVA(:, i) = xtc % pos(:, i)  
      end do

      do i = 1, npos_O
        pos_O(:, i) = xtc % pos(:, npos_PVA+3*i-2)
      end do

    !$omp parallel private(j, xr, dxr)
    !$omp do         
      do i = 1, npos_O
        !write(*,*) i
        dist(i) = 10.0E0
        do j = 1, npos_PVA
          xr = pos_O(:, i) - pos_PVA(:, j)
          xr = xr - box_dim * nint(xr/box_dim)  !wrap distance (periodic boundary condition)
          dxr = sqrt(sum(xr**2.0E0))
          if ( dxr <= dist(i) ) then
            dist(i) = dxr
          end if
        end do
      end do
    !$omp end do
    !$omp end parallel 

     write(*,*) dble(it)*dt, dist(1) 
     write(fname, "('DAT/posO'i8.8'.dat')") int(dble(it)*dt*1000)
     open (17,file=fname)
     if ( mod(it, 50) == 0) then
       write(*,*) int(it*dt)
     end if
     write(17, '(A, 2X, f8.5)') , '#', box_dim
     do i = 1, npos_O
         write (17, '(f8.5,2X,f8.5,2X,f8.5,2X,f7.5)') pos_O(1, i), pos_O(2, i), pos_O(3, i), dist(i)
     end do
     close(17)  

   end do
   

end program example
