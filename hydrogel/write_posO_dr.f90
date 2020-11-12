!  XDR Fortran Interface XTC Example Program with Wrapper
!  2014 (c) James W. Barnett <jbarnet4@tulane.edu>
!  https://github.com/wesbarnett/
!

program example
    implicit none
    character(100), allocatable :: line(:)
    real(kind=8), allocatable :: pos_O(:, :), dr(:) 
    integer :: i, j, it, jt, nt
    character(50) :: fname, ofname
    integer, parameter :: Nist = nmst, Niend = 20000
    integer, parameter :: Nst = 1, Nend = 4000
    real(kind=16), parameter :: dt = 1.0E0
    real(kind=8) :: t

    allocate(line(Nend-Nst+10))
    allocate(pos_O(3, Nend-Nst+10))
    allocate(dr(Nend-Nst+10))

    do i = Nist, Niend
      write(ofname, "('DAT/atoms/log'i7.7'.tt')") int(i)
      open(18, file=ofname, status='replace')
      jt = 1 
 
      do it = Nst, Nend
        !write(*, *) i, it
        t = it * dt
        write(fname, "('DAT/posO'i8.8'.dat')") int(dble(it)*dt*1000)
        open (17, file=fname, status='old')
        do j = 1, i
          !write(*, *) j
          read(17, '()') ! skip Nmax columns
        end do
        read(17, *) pos_O(1, jt), pos_O(2, jt), pos_O(3, jt), dr(jt)
        close(17)  
        !write(*,*) pos_O(1, jt), pos_O(2, jt), pos_O(3, jt), dr(jt)

        !write (*, *) line(jt)
        write(18, '(f7.2, 3X, f7.4, 3X, f7.4, 3X, f7.4, 3X, f6.4)') t, pos_O(1, jt), pos_O(2, jt), pos_O(3, jt), dr(jt)
        jt = jt + 1
      end do

      close(18)
   end do

end program example
