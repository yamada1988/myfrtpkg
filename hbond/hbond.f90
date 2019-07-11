!  XDR Fortran Interface XTC Example Program with Wrapper
!  2014 (c) James W. Barnett <jbarnet4@tulane.edu>
!  https://github.com/wesbarnett/
!

program hbond

    ! 1. Use the xdr interface
    use xdr, only: xtcfile

    implicit none
    integer :: natm
    integer :: no
    double precision ::  dt
    double precision, allocatable, dimension(:) :: box
    integer :: it, pit, pit0, totalstep, skipstep, ptotalstep
    integer :: i,j,k, jo, jh1, jh2
    character outf*128

    type mytype
    real(8), dimension(3):: O
    real(8), dimension(3):: H1
    real(8), dimension(3):: H2
    end type mytype

    type(mytype), allocatable, dimension(:,:) :: tip4p


    ! 2. Declare a variable of type xtcfile
    type(xtcfile) :: xtcf


    ! 3. Initialize it with the names of xtc files you want to read in and write out
    call xtcf % init("../../Production/md.xtc")

    natm = xtcf % NATOMS
    no = natm / 4
    no = 1000
    totalstep = 10000
    skipstep = 100
    ptotalstep = totalstep / skipstep
    allocate(tip4p(0:ptotalstep, 1:no))
    allocate(box(0:ptotalstep))
 
    ! step = 0
    it = 0
    pit = 0
    call xtcf % read
    do i=1, no
      jo = 4*(i-1) + 1
      jh1 = jo + 1
      jh2 = jo + 2 
      tip4p(pit,i) % O = xtcf % pos(:, jo)
      tip4p(pit,i) % H1 = xtcf % pos(:, jh1)
      tip4p(pit,i) % H2 = xtcf % pos(:, jh2)  
    end do
    box(pit) = xtcf % box(1,1)
    call xtcf % read
    dt = xtcf % time
    it = it + 1
    pit = pit + 1

    ! read trajectory
    do while ( xtcf % STAT == 0 )
        if (it > totalstep) exit

        ! skip it % skipstep != 0 frame
        if (mod(it, skipstep) /= 0 ) then
          call xtcf % read
          it = it + 1
          cycle
        end if 

        if (mod(it, 10000) == 0) then
          write(*,'(a,f12.6,a,i0)') " Time (ps): ", xtcf % time, "  Step: ", xtcf % STEP
        end if

        ! get box in it-th step
        box(pit) = xtcf % box(1,1)

        !$omp parallel shared(xtcf, tip4p)
        !$omp private(jo, jh1, jh2, box,)
        !$omp do
        ! get position (j-th O) and box at it-th step
        do i=1, no
          jo = 4*(i-1) + 1
          jh1 = jo + 1
          jh2 = jo + 2
          tip4p(pit,i) % O = xtcf % pos(:, jo)
          tip4p(pit,i) % H1 = xtcf % pos(:, jh1)
          tip4p(pit,i) % H2 = xtcf % pos(:, jh2)
          print *, it, i, tip4p(pit,i)%O
        end do
        !$omp end do
        !$omp end parallel
        
        ! call nextstep
        call xtcf % read
        it = it + 1
        pit = pit + 1
    end do
    ! 5. Close the file
    call xtcf % close
 
    print *, tip4p

end program hbond
