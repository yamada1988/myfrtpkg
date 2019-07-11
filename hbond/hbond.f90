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
    integer :: i,j,k, jo, jh1, jh2, ik
    character outf*128
  
    real(8) :: dist
    real(8), parameter :: dlist = 0.60

    type mytype
    real(8), dimension(3):: O
    real(8), dimension(3):: H1
    real(8), dimension(3):: H2
    end type mytype

    type(mytype), allocatable, dimension(:,:) :: tip4p

    type mytype2
    integer, dimension(50) :: pair = -1
    end type mytype2

    type(mytype2), allocatable, dimension(:) :: DistanceList

    

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
    allocate(DistanceList(1:no))

    call xtcf % read
    it = 0 
    pit = 0
    ! read trajectory
    do while ( xtcf % STAT == 0 )
        if (it > totalstep) exit

        ! skip it % skipstep != 0 frame
        if (mod(it, skipstep) /= 0 ) then
          call xtcf % read
          it = it + 1
          cycle
        end if 

        ! output log
        if (mod(it, 10000) == 0) then
          write(*,'(a,f12.6,a,i0)') " Time (ps): ", xtcf % time, "  Step: ", xtcf % STEP
        end if

        ! get descrete time
        if (it == 1) then
          dt = xtcf % time
        end if

        ! get box in it-th step
        box(pit) = xtcf % box(1,1)

        !$omp parallel shared(xtcf, tip4p)
        !$omp private(jo, jh1, jh2)
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
 


    ! calculate Hydrogen bond pair
    do pit = 1, ptotalstep

      ! make list
      if (mod(pit, 10) == 0) then
        do i=1, no
          ik = 1
          do j=1, no
            dist = sqrt(sum((tip4p(pit,i)%O - tip4p(pit,j)%O)**2))
            if (dist <= dlist ) then
              DistanceList(i)%pair(ik) = j
              ik = ik + 1 
            end if
          end do
          print *, pit, i, DistanceList(i)%pair
        end do
      end if

    end do

end program hbond
