!  XDR Fortran Interface XTC Example Program with Wrapper
!  2014 (c) James W. Barnett <jbarnet4@tulane.edu>
!  https://github.com/wesbarnett/
!

program msd

    ! 1. Use the xdr interface
    use xdr, only: xtcfile

    implicit none
    integer :: natm
    integer :: no, np
    double precision, allocatable, dimension(:,:,:)  ::  pos_o_fly
    double precision, allocatable, dimension(:,:)  ::  msdt
    double precision ::  box, dt
    double precision, allocatable, dimension(:,:)  :: r0, r1, rcount
    integer :: it, pit, pit0, totalstep, skipstep, ptotalstep
    integer :: i,j,k
    character outf*128

    ! 2. Declare a variable of type xtcfile
    type(xtcfile) :: xtcf


    ! 3. Initialize it with the names of xtc files you want to read in
    ! and write out
    call xtcf % init("../../../md01_01_skip.xtc")

    natm = xtcf % NATOMS
    no = natm / 3
    no = 100000
    np = 351000
    totalstep = 4000
    skipstep = 1
    ptotalstep = totalstep / skipstep

    !call xtcf % close
    !call xtcf % init("../../Production/md.xtc")
    allocate(pos_o_fly(0:ptotalstep,no,3))
    allocate(msdt(0:ptotalstep,no))
    allocate(r0(3,no))
    allocate(r1(3,no))
    allocate(rcount(3,no))
    msdt = 0.0E0
    rcount = 0.0E0

    ! step = 0
    it = 0
    pit = 0
    call xtcf % read
    do i=1, no
      j = np + 3*(i-1) + 1
      do k=1, 3
        r1(k, i) = xtcf % pos(k, j)
        pos_o_fly(pit, i, k) = r1(k, i) 
      end do
    end do
    !print *, it, pos_o_fly(it,:)
    r0 = r1
    call xtcf % read
    dt = xtcf % time
    it = it + 1
    pit = pit + 1
    !print *, dt

    ! read trajectory
    do while ( xtcf % STAT == 0 )
        if (it > totalstep) exit

        ! skip it % skipstep != 0 frame
        if (mod(it, skipstep) /= 0 ) then
          call xtcf % read
          it = it + 1
          cycle
        end if 

        if (mod(it, 400) == 0) then
          write(*,'(a,f12.7)') " Time (ps): ", xtcf % time 
        end if

        ! get box in it-th step
        box = xtcf % box(1,1)


        !$omp parallel shared(r1, r0, rcount, box, pos_o_fly, pit)
        !$omp do private(k)
        do i=1, no 
          j = np + 3*(i-1) + 1       
          r1(:, i) = xtcf % pos(:, j) 

          ! count +xk boundary and -xk boundary
          do k=1,3
            if ( r1(k, i)-r0(k, i) > 0.8*box ) then
              rcount(k, i) = rcount(k, i) - 1.0
            else if ( r1(k, i)-r0(k, i) < -0.8*box ) then
              rcount(k, i) = rcount(k, i) + 1.0
            end if
            pos_o_fly(pit, i, k) = r1(k, i) + rcount(k, i)*box
          end do
          ! update position
          r0(:, i) = r1(:, i)
        end do

        !$omp end do
        !$omp end parallel
        !print *, it, pit, pos_o_fly(pit,:)

        
        ! call nextstep
        call xtcf % read
        it = it + 1
        pit = pit + 1
    end do
    ! 5. Close the file
    call xtcf % close
 
    !$omp parallel shared(pos_o_fly, msdt), private(i, pit, pit0)
    !$omp do
    do i = 1, no
      write(outf, '("DAT/msds/msd", i7.7, ".xvg")') i
      open(99, file=outf)
        do pit = 0, ptotalstep
          do pit0 = 0, ptotalstep-pit
            !print *, pit, pit0, pos_o_fly(pit+pit0,:),
            !pos_o_fly(pit0,:)
            msdt(pit, i) = msdt(pit, i) + sum(( pos_o_fly(pit+pit0, i, :) - pos_o_fly(pit0,i,:) )**2)
          end do
          msdt(pit, i) = msdt(pit, i) / dble(ptotalstep - pit+1)
          !print *, pit, msdt(pit)
          write(99, '(f12.3, 4X, f12.7)'), pit*skipstep*dt, msdt(pit, i)
        end do
      close(99)
    end do
    !$omp end do
    !$omp end parallel

end program msd

