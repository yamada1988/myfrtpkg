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
    integer :: it, pit, pit0, ptotalstep
    integer :: i,j,k, jo, jh1, jh2, ik, jk
    character outf*128

 
    integer, parameter :: nmlfu=20
    real(8), parameter :: PI=3.1415927
    real(8) :: dist, cosij, dcosij
    real(8), dimension(3) :: dr, r_hjoi, r_hioi, u_hjoi, u_hioi
    character(len=200) :: inpfile, outfile
    integer :: totalstep, skipstep, logstep, liststep, lmax
    real(8) :: dlist, dhbond, dtheta
    namelist/indata/inpfile,outfile
    namelist/inparam/totalstep,skipstep,logstep,liststep,dlist,lmax,dhbond,dtheta


    type mytype
    real(8), dimension(3):: O
    real(8), dimension(3):: H1
    real(8), dimension(3):: H2
    end type mytype

    type(mytype), allocatable, dimension(:,:) :: tip4p

    type mytype2
    integer, dimension(200) :: pair = -1
    end type mytype2

    type(mytype2), allocatable, dimension(:) :: DistanceList

    ! 2. Declare a variable of type xtcfile
    type(xtcfile) :: xtcf


    ! 2.1 Use namelist
    open(unit=nmlfu, file='namelist')
    read(unit=nmlfu, nml=indata)
    read(unit=nmlfu, nml=inparam)
    close(unit=nmlfu)


    ! 3. Initialize it with the names of xtc files you want to read in and write out
    call xtcf % init(inpfile)

    natm = xtcf % NATOMS
    no = natm / 4
    ptotalstep = totalstep / skipstep
    dcosij = cos(dtheta*PI/180.0)
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
        if (mod(it, logstep) == 0) then
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
          !print *, it, i, tip4p(pit,i)%O
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
      if (mod(pit, liststep) == 0) then
        do i=1, no
          ik = 0
          Distancelist(i)%pair = -1
          do j=1, no
            dr = tip4p(pit,i)%O - tip4p(pit,j)%O
            !print *, i,j, dr
            dr = dr - box(pit)*nint(dr/box(pit))
            !print *, i,j, dr
            dist = sqrt(sum(dr**2))
            if (dist <= dlist .and. 0.010 < dist) then
              ik = ik + 1
              DistanceList(i)%pair(ik) = j
              !print *, pit, i, j, ik, dist
            end if
            !print *, pit, i, ik, j, DistanceList(i)%pair(ik)
          end do
          do jk = 1, lmax
            print *, pit, i, jk, DistanceList(i)%pair(jk)
          end do
        end do
      end if

      ! calc distance for i-DistanceList(i)%pair
      do i=1, no
        do jk=1, lmax
          j = DistanceList(i)%pair(jk)
          if( j < 0 ) exit

          dr = tip4p(pit,i)%O - tip4p(pit,j)%O
          dr = wrap_vector(dr, box(pit))
          dist = sqrt(sum(dr**2))

          if (dist >= dhbond ) cycle
          !print *, i, j, dist 
          r_hjoi = wrap_vector( tip4p(pit,j)%H1-tip4p(pit,i)%O, box(pit) )
          r_hioi = wrap_vector( tip4p(pit,i)%H1-tip4p(pit,i)%O, box(pit) )
          u_hjoi = unit_vector(r_hjoi)
          u_hioi = unit_vector(r_hioi)
          !print *, u_hjoi, u_hioi
          cosij = 0.0E0
          do k=1, 3
            cosij = cosij + u_hjoi(k)*u_hioi(k)
          end do
          if( abs(cosij) < dcosij ) cycle
          print *, i, j, cosij

        end do
      end do
    end do


stop

contains
  ! retrun vector  after PBC procedure
  function wrap_vector(rvec, box_dim) result(rvec_new)
    implicit none
    real(8), dimension(3) :: rvec
    real(8) :: box_dim
    real(8), dimension(3) :: rvec_new

    rvec_new = rvec - box_dim*nint(rvec/box_dim)
    
  end function wrap_vector


  ! return unit vector
  function unit_vector(rvec) result(unit_rvec)
    implicit none
    real(8), dimension(3) :: rvec, unit_rvec
    real(8) :: dr_inverese
     
    dr_inverese = 1.0E0/sqrt(sum(rvec**2))
    unit_rvec = dr_inverese * rvec

  end function unit_vector

  
  ! return theta(rad)


end program hbond
