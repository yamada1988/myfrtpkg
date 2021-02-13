!  Requirement:libxdrfile-1.1.4
!  Download from https://github.com/wesbarnett/libxdrfile
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
    double precision :: box
    integer :: it, pit, pit0, ptotalstep
    integer :: i,j,k, jo, jh1, jh2, ik, jk
    character outf*128


    type mytype
    real(8), dimension(3):: O
    real(8), dimension(3):: H1
    real(8), dimension(3):: H2
    end type mytype

    type(mytype), allocatable, dimension(:) :: tip3p


    type mytype2
    integer, dimension(200) :: pair = -1
    end type mytype2

    type(mytype2), allocatable, dimension(:) :: DistanceList


    real(8), allocatable, dimension(:) :: sm_val
    integer, allocatable, dimension(:) :: sm_indx, sm_jndx


    integer, parameter :: nmlfu=20
    integer, parameter :: hbond_io=17
    real(8), parameter :: PI=3.1415927
    real(8) :: dist, dcosij
    real(8), dimension(2) :: cosij, cosji
    real(8), dimension(3) :: r_ojoi, r_hi1oi, r_hi2oi, r_oioj, r_hj1oj, r_hj2oj
    character(len=200) :: inpfile, outfile
    integer :: totalstep, skipstep, logstep, liststep, lmax, i_hbond
    real(8) :: dlist, dhbond, dtheta
    namelist/indata/inpfile,outfile
    namelist/inparam/totalstep,skipstep,logstep,liststep,dlist,lmax,dhbond,dtheta
    character(len=20) :: hfile = "hbond.sm"
    

    real(8), allocatable, dimension(:) :: binary_val 
    integer, allocatable, dimension(:) :: binary_indx, binary_jndx    
    integer :: ih


    ! 2. Declare a variable of type xtcfile
    type(xtcfile) :: xtcf


    ! 2.1 Use namelist
    open(unit=nmlfu, file='namelist')
    read(unit=nmlfu, nml=indata)
    read(unit=nmlfu, nml=inparam)
    close(unit=nmlfu)

    
    ! 2.1.2 delete hbond.sm
    open(hbond_io, file=hfile)
    close(hbond_io, status='delete')


    ! 3. Initialize it with the names of xtc files you want to read in and write out
    call xtcf % init(inpfile)

    natm = xtcf % NATOMS
    no = natm / 3
    ptotalstep = totalstep / skipstep
    dcosij = cos(dtheta*PI/180.0)
    allocate(tip3p(1:no))
    allocate(DistanceList(1:no))
    allocate(sm_val(10*no))
    allocate(sm_indx(10*no))
    allocate(sm_jndx(10*no)) 


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
        box = xtcf % box(1,1)

        !$omp parallel shared(xtcf, tip3p)
        !$omp private(jo, jh1, jh2)
        !$omp do
        ! get position (j-th O) and box at it-th step
        do i=1, no
          jo = 3*(i-1) + 1
          jh1 = jo + 1
          jh2 = jo + 2
          tip3p(i) % O = xtcf % pos(:, jo)
          tip3p(i) % H1 = xtcf % pos(:, jh1)
          tip3p(i) % H2 = xtcf % pos(:, jh2)
          !print *, it, i, tip3p(i)%O
        end do
        !$omp end do
        !$omp end parallel

        ! Initialize sparse matrix
        sm_val = -1.0E0
        sm_indx = -1
        sm_jndx = -1

        ! make list
        if (mod(pit, liststep) == 0) then
          do i=1, no
            ik = 0
            Distancelist(i)%pair = int(-1)
            do j=1, no
              r_ojoi = tip3p(i)%O - tip3p(j)%O
              r_ojoi = wrap_vector(r_ojoi, box)
              dist = sqrt(sum(r_ojoi**2))
              if (dist <= dlist .and. 0.010 < dist) then
                ik = ik + 1
                DistanceList(i)%pair(ik) = j
              end if
            end do
          end do
        end if


      ! calc hbond-pair for i-DistanceList(i)%pair 
      ! ref:Kumar et.al., J. Chem. Phys.126, 204107, (2007)
      i_hbond = 1
      do i=1, no
        do jk=1, lmax
          j = DistanceList(i)%pair(jk)
          if( j < 0 ) exit

          r_ojoi = tip3p(i)%O - tip3p(j)%O
          r_ojoi = wrap_vector(r_ojoi, box)
          dist = sqrt(sum(r_ojoi**2))

          if (dist >= dhbond ) cycle
          ! calculate angle between oi-oj and oi-hki(k=1,2)
          r_hi1oi = wrap_vector( tip3p(i)%H1-tip3p(i)%O, box )
          r_hi2oi = wrap_vector( tip3p(i)%H2-tip3p(i)%O, box )

          r_ojoi = r_ojoi/dist
          r_hi1oi = unit_vector(r_hi1oi)
          r_hi2oi = unit_vector(r_hi2oi)
          !print *, u_hjoi, u_hioi

          cosij(1) = sum(r_ojoi(:)*r_hi1oi(:))
          cosij(2) = sum(r_ojoi(:)*r_hi2oi(:))

          ! calculate angle between oj-oi and oj-hkj(k=1,2)
          r_hj1oj = wrap_vector( tip3p(j)%H1-tip3p(j)%O, box )
          r_hj2oj = wrap_vector( tip3p(j)%H2-tip3p(j)%O, box )
          r_oioj = -1.0E0*r_ojoi
          r_hj1oj = unit_vector(r_hj1oj)
          r_hj2oj = unit_vector(r_hj2oj)
          cosji(1) = sum(r_oioj(:)*r_hj1oj(:))
          cosji(2) = sum(r_oioj(:)*r_hj2oj(:))


          if ( abs(cosij(1)) > dcosij .or. abs(cosij(2)) > dcosij .or. abs(cosji(1)) > dcosij .or. abs(cosji(2)) > dcosij) then
            sm_val(i_hbond) = 1.0E0
            sm_indx(i_hbond) = int(i)
            sm_jndx(i_hbond) = int(j)
            i_hbond = i_hbond + 1
          end if
        end do
      end do

      ! allocate binary array
      allocate(binary_val(i_hbond-1))
      allocate(binary_indx(i_hbond-1))
      allocate(binary_jndx(i_hbond-1)) 

      do ih=1, i_hbond-1
        binary_val(ih) = sm_val(ih)
        binary_indx(ih) = int(sm_indx(ih))
        binary_jndx(ih) = int(sm_jndx(ih))
      end do
      open(unit = hbond_io, file = hfile, form='UNFORMATTED', action = 'write', position='append')
      write(hbond_io) i_hbond-1
      write(hbond_io) binary_val, binary_indx, binary_jndx
      print *, binary_val(1), binary_indx(1), binary_jndx(1)
      close(hbond_io)

      ! deallocate binary array
      deallocate(binary_val)
      deallocate(binary_indx)
      deallocate(binary_jndx)

      ! call nextstep
      call xtcf % read
      it = it + 1
      pit = pit + 1

    end do
    ! 5. Close the file
    call xtcf % close


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

  

end program hbond
