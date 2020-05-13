! TODO:
! See the Virial calculation page on Gromacs manual.
! periodicity must be considered for the virial calculation.

program stress
    use xdr, only: trrfile
    implicit none
    type(trrfile) :: trr
    integer, parameter :: nsol  = 2720 
    integer, parameter :: nsol_atm = 4
    integer, parameter :: npoly =   50
    integer, parameter :: npoly_atm = 1405
    integer, parameter :: nstep = 10001
    real(kind=8), parameter :: dt = 0.010E0, sfctr = 16.6057788
    real(kind=8), allocatable :: pos(:,:), vel(:,:), force(:,:), sigma_v(:,:,:,:), sigma_k(:,:,:,:), sigma(:,:,:,:), amass(:)
    real(kind=8) :: sigma_buffer_v(3,3) = 0.0E0, sigma_buffer_k(3,3) = 0.0E0 
    real(kind=8) :: box_dim, volume, amass_buffer
    integer :: nhist, ng, ig, npos, i, alpha, beta, ioerr, offset, dummy1

    character(len=2) :: atm
    character(len=30) :: dummy2
    real, parameter :: massH = 1.00794          ! atomic weight (hydrogen) (g/mol)
    real, parameter :: massC = 12.0107          ! atomic weight (carbon)   (g/mol)
    real, parameter :: massO = 15.9994          ! atomic weight (oxygen)   (g/mol)
    real, parameter :: massM = 0

    call trr % init("../MD/md.trr")

    npos = nsol * nsol_atm + npoly * npoly_atm
    offset = nsol * nsol_atm
    allocate(pos(3, npos))
    allocate(vel(3, npos))
    allocate(force(3, npos))
    allocate(sigma_v(2, 3, 3, nstep))
    allocate(sigma_k(2, 3, 3, nstep))
    allocate(sigma(2, 3, 3, nstep))
    allocate(amass(0:npos))
    amass = 0.0E0
    sigma_v = 0.0E0
    sigma_k = 0.0E0
    sigma   = 0.0E0
    ng = 0

    ! read atom mass
    open (71,file='SltInfo',status='old',iostat=ioerr)
     do i = 1, nsol_atm
        read (71,*) dummy1,atm,dummy2
        select case(atm)
        case('H')
        amass(i) = massH
        case('C')
        amass(i) = massC
        case('O')
        amass(i) = massO
        case('M')
        amass(i) = massM
        case default
        stop "Unknown atom type"
        end select
     end do
    close (71)
    amass(0) = amass(nsol_atm)
    do i = nsol_atm+1, nsol*nsol_atm
       amass(i) = amass(mod(i, nsol_atm))
    end do

    open (72,file='MolInfo',status='old',iostat=ioerr)
     do i = offset+1, offset+npoly_atm
        read (72,*) dummy1,atm,dummy2
        select case(atm)
        case('H')
        amass(i) = massH
        case('C')
        amass(i) = massC
        case('O')
        amass(i) = massO
        case('M')
        amass(i) = massM
        case default
        stop "Unknown atom type"
        end select
     end do
    close (72)
    amass_buffer = amass(offset)
    amass(offset) = amass(offset+npoly_atm)
    do i = offset+npoly_atm+1, npos
       amass(i) = amass(offset+mod((i-offset), npoly_atm))
    end do
    amass(offset) = amass_buffer
    !do i = 1, npos
    !   write(*,'(I7,2X,f8.5)') i, amass(i)
    !end do


    call trr % read
    
    do while ( trr % STAT == 0 )
        ng = ng + 1
        box_dim = trr % box(1,1)
        volume = box_dim**(3.0E0) 
        pos = trr % pos(:, :)     !get the position of atom
        force = trr % force(:, :) !get the force of atom
        vel = trr % vel(:, :)     !get the velocity of atom
        do alpha = 1, 3
           do beta = 1, 3
               ! 1st component
               do i = 1, offset
                    sigma_buffer_v(alpha, beta) =  pos(alpha,i)*force(beta,i)/volume
                    sigma_buffer_k(alpha, beta) =  amass(i)*vel(alpha,i)*vel(beta,i)/volume
                    sigma_v(1, alpha, beta, ng) = sigma_v(1, alpha, beta, ng) + sigma_buffer_v(alpha, beta)
                    sigma_k(1, alpha, beta, ng) = sigma_k(1, alpha, beta, ng) + sigma_buffer_k(alpha, beta)
               end do
               sigma_v(1, alpha, beta, ng) = sigma_v(1, alpha, beta, ng) * sfctr  !bar
               sigma_k(1, alpha, beta, ng) = sigma_k(1, alpha, beta, ng) * sfctr  !bar
               sigma(1, alpha, beta, ng) = sigma_v(1, alpha, beta, ng) + sigma_k(1, alpha, beta, ng)
               
               ! 2nd component
               do i = offset+1, npos
                    sigma_buffer_v(alpha, beta) = pos(alpha,i)*force(beta,i)/volume
                    sigma_buffer_k(alpha, beta) = amass(i)*vel(alpha,i)*vel(beta,i)/volume
                    sigma_v(2, alpha, beta, ng) = sigma_v(2, alpha, beta, ng) + sigma_buffer_v(alpha, beta)
                    sigma_k(2, alpha, beta, ng) = sigma_k(2, alpha, beta, ng) + sigma_buffer_k(alpha, beta)
               end do
               sigma_v(2, alpha, beta, ng) = sigma_v(2, alpha, beta, ng) * sfctr  !bar
               sigma_k(2, alpha, beta, ng) = sigma_k(2, alpha, beta, ng) * sfctr  !bar
               sigma(2, alpha, beta, ng) = sigma_v(2, alpha, beta, ng) + sigma_k(2, alpha, beta, ng)
           end do
        end do
        call trr % read
    end do

    ! 5. Close the file
    call trr % close


    ! output results
    open(17, file='sigma-water.xvg',status='replace')
    do i = 1, nstep
      write(17,'(f8.3,2X,f20.9)') (i-1)*dt, sigma(1, 1, 2, i)
    end do
    close(17)

    open(18, file='sigma-pva.xvg',status='replace')
    do i = 1, nstep
      write(18,'(f8.3,2X,f20.9)') (i-1)*dt, sigma(2, 1, 2, i)
    end do
    close(18)

end program stress
