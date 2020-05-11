program stress
    use xdr, only: trrfile
    implicit none
    type(trrfile) :: trr
    integer, parameter :: nsol  = 2720 
    integer, parameter :: nsol_atm = 4
    integer, parameter :: npoly =   50
    integer, parameter :: npoly_atm = 1405
    integer, parameter :: nstep = 501
    real(kind=8), allocatable :: pos(:,:), force(:,:), sigma(:,:,:)
    real(kind=8) :: sigma_buffer(3,3) = 0.0E0
    integer :: nhist, ng, ig, npos, i, alpha, beta

    call trr % init("md.trr")

    npos = nsol * nsol_atm + npoly * npoly_atm
    allocate(pos(3, npos))
    allocate(force(3, npos))
    allocate(sigma(3, 3, nstep))
    sigma = 0.0E0
    ng = 0
    call trr % read
    
    do while ( trr % STAT == 0 )
        ng = ng + 1
        pos = trr % pos(:, :)  !get the position of molecule
        force = trr % force(:, :)
        do i = 1, npos
            do alpha = 1, 3
                do beta = 1, 3
                    sigma_buffer(alpha, beta) = pos(alpha,i)*force(beta,i)
                    sigma(alpha, beta, ng) = sigma(alpha, beta, ng) + sigma_buffer(alpha, beta)
                end do
            end do
        end do
        call trr % read
    end do

    ! 5. Close the file
    call trr % close

    ! output results
    do i = 1, nstep
      write(*,'(f8.3,2X,f20.9)') i*0.1, sigma(1,1,i)
    end do
end program stress
