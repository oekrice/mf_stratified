!*******************************************************************************
PROGRAM main
!*******************************************************************************
! Main fortran file. Use this to initialise parameters, grid, produce/read in initial conditions and then run the code
!*******************************************************************************
    USE shared_data
    USE init, ONLY: initialise
    USE evolve, ONLY: timestep
    USE output, ONLY: save_snap, diagnostics
    IMPLICIT NONE

    ! Put some of the major variables in here - things that can be changed occasionally but not in a series of runs
    cfl  = 0.1
    delta = 1e-4

    smart_velocity = .true.

    ! Import the parameters and set up the grid
    CALL initialise()

    if (hamilton_flag < 0.5) then
        output_directory_root = '/extra/tmp/trcn27/mf2d/'
    else
        output_directory_root = '/nobackup/trcn27/mf2d0/'
    end if

    print*, 'Initial condition set up in Fortran. Running...'
    do n = 0, nt-1  ! Actually run the code

        CALL timestep()  !Does everything except the actual timestep (for diagnostic reasons)

        if (MOD(n, (nt/int(nplots-1))) == 0) then   ! Save a snapshot (prints a message as well)
            CALL save_snap(int(n/(nt/int(nplots-1))))

        end if

        if (MOD(n, (nt/int(ndiags-1))) == 0) then   ! Save a snapshot (prints a message as well)
            CALL diagnostics(int(n/(nt/(ndiags-1))))
            print*, 'Step', n, 'at time', t, 'Open flux:', diag_oflux(int(n/(nt/(ndiags-1))))

        end if

        ax = ax - dt*ex
        ay = ay - dt*ey
        az = az - dt*ez

        t = t + dt

    end do

    CALL save_snap(int(nplots-1.0))
    CALL diagnostics(ndiags-1)

    print*, 'Fortran code completed sucessfully. Carry on.'

END PROGRAM main
