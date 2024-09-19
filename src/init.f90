!*******************************************************************************
MODULE init
!*******************************************************************************
! Initialise the simulation. Read in/compute the initial condition (python code for that I think?), read in the parameters and all will be well
!*******************************************************************************
    USE shared_data
    USE output
    USE pressure
    IMPLICIT NONE

!*******************************************************************************

CONTAINS

SUBROUTINE initialise()

    CALL read_parameters()

    CALL establish_grid()

    CALL initial_condition()

    CALL calculate_timestep()

    CALL set_outflow()

    CALL set_shearing()

    CALL set_friction()

    CALL pressure_function()

END SUBROUTINE initialise

SUBROUTINE read_parameters()
    !Reads in the parameters from the text file variables.txt
    REAL(num), DIMENSION(30):: variables
    CHARACTER(LEN=64):: input_value

    CHARACTER(LEN=64):: parameter_filename
    call get_command_argument(1, input_value)
    read(unit=input_value,fmt=*) run_number

    print*, 'Magnetofriction run number', int(run_number)

    if (run_number < 10) then
      write (parameter_filename, "(A22, A2, I1, A4)") './parameters/variables', '00', int(run_number), '.txt'
    else if (run_number < 100) then
      write (parameter_filename, "(A22, A1, I2, A4)") './parameters/variables', '0', int(run_number), '.txt'
    else
      write (parameter_filename, "(A22, I3, A4)") './parameters/variables', int(run_number), '.txt'
    end if

    OPEN(1, FILE = parameter_filename)
    READ(1, *) variables
    CLOSE(1)

    run_number = variables(1)
    nx = int(variables(2))
    ny = int(variables(21))

    tmax = variables(3)
    nplots = int(variables(4))
    voutfact = variables(5)
    bfact = variables(6)
    shearfact = variables(7)
    eta = variables(8)
    nu0 = variables(9)
    eta0 = variables(10)

    x0 = variables(11)
    x1 = variables(12)
    y0 = variables(13)
    y1 = variables(14)

    a = variables(15)
    b = variables(16)
    deltay = variables(17)
    ystar = variables(18)

    hamilton_flag = int(variables(19))
    ndiags = int(variables(20))

    init_id = nx
    nu0_decay = variables(23)

    decay_type = int(variables(24))

    buoyant_factor = variables(25)

    print*, 'Importing parameters from file', init_id
    print*, 'Output directory', output_directory

    print*, 'eta0 = ', eta0
    print*, 'nu0 = ', nu0

    print*, 'ystar = ', ystar
    print*, 'a = ', a
    print*, 'b = ', b
    print*, 'deltay = ', deltay

    if (decay_type == 0) print*, 'No pressure'
    if (decay_type == 1) print*, 'Exponential pressure decay'
    if (decay_type == 2) print*, 'Smooth tanh pressure decay'
    if (decay_type == 3) print*, 'Sharp tanh pressure decay'


END SUBROUTINE read_parameters

SUBROUTINE establish_grid()
    !Estalish the grid (not particularly difficult for a cartesian box). Try to use the same numbering system as LARE.
    !index 0 is always the cell on the boundary, or just OUTSIDE it. c cells are numbered from -1, s cells from -2
    !Also declare the sizes of the variable arrays (A, B etc.) which should only be multiples of the axes (no subsets, like in the original code)
    INTEGER:: i, j

    dx = (x1-x0)/(nx) ;  dy = (y1-y0)/(ny)
    t = 0.0

    allocate(xs(-2:nx+2)); allocate(ys(-2:ny+2))
    allocate(xc(-1:nx+2)); allocate(yc(-1:ny+2))

    !Allocate global arrays (things that could conceivably be plotted afterwards, I suppose)
    allocate(ax(-1:nx+2,-2:ny+2)); allocate(ay(-2:nx+2,-1:ny+2)); allocate(az(-2:nx+2,-2:ny+2))
    allocate(bx(-2:nx+2,-1:ny+2)); allocate(by(-1:nx+2,-2:ny+2)); allocate(bz(-1:nx+2,-1:ny+2))
    allocate(jx(-1:nx+2,-2:ny+2)); allocate(jy(-2:nx+2,-1:ny+2)); allocate(jz(-2:nx+2,-2:ny+2))
    allocate(ex(-1:nx+2,-2:ny+2)); allocate(ey(-2:nx+2,-1:ny+2)); allocate(ez(-2:nx+2,-2:ny+2))

    allocate(bx1(-2:nx+2,-2:ny+2)); allocate(by1(-2:nx+2,-2:ny+2)); allocate(bz1(-2:nx+2,-2:ny+2))
    allocate(jx1(-2:nx+2,-2:ny+2)); allocate(jy1(-2:nx+2,-2:ny+2)); allocate(jz1(-2:nx+2,-2:ny+2))
    allocate(vx(-2:nx+2,-2:ny+2)); allocate(vy(-2:nx+2,-2:ny+2)); allocate(vz(-2:nx+2,-2:ny+2))
    allocate(ex1(-2:nx+2,-2:ny+2)); allocate(ey1(-2:nx+2,-2:ny+2)); allocate(ez1(-2:nx+2,-2:ny+2))
    allocate(bb1(-2:nx+2,-2:ny+2)); allocate(nu(-2:nx+2,-2:ny+2))

    allocate(fy(0:nx+1,0:ny))
    allocate(vpx(-2:nx+2,-2:ny+2)); allocate(vpy(-2:nx+2,-2:ny+2)); allocate(vpz(-2:nx+2,-2:ny+2))
    allocate(jpz1(-2:nx+2,-2:ny+2))

    allocate(vout_masks(-2:nx+2,-2:ny+2));  allocate(vout_maskc(-1:nx+2,-2:ny+2))
    allocate(buoy_masks(-2:nx+2,-2:ny+2));  allocate(buoy_maskc(-1:nx+2,-2:ny+2))

    ! Establish the actual grid
    do i = -2, nx + 2
        xs(i) = x0 + (real(i)/real(nx))*(x1 - x0)
    end do
    do j = -2, ny + 2
        ys(j) = y0 + (real(j)/real(ny))*(y1 - y0)
    end do

    xc(-1:nx+2) = 0.5*(xs(-2:nx+1) + xs(-1:nx+2))
    yc(-1:ny+2) = 0.5*(ys(-2:ny+1) + ys(-1:ny+2))

    vout_masks = 1.0_num; vout_maskc = 1.0_num

END SUBROUTINE establish_grid

SUBROUTINE initial_condition()

    REAL(num), DIMENSION(0:nx, 0:ny):: az_import
    CHARACTER(LEN =64):: input_directory

    if (run_number < 10) then
    write (input_directory, "(A12, A2, I1, A4)") './inits/init', '00', int(run_number), '.txt'
    else if (run_number < 100) then
    write (input_directory, "(A12, A1, I2, A4)") './inits/init', '0', int(run_number), '.txt'
    else
    write (input_directory, "(A12, I3, A4)") './inits/init', int(run_number), '.txt'
    end if

    open(1,file = input_directory)  !Load in the initial condition in az.
    read(1,*) az_import

    ax(:,:) = 0.0
    ay(:,:) = 0.0
    az(0:nx,0:ny) = -az_import

END SUBROUTINE initial_condition

SUBROUTINE calculate_timestep()
    ! Calculates the timestep based on the parameters and cfl condition. Then adjusts tso there are an integer number of steps in between snapshots, as that's neater
    REAL(num):: dt_ideal
    REAL(num):: plot_dt

    dt_ideal = 1e10 !Maximum time, will get smaller based on the conditions

    if (nu0 > 0) dt_ideal = (min(dx,dy))**2/(nu0*(1.0))
    print*, 'dt due to nu0', cfl*(min(dx,dy))**2/(nu0*(1.0))
    if (voutfact > 0) dt_ideal = min(dt_ideal, dy/voutfact)
    print*, 'dt due to outflow',  cfl*dy/voutfact
    if (eta > 0 ) dt_ideal = min(dt_ideal, min(dx,dy)**2/eta)
    print*, 'dt due to eta', cfl*min(dx,dy)**2/eta
    if (shearfact > 0 ) dt_ideal = min(dt_ideal, dy/shearfact)
    print*, 'dt due to shearing', cfl*dy/shearfact
    if (eta0 > 0) dt_ideal = min(dt_ideal, dx**2/eta0)
    print*, 'dt due to eta0',  cfl*dx**2/eta0
    dt_ideal = cfl*dt_ideal

    !Adjust so there are integer number of timesteps in between plots
    plot_dt = tmax/(ndiags-1)
    nt = int(ndiags-1)*(int((plot_dt-1d-6)/dt_ideal)+1)
    dt = tmax/float(nt)
    print*, 'Final dt', dt, ', total timesteps', nt, ', ', int(nt/(nplots-1)), 'per snapshot'

END SUBROUTINE calculate_timestep

SUBROUTINE set_outflow_old()
    ! Set the outflow arrays vout and voutc, the sizes of which are declared in shared variables
    ! vouts is on gridpoints (like az), voutc is on x ribs (like ax). To allow for upwinding
    ! Need to take care with all this...
    INTEGER:: i, j
    REAL(num):: hfact

    allocate(vouts(-2:nx+2,-2:ny+2))
    allocate(voutc(-1:nx+2,-2:ny+2))

    do j = 0, ny
        hfact = (ys(j) - ys(0))/(ys(ny) - ys(0))   !Distance up the domain
        do i = -2, nx+2
            vouts(i,j) = voutfact*hfact**10
        end do
        do i = -1, nx+2
            voutc(i,j) = voutfact*hfact**10
        end do
    end do

END SUBROUTINE set_outflow_old

SUBROUTINE set_outflow()
    ! Set the outflow arrays vout and voutc, the sizes of which are declared in shared variables
    ! vouts is on gridpoints (like az), voutc is on x ribs (like ax). To allow for upwinding
    ! Need to take care with all this...
    CHARACTER(LEN =64):: input_directory
    INTEGER:: j
    allocate(vouts(-2:nx+2,-2:ny+2))
    allocate(voutc(-1:nx+2,-2:ny+2))

    if (.false.) then  !Import outflow function from initial LARE version
        if (init_id < 10) then
        write (input_directory, "(A10, A2, I1, A4)") './inits/vy', '00', init_id, '.txt'
        else if (init_id < 100) then
        write (input_directory, "(A10, A1, I2, A4)") './inits/vy', '0', init_id, '.txt'
        else
        write (input_directory, "(A10, I3, A4)") './inits/vy', init_id, '.txt'
        end if

        open(1,file = input_directory)
        read(1,*) vouts(-2:nx+2,-2:ny+2)
    else        !Instead define as an explicit function
        do j = -2, ny+2
            vouts(:,j) = (ys(j)/ys(ny))**4
        end do
    end if
    voutc(-1:nx+2,-2:ny+2) = 0.5*(vouts(-2:nx+1,-2:ny+2)+vouts(-1:nx+2,-2:ny+2))

END SUBROUTINE set_outflow

SUBROUTINE set_friction()
    !Sets the nu0 field, which can now vary with height (or some other parameter)
    !To more closely match the MHD results
    INTEGER:: j
    allocate(nu0_field(-2:nx+2,-2:ny+2))

    do j = -2, ny+2
        nu0_field(-2:nx+2,j) = nu0*exp(-ys(j)*nu0_decay)
    end do

END SUBROUTINE set_friction

SUBROUTINE set_shearing()
    ! Set the 1D array for the shearing velocity on the lower boundary. Is added to ex so is aligned with the x grid centres
    INTEGER:: i
    REAL(num) :: ref_speed
    allocate(vz0(-1:nx+2))


    ref_speed = -cos(0.35)*(0.0031 - 0.041*cos(0.35)**2 - 0.031*cos(0.35)**4)

    do i=-1,nx+2
      vz0(i) = -cos(0.35-xc(i))*(0.0031 - 0.041*cos(0.35-xc(i))**2 - 0.031*cos(0.35-xc(i))**4) - ref_speed
    end do

    !vz0(-1:nx+2) = vz0(-1:nx+2)*(1.0_num - sin(pi*xc(-1:nx+2)/(2*xc(0)))**100)
    vz0 = vz0*shearfact

END SUBROUTINE set_shearing

!*******************************************************************************
END MODULE init
!*******************************************************************************
