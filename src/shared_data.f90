!*******************************************************************************
MODULE shared_data
!*******************************************************************************
! Establish the data that needs to be shared among all the bits of the program
!*******************************************************************************

    IMPLICIT NONE
!*******************************************************************************

    INTEGER, PARAMETER:: d = KIND(0.d0) ! precision for floats
    REAL(d), PARAMETER:: PI = 4.0_d * atan(1.0_d)   ! pi
    INTEGER, PARAMETER :: num = KIND(1.D0)

    !Variable parameters that are read in from the variable text file
    REAL(num):: run_number
    INTEGER:: nx
    REAL(num):: tmax
    REAL(num):: nplots
    REAL(num):: voutfact
    REAL(num):: bfact
    REAL(num):: shearfact
    REAL(num):: eta
    REAL(num):: eta0
    REAL(num):: nu0
    INTEGER:: init_id

    !Other parameters hard-coded into the main.f90 file
    REAL(num):: x0
    REAL(num):: x1
    REAL(num):: y0
    REAL(num):: y1
    REAL(num):: cfl
    REAL(num):: delta
    CHARACTER(LEN = 64) :: output_directory_root
    CHARACTER(LEN = 64) :: output_directory
    INTEGER:: hamilton_flag
    INTEGER:: ndiags
    INTEGER:: decay_type
    !Grid data and arrays shared throughout the code
    INTEGER:: ny
    INTEGER:: nt, n
    REAL(num):: dx, dy
    REAL(num):: t, dt
    REAL(num), DIMENSION(:), ALLOCATABLE :: xs, ys, xc, yc
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: ax, ay, az
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: bx, by, bz
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: jx, jy, jz
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: ex, ey, ez

    REAL(num), DIMENSION(:,:), ALLOCATABLE:: vouts, voutc
    REAL(num), DIMENSION(:), ALLOCATABLE:: vz0

    !Temporary arrays that it's probably just quicker to overwrite each time...
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: bx1, by1, bz1
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: jx1, jy1, jz1
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: ex1, ey1, ez1

    REAL(num), DIMENSION(:,:), ALLOCATABLE :: vx, vy, vz
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: bb1, nu, nu0_field

    !Extra arrays for the pressure bits
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: fy
    REAL(num):: a, b, deltay, ystar
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: vpx, vpy, vpz
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: jpz1
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: vout_masks, vout_maskc
    LOGICAL:: smart_velocity
    REAL(num):: nu0_decay
    !Diagnostics
    REAL(num), DIMENSION(:), ALLOCATABLE:: diag_time
    REAL(num), DIMENSION(:), ALLOCATABLE:: diag_oflux
    REAL(num), DIMENSION(:), ALLOCATABLE:: diag_maxj
    REAL(num), DIMENSION(:), ALLOCATABLE:: diag_avgj
    REAL(num), DIMENSION(:), ALLOCATABLE:: diag_energy
    REAL(num), DIMENSION(:), ALLOCATABLE:: diag_maxlorentz
    REAL(num), DIMENSION(:), ALLOCATABLE:: diag_avglorentz
    REAL(num), DIMENSION(:), ALLOCATABLE:: diag_rflux
    REAL(num), DIMENSION(:), ALLOCATABLE:: diag_rheight

!*******************************************************************************
END MODULE shared_data
!*******************************************************************************
