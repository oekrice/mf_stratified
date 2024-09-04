!*******************************************************************************

MODULE pressure
!*******************************************************************************
! Functions for the new pressure part
! Gives a function to add on to the velocity term, which is averaged to grid POINTS
!

!*******************************************************************************
    USE shared_data

    IMPLICIT NONE

!*******************************************************************************
CONTAINS

SUBROUTINE pressure_function()

    !Input the pressure function f(y) here - note y is the vertical coordinate like LARE
    !Should this be positive or negative?!?!
    implicit none
    integer:: j

    do j = 0, ny
        if (decay_type < 0.5) then
            fy(:,j) = a*(1.0_num - b*tanh((ys(j)-ystar)/deltay))
        else
            fy(:,j) = a*exp(-ys(j)/b)
        end if
    end do

END SUBROUTINE pressure_function

SUBROUTINE calculate_pressure()
    !Does the cross product with b, averages to gridpoints and does the softening as for the velocity
    implicit none

    REAL(num), DIMENSION(0:nx,0:ny):: gx2, gy2

    !Can get nu directly from the velocity calculation, so don't need to do this twice. But this function must come after calculate_velocity
    gx2(0:nx,0:ny) = -nu0*(jpz1(0:nx,0:ny)*by1(0:nx,0:ny))/nu(0:nx,0:ny)
    gy2(0:nx,0:ny) = nu0*(jpz1(0:nx,0:ny)*bx1(0:nx,0:ny))/nu(0:nx,0:ny)

    vx(0:nx,0:ny) = vx(0:nx,0:ny) + gx2(0:nx, 0:ny)
    vy(0:nx,0:ny) = vy(0:nx,0:ny) + gy2(0:nx, 0:ny)

    vx(:,0) = 0.0; vy(:,0) = 0.0

END SUBROUTINE calculate_pressure

SUBROUTINE calculate_jp()
    !Calculates jp - the 'pressure current'
    implicit none

    jpz1(0:nx, 0:ny) =  (fy(1:nx+1,0:ny)*by(1:nx+1,0:ny) - fy(0:nx,0:ny)*by(0:nx, 0:ny))/dx

END SUBROUTINE calculate_jp

SUBROUTINE pressure_minimiser()
    !Calculates the functional to be minimised for MHD (normal veclotiy plus the new bit)

    bb1 = bx1**2 + by1**2 + bz1**2 !B squared

    !Cross product:
    vpx(0:nx,0:ny) = nu0*(jy1(0:nx,0:ny)*bz1(0:nx,0:ny) - (jz1(0:nx,0:ny)-jpz1(0:nx,0:ny))*by1(0:nx,0:ny))/nu(0:nx,0:ny)
    vpy(0:nx,0:ny) = nu0*((jz1(0:nx,0:ny)-jpz1(0:nx,0:ny))*bx1(0:nx,0:ny) - jx1(0:nx,0:ny)*bz1(0:nx,0:ny))/nu(0:nx,0:ny)
    vpz(0:nx,0:ny) = nu0*(jx1(0:nx,0:ny)*by1(0:nx,0:ny) - jy1(0:nx,0:ny)*bx1(0:nx,0:ny))/nu(0:nx,0:ny)

    vpx(:,0) = 0.0; vpy(:,0) = 0.0; vpz(:,0) = 0.0  !Ensure zero velocity directly on the lower boundary.

END SUBROUTINE pressure_minimiser



!*******************************************************************************
END MODULE pressure
!*******************************************************************************
