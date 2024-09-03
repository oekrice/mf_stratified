!*******************************************************************************

MODULE evolve
!*******************************************************************************
! Initialise the simulation. Read in/compute the initial condition (python code for that I think?), read in the parameters and all will be well
!*******************************************************************************
    USE shared_data
    USE output
    USE pressure
    IMPLICIT NONE

!*******************************************************************************
CONTAINS

SUBROUTINE timestep()


    !New timestep with a quasi-implicit bit. Using a predictor for the magnfield of the new bit but not the old
    CALL calculate_magnetic()

    CALL calculate_jp()

    CALL calculate_current()

    CALL calculate_velocity()

    if (smart_velocity .and. mod(n,int(0.01_num/dt)) == 0) then
        vout_masks = merge(1.0,0.25,az >= minval(az(0:nx,ny)))
        vout_masks(-1:nx+1,-1:ny+1) = 0.25_num*(vout_masks(0:nx+2,0:ny+2) + vout_masks(0:nx+2,-1:ny+1) + &
        vout_masks(-1:nx+1,0:ny+2) + vout_masks(-1:nx+1,-1:ny+1))
        vout_maskc(-1:nx+2,:) = 0.5_num*(vout_masks(-2:nx+1,:) + vout_masks(-1:nx+2,:))
    end if

    if (mod(n,int(0.01_num/dt)) == 0) then
        buoy_masks = merge(0.0,1.0,az >= minval(az(0:nx,ny)))
        buoy_masks(-1:nx+1,-1:ny+1) = 0.25_num*(buoy_masks(0:nx+2,0:ny+2) + buoy_masks(0:nx+2,-1:ny+1) + &
        buoy_masks(-1:nx+1,0:ny+2) + buoy_masks(-1:nx+1,-1:ny+1))
        buoy_maskc(-1:nx+2,:) = 0.5_num*(buoy_masks(-2:nx+1,:) + buoy_masks(-1:nx+2,:))
    end if


    CALL calculate_electric()


END SUBROUTINE timestep


SUBROUTINE calculate_magnetic()

    !allocate(ax(-1:nx+2,-2:ny+2)); allocate(ay(-2:nx+2,-1:ny+2)); allocate(az(-2:nx+2,-2:ny+2))
    !allocate(bx(-2:nx+2,-1:ny+2)); allocate(by(-1:nx+2,-2:ny+2)); allocate(bz(-1:nx+2,-1:ny+2))

    !Determine the interior points from the vector potential. Some external points (the staggered ones) then determined using the boundary conditions

    ! INTERIOR POINTS
    bx(-2:nx+2, 0:ny+1) =  (az(-2:nx+2, 0:ny+1) - az(-2:nx+2, -1:ny))/dy
    by(0:nx+1, -2:ny+2) = -(az(0:nx+1, -2:ny+2) - az(-1:nx, -2:ny+2))/dx
    bz(0:nx+1, 0:ny+1) = (ay(0:nx+1,0:ny+1) - ay(-1:nx, 0:ny+1))/dx - (ax(0:nx+1,0:ny+1) - ax(0:nx+1, -1:ny))/dy

    ! ZERO CURRENT BOUNDARY CONDITIONS

    bx(-1:nx+1,0)    = bx(-1:nx+1,1)  - dy*(by(0:nx+2,0)  - by(-1:nx+1, 0)) /dx   !Lower boundary: jz = 0 (on grid points)
    bx(-1:nx+1,ny+1) = bx(-1:nx+1,ny) + dy*(by(0:nx+2,ny) - by(-1:nx+1, ny))/dx   !Upper boundary: jz = 0 (on grid points)
    by(0,-1:ny+1)    = by(1,-1:ny+1)  - dx*(bx(0,0:ny+2)  - bx(0, -1:ny+1)) /dy  !Side boundaries
    by(nx+1,-1:ny+1) = by(nx,-1:ny+1) + dx*(bx(nx,0:ny+2) - bx(nx, -1:ny+1))/dy

    ! OTHER BOUNDARY CONDITIONS
    bx(0, :) = 0.0; bx(nx, :) = 0.0   !No flux through the sides
    by(0, :) = by(1, :); by(nx+1,:) = by(nx,:)

    bz(0,:) = bz(1,:); bz(nx+1,:) = bz(nx,:)  !No gradient in the out-of plane field
    bz(:,0) = bz(:,1); bz(:,ny+1) = bz(:,ny)

END SUBROUTINE calculate_magnetic

SUBROUTINE calculate_current()

    !Determine the interior points from the vector potential. Some external points (the staggered ones) then determined using the boundary conditions

    ! ONLY NECESSARY GHOST POINTS - everything else is meaningless and may cause problems
    jx(0:nx+1, 0:ny  ) =   (bz(0:nx+1, 1:ny+1) - bz(0:nx+1, 0:ny))/dy
    jy(0:nx,   0:ny+1) =  -(bz(1:nx+1, 0:ny+1) - bz(0:nx, 0:ny+1))/dx
    jz(0:nx, 0:ny) =  (by(1:nx+1,0:ny) - by(0:nx, 0:ny))/dx - (bx(0:nx,1:ny+1) - bx(0:nx, 0:ny))/dy
END SUBROUTINE calculate_current

SUBROUTINE calculate_velocity()
    !Average magnetic field and current to the gridpoints.
    bx1(0:nx,0:ny) = 0.5*(bx(0:nx,1:ny+1) + bx(0:nx,0:ny))
    by1(0:nx,0:ny) = 0.5*(by(1:nx+1,0:ny) + by(0:nx,0:ny))
    bz1(0:nx,0:ny) = 0.25*(bz(1:nx+1,1:ny+1) + bz(0:nx,1:ny+1) + bz(0:nx,0:ny) + bz(1:nx+1,0:ny))

    jx1(0:nx,0:ny) = 0.5*(jx(1:nx+1,0:ny) + jx(0:nx,0:ny))
    jy1(0:nx,0:ny) = 0.5*(jy(0:nx,1:ny+1) + jy(0:nx,0:ny))
    jz1(0:nx,0:ny) = jz(0:nx, 0:ny) - jpz1(0:nx, 0:ny)
    ! Ensure no currents on the boundaries here (as the averaging takes values that are outside the calculated domain)
    jx1(0,:) = 0.0; jy1(0,:) = 0.0; jz1(0,:) = 0.0
    jx1(nx,:) = 0.0; jy1(nx,:) = 0.0; jz1(nx,:) = 0.0
    jx1(:,0) = 0.0; jy1(:,0) = 0.0; jz1(:,0) = 0.0
    jx1(:,ny) = 0.0; jy1(:,ny) = 0.0; jz1(:,ny) = 0.0

    bb1 = bx1**2 + by1**2 + bz1**2 !B squared

    if (abs(delta) < 1e-10) then !No softening
        nu(0:nx,0:ny) = bb1(0:nx,0:ny)
    else !Softening.
        nu(0:nx,0:ny) = (bb1 (0:nx,0:ny)+ delta*exp(-bb1(0:nx,0:ny)/delta))
    end if
    !Cross product:
    vx(0:nx,0:ny) = nu0_field(0:nx,0:ny)*(jy1(0:nx,0:ny)*bz1(0:nx,0:ny) - jz1(0:nx,0:ny)*by1(0:nx,0:ny))/nu(0:nx,0:ny)
    vy(0:nx,0:ny) = nu0_field(0:nx,0:ny)*(jz1(0:nx,0:ny)*bx1(0:nx,0:ny) - jx1(0:nx,0:ny)*bz1(0:nx,0:ny))/nu(0:nx,0:ny)
    vz(0:nx,0:ny) = nu0_field(0:nx,0:ny)*(jx1(0:nx,0:ny)*by1(0:nx,0:ny) - jy1(0:nx,0:ny)*bx1(0:nx,0:ny))/nu(0:nx,0:ny)

    vx(:,0) = 0.0; vy(:,0) = 0.0; vz(:,0) = 0.0  !Ensure zero velocity directly on the lower boundary.

END SUBROUTINE calculate_velocity

SUBROUTINE calculate_electric()

    ex = 0.0; ey = 0.0; ez = 0.0
    if (eta > 0) then  !There is some coronal diffusion
        ex(1:nx,1:ny) = eta*jx(1:nx,1:ny)
        ey(1:nx,1:ny) = eta*jy(1:nx,1:ny)
        ez(1:nx,1:ny) = eta*(jz(1:nx,1:ny) - jpz1(1:nx,1:ny))
    end if

    if (nu0 > 0) then  !There is some magnetofriction. Use existing velocity field to calculate electric field.
        ex1 = vz*by1 - vy*bz1
        ey1 = vx*bz1 - vz*bx1
        ez1 = vy*bx1 - vx*by1

        !Average to Ribs (interior only):
        ex(1:nx,0:ny) = ex(1:nx,0:ny)  + 0.5*(ex1(0:nx-1,0:ny) + ex1(1:nx,0:ny))
        ey(0:nx,1:ny) = ey(0:nx,1:ny)  + 0.5*(ey1(0:nx,0:ny-1) + ey1(0:nx,1:ny))
        ez(0:nx,0:ny) = ez(0:nx,0:ny) + ez1(0:nx,0:ny)

    end if

    !Add 'buoyant' velocity, assuming that with stratification the rope is heavier than the background and will resist eruption


    !Add outflow
    ex(1:nx, 0:ny) = ex(1:nx, 0:ny) - bz(1:nx, 0:ny)*voutc(1:nx, 0:ny)*vout_maskc(1:nx, 0:ny) !is upwinding here
    ez(0:nx, 0:ny) = ez(0:nx, 0:ny) + bx(0:nx, 0:ny)*vouts(0:nx, 0:ny)*vout_masks(0:nx, 0:ny)

    if (shearfact > 0) then  !Add the shearing on the lower boundary
        ex(1:nx,0) = vz0(1:nx)*by(1:nx, 0)
    end if

    if (eta0 > 0) then  !There is some photospheric diffusion
        ez(0:nx,0) = eta0*(by(1:nx+1,0) - by(0:nx, 0))/dx
    end if


END SUBROUTINE calculate_electric





















!*******************************************************************************
END MODULE evolve
!*******************************************************************************
