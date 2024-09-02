!*******************************************************************************
MODULE output
!*******************************************************************************
! Contains tools for saving and printing data to the correct directories etc. Put separately jsut to stop things being messy
!*******************************************************************************
    USE shared_data
    USE netcdf

    IMPLICIT NONE

!*******************************************************************************

CONTAINS

SUBROUTINE diagnostics(diag_num)
    !Calculates some diagnostics and saves to netcdf file as for the triangle code, which was fairly neat (if i say so myself...). Should make for easy pythonning
    IMPLICIT NONE
    INTEGER:: diag_num
    REAL(num):: diag
    character(len=100):: filename

    integer:: id_1, id_2, id_3, id_4, id_5, id_6, id_7, id_8, id_9, ncid, nd_id
    integer:: j

    real(num), dimension(:,:):: jx0(0:nx-1,0:ny-1),jy0(0:nx-1,0:ny-1),jz0(0:nx-1,0:ny-1) !Current averaged to grid centres
    real(num), dimension(:,:):: j0(0:nx-1,0:ny-1)

    real(num), dimension(:,:):: bx0(0:nx-1,0:ny-1),by0(0:nx-1,0:ny-1),bz0(0:nx-1,0:ny-1) !Current averaged to grid centres
    real(num), dimension(:,:):: b0(0:nx-1,0:ny-1)

    real(num), dimension(:,:):: lx0(0:nx-1,0:ny-1), ly0(0:nx-1,0:ny-1),lz0(0:nx-1,0:ny-1) !Current averaged to grid centres
    real(num), dimension(:,:):: l0(0:nx-1,0:ny-1)

    real(num), dimension(:,:):: az_slice(-2:ny+2)

    !Allocate diagnostic arrays
    if (diag_num == 0) then
        allocate(diag_time(0:ndiags-1))
        allocate(diag_oflux(0:ndiags-1)); allocate(diag_maxj(0:ndiags-1))
        allocate(diag_avgj(0:ndiags-1)); allocate(diag_energy(0:ndiags-1))
        allocate(diag_maxlorentz(0:ndiags-1)); allocate(diag_avglorentz(0:ndiags-1))
        allocate(diag_rflux(0:ndiags-1)); allocate(diag_rheight(0:ndiags-1))
        diag_oflux = 0.0_num; diag_maxj = 0.0_num; diag_avgj = 0.0_num; diag_energy = 0.0_num
        diag_maxlorentz = 0.0_num; diag_avglorentz = 0.0_num; diag_rflux = 0.0_num
        diag_rheight = 0.0_num
    end if

    !TIME
    diag = t
    diag_time(diag_num) = t

    !OPEN FLUX
    diag = sum(abs(by(0:nx-1,ny)))*dx
    diag_oflux(diag_num) = diag

    !CURRENT THINGS
    jx0 = 0.5*(jx(0:nx-1,0:ny-1) + jx(0:nx-1,1:ny))
    jy0 = 0.5*(jy(0:nx-1,0:ny-1) + jy(1:nx,0:ny-1))
    jz0 = 0.25*(jz(0:nx-1,0:ny-1) + jz(1:nx,0:ny-1) + jz(0:nx-1,1:ny) + jz(1:nx,1:ny))

    j0 = jx0**2 + jy0**2 + jz0**2

    diag = maxval(sqrt(j0))
    diag_maxj(diag_num) = diag

    diag = sqrt(sum(j0))/float(nx*ny)
    diag_avgj(diag_num) = diag

    !MAGNETIC FIELD THINGS
    bx0 = 0.5*(bx(0:nx-1,0:ny-1) + bx(1:nx,0:ny-1))
    by0 = 0.5*(by(0:nx-1,0:ny-1) + by(0:nx-1,1:ny))
    bz0 = bz(0:nx-1,0:ny-1)

    b0 = bx0**2 + by0**2 + bz0**2

    diag = 0.5*sum(b0)*dx*dy
    diag_energy(diag_num) = diag

    lx0(0:nx-1,0:ny-1) = (jy0(0:nx-1,0:ny-1)*bz0(0:nx-1,0:ny-1) - jz0(0:nx-1,0:ny-1)*by0(0:nx-1,0:ny-1))
    ly0(0:nx-1,0:ny-1) = (jz0(0:nx-1,0:ny-1)*bx0(0:nx-1,0:ny-1) - jx0(0:nx-1,0:ny-1)*bz0(0:nx-1,0:ny-1))
    lz0(0:nx-1,0:ny-1) = (jx0(0:nx-1,0:ny-1)*by0(0:nx-1,0:ny-1) - jy0(0:nx-1,0:ny-1)*bx0(0:nx-1,0:ny-1))

    l0 = lx0**2 + ly0**2 + lz0**2

    diag = maxval(sqrt(l0))
    diag_maxlorentz(diag_num) = diag

    diag = sqrt(sum(l0))/float(nx*ny)
    diag_avglorentz(diag_num) = diag

    !FLUX ROPE THINGS
    az_slice(-2:ny+2) = -az(nx/2,-2:ny+2)

    diag = 0.0_num
    do j = 0, ny
        diag = maxval(az_slice) - max(az_slice(0), az_slice(ny))
    end do
    diag_rflux(diag_num) = diag

    diag = ys(maxloc(az_slice,1)-2)
    diag_rheight(diag_num) = diag
    !ADMIN

    if (run_number < 10) then
        write (filename, "(A19,I1,A3)") "./diagnostics/run00", int(run_number), ".nc"
    else if (run_number < 100) then
        write (filename, "(A18,I2,A3)") "./diagnostics/run0", int(run_number), ".nc"
    else
        write (filename, "(A17,I3,A3)") "./diagnostics/run", int(run_number), ".nc"
    end if

    !Write to diagnostics file, using netcdf
    call try(nf90_create(trim(filename), nf90_clobber, ncid))
    call try(nf90_def_dim(ncid, 'ndiags', ndiags, nd_id))  !Make up fake dimensions here

    call try(nf90_def_var(ncid, 'time', nf90_double, (/nd_id/), id_1))
    call try(nf90_def_var(ncid, 'openflux', nf90_double, (/nd_id/), id_2))
    call try(nf90_def_var(ncid, 'maxcurrent', nf90_double, (/nd_id/), id_3))
    call try(nf90_def_var(ncid, 'avgcurrent', nf90_double, (/nd_id/), id_4))
    call try(nf90_def_var(ncid, 'energy', nf90_double, (/nd_id/), id_5))
    call try(nf90_def_var(ncid, 'avglorentz', nf90_double, (/nd_id/), id_6))
    call try(nf90_def_var(ncid, 'maxlorentz', nf90_double, (/nd_id/), id_7))
    call try(nf90_def_var(ncid, 'rflux', nf90_double, (/nd_id/), id_8))
    call try(nf90_def_var(ncid, 'rheight', nf90_double, (/nd_id/), id_9))

    call try(nf90_enddef(ncid))

    call try(nf90_put_var(ncid, id_1, diag_time))
    call try(nf90_put_var(ncid, id_2, diag_oflux))
    call try(nf90_put_var(ncid, id_3, diag_maxj))
    call try(nf90_put_var(ncid, id_4, diag_avgj))
    call try(nf90_put_var(ncid, id_5, diag_energy))
    call try(nf90_put_var(ncid, id_6, diag_avglorentz))
    call try(nf90_put_var(ncid, id_7, diag_maxlorentz))
    call try(nf90_put_var(ncid, id_8, diag_rflux))
    call try(nf90_put_var(ncid, id_9, diag_rheight))

    call try(nf90_close(ncid))


END SUBROUTINE diagnostics

subroutine try(status)
! Catch error in reading netcdf fild.
INTEGER, INTENT(IN):: status

if (status /= NF90_noerr) THEN
    PRINT*,TRIM(ADJUSTL(NF90_STRERROR(status)))
end if

end subroutine try

SUBROUTINE save_snap(snap_number)

    !Making modifications so this works with netcdf instead, which would be far more sensible
    IMPLICIT NONE
    INTEGER:: snap_number
    CHARACTER(LEN = 64):: output_filename
    INTEGER:: ncid, vid
    INTEGER:: xs_id, xc_id, ys_id, yc_id
    INTEGER:: ax_id, ay_id, az_id, vy_id, my_id


    if (run_number < 10) then
      write (output_directory, "(A23, A2, I1, A1)") output_directory_root, '00', int(run_number), '/'
    else if (run_number < 100) then
      write (output_directory, "(A23, A1, I2, A1)") output_directory_root, '0', int(run_number), '/'
    else
      write (output_directory, "(A23, I3, A1)") output_directory_root, int(run_number), '/'
    end if

    if (snap_number < 10) then
        write (output_filename, "(A27,A3,I1,A3)") trim(output_directory), "000", snap_number, ".nc"
    else if (snap_number < 100) then
        write (output_filename, "(A27,A2,I2,A3)") trim(output_directory), "00", snap_number, ".nc"
    else if (snap_number < 1000) then
        write (output_filename, "(A27,A1,I3,A3)") trim(output_directory), "0", snap_number, ".nc"
    else if (snap_number < 10000) then
        write (output_filename, "(A27,I4,A3)") trim(output_directory), snap_number, ".nc"
    end if

    print*, output_filename

    call try(nf90_create(trim(output_filename), nf90_clobber, ncid))

    call try(nf90_def_dim(ncid, 'xs', nx+1, xs_id))
    call try(nf90_def_dim(ncid, 'ys', ny+1, ys_id))

    call try(nf90_def_dim(ncid, 'xc', nx, xc_id))
    call try(nf90_def_dim(ncid, 'yc', ny, yc_id))

    call try(nf90_def_var(ncid, 'ax', nf90_double, (/xc_id ,ys_id/), ax_id))
    call try(nf90_def_var(ncid, 'ay', nf90_double, (/xs_id ,yc_id/), ay_id))
    call try(nf90_def_var(ncid, 'az', nf90_double, (/xs_id ,ys_id/), az_id))

    call try(nf90_def_var(ncid, 'vy', nf90_double, (/xs_id ,ys_id/), vy_id))
    call try(nf90_def_var(ncid, 'my', nf90_double, (/xs_id ,ys_id/), my_id))

    call try(nf90_enddef(ncid))
    call try(nf90_close(ncid))

    call try(nf90_open(trim(output_filename), nf90_write, ncid))

    call try(nf90_inq_varid(ncid, 'ax', vid))
    call try(nf90_put_var(ncid, vid, ax(1:nx,0:ny)))

    call try(nf90_inq_varid(ncid, 'ay', vid))
    call try(nf90_put_var(ncid, vid, ay(0:nx,1:ny)))

    call try(nf90_inq_varid(ncid, 'az', vid))
    call try(nf90_put_var(ncid, vid, az(0:nx,0:ny)))

    call try(nf90_inq_varid(ncid, 'vy', vid))
    call try(nf90_put_var(ncid, vid, vouts(0:nx, 0:ny)*vout_masks(0:nx, 0:ny)))

    call try(nf90_inq_varid(ncid, 'my', vid))
    call try(nf90_put_var(ncid, vid, vy(0:nx,0:ny) + vouts(0:nx, 0:ny)*vout_masks(0:nx, 0:ny)))

    call try(nf90_close(ncid))

    print*, 'Saved snapshot number', snap_number, ' at time', t, 'to directory ', output_directory

END SUBROUTINE save_snap

!*******************************************************************************
END MODULE output
!********************************************************************************
