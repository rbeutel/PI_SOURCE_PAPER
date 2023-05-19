!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! - Ariane - (May - 2007)
!! 
!! bruno.blanke@univ-brest.fr and nicolas.grima@univ-brest.fr
!! 
!! This software is a computer program whose purpose is 
!! the computation of 3D streamlines in a given velocity field 
!! (as the output of an Ocean General Circulation Model) and 
!! subsequent water masses analyses.
!! 
!! This software is governed by the CeCILL license under French law and
!! abiding by the rules of distribution of free software.  You can  use, 
!! modify and/ or redistribute the software under the terms of the CeCILL
!! license as circulated by CEA, CNRS and INRIA at the following URL
!! "http://www.cecill.info". 
!! 
!! As a counterpart to the access to the source code and  rights to copy,
!! modify and redistribute granted by the license, users are provided only
!! with a limited warranty  and the software's author,  the holder of the
!! economic rights,  and the successive licensors  have only  limited
!! liability. 
!! 
!! In this respect, the user's attention is drawn to the risks associated
!! with loading,  using,  modifying and/or developing or reproducing the
!! software by the user in light of its specific status of free software,
!! that may mean  that it is complicated to manipulate,  and  that  also
!! therefore means  that it is reserved for developers  and  experienced
!! professionals having in-depth computer knowledge. Users are therefore
!! encouraged to load and test the software's suitability as regards their
!! requirements in conditions enabling the security of their systems and/or 
!! data to be ensured and,  more generally, to use and operate it in the 
!! same conditions as regards security. 
!! 
!! The fact that you are presently reading this means that you have had
!! knowledge of the CeCILL license and that you accept its terms.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!****h* ariane/mod_input_grid
!! NAME
!!   mod_input_grid (mod_input_grid.f90 - Fortran90 module)
!!
!! USAGE
!!   Include 'USE mod_input_grid' in the header of your Fortran 90 source 
!!   code.
!!   Then you'll have access to the subroutine:
!!      - sub_input_grid
!!      - sub_coord_dealloc
!!      - sub_scalef_dealloc
!!      - sub_tmask_dealloc
!!
!! FUNCTION
!!   * Read meshmask file
!! 
!! AUTHOR
!!   * Origin  : Bruno Blanke  (1992)
!!   * F77toF90: Nicolas Grima (April-May 2005)
!! 
!! CREATION DATE
!!   * April-May 2005
!!
!! HISTORY
!!   Date (dd/mm/yyyy/) - Modification(s)
!!
!! RESULT
!!   
!!
!! EXAMPLES
!!   
!!
!! NOTES
!!   ROBODoc header style.
!!
!! TODO
!!   
!!
!! PORTABILITY
!!         Machine-OS    - Fortran90/95 compiler
!!   * i686-pc-linux-gnu -         ifort
!!
!! SEE ALSO
!!   
!!
!! USES
!!   * USE mod_precision
!!   * USE mod_namelist
!!   * USE mod_reducmem
!!   * USE mod_netcdf
!!
!! USED BY
!!   * posini
!!   * trajec
!!
!! SOURCE
!!=========================================================================
MODULE mod_input_grid

  !------------------!
  ! USE ASSOCIAITION !
  !------------------!
  USE mod_precision
  USE mod_memory
  USE mod_namelist
  USE mod_reducmem
  USE mod_netcdf

  !-------------!
  ! DECLARATION !
  !-------------!
  IMPLICIT NONE

  ! 11 (12) variables
  ! We assume that coordinates and scale factors are in the same file.
  ! If not, user have to use "nco" tools to respect this assumption.

  REAL(kind=rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: &
       xx_tt, & ! x coordintes on T grid (tt, ss)
       xx_uu    ! x coordintes on U grid (uu)

  REAL(kind=rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: &
       yy_tt, & ! y coordintes on T grid (tt, ss)
       yy_vv    ! y coordintes on V grid (vv)

  REAL(kind=rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: &
       zz_ww    ! z coordinates on W grid

  REAL(kind=rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: &
       e1f,   & ! scale factor on F Cgrid (for B2C interpolation)
       e2f,   & ! scale factor on F Cgrid (for B2C interpolation)
       e3f,   & ! scale factor on F Cgrid (for B2C interpolation)
       e1t,   & ! scale factor on T Cgrid
       e2t,   & ! scale factor on T Cgrid
       e3t,   & ! scale factor on T Cgrid (also U and V grids)
       e2u,   & ! scale factor on U Cgrid
       e1v      ! scale factor on V Cgrid

  REAL(kind=rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: &
       tmask                              , &   ! mask on T grid
       wmask

  !- This is a ROMS variable to compute U, V and W in the input data module
  REAL(kind = rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: &
       h_roms

  !- This is a MARS variable to compute U, V and W in the input data module
  REAL(kind = rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: &
       h_mars

  !- This is Symphonie bathy used in mkseg0
  REAL(kind = rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: &
       h_symp

  !- This is a NEMO VVL variable to calculate e3t in input data module
  REAL(kind = rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: &
       e3t0
  !- This is a NEMO VVL variable to calculate e3t in input data module
  REAL(kind = rprec), DIMENSION(:,:,:, :), ALLOCATABLE :: &
       totaldepth
  !- This is a NEMO vvl variable marking the last water ceel in depth
  REAL(kind = rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: &
       mbathy
  !! NG 9 nov 2018
  REAL(kind = rprec), DIMENSION(:), ALLOCATABLE :: &
       sc_w , & ! S-coordinate at W-points
       cs_w     ! S-coordinate stretching curves at W-points
  REAL(kind = rprec) :: &
        hc      ! S-coordinate parameter, critical depth (m)


  !!***

CONTAINS
  !=========================================================================
  !!****f* mod_input_grid/sub_input_grid()
  !! NAME
  !!   sub_input_grid()
  !!
  !! FUNCTION
  !!   * Select how to read the grid (ROMS or OPA) and read it.
  !!
  !! AUTHOR
  !!   * Origin  : Nicolas Grima (April-May 2005)
  !! 
  !! CREATION DATE
  !!   * April-May 2005
  !!
  !! HISTORY
  !!   Date (dd/mm/yyyy/) - Modification(s)
  !!
  !! ARGUMENTS
  !!   * No arguments
  !!
  !! TODO
  !!   
  !! USES
  !!   * sub_open_netcdf_file  (mod_netcdf.f90)
  !!   * sub_select_var_dims   (mod_netcdf.f90)
  !!   * sub_read_netcdf_var4d (mod_netcdf.f90)
  !!   * sub_close_netcdf_file (mod_netcdf.f90)
  !!
  !! USED BY
  !!   * trajec
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_input_grid()

    !- Comments -!
    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)'==================='
    WRITE(lun_standard,*)'= INPUT GRID DATA ='
    WRITE(lun_standard,*)'==================='

    IF (key_roms) THEN

      CALL sub_input_grid_roms()

    ELSEIF ( key_mars) THEN

      CALL sub_input_grid_mars()

    ELSEIF ( key_symphonie) THEN

      CALL sub_input_grid_symphonie()

    ELSEIF ( key_B2C_grid) THEN

       IF (B2C_grid_Z_or_Sigma=='Z') THEN

          CALL sub_input_B2C_gridz()

       ELSE

          STOP
          !!TODO:  CALL sub_input_B2CD_grids()

       ENDIF

    ELSE

      CALL sub_input_grid_opa()

    ENDIF

    !--------------------------!
    !- Test tmask z dimension -!
    !--------------------------!
!!NG    IF (key_roms) THEN
!!NG      IF ( SIZE(tmask(:,:,:,:),dim=3) /= kmt-1) THEN
!!NG        WRITE(lun_error,*) 'mod_input_grid: mask have to have &
!!NG             & the same dimensions than kmt-1 !'
!!NG        WRITE(lun_error,*) ' Please check your netcdf grid file...'
!!NG        STOP
!!NG      ENDIF
!!NG    ELSE !- OPA-NEMO -!
!!NG      IF ( SIZE(tmask(:,:,:,:),dim=3) /= kmt) THEN
!!NG        WRITE(lun_error,*) 'mod_input_grid: mask have to have &
!!NG             & the same dimensions than kmt !'
!!NG        WRITE(lun_error,*) ' Please check your netcdf grid file...'
!!NG        STOP
!!NG      ENDIF
!!NG    ENDIF

  END SUBROUTINE sub_input_grid
  !!***

  !=========================================================================
  !!****f* mod_input_grid/sub_input_grid_opa()
  !! NAME
  !!   sub_input_grid()
  !!
  !! FUNCTION
  !!   * Read OPA input grid coordinates, scale factors and mask from one
  !!     netcdf file.
  !!
  !! AUTHOR
  !!   * Origin  : Nicolas Grima (April-May 2005)
  !! 
  !! CREATION DATE
  !!   * April-May 2005
  !!
  !! HISTORY
  !!   Date (dd/mm/yyyy/) - Modification(s)
  !!
  !! ARGUMENTS
  !!   * No arguments
  !!
  !! TODO
  !!   
  !! USES
  !!   * sub_open_netcdf_file  (mod_netcdf.f90)
  !!   * sub_select_var_dims   (mod_netcdf.f90)
  !!   * sub_read_netcdf_var4d (mod_netcdf.f90)
  !!   * sub_close_netcdf_file (mod_netcdf.f90)
  !!
  !! USED BY
  !!   * trajec
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_input_grid_opa()

    !-------------!
    ! Declaration !
    !-------------!
    INTEGER(kind = iprec) :: ncid  ! netcdf file ID
    INTEGER(kind = iprec) :: varid ! netcdf variable ID
    INTEGER(kind = iprec) :: dimx  ! dimension in x (i)
    INTEGER(kind = iprec) :: dimy  ! dimension in y (j)
    INTEGER(kind = iprec) :: dimz  ! dimension in z (k)
    INTEGER(kind = iprec) :: dimt  ! dimension in t (l)

    ! dimsorder: array where netcdf variable dimensions order are stored.
    !            This, because sometime dimensions order could be different 
    !            than (x, y, z, t). It could be (x, y, t, nothing) or 
    !            (z, t, nothing, nothing).
    !            More information: see sub_select_var_dims and 
    !            sub_read_netcdf_var4d in mod_netcdf.f90 file.
    INTEGER(kind = iprec), DIMENSION(4) :: dimsorder

    REAL(kind = rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: e3t2D
    REAL(kind = rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: e3tz

    INTEGER(kind=iprec) :: ii, jj, kk


    !-------------!
    ! Code begins !
    !-------------!

    !---------------------------------------!
    !-- Open NetCdf file with mesh values --!
    !---------------------------------------!
    CALL sub_open_netcdf_file(dir_mesh, fn_mesh, ncid)

    !======================================================================!
    !=========================== COORDINATES ==============================!
    !======================================================================!
    !-----------!
    !-- XX_TT --! (glamt)
    !-----------!
    CALL sub_select_var_dims(ncid, nc_var_xx_tt, varid, dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(xx_tt(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(xx_tt),'r','xx_tt','sub_input_grid_opa')
    CALL sub_read_netcdf_var4d(ncid, varid, xx_tt(:,:,:,:), dimsorder, dims_reg)

    WRITE(lun_standard,*)' - ', TRIM(nc_var_xx_tt),': max ', MAXVAL(xx_tt), ' min ', MINVAL(xx_tt)

    !-----------!
    !-- XX_UU --! (glamu)
    !-----------!
    CALL sub_select_var_dims(ncid, nc_var_xx_uu, varid, dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(xx_uu(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(xx_uu),'r','xx_uu','sub_input_grid_opa')
    CALL sub_read_netcdf_var4d(ncid, varid, xx_uu(:,:,:,:), dimsorder, dims_reg)

    WRITE(lun_standard,*)' - ', TRIM(nc_var_xx_uu),': max ', MAXVAL(xx_uu), ' min ', MINVAL(xx_uu)

    !-----------!
    !-- YY_TT --! (gphit)
    !-----------!
    CALL sub_select_var_dims(ncid, nc_var_yy_tt, varid, dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(yy_tt(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(yy_tt),'r','yy_tt','sub_input_grid_opa')

    CALL sub_read_netcdf_var4d(ncid, varid, yy_tt(:,:,:,:), dimsorder, dims_reg)
    WRITE(lun_standard,*)' - ', TRIM(nc_var_yy_tt),': max ', MAXVAL(yy_tt), ' min ', MINVAL(yy_tt)

    !-----------!
    !-- YY_VV --! (gphiv)
    !-----------!
    CALL sub_select_var_dims(ncid, nc_var_yy_vv, varid, dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(yy_vv(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(yy_vv),'r','yy_vv','sub_input_grid_opa')

    CALL sub_read_netcdf_var4d(ncid, varid, yy_vv(:,:,:,:), dimsorder, dims_reg)
    WRITE(lun_standard,*)' - ', TRIM(nc_var_yy_vv),': max ', MAXVAL(yy_vv), ' min ', MINVAL(yy_vv)

    !-----------!
    !-- ZZ_WW --!
    !-----------!
    CALL sub_select_var_dims(ncid, nc_var_zz_ww, varid, dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(zz_ww(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(zz_ww),'r','zz_ww','sub_input_grid_opa')

    CALL sub_read_netcdf_var4d(ncid, varid, zz_ww(:,:,:,:), dimsorder, dims_reg)
    WRITE(lun_standard,*)' - ', TRIM(nc_var_zz_ww),': max ', MAXVAL(zz_ww), ' min ', MINVAL(zz_ww)

    !======================================================================!
    !========================= SCALE FACTORS ==============================!
    !======================================================================!
    !---------!
    !-- E2U --!
    !---------!
    CALL sub_select_var_dims(ncid, nc_var_e2u, varid, dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(e2u(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(e2u),'r','e2u','sub_input_grid_opa')

    CALL sub_read_netcdf_var4d(ncid, varid, e2u(:,:,:,:), dimsorder, dims_reg)
    WRITE(lun_standard,*)' - ', TRIM(nc_var_e2u),': max ', MAXVAL(e2u), ' min ', MINVAL(e2u)

    !---------!
    !-- E1V --!
    !---------!
    CALL sub_select_var_dims(ncid, nc_var_e1v, varid, dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(e1v(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(e1v),'r','e1v','sub_input_grid_opa')

    CALL sub_read_netcdf_var4d(ncid, varid, e1v(:,:,:,:), dimsorder, dims_reg)
    WRITE(lun_standard,*)' - ', TRIM(nc_var_e1v),': max ', MAXVAL(e1v), ' min ', MINVAL(e1v)

    !---------!
    !-- E1T --!
    !---------!
    CALL sub_select_var_dims(ncid, nc_var_e1t, varid, dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(e1t(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(e1t),'r','e1t','sub_input_grid_opa')

    CALL sub_read_netcdf_var4d(ncid, varid, e1t(:,:,:,:), dimsorder, dims_reg)
    WRITE(lun_standard,*)' - ', TRIM(nc_var_e1t),': max ', MAXVAL(e1t), ' min ', MINVAL(e1t)

    !---------!
    !-- E2T --!
    !---------!
    CALL sub_select_var_dims(ncid, nc_var_e2t, varid, dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(e2t(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(e2t),'r','e2t','sub_input_grid_opa')

    CALL sub_read_netcdf_var4d(ncid, varid, e2t(:,:,:,:), dimsorder, dims_reg)
    WRITE(lun_standard,*)' - ', TRIM(nc_var_e2t),': max ', MAXVAL(e2t), ' min ', MINVAL(e2t)

    !---------!
    !-- E3T --!
    !---------!

    IF (TRIM(mesh_type) == c_def) THEN

       CALL sub_select_var_dims(ncid, nc_var_e3t, varid, dimsorder, dimx, dimy, dimz, dimt)
       ALLOCATE(e3t(1:dimx, 1:dimy, 1:dimz, 1:dimt))
       CALL sub_memory(size(e3t),'r','e3t','sub_input_grid_opa')

       CALL sub_read_netcdf_var4d(ncid, varid, e3t(:,:,:,:), dimsorder, dims_reg)

       write (lun_standard, *) key_vvl, 'key_vvl'
       if (key_vvl) then
          ALLOCATE(e3t0(1:dimx, 1:dimy, 1:dimz, 1:dimt))
          ALLOCATE(totaldepth(1:dimx, 1:dimy, 1, 1:dimt))

          CALL sub_select_var_dims(ncid, nc_var_mbathy, varid, dimsorder, dimx, dimy, dimz, dimt)
          ALLOCATE(mbathy(1:dimx, 1:dimy, 1:dimz, 1:dimt))
          CALL sub_memory(size(mbathy),'r','mbathy','sub_input_grid_opa')

          CALL sub_read_netcdf_var4d(ncid, varid, mbathy(:,:,:,:), dimsorder, dims_reg)

          write (lun_standard, *) "Calculating totaldepth"
          totaldepth(:, :, 1, 1) = 0.d0
          DO jj = 1, dimy
             DO ii = 1, dimx
                      totaldepth(ii, jj, 1, 1) = ZZ_WW(ii, jj, NINT(mbathy(ii, jj, 1, 1)), 1)
             ENDDO
          ENDDO
          write (lun_standard, *) "Total depth max", maxval(totaldepth)
          e3t0 = e3t
       endif

    ELSEIF(TRIM(mesh_type) == 'nemov3') THEN

       CALL sub_select_var_dims(ncid, nc_var_e3t2D, varid, dimsorder, dimx, dimy, dimz, dimt)
       ALLOCATE(e3t2D(1:dimx, 1:dimy, 1:dimz, 1:dimt))
       CALL sub_memory(size(e3t2D),'r','e3t2D','sub_input_grid_opa')

       CALL sub_read_netcdf_var4d(ncid, varid, e3t2D(:,:,:,:), dimsorder, dims_reg)

       CALL sub_select_var_dims(ncid, nc_var_e3tz, varid, dimsorder, dimx, dimy, dimz, dimt)
       ALLOCATE(e3tz(1:dimx, 1:dimy, 1:dimz, 1:dimt))
       CALL sub_memory(size(e3tz),'r','e3tz','sub_input_grid_opa')

       CALL sub_read_netcdf_var4d(ncid, varid, e3tz(:,:,:,:), dimsorder, dims_reg)

       CALL sub_select_var_dims(ncid, nc_var_mbathy, varid, dimsorder, dimx, dimy, dimz, dimt)
       ALLOCATE(mbathy(1:dimx, 1:dimy, 1:dimz, 1:dimt))
       CALL sub_memory(size(mbathy),'r','mbathy','sub_input_grid_opa')

       CALL sub_read_netcdf_var4d(ncid, varid, mbathy(:,:,:,:), dimsorder, dims_reg)

       ALLOCATE(e3t(size(e3t2D,1), size(e3t2D,2), size(e3tz,3),size(e3tz,4)))

       DO jj = 1, size(e3t,2)
           DO ii = 1, size(e3t,1)
              e3t(ii,jj,:,1) = e3tz(1,1,:,1)
           ENDDO
        ENDDO

       DO jj = 1, size(e3t,2)
           DO ii = 1, size(e3t,1)
              IF (NINT(mbathy(ii,jj,1,1)) > 0_iprec) THEN
                   e3t(ii,jj,NINT(mbathy(ii,jj,1,1)),1) = e3t2D(ii,jj,1,1)
              ENDIF
           ENDDO
        ENDDO
  
       DEALLOCATE(e3t2D)
       CALL sub_memory(-size(e3t2D),'r','e3t2D','sub_input_grid_opa')

       DEALLOCATE(e3tz)
       CALL sub_memory(-size(e3tz),'r','e3tz','sub_input_grid_opa')

       DEALLOCATE(mbathy)
       CALL sub_memory(-size(mbathy),'r','mbathy','sub_input_grid_opa')
       
    ENDIF

    WRITE(lun_standard,*)' - ', TRIM(nc_var_e3t),': max ', MAXVAL(e3t), ' min ', MINVAL(e3t)

    !=============================================================!
    !========================= MASK ==============================!
    !=============================================================!
    !-----------!
    !-- TMASK --!
    !-----------!
    CALL sub_select_var_dims(ncid, nc_var_tmask, varid, dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(tmask(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(tmask),'r','tmask','sub_input_grid_opa')

    CALL sub_read_netcdf_var4d(ncid, varid, tmask(:,:,:,:), dimsorder, dims_reg)

    WHERE(tmask(:,:,:,:) == nc_mask_val) tmask(:,:,:,:) = 0._rprec

    WRITE(lun_standard,*)' - ', TRIM(nc_var_tmask),': max ', MAXVAL(tmask), ' min ', MINVAL(tmask)

    IF ( MAXVAL(ABS(tmask(:,:,:,:))) > 1) THEN
!!      WRITE(lun_error,*) 'mod_input_grid: sub_input_grid: problem with tmask values!'
!!      STOP
!! For MITgcm (David Wang)
WRITE(lun_standard,*) 'mod_input_grid: sub_input_grid: we change to 1 the max tmask!'
       WHERE(ABS(tmask(:,:,:,:)) > 1) tmask(:,:,:,:)=1
       WRITE(lun_standard,*)' - ', TRIM(nc_var_tmask),': max ', MAXVAL(tmask), ' min ', MINVAL(tmask)
    ENDIF

    !-------------------------------!
    !----- CLOSE NetCDF file -------!
    !-------------------------------!
    CALL sub_close_netcdf_file(ncid)

    PRINT *, ' '

  END SUBROUTINE sub_input_grid_opa

  !=========================================================================
  !!****f* mod_input_grid/sub_input_tmask_surf_opa()
  !! NAME
  !!   sub_input_tmask_surf_opa()
  !!
  !! FUNCTION
  !!   * Read OPA tmask from a netcdf file.
  !!
  !! AUTHOR
  !!   * Origin  : Nicolas Grima (September 2015)
  !! 
  !! CREATION DATE
  !!   * September 2015
  !!
  !! HISTORY
  !!   Date (dd/mm/yyyy/) - Modification(s)
  !!
  !! ARGUMENTS
  !!   * No arguments
  !!
  !! TODO
  !!   
  !! USES
  !!   * sub_open_netcdf_file  (mod_netcdf.f90)
  !!   * sub_select_var_dims   (mod_netcdf.f90)
  !!   * sub_read_netcdf_var4d (mod_netcdf.f90)
  !!   * sub_close_netcdf_file (mod_netcdf.f90)
  !!
  !! USED BY
  !!   * mkseg0 and mkseg
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_input_tmask_surf_opa()

    !-------------!
    ! Declaration !
    !-------------!
    INTEGER(kind = iprec) :: ncid  ! netcdf file ID
    INTEGER(kind = iprec) :: varid ! netcdf variable ID
    INTEGER(kind = iprec) :: dimx  ! dimension in x (i)
    INTEGER(kind = iprec) :: dimy  ! dimension in y (j)
    INTEGER(kind = iprec) :: dimz  ! dimension in z (k)
    INTEGER(kind = iprec) :: dimt  ! dimension in t (l)

    ! dimsorder: array where netcdf variable dimensions order are stored.
    !            This, because sometime dimensions order could be different 
    !            than (x, y, z, t). It could be (x, y, t, nothing) or 
    !            (z, t, nothing, nothing).
    !            More information: see sub_select_var_dims and 
    !            sub_read_netcdf_var4d in mod_netcdf.f90 file.
    INTEGER(kind = iprec), DIMENSION(4) :: dimsorder

    REAL(kind = rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: e3t2D
    REAL(kind = rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: mbathy
    REAL(kind = rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: e3tz

    INTEGER(kind=iprec) :: ii, jj
    
    !-------------!
    ! Code begins !
    !-------------!

    !---------------------------------------!
    !-- Open NetCdf file with mesh values --!
    !---------------------------------------!
    CALL sub_open_netcdf_file(dir_mesh, fn_mesh, ncid)

    !=============================================================!
    !========================= MASK ==============================!
    !=============================================================!
    !-----------!
    !-- TMASK --!
    !-----------!
    CALL sub_select_var_dims(ncid, nc_var_tmask, varid, dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(tmask(1:dimx, 1:dimy, 1, 1))
    CALL sub_memory(size(tmask),'r','tmask','sub_input_grid_opa')

    CALL sub_read_netcdf_4D_surf(ncid, varid, tmask(:,:,:,:), dimx, dimy)

    WHERE(tmask(:,:,:,:) == nc_mask_val) tmask(:,:,:,:) = 0._rprec

    WRITE(lun_standard,*)' - ', TRIM(nc_var_tmask),': max ', MAXVAL(tmask), ' min ', MINVAL(tmask)

    IF ( MAXVAL(ABS(tmask(:,:,:,:))) > 1) THEN
!!      WRITE(lun_error,*) 'mod_input_grid: sub_input_grid: problem with tmask values!'
!!      STOP
!! For MITgcm (David Wang)
       WRITE(lun_standard,*) 'mod_input_grid: sub_input_grid: we change to 1 the max tmask!'
       WHERE(ABS(tmask(:,:,:,:)) > 1) tmask(:,:,:,:)=1
       WRITE(lun_standard,*)' - ', TRIM(nc_var_tmask),': max ', MAXVAL(tmask), ' min ', MINVAL(tmask)
    ENDIF

    !-------------------------------!
    !----- CLOSE NetCDF file -------!
    !-------------------------------!
    CALL sub_close_netcdf_file(ncid)

    PRINT *, ' '

  END SUBROUTINE sub_input_tmask_surf_opa

  !!***
  !=========================================================================
  !!****f* mod_input_grid/sub_input_grid_roms()
  !! NAME
  !!   sub_input_grid()
  !!
  !! FUNCTION
  !!   * Read ROMS input grid coordinates and mask from one
  !!     netcdf file and compute scale factors.
  !!
  !! AUTHOR
  !!   * Origin  : Nicolas Grima (November 2005)
  !! 
  !! CREATION DATE
  !!   * April-May 2005
  !!
  !! HISTORY
  !!   Date (dd/mm/yyyy/) - Modification(s)
  !!
  !! ARGUMENTS
  !!   * No arguments
  !!
  !! TODO
  !!   
  !! USES
  !!   * sub_open_netcdf_file  (mod_netcdf.f90)
  !!   * sub_select_var_dims   (mod_netcdf.f90)
  !!   * sub_read_netcdf_var4d (mod_netcdf.f90)
  !!   * sub_close_netcdf_file (mod_netcdf.f90)
  !!
  !! USED BY
  !!   * trajec
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_input_grid_roms()

    !-------------!
    ! Declaration !
    !-------------!
    INTEGER(kind = iprec) :: k
    INTEGER(kind = iprec) :: ncid  ! netcdf file ID
    INTEGER(kind = iprec) :: varid ! netcdf variable ID
    INTEGER(kind = iprec) :: dimx  ! dimension in x (i)
    INTEGER(kind = iprec) :: dimy  ! dimension in y (j)
    INTEGER(kind = iprec) :: dimz  ! dimension in z (k)
    INTEGER(kind = iprec) :: dimt  ! dimension in t (l)

    ! dimsorder: array where netcdf variable dimensions order are stored.
    !            This, because sometime dimensions order could be different 
    !            than (x, y, z, t). It could be (x, y, t, nothing) or 
    !            (z, t, nothing, nothing).
    !            More information: see sub_select_var_dims and 
    !            sub_read_netcdf_var4d in mod_netcdf.f90 file.
    INTEGER(kind = iprec), DIMENSION(4) :: dimsorder

    REAL(kind = rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: pn_roms, pm_roms
    REAL(kind = rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: mask_roms

    !-------------!
    ! Code begins !
    !-------------!

    !---------------------------------------!
    !-- Open NetCdf file with mesh values --!
    !---------------------------------------!
    CALL sub_open_netcdf_file(dir_grd_roms, fn_grd_roms, ncid)

    !====================== LON_U ================================!
    !-----------!
    !-- XX_UU --!
    !-----------!
    CALL sub_select_var_dims(ncid, nc_var_lon_u_roms, varid, &
         dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(xx_uu(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(xx_uu),'r','xx_uu','sub_input_grid_roms')

    CALL sub_read_netcdf_var4d(ncid, varid, &
         xx_uu(1:dims_reg(1,3)-1,:,:,:), dimsorder, dims_reg)
    WRITE(lun_standard,*)' - ', TRIM(nc_var_lon_u_roms), &
         ': max ', MAXVAL(xx_uu(1:dims_reg(1,3)-1,:,:,:)), &
         '  min ', MINVAL(xx_uu(1:dims_reg(1,3)-1,:,:,:))

    !====================== LAT_V ================================!
    !-----------!
    !-- YY_VV --!
    !-----------!
    CALL sub_select_var_dims(ncid, nc_var_lat_v_roms, varid, &
         dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(yy_vv(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(yy_vv),'r','yy_vv','sub_input_grid_roms')

    CALL sub_read_netcdf_var4d(ncid, varid, &
         yy_vv(:,1:dims_reg(2,3)-1,:,:), dimsorder, dims_reg)
    WRITE(lun_standard,*)' - ', TRIM(nc_var_lat_v_roms), &
         ': max ', MAXVAL(yy_vv(:,1:dims_reg(2,3)-1,:,:)), &
         '  min ', MINVAL(yy_vv(:,1:dims_reg(2,3)-1,:,:))

    !====================== LON_RHO ================================!
    !-----------!
    !-- XX_TT --!
    !-----------!
    CALL sub_select_var_dims(ncid, nc_var_lon_rho_roms, varid, &
         dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(xx_tt(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(xx_tt),'r','xx_tt','sub_input_grid_roms')

    CALL sub_read_netcdf_var4d(ncid, varid, xx_tt(:,:,:,:), dimsorder, dims_reg)
    WRITE(lun_standard,*)' - ', TRIM(nc_var_lon_rho_roms), &
         ': max ', MAXVAL(xx_tt), ' min ', MINVAL(xx_tt)

    !====================== LAT_RHO ================================!
    !-----------!
    !-- YY_TT --!
    !-----------!
    CALL sub_select_var_dims(ncid, nc_var_lat_rho_roms, varid, &
         dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(yy_tt(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(yy_tt),'r','yy_tt','sub_input_grid_roms')

    CALL sub_read_netcdf_var4d(ncid, varid, yy_tt(:,:,:,:), dimsorder, dims_reg)
    WRITE(lun_standard,*)' - ', TRIM(nc_var_lat_rho_roms), &
         ': max ', MAXVAL(yy_tt), ' min ', MINVAL(yy_tt)

    !=============================================================!
    !========================= MASK ==============================!
    !=============================================================!
    !-----------!
    !-- TMASK --!
    !-----------!
    !==================== READMASK_RHO (2D)=======================!
    CALL sub_select_var_dims(ncid, nc_var_mask_rho_roms, varid, &
         dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(mask_roms(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(mask_roms),'r','mask_roms','sub_input_grid_roms')

    CALL sub_read_netcdf_var4d(ncid, varid, &
         mask_roms(:,:,:,:), dimsorder, dims_reg)
    WRITE(lun_standard,*)' - ', TRIM(nc_var_mask_rho_roms), &
         ': max ', MAXVAL(mask_roms), ' min ', MINVAL(mask_roms)
    !==================== Compute TMASK (3D)======================!
    ALLOCATE(tmask(1:dims_reg(1,3), 1:dims_reg(2,3), 1:dims_reg(3,3), 1))
    CALL sub_memory(size(tmask),'r','tmask','sub_input_grid_roms')

    DO k=1, dims_reg(3,3)
      tmask(:,:,k,1) = mask_roms(:,:,1,1)
    ENDDO
    tmask(:,:,dims_reg(3,3),:) = 0._rprec

    CALL sub_memory(-size(mask_roms),'r','mask_roms','sub_input_grid_roms')
    DEALLOCATE(mask_roms)

    !======================================================================!
    !========================= SCALE FACTORS ==============================!
    !======================================================================!
    !-------------------------------------------!
    !-- E2U -> e2u(i,j)=2./(pn(i,j)+pn(i+1,j)) -!
    !-------------------------------------------!
    !============================= READ PN ================================!
    CALL sub_select_var_dims(ncid, nc_var_pn_roms, varid, &
         dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(pn_roms(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(pn_roms),'r','pn_roms','sub_input_grid_roms')

    CALL sub_read_netcdf_var4d(ncid, varid, &
         pn_roms(:,:,:,:), dimsorder, dims_reg)
    WRITE(lun_standard,*)' - ', TRIM(nc_var_pn_roms), &
         ': max ', MAXVAL(pn_roms), ' min ', MINVAL(pn_roms)
    !=========================== COMPUTE E2T ==============================!
    ALLOCATE(e2t(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(e2t),'r','e2t','sub_input_grid_roms')
    e2t(:,:,:,:) = 1._rprec / pn_roms(:,:,:,:)
    WRITE(lun_standard,*)' - e2t (c): max ', MAXVAL(e2t), ' min ', MINVAL(e2t)
    !=========================== COMPUTE E2U ==============================!
    ALLOCATE(e2u(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(e2u),'r','e2u','sub_input_grid_roms')

    e2u(1:dims_reg(1,3)-1,:,:,:) = 2._rprec / ( pn_roms(1:dims_reg(1,3)-1,:,:,:) + &
         pn_roms(2:dims_reg(1,3),:,:,:) )

    CALL sub_memory(-size(pn_roms),'r','pn_roms','sub_input_grid_roms')
    DEALLOCATE(pn_roms)
    
    e2u(dims_reg(1,3),:,:,:) = e2u(dims_reg(1,3)-1,:,:,:)
    WRITE(lun_standard,*)' - e2u (c): max ', MAXVAL(e2u), ' min ', MINVAL(e2u)

    !-------------------------------------------!
    !-- E1V -> e1v(i,j)=2./(pm(i,j)+pm(i,j+1)) -!
    !-------------------------------------------!
    !============================= READ PM ================================!
    CALL sub_select_var_dims(ncid, nc_var_pm_roms, varid, &
         dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(pm_roms(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(pm_roms),'r','pm_roms','sub_input_grid_roms')

    CALL sub_read_netcdf_var4d(ncid, varid,  &
         pm_roms(:,:,:,:), dimsorder, dims_reg)
    WRITE(lun_standard,*)' - ', TRIM(nc_var_pm_roms), &
         ': max ', MAXVAL(pm_roms), ' min ', MINVAL(pm_roms)
    !=========================== COMPUTE E1T ==============================!
    ALLOCATE(e1t(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(e1t),'r','e1t','sub_input_grid_roms')

    e1t(:,:,:,:) = 1._rprec / pm_roms(:,:,:,:)
    WRITE(lun_standard,*)' - e1t (c): max ', MAXVAL(e1t), ' min ', MINVAL(e1t)
    !=========================== COMPUTE E1V ==============================!
    ALLOCATE(e1v(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(e1v),'r','e1v','sub_input_grid_roms')
    e1v(:,1:dims_reg(2,3)-1,:,:) = 2._rprec / &
         ( pm_roms(:,1:dims_reg(2,3)-1,:,:) + &
         pm_roms(:,2:dims_reg(2,3),:,:) )

    CALL sub_memory(-size(pm_roms),'r','pm_roms','sub_input_grid_roms')
    DEALLOCATE(pm_roms)
    e1v(:,dims_reg(2,3),:,:) =  e1v(:,dims_reg(2,3)-1,:,:)
    WRITE(lun_standard,*)' - e1v (c): max ', MAXVAL(e1v), ' min ', MINVAL(e1v)

    !------------------------------------------!
    !-- E3T is computed in input data module --!
    !------------------------------------------!

    !======================================================================!
    !============================== MORE ==================================!
    !======================================================================!

    !=============================== H ====================================!
    CALL sub_select_var_dims(ncid, nc_var_h_roms, varid, &
         dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(h_roms(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(h_roms),'r','h_roms','sub_input_grid_roms')

    CALL sub_read_netcdf_var4d(ncid, varid,  &
         h_roms(:,:,:,:), dimsorder, dims_reg)
    WRITE(lun_standard,*)' - ', TRIM(nc_var_h_roms), &
         ': max ', MAXVAL(h_roms), ' min ', MINVAL(h_roms)

  END SUBROUTINE sub_input_grid_roms
  !!***
  !!***
  !=========================================================================
  !!****f* mod_input_grid/sub_input_grid_mars()
  !! NAME
  !!   sub_input_grid()
  !!
  !! FUNCTION
  !!   * Read MARS3D input grid coordinates and mask from one
  !!     netcdf file and compute scale factors.
  !!
  !! AUTHOR
  !!   * Origin  : Nicolas Grima (June 2018)
  !! 
  !! CREATION DATE
  !!   * April-May 2005
  !!
  !! HISTORY
  !!   Date (dd/mm/yyyy/) - Modification(s)
  !!
  !! ARGUMENTS
  !!   * No arguments
  !!
  !! TODO
  !!   
  !! USES
  !!   * sub_open_netcdf_file  (mod_netcdf.f90)
  !!   * sub_select_var_dims   (mod_netcdf.f90)
  !!   * sub_read_netcdf_var4d (mod_netcdf.f90)
  !!   * sub_close_netcdf_file (mod_netcdf.f90)
  !!
  !! USED BY
  !!   * trajec
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_input_grid_mars()

    !-------------!
    ! Declaration !
    !-------------!
    INTEGER(kind = iprec) :: k
    INTEGER(kind = iprec) :: ncid  ! netcdf file ID
    INTEGER(kind = iprec) :: varid ! netcdf variable ID
    INTEGER(kind = iprec) :: dimx  ! dimension in x (i)
    INTEGER(kind = iprec) :: dimy  ! dimension in y (j)
    INTEGER(kind = iprec) :: dimz  ! dimension in z (k)
    INTEGER(kind = iprec) :: dimt  ! dimension in t (l)

    ! dimsorder: array where netcdf variable dimensions order are stored.
    !            This, because sometime dimensions order could be different 
    !            than (x, y, z, t). It could be (x, y, t, nothing) or 
    !            (z, t, nothing, nothing).
    !            More information: see sub_select_var_dims and 
    !            sub_read_netcdf_var4d in mod_netcdf.f90 file.
    INTEGER(kind = iprec), DIMENSION(4) :: dimsorder

    REAL(kind = rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: mask_mars

    !-------------!
    ! Code begins !
    !-------------!

    !---------------------------------------!
    !-- Open NetCdf file with mesh values --!
    !---------------------------------------!
    CALL sub_open_netcdf_file(dir_grd_mars, fn_grd_mars, ncid)

    !====================== LON_U ================================!
    !-----------!
    !-- XX_UU --!
    !-----------!
    CALL sub_select_var_dims(ncid, nc_var_lon_u_mars, varid, &
         dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(xx_uu(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(xx_uu),'r','xx_uu','sub_input_grid_mars')

    CALL sub_read_netcdf_var4d(ncid, varid, &
         xx_uu(:,:,:,:), dimsorder, dims_reg)
    WRITE(lun_standard,*)' - ', TRIM(nc_var_lon_u_mars), &
         ': max ', MAXVAL(xx_uu(1:dims_reg(1,3)-1,:,:,:)), &
         '  min ', MINVAL(xx_uu(1:dims_reg(1,3)-1,:,:,:))

    !====================== LAT_V ================================!
    !-----------!
    !-- YY_VV --!
    !-----------!
    CALL sub_select_var_dims(ncid, nc_var_lat_v_mars, varid, &
         dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(yy_vv(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(yy_vv),'r','yy_vv','sub_input_grid_mars')

    CALL sub_read_netcdf_var4d(ncid, varid, &
         yy_vv(:,:,:,:), dimsorder, dims_reg)
    WRITE(lun_standard,*)' - ', TRIM(nc_var_lat_v_mars), &
         ': max ', MAXVAL(yy_vv(:,1:dims_reg(2,3)-1,:,:)), &
         '  min ', MINVAL(yy_vv(:,1:dims_reg(2,3)-1,:,:))

    !====================== LON_T ================================!
    !-----------!
    !-- XX_TT --!
    !-----------!
    CALL sub_select_var_dims(ncid, nc_var_lon_t_mars, varid, &
         dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(xx_tt(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(xx_tt),'r','xx_tt','sub_input_grid_mars')

    CALL sub_read_netcdf_var4d(ncid, varid, xx_tt(:,:,:,:), dimsorder, dims_reg)
    WRITE(lun_standard,*)' - ', TRIM(nc_var_lon_t_mars), &
         ': max ', MAXVAL(xx_tt), ' min ', MINVAL(xx_tt)

    !====================== LAT_T ================================!
    !-----------!
    !-- YY_TT --!
    !-----------!
    CALL sub_select_var_dims(ncid, nc_var_lat_t_mars, varid, &
         dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(yy_tt(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(yy_tt),'r','yy_tt','sub_input_grid_mars')

    CALL sub_read_netcdf_var4d(ncid, varid, yy_tt(:,:,:,:), dimsorder, dims_reg)
    WRITE(lun_standard,*)' - ', TRIM(nc_var_lat_t_mars), &
         ': max ', MAXVAL(yy_tt), ' min ', MINVAL(yy_tt)

    !=============================================================!
    !========================= MASK ==============================!
    !=============================================================!
    !-----------!
    !-- TMASK --!
    !-----------!
    !==================== READMASK_T (2D)=======================!
    CALL sub_select_var_dims(ncid, nc_var_mask_t_mars, varid, &
         dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(mask_mars(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(mask_mars),'r','mask_mars','sub_input_grid_mars')

    CALL sub_read_netcdf_var4d(ncid, varid, &
         mask_mars(:,:,:,:), dimsorder, dims_reg)
    WRITE(lun_standard,*)' - ', TRIM(nc_var_mask_t_mars), &
         ': max ', MAXVAL(mask_mars), ' min ', MINVAL(mask_mars)
    !==================== Compute TMASK (3D)======================!
    ALLOCATE(tmask(1:dims_reg(1,3), 1:dims_reg(2,3), 1:dims_reg(3,3)+1, 1))
    CALL sub_memory(size(tmask),'r','tmask','sub_input_grid_mars')

    DO k=1, dims_reg(3,3)
       tmask(:,:,k,1) = mask_mars(:,:,1,1)                 ! 3D dimension
    ENDDO
    tmask(:,:,dims_reg(3,3)+1,:) = 0._rprec

    CALL sub_memory(-size(mask_mars),'r','mask_mars','sub_input_grid_mars')
    DEALLOCATE(mask_mars)

    !======================================================================!
    !========================= SCALE FACTORS ==============================!
    !======================================================================!

    !=========================== COMPUTE E2T ==============================!
    !- read netcdf variable dimensions !
    CALL sub_select_var_dims(ncid, nc_var_e2t_mars, varid, dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(e2t(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(e2t),'r','e2t','sub_input_grid_mars')
    ! read e2t netcdf variable !
    CALL sub_read_netcdf_var4d(ncid, varid, e2t(:,:,:,:), dimsorder, dims_reg)

    WHERE((e2t > 10000._8).or.(e2t < 0))
       e2t =0._8
    END WHERE

    WRITE(lun_standard,*)' - e2t: max ', MAXVAL(e2t), ' min ', MINVAL(e2t)

    !=========================== COMPUTE E2U ==============================!
    CALL sub_select_var_dims(ncid, nc_var_e2u_mars, varid, dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(e2u(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(e2u),'r','e2u','sub_input_grid_mars')
    ! read e2u netcdf variable !
    CALL sub_read_netcdf_var4d(ncid, varid, e2u(:,:,:,:), dimsorder, dims_reg)

    WHERE((e2u > 10000._8).or.(e2u < 0))
       e2u =0._8
    END WHERE

    WRITE(lun_standard,*)' - e2u: max ', MAXVAL(e2u), ' min ', MINVAL(e2u)

    !=========================== COMPUTE E1T ==============================!
    CALL sub_select_var_dims(ncid, nc_var_e1t_mars, varid, dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(e1t(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(e1t),'r','e1t','sub_input_grid_mars')
    ! read e1t netcdf variable !
    CALL sub_read_netcdf_var4d(ncid, varid, e1t(:,:,:,:), dimsorder, dims_reg)

    WHERE((e1t > 10000._8).or.(e1t < 0))
       e1t =0._8
    END WHERE

    WRITE(lun_standard,*)' - e1t: max ', MAXVAL(e1t), ' min ', MINVAL(e1t)

    !=========================== COMPUTE E1V ==============================!
    CALL sub_select_var_dims(ncid, nc_var_e1v_mars, varid, dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(e1v(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(e1v),'r','e1v','sub_input_grid_mars')
    ! read e1v netcdf variable !
    CALL sub_read_netcdf_var4d(ncid, varid, e1v(:,:,:,:), dimsorder, dims_reg)

    WHERE((e1v > 10000._8).or.(e1v < 0))
       e1v =0._8
    END WHERE

    WRITE(lun_standard,*)' - e1v: max ', MAXVAL(e1v), ' min ', MINVAL(e1v)

    !------------------------------------------!
    !-- E3T is computed in input data module --!
    !------------------------------------------!

    !======================================================================!
    !============================== MORE ==================================!
    !======================================================================!

    !======== Read MARS hc, sc_w and Cs_w ======!
    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' ------ Reading hc, sc_w and Cs_w ------'

    !- Read hc -!
    CALL sub_read_netcdf_varid_ndims(ncid, nc_var_hc_mars, varid)
    CALL sub_read_netcdf_var( ncid, varid, hc)

    !- Read sc_w -!
    ALLOCATE(sc_w(dims_reg(3,3)+1))
    CALL sub_read_netcdf_varid_ndims(ncid, nc_var_sc_w_mars, varid)
    CALL sub_read_netcdf_var1d( ncid, varid, sc_w(dims_reg(3,3)+1:1:-1))

    !!sc_w(:) = [( -(ii-1) * 0.01_8, ii=1,dims_reg(3,3)+1 )]

    !- Read Cs_w -!
    ALLOCATE(cs_w(dims_reg(3,3)+1))
    CALL sub_read_netcdf_varid_ndims(ncid, nc_var_Cs_w_mars, varid)
    CALL sub_read_netcdf_var1d( ncid, varid, Cs_w(dims_reg(3,3)+1:1:-1))

    !!Cs_w(:)=(1.0_8-0.0_8)*SINH(6.0_8*sc_w(:))/SINH(6.0_8) &
    !!     + 0.5_8*0.0_8*(TANH(6.0_8*(sc_w(:)+0.5_8))- &
    !!     TANH(0.5_8*6.0_8))/TANH(6.0_8/2.0_8)

    WRITE(lun_standard,*)' - hc:   ',hc
    WRITE(lun_standard,*)' - SC_W: ',sc_w
    WRITE(lun_standard,*)' - CS_w: ',cs_w

    !============================= Bathy ==================================!
    CALL sub_select_var_dims(ncid, nc_var_bathy_t_mars, varid, &
         dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(h_mars(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(h_mars),'r','h_mars','sub_input_grid_mars')

    CALL sub_read_netcdf_var4d(ncid, varid,  &
         h_mars(:,:,:,:), dimsorder, dims_reg)

    !TODO NG: read HC and replace -999 in bathy by hc
    where ((h_mars(:,:,:,:) < -900._rprec).or.(h_mars(:,:,:,:) > 15000._rprec))
       h_mars(:,:,:,:) = hc
    end where

    WRITE(lun_standard,*)' - ', TRIM(nc_var_bathy_t_mars), &
         ': max ', MAXVAL(h_mars), ' min ', MINVAL(h_mars)

  END SUBROUTINE sub_input_grid_mars
  !!***

  !=========================================================================
  !!****f* mod_input_grid/sub_input_grid_symphonie()
  !! NAME
  !!   sub_input_grid()
  !!
  !! FUNCTION
  !!   * Read SYMPHONIE input grid coordinates from one
  !!     netcdf file and compute scale factors and mask.
  !!
  !! AUTHOR
  !!   * Origin  : Nicolas Grima (April 2007)
  !! 
  !! CREATION DATE
  !!   * April 2007
  !!
  !! HISTORY
  !!   Date (dd/mm/yyyy/) - Modification(s)
  !!
  !! ARGUMENTS
  !!   * No arguments
  !!
  !! TODO
  !!   
  !! USES
  !!   * sub_open_netcdf_file  (mod_netcdf.f90)
  !!   * sub_select_var_dims   (mod_netcdf.f90)
  !!   * sub_read_netcdf_var4d (mod_netcdf.f90)
  !!   * sub_close_netcdf_file (mod_netcdf.f90)
  !!
  !! USED BY
  !!   * trajec
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_input_grid_symphonie()

    !-------------!
    ! Declaration !
    !-------------!
    INTEGER(kind = iprec) :: ncid  ! netcdf file ID
    INTEGER(kind = iprec) :: varid ! netcdf variable ID
    INTEGER(kind = iprec) :: dimx  ! dimension in x (i)
    INTEGER(kind = iprec) :: dimy  ! dimension in y (j)
    INTEGER(kind = iprec) :: dimz  ! dimension in z (k)
    INTEGER(kind = iprec) :: dimt  ! dimension in t (l)

    ! dimsorder: array where netcdf variable dimensions order are stored.
    !            This, because sometime dimensions order could be different 
    !            than (x, y, z, t). It could be (x, y, t, nothing) or 
    !            (z, t, nothing, nothing).
    !            More information: see sub_select_var_dims and 
    !            sub_read_netcdf_var4d in mod_netcdf.f90 file.
    INTEGER(kind = iprec), DIMENSION(4) :: dimsorder
    INTEGER(kind = iprec):: kmax

    REAL(kind = rprec) , DIMENSION(:,:,:,:), ALLOCATABLE :: temp, vlx
!!NG    REAL(kind = rprec) , DIMENSION(:,:,:,:), ALLOCATABLE ::tmask_verif
    !! NG REAL(kind = rprec) , DIMENSION(:,:,:,:), ALLOCATABLE :: zz_tt
    REAL(kind = rprec) :: miss_val, miss_val_zo

    !-------------!
    ! Code begins !
    !-------------!

    !---------------------------------------!
    !-- Open NetCdf file with mesh values --!
    !---------------------------------------!
    CALL sub_open_netcdf_file(dir_grd_symp, fn_grd_symp, ncid)

    !====================== LON_T ================================!
    !-----------!
    !-- XX_TT --!
    !-----------!
    CALL sub_select_var_dims(ncid, nc_var_lon_t_symp, varid, &
         dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(xx_tt(1:dimx, 1:dimy,1,1))
    CALL sub_memory(size(xx_tt),'r','xx_tt','sub_input_grid_symphonie')
    CALL sub_read_netcdf_var4d(ncid, varid, xx_tt(:,:,:,:), dimsorder, dims_reg)
    WRITE(lun_standard,*)' - ', TRIM(nc_var_lon_t_symp), &
         ': max ', MAXVAL(xx_tt), ' min ', MINVAL(xx_tt)

!!NG    !! Save longitude grid for Symphonie data
!!NG    DO j = 1, dims_reg(2,3)
!!NG      WRITE(23,'(340(1x,f0.2))') (xx_tt(i,j,1,1), i= 1, dims_reg(1,3))
!!NG    ENDDO

    !====================== LON_U ================================!
    !-----------!
    !-- XX_UU --!
    !-----------!
    CALL sub_select_var_dims(ncid, nc_var_lon_u_symp, varid, &
         dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(xx_uu(1:dimx, 1:dimy, 1,1))
    CALL sub_memory(size(xx_uu),'r','xx_uu','sub_input_grid_symphonie')
    CALL sub_read_netcdf_var4d(ncid, varid, xx_uu(:,:,:,:), dimsorder, dims_reg)

    xx_uu(:,:,:,:)    = CSHIFT(xx_uu(:,:,:,:), shift=1, dim=1)
    xx_uu(dimx,:,:,:) = 9999._rprec

    WRITE(lun_standard,*)' - ', TRIM(nc_var_lon_u_symp), &
         ': max ', MAXVAL(xx_uu(dimx-1,:,:,:)), &
         '  min ', MINVAL(xx_uu(dimx-1,:,:,:))

    !====================== LAT_T ================================!
    !-----------!
    !-- YY_TT --!
    !-----------!
    CALL sub_select_var_dims(ncid, nc_var_lat_t_symp, varid, &
         dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(yy_tt(1:dimx, 1:dimy, 1, 1))
    CALL sub_memory(size(yy_tt),'r','yy_tt','sub_input_grid_symphonie')
    
    CALL sub_read_netcdf_var4d(ncid, varid, yy_tt(:,:,:,:), dimsorder, dims_reg)
    WRITE(lun_standard,*)' - ', TRIM(nc_var_lat_t_symp), &
         ': max ', MAXVAL(yy_tt), ' min ', MINVAL(yy_tt)

!!NG    !! Save latitude grid for Symphonie data
!!NG    DO j = 1, dims_reg(2,3)
!!NG      WRITE(24,'(340(1x,f0.2))') (yy_tt(i,j,1,1), i= 1, dims_reg(1,3))
!!NG    ENDDO

    !====================== LAT_V ================================!
    !-----------!
    !-- YY_VV --!
    !-----------!
    CALL sub_select_var_dims(ncid, nc_var_lat_v_symp, varid, &
         dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(yy_vv(1:dimx, 1:dimy, 1, 1))
    CALL sub_memory(size(yy_vv),'r','yy_vv','sub_input_grid_symphonie')
    CALL sub_read_netcdf_var4d(ncid, varid, yy_vv(:,:,1:dimz,:), dimsorder, dims_reg)

    yy_vv(:,:,:,:)    = CSHIFT(yy_vv(:,:,:,:), shift=1, dim=2)
    yy_vv(:, dimy,:,:) = 9999._rprec

    WRITE(lun_standard,*)' - ', TRIM(nc_var_lat_v_symp), &
         ': max ', MAXVAL(yy_vv(:, 1:dimy-1,:,:)), &
         '  min ', MINVAL(yy_vv(:, 1:dimy-1,:,:))


    !====================== bathy_m ================================!
    !-------------!
    !-- bathy_m --!
    !-------------!
    CALL sub_select_var_dims(ncid, 'bathy_m', varid, &
         dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(h_symp(1:dimx, 1:dimy, 1, 1))
    CALL sub_memory(size(h_symp),'r','h_symp','sub_input_grid_symphonie')
    CALL sub_read_netcdf_var4d(ncid, varid, h_symp(:,:,:,:), dimsorder, dims_reg)

    WRITE(lun_standard,*)' - bathy_m', &
         ': max ', MAXVAL(h_symp(:,:,:,:)), &
         '  min ', MINVAL(h_symp(:,:,:,:))

    !=============================================================!
    !========================= MASK ==============================!
    !=============================================================!
    !-----------!
    !-- TMASK --! Don't change during time evolution (verified).
    !-----------!
    !==================== READ X Velocity (3D)=======================!
    !! Pas tres propre!!!
    !!==================
    CALL sub_select_var_dims(ncid, TRIM(nc_var_zo), varid, &
         dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(vlx(1:dimx, 1:dimy, 1:dimz, 1))
    CALL sub_memory(size(vlx),'r','vlx','sub_input_grid_symphonie')
    CALL sub_read_netcdf_var4d(ncid, varid, &
         vlx(:,:,dimz-1:1:-1,:), dimsorder, dims_reg)
    !- Read mask value in netcdf file -!
    CALL sub_read_netcdf_att_val( &
         ncid                   , &
         varid                  , &
         TRIM(nc_att_mask_zo)   , &
         miss_val_zo            , &
         0._rprec)

    WHERE(vlx(:,:,:,:) == 0._rprec) vlx(:,:,:,:) = miss_val_zo

    !==================== READ TEMPERATURE (3D)=======================!
    !! Pas tres propre car on suppose que la temperature existe et 
    !! est dans le meme fichier que les coordonnees.
    !!==============================================
    CALL sub_select_var_dims(ncid, nc_var_te, varid, &
         dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(temp(1:dimx, 1:dimy, 1:dimz, 1))
    CALL sub_memory(size(temp),'r','temp','sub_input_grid_symphonie')
    CALL sub_read_netcdf_var4d(ncid, varid, &
         temp(:,:,dimz-1:1:-1,:), dimsorder, dims_reg)
    !- Read mask value in netcdf file -!
    CALL sub_read_netcdf_att_val( &
         ncid                   , &
         varid                  , &
         TRIM(nc_att_mask_te)   , &
         miss_val               , &
         0._rprec)

    temp(:,:,dimz,:) = temp(:,:,dimz-1,:)

    !==================== Compute TMASK (3D)======================!
    ALLOCATE(tmask(1:dims_reg(1,3), 1:dims_reg(2,3), 1:dims_reg(3,3), 1))
    CALL sub_memory(size(tmask),'r','tmask','sub_input_grid_symphonie')
    tmask(:,:,:,:) = 0._rprec

    kmax = dims_reg(3,3)

    WHERE((temp(:,:,kmax:2:-1,1) /= miss_val).AND. &
          (((temp(:,:,kmax:2:-1,1)-temp(:,:,kmax-1:1:-1,1)) > 1.e-11_rprec).OR. &
          (vlx(:,:,kmax:2:-1,1) /= miss_val_zo)) )
      tmask(:,:,kmax:2:-1,1) = 1._rprec
    END WHERE

    WHERE(temp(:,:,1,1) /= miss_val) tmask(:,:,1,1) = 1._rprec

!!NG    WHERE((temp(:,:,kmax:2:-1,dimt) /= miss_val).AND. &
!!NG          (((temp(:,:,kmax:2:-1,dimt)-temp(:,:,kmax-1:1:-1,dimt)) /= 0._rprec).OR. &
!!NG          (vlx(:,:,kmax:2:-1,dimt) /= miss_val_zo)) )
!!NG      tmask_verif(:,:,kmax:2:-1,1) = 1._rprec
!!NG    END WHERE
!!NG
!!NG    WHERE(temp(:,:,1,dimt) /= miss_val) tmask_verif(:,:,1,1) = 1._rprec
!!NG
!!NG    WRITE(0,*)'Verification du Tmask:',MAXVAL(ABS(tmask - tmask_verif))

    CALL sub_memory(-size(vlx),'r','vlx','sub_input_grid_symphonie')
    DEALLOCATE(vlx)
    CALL sub_memory(-size(temp),'r','temp','sub_input_grid_symphonie')
    DEALLOCATE(temp)
    

    WRITE(lun_standard,*)' - Tmask',&
         ': max ', MAXVAL(tmask), ' min ', MINVAL(tmask)

    !========= WMASK ============!
    ALLOCATE(wmask(1:dims_reg(1,3), 1:dims_reg(2,3), 1:dims_reg(3,3), 1))
    CALL sub_memory(size(wmask),'r','wmask','sub_input_grid_symphonie')
    wmask(:,:,:,:) = CSHIFT(tmask(:,:,:,:), shift=-1, dim=3)
    wmask(:,:,1,:) = wmask(:,:,2,:)

    !-----------------!
    !-- READ  ZZ_TT --!
    !-----------------!
!!$    ALLOCATE(zz_tt(1:dims_reg(1,3), 1:dims_reg(2,3), 1:dims_reg(3,3), 1))
!!$    CALL sub_memory(size(zz_tt),'r','zz_tt','sub_input_grid_symphonie')
!!$    zz_tt(:,:,:,:) = 0._rprec
!!$    CALL sub_select_var_dims(ncid, nc_var_depth_t_symp, varid, &
!!$         dimsorder, dimx, dimy, dimz, dimt)
!!$    CALL sub_read_netcdf_var4d(ncid, varid, zz_tt(:,:,dims_reg(3,3)-1:1:-1,:), dimsorder, dims_reg)
!!$
!!$    zz_tt(:,:,:,:) =  zz_tt(:,:,:,:) * tmask(:,:,:,:)
!!$
!!$    WRITE(lun_standard,*)' -  ', TRIM(nc_var_depth_t_symp), &
!!$         ' max ',MAXVAL(zz_tt), ' min ', MINVAL(zz_tt)
!!$
!!$
    ALLOCATE(zz_ww(1:dims_reg(1,3), 1:dims_reg(2,3), 1:dims_reg(3,3), 1))
    CALL sub_memory(size(zz_ww),'r','zz_ww','sub_input_grid_symphonie')
    zz_ww(:,:,:,:) = 0._rprec
!!$
!!$    kmax = dims_reg(3,3)
!!$
!!$    DO k = 2, kmax
!!$      DO j = 1, dims_reg(2,3)
!!$        DO i = 1, dims_reg(1,3)
!!$          zz_ww(i,j,k,1) = 2._rprec * zz_tt(i,j,k-1,1) - zz_ww(i,j,k-1,1)
!!$        ENDDO
!!$      ENDDO
!!$    ENDDO
!!$
!!$    zz_ww(:,:,:,:) =  zz_ww(:,:,:,:) * wmask(:,:,:,:)

!!    WRITE(0,*)MINLOC(zz_ww(:,:,:,:), mask = zz_ww(:,:,:,:) > 0._rprec)

    !! NG: CALL sub_memory(-size(zz_tt),'r','zz_tt','sub_input_grid_symphonie')
    !! NG: DEALLOCATE(zz_tt)
    !! NG: CALL sub_memory(-size(wmask),'r','wmask','sub_input_grid_symphonie')
    !! NG: DEALLOCATE(wmask)

    WRITE(lun_standard,*)' - ZW0: max ', MAXVAL(zz_ww), ' min ', MINVAL(zz_ww)

    !======================================================================!
    !========================= SCALE FACTORS ==============================!
    !======================================================================!

    !=========================== COMPUTE E2T ==============================!
    ALLOCATE(e2t(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(e2t),'r','e2t','sub_input_grid_symphonie')
    e2t(:,:,:,:) = cst_scale_factor * tmask(:,:,:,:)
    WRITE(lun_standard,*)' - e2t (c): max ', MAXVAL(e2t), ' min ', MINVAL(e2t)

    !=========================== COMPUTE E2U ==============================!
    ALLOCATE(e2u(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(e2u),'r','e2u','sub_input_grid_symphonie')
    e2u(:,:,:,:) = cst_scale_factor * tmask(:,:,:,:)
    WRITE(lun_standard,*)' - e2u (c): max ', MAXVAL(e2u), ' min ', MINVAL(e2u)

    !=========================== COMPUTE E1T ==============================!
    ALLOCATE(e1t(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(e1t),'r','e1t','sub_input_grid_symphonie')
    e1t(:,:,:,:) = cst_scale_factor * tmask(:,:,:,:)
    WRITE(lun_standard,*)' - e1t (c): max ', MAXVAL(e1t), ' min ', MINVAL(e1t)

    !=========================== COMPUTE E1V ==============================!
    ALLOCATE(e1v(1:dimx, 1:dimy, 1:dimz, 1:dimt))
    CALL sub_memory(size(e1v),'r','e1v','sub_input_grid_symphonie')
    e1v(:,:,:,:) = cst_scale_factor * tmask(:,:,:,:)
    WRITE(lun_standard,*)' - e1v (c): max ', MAXVAL(e1v), ' min ', MINVAL(e1v)

    !== E3T will be computed in  mod_input_data ==!

    !-------------------------------!
    !----- CLOSE NetCDF file -------!
    !-------------------------------!
    CALL sub_close_netcdf_file(ncid)

    PRINT *, ' '

  END SUBROUTINE sub_input_grid_symphonie
  !!***
  !=========================================================================
  !!****f* mod_input_grid/sub_input_B2C_gridz()
  !! NAME
  !!   sub_input_B2C_gridz()
  !!
  !! FUNCTION
  !!   * Read B-grid Model input grid coordinates, scale factors and mask from one
  !!     netcdf file.
  !!
  !! AUTHOR
  !!   * Origin  : Nicolas Grima (November 2008)
  !! 
  !! CREATION DATE
  !!   * November 2008
  !!
  !! HISTORY
  !!   Date (dd/mm/yyyy/) - Modification(s)
  !!
  !! ARGUMENTS
  !!   * No arguments
  !!
  !! TODO
  !!   
  !! USES
  !!
  !! USED BY
  !!   * trajec and trajec_seq
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_input_B2C_gridz()

    !-------------!
    ! Declaration !
    !-------------!
    INTEGER(kind = iprec) :: ncid  ! netcdf file ID
    INTEGER(kind = iprec) :: varid ! netcdf variable ID
    INTEGER(kind = iprec) :: dimx  ! dimension in x (i)
    INTEGER(kind = iprec) :: dimy  ! dimension in y (j)
    INTEGER(kind = iprec) :: dimz  ! dimension in z (k)
    INTEGER(kind = iprec) :: dimt  ! dimension in t (l)

    ! dimsorder: array where netcdf variable dimensions order are stored.
    !            This, because sometime dimensions order could be different 
    !            than (x, y, z, t). It could be (x, y, t, nothing) or 
    !            (z, t, nothing, nothing).
    !            More information: see sub_select_var_dims and 
    !            sub_read_netcdf_var4d in mod_netcdf.f90 file.
    INTEGER(kind = iprec), DIMENSION(4) :: dimsorder

    REAL(kind=rprec) :: att_mask_val

    !-------------!
    ! Code begins !
    !-------------!

    !---------------------------------------!
    !-- Open NetCdf file with mesh values --!
    !---------------------------------------!
    CALL sub_open_netcdf_file(dir_B2C_grid, file_name_B2C_grid, ncid)

    !======================================================================!
    !=========================== COORDINATES ==============================!
    !======================================================================!
    !-----------!
    !-- XX_TT --!
    !-----------!
    CALL sub_select_var_dims(ncid, nc_var_xx_tt, varid, dimsorder, dimx, dimy, dimz, dimt)
    if (key_add_bottom) THEN
       ALLOCATE(xx_tt(dimx,dimy,dimz+1,dimt))
       xx_tt(:,:,:,:)=0._rprec
    ELSE
       ALLOCATE(xx_tt(dimx,dimy,dimz,dimt))
    ENDIF
    CALL sub_memory(size(xx_tt),'r','xx_tt','sub_input_B2C_gridz')
    CALL sub_read_netcdf_var4d(ncid, varid, xx_tt(:,:,1:dimz,:), dimsorder, dims_reg)
    WRITE(lun_standard,*)' - ', TRIM(nc_var_xx_tt),': max ', MAXVAL(xx_tt), ' min ', MINVAL(xx_tt)

    !-----------!
    !-- XX_UU --!
    !-----------!
    CALL sub_select_var_dims(ncid, nc_var_xx_uu, varid, dimsorder, dimx, dimy, dimz, dimt)
    if (key_add_bottom) THEN
       ALLOCATE(xx_uu(dimx,dimy,dimz+1,dimt))
       xx_uu(:,:,:,:)=0._rprec
    ELSE
       ALLOCATE(xx_uu(dimx,dimy,dimz,dimt))
    ENDIF
    CALL sub_memory(size(xx_uu),'r','xx_uu','sub_input_B2C_gridz')
    CALL sub_read_netcdf_var4d(ncid, varid, xx_uu(:,:,1:dimz,:), dimsorder, dims_reg)
    WRITE(lun_standard,*)' - ', TRIM(nc_var_xx_uu),': max ', MAXVAL(xx_uu), ' min ', MINVAL(xx_uu)

    !-----------!
    !-- YY_TT --!
    !-----------!
    CALL sub_select_var_dims(ncid, nc_var_yy_tt, varid, dimsorder, dimx, dimy, dimz, dimt)
    if (key_add_bottom) THEN
       ALLOCATE(yy_tt(dimx,dimy,dimz+1,dimt))
       yy_tt(:,:,:,:)=0._rprec
    ELSE
       ALLOCATE(yy_tt(dimx,dimy,dimz,dimt))
    ENDIF
    CALL sub_memory(size(yy_tt),'r','yy_tt','sub_input_B2C_gridz')
    CALL sub_read_netcdf_var4d(ncid, varid, yy_tt(:,:,1:dimz,:), dimsorder, dims_reg)
    WRITE(lun_standard,*)' - ', TRIM(nc_var_yy_tt),': max ', MAXVAL(yy_tt), ' min ', MINVAL(yy_tt)

    !-----------!
    !-- YY_VV --!
    !-----------!
    CALL sub_select_var_dims(ncid, nc_var_yy_vv, varid, dimsorder, dimx, dimy, dimz, dimt)
    if (key_add_bottom) THEN
       ALLOCATE(yy_vv(dimx,dimy,dimz+1,dimt))
       yy_vv(:,:,:,:)=0._rprec
    ELSE
       ALLOCATE(yy_vv(dimx,dimy,dimz,dimt))
    ENDIF
    CALL sub_memory(size(yy_vv),'r','yy_vv','sub_input_B2C_gridz')
    CALL sub_read_netcdf_var4d(ncid, varid, yy_vv(:,:,1:dimz,:), dimsorder, dims_reg)
    WRITE(lun_standard,*)' - ', TRIM(nc_var_yy_vv),': max ', MAXVAL(yy_vv), ' min ', MINVAL(yy_vv)

    !-----------!
    !-- ZZ_WW --!
    !-----------!
    CALL sub_select_var_dims(ncid, nc_var_zz_ww, varid, dimsorder, dimx, dimy, dimz, dimt)
    if (key_add_bottom) THEN
       ALLOCATE(zz_ww(dimx,dimy,dimz+1,dimt))
       zz_ww(:,:,:,:)=0._rprec
    ELSE
       ALLOCATE(zz_ww(dimx,dimy,dimz,dimt))
    ENDIF
    CALL sub_memory(size(zz_ww),'r','zz_ww','sub_input_B2C_gridz')
    CALL sub_read_netcdf_var4d(ncid, varid, zz_ww(:,:,1:dimz,:), dimsorder, dims_reg)
    WRITE(lun_standard,*)' - ', TRIM(nc_var_zz_ww),': max ', MAXVAL(zz_ww), ' min ', MINVAL(zz_ww)

    !======================================================================!
    !========================= SCALE FACTORS ==============================!
    !======================================================================!
    !---------!
    !-- E1F --!
    !---------!
    CALL sub_select_var_dims(ncid, nc_var_e1f, varid, dimsorder, dimx, dimy, dimz, dimt)
    ALLOCATE(e1f(dimx,dimy,dimz,dimt))
    CALL sub_memory(size(e1f),'r','e1f','sub_input_B2C_gridz')
    CALL sub_read_netcdf_var4d(ncid, varid, e1f(:,:,1:dimz,:), dimsorder, dims_reg)
    WRITE(lun_standard,*)' - ', TRIM(nc_var_e1f),': max ', MAXVAL(e1f), ' min ', MINVAL(e1f)

    !---------!
    !-- E2F --!
    !---------!
    CALL sub_select_var_dims(ncid, nc_var_e2f, varid, dimsorder, dimx, dimy, dimz, dimt)
    if (key_add_bottom) THEN
       ALLOCATE(e2f(dimx,dimy,dimz+1,dimt))
       e2f(:,:,:,:)=0._rprec
    ELSE
       ALLOCATE(e2f(dimx,dimy,dimz,dimt))
    ENDIF
    CALL sub_memory(size(e2f),'r','e2f','sub_input_B2C_gridz')
    CALL sub_read_netcdf_var4d(ncid, varid, e2f(:,:,1:dimz,:), dimsorder, dims_reg)
    WRITE(lun_standard,*)' - ', TRIM(nc_var_e2f),': max ', MAXVAL(e2f), ' min ', MINVAL(e2f)

    !---------!
    !-- E3F --!
    !---------!
    CALL sub_select_var_dims(ncid, nc_var_e3f, varid, dimsorder, dimx, dimy, dimz, dimt)
    if (key_add_bottom) THEN
       ALLOCATE(e3f(dimx, dimy, dimz+1, dimt))
       e3f(:,:,:,:)=1.e+10_rprec
    ELSE
       ALLOCATE(e3f(dimx, dimy, dimz, dimt))
       
    ENDIF
    CALL sub_memory(size(e3f),'r','e3f','sub_input_B2C_gridz')
    CALL sub_read_netcdf_var4d(ncid, varid, e3f(:,:,1:dimz,:), dimsorder, dims_reg)

    !! NG: DEBUG MOM
    !! e3f(:,:,1,:)=10._rprec
    !! NG: DEBUG MOM

    !- Read mask value in netcdf file -!
    CALL sub_read_netcdf_att_val( &
         ncid                   , &
         varid                  , &
         'missing_value'        , &
         att_mask_val           , &
         nc_mask_val              )

    WHERE(e3f(:,:,:,:) == att_mask_val) e3f(:,:,:,:) = 1.e+10_rprec

    WRITE(lun_standard,*)' - ', TRIM(nc_var_e3f),': max ', MAXVAL(e3f), ' min ', MINVAL(e3f)

    !---------!
    !-- E2U --!
    !---------!
    CALL sub_select_var_dims(ncid, nc_var_e2u, varid, dimsorder, dimx, dimy, dimz, dimt)
    if (key_add_bottom) THEN
       ALLOCATE(e2u(dimx,dimy,dimz+1,dimt))
       e2u(:,:,:,:)=0._rprec
    ELSE
       ALLOCATE(e2u(dimx,dimy,dimz,dimt))
    ENDIF
    CALL sub_memory(size(e2u),'r','e2u','sub_input_B2C_gridz')
    CALL sub_read_netcdf_var4d(ncid, varid, e2u(:,:,1:dimz,:), dimsorder, dims_reg)
    WRITE(lun_standard,*)' - ', TRIM(nc_var_e2u),': max ', MAXVAL(e2u), ' min ', MINVAL(e2u)

    !---------!
    !-- E1V --!
    !---------!
    CALL sub_select_var_dims(ncid, nc_var_e1v, varid, dimsorder, dimx, dimy, dimz, dimt)
    if (key_add_bottom) THEN
       ALLOCATE(e1v(dimx,dimy,dimz+1,dimt))
       e1v(:,:,:,:)=0._rprec
    ELSE
       ALLOCATE(e1v(dimx,dimy,dimz,dimt))
    ENDIF
    CALL sub_memory(size(e1v),'r','e1v','sub_input_B2C_gridz')
    CALL sub_read_netcdf_var4d(ncid, varid, e1v(:,:,1:dimz,:), dimsorder, dims_reg)
    WRITE(lun_standard,*)' - ', TRIM(nc_var_e1v),': max ', MAXVAL(e1v), ' min ', MINVAL(e1v)

    !---------!
    !-- E1T --!
    !---------!
    CALL sub_select_var_dims(ncid, nc_var_e1t, varid, dimsorder, dimx, dimy, dimz, dimt)
    if (key_add_bottom) THEN
       ALLOCATE(e1t(dimx,dimy,dimz+1,dimt))
       e1t(:,:,:,:)=0._rprec
    ELSE
       ALLOCATE(e1t(dimx,dimy,dimz,dimt))
    ENDIF
    CALL sub_memory(size(e1t),'r','e1t','sub_input_B2C_gridz')
    CALL sub_read_netcdf_var4d(ncid, varid, e1t(:,:,1:dimz,:), dimsorder, dims_reg)
    WRITE(lun_standard,*)' - ', TRIM(nc_var_e1t),': max ', MAXVAL(e1t), ' min ', MINVAL(e1t)

    !---------!
    !-- E2T --!
    !---------!
    CALL sub_select_var_dims(ncid, nc_var_e2t, varid, dimsorder, dimx, dimy, dimz, dimt)
    if (key_add_bottom) THEN
       ALLOCATE(e2t(dimx,dimy,dimz+1,dimt))
       e2t(:,:,:,:)=0._rprec
    ELSE
       ALLOCATE(e2t(dimx,dimy,dimz,dimt))
    ENDIF
    CALL sub_memory(size(e2t),'r','e2t','sub_input_B2C_gridz')
    CALL sub_read_netcdf_var4d(ncid, varid, e2t(:,:,1:dimz,:), dimsorder, dims_reg)
    WRITE(lun_standard,*)' - ', TRIM(nc_var_e2t),': max ', MAXVAL(e2t), ' min ', MINVAL(e2t)

    !=============================================================!
    !========================= MASK ==============================!
    !=============================================================!
    !-----------!
    !-- TMASK --!
    !-----------!
    CALL sub_select_var_dims(ncid, nc_var_tmask, varid, dimsorder, dimx, dimy, dimz, dimt)
    if (key_add_bottom) THEN
       ALLOCATE(tmask(dimx,dimy,dimz+1,dimt))
       tmask(:,:,:,:)=0._rprec
    ELSE
       ALLOCATE(tmask(dimx,dimy,dimz,dimt))
    ENDIF
    CALL sub_memory(size(tmask),'r','tmask','sub_input_B2C_gridz')
    CALL sub_read_netcdf_var4d(ncid, varid, tmask(:,:,1:dimz,:), &
         dimsorder, dims_reg)

    !- Read mask value in netcdf file -!
    CALL sub_read_netcdf_att_val( &
         ncid                   , &
         varid                  , &
         'missing_value'        , &
         att_mask_val           , &
         nc_mask_val              )

    WHERE(tmask(:,:,:,:) /= att_mask_val) tmask(:,:,:,:) = 1._rprec

    WHERE(tmask(:,:,:,:) == att_mask_val) tmask(:,:,:,:) = 0._rprec

    IF ( SUM(ABS(tmask(:,:,:,:))) == 0) THEN
      WRITE(lun_error,*) 'mod_input_grid: sub_input_grid: problem with tmask values!'
      STOP
    ENDIF


    !---------!
    !-- E3T --!
    !---------!
    CALL sub_select_var_dims(ncid, nc_var_e3t, varid, dimsorder, dimx, dimy, dimz, dimt)
    if (key_add_bottom) THEN
       ALLOCATE(e3t(dimx, dimy, dimz+1, dimt))
       e3t(:,:,:,:)=0._rprec
    ELSE
       ALLOCATE(e3t(dimx, dimy, dimz, dimt))
    ENDIF
    CALL sub_memory(size(e3t),'r','e3t','sub_input_B2C_gridz')
    CALL sub_read_netcdf_var4d(ncid, varid, e3t(:,:,1:dimz,:), dimsorder, dims_reg)

    WHERE(e3t(:,:,:,:) == att_mask_val) e3t(:,:,:,:) = 1.e+10_rprec
    
    WRITE(lun_standard,*)' - ', TRIM(nc_var_e3t),': max ', MAXVAL(e3t), ' min ', MINVAL(e3t)


    !- Tmask statistics -!

    WRITE(lun_standard,*)' - ', TRIM(nc_var_tmask),': max ', MAXVAL(tmask), ' min ', MINVAL(tmask)


    !-------------------------------!
    !----- CLOSE NetCDF file -------!
    !-------------------------------!
    CALL sub_close_netcdf_file(ncid)

    PRINT *, ' '

    !! NG: for MOM
    !! NG: zz_ww(:,:,:,:)=CSHIFT(zz_ww(:,:,:,:),shift=-1,dim=3)
    !! NG: zz_ww(:,:,1,:)=0._rprec
    !! NG: WRITE(lun_standard,*)' DEBUG --- zz_ww MOM: ',zz_ww
    !! NG: PRINT *, ' '
    !!

  END SUBROUTINE sub_input_B2C_gridz
  !!***
  !=========================================================================
  !!****f* mod_input_grid/sub_coord_dealloc()
  !! NAME
  !!   sub_coor_dealloc()
  !!
  !! FUNCTION
  !!   * Deallocate memory for coordinate arrays
  !!
  !! AUTHOR
  !!   * Origin  : Nicolas Grima (April-May 2005)
  !! 
  !! CREATION DATE
  !!   * April-May 2005
  !!
  !! HISTORY
  !!   Date (dd/mm/yyyy/) - Modification(s)
  !!
  !! ARGUMENTS
  !!   * No argument
  !!
  !! TODO
  !!   
  !!
  !! USED BY
  !!   * trajec
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_coord_dealloc()

    !-------------!
    ! Code begins !
    !-------------!
    !- Deallocate arrays memory -!
    IF (ALLOCATED(xx_tt)) THEN
       CALL sub_memory(-size(xx_tt),'r','xx_tt','sub_coord_dealloc')
       DEALLOCATE(xx_tt)
    ENDIF
    IF (ALLOCATED(xx_uu)) THEN
       CALL sub_memory(-size(xx_uu),'r','xx_uu','sub_coord_dealloc')
       DEALLOCATE(xx_uu)
    ENDIF
    IF (ALLOCATED(yy_tt)) THEN
       CALL sub_memory(-size(yy_tt),'r','yy_tt','sub_coord_dealloc')
       DEALLOCATE(yy_tt)
    ENDIF
    IF (ALLOCATED(yy_vv)) THEN
       CALL sub_memory(-size(yy_vv),'r','yy_vv','sub_coord_dealloc')
       DEALLOCATE(yy_vv)
    ENDIF
    IF (ALLOCATED(zz_ww)) THEN
       CALL sub_memory(-size(zz_ww),'r','zz_ww','sub_coord_dealloc')
       DEALLOCATE(zz_ww)
    ENDIF

  END SUBROUTINE sub_coord_dealloc
  !!***
  !=========================================================================
  !!****f* mod_input_grid/sub_scalef_dealloc()
  !! NAME
  !!   sub_scalef_dealloc()
  !!
  !! FUNCTION
  !!   * Deallocate memory of scale factor arrays
  !!
  !! AUTHOR
  !!   * Origin  : Nicolas Grima (April-May 2005)
  !! 
  !! CREATION DATE
  !!   * April-May 2005
  !!
  !! HISTORY
  !!   Date (dd/mm/yyyy/) - Modification(s)
  !!
  !! ARGUMENTS
  !!   * No argument
  !!
  !! TODO
  !!   
  !!
  !! USED BY
  !!   * trajec
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_scalef_dealloc()

    !-------------!
    ! Code begins !
    !-------------!
    !- Deallocate arrays memory -!
    IF (ALLOCATED(e1t)) THEN
       CALL sub_memory(-size(e1t),'r','e1t','sub_scalef_dealloc')
       DEALLOCATE(e1t)
    ENDIF
    IF (ALLOCATED(e2t)) THEN
       CALL sub_memory(-size(e2t),'r','e2t','sub_scalef_dealloc')
       DEALLOCATE(e2t)
    ENDIF
    IF (ALLOCATED(e2u)) THEN
       CALL sub_memory(-size(e2u),'r','e2u','sub_scalef_dealloc')
       DEALLOCATE(e2u)
    ENDIF
    IF (ALLOCATED(e1v)) THEN
       CALL sub_memory(-size(e1v),'r','e1v','sub_scalef_dealloc')
       DEALLOCATE(e1v)
    ENDIF
    IF (ALLOCATED(e3t)) THEN
       CALL sub_memory(-size(e3t),'r','e3t','sub_scalef_dealloc')
       DEALLOCATE(e3t)
    ENDIF

  END SUBROUTINE sub_scalef_dealloc

  !!***
  !=========================================================================
  !!****f* mod_input_grid/sub_h_roms_dealloc()
  !! NAME
  !!   sub_h_roms_dealloc()
  !!
  !! FUNCTION
  !!   * Deallocate memory of h_roms arrays
  !!
  !! AUTHOR
  !!   * Origin  : Nicolas Grima (November 2006)
  !! 
  !! CREATION DATE
  !!   * November 2006
  !!
  !! HISTORY
  !!   Date (dd/mm/yyyy/) - Modification(s)
  !!
  !! ARGUMENTS
  !!   * No argument
  !!
  !! TODO
  !!   
  !!
  !! USED BY
  !!   * trajec and trajec_seq
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_h_roms_dealloc()

    !-------------!
    ! Code begins !
    !-------------!
    !- Deallocate arrays memory -!
    IF (ALLOCATED(h_roms)) THEN
       CALL sub_memory(-size(h_roms),'r','h_roms','sub_h_roms_dealloc')
       DEALLOCATE(h_roms)
    ENDIF

  END SUBROUTINE sub_h_roms_dealloc
  !!***
  !=========================================================================
  !!****f* mod_input_grid/sub_tmask_dealloc()
  !! NAME
  !!   sub_tmask_dealloc()
  !!
  !! FUNCTION
  !!   Deallocate memory of T mask array
  !!
  !! AUTHOR
  !!   * Origin  : Nicolas Grima (April-May 2005)
  !! 
  !! CREATION DATE
  !!   * April-May 2005
  !!
  !! HISTORY
  !!   Date (dd/mm/yyyy/) - Modification(s)
  !!
  !! ARGUMENTS
  !!   * No argument
  !!
  !! TODO
  !!   
  !!
  !! USED BY
  !!   * trajec
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_tmask_dealloc()

    !-------------!
    ! Code begins !
    !-------------!
    !- Deallocate array memory -!
    CALL sub_memory(-size(tmask),'r','tmask','sub_tmask_dealloc')
    DEALLOCATE(tmask)

  END SUBROUTINE sub_tmask_dealloc
  !!***
END MODULE mod_input_grid