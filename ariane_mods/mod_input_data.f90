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
!!****h* ariane/mod_input_data
!! NAME
!!   mod_input_data (mod_input_data.f90 - Fortran90 module)
!!
!! USAGE
!!   Include 'USE mod_input_data' in the header of your Fortran 90 source 
!!   code.
!!   Then you'll have access to the subroutine:
!!      - sub_input_data
!!      - sub_input_data_main
!!      - sub_transp_alloc
!!      - sub_transp_dealloc
!!      - sub_tracer_alloc
!!      - sub_tracer_dealloc
!!
!! FUNCTION
!!   Read input data: tracers and currents (compute transport).
!!   File format is netcdf (netcdf version 3.6.0 or newer is required).
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
!! RESULT
!!   
!!
!! EXAMPLES
!!   * USE mod_input_data
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
!!   * i686-pc-linux-gnu -          g95
!!
!! SEE ALSO
!!   
!!
!! USES
!!   * USE mod_precision
!!   * USE mod_namelist
!!   * USE mod_input_grid
!!   * USE mod_w_comput
!!   * USE mod_rhostp
!!   * USE mod_netcdf
!!   * USE reducmem
!!
!! USED BY
!!   * posini
!!   * trajec
!!
!! SOURCE
!!=========================================================================
MODULE mod_input_data

  !------------------!
  ! USE ASSOCIAITION !
  !------------------!
  USE mod_precision
  USE mod_namelist
  USE mod_input_grid
  USE mod_seq
  USE mod_w_comput
  USE mod_rhostp
  !! USE mod_netcdf
  USE mod_netcdf_write_fast
  USE mod_reducmem
  USE netcdf
  USE mod_B2C_grid_interpolation
  USE mod_B2C_grid_save
  USE mod_input_data_main

  !-------------!
  ! DECLARATION !
  !-------------!
  IMPLICIT NONE

  !- Global variables -!
  REAL(kind=rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: &
       uu      , & ! Zonal Transport
       vv      , & ! Meridional transport
       ww          ! Vertical transport

  REAL(kind=rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: &
       epr       ! + Evaporation - Precipitaiton - Runoff

  REAL(kind=rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: &
       tt      , & ! Temperature
       ss      , & ! Salinity
       rr          ! Density

  REAL(kind=rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: &
       ssh_vvl         ! sea surface height (for VVL)

  !! NG: 05/08/2008 : global variables for ROMS
  !! If ROMS is not selected these arrays don't take memory...
  REAL(kind=rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: &
       zeta_a, &
       zw0 , &
       e3u , &
       e3v , &
       zz_tt, &
       sses

  INTEGER(kind = iprec), DIMENSION(:,:,:,:), ALLOCATABLE :: &
       uBgridmask, &
       vBgridmask

CONTAINS
  !!***
  !=========================================================================
  !!****f* mod_input_data/sub_input_data()
  !! NAME
  !!   sub_input_data()
  !!
  !! FUNCTION
  !!   * loop on current and tracer variables: U, V, [W, T, S, R].
  !!   * call read netcdf.
  !!   * compute transport from currents.
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
  !!   * no arguments
  !!
  !! TODO
  !!   
  !!
  !! USED BY
  !!   * trajec
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_input_data()

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)'========================================================'
    WRITE(lun_standard,*)'= READ INPUT GEOPHYSICAL DATA AND STORE THEM IN MEMORY ='
    WRITE(lun_standard,*)'========================================================'

    ! allocate uu, vv, ww
    CALL sub_transp_alloc(dims_reg(1,3),dims_reg(2,3),dims_reg(3,3),dims_reg(4,3))
    uu(:,:,:,:) = 0._rprec
    vv(:,:,:,:) = 0._rprec
    ww(:,:,:,:) = 0._rprec

    IF (key_alltracers) THEN
      ! allocate tt, ss, rr
      CALL sub_tracer_alloc(dims_reg(1,3),dims_reg(2,3),dims_reg(3,3),dims_reg(4,3))
      tt(:,:,:,:) = 0._rprec
      ss(:,:,:,:) = 0._rprec
      rr(:,:,:,:) = 0._rprec
    ENDIF

    if (key_vvl) then
       ! allocate ssh
       CALL sub_ssh_alloc(dims_reg(1,3),dims_reg(2,3))
       ssh_vvl(:,:,1,1) = 0._rprec
    endif

    IF (key_roms) THEN

      CALL sub_input_data_roms()

    ELSEIF (key_mars) THEN

      !!- TODO - CALL sub_input_data_mars()
      STOP

    ELSEIF (key_symphonie) THEN

      CALL sub_input_data_symphonie()

    ELSEIF (key_B2C_grid) THEN

      IF (B2C_grid_Z_or_Sigma=='Z') THEN
        CALL sub_input_data_B2C_gridz()
      ELSE
        STOP
        !!TODO:  CALL sub_input_data_B2CD_grids()
      ENDIF

    ELSE

      CALL sub_input_data_opa()

    ENDIF

  END SUBROUTINE sub_input_data
  !!***
  !!***
  !=========================================================================
  !!****f* mod_input_data_seq/sub_input_data_seq()
  !! NAME
  !!   sub_input_data_seq()
  !!
  !! FUNCTION
  !!   * loop on current and tracer variables: U, V, [W, T, S, R].
  !!   * call read netcdf.
  !!   * compute transport from currents.
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
  !!   * no arguments
  !!
  !! TODO
  !!   
  !!
  !! USED BY
  !!   * trajec
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_input_data_alloc_init_seq()

    INTEGER(kind = iprec) :: alloc_size

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)'====================================================='
    WRITE(lun_standard,*)'= ALLOCATE INPUT GEOPHYSICAL DATA (sequential mode) ='
    WRITE(lun_standard,*)'====================================================='

    ! allocate uu, vv, ww
    CALL sub_transp_alloc(dims_reg(1,3),dims_reg(2,3),dims_reg(3,3),it_ind)
    uu(:,:,:,:) = 0._rprec
    vv(:,:,:,:) = 0._rprec
    ww(:,:,:,:) = 0._rprec

    IF (key_alltracers) THEN
      ! allocate tt, ss, rr
      CALL sub_tracer_alloc(dims_reg(1,3),dims_reg(2,3),dims_reg(3,3),it_ind)
      tt(:,:,:,:) = 0._rprec
      ss(:,:,:,:) = 0._rprec
      rr(:,:,:,:) = 0._rprec
    ENDIF

    if (key_vvl) then
       ! allocate ssh
       CALL sub_ssh_alloc(dims_reg(1,3),dims_reg(2,3))
       ssh_vvl(:,:,1,1) = 0._rprec
    endif
    !!NG: This is done if key_interp_temporal = .TRUE.
    !!NG: to read a first input data.
    !!NG: the last  input data time step in forward mode.
    !!NG: the first input data time step in backward mode.
    IF (key_interp_temporal) THEN

      IF (key_alltracers) THEN
        IF ((key_roms).OR.(key_mars).OR.(key_symphonie))THEN
          alloc_size = 7
        ELSE
          IF (TRIM(w_surf_option) == 'E-P-R') THEN
             alloc_size = 7
          elseif (key_vvl) then
                write (*,*) 'key_vvl and key_alltracers not written'
                stop
          ELSE
            alloc_size = 6
          ENDIF
        ENDIF
      ELSE
        IF (TRIM(w_surf_option) == 'E-P-R') THEN
           alloc_size = 4
        elseif (key_vvl) then
             alloc_size = 4
        ELSE
          alloc_size = 3
        ENDIF
      ENDIF

      CALL sub_seq_alloc(alloc_size)
      !!NG      CALL sub_seq_interp_init()

      IF (TRIM(forback) == 'forward' ) THEN
        forback = 'backward'
      ELSE
        forback = 'forward'
      ENDIF

      CALL sub_seq_init()

      IF (key_roms) THEN
        CALL sub_input_data_seq_roms( &
             ncids(:)               , &
             varids(:)              , &
             new_file(:)            , &
             ind_file(:)            , &
             ind_time(:)            , &
             ind_time_size(:)       , &
             sdimsorders(:,:)        )

      ELSEIF (key_mars) THEN
        CALL sub_input_data_seq_mars( &
             ncids(:)               , &
             varids(:)              , &
             new_file(:)            , &
             ind_file(:)            , &
             ind_time(:)            , &
             ind_time_size(:)       , &
             sdimsorders(:,:)        )

      ELSEIF ( key_symphonie) THEN
        CALL sub_input_data_seq_symphonie( &
             ncids(:)                    , &
             varids(:)                   , &
             new_file(:)                 , &
             ind_file(:)                 , &
             ind_time(:)                 , &
             ind_time_size(:)            , &
             sdimsorders(:,:)            )

      ELSEIF (key_B2C_grid) THEN

        WRITE(lun_error,*)''
        WRITE(lun_error,*)'STOP === key_B2C_grid is not available === STOP'
        STOP

      ELSE
        CALL sub_input_data_seq_opa( &
             ncids(:)              , &
             varids(:)             , &
             new_file(:)           , &
             ind_file(:)           , &
             ind_time(:)           , &
             ind_time_size(:)      , &
             sdimsorders(:,:)       )
      ENDIF ! (key_roms)

      IF (TRIM(forback) == 'forward' ) THEN
        forback = 'backward'
      ELSE
        forback = 'forward'
      ENDIF

    ENDIF ! interpolation

  END SUBROUTINE sub_input_data_alloc_init_seq
  !!***

  !=========================================================================
  !!****f* mod_input_data/sub_input_data_opa()
  !! NAME
  !!   sub_input_data_opa()
  !!
  !! FUNCTION
  !!   * read OPA tracer variables.
  !!   * loop on current and tracer variables: U, V, [W, T, S, R].
  !!   * call read netcdf.
  !!   * compute transport from currents.
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
  !!   * no arguments
  !!
  !! TODO
  !!   
  !!
  !! USED BY
  !!   * trajec
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_input_data_opa()

    ! Variable number: U, V, [W, T, S, R]
    CHARACTER(len=4), PARAMETER :: c_none='NONE'

    INTEGER(kind  = iprec) :: kk, ll
    INTEGER(kind  = iprec) :: i,j,k,l

    !-------------!
    ! Code begins !
    !-------------!

     !======================= SEA SURFACE HEIGHT ===============!

    if (key_vvl) then
       CALL sub_input_data_main(c_dir_ssh,c_prefix_ssh,ind0_ssh,indn_ssh, &
            maxsize_ssh,c_suffix_ssh,nc_var_ssh,c_none,nc_att_mask_ssh,ssh_vvl(:,:,:,:))

       if (.not. key_sequential) then
          write (*,*) "key_vvl currently requires key_sequential"
          stop
       endif

       e3t(:, :, kk, 1) = e3t0(:, :, kk, 1) * (1 + ssh_vvl(:, :, 1, 1)/totaldepth(:, :,1,1))
    endif


    !======================= ZONAL TRANSPORT ==================!
    !-- Read Zonal Current --!
    CALL sub_input_data_main(c_dir_zo,c_prefix_zo,ind0_zo,indn_zo,maxsize_zo, &
         c_suffix_zo, nc_var_zo,nc_var_eivu,nc_att_mask_zo,uu(:,:,:,:))

    !-- Compute transport --!
    IF (key_partialsteps) THEN

      DO ll = 1, lmt
        DO kk = 1, dims_reg(3,3)
          ! partial steps case e3t is a 3D array
          uu(:,:,kk,ll) = uu(:,:,kk,ll) * e2u(:,:,1,1) * &
               MIN(e3t(:,:,kk,1), CSHIFT(e3t(:,:,kk,1),shift=1,dim=1))
        ENDDO
      ENDDO

    ELSE

      DO ll = 1, lmt
        DO kk = 1, dims_reg(3,3)
          ! e3t is a vector (no partial steps)
          uu(:,:,kk,ll) = uu(:,:,kk,ll) * e2u(:,:,1,1) * e3t(1,1,kk,1)
        ENDDO
      ENDDO
    ENDIF

    WRITE(lun_standard,*)' - Transport: max ', MAXVAL(uu), ' min ', MINVAL(uu)

    !============= MERIDIONAL TRANSPORT =======================!

    !-- Read Meridional Current --!
    CALL sub_input_data_main(c_dir_me,c_prefix_me,ind0_me,indn_me,maxsize_me, &
         c_suffix_me, nc_var_me,nc_var_eivv,nc_att_mask_me,vv(:,:,:,:))

    !-- Compute transport --
    IF (key_partialsteps) THEN

      DO ll = 1, lmt
        DO kk = 1, dims_reg(3,3)
          ! partial steps case e3t is a 3D array
          vv(:,:,kk,ll) = vv(:,:,kk,ll) * e1v(:,:,1,1) * &
               MIN(e3t(:,:,kk,1), CSHIFT(e3t(:,:,kk,1),shift=1,dim=2))
        ENDDO
      ENDDO

    ELSE

      DO ll = 1, lmt
        DO kk = 1, dims_reg(3,3)
          ! e3t is a vector (no partial steps)
          vv(:,:,kk,ll) = vv(:,:,kk,ll) * e1v(:,:,1,1) * e3t(1,1,kk,1)
        ENDDO
      ENDDO

    ENDIF

    WRITE(lun_standard,*)' - Transport: max ', MAXVAL(vv), ' min ', MINVAL(vv)

    !=============== VERTICAL TRANSPORT =======================!
    IF (key_computew) THEN

      !-- compute w transport from U and V transport --!
      CALL sub_w_comput(uu,vv,ww)

    ELSE

      !-- Read Vertical Current --!
      CALL sub_input_data_main(c_dir_ve,c_prefix_ve,ind0_ve,indn_ve,maxsize_ve, &
           c_suffix_ve, nc_var_ve,nc_var_eivw,nc_att_mask_ve,ww(:,:,:,:))

      !-- Compute transport --
      DO ll = 1, lmt
        DO kk = 1, dims_reg(3,3)
          ww(:,:,kk,ll) = ww(:,:,kk,ll) * e1t(:,:,1,1) * e2t(:,:,1,1)
        ENDDO
      ENDDO

    ENDIF

    !!NG: 24/07/2009 
    IF (TRIM(w_surf_option) == 'E-P-R') THEN
      !!NG: 24/07/2009 We force here W="E-P-R" at the surface.
      !!NG: 24/07/2009 Physically is the better solution,
      !!NG: 24/07/2009 but in this case the non-divergence criteria
      !!NG: 24/07/2009 is not at all respected at the sea surface. 
      !!NG: 24/07/2009 This option has to be activated
      !!NG: 24/07/2009 if you know what you are doing.

      !-- Read E-R-R --!
      CALL sub_input_data_main(c_dir_ep,c_prefix_ep,ind0_ep,indn_ep,maxsize_ep, &
           c_suffix_ep, nc_var_ep,c_none,nc_att_mask_ep,epr(:,:,:,:))

      DO ll = 1, lmt
        ww(:,:,1,ll) = epr(:,:,1,ll) * epr_coef * e1t(:,:,1,1) * e2t(:,:,1,1)
      ENDDO

    ELSEIF (TRIM(w_surf_option) == 'zero') THEN
      !!NG: 24/07/2009 We force here W=0 at the surface.
      !!NG: 24/07/2009 Is the same case that the precedent, but you don't
      !!NG: 24/07/2009 have the "E-P-R" fields available.
      !!NG: 24/07/2009 This option has to be activated
      !!NG: 24/07/2009 if you know what you are doing.
      ww(:,:,1,:) = 0._rprec
    ELSE
      !!NG: 24/07/2009 Nothing is changed compared to the precedent version
      !!NG: 24/07/2009 of Ariane. The non-divergence criteria is respected,
      !!NG: 24/07/2009 but a lot of particles can be intercepted at the surface.
      !!NG: 24/07/2009 To reduce the number of particles intercepted at the sea
      !!NG: 24/07/2009 surface, please activate one the two options above mentioned.
    ENDIF

    WRITE(lun_standard,*)' - Transport: max ', MAXVAL(ww), ' min ', MINVAL(ww)

    IF (key_alltracers) THEN

      !==================== TEMPERATURE =========================!
      !-- Read Temperature --!
      CALL sub_input_data_main(c_dir_te,c_prefix_te,ind0_te,indn_te,maxsize_te, &
           c_suffix_te, nc_var_te,c_none,nc_att_mask_te,tt(:,:,:,:))

      !====================== SALINITY ===========================!
      !-- Read Salinity --!
      CALL sub_input_data_main(c_dir_sa,c_prefix_sa,ind0_sa,indn_sa,maxsize_sa, &
           c_suffix_sa, nc_var_sa,c_none,nc_att_mask_sa,ss(:,:,:,:))

      !======================= DENSITY ===========================!
      IF (key_computesigma) THEN

        !-- Compute Density from Temperature and Salinity --!
        DO l = 1, lmt
          DO k = 1,  dims_reg(3,3)
            DO j = 1,  dims_reg(2,3)
              DO i = 1,  dims_reg(1,3)
                rr(i,j,k,l)=rhostp(1000._rprec,ss(i,j,k,l),tt(i,j,k,l),-ABS(zsigma))
              END DO
            END DO
          END DO
        ENDDO

      ELSE

        !-- Read Density --!
        CALL sub_input_data_main(c_dir_de,c_prefix_de,ind0_de,indn_de,maxsize_de, &
             c_suffix_de, nc_var_de,c_none,nc_att_mask_de,rr(:,:,:,:))

      ENDIF

    ENDIF

  END SUBROUTINE sub_input_data_opa
  !!***
  !=========================================================================
  !!****f* mod_input_data_seq/sub_input_data_seq_opa()
  !! NAME
  !!   sub_input_data_seq_opa()
  !!
  !! FUNCTION
  !!   * read OPA tracer variables.
  !!   * loop on current and tracer variables: U, V, [W, T, S, R].
  !!   * call read netcdf.
  !!   * compute transport from currents.
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
  !!   * no arguments
  !!
  !! TODO
  !!   
  !!
  !! USED BY
  !!   * trajec
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_input_data_seq_opa( &
       ncids                       , &
       varids                      , &
       new_file                    , &
       ind_file                    , &
       ind_time                    , &
       ind_time_size               , &
       dimsorders                  )


    INTEGER(kind = iprec) , DIMENSION(:), INTENT(inout) :: ncids
    INTEGER(kind = iprec) , DIMENSION(:), INTENT(inout) :: varids
    LOGICAL               , DIMENSION(:), INTENT(inout) :: new_file
    INTEGER(kind = iprec) , DIMENSION(:), INTENT(inout) :: ind_file
    INTEGER(kind = iprec) , DIMENSION(:), INTENT(inout) :: ind_time
    INTEGER(kind = iprec) , DIMENSION(:), INTENT(inout) :: ind_time_size
    INTEGER(kind = iprec) , DIMENSION(:,:), INTENT(inout) :: dimsorders

    ! Variable number: U, V, [W, T, S, R]
    CHARACTER(len=4), PARAMETER :: c_none='NONE'

    INTEGER(kind  = iprec) :: kk, ll
    INTEGER(kind  = iprec) :: i,j,k,l

    INTEGER(kind  = iprec) :: ind_ep, ind_ssh


    !-------------!
    ! Code begins !
    !-------------!

    !======================= SEA SURFACE HEIGHT ===============!
    if (key_vvl) then
       if (key_alltracers) then
          ind_ssh = 7
       else
          ind_ssh = 4
       endif
       CALL sub_input_data_seq_main(                      &
         ncids(ind_ssh), varids(ind_ssh)                            , &
         new_file(ind_ssh), ind_file(ind_ssh)                       , &
         ind_time(ind_ssh), ind_time_size(ind_ssh)                  , &
         dimsorders(:,ind_ssh)                                , &
         c_dir_ssh, c_prefix_ssh, maxsize_ssh, c_suffix_ssh , &
         nc_var_ssh, c_none, nc_att_mask_ssh,ssh_vvl(:,:,:,1:1))

         do kk = 1, dims_reg(3,3)
            e3t(:, :, kk, 1) = e3t0(:, :, kk, 1) * (1 + ssh_vvl(:, :, 1, 1)/totaldepth(:, :,1,1))
         enddo
    endif

    !======================= ZONAL TRANSPORT ==================!
    IF (key_interp_temporal) THEN
      uu(:,:,:,2) = uu(:,:,:,1) ! We are in temporal interpolation case
    ENDIF
    !-- Read Zonal Current --!
    CALL sub_input_data_seq_main(                         &
         ncids(1), varids(1)                            , &
         new_file(1), ind_file(1)                       , &
         ind_time(1), ind_time_size(1)                  , &
         dimsorders(:,1)                                , &
         c_dir_zo, c_prefix_zo, maxsize_zo, c_suffix_zo , &
         nc_var_zo, nc_var_eivu, nc_att_mask_zo,uu(:,:,:,1:1))

    !-- Compute transport --!
    IF (key_partialsteps) THEN

      DO ll = 1, 1
        DO kk = 1, dims_reg(3,3)
          ! partial steps case e3t is a 3D array
          uu(:,:,kk,ll) = uu(:,:,kk,ll) * e2u(:,:,1,1) * &
               MIN(e3t(:,:,kk,1), CSHIFT(e3t(:,:,kk,1),shift=1,dim=1))
        ENDDO
      ENDDO

    ELSE

      DO ll = 1, 1
        DO kk = 1, dims_reg(3,3)
          ! e3t is a vector (no partial steps)
          uu(:,:,kk,ll) = uu(:,:,kk,ll) * e2u(:,:,1,1) * e3t(1,1,kk,1)
        ENDDO
      ENDDO
    ENDIF

    IF (id_comments) THEN
      WRITE(lun_standard,*)' - Transport U: max ', MAXVAL(uu), ' min ', MINVAL(uu)
      ! lets look at the transport at a specific location RB
      WRITE(lun_standard,*)' - Transport: uu(625+1,250,16,1)', uu(95+1,52,16,1) !RB 2022/07/19
      !WRITE(lun_standard,*)' - shaaaape ',  SHAPE(uu) !shape suggests that uu is on the REDUCED grid
    ENDIF

    !============= MERIDIONAL TRANSPORT =======================!
    IF (key_interp_temporal) THEN
      vv(:,:,:,2) = vv(:,:,:,1) ! We are in temporal interpolation case
    ENDIF
    !-- Read Meridional Current --!
    CALL sub_input_data_seq_main(                         &
         ncids(2), varids(2)                            , &
         new_file(2), ind_file(2)                       , &
         ind_time(2), ind_time_size(2)                  , &
         dimsorders(:,2)                                , &
         c_dir_me,c_prefix_me, maxsize_me, c_suffix_me  , &
         nc_var_me,nc_var_eivv,nc_att_mask_me,vv(:,:,:,1:1))

    !-- Compute transport --
    IF (key_partialsteps) THEN

      DO ll = 1, 1
        DO kk = 1, dims_reg(3,3)
          ! partial steps case e3t is a 3D array
          vv(:,:,kk,ll) = vv(:,:,kk,ll) * e1v(:,:,1,1) * &
               MIN(e3t(:,:,kk,1), CSHIFT(e3t(:,:,kk,1),shift=1,dim=2))
        ENDDO
      ENDDO

    ELSE

      DO ll = 1, 1
        DO kk = 1, dims_reg(3,3)
          ! e3t is a vector (no partial steps)
          vv(:,:,kk,ll) = vv(:,:,kk,ll) * e1v(:,:,1,1) * e3t(1,1,kk,1)
        ENDDO
      ENDDO

    ENDIF

    IF (id_comments) THEN
      WRITE(lun_standard,*)' - Transport V: max ', MAXVAL(vv), ' min ', MINVAL(vv)
      ! lets look at the transport at a specific location RB (reduced grid)
      WRITE(lun_standard,*)' - Transport: vv(625+1,250,16,1)', vv(95+1,52,16,1) !RB 2022/07/19 
      !WRITE(lun_standard,*)' - Transport V problem: (642,274,16,ll) ', vv(3,32,16,1)
    ENDIF

    !=============== VERTICAL TRANSPORT =======================!
    IF (key_interp_temporal) THEN
      ww(:,:,:,2) = ww(:,:,:,1) ! We are in temporal interpolation case
    ENDIF

    WRITE(lun_standard,*)''

    IF (key_computew) THEN
      WRITE(lun_standard,*)'  - Vertical Transport is computing'
      !-- compute w transport from U and V transport --!
      CALL sub_w_comput(uu,vv,ww)

    ELSE

      !-- Read Vertical Current --!
      CALL sub_input_data_seq_main(                       &
           ncids(3), varids(3)                            , &
           new_file(3), ind_file(3)                       , &
           ind_time(3), ind_time_size(3)                  , &
           dimsorders(:,3)                                , &
           c_dir_ve,c_prefix_ve, maxsize_ve, c_suffix_ve  , &
           nc_var_ve,nc_var_eivw,nc_att_mask_ve,ww(:,:,:,1:1))

      !-- Compute transport --
      DO ll = 1, 1
        DO kk = 1, dims_reg(3,3)
          ww(:,:,kk,ll) = ww(:,:,kk,ll) * e1t(:,:,1,1) * e2t(:,:,1,1)
        ENDDO
      ENDDO

    ENDIF

    !!NG: 24/07/2009 
    IF (TRIM(w_surf_option) == 'E-P-R') THEN
      !!NG: 24/07/2009 We force here W="E-P-R" at the surface.
      !!NG: 24/07/2009 Physically is the better solution,
      !!NG: 24/07/2009 but in this case the non-divergence criteria
      !!NG: 24/07/2009 is not at all respected at the sea surface. 
      !!NG: 24/07/2009 This option has to be activated
      !!NG: 24/07/2009 if you know what you are doing.
      IF (key_alltracers) THEN
        ind_ep=7
      ELSE
        ind_ep=4
      ENDIF

      !-- Read E-R-R --!
      CALL sub_input_data_seq_main(                          &
           ncids(ind_ep), varids(ind_ep)                   , &
           new_file(ind_ep), ind_file(ind_ep)              , &
           ind_time(ind_ep), ind_time_size(ind_ep)         , &
           dimsorders(:,4)                                 , &
           c_dir_ep,c_prefix_ep, maxsize_ep, c_suffix_ep   , &
           nc_var_ep,c_none,nc_att_mask_ep,epr(:,:,:,1:1))

      ww(:,:,1,1) = epr(:,:,1,1) * epr_coef * e1t(:,:,1,1) * e2t(:,:,1,1)

    ELSEIF (TRIM(w_surf_option) == 'zero') THEN
      !!NG: 24/07/2009 We force here W=0 at the surface.
      !!NG: 24/07/2009 Is the same case that the precedent, but you don't
      !!NG: 24/07/2009 have the "E-P-R" fields available.
      !!NG: 24/07/2009 This option has to be activated
      !!NG: 24/07/2009 if you know what you are doing.
      ww(:,:,1,1) =0._rprec
    ELSE
      !!NG: 24/07/2009 Nothing is changed compared to the precedent version
      !!NG: 24/07/2009 of Ariane. The non-divergence criteria is respected,
      !!NG: 24/07/2009 but a lot of particles can be intercepted at the surface.
      !!NG: 24/07/2009 To reduce the number of particles intercepted at the sea
      !!NG: 24/07/2009 surface, please activate one the two options above mentioned.
    ENDIF

    IF (id_comments) THEN
      WRITE(lun_standard,*)' - Vertical Transport: max ', MAXVAL(ww), ' min ', MINVAL(ww)
      ! lets look at the transport at a specific location RB Mar 1, 2022
      WRITE(lun_standard,*)' - Transport: ww(625+1,250,16,1)', ww(95+1,52,16,1) !RB 2022/07/19
      !WRITE(lun_standard,*)' - Transport W problem: (642,274,16,ll) ', ww(3,32,16,1)
    ENDIF

    IF (key_alltracers) THEN

      !==================== TEMPERATURE =========================!
      IF (key_interp_temporal) THEN
        tt(:,:,:,2) = tt(:,:,:,1) ! We are in temporal interpolation case
      ENDIF
      !-- Read Temperature --!
      CALL sub_input_data_seq_main(                         &
           ncids(4), varids(4)                            , &
           new_file(4), ind_file(4)                       , &
           ind_time(4), ind_time_size(4)                  , &
           dimsorders(:,4)                                , &
           c_dir_te, c_prefix_te, maxsize_te, c_suffix_te , &
           nc_var_te,c_none,nc_att_mask_te,tt(:,:,:,1:1))

      !====================== SALINITY ===========================!
      IF (key_interp_temporal) THEN
        ss(:,:,:,2) = ss(:,:,:,1) ! We are in temporal interpolation case
      ENDIF
      !-- Read Salinity --!
      CALL sub_input_data_seq_main(                       &
           ncids(5), varids(5)                            , &
           new_file(5), ind_file(5)                       , &
           ind_time(5), ind_time_size(5)                  , &
           dimsorders(:,5)                                , &
           c_dir_sa,c_prefix_sa, maxsize_sa, c_suffix_sa  , &
           nc_var_sa,c_none,nc_att_mask_sa,ss(:,:,:,1:1))

      !======================= DENSITY ===========================!
      IF (key_interp_temporal) THEN
        rr(:,:,:,2) = rr(:,:,:,1) ! We are in temporal interpolation case
      ENDIF
      IF (key_computesigma) THEN

        WRITE(lun_standard,*)''
        WRITE(lun_standard,*)'------ Density is computing -------'
        !-- Compute Density from Temperature and Salinity --!
        DO l = 1, 1
          DO k = 1,  dims_reg(3,3)
            DO j = 1,  dims_reg(2,3)
              DO i = 1,  dims_reg(1,3)
                rr(i,j,k,l)=rhostp(1000._rprec,ss(i,j,k,l),tt(i,j,k,l),-ABS(zsigma))
              END DO
            END DO

          END DO
        ENDDO

        !- Comments -!
        IF (id_comments) THEN
          WRITE(lun_standard,*)' - Density: max ', MAXVAL(rr), &
               ' min ', MINVAL(rr)
        ENDIF

      ELSE

        !-- Read Density --!
        CALL sub_input_data_seq_main(                     &
             ncids(6), varids(6)                            , &
             new_file(6), ind_file(6)                       , &
             ind_time(6), ind_time_size(6)                  , &
             dimsorders(:,6)                                , &
             c_dir_de,c_prefix_de,maxsize_de, c_suffix_de   , &
             nc_var_de,c_none,nc_att_mask_de,rr(:,:,:,1:1))

      ENDIF

    ENDIF

  END SUBROUTINE sub_input_data_seq_opa
  !!***
  !=========================================================================
  !!****f* mod_input_data/sub_input_data_B2C_gridz()
  !! NAME
  !!   sub_input_data_B2C_gridz()
  !!
  !! FUNCTION
  !!   * read B-grid Model tracer variables. Depth levels are on Z (Not Sigma)
  !!   * loop on current and tracer variables: U, V, [W, T, S, R].
  !!   * call read netcdf.
  !!   * compute transport from currents.
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
  !!   * no arguments
  !!
  !! TODO
  !!   
  !!
  !! USED BY
  !!   * trajec
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_input_data_B2C_gridz()
    
    ! Variable number: U, V, [W, T, S, R]
    CHARACTER(len=4), PARAMETER :: c_none='NONE'

    INTEGER(kind  = iprec) :: kk, ll
    INTEGER(kind  = iprec) :: i,j,k,l

    INTEGER(kind = iprec) :: nb_i, nb_j, nb_k

    !-------------!
    ! Code begins !
    !-------------!
    nb_i=SIZE(uu,dim=1)
    nb_j=SIZE(uu,dim=2)
    nb_k=SIZE(uu,dim=3)

    !======================= ZONAL TRANSPORT ==================!
    !-- Read Zonal Current on B-grid --!
    ALLOCATE(uBgridmask(nb_i,nb_j,nb_k,1))
    CALL sub_memory(SIZE(uBgridmask),'i','uBgridmask','sub_input_data_B2C_gridz')

    CALL sub_input_data_main(c_dir_zo,c_prefix_zo,ind0_zo,indn_zo,maxsize_zo, &
         c_suffix_zo, nc_var_zo,nc_var_eivu,nc_att_mask_zo,uu(:,:,:,:), &
         uBgridmask(:,:,:,1:1))

    !-- Interpolate Zonal Current from B-grid to C-grid --!
    !! CALL sub_B2C_grid_interpolation(uu, e2f, e2u, 'u', uBgridmask)
    CALL sub_B2C_grid_interpolation(uu, e2f, e2u, e3f,e3t,'u', key_partialsteps)

    !-- Save U interpolated from B grid to C grid --!
!!$    IF (key_B2C_save_data) THEN
!!$      CALL sub_B2C_grid_save_U(uu,uBgridmask)
!!$    ENDIF

    !-- Compute transport --!
    IF (key_partialsteps) THEN

      DO ll = 1, lmt
        DO kk = 1, dims_reg(3,3)
          ! partial steps case e3t is a 3D array
          uu(:,:,kk,ll) = uu(:,:,kk,ll) * e2u(:,:,1,1) * &
               MIN(e3t(:,:,kk,1), CSHIFT(e3t(:,:,kk,1),shift=1,dim=1))
        ENDDO
      ENDDO

    ELSE

      DO ll = 1, lmt
        DO kk = 1, dims_reg(3,3)
          ! e3t is a vector (no partial steps)
          uu(:,:,kk,ll) = uu(:,:,kk,ll) * e2u(:,:,1,1) * e3t(1,1,kk,1)
        ENDDO
      ENDDO
    ENDIF

    WRITE(lun_standard,*)' - Transport U: max ', MAXVAL(uu), ' min ', MINVAL(uu)

    CALL sub_memory(-SIZE(uBgridmask),'i','uBgridmask','sub_input_data_B2C_gridz')
    DEALLOCATE(uBgridmask)

    !============= MERIDIONAL TRANSPORT =======================!

    !-- Read Meridional Current --!
    ALLOCATE(vBgridmask(nb_i,nb_j,nb_k,1))
    CALL sub_memory(SIZE(vBgridmask),'i','vBgridmask','sub_input_data_B2C_gridz')

    CALL sub_input_data_main(c_dir_me,c_prefix_me,ind0_me,indn_me,maxsize_me, &
         c_suffix_me, nc_var_me,nc_var_eivv,nc_att_mask_me,vv(:,:,:,:), &
         vBgridmask(:,:,:,1:1))

    !-- Interpolate Zonal Current from B-grid to C-grid --!
    !! CALL sub_B2C_grid_interpolation(vv, e1f, e1v,'v', vBgridmask)
    CALL sub_B2C_grid_interpolation(vv, e1f, e1v, e3f,e3t,'v', key_partialsteps)

    !-- Save V interpolated from B grid to C grid --!
!!$    IF (key_B2C_save_data) THEN
!!$      CALL sub_B2C_grid_save_V(vv,vBgridmask)
!!$    ENDIF


    !-- Compute transport --
    IF (key_partialsteps) THEN

      DO ll = 1, lmt
        DO kk = 1, dims_reg(3,3)
          ! partial steps case e3t is a 3D array
          vv(:,:,kk,ll) = vv(:,:,kk,ll) * e1v(:,:,1,1) * &
               MIN(e3t(:,:,kk,1), CSHIFT(e3t(:,:,kk,1),shift=1,dim=2))
        ENDDO
      ENDDO

    ELSE

      DO ll = 1, lmt
        DO kk = 1, dims_reg(3,3)
          ! e3t is a vector (no partial steps)
          vv(:,:,kk,ll) = vv(:,:,kk,ll) * e1v(:,:,1,1) * e3t(1,1,kk,1)
        ENDDO
      ENDDO

    ENDIF

    WRITE(lun_standard,*)' - Transport V: max ', MAXVAL(vv), ' min ', MINVAL(vv)

    CALL sub_memory(-SIZE(vBgridmask),'i','vBgridmask','sub_input_data_B2C_gridz')
    DEALLOCATE(vBgridmask)

    !=============== VERTICAL TRANSPORT =======================!
    IF (key_read_w) THEN

      !-- Read Vertical Current --!
      CALL sub_input_data_main(c_dir_ve,c_prefix_ve,ind0_ve,indn_ve,maxsize_ve, &
           c_suffix_ve, nc_var_ve,nc_var_eivw,nc_att_mask_ve,ww(:,:,:,:))


      !-- Save V interpolated on from B grid to C grid --!
      IF (key_B2C_save_data) THEN
        CALL sub_B2C_grid_save_W(ww)
      ENDIF

      !-- Compute transport --
      DO ll = 1, lmt
        DO kk = 1, dims_reg(3,3)
          ww(:,:,kk,ll) = ww(:,:,kk,ll) * e1t(:,:,1,1) * e2t(:,:,1,1)
        ENDDO
      ENDDO

      WRITE(lun_standard,*)''
      WRITE(lun_standard,*)' - vertical current is read -'

    ELSE

      !-- compute w transport from U and V transport --!
      WRITE(lun_standard,*)''
      WRITE(lun_standard,*)' - vertical current is computed -'
      CALL sub_w_comput(uu,vv,ww)

    ENDIF


    WRITE(lun_standard,*)' - Transport W: max ', MAXVAL(ww), ' min ', MINVAL(ww)

    IF (key_alltracers) THEN

      !==================== TEMPERATURE =========================!
      !-- Read Temperature --!
      CALL sub_input_data_main(c_dir_te,c_prefix_te,ind0_te,indn_te,maxsize_te, &
           c_suffix_te, nc_var_te,c_none,nc_att_mask_te,tt(:,:,:,:))

      !====================== SALINITY ===========================!
      !-- Read Salinity --!
      CALL sub_input_data_main(c_dir_sa,c_prefix_sa,ind0_sa,indn_sa,maxsize_sa, &
           c_suffix_sa, nc_var_sa,c_none,nc_att_mask_sa,ss(:,:,:,:))

      !======================= DENSITY ===========================!
      IF (key_computesigma) THEN

        !-- Compute Density from Temperature and Salinity --!
        DO l = 1, lmt
          DO k = 1,  dims_reg(3,3)
            DO j = 1,  dims_reg(2,3)
              DO i = 1,  dims_reg(1,3)
                rr(i,j,k,l)=rhostp(1000._rprec,ss(i,j,k,l),tt(i,j,k,l),-ABS(zsigma))
              END DO
            END DO
          END DO
        ENDDO

      ELSE

        !-- Read Density --!
        CALL sub_input_data_main(c_dir_de,c_prefix_de,ind0_de,indn_de,maxsize_de, &
             c_suffix_de, nc_var_de,c_none,nc_att_mask_de,rr(:,:,:,:))

      ENDIF

    ENDIF

  END SUBROUTINE sub_input_data_B2C_gridz
  !!***
  !=========================================================================
  !!****f* mod_input_data_seq/sub_input_data_B2C_gridz_seq()
  !! NAME
  !!   ssub_input_data_B2C_gridz_seq()
  !!
  !! FUNCTION
  !!
  !! AUTHOR
  !!   * Origin  : Nicolas Grima (April 2009)
  !! 
  !! CREATION DATE
  !!   * April 2009
  !!
  !! HISTORY
  !!   Date (dd/mm/yyyy/) - Modification(s)
  !!
  !! ARGUMENTS
  !!   *
  !!
  !! TODO
  !!   
  !!
  !! USED BY
  !!   * posini_seq
  !!   * trajec_seq
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_input_data_B2C_gridz_seq( &
       ncids                             , &
       varids                            , &
       new_file                          , &
       ind_file                          , &
       ind_time                          , &
       ind_time_size                     , &
       dimsorders                        , &
       first_call                        , &
       lref                              )


    INTEGER(kind = iprec) , DIMENSION(:), INTENT(inout) :: ncids
    INTEGER(kind = iprec) , DIMENSION(:), INTENT(inout) :: varids
    LOGICAL               , DIMENSION(:), INTENT(inout) :: new_file
    INTEGER(kind = iprec) , DIMENSION(:), INTENT(inout) :: ind_file
    INTEGER(kind = iprec) , DIMENSION(:), INTENT(inout) :: ind_time
    INTEGER(kind = iprec) , DIMENSION(:), INTENT(inout) :: ind_time_size
    INTEGER(kind = iprec) , DIMENSION(:,:), INTENT(inout) :: dimsorders
    LOGICAL               , OPTIONAL    , INTENT(in)    :: first_call
    INTEGER(kind = iprec) , OPTIONAL      , INTENT(in)  :: lref

    ! Variable number: U, V, [W, T, S, R]
    CHARACTER(len=4), PARAMETER :: c_none='NONE'

    INTEGER(kind  = iprec) :: kk, ll
    INTEGER(kind  = iprec) :: i,j,k,l
    INTEGER(kind  = iprec) :: nb_i, nb_j, nb_k

    ! To save data interpolated on C grid
    LOGICAL                :: new_file_tmp
    LOGICAL                :: close_file_tmp
    INTEGER(kind  = iprec) :: ind_file_tmp
    INTEGER(kind  = iprec) :: dimt_tmp

    INTEGER(kind  = iprec) :: iter

    !-------------!
    ! Code begins !
    !-------------!
    IF (PRESENT(lref))THEN
      iter = lref
    ELSE
      iter = iZero
    ENDIF

    nb_i=SIZE(uu,dim=1)
    nb_j=SIZE(uu,dim=2)
    nb_k=SIZE(uu,dim=3)

    !======================= ZONAL TRANSPORT ==================!
    IF (key_interp_temporal) THEN
      uu(:,:,:,2) = uu(:,:,:,1) ! We are in temporal interpolation case
    ENDIF

    !-- Read Zonal Current on B-grid --!
    IF (.NOT.ALLOCATED(uBgridmask)) THEN
      ALLOCATE(uBgridmask(nb_i,nb_j,nb_k,1))
      CALL sub_memory(SIZE(uBgridmask),'i','uBgridmask','sub_input_data_B2C_gridz_seq')
    ENDIF

    new_file_tmp = new_file(1)

    !-- Read Zonal Current --!
    CALL sub_input_data_seq_main(                         &
         ncids(1), varids(1)                            , &
         new_file(1), ind_file(1)                       , &
         ind_time(1), ind_time_size(1)                  , &
         dimsorders(:,1)                                , &
         c_dir_zo, c_prefix_zo, maxsize_zo, c_suffix_zo , &
         nc_var_zo, nc_var_eivu, nc_att_mask_zo,uu(:,:,:,1:1), &
         uBgridmask(:,:,:,1:1))

    ind_file_tmp   = ind_file(1)
    close_file_tmp = new_file(1)
    dimt_tmp       = ind_time_size(1)

    IF (key_B2C_save_data) THEN
      CALL write_netcdf ( &
           uu(:,:,:,1)  , &
           'uuBgrid'    , &
           'current.nc' , &
           iter         )
    ENDIF

    !-- Interpolate Zonal Current from B-grid to C-grid --!
    !! CALL sub_B2C_grid_interpolation(uu, e2f, e2u, 'u', uBgridmask)
    CALL sub_B2C_grid_interpolation(uu, e2f, e2u, e3f,e3t,'u', key_partialsteps)

    IF (key_B2C_save_data) THEN
      CALL write_netcdf ( &
           uu(:,:,:,1)  , &
           'uuCgrid'    , &
           'current.nc' , &
           iter         )
    ENDIF

    !-- Compute transport --!
    IF (key_partialsteps) THEN

      DO ll = 1, 1
        DO kk = 1, dims_reg(3,3)
          ! partial steps case e3t is a 3D array
          uu(:,:,kk,ll) = uu(:,:,kk,ll) * e2u(:,:,1,1) * &
               MIN(e3t(:,:,kk,1), CSHIFT(e3t(:,:,kk,1),shift=1,dim=1))
        ENDDO
      ENDDO

    ELSE

      DO ll = 1, 1
        DO kk = 1, dims_reg(3,3)
          ! e3t is a vector (no partial steps)
          uu(:,:,kk,ll) = uu(:,:,kk,ll) * e2u(:,:,1,1) * e3t(1,1,kk,1)
        ENDDO
      ENDDO
    ENDIF

    IF (id_comments) THEN
      WRITE(lun_standard,*)' - Transport U: max ', MAXVAL(uu), ' min ', MINVAL(uu)
    ENDIF

    IF (key_B2C_save_data) THEN
      CALL write_netcdf  ( &
           uu(:,:,:,1)   , &
           'uuCgrid'     , &
           'transport.nc', &
           iter     )
    ENDIF

    !============= MERIDIONAL TRANSPORT =======================!
    IF (.NOT.ALLOCATED(vBgridmask)) THEN
      ALLOCATE(vBgridmask(nb_i,nb_j,nb_k,1))
      CALL sub_memory(SIZE(vBgridmask),'i','vBgridmask','sub_input_data_B2C_gridz_seq')
    ENDIF

    IF (key_interp_temporal) THEN
      vv(:,:,:,2) = vv(:,:,:,1) ! We are in temporal interpolation case
    ENDIF

    new_file_tmp = new_file(2)

    !-- Read Meridional Current --!
    CALL sub_input_data_seq_main(                         &
         ncids(2), varids(2)                            , &
         new_file(2), ind_file(2)                       , &
         ind_time(2), ind_time_size(2)                  , &
         dimsorders(:,2)                                , &
         c_dir_me,c_prefix_me, maxsize_me, c_suffix_me  , &
         nc_var_me,nc_var_eivv,nc_att_mask_me,vv(:,:,:,1:1), &
         vBgridmask(:,:,:,1:1))

    ind_file_tmp   = ind_file(2)
    close_file_tmp = new_file(2)
    dimt_tmp       = ind_time_size(2)

    IF (key_B2C_save_data) THEN
      CALL write_netcdf ( &
           vv(:,:,:,1)  , &
           'vvBgrid'    , &
           'current.nc' , &
           iter         )
    ENDIF

    !-- Interpolate Zonal Current from B-grid to C-grid --!
    !! CALL sub_B2C_grid_interpolation(vv, e1f, e1v,'v', vBgridmask)
    CALL sub_B2C_grid_interpolation(vv, e1f, e1v, e3f,e3t,'v', key_partialsteps)

    IF (key_B2C_save_data) THEN
      CALL write_netcdf ( &
           vv(:,:,:,1)  , &
           'vvCgrid'    , &
           'current.nc' , &
           iter         )
    ENDIF

    !-- Save V interpolated from B grid to C grid --!
!!$    IF (key_B2C_save_data) THEN
!!$      IF (PRESENT(first_call)) THEN
!!$        IF (first_call) THEN
!!$          CALL sub_B2C_grid_save_V_seq(&
!!$               vv                    , &
!!$               new_file_tmp          , &
!!$               ind_file_tmp          , &
!!$               maxsize_zo            , &
!!$               dimt_tmp              , &
!!$               close_file_tmp          &
!!$               )
!!$        ENDIF
!!$      ENDIF
!!$    ENDIF

    !-- Compute transport --
    IF (key_partialsteps) THEN

      DO ll = 1, 1
        DO kk = 1, dims_reg(3,3)
          ! partial steps case e3t is a 3D array
          vv(:,:,kk,ll) = vv(:,:,kk,ll) * e1v(:,:,1,1) * &
               MIN(e3t(:,:,kk,1), CSHIFT(e3t(:,:,kk,1),shift=1,dim=2))
        ENDDO
      ENDDO

    ELSE

      DO ll = 1, 1
        DO kk = 1, dims_reg(3,3)
          ! e3t is a vector (no partial steps)
          vv(:,:,kk,ll) = vv(:,:,kk,ll) * e1v(:,:,1,1) * e3t(1,1,kk,1)
        ENDDO
      ENDDO

    ENDIF

    IF (id_comments) THEN
      WRITE(lun_standard,*)' - Transport V: max ', MAXVAL(vv), ' min ', MINVAL(vv)
    ENDIF

    IF (key_B2C_save_data) THEN
      CALL write_netcdf  ( &
           vv(:,:,:,1)   , &
           'vvCgrid'     , &
           'transport.nc', &
           iter     )

    ENDIF

    !=============== VERTICAL TRANSPORT =======================!
    IF (key_interp_temporal) THEN
      ww(:,:,:,2) = ww(:,:,:,1) ! We are in temporal interpolation case
    ENDIF

    IF (key_read_w) THEN

      new_file_tmp = new_file(3)

      !-- Read Vertical Current --!
      CALL sub_input_data_seq_main(                       &
           ncids(3), varids(3)                            , &
           new_file(3), ind_file(3)                       , &
           ind_time(3), ind_time_size(3)                  , &
           dimsorders(:,3)                                , &
           c_dir_ve,c_prefix_ve, maxsize_ve, c_suffix_ve  , &
           nc_var_ve,nc_var_eivw,nc_att_mask_ve,ww(:,:,:,1:1))

      ind_file_tmp   = ind_file(3)
      close_file_tmp = new_file(3)
      dimt_tmp       = ind_time_size(3)

      !-- Save W interpolated from B grid to C grid --!
      IF (key_B2C_save_data) THEN
        IF (PRESENT(first_call)) THEN
          IF (first_call) THEN
            CALL sub_B2C_grid_save_W_seq(&
                 ww                    , &
                 new_file_tmp          , &
                 ind_file_tmp          , &
                 maxsize_zo            , &
                 dimt_tmp              , &
                 close_file_tmp          &
                 )
          ENDIF
        ENDIF
      ENDIF

      !-- Compute transport --
      DO ll = 1, 1
        DO kk = 1, dims_reg(3,3)
          ww(:,:,kk,ll) = ww(:,:,kk,ll) * e1t(:,:,1,1) * e2t(:,:,1,1)
        ENDDO
      ENDDO

    ELSE

      !-- compute w transport from U and V transport --!
      WRITE(lun_standard,*)''
      WRITE(lun_standard,*)' - vertical current is computed -'
      CALL sub_w_comput(uu,vv,ww)

    ENDIF

    IF (id_comments) THEN
      WRITE(lun_standard,*)' - Transport W: max ', MAXVAL(ww), ' min ', MINVAL(ww)
    ENDIF

    IF (key_alltracers) THEN

      !==================== TEMPERATURE =========================!
      IF (key_interp_temporal) THEN
        tt(:,:,:,2) = tt(:,:,:,1) ! We are in temporal interpolation case
      ENDIF
      !-- Read Temperature --!
      CALL sub_input_data_seq_main(                       &
           ncids(4), varids(4)                            , &
           new_file(4), ind_file(4)                       , &
           ind_time(4), ind_time_size(4)                  , &
           dimsorders(:,3)                                , &
           c_dir_te, c_prefix_te, maxsize_te, c_suffix_te , &
           nc_var_te,c_none,nc_att_mask_te,tt(:,:,:,1:1))

      !====================== SALINITY ===========================!
      IF (key_interp_temporal) THEN
        ss(:,:,:,2) = ss(:,:,:,1) ! We are in temporal interpolation case
      ENDIF
      !-- Read Salinity --!
      CALL sub_input_data_seq_main(                       &
           ncids(5), varids(5)                            , &
           new_file(5), ind_file(5)                       , &
           ind_time(5), ind_time_size(5)                  , &
           dimsorders(:,5)                                , &
           c_dir_sa,c_prefix_sa, maxsize_sa, c_suffix_sa  , &
           nc_var_sa,c_none,nc_att_mask_sa,ss(:,:,:,1:1))

      !======================= DENSITY ===========================!
      IF (key_interp_temporal) THEN
        rr(:,:,:,2) = rr(:,:,:,1) ! We are in temporal interpolation case
      ENDIF
      IF (key_computesigma) THEN

        !-- Compute Density from Temperature and Salinity --!
        DO l = 1, 1
          DO k = 1,  dims_reg(3,3)
            DO j = 1,  dims_reg(2,3)
              DO i = 1,  dims_reg(1,3)
                rr(i,j,k,l)=rhostp(1000._rprec,ss(i,j,k,l),tt(i,j,k,l),-ABS(zsigma))
              END DO
            END DO

          END DO
        ENDDO

      ELSE

        !-- Read Density --!
        CALL sub_input_data_seq_main(                     &
             ncids(6), varids(6)                            , &
             new_file(6), ind_file(6)                       , &
             ind_time(6), ind_time_size(6)                  , &
             dimsorders(:,6)                                , &
             c_dir_de,c_prefix_de,maxsize_de, c_suffix_de   , &
             nc_var_de,c_none,nc_att_mask_de,rr(:,:,:,1:1))

      ENDIF

    ENDIF

  END SUBROUTINE sub_input_data_B2C_gridz_seq
  !!***

  !=========================================================================
  !!****f* mod_input_data/sub_input_data_roms()
  !! NAME
  !!   sub_input_data_roms()
  !!
  !! FUNCTION
  !!   Data used in Ariane are on C grid and they are based on 
  !!   OPA conventions. ROMS data are on a C grid but don't respect
  !!   all OPA conventions.
  !!   The vertical axis is reversed.
  !!   U is defined only on (imt-1, jmt  , kmt-1)
  !!   V is defined only on (imt  , jmt-1, kmt-1)
  !!   W must be computed
  !!   T and S are defined on (imt, jmt, kmt-1)
  !!
  !!   * read ROMS tracer variables.
  !!   * loop on current and tracer variables: U, V, [W, T, S, R].
  !!   * call read netcdf.
  !!   * compute scale factor e3t.
  !!   * compute transport from currents.
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
  !!   * no arguments
  !!
  !! TODO
  !!   
  !!
  !! USED BY
  !!   * trajec
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_input_data_roms()

    INTEGER(kind=iprec) :: i,j,k,l
    CHARACTER(len=4), PARAMETER :: c_none='NONE'

    !!NG: 05/08/2008 REAL(kind=rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: &
    !!NG: 05/08/2008      zeta, &
    !!NG: 05/08/2008      zw0 , &
    !!NG: 05/08/2008      e3u , &
    !!NG: 05/08/2008      e3v

    INTEGER(kind=iprec) :: ncid
    INTEGER(kind=iprec) :: varid
    INTEGER(kind = ishort) :: is_err !cjmp!

    LOGICAL :: key_ijt

    REAL(kind = rprec) :: hc ! S-coordinate parameter, critical depth (m)
    REAL(kind = rprec), DIMENSION(:), ALLOCATABLE :: &
         sc_w , & ! S-coordinate at W-points
         cs_w     ! S-coordinate stretching curves at W-points
    CHARACTER(len = 128)  :: c_filename   ! netcdf file name             !cjmp!

    CHARACTER(len=11 ), PARAMETER :: nc_glbatt_Vstretching='Vstretching' !cjmp!
    CHARACTER(len=10 ), PARAMETER :: nc_glbatt_Vtransform='Vtransform'   !cjmp!
    INTEGER(kind = iprec) :: Vstretching                                 !cjmp!
    INTEGER(kind = iprec) :: Vtransform                                  !cjmp!



    !======================= ZONAL CURRENT ==================!
    !- Read U on ROMS grid -!
    CALL sub_input_data_main(c_dir_zo,c_prefix_zo, &
         ind0_zo,indn_zo,maxsize_zo, &
         c_suffix_zo, nc_var_zo,nc_var_eivu,nc_att_mask_zo, &
         uu(1:dims_reg(1,3)-1,:,dims_reg(3,3)-1:1:-1,:))

    !============= MERIDIONAL CURRENT =======================!
    !- Read V on ROMS grid -!
    CALL sub_input_data_main(c_dir_me,c_prefix_me,&
         ind0_me,indn_me,maxsize_me, &
         c_suffix_me, nc_var_me,nc_var_eivv,nc_att_mask_me,&
         vv(:,1:dims_reg(2,3)-1,dims_reg(3,3)-1:1:-1,:))

    IF (key_alltracers) THEN
      !==================== TEMPERATURE =========================!
      !- Read T on ROMS grid -!
      CALL sub_input_data_main(c_dir_te,c_prefix_te,&
           ind0_te,indn_te,maxsize_te, &
           c_suffix_te, nc_var_te,c_none,nc_att_mask_te,&
           tt(:,:,dims_reg(3,3)-1:1:-1,:))

      !====================== SALINITY ===========================!
      !- Read S on ROMS grid -!
      CALL sub_input_data_main(c_dir_sa,c_prefix_sa, &
           ind0_sa,indn_sa,maxsize_sa, &
           c_suffix_sa, nc_var_sa,c_none,nc_att_mask_sa, &
           ss(:,:,dims_reg(3,3)-1:1:-1,:))

      !====================== COMPUTE RR =========================!
      WRITE(lun_standard,*)''
      WRITE(lun_standard,*)'------ Density is computing -------'
      !-- Compute Density from Temperature and Salinity --!
      DO l = 1, lmt
        DO k = 1,  dims_reg(3,3)
          DO j = 1,  dims_reg(2,3)
            DO i = 1,  dims_reg(1,3)
              rr(i,j,k,l)=rhostp(1000._rprec,ss(i,j,k,l),tt(i,j,k,l),-ABS(zsigma))
            END DO
          END DO
        END DO
      ENDDO
      !- Comments -!
      WRITE(lun_standard,*)' - Density: max ', MAXVAL(rr), ' min ', MINVAL(rr)

    ENDIF

    !====================== ZETA 2D + Time =======================!
    !- Read Zeta on ROMS grid -!
    ALLOCATE(zeta_a(dims_reg(1,3),dims_reg(2,3), 1, dims_reg(4,3)))
    CALL sub_memory(SIZE(zeta_a),'r','zeta_a','sub_input_data_roms')

    key_ijt = .TRUE.

    CALL sub_input_data_main(c_dir_ze,c_prefix_ze, &
         ind0_ze,indn_ze,maxsize_ze, &
         c_suffix_ze, nc_var_ze,c_none,nc_att_mask_ze, &
         zeta_a(:,:,:,:), key_ijt = key_ijt)

    !========================= ZW0 ==============================!
    !- Read global attributes hc (scalar) sc_w and cs_w (vectors) -!
    !- in AVG files.
    !- Build file name from temperature file -!
    CALL sub_build_filename(ind0_ze, maxsize_ze, &
         c_prefix_ze, c_suffix_ze, c_filename)

    !======== Read ROMS global attributes hc, sc_w and Cs_w ======!
    ALLOCATE(sc_w(dims_reg(3,3)))
    CALL sub_memory(SIZE(sc_w),'r','sc_w','sub_input_data_roms')

    ALLOCATE(cs_w(dims_reg(3,3)))
    CALL sub_memory(SIZE(cs_w),'r','cs_w','sub_input_data_roms')

    !! NG : 29/09/2009
    IF (roms_global_attribute) THEN 
      CALL sub_read_netcdf_global_att(                   &
           dir_glbatt, fn_glbatt                       , &
           nc_glbatt_hc, nc_glbatt_sc_w, nc_glbatt_Cs_w, &
           hc, sc_w, Cs_w)
    ELSE

      !- Open NetCDF file -!
      CALL sub_open_netcdf_file(dir_glbatt, fn_glbatt, ncid)

      !- Read hc -!
      CALL sub_read_netcdf_varid_ndims(ncid, nc_glbatt_hc, varid)
      CALL sub_read_netcdf_var( ncid, varid, hc)

      !- Read sc_w -!

      !cjmp modified to calculate sc_w (called s_w in ROMS 3)
      !analytically, since sc_w does not vary by coordinate
      !transform.  NOTE WELL, THE VERTICAL INDEX IS REVERSED BETWEEN
      !ROMS AND ARIANE

      !cjmp_old CALL sub_read_netcdf_varid_ndims(ncid, nc_glbatt_sc_w, varid)
      !cjmp_old CALL sub_read_netcdf_var1d( ncid, varid, sc_w(dims_reg(3,3):1:-1))

      DO k=0,dims_reg(3,3)-1
        sc_w(dims_reg(3,3)-k)=(1.0*k-(dims_reg(3,3)-1))/(dims_reg(3,3)-1)
      ENDDO
      !cjmp end modification

      !- Read Cs_w -!
      CALL sub_read_netcdf_varid_ndims(ncid, nc_glbatt_Cs_w, varid)
      CALL sub_read_netcdf_var1d( ncid, varid, Cs_w(dims_reg(3,3):1:-1))

      !- Close netcdf file -!
      CALL sub_close_netcdf_file(ncid)

    ENDIF


    !cjmp modified to attempt to read Vtransform and Vstretching from ROMS !
    !cjmp history files.  If they are not present in file, assume Vtransform=1  !
    !cjmp and Vstretching=1 ; Note, this requires checking if the variable exists !
    !cjmp before trying to read it, and also reading it directly, for the routines
    !cjmp in mod_netcdf assume the variable has type real. !
    !- Open NetCDF file -!
    !!NG: 16 nov 2018 add Vtransform 2 for Pierre Amael
    !! vtransform=1: z_unperturbated (s) = hc * (sc(s) - Cs(s)) + Cs(s) * h
    !! vtransform=2: z_unperturbated (s) = h * (sc(s) * hc + Cs(s) * h) / (h + hc)
    !!NG: 16 nov 2018
    CALL sub_open_netcdf_file(dir_glbatt, fn_glbatt, ncid)

    !now read Vtransform
    !first, check if it exists
    is_err     = nf90_inq_varid( &
         ncid  = ncid          , &
         name  = TRIM(nc_glbatt_Vtransform)  , &
         varid = varid           &
         )
    !- Test if Vtransform exists -!
    IF (is_err /= nf90_noerr) THEN
      WRITE(lun_standard,*)' Vtransform does not exits, setting to 1 '
      Vtransform=1
    ELSE
      is_err=nf90_get_var(ncid=ncid,varid=varid,values=Vtransform)
      IF (is_err /= nf90_noerr) THEN
        WRITE(lun_error,*) 'Error reading Vtransform'
      ELSE
        WRITE(lun_standard,*) 'Vtransform is ',Vtransform
      ENDIF
    ENDIF

    !now read Vstretching
    !first, check if it exists
    is_err     = nf90_inq_varid( &
         ncid  = ncid          , &
         name  = TRIM(nc_glbatt_Vstretching)  , &
         varid = varid           &
         )
    !- Test if Vstretching exists -!
    IF (is_err /= nf90_noerr) THEN
      WRITE(lun_standard,*)' Vstretching does not exits, setting to 1 '
      Vstretching=1
    ELSE
      is_err=nf90_get_var(ncid=ncid,varid=varid,values=Vstretching)
      IF (is_err /= nf90_noerr) THEN
        WRITE(lun_error,*) 'Error reading Vstretching'
      ELSE
        WRITE(lun_standard,*) 'Vstretching is ',Vstretching
      ENDIF
    ENDIF

    !- Close netcdf file -!
    CALL sub_close_netcdf_file(ncid)

    !cjmp check if Vstretching or Vtransform .neq. 1; if so, stop
    IF (Vstretching /= 1) THEN
      WRITE(lun_error,*) 'Currently, ARIANE only supports Vstretching=1'
      STOP
    ENDIF
    IF ((Vtransform < 1).or.(Vtransform > 2)) THEN
      WRITE(lun_error,*) 'Currently, ARIANE only supports Vtransform=1 or 2'
      STOP
    ENDIF
    !cjmp end modification

    !============== ZW0 (3D - constant in time) =================!
    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' ------ ZW0 is computing ------'
    WRITE(lun_standard,*)' - hc: '  , hc
    WRITE(lun_standard,*)' - SC_W: ', sc_w
    WRITE(lun_standard,*)' - CS_w: ', cs_w
    ALLOCATE(zw0(dims_reg(1,3), dims_reg(2,3), dims_reg(3,3), 1))
    CALL sub_memory(SIZE(zw0),'r','zw0','sub_input_data_roms')

    !!NG: 16 nov 2018 add Vtransform 2 for Pierre Amael
    !! vtransform=1: z_unperturbated (s) = hc * (sc(s) - Cs(s)) + Cs(s) * h
    !! vtransform=2: z_unperturbated (s) = h * (sc(s) * hc + Cs(s) * h) / (h + hc)
    !!NG: 16 nov 2018
    IF (Vtransform ==1) THEN
      DO k = 1, dims_reg(3,3)
        DO j = 1, dims_reg(2,3)
          DO i = 1, dims_reg(1,3)
            zw0(i,j,k,1) = hc * sc_w(k) + &
                 (h_roms(i,j,1,1) - hc) * cs_w(k)
          ENDDO
        ENDDO
      ENDDO
    ELSEIF (Vtransform ==2) THEN
      DO k = 1, dims_reg(3,3)
        DO j = 1, dims_reg(2,3)
          DO i = 1, dims_reg(1,3)
            zw0(i,j,k,1) = h_roms(i,j,1,1) * &
                 ( sc_w(k) * hc + cs_w(k) * h_roms(i,j,1,1)) / &
                 (h_roms(i,j,1,1) + hc)
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    WRITE(lun_standard,*)' - ZW0 : max ', MAXVAL(zw0), ' min ', MINVAL(zw0)

    CALL sub_memory(-SIZE(sc_w),'r','sc_w','sub_input_data_roms')
    DEALLOCATE(sc_w)

    CALL sub_memory(-SIZE(cs_w),'r','cs_w','sub_input_data_roms')
    DEALLOCATE(cs_w)


    !==================== ZZ_WW (3D + time) ========================!
    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' ------ ZZ_WW is computing ------'

    ALLOCATE(zz_ww(dims_reg(1,3), dims_reg(2,3), dims_reg(3,3), dims_reg(4,3)))
    CALL sub_memory(SIZE(zz_ww),'r','zz_ww','sub_input_data_roms')

    DO l = 1, dims_reg(4,3)
      DO k = 1, dims_reg(3,3)
        DO j = 1, dims_reg(2,3) 
          DO i = 1, dims_reg(1,3)
            zz_ww(i,j,k,l) = &
                 zeta_a(i,j,1,l) + zw0(i,j,k,1) * ( 1._rprec + zeta_a(i,j,1,l) / h_roms(i,j,1,1))
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    WRITE(lun_standard,*)' - ZZ_WW : max ', MAXVAL(zz_ww), ' min ', MINVAL(zz_ww)

    CALL sub_memory(-SIZE(zw0),'r','zw0','sub_input_data_roms')
    DEALLOCATE(zw0)

    !=================== E3T (3D + time) ========================!
    !- Compute E3T from ZZ_WW -!
    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' ------ E3T scale factor is computing ------'
    ALLOCATE(e3t(dims_reg(1,3), dims_reg(2,3), dims_reg(3,3), dims_reg(4,3)))
    CALL sub_memory(SIZE(e3t),'r','e3t','sub_input_data_roms')

    e3t(:,:,1:dims_reg(3,3)-1,:)= &
         zz_ww(:,:,1:dims_reg(3,3)-1,:) - zz_ww(:,:,2:dims_reg(3,3),:)

    e3t(:,:,dims_reg(3,3),:) = e3t(:,:,dims_reg(3,3)-1,:)
    WRITE(lun_standard,*)' - E3T : max ', MAXVAL(e3t), ' min ', MINVAL(e3t)

    !========================= E3U ==============================!
    ALLOCATE(e3u(dims_reg(1,3), dims_reg(2,3), dims_reg(3,3), dims_reg(4,3)))
    CALL sub_memory(SIZE(e3u),'r','e3u','sub_input_data_roms')

    e3u(1:dims_reg(1,3)-1,:,:,:) = &
         .5_rprec *( e3t(1:dims_reg(1,3)-1,:,:,:) + e3t(2:dims_reg(1,3),:,:,:) )

    e3u(dims_reg(1,3),:,:,:) = e3u(dims_reg(1,3)-1,:,:,:)

    !========================= E3V ==============================!
    ALLOCATE(e3v(dims_reg(1,3), dims_reg(2,3), dims_reg(3,3), dims_reg(4,3)))
    CALL sub_memory(SIZE(e3v),'r','e3v','sub_input_data_roms')

    e3v(:,1:dims_reg(2,3)-1,:,:) = &
         .5_rprec *( e3t(:,1:dims_reg(2,3)-1,:,:) + e3t(:,2:dims_reg(2,3),:,:) )

    e3v(:,dims_reg(2,3),:,:) = e3v(:,dims_reg(2,3)-1,:,:)

    !====================== TRANSPORT ===========================!
    !=== U ===!
    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' ------ Zonal transport is computing ------'
    DO l = 1, dims_reg(4,3)
      DO k = 1, dims_reg(3,3)
        DO j = 1, dims_reg(2,3)
          DO i = 1, dims_reg(1,3)
            uu(i,j,k,l) =  uu(i,j,k,l) * e3u(i,j,k,l) * e2u(i,j,1,1)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    WRITE(lun_standard,*)' - U Transport: max ', MAXVAL(uu), ' min ', MINVAL(uu)

    CALL sub_memory(-SIZE(e3u),'r','e3u','sub_input_data_roms')
    DEALLOCATE(e3u)

    !=== V ===!
    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' ------ Meridional transport is computing ------'
    DO l = 1, dims_reg(4,3)
      DO k = 1, dims_reg(3,3)
        DO j = 1, dims_reg(2,3)
          DO i = 1, dims_reg(1,3)
            vv(i,j,k,l) =  vv(i,j,k,l) * e3v(i,j,k,l) * e1v(i,j,1,1)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    WRITE(lun_standard,*)' - V Transport: max ', MAXVAL(vv), ' min ', MINVAL(vv)

    CALL sub_memory(-SIZE(e3v),'r','e3v','sub_input_data_roms')
    DEALLOCATE(e3v)


    !=== W ===!
    !- Must be no divergent -!
    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' ------ Vertical transport is computing ------'
    ww(:,:,dims_reg(3,3),:) = 0._rprec

    DO l = 1, dims_reg(4,3)
      DO k = dims_reg(3,3)-1, 1, -1
        DO j = 2, dims_reg(2,3)
          DO i = 2, dims_reg(1,3)
            ww(i,j,k,l) = ww(i,j,k+1,l) &
                 + uu(i-1,   j, k, l) - uu(i, j, k, l) &
                 + vv(  i, j-1, k, l) - vv(i, j, k, l)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO l = 1, dims_reg(4,3)
      DO k = 2, dims_reg(3,3)-1
        DO j = 1, dims_reg(2,3)
          DO i = 1, dims_reg(1,3)
            ww(i,j,k,l) = &
                 +   ww(i,j,k,l) &
                 - ( ww(i,j,1,l) / ( zz_ww(i,j,1,l) - zz_ww(i,j,dims_reg(3,3),l) ) ) &
                 *                 ( zz_ww(i,j,k,l) - zz_ww(i,j,dims_reg(3,3),l) )
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ww(:,:, 1,:) = 0._rprec
    ww(:,:,dims_reg(3,3),:) = 0._rprec

    WRITE(lun_standard,*)' - W Transport: max ', MAXVAL(ww), ' min ', MINVAL(ww)

  END SUBROUTINE sub_input_data_roms
  !!***
  !=========================================================================
  !!****f* mod_input_data_seq/sub_input_data_seq_roms()
  !! NAME
  !!   sub_input_data_seq_roms()
  !!
  !! FUNCTION
  !!   Data used in Ariane are on C grid and they are based on 
  !!   OPA conventions. ROMS data are on a C grid but don't respect
  !!   all OPA conventions.
  !!   The vertical axis is reversed.
  !!   U is defined only on (imt-1, jmt  , kmt-1)
  !!   V is defined only on (imt  , jmt-1, kmt-1)
  !!   W must be computed
  !!   T and S are defined on (imt, jmt, kmt-1)
  !!
  !!   * read ROMS tracer variables.
  !!   * loop on current and tracer variables: U, V, [W, T, S, R].
  !!   * call read netcdf.
  !!   * compute scale factor e3t.
  !!   * compute transport from currents.
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
  !!   * no arguments
  !!
  !! TODO
  !!   
  !!
  !! USED BY
  !!   * trajec
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_input_data_seq_roms(&
       ncids                       , &
       varids                      , &
       new_file                    , &
       ind_file                    , &
       ind_time                    , &
       ind_time_size               , &
       dimsorders                  )


    INTEGER(kind = iprec) , DIMENSION(:), INTENT(inout) :: ncids
    INTEGER(kind = iprec) , DIMENSION(:), INTENT(inout) :: varids
    LOGICAL               , DIMENSION(:), INTENT(inout) :: new_file
    INTEGER(kind = iprec) , DIMENSION(:), INTENT(inout) :: ind_file
    INTEGER(kind = iprec) , DIMENSION(:), INTENT(inout) :: ind_time
    INTEGER(kind = iprec) , DIMENSION(:), INTENT(inout) :: ind_time_size
    INTEGER(kind = iprec) , DIMENSION(:,:), INTENT(inout) :: dimsorders

    INTEGER(kind=iprec) :: i,j,k,l
    INTEGER(kind=iprec) :: ncid
    INTEGER(kind=iprec) :: varid
    INTEGER(kind = ishort) :: is_err      !cjmp!
    INTEGER(kind = iprec) :: Vstretching  !cjmp!
    INTEGER(kind = iprec) :: Vtransform   !cjmp!
    CHARACTER(len=4), PARAMETER :: c_none='NONE'
    CHARACTER(len=11 ), PARAMETER :: nc_glbatt_Vstretching='Vstretching' !cjmp!
    CHARACTER(len=10 ), PARAMETER :: nc_glbatt_Vtransform='Vtransform'   !cjmp!

    LOGICAL :: key_ijt

    !!NG: 05/08/2008 REAL(kind=rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: &
    !!NG: 05/08/2008         zeta, &
    !!NG: 05/08/2008         zw0 , &
    !!NG: 05/08/2008         e3u , &
    !!NG: 05/08/2008         e3v

    REAL(kind = rprec) :: hc ! S-coordinate parameter, critical depth (m)
    REAL(kind = rprec), DIMENSION(:), ALLOCATABLE :: &
         sc_w , & ! S-coordinate at W-points
         cs_w     ! S-coordinate stretching curves at W-points
    !!NG: 05/08/2008 CHARACTER(len = 128)  :: c_filename   ! netcdf file name

    !======================= ZONAL CURRENT ==================!
    IF (key_interp_temporal) THEN
      uu(:,:,:,2) = uu(:,:,:,1) ! We are in temporal interpolation case
    ENDIF
    !- Read U on ROMS grid -!
    CALL sub_input_data_seq_main(                         &
         ncids(1), varids(1)                            , &
         new_file(1), ind_file(1)                       , &
         ind_time(1), ind_time_size(1)                  , &
         dimsorders(:,1)                                , &
         c_dir_zo,c_prefix_zo, maxsize_zo, c_suffix_zo  , &
         nc_var_zo,nc_var_eivu,nc_att_mask_zo           , &
         uu(1:dims_reg(1,3)-1,:,dims_reg(3,3)-1:1:-1,1:1))

    !============= MERIDIONAL CURRENT =======================!
    IF (key_interp_temporal) THEN
      vv(:,:,:,2) = vv(:,:,:,1) ! We are in temporal interpolation case
    ENDIF
    !- Read V on ROMS grid -!
    CALL sub_input_data_seq_main(                         &
         ncids(2), varids(2)                            , &
         new_file(2), ind_file(2)                       , &
         ind_time(2), ind_time_size(2)                  , &
         dimsorders(:,2)                                , &
         c_dir_me,c_prefix_me,maxsize_me, c_suffix_me   , &
         nc_var_me,nc_var_eivv,nc_att_mask_me           , &
         vv(:,1:dims_reg(2,3)-1,dims_reg(3,3)-1:1:-1,1:1))

    IF (key_alltracers) THEN
      !==================== TEMPERATURE =========================!
      IF (key_interp_temporal) THEN
        tt(:,:,:,2) = tt(:,:,:,1) ! We are in temporal interpolation case
      ENDIF
      !- Read T on ROMS grid -!
      CALL sub_input_data_seq_main(                       &
           ncids(4), varids(4)                            , &
           new_file(4), ind_file(4)                       , &
           ind_time(4), ind_time_size(4)                  , &
           dimsorders(:,4)                                , &
           c_dir_te,c_prefix_te,maxsize_te, c_suffix_te   , &
           nc_var_te,c_none,nc_att_mask_te                , &
           tt(:,:,dims_reg(3,3)-1:1:-1,1:1))

      !====================== SALINITY ===========================!
      IF (key_interp_temporal) THEN
        ss(:,:,:,2) = ss(:,:,:,1) ! We are in temporal interpolation case
      ENDIF
      !- Read S on ROMS grid -!
      CALL sub_input_data_seq_main(                       &
           ncids(5), varids(5)                            , &
           new_file(5), ind_file(5)                       , &
           ind_time(5), ind_time_size(5)                  , &
           dimsorders(:,5)                                , &
           c_dir_sa,c_prefix_sa, maxsize_sa, c_suffix_sa  , &
           nc_var_sa,c_none,nc_att_mask_sa                , &
           ss(:,:,dims_reg(3,3)-1:1:-1,1:1))

      !====================== COMPUTE RR =========================!
      IF (key_interp_temporal) THEN
        rr(:,:,:,2) = rr(:,:,:,1) ! We are in temporal interpolation case
      ENDIF
      WRITE(lun_standard,*)''
      WRITE(lun_standard,*)'------ Density is computing -------'
      !-- Compute Density from Temperature and Salinity --!
      DO l = 1, 1
        DO k = 1,  dims_reg(3,3)
          DO j = 1,  dims_reg(2,3)
            DO i = 1,  dims_reg(1,3)
              rr(i,j,k,l)=rhostp(1000._rprec,ss(i,j,k,l),tt(i,j,k,l),-ABS(zsigma))
            END DO
          END DO
        END DO
      ENDDO
      !- Comments -!
      IF (id_comments) THEN
        WRITE(lun_standard,*)' - Density: max ', MAXVAL(rr), ' min ', MINVAL(rr)
      ENDIF

      !====================== ZETA 2D + Time =======================!
      !- Read Zeta on ROMS grid -!
      IF (.NOT.ALLOCATED(zeta_a)) THEN
        ALLOCATE(zeta_a(dims_reg(1,3),dims_reg(2,3), 1, 1))
        CALL sub_memory(SIZE(zeta_a),'r','zeta_a','sub_input_data_seq_roms')
      ENDIF

      key_ijt = .TRUE.

      CALL sub_input_data_seq_main(                         &
           ncids(7), varids(7)                            , &
           new_file(7), ind_file(7)                       , &
           ind_time(7), ind_time_size(7)                  , &
           dimsorders(:,7)                                , &
           c_dir_ze,c_prefix_ze, maxsize_ze, c_suffix_ze  , &
           nc_var_ze,c_none,nc_att_mask_ze, &
           zeta_a(:,:,:,:), key_ijt = key_ijt)

    ELSE ! key_alltracers

      !====================== ZETA 2D + Time =======================!
      !- Read Zeta on ROMS grid -!
      IF (.NOT.ALLOCATED(zeta_a)) THEN
        ALLOCATE(zeta_a(dims_reg(1,3),dims_reg(2,3), 1, 1))
        CALL sub_memory(SIZE(zeta_a),'r','zeta_a','sub_input_data_seq_roms')
      ENDIF

      key_ijt = .TRUE.

      CALL sub_input_data_seq_main(                         &
           ncids(3), varids(3)                            , &
           new_file(3), ind_file(3)                       , &
           ind_time(3), ind_time_size(3)                  , &
           dimsorders(:,3)                                , &
           c_dir_ze,c_prefix_ze, maxsize_ze, c_suffix_ze  , &
           nc_var_ze,c_none,nc_att_mask_ze, &
           zeta_a(:,:,:,:), key_ijt = key_ijt)

    ENDIF ! END key_alltracers

    !! NG: 05/08/2008: We assume that ZW0 is constant in time
    !! NG: 05/08/2008: It will be read only one time...

    IF (.NOT.ALLOCATED(zw0)) THEN

      !! NG: 05/08/2008!========================= ZW0 ==============================!
      !! NG: 05/08/2008!- Read global attributes hc (scalar) sc_w and cs_w (vectors) -!
      !! NG: 05/08/2008!- in AVG files.
      !! NG: 05/08/2008!- Build file name from temperature file -!
      !! NG: 05/08/2008 CALL sub_build_filename(ind0_te, maxsize_te, &
      !! NG: 05/08/2008     c_prefix_te, c_suffix_te, c_filename)

      !======== Read ROMS global attributes hc, sc_w and Cs_w ======!
      IF (.NOT.ALLOCATED(sc_w)) THEN
        ALLOCATE(sc_w(dims_reg(3,3)))
        !!NG 1dec2009 write(*,*)' sc_w: ',size(sc_w)
        CALL sub_memory(SIZE(sc_w),'r','sc_w','sub_input_data_seq_roms')
      ENDIF

      IF (.NOT.ALLOCATED(cs_w)) THEN
        ALLOCATE(cs_w(dims_reg(3,3)))
        !!NG 1dec2009 write(*,*)' cs_w: ',size(cs_w)
        CALL sub_memory(SIZE(cs_w),'r','cs_w','sub_input_data_seq_roms')
      ENDIF

      !! NG : 29/09/2009
      IF (roms_global_attribute) THEN 

        CALL sub_read_netcdf_global_att(                   &
             dir_glbatt, fn_glbatt                       , &
             nc_glbatt_hc, nc_glbatt_sc_w, nc_glbatt_Cs_w, &
             hc, sc_w, Cs_w)

      ELSE

        !- Open NetCDF file -!
        CALL sub_open_netcdf_file(dir_glbatt, fn_glbatt, ncid)

        !- Read hc -!
        CALL sub_read_netcdf_varid_ndims(ncid, nc_glbatt_hc, varid)
        CALL sub_read_netcdf_var( ncid, varid, hc)

        !- Read sc_w -!

        !cjmp modified to calculate sc_w (called s_w in ROMS 3)
        !analytically, since sc_w does not vary by coordinate
        !transform.  NOTE WELL, THE VERTICAL INDEX IS REVERSED
        !BETWEEN ROMS AND ARIANE

        !cjmp_old CALL sub_read_netcdf_varid_ndims(ncid, nc_glbatt_sc_w, varid)
        !cjmp_old CALL sub_read_netcdf_var1d( ncid, varid, sc_w(dims_reg(3,3):1:-1))
        DO k=0,dims_reg(3,3)-1
          sc_w(dims_reg(3,3)-k)=(1.0*k-(dims_reg(3,3)-1))/(dims_reg(3,3)-1)
        ENDDO
        !cjmp end modification

        !- Read Cs_w -!
        CALL sub_read_netcdf_varid_ndims(ncid, nc_glbatt_Cs_w, varid)
        CALL sub_read_netcdf_var1d( ncid, varid, Cs_w(dims_reg(3,3):1:-1))

        !- Close netcdf file -!
        CALL sub_close_netcdf_file(ncid)

      ENDIF

      !cjmp modified to attempt to read Vtransform and Vstretching from ROMS !
      !cjmp history files.  If they are not present in file, assume Vtransform=1  !
      !cjmp and Vstretching=1 ; Note, this requires checking if the variable exists !
      !cjmp before trying to read it, and also reading it directly, for the routines
      !cjmp in mod_netcdf assume the variable has type real. !
      !- Open NetCDF file -!
      !!NG: 16 nov 2018 add Vtransform 2 for Pierre Amael
      !! vtransform=1: z_unperturbated (s) = hc * (sc(s) - Cs(s)) + Cs(s) * h
      !! vtransform=2: z_unperturbated (s) = h * (sc(s) * hc + Cs(s) * h) / (h + hc)
      !!NG: 16 nov 2018
      CALL sub_open_netcdf_file(dir_glbatt, fn_glbatt, ncid)

      !now read Vtransform
      !first, check if it exists
      is_err     = nf90_inq_varid( &
           ncid  = ncid          , &
           name  = TRIM(nc_glbatt_Vtransform)  , &
           varid = varid           &
           )
      !- Test if Vtransform exists -!
      IF (is_err /= nf90_noerr) THEN
        WRITE(lun_standard,*)' Vtransform does not exits, setting to 1 '
        Vtransform=1
      ELSE
        is_err=nf90_get_var(ncid=ncid,varid=varid,values=Vtransform)
        IF (is_err /= nf90_noerr) THEN
          WRITE(lun_error,*) 'Error reading Vtransform'
        ELSE
          WRITE(lun_standard,*) 'Vtransform is ',Vtransform
        ENDIF
      ENDIF

      !now read Vstretching
      !first, check if it exists
      is_err     = nf90_inq_varid( &
           ncid  = ncid          , &
           name  = TRIM(nc_glbatt_Vstretching)  , &
           varid = varid           &
           )
      !- Test if Vstretching exists -!
      IF (is_err /= nf90_noerr) THEN
        WRITE(lun_standard,*)' Vstretching does not exits, setting to 1 '
        Vstretching=1
      ELSE
        is_err=nf90_get_var(ncid=ncid,varid=varid,values=Vstretching)
        IF (is_err /= nf90_noerr) THEN
          WRITE(lun_error,*) 'Error reading Vstretching'
        ELSE
          WRITE(lun_standard,*) 'Vstretching is ',Vstretching
        ENDIF
      ENDIF

      !- Close netcdf file -!
      CALL sub_close_netcdf_file(ncid)

      !cjmp check if Vstretching or Vtransform .neq. 1; if so, stop
      IF (Vstretching /= 1) THEN
        WRITE(lun_error,*) 'Currently, ARIANE only supports Vstretching=1'
        STOP
      ENDIF
      IF ((Vtransform > 2).or.(Vtransform < 1)) THEN
        WRITE(lun_error,*) 'Currently, ARIANE only supports Vtransform=1 or 2'
        STOP
      ENDIF
      !cjmp end modification

      !============== ZW0 (3D - constant in time) =================!
      WRITE(lun_standard,*)''
      WRITE(lun_standard,*)' ------ ZW0 is computing ------'
      WRITE(lun_standard,*)' - hc: ',hc
      WRITE(lun_standard,*)' - SC_W: ',sc_w
      WRITE(lun_standard,*)' - CS_w: ',cs_w

      IF (.NOT.ALLOCATED(zw0)) THEN
        ALLOCATE(zw0(dims_reg(1,3), dims_reg(2,3), dims_reg(3,3), 1))
        CALL sub_memory(SIZE(zw0),'r','zw0','sub_input_data_seq_roms')
      ENDIF

      !!NG: 16 nov 2018 add Vtransform 2 for Pierre Amael
      !! vtransform=1: z_unperturbated (s) = hc * (sc(s) - Cs(s)) + Cs(s) * h
      !! vtransform=2: z_unperturbated (s) = h * (sc(s) * hc + Cs(s) * h) / (h + hc)
      !!NG: 16 nov 2018
      IF (Vtransform ==1) THEN
        DO k = 1, dims_reg(3,3)
          DO j = 1, dims_reg(2,3)
            DO i = 1, dims_reg(1,3)
              zw0(i,j,k,1) = hc * sc_w(k) + &
                   (h_roms(i,j,1,1) - hc) * cs_w(k)
            ENDDO
          ENDDO
        ENDDO
      ELSEIF (Vtransform ==2) THEN
        DO k = 1, dims_reg(3,3)
          DO j = 1, dims_reg(2,3)
            DO i = 1, dims_reg(1,3)
              zw0(i,j,k,1) = h_roms(i,j,1,1) * &
                   ( sc_w(k) * hc + cs_w(k) * h_roms(i,j,1,1)) / &
                   (h_roms(i,j,1,1) + hc)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      !!    DEALLOCATE(sc_w)
      !!    DEALLOCATE(cs_w)
      !! NG:10/07/1009
      CALL sub_memory(-SIZE(sc_w),'r','sc_w','sub_input_approx_zz_ww_roms')
      DEALLOCATE(sc_w)
      CALL sub_memory(-SIZE(cs_w),'r','cs_w','sub_input_approx_zz_ww_roms')
      DEALLOCATE(cs_w)
      !! NG:10/07/1009

    ELSE

      WRITE(lun_standard,*)''
      WRITE(lun_standard,*) &
           '  ------  We assume that ZW0 is constant in time. ------'

    ENDIF

    IF (id_comments) THEN
      WRITE(lun_standard,*)' - ZW0 : max ', MAXVAL(zw0), ' min ', MINVAL(zw0)
    ENDIF

    !==================== ZZ_WW (3D + time) ========================!
    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' ------ ZZ_WW is computing ------'

    IF (.NOT.ALLOCATED(zz_ww)) THEN
      ALLOCATE(zz_ww(dims_reg(1,3), dims_reg(2,3), dims_reg(3,3),1))
      CALL sub_memory(SIZE(zz_ww),'r','zz_ww','sub_input_data_seq_roms')
    ENDIF

    DO l = 1, 1
      DO k = 1, dims_reg(3,3)
        DO j = 1, dims_reg(2,3)
          DO i = 1, dims_reg(1,3)
            zz_ww(i,j,k,l) = &
                 zeta_a(i,j,1,l) + zw0(i,j,k,1) * ( 1. + zeta_a(i,j,1,l) / h_roms(i,j,1,1))
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    IF (id_comments) THEN
      WRITE(lun_standard,*)' - ZZ_WW : max ', MAXVAL(zz_ww), ' min ', MINVAL(zz_ww)
    ENDIF
    !!    DEALLOCATE(zw0)

    !=================== E3T (3D + time) ========================!
    !- Compute E3T from ZZ_WW -!
    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' ------ E3T scale factor is computing ------'

    IF (.NOT.ALLOCATED(e3t)) THEN
      ALLOCATE(e3t(dims_reg(1,3), dims_reg(2,3), dims_reg(3,3), 1))
      CALL sub_memory(SIZE(e3t),'r','e3t','sub_input_data_seq_roms')
    ENDIF

    e3t(:,:,1:dims_reg(3,3)-1,:)= &
         zz_ww(:,:,1:dims_reg(3,3)-1,:) - zz_ww(:,:,2:dims_reg(3,3),:)

    e3t(:,:,dims_reg(3,3),:) = e3t(:,:,dims_reg(3,3)-1,:)

    IF (id_comments) THEN
      WRITE(lun_standard,*)' - E3T : max ', MAXVAL(e3t), ' min ', MINVAL(e3t)
    ENDIF
    !========================= E3U ==============================!
    IF (.NOT.ALLOCATED(e3u)) THEN
      ALLOCATE(e3u(dims_reg(1,3), dims_reg(2,3), dims_reg(3,3), 1))
      CALL sub_memory(SIZE(e3u),'r','e3u','sub_input_data_seq_roms')
    ENDIF

    e3u(1:dims_reg(1,3)-1,:,:,:) = &
         .5_rprec *( e3t(1:dims_reg(1,3)-1,:,:,:) + e3t(2:dims_reg(1,3),:,:,:) )

    e3u(dims_reg(1,3),:,:,:) = e3u(dims_reg(1,3)-1,:,:,:)

    !========================= E3V ==============================!
    IF (.NOT.ALLOCATED(e3v)) THEN
      ALLOCATE(e3v(dims_reg(1,3), dims_reg(2,3), dims_reg(3,3), 1))
      CALL sub_memory(SIZE(e3v),'r','e3u','sub_input_data_seq_roms')
    ENDIF

    e3v(:,1:dims_reg(2,3)-1,:,:) = &
         .5_rprec *( e3t(:,1:dims_reg(2,3)-1,:,:) + e3t(:,2:dims_reg(2,3),:,:) )

    e3v(:,dims_reg(2,3),:,:) = e3v(:,dims_reg(2,3)-1,:,:)

    !====================== TRANSPORT ===========================!
    !=== U ===!
    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' ------ Zonal transport is computing ------'
    DO l = 1, 1
      DO k = 1, dims_reg(3,3)
        DO j = 1, dims_reg(2,3)
          DO i = 1, dims_reg(1,3)
            uu(i,j,k,l) =  uu(i,j,k,l) * e3u(i,j,k,l) * e2u(i,j,1,1)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    IF (id_comments) THEN
      WRITE(lun_standard,*)' - U Transport: max ', MAXVAL(uu), ' min ', MINVAL(uu)
    ENDIF
    !!    DEALLOCATE(e3u)

    !=== V ===!
    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' ------ Meridional transport is computing ------'
    DO l = 1, 1
      DO k = 1, dims_reg(3,3)
        DO j = 1, dims_reg(2,3)
          DO i = 1, dims_reg(1,3)
            vv(i,j,k,l) =  vv(i,j,k,l) * e3v(i,j,k,l) * e1v(i,j,1,1)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    IF (id_comments) THEN
      WRITE(lun_standard,*)' - V Transport: max ', MAXVAL(vv), ' min ', MINVAL(vv)
    ENDIF
    !!    DEALLOCATE(e3v)

    !=== W ===!
    IF (key_interp_temporal) THEN
      ww(:,:,:,2) = ww(:,:,:,1) ! We are in temporal interpolation case
    ENDIF
    !- Must be no divergent -!
    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' ------ Vertical transport is computing ------'
    ww(:,:,dims_reg(3,3),1:1) = 0._rprec

    DO l = 1, 1
      DO k = dims_reg(3,3)-1, 1, -1
        DO j = 2, dims_reg(2,3)
          DO i = 2, dims_reg(1,3)
            ww(i,j,k,l) = ww(i,j,k+1,l) &
                 + uu(i-1,   j, k, l) - uu(i, j, k, l) &
                 + vv(  i, j-1, k, l) - vv(i, j, k, l)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO l = 1, 1
      DO k = 2, dims_reg(3,3)-1
        DO j = 1, dims_reg(2,3)
          DO i = 1, dims_reg(1,3)
            ww(i,j,k,l) = &
                 +   ww(i,j,k,l) &
                 - ( ww(i,j,1,l) / ( zz_ww(i,j,1,l) - zz_ww(i,j,dims_reg(3,3),l) ) ) &
                 *                 ( zz_ww(i,j,k,l) - zz_ww(i,j,dims_reg(3,3),l) )
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ww(:,:, 1,:) = 0._rprec
    ww(:,:,dims_reg(3,3),:) = 0._rprec

    IF (id_comments) THEN
      WRITE(lun_standard,*)' - W Transport: max ', MAXVAL(ww), ' min ', MINVAL(ww)
    ENDIF

  END SUBROUTINE sub_input_data_seq_roms
  !!***
  !=========================================================================
  !!****f* mod_input_data/sub_input_approx_zz_ww_roms()
  !! NAME
  !!  sub_input_approx_zz_ww_roms ()
  !!
  !! FUNCTION
  !!
  !! AUTHOR
  !!   * Origin  : Nicolas Grima (February 2008)
  !! 
  !! CREATION DATE
  !!   * February 2008
  !!
  !! HISTORY
  !!   Date (dd/mm/yyyy/) - Modification(s)
  !!
  !! ARGUMENTS
  !!   * no arguments
  !!
  !! TODO
  !!   
  !!
  !! USED BY
  !!   * trajec_seq
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_input_approx_zz_ww_roms()

    !! WE ASSUME HERE THAT ZETA=0. !!

    REAL(kind=rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: zw0

    REAL(kind = rprec) :: hc ! S-coordinate parameter, critical depth (m)
    REAL(kind = rprec), DIMENSION(:), ALLOCATABLE :: &
         sc_w , & ! S-coordinate at W-points
         cs_w     ! S-coordinate stretching curves at W-points
    !!NG: 16 march 2011 CHARACTER(len = 128)  :: c_filename   ! netcdf file name

    INTEGER(kind=iprec) :: i,j,k,l
    INTEGER(kind=iprec) :: ncid
    INTEGER(kind=iprec) :: varid

    !!NG: 16 march 2011 !========================= ZW0 ==============================!
    !!NG: 16 march 2011 !- Read global attributes hc (scalar) sc_w and cs_w (vectors) -!
    !!NG: 16 march 2011 !- in AVG files.
    !!NG: 16 march 2011 !- Build file name from temperature file -!
    !!NG: 16 march 2011 CALL sub_build_filename(ind0_te, maxsize_te, &
    !!NG: 16 march 2011      c_prefix_te, c_suffix_te, c_filename)

    !======== Read ROMS global attributes hc, sc_w and Cs_w ======!
    IF (.NOT.ALLOCATED(sc_w)) THEN
      ALLOCATE(sc_w(dims_reg(3,3)))
      CALL sub_memory(SIZE(sc_w),'r','sc_w','sub_input_approx_zz_ww_roms')
    ENDIF

    IF (.NOT.ALLOCATED(cs_w)) THEN
      ALLOCATE(cs_w(dims_reg(3,3)))
      CALL sub_memory(SIZE(cs_w),'r','cs_w','sub_input_approx_zz_ww_roms')
    ENDIF

    !! NG : 29/09/2009
    IF (roms_global_attribute) THEN 
      CALL sub_read_netcdf_global_att(                   &
           dir_glbatt, fn_glbatt                       , &
           nc_glbatt_hc, nc_glbatt_sc_w, nc_glbatt_Cs_w, &
           hc, sc_w, Cs_w)
    ELSE

      !- Open NetCDF file -!
      CALL sub_open_netcdf_file(dir_glbatt, fn_glbatt, ncid)

      !- Read hc -!
      CALL sub_read_netcdf_varid_ndims(ncid, nc_glbatt_hc, varid)
      CALL sub_read_netcdf_var( ncid, varid, hc)

      !- Read sc_w -!
      CALL sub_read_netcdf_varid_ndims(ncid, nc_glbatt_sc_w, varid)
      CALL sub_read_netcdf_var1d( ncid, varid, sc_w)

      !- Read Cs_w -!
      CALL sub_read_netcdf_varid_ndims(ncid, nc_glbatt_Cs_w, varid)
      CALL sub_read_netcdf_var1d( ncid, varid, Cs_w)

      !- Close netcdf file -!
      CALL sub_close_netcdf_file(ncid)

    ENDIF

    !============== ZW0 (3D - constant in time) =================!
    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' ------ ZW0 is computing ------'
    WRITE(lun_standard,*)' - hc: ',hc
    WRITE(lun_standard,*)' - SC_W: ',sc_w
    WRITE(lun_standard,*)' - CS_w: ',cs_w

    IF (.NOT.ALLOCATED(zw0)) THEN
      ALLOCATE(zw0(dims_reg(1,3), dims_reg(2,3), dims_reg(3,3), 1))
      CALL sub_memory(SIZE(zw0),'r','zw0','sub_input_approx_zz_ww_roms')
    ENDIF

    DO k = 1, dims_reg(3,3)
      DO j = 1, dims_reg(2,3)
        DO i = 1, dims_reg(1,3)
          zw0(i,j,k,1) = hc * sc_w(k) + &
               (h_roms(i,j,1,1) - hc) * cs_w(k)
        ENDDO
      ENDDO
    ENDDO

    IF (id_comments) THEN
      WRITE(lun_standard,*)' - ZW0 : max ', MAXVAL(zw0), ' min ', MINVAL(zw0)
    ENDIF

    !==================== ZZ_WW (3D + time) ========================!
    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' ------ ZZ_WW is computing ------'

    IF (.NOT.ALLOCATED(zz_ww)) THEN
      ALLOCATE(zz_ww(dims_reg(1,3), dims_reg(2,3), dims_reg(3,3),1))
      CALL sub_memory(SIZE(zz_ww),'r','zz_ww','sub_input_approx_zz_ww_roms')
    ENDIF

    DO l = 1, 1
      DO k = 1, dims_reg(3,3)
        DO j = 1, dims_reg(2,3)
          DO i = 1, dims_reg(1,3)
            zz_ww(i,j,k,l) = zw0(i,j,k,1)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    IF (id_comments) THEN
      WRITE(lun_standard,*)' - ZZ_WW : max ', MAXVAL(zz_ww), ' min ', MINVAL(zz_ww)
    ENDIF

    CALL sub_memory(-SIZE(zw0),'r','zw0','sub_input_approx_zz_ww_roms')
    DEALLOCATE(zw0)
    CALL sub_memory(-SIZE(sc_w),'r','sc_w','sub_input_approx_zz_ww_roms')
    DEALLOCATE(sc_w)
    CALL sub_memory(-SIZE(cs_w),'r','cs_w','sub_input_approx_zz_ww_roms')
    DEALLOCATE(cs_w)


  END SUBROUTINE sub_input_approx_zz_ww_roms


  !!***
  !=========================================================================
  !!****f* mod_input_data_seq/sub_input_data_seq_mars()
  !! NAME
  !!   sub_input_data_seq_roms()
  !!
  !! FUNCTION
  !!   Data used in Ariane are on C grid and they are based on 
  !!   OPA conventions.
  !!
  !!   * read MARS tracer variables.
  !!   * loop on current and tracer variables: U, V, [W, T, S, R].
  !!   * call read netcdf.
  !!   * compute scale factor e3t.
  !!   * compute transport from currents.
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
  !!   * no arguments
  !!
  !! TODO
  !!   
  !!
  !! USED BY
  !!   * trajec
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_input_data_seq_mars(&
       ncids                       , &
       varids                      , &
       new_file                    , &
       ind_file                    , &
       ind_time                    , &
       ind_time_size               , &
       dimsorders                  , &
       lref)


    INTEGER(kind = iprec) , DIMENSION(:), INTENT(inout) :: ncids
    INTEGER(kind = iprec) , DIMENSION(:), INTENT(inout) :: varids
    LOGICAL               , DIMENSION(:), INTENT(inout) :: new_file
    INTEGER(kind = iprec) , DIMENSION(:), INTENT(inout) :: ind_file
    INTEGER(kind = iprec) , DIMENSION(:), INTENT(inout) :: ind_time
    INTEGER(kind = iprec) , DIMENSION(:), INTENT(inout) :: ind_time_size
    INTEGER(kind = iprec) , DIMENSION(:,:), INTENT(inout) :: dimsorders
    INTEGER(kind = iprec) , OPTIONAL      , INTENT(in)  :: lref

    INTEGER(kind=iprec) :: i,j,k,l, ii
    INTEGER(kind=iprec) :: ncid
    INTEGER(kind=iprec) :: varid
    INTEGER(kind=iprec) :: is_err
    INTEGER(kind=iprec) :: iter
    CHARACTER(len=4), PARAMETER :: c_none='NONE'

    LOGICAL :: key_ijt

    IF (PRESENT(lref)) THEN
      iter = lref
    ELSE
      iter = 0
    ENDIF

    !======================= ZONAL CURRENT ==================!
    IF (key_interp_temporal) THEN
      uu(:,:,:,2) = uu(:,:,:,1) ! We are in temporal interpolation case
    ENDIF
    !- Read U on MARS grid -!
    CALL sub_input_data_seq_main(                         &
         ncids(1), varids(1)                            , &
         new_file(1), ind_file(1)                       , &
         ind_time(1), ind_time_size(1)                  , &
         dimsorders(:,1)                                , &
         c_dir_zo,c_prefix_zo, maxsize_zo, c_suffix_zo  , &
         nc_var_zo,nc_var_eivu,nc_att_mask_zo           , &
         uu(:,:,dims_reg(3,3):1:-1,1:1))  ! reverse z dimension !

    !============= MERIDIONAL CURRENT =======================!
    IF (key_interp_temporal) THEN
      vv(:,:,:,2) = vv(:,:,:,1) ! We are in temporal interpolation case
    ENDIF
    !- Read V on MARS grid -!
    CALL sub_input_data_seq_main(                         &
         ncids(2), varids(2)                            , &
         new_file(2), ind_file(2)                       , &
         ind_time(2), ind_time_size(2)                  , &
         dimsorders(:,2)                                , &
         c_dir_me,c_prefix_me,maxsize_me, c_suffix_me   , &
         nc_var_me,nc_var_eivv,nc_att_mask_me           , &
         vv(:,:,dims_reg(3,3):1:-1,1:1))  ! reverse z dimension !

    IF (key_alltracers) THEN
      !==================== TEMPERATURE =========================!
      IF (key_interp_temporal) THEN
        tt(:,:,:,2) = tt(:,:,:,1) ! We are in temporal interpolation case
      ENDIF
      !- Read T on MARS grid -!
      CALL sub_input_data_seq_main(                      &
           ncids(4), varids(4)                         , &
           new_file(4), ind_file(4)                    , &
           ind_time(4), ind_time_size(4)               , &
           dimsorders(:,4)                             , &
           c_dir_te,c_prefix_te,maxsize_te, c_suffix_te, &
           nc_var_te,c_none,nc_att_mask_te             , &
           tt(:,:,dims_reg(3,3):1:-1,1:1))               ! reverse z dimension !

      !====================== SALINITY ===========================!
      IF (key_interp_temporal) THEN
        ss(:,:,:,2) = ss(:,:,:,1) ! We are in temporal interpolation case
      ENDIF
      !- Read S on MARS grid -!
      CALL sub_input_data_seq_main(                       &
           ncids(5), varids(5)                          , &
           new_file(5), ind_file(5)                     , &
           ind_time(5), ind_time_size(5)                , &
           dimsorders(:,5)                              , &
           c_dir_sa,c_prefix_sa, maxsize_sa, c_suffix_sa, &
           nc_var_sa,c_none,nc_att_mask_sa              , &
           ss(:,:,dims_reg(3,3):1:-1,1:1))                ! reverse z dimension !

      !====================== COMPUTE RR =========================!
      IF (key_interp_temporal) THEN
        rr(:,:,:,2) = rr(:,:,:,1) ! We are in temporal interpolation case
      ENDIF
      WRITE(lun_standard,*)''
      WRITE(lun_standard,*)'------ Density is computing -------'
      !-- Compute Density from Temperature and Salinity --!
      DO l = 1, 1
        DO k = 1,  dims_reg(3,3)
          DO j = 1,  dims_reg(2,3)
            DO i = 1,  dims_reg(1,3)
              rr(i,j,k,l)=rhostp(1000._rprec,ss(i,j,k,l),tt(i,j,k,l),-ABS(zsigma))
            END DO
          END DO
        END DO
      ENDDO
      !- Comments -!
      IF (id_comments) THEN
        WRITE(lun_standard,*)' - Density: max ', MAXVAL(rr), ' min ', MINVAL(rr)
      ENDIF

      !====================== ZETA 2D + Time =======================!
      !- Read Zeta on MARS grid -!
      IF (.NOT.ALLOCATED(zeta_a)) THEN
        ALLOCATE(zeta_a(dims_reg(1,3),dims_reg(2,3), 1, 1))
        CALL sub_memory(SIZE(zeta_a),'r','zeta_a','sub_input_data_seq_mars')
      ENDIF

      key_ijt = .TRUE.

      CALL sub_input_data_seq_main(                         &
           ncids(7), varids(7)                            , &
           new_file(7), ind_file(7)                       , &
           ind_time(7), ind_time_size(7)                  , &
           dimsorders(:,7)                                , &
           c_dir_ze,c_prefix_ze, maxsize_ze, c_suffix_ze  , &
           nc_var_ze,c_none,nc_att_mask_ze, &
           zeta_a(:,:,:,:), key_ijt = key_ijt)               ! key_ijt !!!!

    ELSE ! key_alltracers

      !====================== ZETA 2D + Time =======================!
      !- Read Zeta on MARS grid -!
      IF (.NOT.ALLOCATED(zeta_a)) THEN
        ALLOCATE(zeta_a(dims_reg(1,3),dims_reg(2,3), 1, 1))
        CALL sub_memory(SIZE(zeta_a),'r','zeta_a','sub_input_data_seq_mars')
      ENDIF

      key_ijt = .TRUE.

      CALL sub_input_data_seq_main(                         &
           ncids(3), varids(3)                            , &
           new_file(3), ind_file(3)                       , &
           ind_time(3), ind_time_size(3)                  , &
           dimsorders(:,3)                                , &
           c_dir_ze,c_prefix_ze, maxsize_ze, c_suffix_ze  , &
           nc_var_ze,c_none,nc_att_mask_ze, &
           zeta_a(:,:,:,:), key_ijt = key_ijt)              ! key_ijt !!!!

    ENDIF ! END key_alltracers

    !==================== ZZ_WW (3D + time) ========================!
    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' ------ ZZ_WW is computing ------'

    IF (.NOT.ALLOCATED(zz_ww)) THEN
      ALLOCATE(zz_ww(dims_reg(1,3), dims_reg(2,3), dims_reg(3,3)+1,1))
      CALL sub_memory(SIZE(zz_ww),'r','zz_ww','sub_input_data_seq_mars')
    ENDIF

    DO l = 1, 1
      DO k = 1, dims_reg(3,3)+1
        DO j = 1, dims_reg(2,3)
          DO i = 1, dims_reg(1,3)
            zz_ww(i,j,k,l) = &
                 zeta_a(i,j,1,l) * (1 + sc_w(k)) + &
                 hc * sc_w(k) + &
                 (h_mars(i,j,1,1) - hc) * cs_w(k)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    IF (id_comments) THEN
      WRITE(lun_standard,*)' - ZZ_WW : max ', MAXVAL(zz_ww), ' min ', MINVAL(zz_ww)
      !WRITE(lun_standard,*)' - zz_ww(imt/2,jmt/2,:):', zz_ww(dims_reg(1,3)/2,dims_reg(2,3)/2,:,1)
    ENDIF

    !=================== E3T (3D + time) ========================!
    !- Compute E3T from ZZ_WW -!
    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' ------ E3T scale factor is computing ------'

    IF (.NOT.ALLOCATED(e3t)) THEN
      ALLOCATE(e3t(dims_reg(1,3), dims_reg(2,3), dims_reg(3,3)+1, 1))
      CALL sub_memory(SIZE(e3t),'r','e3t','sub_input_data_seq_mars')
    ENDIF

    e3t(:,:,1:dims_reg(3,3),:)= &
         zz_ww(:,:,1:dims_reg(3,3),:) - zz_ww(:,:,2:dims_reg(3,3)+1,:)

    e3t(:,:,dims_reg(3,3)+1,:) = e3t(:,:,dims_reg(3,3),:)

    IF (id_comments) THEN
      WRITE(lun_standard,*)' - E3T : max ', MAXVAL(e3t), ' min ', MINVAL(e3t)
    ENDIF

    IF (key_write_transport) THEN

      CALL write_netcdf     ( &
           h_mars(:,:,1,1)  , &
           'bathy'          , &
           'scalefactors.nc', &
           iter             )

      CALL write_netcdf     ( &
           zz_ww(:,:,:,1)   , &
           'zz_ww'          , &
           'scalefactors.nc', &
           iter             )

      CALL write_netcdf     ( &
           e3t(:,:,:,1)     , &
           'e3t'            , &
           'scalefactors.nc', &
           iter             )
    ENDIF

    !========================= E3U ==============================!
    IF (.NOT.ALLOCATED(e3u)) THEN
      ALLOCATE(e3u(dims_reg(1,3), dims_reg(2,3), dims_reg(3,3)+1, 1))
      CALL sub_memory(SIZE(e3u),'r','e3u','sub_input_data_seq_mars')
    ENDIF

    e3u(1:dims_reg(1,3)-1,:,:,:) = &
         .5_rprec *( e3t(1:dims_reg(1,3)-1,:,:,:) + e3t(2:dims_reg(1,3),:,:,:) )

    e3u(dims_reg(1,3),:,:,:) = e3u(dims_reg(1,3)-1,:,:,:)

    IF (id_comments) THEN
      WRITE(lun_standard,*)' - E3U : max ', MAXVAL(e3u), ' min ', MINVAL(e3u)
    ENDIF

    !========================= E3V ==============================!
    IF (.NOT.ALLOCATED(e3v)) THEN
      ALLOCATE(e3v(dims_reg(1,3), dims_reg(2,3), dims_reg(3,3)+1, 1))
      CALL sub_memory(SIZE(e3v),'r','e3u','sub_input_data_seq_mars')
    ENDIF

    e3v(:,1:dims_reg(2,3)-1,:,:) = &
         .5_rprec *( e3t(:,1:dims_reg(2,3)-1,:,:) + e3t(:,2:dims_reg(2,3),:,:) )

    e3v(:,dims_reg(2,3),:,:) = e3v(:,dims_reg(2,3)-1,:,:)

    IF (id_comments) THEN
      WRITE(lun_standard,*)' - E3V : max ', MAXVAL(e3v), ' min ', MINVAL(e3v)
    ENDIF

    IF (key_write_transport) THEN
      CALL write_netcdf  ( &
           uu(:,:,:,1)   , &
           'uu'          , &
           'velocity.nc', &
           iter     )
      CALL write_netcdf      ( &
           e3u(:,:,:,1)     , &
           'e3u'             , &
           'scalefactors.nc', &
           iter            )

      CALL write_netcdf      ( &
           e2u(:,:,1,1)     , &
           'dy_uu'          , &
           'scalefactors.nc', &
           iter            )

    ENDIF

    !====================== TRANSPORT ===========================!
    !=== U ===!
    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' ------ Zonal transport is computing ------'
    DO l = 1, 1
      DO k = 1, dims_reg(3,3)
        DO j = 1, dims_reg(2,3)
          DO i = 1, dims_reg(1,3)
            uu(i,j,k,l) =  uu(i,j,k,l) * e3u(i,j,k,l) * e2u(i,j,1,1)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    IF (id_comments) THEN
      WRITE(lun_standard,*)' - U Transport: max ', MAXVAL(uu), ' min ', MINVAL(uu)
    ENDIF
    !!    DEALLOCATE(e3u)

    IF (key_write_transport) THEN
      CALL write_netcdf  ( &
           uu(:,:,:,1)   , &
           'uu'          , &
           'transport.nc', &
           iter     )

      CALL write_netcdf  ( &
           vv(:,:,:,1)  , &
           'v'          , &
           'velocity.nc', &
           iter     )

      CALL write_netcdf     ( &
           e3v(:,:,:,1)     , &
           'e3v'            , &
           'scalefactors.nc', &
           iter             )

      CALL write_netcdf      ( &
           e1v(:,:,1,1)     , &
           'dx_vv'          , &
           'scalefactors.nc', &
           iter            )
    ENDIF

    !=== V ===!
    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' ------ Meridional transport is computing ------'
    DO l = 1, 1
      DO k = 1, dims_reg(3,3)
        DO j = 1, dims_reg(2,3)
          DO i = 1, dims_reg(1,3)
            vv(i,j,k,l) =  vv(i,j,k,l) * e3v(i,j,k,l) * e1v(i,j,1,1)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    IF (id_comments) THEN
      WRITE(lun_standard,*)' - V Transport: max ', MAXVAL(vv), ' min ', MINVAL(vv)
    ENDIF
    !!    DEALLOCATE(e3v)

    IF (key_write_transport) THEN
      CALL write_netcdf  ( &
           vv(:,:,:,1)   , &
           'vv'          , &
           'transport.nc', &
           iter     )
    ENDIF

    !=== W ===!
    IF (key_interp_temporal) THEN
      ww(:,:,:,2) = ww(:,:,:,1) ! We are in temporal interpolation case
    ENDIF
    !- Must be no divergent -!
    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' ------ Vertical transport is computing ------'
    ww(:,:,dims_reg(3,3)+1,1:1) = 0._rprec

    DO l = 1, 1
      DO k = dims_reg(3,3), 1, -1
        DO j = 2, dims_reg(2,3)
          DO i = 2, dims_reg(1,3)
            ww(i,j,k,l) = ww(i,j,k+1,l) &
                 + uu(i-1,   j, k, l) - uu(i, j, k, l) &
                 + vv(  i, j-1, k, l) - vv(i, j, k, l)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO l = 1, 1
      DO k = 2, dims_reg(3,3)
        DO j = 1, dims_reg(2,3)
          DO i = 1, dims_reg(1,3)
            ww(i,j,k,l) = &
                 +   ww(i,j,k,l) &
                 - ( ww(i,j,1,l) / ( zz_ww(i,j,1,l) - zz_ww(i,j,dims_reg(3,3)+1,l) ) ) &
                 *                 ( zz_ww(i,j,k,l) - zz_ww(i,j,dims_reg(3,3)+1,l) )
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ww(:,:, 1,:) = 0._rprec
    ww(:,:,dims_reg(3,3)+1,:) = 0._rprec

    IF (id_comments) THEN
      WRITE(lun_standard,*)' - W Transport: max ', MAXVAL(ww), ' min ', MINVAL(ww)
      !WRITE(lun_standard,*)' - ww(imt/2,jmt/2,:):', ww(dims_reg(1,3)/2,dims_reg(2,3)/2,:,1)
    ENDIF

    IF (key_write_transport) THEN
      CALL write_netcdf  ( &
           ww(:,:,:,1)   , &
           'ww'          , &
           'transport.nc', &
           iter     )
    ENDIF

  END SUBROUTINE sub_input_data_seq_mars

  !!***
  !=========================================================================
  !!****f* mod_input_data/sub_input_approx_zz_ww_mars()
  !! NAME
  !!  sub_input_approx_zz_ww_roms ()
  !!
  !! FUNCTION
  !!
  !! AUTHOR
  !!   * Origin  : Nicolas Grima (November 2018)
  !! 
  !! CREATION DATE
  !!   * November 2018
  !!
  !! HISTORY
  !!   Date (dd/mm/yyyy/) - Modification(s)
  !!
  !! ARGUMENTS
  !!   * no arguments
  !!
  !! TODO
  !!   
  !!
  !! USED BY
  !!   * trajec_seq
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_input_approx_zz_ww_mars()

    !!NG: 16nov2018 NOT TESTED !! NOT TESTED !! NOT TESTED !! NOT TESTED

    REAL(kind=rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: zw0

    INTEGER(kind=iprec) :: i,j,k,l
    INTEGER(kind=iprec) :: ncid
    INTEGER(kind=iprec) :: varid

    !==================== ZZ_WW (3D + time) ========================!
    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' ------ ZZ_WW is computing ------'

    IF (.NOT.ALLOCATED(zz_ww)) THEN
      ALLOCATE(zz_ww(dims_reg(1,3), dims_reg(2,3), dims_reg(3,3)+1,1))
      CALL sub_memory(SIZE(zz_ww),'r','zz_ww','sub_input_approx_zz_ww_roms')
    ENDIF

    DO l = 1, 1
      DO k = 1, dims_reg(3,3)+1
        DO j = 1, dims_reg(2,3)
          DO i = 1, dims_reg(1,3)
            zz_ww(i,j,k,l) = &
                 hc * sc_w(k) + &
                 (h_mars(i,j,1,1) - hc) * cs_w(k)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    IF (id_comments) THEN
      WRITE(lun_standard,*)' - ZZ_WW(approx) : max ', MAXVAL(zz_ww), ' min ', MINVAL(zz_ww)
    ENDIF

    !!NG: 16nov2018 NOT TESTED !! NOT TESTED !! NOT TESTED !! NOT TESTED

  END SUBROUTINE sub_input_approx_zz_ww_mars

  !!***
  !=========================================================================
  !!****f* mod_input_data/sub_input_data_symphonie()
  !! NAME
  !!   sub_input_data_symphonie()
  !!
  !! FUNCTION
  !!   Data used in Ariane are on C grid and they are based on 
  !!   OPA conventions. SYMPHONIE data are on a C grid but don't respect
  !!   all OPA conventions.
  !!   The vertical axis is reversed.
  !!   U is defined only on (imt-1, jmt  , kmt-1)
  !!   V is defined only on (imt  , jmt-1, kmt-1)
  !!   W must be computed
  !!   T and S are defined on (imt, jmt, kmt-1)
  !!
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
  !!   * no arguments
  !!
  !! TODO
  !!   
  !!
  !! USED BY
  !!   * trajec
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_input_data_symphonie()

    INTEGER(kind=iprec) :: i,j,k,l
    CHARACTER(len=4), PARAMETER :: c_none='NONE'

    INTEGER(kind = iprec) :: kmax

    REAL(kind=rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: &
         zw0, &
         e3u  , &
         e3v

    REAL(kind=rprec), DIMENSION(:,:), ALLOCATABLE :: depth_max

    kmax =  dims_reg(3,3)


    !======================= ZONAL CURRENT ==================!
    !- Read U on SYMPHONIE grid -!
    CALL sub_input_data_main(c_dir_zo,c_prefix_zo, &
         ind0_zo,indn_zo,maxsize_zo, &
         c_suffix_zo, nc_var_zo,nc_var_eivu,nc_att_mask_zo, &
         uu(:,:,dims_reg(3,3)-1:1:-1,:))

    uu(:,:,:,:) = CSHIFT(uu(:,:,:,:), shift=1, dim=1)
    uu(dims_reg(1,3),:,:,:) = 0._rprec

    !============= MERIDIONAL CURRENT =======================!
    !- Read V on SYMPHONIE grid -!
    CALL sub_input_data_main(c_dir_me,c_prefix_me,&
         ind0_me,indn_me,maxsize_me, &
         c_suffix_me, nc_var_me,nc_var_eivv,nc_att_mask_me,&
         vv(:,:,dims_reg(3,3)-1:1:-1,:))

    vv(:,:,:,:) = CSHIFT(vv(:,:,:,:), shift=1, dim=2)
    vv(:,dims_reg(2,3),:,:) = 0._rprec

    IF (key_alltracers) THEN
      !==================== TEMPERATURE =========================!
      !- Read T on SYMPHONIE grid -!
      CALL sub_input_data_main(c_dir_te,c_prefix_te,&
           ind0_te,indn_te,maxsize_te, &
           c_suffix_te, nc_var_te,c_none,nc_att_mask_te,&
           tt(:,:,dims_reg(3,3)-1:1:-1,:))

      !====================== SALINITY ===========================!
      !- Read S on SYMPHONIE grid -!
      CALL sub_input_data_main(c_dir_sa,c_prefix_sa, &
           ind0_sa,indn_sa,maxsize_sa, &
           c_suffix_sa, nc_var_sa,c_none,nc_att_mask_sa, &
           ss(:,:,dims_reg(3,3)-1:1:-1,:))



      IF (key_computesigma) THEN

        !====================== COMPUTE RR =========================!
        WRITE(lun_standard,*)''
        WRITE(lun_standard,*)'------ Density is computing -------'
        !-- Compute Density from Temperature and Salinity --!
        DO l = 1, lmt
          DO k = 1,  dims_reg(3,3)-1
            DO j = 1,  dims_reg(2,3)
              DO i = 1,  dims_reg(1,3)
                rr(i,j,k,l)=rhostp(1000._rprec,ss(i,j,k,l),tt(i,j,k,l),-ABS(zsigma))
              END DO
            END DO
          END DO
        ENDDO

      ELSE

        !-- Read Density --!
        CALL sub_input_data_main(                           &
             c_dir_de,c_prefix_de                         , &
             ind0_de,indn_de, maxsize_de                  , &
             c_suffix_de, nc_var_de,c_none,nc_att_mask_de , &
             rr(:,:,dims_reg(3,3)-1:1:-1,1:1) )

      ENDIF
    ENDIF

    !- Comments -!
    WRITE(lun_standard,*)' - Density: max ', MAXVAL(rr), ' min ', MINVAL(rr)

    !====================== READ ZZ_TT =======================!
    !- Read ZZ_TT on SYMPHONIE grid -!
    !- surface is at indice kmax and deep at 0.

    ALLOCATE(zw0(dims_reg(1,3),dims_reg(2,3),dims_reg(3,3), 1))
    CALL sub_memory(SIZE(zw0),'r','zw0','sub_input_data_symphonie')

    zw0(:,:,:,:) = zz_ww(:,:,:,:)

    !====================== READ SSE 2D + Time =======================!
    !- Read SSE on SYMPHONIE grid -!
    !! NG: ALLOCATE(sse(dims_reg(1,3),dims_reg(2,3), 1, dims_reg(4,3)))
    !! NG: call sub_memory(size(sse),'r','sse','sub_input_data_symphonie')

    CALL sub_input_data_main(c_dir_sse,c_prefix_sse, &
         ind0_sse,indn_sse,maxsize_sse, &
         c_suffix_sse, nc_var_sse,c_none,nc_att_mask_sse, &
         sses(:,:,:,:))

    !==================== COMPUTE ZZ_WW (3D + time) =====================!
    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' ------ ZZ_WW is computing ------'
    ALLOCATE(zz_ww(dims_reg(1,3), dims_reg(2,3), dims_reg(3,3), dims_reg(4,3)))
    CALL sub_memory(SIZE(zz_ww),'r','zz_ww','sub_input_data_symphonie')

    zz_ww(:,:,:,:)=0._rprec

    ALLOCATE(depth_max(dims_reg(1,3),dims_reg(2,3)))
    CALL sub_memory(SIZE(depth_max),'r','depth_max','sub_input_data_symphonie')

    depth_max(:,:) = MAXVAL(ABS(zw0(:,:,1:kmax,1)), dim = 3)

    DO l = 1, dims_reg(4,3)
      DO k = 1, kmax
        DO j = 1, dims_reg(2,3)
          DO i = 1, dims_reg(1,3)
            IF (depth_max(i,j) /= 0._RPREC) THEN
              zz_ww(i,j,k,l) = sses(i,j,1,l) + zw0(i,j,k,1) * &
                   (1 + ( sses(i,j,1,l)/depth_max(i,j) ) )
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    CALL sub_memory(-SIZE(depth_max),'r','depth_max','sub_input_data_symphonie')
    DEALLOCATE(depth_max)

    WRITE(lun_standard,*)' - ZZ_WW : max ', MAXVAL(zz_ww), ' min ', MINVAL(zz_ww)

    !=================== E3T (3D + time) ========================!
    !- Compute E3T from ZZ_WW -!
    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' ------ E3T scale factor is computing ------'
    ALLOCATE(e3t(dims_reg(1,3), dims_reg(2,3), dims_reg(3,3), dims_reg(4,3)))
    CALL sub_memory(SIZE(e3t),'r','e3t','sub_input_data_symphonie')

    e3t(:,:,1:dims_reg(3,3)-1,:)= &
         zz_ww(:,:,1:dims_reg(3,3)-1,:) - zz_ww(:,:,2:dims_reg(3,3),:)

    e3t(:,:,dims_reg(3,3),:) = e3t(:,:,dims_reg(3,3)-1,:)

    e3t(:,:,:,:) = e3t(:,:,:,:) * tmask(:,:,:,:)

    WRITE(lun_standard,*)' - E3T : max ', MAXVAL(e3t), ' min ', MINVAL(e3t)

    !========================= E3U ==============================!
    ALLOCATE(e3u(dims_reg(1,3), dims_reg(2,3), dims_reg(3,3), dims_reg(4,3)))
    CALL sub_memory(SIZE(e3u),'r','e3u','sub_input_data_symphonie')

    e3u(1:dims_reg(1,3)-1,:,:,:) = &
         .5_rprec *( e3t(1:dims_reg(1,3)-1,:,:,:) + e3t(2:dims_reg(1,3),:,:,:) )

    e3u(dims_reg(1,3),:,:,:) = e3u(dims_reg(1,3)-1,:,:,:)

    !========================= E3V ==============================!
    ALLOCATE(e3v(dims_reg(1,3), dims_reg(2,3), dims_reg(3,3), dims_reg(4,3)))
    CALL sub_memory(SIZE(e3u),'r','e3u','sub_input_data_symphonie')

    e3v(:,1:dims_reg(2,3)-1,:,:) = &
         .5_rprec *( e3t(:,1:dims_reg(2,3)-1,:,:) + e3t(:,2:dims_reg(2,3),:,:) )

    e3v(:,dims_reg(2,3),:,:) = e3v(:,dims_reg(2,3)-1,:,:)

    !====================== TRANSPORT ===========================!
    !=== U ===!
    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' ------ Zonal transport is computing ------'
    DO l = 1, dims_reg(4,3)
      DO k = 1, dims_reg(3,3)
        DO j = 1, dims_reg(2,3)
          DO i = 1, dims_reg(1,3)
            uu(i,j,k,l) =  uu(i,j,k,l) * e3u(i,j,k,l) * e2u(i,j,1,1)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    WRITE(lun_standard,*)' - U Transport: max ', MAXVAL(uu), ' min ', MINVAL(uu)

    CALL sub_memory(-SIZE(e3u),'r','e3u','sub_input_data_symphonie')
    DEALLOCATE(e3u)


    !=== V ===!
    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' ------ Meridional transport is computing ------'
    DO l = 1, dims_reg(4,3)
      DO k = 1, dims_reg(3,3)
        DO j = 1, dims_reg(2,3)
          DO i = 1, dims_reg(1,3)
            vv(i,j,k,l) =  vv(i,j,k,l) * e3v(i,j,k,l) * e1v(i,j,1,1)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    WRITE(lun_standard,*)' - V Transport: max ', MAXVAL(vv), ' min ', MINVAL(vv)

    CALL sub_memory(-SIZE(e3v),'r','e3v','sub_input_data_symphonie')
    DEALLOCATE(e3v)


    !=== W ===!
    !- Must be no divergent -!
    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' ------ Vertical transport is computing ------'
    ww(:,:,dims_reg(3,3),:) = 0._rprec

    DO l = 1, dims_reg(4,3)
      DO k = dims_reg(3,3)-1, 1, -1
        DO j = 2, dims_reg(2,3)-1
          DO i = 2, dims_reg(1,3)-1
            ww(i,j,k,l) = ww(i,j,k+1,l) &
                 + uu(i-1,   j, k, l) - uu(i, j, k, l) &
                 + vv(  i, j-1, k, l) - vv(i, j, k, l)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO l = 1, dims_reg(4,3)
      DO k = 2, dims_reg(3,3)-1
        DO j = 1, dims_reg(2,3)
          DO i = 1, dims_reg(1,3)
            IF (zz_ww(i,j,1,l) - zz_ww(i,j,dims_reg(3,3),l) /= 0._rprec) THEN
              ww(i,j,k,l) = &
                   +   ww(i,j,k,l) &
                   - ( ww(i,j,1,l) / ( zz_ww(i,j,1,l) - zz_ww(i,j,dims_reg(3,3),l) ) ) &
                   *                 ( zz_ww(i,j,k,l) - zz_ww(i,j,dims_reg(3,3),l) )
            ELSE
              ww(i,j,k,l) = 0._rprec
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ww(:,:, 1,:) = 0._rprec
    ww(:,:,dims_reg(3,3),:) = 0._rprec

    WRITE(lun_standard,*)' - W Transport: max ', MAXVAL(ww), ' min ', MINVAL(ww)

    zz_ww(:,:,:,:) = zw0(:,:,:,:)

    CALL sub_memory(-SIZE(zw0),'r','zw0','sub_input_data_symphonie')
    DEALLOCATE(zw0)

    !! NG: call sub_memory(-size(sses),'r','sse','sub_input_data_symphonie')
    !! NG: DEALLOCATE(sse)


  END SUBROUTINE sub_input_data_symphonie
  !!***
  !=========================================================================
  !!****f* mod_input_data/sub_input_data_seq_symphonie()
  !! NAME
  !!   sub_input_data_seq_symphonie()
  !!
  !! FUNCTION
  !!   Data used in Ariane are on C grid and they are based on 
  !!   OPA conventions. SYMPHONIE data are on a C grid but don't respect
  !!   all OPA conventions.
  !!   The vertical axis is reversed.
  !!   U is defined only on (imt-1, jmt  , kmt-1)
  !!   V is defined only on (imt  , jmt-1, kmt-1)
  !!   W must be computed
  !!   T and S are defined on (imt, jmt, kmt-1)
  !!
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
  !!   * no arguments
  !!
  !! TODO
  !!   
  !!
  !! USED BY
  !!   * trajec
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_input_data_seq_symphonie( &
       ncids                             , &
       varids                            , &
       new_file                          , &
       ind_file                          , &
       ind_time                          , &
       ind_time_size                     , &
       dimsorders                        )

    INTEGER(kind = iprec) , DIMENSION(:), INTENT(inout) :: ncids
    INTEGER(kind = iprec) , DIMENSION(:), INTENT(inout) :: varids
    LOGICAL               , DIMENSION(:), INTENT(inout) :: new_file
    INTEGER(kind = iprec) , DIMENSION(:), INTENT(inout) :: ind_file
    INTEGER(kind = iprec) , DIMENSION(:), INTENT(inout) :: ind_time
    INTEGER(kind = iprec) , DIMENSION(:), INTENT(inout) :: ind_time_size
    INTEGER(kind = iprec) , DIMENSION(:,:), INTENT(inout) :: dimsorders

    INTEGER(kind=iprec) :: i,j,k,l
    CHARACTER(len=4), PARAMETER :: c_none='NONE'
    CHARACTER(len = 128)  :: c_filename

    INTEGER(kind = iprec) :: kmax

    REAL(kind=rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: &
         zw0, &
         e3u  , &
         e3v

    REAL(kind=rprec), DIMENSION(:,:), ALLOCATABLE :: depth_max
    INTEGER(kind = iprec) , DIMENSION(:,:), ALLOCATABLE :: ind_depth_max
    REAL(kind=rprec):: miss_val

    INTEGER(kind = iprec) :: ncid, varid
    INTEGER(kind = iprec) :: dimx  ! dimension in x (i)
    INTEGER(kind = iprec) :: dimy  ! dimension in y (j)
    INTEGER(kind = iprec) :: dimz  ! dimension in z (k)
    INTEGER(kind = iprec) :: dimt  ! dimension in t (l)
    INTEGER(kind = iprec), DIMENSION(4) :: dimsorder

    kmax =  dims_reg(3,3)

    !=== ZZ_TT has no time dimension = problem !!!! ====!!

    CALL sub_build_filename( &
         ind_file(1), maxsize_zt,c_prefix_zt, c_suffix_zt, c_filename)
    CALL sub_open_netcdf_file(c_dir_zt, c_filename, ncid)
    CALL sub_select_var_dims(ncid, nc_var_zt, varid, &
         dimsorder, dimx, dimy, dimz, dimt)
    CALL sub_read_netcdf_var4d(ncid, varid, zz_tt(:,:,dims_reg(3,3)-1:1:-1,:), dimsorder, dims_reg)
    zz_tt(:,:,:,:) =  zz_tt(:,:,:,:) * tmask(:,:,:,:)

    WRITE(lun_standard,*)' -  ', TRIM(nc_var_depth_t_symp), &
         ' max ',MAXVAL(zz_tt), ' min ', MINVAL(zz_tt)

    !======================= ZONAL CURRENT ==================!
    !- Read U on SYMPHONIE grid -!
    CALL sub_input_data_seq_main(                         &
         ncids(1), varids(1)                            , &
         new_file(1), ind_file(1)                       , &
         ind_time(1), ind_time_size(1)                  , &
         dimsorders(:,1)                                , &
         c_dir_zo,c_prefix_zo, maxsize_zo, c_suffix_zo  , &
         nc_var_zo,nc_var_eivu,nc_att_mask_zo           , &
         uu(:,:,dims_reg(3,3)-1:1:-1,1:1))

    uu(:,:,:,:) = CSHIFT(uu(:,:,:,:), shift=1, dim=1)
    uu(dims_reg(1,3),:,:,:) = 0._rprec

    !============= MERIDIONAL CURRENT =======================!
    !- Read V on SYMPHONIE grid -!
    CALL sub_input_data_seq_main(                         &
         ncids(2), varids(2)                            , &
         new_file(2), ind_file(2)                       , &
         ind_time(2), ind_time_size(2)                  , &
         dimsorders(:,2)                                , &
         c_dir_me,c_prefix_me,maxsize_me, c_suffix_me   , &
         nc_var_me,nc_var_eivv,nc_att_mask_me           , &
         vv(:,:,dims_reg(3,3)-1:1:-1,1:1))

    vv(:,:,:,:) = CSHIFT(vv(:,:,:,:), shift=1, dim=2)
    vv(:,dims_reg(2,3),:,:) = 0._rprec

    IF (key_alltracers) THEN
      !==================== TEMPERATURE =========================!
      !- Read T on SYMPHONIE grid -!
      CALL sub_input_data_seq_main(                       &
           ncids(4), varids(4)                          , &
           new_file(4), ind_file(4)                     , &
           ind_time(4), ind_time_size(4)                , &
           dimsorders(:,4)                              , &
           c_dir_te,c_prefix_te,maxsize_te, c_suffix_te , &
           nc_var_te,c_none,nc_att_mask_te              , &
           tt(:,:,dims_reg(3,3)-1:1:-1,1:1))

      !- Read mask value in netcdf file -!
      CALL sub_read_netcdf_att_val(ncids(4), varids(4), &
           TRIM(nc_att_mask_te), miss_val,0._rprec)

      !====================== SALINITY ===========================!
      !- Read S on SYMPHONIE grid -!
      CALL sub_input_data_seq_main(                       &
           ncids(5), varids(5)                          , &
           new_file(5), ind_file(5)                     , &
           ind_time(5), ind_time_size(5)                , &
           dimsorders(:,5)                              , &
           c_dir_sa,c_prefix_sa, maxsize_sa, c_suffix_sa, &
           nc_var_sa,c_none,nc_att_mask_sa              , &
           ss(:,:,dims_reg(3,3)-1:1:-1,1:1))


      IF (key_computesigma) THEN
        !====================== COMPUTE RR =========================!
        WRITE(lun_standard,*)''
        WRITE(lun_standard,*)'------ Density is computing -------'
        !-- Compute Density from Temperature and Salinity --!
        DO l = 1, 1
          DO k = 1,  dims_reg(3,3)-1
            DO j = 1,  dims_reg(2,3)
              DO i = 1,  dims_reg(1,3)
                rr(i,j,k,l)=rhostp(1000._rprec,ss(i,j,k,l),tt(i,j,k,l),-ABS(zsigma))
              END DO
            END DO
          END DO
        ENDDO

      ELSE

        !-- Read Density --!
        CALL sub_input_data_seq_main(                         &
             ncids(6), varids(6)                            , &
             new_file(6), ind_file(6)                       , &
             ind_time(6), ind_time_size(6)                  , &
             dimsorders(:,6)                                , &
             c_dir_de,c_prefix_de,maxsize_de, c_suffix_de   , &
             nc_var_de,c_none,nc_att_mask_de                , &
             rr(:,:,dims_reg(3,3)-1:1:-1,1:1))

      ENDIF

      !- Comments -!
      IF (id_comments) THEN
        WRITE(lun_standard,*)' - Density: max ', MAXVAL(rr), ' min ', MINVAL(rr)
      ENDIF
    ENDIF

    !====================== READ ZZ_TT =======================!
    !- Read ZZ_TT on SYMPHONIE grid -!
    !- surface is at indice kmax and deep at 0.

    !! Problem here because the ZZ_TT hasno time dimension in the output file
!!$    CALL sub_input_data_seq_main(                               &
!!$         ncids(8), varids(8)                                  , &
!!$         new_file(8), ind_file(8)                             , &
!!$         ind_time(8), ind_time_size(8)                        , &
!!$         dimsorders(:,8)                                      , &
!!$         c_dir_zt , c_prefix_zt, maxsize_zt, c_suffix_zt      , &
!!$         nc_var_zt, c_none, nc_att_mask_zt                    , &
!!$         zz_tt(:,:,dims_reg(3,3)-1:1:-1,:))
!!$
!!$    zz_tt(:,:,:,:) =  zz_tt(:,:,:,:) * tmask(:,:,:,:)

    ALLOCATE(zw0(dims_reg(1,3),dims_reg(2,3),dims_reg(3,3), 1))
    CALL sub_memory(SIZE(zw0),'r','zw0','sub_input_data_seq_symphonie')
    zw0(:,:,1,:) = 0._rprec

    DO k = 2, kmax
      DO j = 1, dims_reg(2,3)
        DO i = 1, dims_reg(1,3)
          zw0(i,j,k,1) = 2._rprec * zz_tt(i,j,k-1,1) - zw0(i,j,k-1,1)
        ENDDO
      ENDDO
    ENDDO

    zw0(:,:,:,:) =  zw0(:,:,:,:) * wmask(:,:,:,:)

    !====================== READ SSE 2D + Time =======================!
    !- Read SSE on SYMPHONIE grid -!
    !!  NG: ALLOCATE(sse(dims_reg(1,3),dims_reg(2,3), 1, 1))
    !!  NG: call sub_memory(size(sse),'r','sse','sub_input_data_seq_symphonie')

    CALL sub_input_data_seq_main(                               &
         ncids(7), varids(7)                                  , &
         new_file(7), ind_file(7)                             , &
         ind_time(7), ind_time_size(7)                        , &
         dimsorders(:,7)                                      , &
         c_dir_sse , c_prefix_sse, maxsize_sse, c_suffix_sse  , &
         nc_var_sse, c_none,nc_att_mask_sse, &
         sses(:,:,:,:))

    !==================== COMPUTE ZZ_WW (3D + time) =====================!
    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' ------ ZZ_WW is computing ------'

    zz_ww(:,:,:,1) = 0._rprec

    ALLOCATE(depth_max(dims_reg(1,3),dims_reg(2,3)))
    CALL sub_memory(SIZE(depth_max),'r','depth_max','sub_input_data_seq_symphonie')

    depth_max(:,:) = MAXVAL(ABS(zw0(:,:,1:kmax,1)), dim = 3)

    !!NG !! Save bathymetry for Symphonie data.
    !!NG  IF (ind_time(1) == 2) THEN
    !!NG    DO j = 1, dims_reg(2,3)
    !!NG      WRITE(22,'(340(1x,f0.2))') (depth_max(i,j), i= 1, dims_reg(1,3))
    !!NG    ENDDO
    !!NG  ENDIF

    ALLOCATE(ind_depth_max(dims_reg(1,3),dims_reg(2,3)))
    CALL sub_memory(SIZE(ind_depth_max),'i','ind_depth_max','sub_input_data_seq_symphonie')

    ind_depth_max(:,:) = MAXLOC(ABS(zw0(:,:,1:kmax,1)), dim = 3)

    DO k = 1, kmax
      DO j = 1, dims_reg(2,3)
        DO i = 1, dims_reg(1,3)
          IF (depth_max(i,j) /= 0._rprec) THEN
            zz_ww(i,j,k,1) = sses(i,j,1,1) + zw0(i,j,k,1) * &
                 (1 + ( sses(i,j,1,1)/depth_max(i,j) ) )
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    WRITE(lun_standard,*)' - ZZ_WW : max ', MAXVAL(zz_ww), ' min ', MINVAL(zz_ww)

    !! NG: write(*,*)' depth_max(158+1,70+1) :', depth_max(158+1,70+1)
    !! NG: write(*,*)' sse(158+1,70+1,:,1) :',   sse(158+1,70+1,:,1)
    !! NG: write(*,*)' zw0(158+1,70+1,:,1) :',   zw0(158+1,70+1,:,1)
    !! NG: write(*,*)' zz_ww(158+1,70+1,:,1) :', zz_ww(158+1,70+1,:,1)


    !=================== E3T (3D + time) ========================!
    !- Compute E3T from ZZ_WW -!
    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' ------ E3T scale factor is computing ------'
    IF (.NOT.ALLOCATED(e3t)) THEN
      ALLOCATE(e3t(dims_reg(1,3), dims_reg(2,3), dims_reg(3,3), 1))
      CALL sub_memory(SIZE(e3t),'r','e3t','sub_input_data_seq_symphonie')
    ENDIF

    e3t(:,:,1:dims_reg(3,3)-1,:)= &
         zz_ww(:,:,1:dims_reg(3,3)-1,:) - zz_ww(:,:,2:dims_reg(3,3),:)


    e3t(:,:,:,:) = e3t(:,:,:,:) * tmask(:,:,:,:)

    !!  WRITE(0,*)MINLOC(e3t(:,:,:,:), mask = e3t(:,:,:,:) < 0._rprec)

    WRITE(lun_standard,*)' - E3T : max ', MAXVAL(e3t), ' min ', MINVAL(e3t)

    !========================= E3U ==============================!
    IF (.NOT.ALLOCATED(e3u)) THEN
      ALLOCATE(e3u(dims_reg(1,3), dims_reg(2,3), dims_reg(3,3), 1))
      CALL sub_memory(SIZE(e3u),'r','e3u','sub_input_data_seq_symphonie')
    ENDIF
    e3u(1:dims_reg(1,3)-1,:,:,:) = &
         .5_rprec *( e3t(1:dims_reg(1,3)-1,:,:,:) + e3t(2:dims_reg(1,3),:,:,:) )

    e3u(dims_reg(1,3),:,:,:) = e3u(dims_reg(1,3)-1,:,:,:)

    !========================= E3V ==============================!
    IF (.NOT.ALLOCATED(e3v)) THEN
      ALLOCATE(e3v(dims_reg(1,3), dims_reg(2,3), dims_reg(3,3), 1))
      CALL sub_memory(SIZE(e3v),'r','e3v','sub_input_data_seq_symphonie')
    ENDIF
    e3v(:,1:dims_reg(2,3)-1,:,:) = &
         .5_rprec *( e3t(:,1:dims_reg(2,3)-1,:,:) + e3t(:,2:dims_reg(2,3),:,:) )

    e3v(:,dims_reg(2,3),:,:) = e3v(:,dims_reg(2,3)-1,:,:)

    !====================== TRANSPORT ===========================!
    !=== U ===!
    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' ------ Zonal transport is computing ------'
    DO l = 1, 1
      DO k = 1, dims_reg(3,3)
        DO j = 1, dims_reg(2,3)
          DO i = 1, dims_reg(1,3)
            uu(i,j,k,l) =  uu(i,j,k,l) * e3u(i,j,k,l) * e2u(i,j,1,1)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    WRITE(lun_standard,*)' - U Transport: max ', MAXVAL(uu), ' min ', MINVAL(uu)

    CALL sub_memory(-SIZE(e3u),'r','e3u','sub_input_data_seq_symphonie')
    DEALLOCATE(e3u)


    !=== V ===!
    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' ------ Meridional transport is computing ------'
    DO l = 1, 1
      DO k = 1, dims_reg(3,3)
        DO j = 1, dims_reg(2,3)
          DO i = 1, dims_reg(1,3)
            vv(i,j,k,l) =  vv(i,j,k,l) * e3v(i,j,k,l) * e1v(i,j,1,1)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    WRITE(lun_standard,*)' - V Transport: max ', MAXVAL(vv), ' min ', MINVAL(vv)

    CALL sub_memory(-SIZE(e3v),'r','e3v','sub_input_data_seq_symphonie')
    DEALLOCATE(e3v)


    !=== W ===!
    !- Must be no divergent -!
    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' ------ Vertical transport is computing ------'
    ww(:,:,dims_reg(3,3),:) = 0._rprec

    DO l = 1,1
      DO k = dims_reg(3,3)-1, 1, -1
        DO j = 2, dims_reg(2,3)-1 
          DO i = 2, dims_reg(1,3)-1
            ww(i,j,k,l) = ww(i,j,k+1,l) &
                 + uu(i-1,   j, k, l) - uu(i, j, k, l) &
                 + vv(  i, j-1, k, l) - vv(i, j, k, l)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO l = 1, 1
      DO j = 1, dims_reg(2,3)
        DO i = 1, dims_reg(1,3)
          DO k = 2, ind_depth_max(i,j) - 1
            ww(i,j,k,l) = &
                 +   ww(i,j,k,l) &
                 - ( ww(i,j,1,l) / ( zz_ww(i,j,1,l) - depth_max(i,j) ) ) &
                 *                 ( zz_ww(i,j,k,l) - depth_max(i,j) )
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ww(:,:, 1,:) = 0._rprec
    ww(:,:,dims_reg(3,3),:) = 0._rprec

    WRITE(lun_standard,*)' - W Transport: max ', MAXVAL(ww), ' min ', MINVAL(ww)

    zz_ww(:,:,:,:) = zw0(:,:,:,:)

    CALL sub_memory(-SIZE(depth_max),'r','depth_max','sub_input_data_seq_symphonie')
    DEALLOCATE(depth_max)
    CALL sub_memory(-SIZE(ind_depth_max),'i','ind_depth_max','sub_input_data_seq_symphonie')
    DEALLOCATE(ind_depth_max)
    CALL sub_memory(-SIZE(zw0),'r','zw0','sub_input_data_seq_symphonie')
    DEALLOCATE(zw0)
    !! NG: call sub_memory(-size(sses),'r','sse','sub_input_data_seq_symphonie')
    !! NG: DEALLOCATE(sses)


  END SUBROUTINE sub_input_data_seq_symphonie
  !!***


  !!****f* mod_input_data/sub_transp_alloc()
  !! NAME
  !!   sub_transp_alloc()
  !!
  !! FUNCTION
  !!   Dynamic allocation of the transport arrays (uu, vv, ww)
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
  !!   * dim1: dimension 1, or x, or i
  !!   * dim1: dimension 2, or y, or j
  !!   * dim1: dimension 3, or z, or k
  !!   * dim1: dimension 4, or t, or l
  !!
  !! TODO
  !!   
  !!
  !! USED BY
  !!   * mod_input_data (Private)
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_transp_alloc(dim1, dim2, dim3, dim4)

    !-------------!
    ! Declaration !
    !-------------!
    !- arguments -!
    INTEGER(kind=iprec), INTENT(in) :: dim1, dim2, dim3, dim4

    !-------------!
    ! Code begins !
    !-------------!
    !- Dynamic allocation -!
    ALLOCATE(uu(dim1, dim2, dim3, dim4))
    CALL sub_memory(SIZE(uu),'r','uu','sub_transp_alloc')
    ALLOCATE(vv(dim1, dim2, dim3, dim4))
    CALL sub_memory(SIZE(vv),'r','vv','sub_transp_alloc')
    IF (key_mars) THEN
      ALLOCATE(ww(dim1, dim2, dim3+1, dim4))
    ELSE
      ALLOCATE(ww(dim1, dim2, dim3, dim4))
    ENDIF
    CALL sub_memory(SIZE(ww),'r','ww','sub_transp_alloc')

    IF (TRIM(w_surf_option) == 'E-P-R') THEN
      ALLOCATE(epr(dim1, dim2, 1, dim4))
      CALL sub_memory(SIZE(epr),'r','epr','sub_transp_alloc')
    ENDIF

    IF (key_symphonie) THEN
      ALLOCATE(zz_tt(dim1, dim2, dim3, dim4))
      CALL sub_memory(SIZE(zz_tt),'r','zz_tt','sub_transp_alloc')
      ALLOCATE(sses(dim1, dim2, dim3, dim4))
      CALL sub_memory(SIZE(sses),'r','sse','sub_transp_alloc')
    ENDIF

  END SUBROUTINE sub_transp_alloc
  !!***
  !!****f* mod_input_data/sub_transp_dealloc()
  !! NAME
  !!   sub_transp_dealloc()
  !!
  !! FUNCTION
  !!   Deallocate memory of the transport arrays (uu, vv, ww)
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
  SUBROUTINE sub_transp_dealloc()

    !-------------!
    ! Code begins !
    !-------------!
    !- Deallocate arrays memory -!
    IF (ALLOCATED(uu)) THEN
      CALL sub_memory(-SIZE(uu),'r','uu','sub_transp_dealloc')
      DEALLOCATE(uu)
    ENDIF
    IF (ALLOCATED(vv)) THEN
      CALL sub_memory(-SIZE(vv),'r','vv','sub_transp_dealloc')
      DEALLOCATE(vv)
    ENDIF
    IF (ALLOCATED(ww)) THEN
      CALL sub_memory(-SIZE(ww),'r','ww','sub_transp_dealloc')
      DEALLOCATE(ww)
    ENDIF

    IF (ALLOCATED(epr)) THEN
      CALL sub_memory(-SIZE(epr),'r','epr','sub_transp_dealloc')
      DEALLOCATE(epr)
    ENDIF

    IF (ALLOCATED(zz_tt)) THEN
      CALL sub_memory(-SIZE(zz_tt),'r','zz_tt','sub_transp_dealloc')
      DEALLOCATE(zz_tt)
      CALL sub_memory(SIZE(sses),'r','sse','sub_transp_alloc')
      DEALLOCATE(sses)
    ENDIF


  END SUBROUTINE sub_transp_dealloc
  !!***
  !=========================================================================
  !!****f* mod_input_data/sub_tracer_alloc()
  !! NAME
  !!   sub_tracer_alloc()
  !!
  !! FUNCTION
  !!   Dynamic allocation of the tracer arrays (tt, ss, rr)
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
  !!   * dim1: dimension 1, or x, or i
  !!   * dim1: dimension 2, or y, or j
  !!   * dim1: dimension 3, or z, or k
  !!   * dim1: dimension 4, or t, or l
  !!
  !! TODO
  !!   
  !!
  !! USED BY
  !!   * mod_input_data (Private)
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_tracer_alloc(dim1, dim2, dim3, dim4)

    !-------------!
    ! Declaration !
    !-------------!
    !- arguments -!
    INTEGER(kind=iprec), INTENT(in) :: dim1, dim2, dim3, dim4

    !-------------!
    ! Code begins !
    !-------------!
    !- Dynamic allocation -!
    ALLOCATE(tt(dim1, dim2, dim3, dim4))
    CALL sub_memory(SIZE(tt),'r','tt','sub_tracer_alloc')

    ALLOCATE(ss(dim1, dim2, dim3, dim4))
    CALL sub_memory(SIZE(ss),'r','ss','sub_tracer_alloc')

    ALLOCATE(rr(dim1, dim2, dim3, dim4))
    CALL sub_memory(SIZE(rr),'r','rr','sub_tracer_alloc')

  END SUBROUTINE sub_tracer_alloc
  !!***

  !!****f* mod_input_data/sub_tracer_dealloc()
  !! NAME
  !!   sub_tracer_dealloc()
  !!
  !! FUNCTION
  !!   Deallocate memory of the tracer arrays (tt, ss, rr)
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
  SUBROUTINE sub_tracer_dealloc()

    !-------------!
    ! Code begins !
    !-------------!
    !- Deallocate arrays -!
    IF (ALLOCATED(tt)) THEN
      CALL sub_memory(-SIZE(tt),'r','tt','sub_tracer_dealloc')
      DEALLOCATE(tt)
    ENDIF
    IF (ALLOCATED(ss)) THEN
      CALL sub_memory(-SIZE(ss),'r','ss','sub_tracer_dealloc')
      DEALLOCATE(ss)
    ENDIF
    IF (ALLOCATED(rr)) THEN
      CALL sub_memory(-SIZE(rr),'r','rr','sub_tracer_dealloc')
      DEALLOCATE(rr)
    ENDIF

  END SUBROUTINE sub_tracer_dealloc


  !=========================================================================
  !!****f* mod_input_data/sub_ssh_alloc()
  !! NAME
  !!   sub_ssh_alloc()
  !!
  !! FUNCTION
  !!   Dynamic allocation of the ssh arrays
  !!
  !! AUTHOR
  !!   * Origin  : Susan Allen Jan 2018, adapted from sub_tracers_alloc()
  !!
  !! CREATION DATE
  !!   * Jan 2018
  !!
  !! HISTORY
  !!   Date (dd/mm/yyyy/) - Modification(s)
  !!
  !! ARGUMENTS
  !!   * dim1: dimension 1, or x, or i
  !!   * dim1: dimension 2, or y, or j
  !!
  !!
  !! TODO
  !!
  !!
  !! USED BY
  !!   * mod_input_data (Private)
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_ssh_alloc(dim1, dim2)

    !-------------!
    ! Declaration !
    !-------------!
    !- arguments -!
    INTEGER(kind=iprec), INTENT(in) :: dim1, dim2

    !-------------!
    ! Code begins !
    !-------------!
    !- Dynamic allocation -!
    ALLOCATE(ssh_vvl(dim1, dim2, 1, 1))
    CALL sub_memory(size(ssh_vvl),'r','ssh','sub_ssh_alloc')

  END SUBROUTINE sub_ssh_alloc
  !!***

  !!****f* mod_input_data/sub_ssh_dealloc()
  !! NAME
  !!   sub_ssh_dealloc()
  !!
  !! FUNCTION
  !!   Deallocate memory of the ssh array
  !!
  !! AUTHOR
  !!   * Origin  : Susan Allen Jan 2018 adapted from sub_tracers_dealloc()
  !!
  !! CREATION DATE
  !!   * Jan 2018
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
  !!   * trajec_seq
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_ssh_dealloc()

    !-------------!
    ! Code begins !
    !-------------!
    !- Deallocate arrays -!
    IF (ALLOCATED(ssh_vvl)) THEN
       CALL sub_memory(-size(ssh_vvl),'r','ssh','sub_ssh_dealloc')
       DEALLOCATE(ssh_vvl)
    ENDIF

  END SUBROUTINE sub_ssh_dealloc


  !!***
  !!****f* mod_input_data/sub_tracer_dealloc()
  !! NAME
  !!   sub_tracer_dealloc()
  !!
  !! FUNCTION
  !!   Deallocate memory of the tracer arrays (tt, ss, rr)
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
  SUBROUTINE sub_input_data_dealloc_mem()

    IF (ALLOCATED(zeta_a)) THEN
      CALL sub_memory(-SIZE(zeta_a),'r','zeta_a','sub_input_data_dealloc_mem')
      DEALLOCATE(zeta_a)
    END IF

    IF (ALLOCATED(zw0)) THEN
      CALL sub_memory(-SIZE(zw0),'r','zw0','sub_input_data_dealloc_mem')
      DEALLOCATE(zw0)
    END IF

    IF (ALLOCATED(e3u)) THEN
      CALL sub_memory(-SIZE(e3u),'r','e3u','sub_input_data_dealloc_mem')
      DEALLOCATE(e3u)
    END IF

    IF (ALLOCATED(e3v)) THEN
      CALL sub_memory(-SIZE(e3v),'r','e3v','sub_input_data_dealloc_mem')
      DEALLOCATE(e3v)
    END IF

    IF (ALLOCATED(tmp_array)) THEN
      CALL sub_memory(-SIZE(tmp_array),'r','tmp_array','sub_input_data_dealloc_mem')
      DEALLOCATE(tmp_array)
    END IF

  END SUBROUTINE sub_input_data_dealloc_mem


END MODULE mod_input_data
