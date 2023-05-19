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
MODULE mod_seq

  !------------------!
  ! USE ASSOCIAITION !
  !------------------!
  USE mod_precision
  USE mod_memory
  USE mod_namelist
  USE mod_netcdf

  !-------------!
  ! DECLARATION !
  !-------------!
  IMPLICIT NONE

  LOGICAL, DIMENSION(:), ALLOCATABLE :: new_file
  !
  INTEGER (kind = iprec), DIMENSION(:), ALLOCATABLE :: &
       ncids        , & ! Netcdf file IDs
       varids       , & ! Netcdf file variable IDs
       ind_file     , & ! File indices
       ind_time     , & ! Time indices
       ind_time_size    ! Time size

  INTEGER (kind = iprec), DIMENSION(:,:), ALLOCATABLE :: &
       sdimsorders

CONTAINS

  !!=========================================================================
  SUBROUTINE sub_seq_alloc(alloc_size)

    INTEGER(kind = iprec), INTENT(in) :: alloc_size

    write (*,*) 'alloc_size is ', alloc_size

    IF (.NOT.ALLOCATED(ncids)) THEN
       ALLOCATE(ncids(alloc_size))
       CALL sub_memory(size(ncids),'i','ncids','sub_seq_alloc')
    ENDIF 
    IF (.NOT.ALLOCATED(varids))  THEN
       ALLOCATE(varids(alloc_size))
       CALL sub_memory(size(varids),'i','varids','sub_seq_alloc')
    ENDIF
    IF (.NOT.ALLOCATED(new_file))THEN
       ALLOCATE(new_file(alloc_size))
    ENDIF
    IF (.NOT.ALLOCATED(ind_file)) THEN
       ALLOCATE(ind_file(alloc_size))
       CALL sub_memory(size(ind_file),'i','ind_file','sub_seq_alloc')
    ENDIF
    IF (.NOT.ALLOCATED(ind_time)) THEN
       ALLOCATE(ind_time(alloc_size))
       CALL sub_memory(size(ind_time),'i','ind_time','sub_seq_alloc')
    ENDIF
    IF (.NOT.ALLOCATED(ind_time_size)) THEN
       ALLOCATE(ind_time_size(alloc_size))
       CALL sub_memory(size(ind_time_size),'i','ind_time_size','sub_seq_alloc')
    ENDIF
    IF (.NOT.ALLOCATED(sdimsorders  )) THEN
       ALLOCATE(sdimsorders(4,alloc_size)) 
       CALL sub_memory(size(sdimsorders),'i','sdimsorders','sub_seq_alloc')
    ENDIF

  END SUBROUTINE sub_seq_alloc

  !!=========================================================================
  SUBROUTINE sub_seq_dealloc()


    IF (ALLOCATED(ncids)) THEN
       CALL sub_memory(-size(ncids),'i','ncids','sub_seq_dealloc')
       DEALLOCATE(ncids)
    ENDIF
    IF (ALLOCATED(varids)) THEN
       CALL sub_memory(-size(varids),'i','varids','sub_seq_dealloc')
       DEALLOCATE(varids)
    ENDIF
    IF (ALLOCATED(new_file)) THEN
       DEALLOCATE(new_file)
    ENDIF
    IF (ALLOCATED(ind_file)) THEN
       CALL sub_memory(-size(ind_file),'i','ind_file','sub_seq_dealloc')
       DEALLOCATE(ind_file)
    ENDIF
    IF (ALLOCATED(ind_time)) THEN
       CALL sub_memory(-size(ind_time),'i','ind_time','sub_seq_dealloc')
       DEALLOCATE(ind_time)
    ENDIF
    IF (ALLOCATED(ind_time_size)) THEN
       CALL sub_memory(-size(ind_time_size),'i','ind_time_size','sub_seq_dealloc')
       DEALLOCATE(ind_time_size)
    ENDIF
    IF (ALLOCATED(sdimsorders)) THEN
       CALL sub_memory(-size(sdimsorders),'i','sdimsorders','sub_seq_dealloc')
       DEALLOCATE(sdimsorders)
    ENDIF 

  END SUBROUTINE sub_seq_dealloc

  !!=========================================================================
  SUBROUTINE sub_seq_init(i_time)

    INTEGER(kind = iprec), OPTIONAL, INTENT(in) :: i_time

    INTEGER(kind = iprec) :: &
         ind_zo, & !
         ind_me, & !
         ind_ve, & !
         ind_te, & !
         ind_sa, & !
         ind_rr, & !
         ind_se, & !
         ind_ssh, & !
         ind_ep

    INTEGER(kind = iprec), DIMENSION(:), ALLOCATABLE :: &
         ind_tt

    ncids(:)         = 0
    varids(:)        = 0
    ind_file(:)      = 0
    ind_time(:)      = 0
    new_file(:)      = .TRUE.
    ind_time_size(:) = 0
    sdimsorders(:,:) = 0

    IF (PRESENT(i_time)) THEN


       IF (key_alltracers) THEN

          IF (key_roms.or.key_mars) THEN

             ALLOCATE(ind_tt(7))
             CALL sub_memory(size(ind_tt),'i','ind_tt','sub_seq_init')

             CALL sub_search_forward_inds(ind_min=i_time,             &
                  ind_zo = ind_zo , ind_me = ind_me, ind_ve = ind_ve, & 
                  ind_te = ind_te, ind_sa = ind_sa,  ind_rr = ind_rr, &
                  ind_se = ind_se, ind_tt = ind_tt)

             ind_file(4) = ind_te
             ind_file(5) = ind_sa
             ind_file(6) = ind_rr
             ind_file(7) = ind_se  ! ind_se = ind_ze

          ELSEIF (key_symphonie) THEN
             ALLOCATE(ind_tt(8))
             CALL sub_memory(size(ind_tt),'i','ind_tt','sub_seq_init')

             CALL sub_search_forward_inds(ind_min=i_time,             &
                  ind_zo = ind_zo , ind_me = ind_me, ind_ve = ind_ve, &
                  ind_te = ind_te, ind_sa = ind_sa,  ind_rr = ind_rr, &
                  ind_se = ind_se, ind_tt = ind_tt)

             ind_file(4) = ind_te
             ind_file(5) = ind_sa
             ind_file(6) = ind_rr
             ind_file(7) = ind_se
             ind_file(8) = ind_te
             ind_tt(8) = 1   !! to include in sub_search_forward_inds

          ELSE
             IF (TRIM(w_surf_option) == 'E-P-R') THEN
                ALLOCATE(ind_tt(7))
                CALL sub_memory(size(ind_tt),'i','ind_tt','sub_seq_init')

                CALL sub_search_forward_inds(ind_min=i_time,             &
                     ind_zo = ind_zo , ind_me = ind_me, ind_ve = ind_ve, &
                     ind_te = ind_te, ind_sa = ind_sa,  ind_rr = ind_rr, &
                     ind_tt = ind_tt, ind_ep = ind_ep)

                ind_file(7) = ind_ep
             ELSEIF (key_vvl) THEN
                ALLOCATE(ind_tt(7))
                CALL sub_memory(size(ind_tt),'i','ind_tt','sub_seq_init')

                CALL sub_search_forward_inds(ind_min=i_time,             &
                     ind_zo = ind_zo , ind_me = ind_me, ind_ve = ind_ve, &
                     ind_te = ind_te, ind_sa = ind_sa,  ind_rr = ind_rr, &
                     ind_tt = ind_tt, ind_ssh = ind_ssh)

                ind_file(7) = ind_ssh
             ELSE
                ALLOCATE(ind_tt(6))
                CALL sub_memory(size(ind_tt),'i','ind_tt','sub_seq_init')

                CALL sub_search_forward_inds(ind_min=i_time,             &
                     ind_zo = ind_zo , ind_me = ind_me, ind_ve = ind_ve, &
                     ind_te = ind_te, ind_sa = ind_sa,  ind_rr = ind_rr, &
                     ind_tt = ind_tt)
             ENDIF
             ind_file(4) = ind_te
             ind_file(5) = ind_sa
             ind_file(6) = ind_rr
          ENDIF

          ind_file(1) = ind_zo
          ind_file(2) = ind_me
          ind_file(3) = ind_ve
          ind_time(:) = ind_tt(:)


       ELSE !! all_tracers

          !!NG: 25 july 2011 (modifications in ROMS, add se)
          IF (key_roms.or.key_mars) THEN
             ALLOCATE(ind_tt(3))
             CALL sub_memory(size(ind_tt),'i','ind_tt','sub_seq_init')

             CALL sub_search_forward_inds(ind_min=i_time,             &
                  ind_zo = ind_zo , ind_me = ind_me, ind_ve = ind_ve, &
                  ind_se = ind_se , ind_tt = ind_tt )

             ind_file(1) = ind_zo
             ind_file(2) = ind_me
             ind_file(3) = ind_se
             ind_time(:) = ind_tt(:)

          ELSEIF (key_symphonie) THEN
             ALLOCATE(ind_tt(4))
             CALL sub_memory(size(ind_tt),'i','ind_tt','sub_seq_init')

             CALL sub_search_forward_inds(ind_min=i_time            , &
                  ind_zo = ind_zo , ind_me = ind_me, ind_ve = ind_ve, &
                  ind_tt = ind_tt)

             ind_tt(4)=1 !! To include in sub_search_forward_inds feb 2010 !!

             ind_file(1) = ind_zo
             ind_file(2) = ind_me
             ind_file(3) = ind_ve
             ind_time(:) = ind_tt(:)

          ELSEIF (key_vvl) THEN
             ALLOCATE(ind_tt(4))
             CALL sub_memory(size(ind_tt),'i','ind_tt','sub_seq_init')

             CALL sub_search_forward_inds(ind_min=i_time,             &
                  ind_zo = ind_zo , ind_me = ind_me, ind_ve = ind_ve, &
                  ind_ssh = ind_ssh , ind_tt = ind_tt )

             ind_file(1) = ind_zo
             ind_file(2) = ind_me
             ind_file(3) = ind_ve
             ind_file(4) = ind_ssh
             ind_time(:) = ind_tt(:)

          ELSE
             ALLOCATE(ind_tt(3))
             CALL sub_memory(size(ind_tt),'i','ind_tt','sub_seq_init')

             CALL sub_search_forward_inds(ind_min=i_time            , &
                  ind_zo = ind_zo , ind_me = ind_me, ind_ve = ind_ve, &
                  ind_tt = ind_tt)

             ind_file(1) = ind_zo
             ind_file(2) = ind_me
             ind_file(3) = ind_ve
             ind_time(:) = ind_tt(:)

          ENDIF

          !!NG: 25 july 2011

       ENDIF



    ELSE !! i_time

       !! NG : debug  write(*,*)'(++) NO i_time : '

       IF (TRIM(forback) == 'forward' ) THEN

          !! NG : debug write(*,*)'(++) forward'

          ind_file(1) = ind0_zo
          ind_file(2) = ind0_me
          ind_file(3) = ind0_ve

          IF (key_alltracers) THEN
             ind_file(4) = ind0_te
             ind_file(5) = ind0_sa
             IF (.NOT.key_computesigma) THEN
                ind_file(6) = ind0_de
             ENDIF
             IF (TRIM(w_surf_option) == 'E-P-R') THEN
                ind_file(7) = ind0_ep
             ELSEIF (key_vvl) THEN
                ind_file(7) = ind0_ssh
             ENDIF
             IF (key_roms.or.key_mars) THEN
                ind_file(7) = ind0_ze
             ELSEIF(key_symphonie) THEN
                ind_file(7) = ind0_sse
                ind_file(8) = ind0_zt
             ENDIF
          ELSEIF(TRIM(w_surf_option) == 'E-P-R') THEN
             ind_file(4) = ind0_ep
          ELSEIF (key_vvl) THEN
             ind_file(4) = ind0_ssh
          ELSEIF (key_roms.or.key_mars) THEN
             ind_file(3) = ind0_ze
          ELSEIF(key_symphonie) THEN
             ind_file(4) = ind0_zt
          ENDIF
          ind_time(:) = 1

       ELSE

          ind_file(1) = indn_zo
          ind_file(2) = indn_me
          ind_file(3) = indn_ve
          IF (key_alltracers) THEN
             ind_file(4) = indn_te
             ind_file(5) = indn_sa
             IF (.NOT.key_computesigma) THEN
                ind_file(6) = indn_de
             ENDIF
             IF(TRIM(w_surf_option) == 'E-P-R') THEN
                ind_file(7) = indn_ep
             ELSEIF (key_vvl) THEN
                ind_file(7) = indn_ssh
             ENDIF
             IF (key_roms.or.key_mars) THEN
                ind_file(7) = indn_ze
             ELSEIF(key_symphonie) THEN
                ind_file(7) = indn_sse
                ind_file(8) = indn_zt
             ENDIF
          ELSEIF(TRIM(w_surf_option) == 'E-P-R') THEN
             ind_file(4) = indn_ep
          ELSEIF (key_vvl) THEN
             ind_file(4) = indn_ssh
          ELSEIF (key_roms.or.key_mars) THEN
             ind_file(3) = indn_ze
          ELSEIF(key_symphonie) THEN
             ind_file(4) = indn_zt
          ENDIF
          ind_time(:) = 0 !! Because we don't know how many records are available in the Netcdf File
       ENDIF

    ENDIF !! i_time

    IF (ALLOCATED(ind_tt)) THEN
       CALL sub_memory(-size(ind_tt),'i','ind_tt','sub_seq_init')
       DEALLOCATE(ind_tt)
    END IF

  END SUBROUTINE sub_seq_init

  !!=========================================================================
  SUBROUTINE sub_search_forward_inds(  &
       ind_min, &
       ind_zo, ind_me, ind_ve, &
       ind_te, ind_sa, ind_rr, &
       ind_se, ind_ssh, ind_tt, ind_ep)

    INTEGER(kind = iprec), INTENT(in)                :: &
         ind_min
    INTEGER(kind = iprec), INTENT(out)               :: &
         ind_zo, ind_me, ind_ve
    INTEGER(kind = iprec), OPTIONAL, INTENT(out)     :: &
         ind_te, ind_sa, ind_rr, ind_se, ind_ssh, ind_ep
    INTEGER(kind = iprec), DIMENSION(:), INTENT(out) :: &
         ind_tt

    INTEGER(kind = iprec) :: ind_w
    INTEGER(kind = iprec) :: tmp_ind

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)'====================================================='
    WRITE(lun_standard,*)'= Searching in Netcdf files the time limit indices  ='
    WRITE(lun_standard,*)'====================================================='

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)'- Zonal current file -'
    CALL sub_search_forward_ind(ind_min, ind_zo, ind_tt(1), &
         ind0_zo, c_dir_zo, c_prefix_zo, maxsize_zo, c_suffix_zo, &
         nc_var_zo)
    WRITE(lun_standard,*)'    - file name indice: ', ind_zo
    WRITE(lun_standard,*)'    - file time indice: ', ind_tt(1)

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)'- Meridional current file -'
    CALL sub_search_forward_ind(ind_min, ind_me, ind_tt(2), &
         ind0_me, c_dir_me, c_prefix_me, maxsize_me, c_suffix_me, &
         nc_var_me)
    WRITE(lun_standard,*)'    - file name indice: ', ind_me
    WRITE(lun_standard,*)'    - file time indice: ', ind_tt(2)

    IF ((key_computew).OR.(.NOT.key_read_w).OR.(key_roms).OR.(key_mars).OR.(key_symphonie)) THEN
       ind_ve    = iZero
       ind_tt(3) = iZero
    ELSE
       WRITE(lun_standard,*)''
       WRITE(lun_standard,*)'- Vertical current file -'
       CALL sub_search_forward_ind(ind_min, ind_ve, ind_tt(3)   , &
            ind0_ve, c_dir_ve,c_prefix_ve, maxsize_ve, c_suffix_ve, &
            nc_var_ve)
       WRITE(lun_standard,*)'    - file name indice: ', ind_ve
       WRITE(lun_standard,*)'    - file time indice: ', ind_tt(3)
    ENDIF

    IF (PRESENT(ind_te)) THEN
       WRITE(lun_standard,*)''
       WRITE(lun_standard,*)'- Temperature file -'
       CALL sub_search_forward_ind(ind_min, ind_te, ind_tt(4)   , &
            ind0_te, c_dir_te,c_prefix_te, maxsize_te, c_suffix_te, &
            nc_var_te)
       WRITE(lun_standard,*)'    - file name indice: ', ind_te
       WRITE(lun_standard,*)'    - file time indice: ', ind_tt(4)
    ENDIF

    IF (PRESENT(ind_sa)) THEN
       WRITE(lun_standard,*)''
       WRITE(lun_standard,*)'- Salinity file -'
       CALL sub_search_forward_ind(ind_min, ind_sa, ind_tt(5)   , &
            ind0_sa, c_dir_sa,c_prefix_sa, maxsize_sa, c_suffix_sa, &
            nc_var_sa)
       WRITE(lun_standard,*)'    - file name indice: ', ind_sa
       WRITE(lun_standard,*)'    - file time indice: ', ind_tt(5)
    ENDIF

    IF (PRESENT(ind_ep)) THEN
       WRITE(lun_standard,*)''
       WRITE(lun_standard,*)'- E-P-R file -'
       IF (key_alltracers) THEN
          ind_w=7
       ELSE
          ind_w=4
       ENDIF
       CALL sub_search_forward_ind(ind_min, ind_ep, ind_tt(ind_w) , &
            ind0_ep, c_dir_ep,c_prefix_ep, maxsize_ep, c_suffix_ep, &
            nc_var_ep)
       WRITE(lun_standard,*)'    - file name indice: ', ind_ep
       WRITE(lun_standard,*)'    - file time indice: ', ind_tt(ind_w)
    ENDIF

    IF (PRESENT(ind_ssh)) THEN
       WRITE(lun_standard,*)''
       WRITE(lun_standard,*)'- SSH file -'
       IF (key_alltracers) THEN
          ind_ssh=7
       ELSE
          ind_ssh=4
       ENDIF
       CALL sub_search_forward_ind(ind_min, ind_ssh, ind_tt(ind_ssh) , &
            ind0_ssh, c_dir_ssh,c_prefix_ssh, maxsize_ssh, c_suffix_ssh, &
            nc_var_ssh)
       WRITE(lun_standard,*)'    - file name indice: ', ind_ssh
       IF (key_alltracers) THEN
          WRITE(lun_standard,*)'    - file time indice: ', ind_tt(7)
       ELSE
          WRITE(lun_standard,*)'    - file time indice: ', ind_tt(4)
       ENDIF
    ENDIF
    

    IF (PRESENT(ind_rr)) THEN
       IF (key_computesigma) THEN
          ind_rr    = 0
          ind_tt(6) = 0
       ELSE
          WRITE(lun_standard,*)''
          WRITE(lun_standard,*)'- Density file-'
          CALL sub_search_forward_ind(ind_min, ind_rr, ind_tt(6)      , &
               ind0_de, c_dir_de, c_prefix_de, maxsize_de, c_suffix_de, &
               nc_var_de)
          WRITE(lun_standard,*)'    - file name indice: ', ind_rr
          WRITE(lun_standard,*)'    - file time indice: ', ind_tt(6)
       ENDIF
    ENDIF

    IF (PRESENT(ind_se)) THEN

       !!NG: 25 july 2011
       IF (key_alltracers) THEN 
          tmp_ind = iSeven
       ELSE
          tmp_ind = iThree
       END IF
       !!NG: 25 july 2011

       IF (key_roms.or.key_mars) THEN
          WRITE(lun_standard,*)''
          WRITE(lun_standard,*)'- Zeta file (roms)-'
          CALL sub_search_forward_ind(ind_min, ind_se, ind_tt(tmp_ind)      , &
               ind0_ze, c_dir_ze, c_prefix_ze, maxsize_ze, c_suffix_ze, &
               nc_var_ze)
       ELSEIF(key_symphonie) THEN
          WRITE(lun_standard,*)''
          WRITE(lun_standard,*)'- SSE file (symphonie)-'
          CALL sub_search_forward_ind(ind_min, ind_se, ind_tt(tmp_ind)           , &
               ind0_sse, c_dir_sse, c_prefix_sse, maxsize_sse, c_suffix_sse, &
               nc_var_sse)
       ENDIF
       WRITE(lun_standard,*)'    - file name indice: ', ind_se
       WRITE(lun_standard,*)'    - file time indice: ', ind_tt(tmp_ind)
    ENDIF

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)'================================'
    WRITE(lun_standard,*)'= Indice searches are finished ='
    WRITE(lun_standard,*)'================================'

  END SUBROUTINE sub_search_forward_inds

  !!=========================================================================
  SUBROUTINE sub_search_forward_ind(        &
       ind_min,                             &
       ind, ind_t,                          &
       indf, cdir, prefix, ndigits, suffix, &
       ncvar )

    INTEGER(kind = iprec), INTENT(in)  :: &
         ind_min
    INTEGER(kind = iprec), INTENT(out) :: &
         ind, & !
         ind_t  !
    INTEGER(kind = iprec), INTENT(in) :: &
         ndigits, indf
    CHARACTER(len = *)   , INTENT(in)  :: &
         cdir, prefix, suffix, ncvar

    INTEGER(kind = iprec) :: ncid
    INTEGER(kind = iprec) :: varid
    INTEGER(kind = iprec) :: dimx         ! dimension in x (i)
    INTEGER(kind = iprec) :: dimy         ! dimension in y (j)
    INTEGER(kind = iprec) :: dimz         ! dimension in z (k)
    INTEGER(kind = iprec) :: dimt         ! dimension in t (l)
    INTEGER(kind = iprec) :: nb_read
    INTEGER(kind = iprec), DIMENSION(4) :: dimsorder 
    CHARACTER(len = 128)  :: c_filename   ! file name

    nb_read = 0
    ind     = indf - 1

    DO WHILE (nb_read < ind_min)

      ind = ind + 1

      !- Build file name -!
      CALL sub_build_filename(ind, ndigits, prefix, suffix, c_filename)

      !- Open Netcdf file -!
      CALL sub_open_netcdf_file(cdir, c_filename, ncid)

      !- Read variable dimensions -!
      CALL sub_select_var_dims(ncid, ncvar, varid, dimsorder, &
           dimx, dimy, dimz, dimt)

      nb_read= nb_read + dimt

      CALL sub_close_netcdf_file(ncid)

    ENDDO

    ind_t = dimt - (nb_read - ind_min)

  END SUBROUTINE sub_search_forward_ind

  !!=========================================================================
  SUBROUTINE sub_search_backward_inds(  &
       ind_max, &
       ind_zo, ind_me, ind_ve, &
       ind_te, ind_sa, ind_rr, &
       ind_se, ind_tt)

    INTEGER(kind = iprec), INTENT(in)                :: &
         ind_max
    INTEGER(kind = iprec), INTENT(out)               :: &
         ind_zo, ind_me, ind_ve
    INTEGER(kind = iprec), OPTIONAL, INTENT(out)     :: &
         ind_te, ind_sa, ind_rr, ind_se
    INTEGER(kind = iprec), DIMENSION(:), INTENT(out) :: &
         ind_tt

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)'==========================================================='
    WRITE(lun_standard,*)'= Searching in Netcdf files indices corresponding to lmax ='
    WRITE(lun_standard,*)'==========================================================='

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)'- Zonal current file -'
    CALL sub_search_backward_ind(ind_max, ind_zo, ind_tt(1), &
         indn_zo, c_dir_zo, c_prefix_zo, maxsize_zo, c_suffix_zo, &
         nc_var_zo)
    WRITE(lun_standard,*)'    - file name indice: ', ind_zo
    WRITE(lun_standard,*)'    - file time indice: ', ind_tt(1)

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)'- Meridional current file -'
    CALL sub_search_backward_ind(ind_max, ind_me, ind_tt(2), &
         indn_me, c_dir_me, c_prefix_me, maxsize_me, c_suffix_me, &
         nc_var_me)
    WRITE(lun_standard,*)'    - file name indice: ', ind_me
    WRITE(lun_standard,*)'    - file time indice: ', ind_tt(2)

    IF ((key_computew).OR.(.NOT.key_read_w).OR.(key_roms).OR.(key_mars).OR.(key_symphonie)) THEN
       ind_ve    = 0
       ind_tt(3) = 0
    ELSE
       WRITE(lun_standard,*)''
       WRITE(lun_standard,*)'- Vertical current file -'
       CALL sub_search_backward_ind(ind_max, ind_ve, ind_tt(3)   , &
            indn_ve, c_dir_ve,c_prefix_ve, maxsize_ve, c_suffix_ve, &
            nc_var_ve)
       WRITE(lun_standard,*)'    - file name indice: ', ind_ve
       WRITE(lun_standard,*)'    - file time indice: ', ind_tt(3)
    ENDIF

    IF (PRESENT(ind_te)) THEN
       WRITE(lun_standard,*)''
       WRITE(lun_standard,*)'- Temperature file -'
       CALL sub_search_backward_ind(ind_max, ind_te, ind_tt(4)   , &
            indn_te, c_dir_te,c_prefix_te, maxsize_te, c_suffix_te, &
            nc_var_te)
       WRITE(lun_standard,*)'    - file name indice: ', ind_te
       WRITE(lun_standard,*)'    - file time indice: ', ind_tt(4)
    ENDIF
    IF (PRESENT(ind_sa)) THEN
       WRITE(lun_standard,*)''
       WRITE(lun_standard,*)'- Salinity file -'
       CALL sub_search_backward_ind(ind_max, ind_sa, ind_tt(5)   , &
            indn_sa, c_dir_sa,c_prefix_sa, maxsize_sa, c_suffix_sa, &
            nc_var_sa)
       WRITE(lun_standard,*)'    - file name indice: ', ind_sa
       WRITE(lun_standard,*)'    - file time indice: ', ind_tt(5)
    ENDIF

    IF (PRESENT(ind_rr)) THEN
       IF (key_computesigma) THEN
          ind_rr    = 0
          ind_tt(6) = 0
       ELSE
          WRITE(lun_standard,*)''
          WRITE(lun_standard,*)'- Density file -'

          CALL sub_search_backward_ind(ind_max, ind_rr, ind_tt(6)      , &
               indn_de, c_dir_de, c_prefix_de, maxsize_de, c_suffix_de, &
               nc_var_de)

          WRITE(lun_standard,*)'    - file name indice: ', ind_rr
          WRITE(lun_standard,*)'    - file time indice: ', ind_tt(6)
       ENDIF
    ENDIF

    IF (PRESENT(ind_se)) THEN
       IF (key_roms.or.key_mars) THEN
          WRITE(lun_standard,*)''
          WRITE(lun_standard,*)'- Zeta file (roms/mars)-'
          CALL sub_search_forward_ind(ind_max, ind_se, ind_tt(7)      , &
               ind0_ze, c_dir_ze, c_prefix_ze, maxsize_ze, c_suffix_ze, &
               nc_var_ze)
       ELSEIF(key_symphonie) THEN
          WRITE(lun_standard,*)''
          WRITE(lun_standard,*)'- SSE file (symphonie)-'
          CALL sub_search_forward_ind(ind_max, ind_se, ind_tt(7)           , &
               ind0_sse, c_dir_sse, c_prefix_sse, maxsize_sse, c_suffix_sse, &
               nc_var_sse)
       ENDIF
       WRITE(lun_standard,*)'    - file name indice: ', ind_se
       WRITE(lun_standard,*)'    - file time indice: ', ind_tt(7)
    ENDIF


    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)'================================'
    WRITE(lun_standard,*)'= Indice searches are finished ='
    WRITE(lun_standard,*)'================================'

  END SUBROUTINE sub_search_backward_inds

  !!=========================================================================
  SUBROUTINE sub_search_backward_ind( &
       ind_max,                        &
       ind, ind_t,                    &
       indf, cdir, prefix, ndigits, suffix, &
       ncvar )

    INTEGER(kind = iprec), INTENT(in) :: &
         ind_max
    INTEGER(kind = iprec), INTENT(out) :: &
         ind, & !
         ind_t  !
    INTEGER(kind = iprec), INTENT(in) :: &
         ndigits,indf
    CHARACTER(len = *)   , INTENT(in)  :: &
         cdir, prefix, suffix, ncvar

    INTEGER(kind = iprec) :: ncid
    INTEGER(kind = iprec) :: varid
    INTEGER(kind = iprec) :: dimx         ! dimension in x (i)
    INTEGER(kind = iprec) :: dimy         ! dimension in y (j)
    INTEGER(kind = iprec) :: dimz         ! dimension in z (k)
    INTEGER(kind = iprec) :: dimt         ! dimension in t (l)
    INTEGER(kind = iprec) :: nb_read
    INTEGER(kind = iprec), DIMENSION(4) :: dimsorder 
    CHARACTER(len = 128)  :: c_filename   ! file name

    nb_read = lmt
    ind     = indf

    DO WHILE (nb_read > ind_max)

      !- Build file name -!
      CALL sub_build_filename(ind, ndigits, prefix, suffix, c_filename)

      !- Open Netcdf file -!
      CALL sub_open_netcdf_file(cdir, c_filename, ncid)

      !- Read variable dimensions -!
      CALL sub_select_var_dims(ncid, ncvar, varid, dimsorder, &
           dimx, dimy, dimz, dimt)

      nb_read = nb_read - dimt

      CALL sub_close_netcdf_file(ncid)

      ind = ind - 1

    ENDDO

    ind_t = ind_max - nb_read + 1

  END SUBROUTINE sub_search_backward_ind

END MODULE mod_seq