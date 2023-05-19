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
!!****f* ariane/posini_seq()
!! NAME
!!   posini_seq() (posini-seq.f90)
!!
!! USAGE
!!   
!!
!! FUNCTION
!!   
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
!!   * i686-pc-linux-gnu -          g95
!!
!! SEE ALSO
!!
!! TODO
!!   * Should be transform to a module.
!!   * Should be migrate to F90 !!!
!!   * 3 similar loops, try to find a common base to reduce writing.
!!   * more modularity.
!!   * more comments.
!!
!! USES
!!   * USE mod_precision
!!   * USE mod_namelist
!!   * USE mod_input_grid
!!   * USE mod_input_data
!!   * USE mod_posin
!!   * USE mod_flags
!!   * USE mod_stats
!!   * USE mod_stati
!!   * USE mod_txt
!!   * USE mod_fx
!!   * USE mod_fy
!!   * USE mod_fz
!!   * USE mod_criter0
!!   * USE mod_sigma
!!   * USE mod_zinter
!!
!! USED BY
!!   * trajec
!!
!! SOURCE
!!=========================================================================
SUBROUTINE posini_seq(trmin, nb_traj, nb_sect, trans_total, pos_nfnt)
  !-----------------------------------------------------------------------
  ! automatic positioning: particles are gathered where the transport is
  ! the largest, giving an equivalent individual tranport to all
  ! particles (not to exceed max_transport)
  !------------------------------------------------------------------------

  !------------------!
  ! USE ASSOCIAITION !
  !------------------!
  USE mod_precision
  USE mod_namelist
  USE mod_input_grid
  USE mod_seq
  USE mod_input_data
  USE mod_init_particules
  USE mod_posin
  USE mod_flags
  USE mod_stats
  USE mod_stati
  USE mod_txt
  USE mod_fx
  USE mod_fy
  USE mod_fz
  USE mod_criter0
  USE mod_sigma
  USE mod_zinter
  USE mod_reducmem
  USE mod_lun

  !-------------!
  ! Declaration !
  !-------------!
  IMPLICIT NONE

  !- arguments -!
  INTEGER(kind=iprec), INTENT(out) :: nb_traj, nb_sect, pos_nfnt
  REAL(kind=rprec)   , INTENT(in)  :: trmin
  REAL(kind=rprec)   , INTENT(out) :: trans_total

  INTEGER(kind=iprec) :: &
       born1, & !
       born2, &
       alloc_size, & !
       n         , & !
       i         , & !
       k         , & !
       j         , & !
       l         , & !
       it0       , & !
       iu0       , & !
       it00      , & !
       n1        , & !
       n2        , & !
       nk        , & !
       nl        , & !
       j1        , & !
       k1        , & !
       l1        , & !
       jt0       , & !
       jv0       , & !
       jt00      , & !
       i1        , & !
       kt0       , & !
       kw0       , & !
       kt00      , & !
       is,js,ks  , & !
       it0s,jt0s,kt0s , & !
       iu0s,jv0s,kw0s 
  !
  INTEGER(kind=iprec), DIMENSION(:), ALLOCATABLE :: &
       num , & !
       it1, & !
       it2, & !
       jt1, & !
       jt2, & !
       kt1, & !
       kt2, & !
       ndir

  REAL(kind=rprec) :: u  ! zonal current
  REAL(kind=rprec) :: v  ! meridional current
  REAL(kind=rprec) :: w  ! vertical current

  REAL(kind=rprec) :: init_trans ! transport

  REAL(kind=rprec) :: prof ! depth

  LOGICAL :: switch=.FALSE.

  !- Comments -!
  WRITE(lun_standard,*)''
  WRITE(lun_standard,*)'  ---------------------------'
  WRITE(lun_standard,*)'  = ENTER >>>> POSINI (SEQ) ='
  WRITE(lun_standard,*)'  ---------------------------'

  !-------------!
  ! Code begins !
  !-------------!

  !- Dynamic Allocations -!

  CALL sub_txt_alloc()
  CALL sub_txt_init()

  CALL sub_flags_alloc(imt, jmt, kmt)
  CALL sub_flags_init(-1)

  CALL sub_stati_alloc()
  CALL sub_stati_init(0)

  CALL sub_stats_alloc() 
  CALL sub_stats_init()

  ALLOCATE(num(maxsegm)) ; num(:)  = 0
  CALL sub_memory(SIZE(num),'i','num','posini_seq')
  ALLOCATE(it1(maxsegm)) ; it1(:)  = 0
  CALL sub_memory(SIZE(it1),'i','it1','posini_seq')
  ALLOCATE(it2(maxsegm)) ; it2(:)  = 0
  CALL sub_memory(SIZE(it2),'i','it2','posini_seq')
  ALLOCATE(jt1(maxsegm)) ; jt1(:)  = 0
  CALL sub_memory(SIZE(jt1),'i','jt1','posini_seq')
  ALLOCATE(jt2(maxsegm)) ; jt2(:)  = 0
  CALL sub_memory(SIZE(jt2),'i','jt2','posini_seq')
  ALLOCATE(kt1(maxsegm)) ; kt1(:)  = 0
  CALL sub_memory(SIZE(kt1),'i','kt1','posini_seq')
  ALLOCATE(kt2(maxsegm)) ; kt2(:)  = 0
  CALL sub_memory(SIZE(kt2),'i','kt2','posini_seq')
  ALLOCATE(ndir(maxsegm)); ndir(:) = 0
  CALL sub_memory(SIZE(ndir),'i','ndir','posini_seq')

  !- Initializations -!
  nb_traj     = 0 
  pos_nfnt    = 0
  nb_sect     = 0
  trans_total = 0._rprec

  !- Open sections.txt file -!
  OPEN(unit=lun_dummy, file='sections.txt', action="read", &
       POSITION='rewind')

  WRITE(lun_standard,*)'--- sections.txt is opened: '

  !- DoLoop on segment number -!
  DO n = 1, maxsegm

    READ(lun_dummy,*) &
         num(n),it1(n),it2(n),jt1(n),jt2(n),kt1(n),kt2(n),sname

    WRITE(lun_standard,*)'    --- Reading section: ', sname

    ndir(n) = 0

    IF (key_ascii_outputs) THEN
      WRITE(lun_output,*)         &
           num(n),' '           , &
           it1(n),' ',it2(n),' ', &
           jt1(n),' ',jt2(n),' ', &
           kt1(n),' ',kt2(n),' ',sname
    ENDIF

    IF (num(n) < 0) THEN
      WRITE(lun_error,*)'    transparent sections are not allowed'
      STOP
    ENDIF

    IF (ABS(num(n)) > maxsect) THEN
      WRITE(lun_error,*)'    problem with segment #',n
      WRITE(lun_error,*)'      maximum number of sections is: ',maxsect
      STOP
    ENDIF

    IF ( (it1(n) /= it2(n)).AND. &
         (jt1(n) /= jt2(n)).AND. &
         (kt1(n) /= kt2(n))) THEN
      WRITE(lun_error,*)'    problem with segment #',n
      WRITE(lun_error,*)'      one pair of indices must be identical'
      STOP
    ENDIF

    IF ( ( (it1(n) == it2(n)).AND.(jt1(n) == jt2(n)) ).OR. &
         ( (it1(n) == it2(n)).AND.(kt1(n) == kt2(n)) ).OR. &
         ( (jt1(n) == jt2(n)).AND.(kt1(n) == kt2(n)) )  ) THEN
      WRITE(lun_error,*)'    problem with segment #',n
      WRITE(lun_error,*)'      only one pair of indices can be identical'
      STOP
    ENDIF

    secname(num(n)) = sname

    IF (num(n) > nb_sect) nb_sect = num(n)


    !-----------------------!
    !- Latitudinal section -!
    !-----------------------!
    IF (it1(n) == it2(n)) THEN

      ndir(n) = 1
      i       = it1(n)

      IF (i < 0) THEN
        i       = -i
        ndir(n) = -1
      ENDIF

      DO k = kt1(n), kt2(n)
        DO j = jt1(n), jt2(n)
          mtfin(i,j,k) = num(n)
        END DO
      END DO

    ENDIF

    IF (it1(n) < 0) it1(n) = -it1(n)

    IF (it2(n) < 0) it2(n) = -it2(n)

    !------------------------!
    !- Longitudinal section -!
    !------------------------!
    IF (jt1(n) == jt2(n)) THEN

      ndir(n) = 2
      j       = jt1(n)

      IF (j < 0) THEN
        j       = -j
        ndir(n) = -2
      ENDIF

      DO k = kt1(n), kt2(n)
        DO i = it1(n), it2(n)
          mtfin(i,j,k) = num(n)
        END DO
      END DO
    ENDIF

    IF (jt1(n) < 0) jt1(n) = -jt1(n)

    IF (jt2(n) < 0) jt2(n) = -jt2(n)


    !--------------------!
    !- Vertical section -!
    !--------------------!
    IF (kt1(n) == kt2(n)) THEN

      ndir(n) = 3
      k       = kt1(n)

      IF (k < 0) THEN
        k       = -k
        ndir(n) = -3
      ENDIF

      DO j = jt1(n), jt2(n)
        DO i = it1(n), it2(n)
          mtfin(i,j,k) = num(n)
        END DO
      END DO

    ENDIF

    IF (kt1(n) < 0) kt1(n) = -kt1(n)

    IF (kt2(n) < 0) kt2(n) = -kt2(n)

  END DO

  secname(0)='meanders'

  !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
  !!===============================================================================!!
  !!_______________________________________________________________________________!!

  IF (TRIM(bin) == 'nobin') THEN

    IF (key_ascii_outputs) THEN
      WRITE(lun_output,*)'    initial positions are automatically computed'
    ENDIF

    !!NG: to have the results in  same order than the old version of Ariane
    IF (TRIM(forback) == 'backward') THEN
      forback='forward'
      switch=.TRUE.
    ELSE
      switch=.FALSE.
    ENDIF

    IF (key_alltracers) THEN
      IF (key_roms.OR.key_mars) THEN
        alloc_size = 7
      ELSEIF(key_symphonie) THEN
        alloc_size = 8
      ELSE
        IF (TRIM(w_surf_option) == 'E-P-R') THEN
           alloc_size = 7
        ELSEIF (key_vvl) THEN
           alloc_size = 7
        ELSE
          alloc_size = 6
        ENDIF
      ENDIF
    ELSE
      IF (key_roms.OR.key_mars)THEN
        alloc_size = 3
      ELSEIF(key_symphonie) THEN
        alloc_size = 4
      ELSE
        IF (TRIM(w_surf_option) == 'E-P-R') THEN
           alloc_size = 4
        ELSEIF (key_vvl) THEN
           alloc_size = 4
        ELSE
          alloc_size = 3
        ENDIF
      ENDIF
    ENDIF

    CALL sub_seq_alloc(alloc_size)

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)'posini_seq: call sub_seq_init: start'

    IF ( (lmin /= 1).AND.(TRIM(forback) == 'forward')) THEN
      WRITE(lun_standard,*)'posini_seq: 1', lmin
      CALL sub_seq_init(lmin)
    ELSEIF (TRIM(forback) == 'backward') THEN  !! Never used !!
      WRITE(lun_standard,*)'posini_seq: 2', lmax
      CALL sub_seq_init(lmax)
    ELSE
      WRITE(lun_standard,*)'posini_seq: 3'
      CALL sub_seq_init()
    ENDIF
    WRITE(lun_standard,*)'posini_seq: call sub_seq_init: end', lmin, lmax

    DO l = lmin, lmax

      !! READ INPUT DATA SEQUENTIALLY !!
      IF (key_roms) THEN
        CALL sub_input_data_seq_roms( &
             ncids(:)               , &
             varids(:)              , &
             new_file(:)            , &
             ind_file(:)            , &
             ind_time(:)            , &
             ind_time_size(:)       , &
             sdimsorders(:,:)        )
      ELSEIF(key_symphonie) THEN
        CALL sub_input_data_seq_symphonie( &
             ncids(:)                    , &
             varids(:)                   , &
             new_file(:)                 , &
             ind_file(:)                 , &
             ind_time(:)                 , &
             ind_time_size(:)            , &
             sdimsorders(:,:)              )

      ELSEIF (key_mars) THEN
        CALL sub_input_data_seq_mars( &
             ncids(:)               , &
             varids(:)              , &
             new_file(:)            , &
             ind_file(:)            , &
             ind_time(:)            , &
             ind_time_size(:)       , &
             sdimsorders(:,:)       )

      ELSEIF (key_B2C_grid) THEN

        CALL sub_input_data_B2C_gridz_seq( &
             ncids(:)              , &
             varids(:)             , &
             new_file(:)           , &
             ind_file(:)           , &
             ind_time(:)           , &
             ind_time_size(:)      , &
             sdimsorders(:,:)       )

        !!NG: 16/04/2009 WRITE(lun_error,*)''
        !!NG: 16/04/2009 WRITE(lun_error,*)'STOP === key_B2C_grid in sequential mode is not available === STOP'
        !!NG: 16/04/2009 STOP

      ELSE
        CALL sub_input_data_seq_opa( &
             ncids(:)              , &
             varids(:)             , &
             new_file(:)           , &
             ind_file(:)           , &
             ind_time(:)           , &
             ind_time_size(:)      , &
             sdimsorders(:,:)       )
      ENDIF

      !!NG: Bug fixed 27/12/2006
      !!NG:    IF(TRIM(forback) == 'backward') THEN
      IF (switch) THEN
        uu(:,:,:,:) = -uu(:,:,:,:)
        vv(:,:,:,:) = -vv(:,:,:,:)
        ww(:,:,:,:) = -ww(:,:,:,:)
      ENDIF

      DO n = 1, maxsegm

        IF (num(n) == 1) THEN  !! only initial section

          IF ((ndir(n) == 1).OR.(ndir(n) == -1)) THEN
            it0  = it1(n)
            iu0  = it0
            it00 = it0 + 1

            IF (ndir(n) < 0) THEN
              iu0  = it0 - 1
              it00 = it0 - 1
            ENDIF

            DO k = kt1(n), kt2(n)
              DO j = jt1(n), jt2(n)

                CALL sub_reducmem_shift_or_not_ind(it0, j, k, it0s, js, ks)

                IF ( (tmask(it0s,js,ks,1) > 0.5).AND.(mtfin(it00,j,k) <= 0) ) THEN

                  !!NG: bug 14/04/2009 IF (key_nointerpolstats) THEN
                  !!NG: bug 14/04/2009   IF(key_alltracers) THEN
                  !!NG: bug 14/04/2009      init_temp(nb_traj) = tt(it0s,js,ks,1)
                  !!NG: bug 14/04/2009      init_salt(nb_traj) = ss(it0s,js,ks,1)
                  !!NG: bug 14/04/2009      init_dens(nb_traj) = rr(it0s,js,ks,1)
                  !!NG: bug 14/04/2009   ENDIF ! key_alltracers
                  !!NG: bug 14/04/2009 ENDIF ! key_nointerpolstats

                  CALL sub_reducmem_shift_or_not_ind(iu0, j, k, iu0s, js, ks)

                  u = uu(iu0s, js, ks, 1)

                  !----------------------------------
                  ! test on the sign of the velocity
                  !----------------------------------
                  IF ((ABS(u) > trmin).AND.(REAL(ndir(n),kind=rprec) * u > 0._rprec)) THEN

                    IF (key_2D) THEN
                      n1 = 1 + INT( SQRT(ABS(u)/max_transport) )
                      n2 = 1
                      nk = n2
                      nl = n1
                    ELSE
                      IF (lmt == 1) THEN
                        n1 = 1 + INT( SQRT(ABS(u)/max_transport) )
                        n2 = 1
                      ELSE
                        n1 = 1 + INT( (ABS(u)/max_transport)**(1._rprec/3._rprec) )
                        n2 = n1
                      ENDIF
                      nk = n1
                      nl = n2
                    ENDIF

                    init_trans=ABS(u)/REAL(n1,kind=rprec)/REAL(n1,kind=rprec)/REAL(n2,kind=rprec)

                    DO j1 = 1, n1
                      DO k1 = 1, nk
                        DO l1 = 1, nl

                          nb_traj = nb_traj + 1

                          IF (nb_traj > nmax) THEN
                            WRITE(lun_error,*)'too many particles, maximum is: ',nmax
                            STOP
                          ENDIF

                          tfi(nb_traj) = REAL(iu0,kind=rprec)

                          tfj(nb_traj) = REAL(j,kind=rprec) - 1._rprec - &
                               .5_rprec/REAL(n1,kind=rprec) + &
                               REAL(j1,kind=rprec)/REAL(n1,kind=rprec)

                          tfk(nb_traj) = REAL(k,kind=rprec)            - &
                               .5_rprec/REAL(nk,kind=rprec) + &
                               REAL(k1,kind=rprec)/REAL(nk,kind=rprec)

                          tfl(nb_traj) = REAL(l,kind=rprec) - .5_rprec - &
                               .5_rprec/REAL(nl,kind=rprec) + &
                               REAL(l1,kind=rprec)/REAL(nl,kind=rprec)

                          ttr(nb_traj) = init_trans

                          IF (criter0(tfi(nb_traj),tfj(nb_traj),-tfk(nb_traj), &
                               tfl(nb_traj),it0,j,k,l)) THEN

                            trans_total = trans_total + init_trans

                            isn(-1) = isn(-1) + 1
                            sn(-1)  =  sn(-1) + init_trans

                            IF (key_alltracers) THEN

                              IF (.NOT.key_nointerpolstats) THEN
                                init_temp(nb_traj) = zinter(tt, &
                                     tfi(nb_traj),tfj(nb_traj),-tfk(nb_traj), &
                                     1._rprec + (tfl(nb_traj) - INT(tfl(nb_traj))) )

                                init_salt(nb_traj) = zinter(ss, &
                                     tfi(nb_traj),tfj(nb_traj),-tfk(nb_traj), &
                                     1._rprec + (tfl(nb_traj) - INT(tfl(nb_traj))) )

                                IF (key_approximatesigma) THEN
                                  init_dens(nb_traj) = zinter(rr, & 
                                       tfi(nb_traj),tfj(nb_traj),-tfk(nb_traj), &
                                       1._rprec + (tfl(nb_traj) - INT(tfl(nb_traj))) )

                                ELSE ! key_approximatesigma
                                  init_dens(nb_traj) = sigma(zsigma, init_salt(nb_traj), init_temp(nb_traj))

                                ENDIF ! key_approximatesigma

                              ELSE                                     !!NG: bug correction 14/04/2009

                                init_temp(nb_traj) = tt(it0s,js,ks,1) !!NG: bug correction 14/04/2009
                                init_salt(nb_traj) = ss(it0s,js,ks,1) !!NG: bug correction 14/04/2009
                                init_dens(nb_traj) = rr(it0s,js,ks,1) !!NG: bug correction 14/04/2009

                              ENDIF ! key_nointerpolstats

                              IF (init_temp(nb_traj) <= s0min(5,-1)) s0min(5,-1) = init_temp(nb_traj)
                              IF (init_salt(nb_traj) <= s0min(6,-1)) s0min(6,-1) = init_salt(nb_traj)
                              IF (init_dens(nb_traj) <= s0min(7,-1)) s0min(7,-1) = init_dens(nb_traj)
                              IF (init_temp(nb_traj) >= s0max(5,-1)) s0max(5,-1) = init_temp(nb_traj)
                              IF (init_salt(nb_traj) >= s0max(6,-1)) s0max(6,-1) = init_salt(nb_traj)
                              IF (init_dens(nb_traj) >= s0max(7,-1)) s0max(7,-1) = init_dens(nb_traj)

                              s0x(5,-1)  = s0x(5,-1)  + init_temp(nb_traj)        * init_trans
                              s0x(6,-1)  = s0x(6,-1)  + init_salt(nb_traj)        * init_trans
                              s0x(7,-1)  = s0x(7,-1)  + init_dens(nb_traj)        * init_trans
                              s0x2(5,-1) = s0x2(5,-1) + init_temp(nb_traj) * init_temp(nb_traj) * init_trans
                              s0x2(6,-1) = s0x2(6,-1) + init_salt(nb_traj) * init_salt(nb_traj) * init_trans
                              s0x2(7,-1) = s0x2(7,-1) + init_dens(nb_traj) * init_dens(nb_traj) * init_trans

                            ENDIF !  key_alltracers

                            IF (key_roms.OR.key_mars.OR.key_symphonie) THEN
                              prof = &
                                   fz(gi=tfi(nb_traj), gj=tfj(nb_traj), gk=-tfk(nb_traj), il=1)
                            ELSE
                              prof = fz(gi=tfi(nb_traj), gj=tfj(nb_traj), gk=-tfk(nb_traj))
                            ENDIF

                            IF (fx(tfi(nb_traj),tfj(nb_traj)) <= s0min(1,-1)) &
                                 s0min(1,-1) = fx(tfi(nb_traj),tfj(nb_traj))

                            IF (fy(tfi(nb_traj),tfj(nb_traj)) <= s0min(2,-1)) &
                                 s0min(2,-1) = fy(tfi(nb_traj),tfj(nb_traj))

                            IF (-prof <= s0min(3,-1)) s0min(3,-1) = -prof

                            IF (fx(tfi(nb_traj),tfj(nb_traj)) >= s0max(1,-1)) &
                                 s0max(1,-1) = fx(tfi(nb_traj),tfj(nb_traj))

                            IF (fy(tfi(nb_traj),tfj(nb_traj)) >= s0max(2,-1)) &
                                 s0max(2,-1) = fy(tfi(nb_traj),tfj(nb_traj))

                            IF (-prof >= s0max(3,-1)) s0max(3,-1) = -prof

                            s0x(1,-1) = s0x(1,-1) + fx(tfi(nb_traj),tfj(nb_traj)) * init_trans

                            s0x(2,-1) = s0x(2,-1) + fy(tfi(nb_traj),tfj(nb_traj)) * init_trans

                            s0x(3,-1) = s0x(3,-1) - prof * init_trans

                            s0x2(1,-1) = s0x2(1,-1) + &
                                 fx(tfi(nb_traj),tfj(nb_traj)) * &
                                 fx(tfi(nb_traj),tfj(nb_traj)) * init_trans

                            s0x2(2,-1) = s0x2(2,-1) + &
                                 fy(tfi(nb_traj),tfj(nb_traj)) * &
                                 fy(tfi(nb_traj),tfj(nb_traj)) * init_trans

                            s0x2(3,-1) = s0x2(3,-1) + prof * prof * init_trans

                          ELSE

                            nb_traj = nb_traj - 1

                          ENDIF

                        END DO
                      END DO
                    END DO

                    ! fin de test sur signe de la vitesse
                  ENDIF

                  ! fin de test sur masque/critere
                ENDIF

                ! fins de boucle sur j,k,l
              END DO

            END DO

            ! fin test it1=it2
          ENDIF

          !!======================================================================
          IF ((ndir(n) == 2).OR.(ndir(n) == -2)) THEN

            jt0  = jt1(n)
            jv0  = jt0
            jt00 = jt0 + 1

            IF (ndir(n) < 0) THEN
              jv0  = jt0 - 1
              jt00 = jt0 - 1
            ENDIF

            DO k = kt1(n), kt2(n)

              IF (key_periodic) THEN
                born1 = MAX(2,it1(n))
                born2 = MIN(imt-1,it2(n))

              ELSE ! key_periodic
                born1 = it1(n)
                born2 = it2(n)

              ENDIF ! key_periodic

              DO i = born1, born2

                CALL sub_reducmem_shift_or_not_ind(i, jt0, k, is, jt0s, ks)

                IF ( (tmask(is,jt0s,ks,1) > 0.5).AND.(mtfin(i,jt00,k) <= 0) ) THEN

                  !!NG: bug 14/04/2009 IF (key_nointerpolstats) THEN
                  !!NG: bug 14/04/2009   IF (key_alltracers) THEN
                  !!NG: bug 14/04/2009      init_temp(nb_traj) = tt(is,jt0s,ks,1)
                  !!NG: bug 14/04/2009      init_salt(nb_traj) = ss(is,jt0s,ks,1)
                  !!NG: bug 14/04/2009      init_dens(nb_traj) = rr(is,jt0s,ks,1)
                  !!NG: bug 14/04/2009   ENDIF
                  !!NG: bug 14/04/2009 ENDIF

                  CALL sub_reducmem_shift_or_not_ind(i, jv0, k, is, jv0s, ks)

                  v = vv(is,jv0s,ks,1)

                  !------------------------------------!
                  !- test on the sign of the velocity -!
                  !------------------------------------!
                  IF ((ABS(v) > trmin).AND.(REAL(ndir(n),kind=rprec) * v > 0._rprec)) THEN
                    !
                    IF (key_2D) THEN
                      n1 = 1 + INT( SQRT(ABS(v)/max_transport) )
                      n2 = 1
                      nk = n2
                      nl = n1
                    ELSE
                      IF (lmt == 1) THEN
                        n1 = 1 + INT( SQRT(ABS(v)/max_transport) )
                        n2 = 1
                      ELSE
                        n1=  1 + INT( (ABS(v)/max_transport)**(1._rprec/3._rprec) )
                        n2= n1
                      ENDIF
                      nk = n1
                      nl = n2
                    ENDIF

                    init_trans = ABS(v)/REAL(n1,kind=rprec)/REAL(n1,kind=rprec)/REAL(n2,kind=rprec)

                    DO i1 = 1, n1
                      DO k1 = 1, nk
                        DO l1 = 1, nl

                          nb_traj = nb_traj + 1

                          IF (nb_traj > nmax) THEN
                            WRITE(lun_error,*)'too many particles, maximum is: ',nmax
                            STOP
                          ENDIF

                          tfi(nb_traj) = REAL(i,kind=rprec) -1._rprec - &
                               .5_rprec/REAL(n1,kind=rprec)           + &
                               REAL(i1,kind=rprec)/REAL(n1,kind=rprec)

                          tfj(nb_traj) = REAL(jv0,kind=rprec)

                          tfk(nb_traj) = REAL(k,kind=rprec)           - &
                               .5_rprec/REAL(nk,kind=rprec)           + &
                               REAL(k1,kind=rprec)/REAL(nk,kind=rprec)

                          tfl(nb_traj) = REAL(l,kind=rprec) -.5_rprec - &
                               .5_rprec/REAL(nl,kind=rprec)           + &
                               REAL(l1,kind=rprec)/REAL(nl,kind=rprec)

                          ttr(nb_traj) = init_trans

                          IF (criter0(tfi(nb_traj),tfj(nb_traj),-tfk(nb_traj), &
                               tfl(nb_traj),i,jt0,k,l)) THEN

                            trans_total = trans_total + init_trans

                            isn(-1) = isn(-1) + 1
                            sn(-1)  =  sn(-1) + init_trans

                            IF (key_alltracers) THEN

                              IF (.NOT.key_nointerpolstats) THEN
                                init_temp(nb_traj) = zinter(tt, &
                                     tfi(nb_traj),tfj(nb_traj),-tfk(nb_traj), &
                                     1._rprec + (tfl(nb_traj) - INT(tfl(nb_traj))) )

                                init_salt(nb_traj) = zinter(ss, &
                                     tfi(nb_traj),tfj(nb_traj),-tfk(nb_traj), &
                                     1._rprec + (tfl(nb_traj) - INT(tfl(nb_traj))) )

                                IF (key_approximatesigma) THEN
                                  init_dens(nb_traj) = zinter(rr, &
                                       tfi(nb_traj),tfj(nb_traj),-tfk(nb_traj), &
                                       1._rprec + (tfl(nb_traj) - INT(tfl(nb_traj))) )

                                ELSE
                                  init_dens(nb_traj) = sigma(zsigma, init_salt(nb_traj), init_temp(nb_traj))

                                ENDIF ! key_approximatesigma

                              ELSE                                     !!NG: bug correction 14/04/2009

                                init_temp(nb_traj) = tt(is,jt0s,ks,1) !!NG: bug correction 14/04/2009
                                init_salt(nb_traj) = ss(is,jt0s,ks,1) !!NG: bug correction 14/04/2009
                                init_dens(nb_traj) = rr(is,jt0s,ks,1) !!NG: bug correction 14/04/2009

                              ENDIF !  key_nointerpolstats

                              IF (init_temp(nb_traj) <= s0min(5,-1)) s0min(5,-1) = init_temp(nb_traj)
                              IF (init_salt(nb_traj) <= s0min(6,-1)) s0min(6,-1) = init_salt(nb_traj)
                              IF (init_dens(nb_traj) <= s0min(7,-1)) s0min(7,-1) = init_dens(nb_traj)
                              IF (init_temp(nb_traj) >= s0max(5,-1)) s0max(5,-1) = init_temp(nb_traj)
                              IF (init_salt(nb_traj) >= s0max(6,-1)) s0max(6,-1) = init_salt(nb_traj)
                              IF (init_dens(nb_traj) >= s0max(7,-1)) s0max(7,-1) = init_dens(nb_traj)

                              s0x(5,-1)  = s0x(5,-1)  + init_temp(nb_traj)        * init_trans
                              s0x(6,-1)  = s0x(6,-1)  + init_salt(nb_traj)        * init_trans
                              s0x(7,-1)  = s0x(7,-1)  + init_dens(nb_traj)        * init_trans
                              s0x2(5,-1) = s0x2(5,-1) + init_temp(nb_traj) * init_temp(nb_traj) * init_trans
                              s0x2(6,-1) = s0x2(6,-1) + init_salt(nb_traj) * init_salt(nb_traj) * init_trans
                              s0x2(7,-1) = s0x2(7,-1) + init_dens(nb_traj) * init_dens(nb_traj) * init_trans

                            ENDIF ! key_alltracers

                            IF (key_roms.OR.key_mars.OR.key_symphonie) THEN

                              prof = &
                                   fz(gi=tfi(nb_traj), gj=tfj(nb_traj), gk=-tfk(nb_traj), il=1)

                            ELSE
                              prof = fz(gi=tfi(nb_traj), gj=tfj(nb_traj), gk=-tfk(nb_traj))
                            ENDIF

                            IF (fx(tfi(nb_traj),tfj(nb_traj)) <= s0min(1,-1)) &
                                 s0min(1,-1) = fx(tfi(nb_traj),tfj(nb_traj))

                            IF (fy(tfi(nb_traj),tfj(nb_traj)) <= s0min(2,-1)) &
                                 s0min(2,-1) = fy(tfi(nb_traj),tfj(nb_traj))

                            IF (-prof <= s0min(3,-1)) s0min(3,-1) = -prof

                            IF (fx(tfi(nb_traj),tfj(nb_traj)) >= s0max(1,-1)) &
                                 s0max(1,-1)=fx(tfi(nb_traj),tfj(nb_traj))

                            IF (fy(tfi(nb_traj),tfj(nb_traj)) >= s0max(2,-1)) &
                                 s0max(2,-1) = fy(tfi(nb_traj),tfj(nb_traj))

                            IF (-prof >= s0max(3,-1)) s0max(3,-1) = -prof

                            s0x(1,-1) = s0x(1,-1) + fx(tfi(nb_traj),tfj(nb_traj)) * init_trans

                            s0x(2,-1) = s0x(2,-1) + fy(tfi(nb_traj),tfj(nb_traj)) * init_trans

                            s0x(3,-1) = s0x(3,-1) - prof * init_trans

                            s0x2(1,-1) = s0x2(1,-1) + &
                                 fx(tfi(nb_traj),tfj(nb_traj)) * &
                                 fx(tfi(nb_traj),tfj(nb_traj)) * init_trans

                            s0x2(2,-1) = s0x2(2,-1) + &
                                 fy(tfi(nb_traj),tfj(nb_traj)) * &
                                 fy(tfi(nb_traj),tfj(nb_traj)) * init_trans

                            s0x2(3,-1) = s0x2(3,-1) + prof * prof * init_trans

                          ELSE

                            nb_traj = nb_traj - 1

                          ENDIF

                        END DO
                      END DO
                    END DO
                    !
                    ! fin de test sur signe de la vitesse
                  ENDIF
                  !
                  ! fin de test sur masque/critere
                ENDIF
                !
                ! fins de boucle sur i,k,l
              END DO
            END DO
            !
            ! fin test jt1=jt2
          ENDIF

          !!======================================================================
          IF ((ndir(n) == 3).OR.(ndir(n) == -3)) THEN

            kt0  = kt1(n)
            kw0  = kt0 + 1
            kt00 = kt0 + 1

            IF (ndir(n) < 0) THEN
              kw0  = kt0
              kt00 = kt0 - 1
            ENDIF

            IF (kt0 == 0) THEN
              kt0 = 1
            ENDIF

            DO j = jt1(n), jt2(n)
              DO i = it1(n), it2(n)

                CALL sub_reducmem_shift_or_not_ind(i, j, kt0, is, js, kt0s)

                IF ( (tmask(is,js,kt0s,1) > 0.5).AND.(mtfin(i,j,kt00) <= 0) ) THEN

                  !!NG: bug 14/04/2009 IF (key_nointerpolstats) THEN
                  !!NG: bug 14/04/2009    IF (key_alltracers) THEN
                  !!NG: bug 14/04/2009       init_temp(nb_traj) = tt(is,js,kt0s,1)
                  !!NG: bug 14/04/2009       init_salt(nb_traj) = ss(is,js,kt0s,1)
                  !!NG: bug 14/04/2009       init_dens(nb_traj) = rr(is,js,kt0s,1)
                  !!NG: bug 14/04/2009    ENDIF
                  !!NG: bug 14/04/2009 ENDIF

                  CALL sub_reducmem_shift_or_not_ind(i, j, kw0, is, js, kw0s)

                  w=ww(is,js,kw0s,1)

                  !------------------------------------!
                  !- test on the sign of the velocity -!
                  !------------------------------------!
                  IF ((ABS(w) > trmin).AND.(REAL(ndir(n),kind=rprec) * w < 0._rprec)) THEN

                    IF (lmt == 1) THEN
                      n1 = 1 + INT( SQRT(ABS(w)/max_transport) )
                      n2 = 1
                    ELSE
                      n1 = 1 + INT( (ABS(w)/max_transport)**(1._rprec/3._rprec) )
                      n2 = n1
                    ENDIF

                    init_trans = ABS(w)/REAL(n1,kind=rprec)/REAL(n1,kind=rprec)/REAL(n2,kind=rprec)

                    DO j1 = 1, n1
                      DO i1 = 1, n1
                        DO l1 = 1, n2

                          nb_traj = nb_traj + 1

                          IF (nb_traj > nmax) THEN
                            WRITE(lun_error,*)'too many particles, maximum is: ', nmax
                            STOP
                          ENDIF

                          tfi(nb_traj) = REAL(i,kind=rprec)   -1._rprec - &
                               .5_rprec/REAL(n1,kind=rprec)             + &
                               REAL(i1,kind=rprec)/REAL(n1,kind=rprec)

                          tfj(nb_traj) = REAL(j,kind=rprec)   -1._rprec - &
                               .5_rprec/REAL(n1,kind=rprec)             + &
                               REAL(j1,kind=rprec)/REAL(n1,kind=rprec)

                          tfk(nb_traj) = REAL(kw0,kind=rprec)

                          tfl(nb_traj) = REAL(l,kind=rprec)   -.5_rprec - &
                               .5_rprec/REAL(n2,kind=rprec)             + &
                               REAL(l1,kind=rprec)/REAL(n2,kind=rprec)

                          ttr(nb_traj) = init_trans
                          !
                          IF (criter0(tfi(nb_traj),tfj(nb_traj),-tfk(nb_traj), &
                               tfl(nb_traj),i,j,kt0,l)) THEN
                            !
                            trans_total = trans_total + init_trans

                            isn(-1) = isn(-1) + 1
                            sn(-1)  =  sn(-1) + init_trans

                            IF (key_alltracers) THEN

                              IF (.NOT.key_nointerpolstats) THEN
                                init_temp(nb_traj) = zinter(tt, &
                                     tfi(nb_traj),tfj(nb_traj),-tfk(nb_traj), &
                                     1._rprec + (tfl(nb_traj) - INT(tfl(nb_traj))) )
                                init_salt(nb_traj) = zinter(ss, &
                                     tfi(nb_traj),tfj(nb_traj),-tfk(nb_traj), &
                                     1._rprec + (tfl(nb_traj) - INT(tfl(nb_traj))) )

                                IF (key_approximatesigma) THEN
                                  init_dens(nb_traj) = zinter(rr, &
                                       tfi(nb_traj),tfj(nb_traj),-tfk(nb_traj), &
                                       1._rprec + (tfl(nb_traj) - INT(tfl(nb_traj))) )

                                ELSE
                                  init_dens(nb_traj) = sigma(zsigma, init_salt(nb_traj), init_temp(nb_traj))

                                ENDIF

                              ELSE                                     !!NG: bug correction 14/04/2009

                                init_temp(nb_traj) = tt(is,js,kt0s,1) !!NG: bug correction 14/04/2009
                                init_salt(nb_traj) = ss(is,js,kt0s,1) !!NG: bug correction 14/04/2009
                                init_dens(nb_traj) = rr(is,js,kt0s,1) !!NG: bug correction 14/04/2009

                              ENDIF

                              IF (init_temp(nb_traj) <= s0min(5,-1)) s0min(5,-1) = init_temp(nb_traj)
                              IF (init_salt(nb_traj) <= s0min(6,-1)) s0min(6,-1) = init_salt(nb_traj)
                              IF (init_dens(nb_traj) <= s0min(7,-1)) s0min(7,-1) = init_dens(nb_traj)
                              IF (init_temp(nb_traj) >= s0max(5,-1)) s0max(5,-1) = init_temp(nb_traj)
                              IF (init_salt(nb_traj) >= s0max(6,-1)) s0max(6,-1) = init_salt(nb_traj)
                              IF (init_dens(nb_traj) >= s0max(7,-1)) s0max(7,-1) = init_dens(nb_traj)

                              s0x(5,-1)  = s0x(5,-1)  + init_temp(nb_traj)        * init_trans
                              s0x(6,-1)  = s0x(6,-1)  + init_salt(nb_traj)        * init_trans
                              s0x(7,-1)  = s0x(7,-1)  + init_dens(nb_traj)        * init_trans
                              s0x2(5,-1) = s0x2(5,-1) + init_temp(nb_traj) * init_temp(nb_traj) * init_trans
                              s0x2(6,-1) = s0x2(6,-1) + init_salt(nb_traj) * init_salt(nb_traj) * init_trans
                              s0x2(7,-1) = s0x2(7,-1) + init_dens(nb_traj) * init_dens(nb_traj) * init_trans


                            ENDIF

                            IF (key_roms.OR.key_mars.OR.key_symphonie) THEN
                              prof = &
                                   fz(gi=tfi(nb_traj), gj=tfj(nb_traj), gk=-tfk(nb_traj), il=1)
                            ELSE
                              prof = fz(gi=tfi(nb_traj), gj=tfj(nb_traj), gk=-tfk(nb_traj))
                            ENDIF

                            IF (fx(tfi(nb_traj),tfj(nb_traj)) <= s0min(1,-1)) &
                                 s0min(1,-1) = fx(tfi(nb_traj),tfj(nb_traj))

                            IF (fy(tfi(nb_traj),tfj(nb_traj)) <= s0min(2,-1)) &
                                 s0min(2,-1) = fy(tfi(nb_traj),tfj(nb_traj))

                            IF (-prof <= s0min(3,-1)) s0min(3,-1) = -prof

                            IF (fx(tfi(nb_traj),tfj(nb_traj)) >= s0max(1,-1)) &
                                 s0max(1,-1) = fx(tfi(nb_traj),tfj(nb_traj))

                            IF (fy(tfi(nb_traj),tfj(nb_traj)) >= s0max(2,-1)) &
                                 s0max(2,-1) = fy(tfi(nb_traj),tfj(nb_traj))

                            IF (-prof >= s0max(3,-1)) s0max(3,-1) = -prof

                            s0x(1,-1)  = s0x(1,-1) + fx(tfi(nb_traj),tfj(nb_traj))*init_trans

                            s0x(2,-1)  = s0x(2,-1) + fy(tfi(nb_traj),tfj(nb_traj))*init_trans

                            s0x(3,-1)  = s0x(3,-1) - prof * init_trans

                            s0x2(1,-1) = s0x2(1,-1) + &
                                 fx(tfi(nb_traj),tfj(nb_traj)) * &
                                 fx(tfi(nb_traj),tfj(nb_traj)) * init_trans

                            s0x2(2,-1) = s0x2(2,-1) + &
                                 fy(tfi(nb_traj),tfj(nb_traj)) * &
                                 fy(tfi(nb_traj),tfj(nb_traj)) * init_trans

                            s0x2(3,-1) = s0x2(3,-1) + prof * prof * init_trans

                          ELSE

                            nb_traj = nb_traj - 1

                          ENDIF

                        END DO
                      END DO
                    END DO
                    !
                    ! fin de test sur signe de la vitesse
                  ENDIF
                  !
                  ! fin de test sur masque/critere
                ENDIF
                !
                ! fins de boucle sur i,j,l
              END DO
            END DO
            !
            ! fin test kt1=kt2
          ENDIF
          !
          ! fin test num(n)
        ENDIF
        !
        ! fin du positionnement
      END DO
    END DO

    IF (switch) THEN
      forback='backward'
    ENDIF

  ELSE

    WRITE(lun_standard,*)'     - in Posini just to initialize some arrays...'
    WRITE(lun_standard,*)'       bin = ', TRIM(bin)

  ENDIF

  IF (ALLOCATED(num)) THEN
    CALL sub_memory(-SIZE(num),'i','num','posini_seq')
    DEALLOCATE(num)
  END IF
  IF (ALLOCATED(it1)) THEN
    CALL sub_memory(-SIZE(it1),'i','it1','posini_seq')
    DEALLOCATE(it1)
  END IF
  IF (ALLOCATED(it2)) THEN
    CALL sub_memory(-SIZE(it2),'i','it2','posini_seq')
    DEALLOCATE(it2)
  END IF
  IF (ALLOCATED(jt1)) THEN
    CALL sub_memory(-SIZE(jt1),'i','jt1','posini_seq')
    DEALLOCATE(jt1)
  END IF
  IF (ALLOCATED(jt2)) THEN
    CALL sub_memory(-SIZE(jt2),'i','jt2','posini_seq')
    DEALLOCATE(jt2)
  END IF
  IF (ALLOCATED(kt1)) THEN
    CALL sub_memory(-SIZE(kt1),'i','kt1','posini_seq')
    DEALLOCATE(kt1)
  END IF
  IF (ALLOCATED(kt2)) THEN
    CALL sub_memory(-SIZE(kt2),'i','kt2','posini_seq')
    DEALLOCATE(kt2)
  END IF
  IF ( ALLOCATED(ndir)) THEN
    CALL sub_memory(-SIZE(ndir),'i','ndir','posini_seq')
    DEALLOCATE(ndir)
  END IF

  WRITE(lun_standard,*)''
  WRITE(lun_standard,*)'  ----------------------------'
  WRITE(lun_standard,*)'  =  POSINI (SEQ)  >>>> EXIT ='
  WRITE(lun_standard,*)'  ----------------------------'

  RETURN

END SUBROUTINE posini_seq
