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
!!****f* ariane/posini()
!! NAME
!!   posini() (posini.f90)
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
!!   * more modularity
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
SUBROUTINE posini(trmin, nb_traj,nb_sect,trans_total,pos_nfnt)
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
  USE mod_input_data
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
  USE mod_save_netcdf

  !-------------!
  ! Declaration !
  !-------------!
  IMPLICIT NONE

  !- arguments -!
  INTEGER(kind=iprec), INTENT(out) :: nb_traj, nb_sect, pos_nfnt
  REAL(kind=rprec)   , INTENT(in)  :: trmin
  REAL(kind=rprec)   , INTENT(out) :: trans_total

  !- local variables -!
  INTEGER(kind=iprec), DIMENSION(:), ALLOCATABLE :: &
       it1, & !
       it2, & !
       jt1, & ! 
       jt2, & !
       kt1, & !
       kt2, & !
       num, & !
       ndir   !

  INTEGER(kind=iprec) :: n     !
  INTEGER(kind=iprec) :: i,j,k
  INTEGER(kind=iprec) :: is,js,ks
  INTEGER(kind=iprec) :: it0s,jt0s,kt0s
  INTEGER(kind=iprec) :: iu0s,jv0s,kw0s
  INTEGER(kind=iprec) :: l     !
  INTEGER(kind=iprec) :: iu0   !
  INTEGER(kind=iprec) :: it0   !
  INTEGER(kind=iprec) :: it00  !
  INTEGER(kind=iprec) :: n1    !
  INTEGER(kind=iprec) :: n2    !
  INTEGER(kind=iprec) :: i1    !
  INTEGER(kind=iprec) :: j1    !
  INTEGER(kind=iprec) :: k1    !
  INTEGER(kind=iprec) :: l1    !
  INTEGER(kind=iprec) :: jt0   !
  INTEGER(kind=iprec) :: jt00  !
  INTEGER(kind=iprec) :: ilim1 !
  INTEGER(kind=iprec) :: ilim2 !
  INTEGER(kind=iprec) :: jv0   !
  INTEGER(kind=iprec) :: kt0   !
  INTEGER(kind=iprec) :: kt00  !
  INTEGER(kind=iprec) :: kw0   !

  REAL(kind=rprec) :: u    ! zonal current
  REAL(kind=rprec) :: v    ! meridional current
  REAL(kind=rprec) :: w    ! vertical current
  REAL(kind=rprec) :: init_trans
  REAL(kind=rprec) :: prof ! W depth


  !- Comments -!
  WRITE(lun_standard,*)''
  WRITE(lun_standard,*)'  ---------------------'
  WRITE(lun_standard,*)'  = ENTER >>>> POSINI ='
  WRITE(lun_standard,*)'  ---------------------'

  !-------------!
  ! Code begins !
  !-------------!
  !- Dynamic Allocations -!

  CALL sub_txt_alloc()
  CALL sub_txt_init()

  CALL sub_flags_alloc(imt, jmt, kmt)
  CALL sub_flags_init(0)

  CALL sub_stati_alloc()
  CALL sub_stati_init(0)

  CALL sub_stats_alloc() 
  CALL sub_stats_init()

  ALLOCATE(it1(maxsegm))
  CALL sub_memory(size(it1),'i','it1','posini')
  ALLOCATE(it2(maxsegm))
  CALL sub_memory(size(it2),'i','it2','posini')
  ALLOCATE(jt1(maxsegm))
  CALL sub_memory(size(jt1),'i','jt1','posini')
  ALLOCATE(jt2(maxsegm))
  CALL sub_memory(size(jt2),'i','jt2','posini')
  ALLOCATE(kt1(maxsegm))
  CALL sub_memory(size(kt1),'i','kt1','posini')
  ALLOCATE(kt2(maxsegm))
  CALL sub_memory(size(kt2),'i','kt2','posini')
  ALLOCATE(num(maxsegm))
  CALL sub_memory(size(num),'i','num','posini')
  ALLOCATE(ndir(maxsegm))
  CALL sub_memory(size(ndir),'i','ndir','posini')

  !- Initializations -!
  nb_traj = 0
  pos_nfnt  = 0
  nb_sect = 0
  trans_total = 0._rprec

  !- Open sections.txt file -!
  OPEN(unit=lun_dummy, file='sections.txt', action="read", &
       POSITION='rewind')

  !- DoLoop on segment number -!
  DO n = 1, maxsegm

     READ(lun_dummy,*) num(n),it1(n),it2(n),jt1(n),jt2(n),kt1(n),kt2(n),sname

     ndir(n) = 0

     WRITE(lun_output,*)num(n),' ',it1(n),' ',it2(n),' ',jt1(n),' ',jt2(n), &
          ' ',kt1(n),' ',kt2(n),' ',sname

     IF (num(n) < 0) THEN
        WRITE(lun_standard,*)''
        WRITE(lun_standard,*)'transparent section associated with segment #',n
        pos_nfnt = max0(pos_nfnt,-num(n))
     ENDIF

     IF (ABS(num(n)) > maxsect) THEN
        WRITE(lun_error,*)'Problem with segment #'         , n
        WRITE(lun_error,*)'Maximum number of sections is: ', maxsect
        STOP
     ENDIF

     IF ((it1(n).NE.it2(n)).AND. &
          (jt1(n).NE.jt2(n)).AND. &
          (kt1(n).NE.kt2(n))) THEN
        WRITE(lun_error,*)'Problem with segment #',n
        WRITE(lun_error,*)'One pair of indices must be identical.'
        STOP
     ENDIF

     IF ( ( (it1(n) == it2(n)).AND.(jt1(n) == jt2(n)) ).OR. &
          ( (it1(n) == it2(n)).AND.(kt1(n) == kt2(n)) ).OR. &
          ( (jt1(n) == jt2(n)).AND.(kt1(n) == kt2(n)) )  ) THEN
        WRITE(lun_error,*)'Problem with segment #',n
        WRITE(lun_error,*)'Only one pair of indices can be identical'
        STOP
     ENDIF

     IF (num(n) > 0) THEN
        secname(num(n)) = sname
        IF (num(n) > nb_sect) nb_sect = num(n)
     ELSE
        fntname(-num(n)) = sname
     ENDIF

     !-----------------------!
     !- Latitudinal section -!
     !-----------------------!
     IF (it1(n) == it2(n)) THEN

        ndir(n) = 1
        i = it1(n)

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

     WRITE(lun_standard,*)'-- initial positions are automatically computed --'
     WRITE(lun_standard,*)'                 bin = ', TRIM(bin)

     !----------------!
     !- Loop on time -!
     !----------------!
     DO l = lmin, lmax

        !--------------------!
        !- Loop on segments -!
        !--------------------!
        DO n = 1, maxsegm

           !------------------------------------- --!
           !-    Test if the segment number = 1    -!
           !- Segment where particules will start. -!
           !------------------------------------- --!
           IF (num(n) == 1) THEN

              !!======================================================================
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

                       IF ( (tmask(it0s,js,ks,1) >  0.5).AND. &
                            (mtfin(it00,j,k)     <= 0  ) )    THEN

                          !!NG: bug 14/04/2009 IF (key_nointerpolstats) THEN
                          !!NG: bug 14/04/2009   IF(key_alltracers) THEN
                          !!NG: bug 14/04/2009      init_temp(nb_traj) = tt(it0s,js,ks,1)
                          !!NG: bug 14/04/2009      init_salt(nb_traj) = ss(it0s,js,ks,1)
                          !!NG: bug 14/04/2009      init_dens(nb_traj) = rr(it0s,js,ks,1)
                          !!NG: bug 14/04/2009   ENDIF ! key_alltracers
                          !!NG: bug 14/04/2009 ENDIF ! key_nointerpolstats

                          CALL sub_reducmem_shift_or_not_ind(iu0, j, k, iu0s, js, ks)

                          u = uu(iu0s, js, ks, l)

                          !----------------------------------
                          ! test on the sign of the velocity
                          !----------------------------------
                          IF ((ABS(u) > trmin).AND.(REAL(ndir(n), kind = rprec)*u > 0._rprec)) THEN

                             IF (lmt == 1) THEN
                                n1 = 1 + INT( SQRT(ABS(u)/max_transport) )
                                n2 = 1
                             ELSE
                                n1 = 1 + INT( (ABS(u)/max_transport)**(1._rprec/3._rprec) )
                                n2 = n1
                             ENDIF

                             init_trans = ABS(u)/REAL(n1,kind=rprec)/REAL(n1,kind=rprec)/REAL(n2,kind=rprec)

                             DO j1 = 1, n1
                                DO k1 = 1, n1
                                   DO l1 = 1, n2

                                      nb_traj = nb_traj + 1

                                      IF (nb_traj > nmax) THEN
                                         WRITE(lun_error,*)'too many particles, maximum is: ', nmax
                                         STOP
                                      ENDIF

                                      tfi(nb_traj) = REAL(iu0, kind=rprec)

                                      tfj(nb_traj) = REAL(j, kind=rprec) - 1._rprec - &
                                           .5_rprec/REAL(n1, kind=rprec)            + &
                                           REAL(j1, kind=rprec)/REAL(n1, kind=rprec)

                                      tfk(nb_traj) = REAL(k, kind=rprec)            - &
                                           .5_rprec/REAL(n1, kind=rprec)            + &
                                           REAL(k1, kind=rprec)/REAL(n1, kind=rprec)

                                      tfl(nb_traj) = REAL(l, kind=rprec) - .5_rprec - &
                                           .5_rprec/REAL(n2, kind=rprec)            + &
                                           REAL(l1, kind=rprec)/REAL(n2, kind=rprec)

                                      ttr(nb_traj) = init_trans

                                      IF (criter0(tfi(nb_traj),tfj(nb_traj),-tfk(nb_traj), &
                                           tfl(nb_traj),it0,j,k,l)) THEN

                                         trans_total = trans_total + init_trans

                                         isn(-1) = isn(-1) + 1
                                         sn(-1)  =  sn(-1) + init_trans

                                         IF (key_alltracers) THEN

                                            IF (.NOT.key_nointerpolstats) THEN
                                               init_temp(nb_traj) = zinter(tt, &
                                                    tfi(nb_traj),tfj(nb_traj),-tfk(nb_traj),tfl(nb_traj))
                                               init_salt(nb_traj) = zinter(ss, &
                                                    tfi(nb_traj),tfj(nb_traj),-tfk(nb_traj),tfl(nb_traj))
                                               IF (key_approximatesigma) THEN
                                                  init_dens(nb_traj) = zinter(rr, &
                                                       tfi(nb_traj),tfj(nb_traj),-tfk(nb_traj),tfl(nb_traj))
                                               ELSE
                                                  init_dens(nb_traj) = sigma(zsigma, init_salt(nb_traj), init_temp(nb_traj))
                                               ENDIF

                                            ELSE                                    !!NG: bug correction 14/04/2009

                                               init_temp(nb_traj) = tt(it0s,js,ks,1) !!NG: bug correction 14/04/2009
                                               init_salt(nb_traj) = ss(it0s,js,ks,1) !!NG: bug correction 14/04/2009
                                               init_dens(nb_traj) = rr(it0s,js,ks,1) !!NG: bug correction 14/04/2009
                                            ENDIF

                                            IF (init_temp(nb_traj) <= s0min(5,-1)) s0min(5,-1) = init_temp(nb_traj)
                                            IF (init_salt(nb_traj) <= s0min(6,-1)) s0min(6,-1) = init_salt(nb_traj)
                                            IF (init_dens(nb_traj) <= s0min(7,-1)) s0min(7,-1) = init_dens(nb_traj)
                                            IF (init_temp(nb_traj) >= s0max(5,-1)) s0max(5,-1) = init_temp(nb_traj)
                                            IF (init_salt(nb_traj) >= s0max(6,-1)) s0max(6,-1) = init_salt(nb_traj)
                                            IF (init_dens(nb_traj) >= s0max(7,-1)) s0max(7,-1) = init_dens(nb_traj)

                                            s0x(5,-1)  = s0x(5,-1)  + init_temp(nb_traj)     * init_trans
                                            s0x(6,-1)  = s0x(6,-1)  + init_salt(nb_traj)     * init_trans
                                            s0x(7,-1)  = s0x(7,-1)  + init_dens(nb_traj)     * init_trans
                                            s0x2(5,-1) = s0x2(5,-1) + init_temp(nb_traj) * init_temp(nb_traj) * init_trans
                                            s0x2(6,-1) = s0x2(6,-1) + init_salt(nb_traj) * init_salt(nb_traj) * init_trans
                                            s0x2(7,-1) = s0x2(7,-1) + init_dens(nb_traj) * init_dens(nb_traj) * init_trans

                                         ENDIF

                                         IF (key_roms.OR.key_symphonie) THEN

                                            prof = &
                                                 fz(gi=tfi(nb_traj), gj=tfj(nb_traj), gk=-tfk(nb_traj), il=l)

                                         ELSE
                                            prof = fz(gi=tfi(nb_traj), gj=tfj(nb_traj), gk=-tfk(nb_traj))
                                         ENDIF

                                         IF (fx(tfi(nb_traj),tfj(nb_traj)) <= s0min(1,-1)) &
                                              s0min(1,-1) = fx(tfi(nb_traj),tfj(nb_traj))

                                         IF (fy(tfi(nb_traj),tfj(nb_traj)) <= s0min(2,-1)) &
                                              s0min(2,-1) = fy(tfi(nb_traj),tfj(nb_traj))

                                         IF (-prof <= s0min(3,-1)) &
                                              s0min(3,-1) = -prof

                                         IF (fx(tfi(nb_traj),tfj(nb_traj)) >= s0max(1,-1)) &
                                              s0max(1,-1) = fx(tfi(nb_traj),tfj(nb_traj))

                                         IF (fy(tfi(nb_traj),tfj(nb_traj)) >= s0max(2,-1)) &
                                              s0max(2,-1) = fy(tfi(nb_traj),tfj(nb_traj))

                                         IF (-prof >= s0max(3,-1)) s0max(3,-1) = -prof

                                         s0x(1,-1)  = s0x(1,-1)  + fx(tfi(nb_traj),tfj(nb_traj)) * init_trans

                                         s0x(2,-1)  = s0x(2,-1)  + fy(tfi(nb_traj),tfj(nb_traj)) * init_trans

                                         s0x(3,-1)  = s0x(3,-1)  - prof                          * init_trans

                                         s0x2(1,-1) = s0x2(1,-1) + fx(tfi(nb_traj),tfj(nb_traj)) * &
                                              fx(tfi(nb_traj),tfj(nb_traj)) * init_trans

                                         s0x2(2,-1)  = s0x2(2,-1)+ fy(tfi(nb_traj),tfj(nb_traj)) * &
                                              fy(tfi(nb_traj),tfj(nb_traj)) * init_trans

                                         s0x2(3,-1)  = s0x2(3,-1)+ prof * prof * init_trans

                                      ELSE

                                         nb_traj=nb_traj-1

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
                       ! fins de boucle sur j,k,l
                    END DO
                 END DO
                 !!NG:          END DO
                 !
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

                 DO k = kt1(n),kt2(n)

                    IF (key_periodic) THEN
                       ilim1 = MAX(2,it1(n))
                       ilim2 = MIN(imt-1,it2(n))
                    ELSE
                       ilim1 = it1(n)
                       ilim2 = it2(n)
                    ENDIF

                    DO i = ilim1, ilim2

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

                          v = vv(is,jv0s,ks,l)

                          !----------------------------------
                          ! test on the sign of the velocity
                          !----------------------------------
                          IF ((ABS(v) > trmin).AND.(REAL(ndir(n),kind=rprec)*v > 0._rprec)) THEN

                             IF (lmt == 1) THEN
                                n1 = 1 + INT( SQRT(ABS(v)/max_transport) )
                                n2 = 1
                             ELSE
                                n1 = 1 + INT( (ABS(v)/max_transport)**(1._rprec/3._rprec) )
                                n2 = n1
                             ENDIF

                             init_trans = ABS(v)/REAL(n1,kind=rprec)/REAL(n1,kind=rprec)/REAL(n2,kind=rprec)

                             DO i1 = 1,n1
                                DO k1 = 1,n1
                                   DO l1 = 1,n2

                                      nb_traj = nb_traj + 1

                                      IF (nb_traj > nmax) THEN
                                         WRITE(lun_error,*)'too many particles, maximum is: ',nmax
                                         STOP
                                      ENDIF

                                      tfi(nb_traj) = REAL(i,kind=rprec)   -1._rprec - &
                                           .5_rprec/REAL(n1,kind=rprec)             + &
                                           REAL(i1,kind=rprec)/REAL(n1,kind=rprec)

                                      tfj(nb_traj) = REAL(jv0,kind=rprec)

                                      tfk(nb_traj) = REAL(k,kind=rprec)             - &
                                           .5_rprec/REAL(n1,kind=rprec)             + &
                                           REAL(k1,kind=rprec)/REAL(n1,kind=rprec)

                                      tfl(nb_traj) = REAL(l,kind=rprec)   -.5_rprec - &
                                           .5_rprec/REAL(n2,kind=rprec)             + &
                                           REAL(l1,kind=rprec)/REAL(n2,kind=rprec)

                                      ttr(nb_traj) = init_trans

                                      IF (criter0(tfi(nb_traj),tfj(nb_traj),-tfk(nb_traj), &
                                           tfl(nb_traj),i,jt0,k,l)) THEN

                                         trans_total = trans_total + init_trans

                                         isn(-1) = isn(-1) + 1
                                         sn(-1)  = sn(-1)  + init_trans

                                         IF (key_alltracers) THEN

                                            IF (.NOT.key_nointerpolstats) THEN
                                               init_temp(nb_traj) = zinter(tt, &
                                                    tfi(nb_traj),tfj(nb_traj),-tfk(nb_traj),tfl(nb_traj))
                                               init_salt(nb_traj) = zinter(ss, &
                                                    tfi(nb_traj),tfj(nb_traj),-tfk(nb_traj),tfl(nb_traj))
                                               IF (key_approximatesigma) THEN
                                                  init_dens(nb_traj) = zinter(rr, &
                                                       tfi(nb_traj),tfj(nb_traj),-tfk(nb_traj),tfl(nb_traj))
                                               ELSE
                                                  init_dens(nb_traj)  = sigma(zsigma, init_salt(nb_traj), init_temp(nb_traj))
                                               ENDIF

                                            ELSE                                    !!NG: bug correction 14/04/2009

                                               init_temp(nb_traj) = tt(is,jt0s,ks,1) !!NG: bug correction 14/04/2009
                                               init_salt(nb_traj) = ss(is,jt0s,ks,1) !!NG: bug correction 14/04/2009
                                               init_dens(nb_traj) = rr(is,jt0s,ks,1) !!NG: bug correction 14/04/2009
                                            ENDIF

                                            IF (init_temp(nb_traj) <= s0min(5,-1)) s0min(5,-1) = init_temp(nb_traj)
                                            IF (init_salt(nb_traj) <= s0min(6,-1)) s0min(6,-1) = init_salt(nb_traj)
                                            IF (init_dens(nb_traj) <= s0min(7,-1)) s0min(7,-1) = init_dens(nb_traj)
                                            IF (init_temp(nb_traj) >= s0max(5,-1)) s0max(5,-1) = init_temp(nb_traj)
                                            IF (init_salt(nb_traj) >= s0max(6,-1)) s0max(6,-1) = init_salt(nb_traj)
                                            IF (init_dens(nb_traj) >= s0max(7,-1)) s0max(7,-1) = init_dens(nb_traj)

                                            s0x(5,-1)  = s0x(5,-1)  + init_temp(nb_traj)     * init_trans
                                            s0x(6,-1)  = s0x(6,-1)  + init_salt(nb_traj)     * init_trans
                                            s0x(7,-1)  = s0x(7,-1)  + init_dens(nb_traj)     * init_trans
                                            s0x2(5,-1) = s0x2(5,-1) + init_temp(nb_traj) * init_temp(nb_traj) * init_trans
                                            s0x2(6,-1) = s0x2(6,-1) + init_salt(nb_traj) * init_salt(nb_traj) * init_trans
                                            s0x2(7,-1) = s0x2(7,-1) + init_dens(nb_traj) * init_dens(nb_traj) * init_trans

                                         ENDIF

                                         IF (key_roms.OR.key_symphonie) THEN
                                            prof = &
                                                 fz(gi =  tfi(nb_traj), &
                                                 &  gj =  tfj(nb_traj), &
                                                 &  gk = -tfk(nb_traj), &
                                                 &  il =    l           )
                                         ELSE
                                            prof = fz(gi =  tfi(nb_traj), gj =  tfj(nb_traj), gk=-tfk(nb_traj))
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

                                         s0x(1,-1)  = s0x(1,-1)  + fx(tfi(nb_traj),tfj(nb_traj)) * init_trans

                                         s0x(2,-1)  = s0x(2,-1)  + fy(tfi(nb_traj),tfj(nb_traj)) * init_trans

                                         s0x(3,-1)  = s0x(3,-1)  - prof * init_trans

                                         s0x2(1,-1) = s0x2(1,-1) + fx(tfi(nb_traj),tfj(nb_traj)) * &
                                              fx(tfi(nb_traj),tfj(nb_traj)) * init_trans

                                         s0x2(2,-1) = s0x2(2,-1) + fy(tfi(nb_traj),tfj(nb_traj)) * &
                                              fy(tfi(nb_traj),tfj(nb_traj)) * init_trans

                                         s0x2(3,-1) = s0x2(3,-1) + prof * prof *init_trans

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

                          w = ww(is,js,kw0s,l)

                          !----------------------------------
                          ! test on the sign of the velocity
                          !----------------------------------
                          IF ((ABS(w) > trmin).AND.(REAL(ndir(n),kind=rprec)*w < 0._rprec)) THEN

                             IF (lmt == 1) THEN
                                n1 = 1 + INT( SQRT(ABS(w)/max_transport) )
                                n2 = 1
                             ELSE
                                n1 = 1 + INT( (ABS(w)/max_transport)**(1._rprec/3._rprec) )
                                n2 = n1
                             ENDIF

                             init_trans = ABS(w)/REAL(n1,kind=rprec)/REAL(n1,kind=rprec)/REAL(n2,kind=rprec)

                             DO j1 = 1,n1
                                DO i1 = 1,n1
                                   DO l1 = 1,n2

                                      nb_traj = nb_traj+1

                                      IF (nb_traj > nmax) THEN
                                         WRITE(lun_error,*)'too many particles, maximum is: ',nmax
                                         STOP
                                      ENDIF

                                      tfi(nb_traj) = REAL(i,kind=rprec)   - 1._rprec - &
                                           .5_rprec/REAL(n1,kind=rprec)              + &
                                           REAL(i1,kind=rprec)/REAL(n1,kind=rprec)

                                      tfj(nb_traj) = REAL(j,kind=rprec)   - 1._rprec - &
                                           .5_rprec/REAL(n1,kind=rprec)              + &
                                           REAL(j1,kind=rprec)/REAL(n1,kind=rprec)

                                      tfk(nb_traj) = REAL(kw0,kind=rprec)

                                      tfl(nb_traj) = REAL(l,kind=rprec)   - .5_rprec - &
                                           .5_rprec/REAL(n2,kind=rprec)              + &
                                           REAL(l1,kind=rprec)/REAL(n2,kind=rprec)

                                      ttr(nb_traj) = init_trans

                                      IF (criter0(tfi(nb_traj),tfj(nb_traj),-tfk(nb_traj), &
                                           tfl(nb_traj),i,j,kt0,l)) THEN

                                         trans_total = trans_total + init_trans

                                         isn(-1) = isn(-1) + 1
                                         sn(-1)  =  sn(-1) + init_trans

                                         IF (key_alltracers) THEN

                                            IF (.NOT.key_nointerpolstats) THEN
                                               init_temp(nb_traj) = zinter(tt, &
                                                    tfi(nb_traj),tfj(nb_traj),-tfk(nb_traj),tfl(nb_traj))
                                               init_salt(nb_traj) = zinter(ss, &
                                                    tfi(nb_traj),tfj(nb_traj),-tfk(nb_traj),tfl(nb_traj))
                                               IF (key_approximatesigma) THEN
                                                  init_dens(nb_traj) = zinter(rr, &
                                                       tfi(nb_traj),tfj(nb_traj),-tfk(nb_traj),tfl(nb_traj))
                                               ELSE
                                                  init_dens(nb_traj) = sigma(zsigma, init_salt(nb_traj), init_temp(nb_traj))
                                               ENDIF
                                            ELSE                                    !!NG: bug correction 14/04/2009

                                               init_temp(nb_traj) = tt(is,js,kt0s,1) !!NG: bug correction 14/04/2009
                                               init_salt(nb_traj) = ss(is,js,kt0s,1) !!NG: bug correction 14/04/2009
                                               init_dens(nb_traj) = rr(is,js,kt0s,1) !!NG: bug correction 14/04/2009
                                            ENDIF

                                            IF (init_temp(nb_traj) <= s0min(5,-1)) &
                                                 s0min(5,-1) = init_temp(nb_traj)

                                            IF (init_salt(nb_traj) <= s0min(6,-1)) &
                                                 s0min(6,-1) = init_salt(nb_traj)

                                            IF (init_dens(nb_traj) <= s0min(7,-1)) &
                                                 s0min(7,-1) = init_dens(nb_traj)

                                            IF (init_temp(nb_traj) >= s0max(5,-1)) &
                                                 s0max(5,-1) = init_temp(nb_traj)

                                            IF (init_salt(nb_traj) >= s0max(6,-1)) &
                                                 s0max(6,-1) = init_salt(nb_traj)

                                            IF (init_dens(nb_traj) >= s0max(7,-1)) &
                                                 s0max(7,-1) = init_dens(nb_traj)

                                            s0x(5,-1)  = s0x(5,-1)  + init_temp(nb_traj) * init_trans
                                            s0x(6,-1)  = s0x(6,-1)  + init_salt(nb_traj) * init_trans
                                            s0x(7,-1)  = s0x(7,-1)  + init_dens(nb_traj) * init_trans
                                            s0x2(5,-1) = s0x2(5,-1) + init_temp(nb_traj) * init_temp(nb_traj) * init_trans
                                            s0x2(6,-1) = s0x2(6,-1) + init_salt(nb_traj) * init_salt(nb_traj) * init_trans
                                            s0x2(7,-1) = s0x2(7,-1) + init_dens(nb_traj) * init_dens(nb_traj) * init_trans

                                         ENDIF

                                         IF (key_roms.OR.key_symphonie) THEN
                                            prof= &
                                                 fz(gi=tfi(nb_traj), gj=tfj(nb_traj), gk=-tfk(nb_traj), il=l)
                                         ELSE
                                            prof=fz(gi=tfi(nb_traj), gj=tfj(nb_traj), gk=-tfk(nb_traj))
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

                                         IF (-prof>= s0max(3,-1)) s0max(3,-1) = -prof

                                         s0x(1,-1)  = s0x(1,-1) + fx(tfi(nb_traj),tfj(nb_traj)) * init_trans

                                         s0x(2,-1)  = s0x(2,-1) + fy(tfi(nb_traj),tfj(nb_traj)) * init_trans

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
        ENDDO

     ENDDO

  ELSE
     WRITE(lun_standard,*) '-- posini: bin = ', TRIM(bin), ' --'
  ENDIF

  !- Deallocate dynamic memory array
  CALL sub_memory(-size(it1),'i','it1','posini')
  DEALLOCATE(it1)
  CALL sub_memory(-size(it2),'i','it2','posini')
  DEALLOCATE(it2)
  CALL sub_memory(-size(jt1),'i','jt1','posini')
  DEALLOCATE(jt1)
  CALL sub_memory(-size(jt2),'i','jt2','posini')
  DEALLOCATE(jt2)
  CALL sub_memory(-size(kt1),'i','kt1','posini')
  DEALLOCATE(kt1)
  CALL sub_memory(-size(kt2),'i','kt2','posini')
  DEALLOCATE(kt2)
  CALL sub_memory(-size(num),'i','num','posini')
  DEALLOCATE(num)
  CALL sub_memory(-size(ndir),'i','ndir','posini')
  DEALLOCATE(ndir)

  WRITE(lun_standard,*)''
  WRITE(lun_standard,*)'  ----------------------'
  WRITE(lun_standard,*)'  =  POSINI  >>>> EXIT ='
  WRITE(lun_standard,*)'  ----------------------'

  RETURN
END SUBROUTINE posini
!!***
