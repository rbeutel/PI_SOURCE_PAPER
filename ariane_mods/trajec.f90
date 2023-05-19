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
!!****f* ariane/trajec
!! NAME
!!   trajec (trajec.f90 - Main program - Fortran 90)
!!
!! USAGE
!!   
!!
!! FUNCTION
!! --- Particle trajectories from the analytical computation of 3D
!!     streamlines in a 3D velocity fiels.
!!
!! --- Qualitative or quantitative diagnostics
!!
!! --- Computations done over periods equal to the available sampling
!!     for the velocity field, and over which the velocity is assumed
!!     constant
!!
!! --- The program works with the OPA tensorial formalism
!!
!! --- input: simplified meshmask, 3D TRANSPORT field (velocity x
!!     section), parameters (on various files)
!!
!! --- output: individual trajectories or quantitative outputs
!!
!! --- no-vectorized version
!!
!! --- some comments (yes !) in English
!!
!! --- copyright: Bruno Blanke (Bruno.Blanke@univ-brest.fr)
!! first version: Summer 1992
!! this version: March 2007
!! http://www.ifremer.fr/lpo/blanke/ARIANE
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
!!
!! (see mod_lun)
!! current input files:
!! -------------------

!! current output files:
!! ---------------------
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
!! SOURCE
!!=========================================================================
SUBROUTINE trajec

  USE mod_precision
  USE mod_cst
  USE mod_namelist
  USE mod_input_grid
  USE mod_input_data
  USE mod_output_data
  USE mod_trajec_subs
  USE mod_init_particules
  USE mod_posin
  USE mod_flags
  USE mod_stats
  USE mod_stati
  USE mod_txt
  USE mod_lun
  USE mod_fx
  USE mod_fy
  USE mod_fz
  USE mod_zinter
  USE mod_sigma
  USE mod_quant
  USE mod_reducmem
  USE mod_orca
  USE mod_save_netcdf

  !--------------!
  ! DECLARATIONS !
  !--------------!
  IMPLICIT NONE

  INTEGER(kind=iprec) :: &
       nsect2  =  0    , &   
       ind_traj = 0    , &
       nfinold =  0    , & !
       ibuoy   =  0    , & !
       iswap   =  0    , & !
       iperfl  =  0    , & !
       jperfl  =  0    , & !
       i_mod   =  1    , & ! ntraj / 10 
       inext           , & !
       i, j,n          , & !
       is, js, ks      , & !
       nis, njs, nks   , & !
       i1, i2          , & !
       i1s, i2s        , & !
       j1, j2          , & !
       j1s, j2s        , & !
       k1, k2          , & !
       k1s, k2s        , & !
       lu              , & !
       imin_tmp        , & !
       imax_tmp            !

  INTEGER(kind=iprec), DIMENSION(:), ALLOCATABLE ::  &
       i0    , & !
       j0    , & !
       k0    , & !
       l0    , & !  
       nfin

  REAL(kind = qprec) :: ttt

  REAL(kind=rprec) ::     & !
       subtime =   -1._rprec  , & !
       tfic    =    0._rprec  , & !
       tlap    =    0._rprec  , & !
       tfin    =    0._rprec  , & !
       ft               , & !
       rrr              , & !
       prof             , & !
       age              , & !
       tnext            , & !
       tswap            , & !
       x0, y0, z0       , & !
       u, v, w          , & !
       du, dv, dw       , & !
       tx, ty, tz       , & !     
       tfac             , & !
       t, s, r          , & !
       xt, yt, zt, st   , & !
       trtot1           , & !
       r_lmt            , & !
       dumm

  REAL(kind=rprec), DIMENSION(:), ALLOCATABLE :: & !
       trueage, & !
       agenew , & !
       gi     , & !
       gj     , & !
       gk     , & !
       gl     , & !
       hi     , & !
       hj     , & !
       hk     , & !
       hl     , & !
       tf     , & !
       sf     , & !
       rf     , & !
       init_depth , & ! !! NG: 15_09_2008
       final_depth !! NG: 15_09_2008 

  REAL(kind=rprec), DIMENSION(:,:,:), ALLOCATABLE :: & !
       array_qual_traj

  !-------------------------------------------!
  !- DYNAMIC ALLOCATIONS AND INITIALIZATIONS -!
  !-------------------------------------------!
  CALL sub_posin_alloc(nmax)

  ALLOCATE(i0(1))      ; i0(:)   = 0
  CALL sub_memory(size(i0),'i','i0','trajec')
  ALLOCATE(j0(1))      ; j0(:)   = 0
  CALL sub_memory(size(j0),'i','j0','trajec')
  ALLOCATE(k0(1))      ; k0(:)   = 0
  CALL sub_memory(size(k0),'i','k0','trajec')
  ALLOCATE(l0(1))      ; l0(:)   = 0
  CALL sub_memory(size(l0),'i','l0','trajec')
  ALLOCATE(nfin(nmax)) ; nfin(:) = -1
  CALL sub_memory(size(nfin),'i','nfin','trajec')

  ALLOCATE(trueage(nmax)); trueage(:)= 0._rprec
  CALL sub_memory(size(trueage),'r','trueage','trajec')
  ALLOCATE(agenew(nmax)); agenew(:)= 0._rprec
  CALL sub_memory(size(agenew),'r','agenew','trajec')
  tage(:) = 0._rprec ! define in mod_posin.f90

  ALLOCATE(gi(1))    ; gi(:) = 0._rprec
  CALL sub_memory(size(gi),'r','gi','trajec')
  ALLOCATE(gj(1))    ; gj(:) = 0._rprec
  CALL sub_memory(size(gj),'r','gj','trajec')
  ALLOCATE(gk(1))    ; gk(:) = 0._rprec
  CALL sub_memory(size(gk),'r','gk','trajec')
  ALLOCATE(gl(1))    ; gl(:) = 0._rprec
  CALL sub_memory(size(gl),'r','gl','trajec')
  ALLOCATE(hi(nmax)) ; hi(:) = 0._rprec
  CALL sub_memory(size(hi),'r','hi','trajec')
  ALLOCATE(hj(nmax)) ; hj(:) = 0._rprec
  CALL sub_memory(size(hj),'r','hj','trajec')
  ALLOCATE(hk(nmax)) ; hk(:) = 0._rprec
  CALL sub_memory(size(hk),'r','hk','trajec')
  ALLOCATE(hl(nmax)) ; hl(:) = 0._rprec
  CALL sub_memory(size(hl),'r','hl','trajec')
  ALLOCATE(tf(nmax)) ; tf(:) = 0._rprec
  CALL sub_memory(size(tf),'r','tf','trajec')
  ALLOCATE(sf(nmax)) ; sf(:) = 0._rprec
  CALL sub_memory(size(sf),'r','sf','trajec')
  ALLOCATE(rf(nmax)) ; rf(:) = 0._rprec
  CALL sub_memory(size(rf),'r','rf','trajec')

  !! NG: 15_09_2008
  ALLOCATE(init_depth(nmax))  ; init_depth(:)  = 0._rprec
  CALL sub_memory(size(init_depth),'i','init_depth','trajec')
  ALLOCATE(final_depth(nmax)) ; final_depth(:) = 0._rprec
  CALL sub_memory(size(final_depth),'i','final_depth','trajec')

  !======================!
  !- Read region limits -!
  !======================!
  CALL sub_reducmem_read_reg_lim()

  !==========================!
  ! Open ASCII output files -! (if key_ascii_outputs=True)
  !==========================!
  CALL sub_output_data_open_files()

  !-------------------!
  !- READ INPUT MESH -!
  !-------------------!---------------------------------------------------!
  !-- Allocate and read coordinates and scale factors 
  !-- xx_tt, xx_uu, xx_vv, xx_ff
  !-- yy_tt, yy_uu, yy_vv, yy_ff
  !-- zz_tt, zz_ww
  !-- e2u
  !-- e1v
  !-- e1t, e2t, e3t
  !-- tmask
  !-----------------------------------------------------------------------!
  CALL sub_input_grid()

  !==========================
  !- QUALITATIVE VARIABLES -!
  !=================================================================
  ! sampling time (in seconds) for the available transport field
  !=================================================================
  tfic = tunit   * REAL(ntfic, kind=rprec)

  !=================================================================
  ! sampling time (in seconds) between two successive output positions
  !=================================================================
  tlap = ABS(delta_t) * REAL(frequency, kind=rprec)

  !=================================================================
  ! maximum integration time (im seconds) for individual trajectories
  !=================================================================
  tfin = tlap    * REAL(ABS(nb_output), kind=rprec)

  !==============================================
  ! mask is written (for graphical trajectories)
  !==============================================
  IF ((TRIM(mode) == 'qualitative').AND.(mask)) THEN
     DO j = 1, jmt
        DO i = 1, imt

           !!NG: 20 may 2010 IF (tmask(i,j,1,1) == 0._rprec) THEN
           IF (tmask(i,j,1,1) <= 0.5_rprec) THEN

              IF (key_alltracers) THEN
                 WRITE(lun_traj,7055)0,xx_tt(i,j,1,1),yy_tt(i,j,1,1),0.,0.,0.,0.,0.
              ELSE
                 WRITE(lun_traj,7055)0,xx_tt(i,j,1,1),yy_tt(i,j,1,1),0.,0.
              ENDIF

           ENDIF

        END DO
     END DO
  ENDIF

  !!NG: BEGIN DIRECT DATA ACCESS
  !-------------------!
  !- READ INPUT DATA -!
  !-------------------!---------------------------------------------------!
  !-- Allocate and read (or compute) uu, vv, ww, [tt], [ss] and [rr].
  !-- and store them in the memory (RAM) of the computer.
  !-----------------------------------------------------------------------!
  CALL sub_input_data() 

  !=================================================================
  ! opposite transports (-1) if negative time step
  !=================================================================
  IF (TRIM(forback) == 'backward') THEN
     delta_t = - delta_t
     uu(:,:,:,:) = -uu(:,:,:,:)
     vv(:,:,:,:) = -vv(:,:,:,:)
     ww(:,:,:,:) = -ww(:,:,:,:)
  ENDIF
  !!NG: END DIRECT DATA ACCESS

  !=================================================================
  ! we initialize the *new* particle position 
  !=================================================================
  hi(:) = 0._rprec
  hj(:) = 0._rprec
  hk(:) = 0._rprec
  hl(:) = 0._rprec

  !=======================================!
  !- Read or Compute particule positions -!
  !=======================================!
  CALL sub_init_particules_positions()

  !=======================!
  !- NetCDF data outputs -!
  !=======================!
  !---------------------------!
  !- Open Netcdf output file -!
  !---------------------------!
  CALL sub_save_netcdf_data_init()

  !------------------------------------------------!
  !- Save Initial Positions in netcdf output file -!
  !------------------------------------------------!
  CALL sub_save_netcdf_init_pos()

  !- Deallocate init_temp, init_salt, init_dens (mod_posin) -!
  CALL sub_posin_dealloc_init_tracers()

  IF (TRIM(mode) == 'quantitative') THEN
     CALL sub_quant_alloc(imt, jmt, kmt)
     CALL sub_quant_init() !  initialize a few arrays

     !-------------------------------------------------------------------------
     ! we do not yet know whether a final criterion (criter1) will be satisfied
     !-------------------------------------------------------------------------
     icrit1 = 0

     secname(nsect+1) = 'Criter1'

  ELSEIF (TRIM(mode) == 'qualitative') THEN

     ALLOCATE(array_qual_traj(ntraj,nb_output+1,nstat))
     CALL sub_memory(size(array_qual_traj),'r','array_qual_traj','trajec')

     array_qual_traj(:,:,:)=mask_value

  ELSE
     STOP

  ENDIF

  !==================================================================!
  !==================================================================!
  !          General loop on the number of trajectories.             !
  !==================================================================!
  !==================================================================!

  WRITE(lun_standard,*)''
  WRITE(lun_standard,*)'============================================='
  WRITE(lun_standard,1122) ntraj
  WRITE(lun_standard,*)'============================================='
1122 FORMAT (' = ARIANE IS COMPUTING ',I0.10,' PARTICLES =')

  IF ( ntraj < 20 ) THEN
     i_mod = 1
  ELSEIF ( ntraj < 100) THEN
     i_mod = 10
  ELSEIF ( ntraj < 1000) THEN
     i_mod = 100
  ELSEIF ( ntraj < 10000) THEN
     i_mod = 1000
  ELSEIF ( ntraj < 100000) THEN
     i_mod = 10000
  ELSEIF ( ntraj < 1000000) THEN
     i_mod = 100000
  ELSEIF ( ntraj < 10000000) THEN
     i_mod = 1000000
  ELSEIF ( ntraj < 100000000) THEN
     i_mod = 10000000
  ELSEIF ( ntraj < 1000000000) THEN
     i_mod = 100000000
  ELSE
     i_mod = 1000000000
  ENDIF

  !=====================================================================!
  ! DL0 DL0 DL0 DL0 DL0 DL0 DL0 DL0 DL0 DL0 DL0 DL0 DL0 DL0 DL0 DL0 DL0 !
  !=====================================================================!
  DL0:DO n = 1, ntraj

     !------------!
     !- Comments -!
     !------------!
     !!NG    IF(n == ntraj) THEN
     !!NG      WRITE(lun_standard,'(I0,A)', advance='yes') n, '.'
     !!NG    ELSEIF  ((n == 1).OR.(MOD(n, i_mod) == 0)) THEN
     !!NG      WRITE(lun_standard,'(I0,3A)', advance='no') n, ' - '
     !!NG    ENDIF

     IF  ((n == 1).OR.(MOD(n, i_mod) == 0).OR.(n == ntraj)) THEN
        WRITE(lun_standard,'(10X,I0)') n
     ENDIF

     !====================================================================
     ! Few initializations (transport arrays) [QUANT]
     !====================================================================
     IF (TRIM(mode) == 'quantitative') THEN

        nfinold = 0

        IF (.NOT.key_eco) THEN
           !----------------------!
           !- Take a lot of time -!
           !----------------------!
           CALL sub_quant_inittmp() ! Inititalize arrays
        ENDIF

     ENDIF

     !====================================================================
     ! time index must be within [0.5, lmt+.5[; otherwise time translation
     !NG time index must be within ]0.5, lmt+.5]; otherwise time translation
     !====================================================================
     !NG !!! To work "lmt" must be greater than zero !!! (else infinit loop)
     !NG !!! "lmt" must be tested before to use it.  !!!
     !NG !!! If "lmt" is tested, GOTO could be removed. !!!
100  IF (tfl(n) < 0.5_rprec) THEN
        tfl(n) = tfl(n) + REAL(lmt, kind=rprec)
        GOTO 100
     ENDIF
     IF (tfl(n) >= (REAL(lmt,kind=rprec)+0.5_rprec)) THEN
        tfl(n) = tfl(n) - REAL(lmt, kind=rprec)
        GOTO 100
     ENDIF

     !====================================================================
     ! the program is able to deal with constant-depth particles
     !  (qualitative experiments only) [QUALI]
     !====================================================================
     ibuoy = 0

     IF (tfk(n) <= 0._rprec) THEN

        IF (TRIM(mode) == 'quantitative') THEN 
           WRITE(lun_error,*)''
           WRITE(lun_error,*)'---    ERROR    ---'
           WRITE(lun_error,*)'Isobaric particles allowed in QUALITATIVE experiments ONLY '
           WRITE(lun_error,*)'Particule number n      = ', n
           WRITE(lun_error,*)'                 tfk(n) = ', tfk(n)
           WRITE(lun_error,*)'--- STOP ARIANE ---'
           STOP
        ENDIF

        ibuoy = 1

        tfk(n) = -tfk(n)

     ENDIF

     !====================================================================
     ! time (tfl(n)) and position (VELOCITY grid: fi, tfj(n), tfk(n))
     ! indices used to  detect the relevant embedding temperature
     !  gridcell (i, j, k, l).
     !====================================================================
     i0(1) =  INT(tfi(n), kind = iprec) + 1_iprec
     j0(1) =  INT(tfj(n), kind = iprec) + 1_iprec
     k0(1) =  INT(tfk(n), kind = iprec)
     l0(1) = NINT(tfl(n), kind = iprec)

     !===========================!
     !- Periodicity in OPA-ORCA -!
     !===========================!
     IF (.NOT.(key_roms.OR.key_symphonie)) THEN

        CALL sub_orca_north_pole(i0(1), j0(1))

        CALL sub_orca_east_west_periodic(i0(1))

     ENDIF

     !====================================================================
     ! test of the initial position
     !====================================================================
     CALL sub_reducmem_shift_or_not_ind(i0(1),j0(1),k0(1),is,js,ks)
     write (*,*) i0(1), j0(1), k0(1), is, js, ks, e3t(is,js,ks,1)
     IF (tmask(is,js,ks,1) <= .5_rprec) THEN

        WRITE(lun_error,7000)'   false start', &
             n,i0(1),j0(1),k0(1),l0(1),tfi(n),tfj(n),tfk(n),tfl(n)

        CYCLE DL0

     ENDIF

     !====================================================================
     ! tfi(n), tfj(n), tfk(n) and tfl(n) refer to the initial position
     ! gi, gj, gk and gl refer to the current position
     ! BEWARE: we use a reverse vertical axis for easier computations
     !====================================================================
     gi(1) =  tfi(n)
     gj(1) =  tfj(n)
     gk(1) = -tfk(n)
     gl(1) =  tfl(n)
     !! write(*,*)'1', n, -tfk(n), gk(1)

     !====================================================================
     ! age initialization
     !====================================================================
     age = tage(n)

     !====================================================================
     ! next output instant (for positions in qualitative experiments)
     !====================================================================
     inext = 1
     tnext = tlap

     !====================================================================
     ! next time swap instant (between two samples of velocity field)
     !====================================================================
     !NG "l" is an integer => REAL(l, kind=rprec)
     IF (TRIM(forback) == 'forward') THEN
        tswap = tfic * ABS( tfl(n)-(l0(1)+.5_rprec) ) + tage(n)
     ELSE
        tswap = tfic * ABS( tfl(n)-(l0(1)-.5_rprec) ) + tage(n)
     ENDIF

     !!BBL - march 2007
     !!NG: option to follow an experience interrupted. Age is stored
     !!NG: in the final position file.
     !!NG: Dont USE this option in backward mode !!!
     IF (key_read_age) THEN
        dumm = tswap / tfic
        tswap = NINT(tswap-NINT(dumm)*tfic)+NINT(dumm)*tfic
     ENDIF

     !!NG: to be compatible to sequential version - 09/01/2007
     !!NG: we had also the key_nointerpolstats
     !===================================================================
     ! Compute diagnostics for initial tracers values
     !===================================================================
     CALL sub_reducmem_shift_or_not_ind(i0(1),j0(1),k0(1),is,js,ks)

     IF (key_alltracers) THEN
        IF (key_nointerpolstats) THEN
           tf(n) = tt(is,js,ks,l0(1))
           sf(n) = ss(is,js,ks,l0(1))
           rf(n) = rr(is,js,ks,l0(1))
        ELSE
           tf(n)=zinter(tt,tfi(n),tfj(n),-tfk(n),tfl(n))
           sf(n)=zinter(ss,tfi(n),tfj(n),-tfk(n),tfl(n))
           IF (key_approximatesigma) THEN
              rrr=zinter(rr,tfi(n),tfj(n),-tfk(n),tfl(n))
           ELSE
              rrr=sigma(zsigma,sf(n),tf(n))
           ENDIF
        ENDIF
     ENDIF

     !!NG : 09/01/2007

     !====================================================================
     ! the transport on the initial section is stored in the relevant
     !  transport arrays (quantitative experiment) [QUANT]
     !====================================================================
     IF (TRIM(mode) == 'quantitative') THEN

        IF (key_roms.OR.key_symphonie) THEN
           prof=fz(gi=tfi(n), gj=tfj(n), gk=-tfk(n), il=l0(1))
        ELSE
           prof=fz(gi=tfi(n), gj=tfj(n), gk=-tfk(n))
        ENDIF

        !! NG: 15_09_2008
        !! NG: Here we stored the depth of the particle at the intial position.
        !! NG: Because with ROMS or symphonie the zz_ww in mod_fz is different at
        !! NG: time step.

        init_depth(n)=prof

        CALL sub_reducmem_shift_or_not_ind(NINT(tfi(n)),NINT(tfj(n)),NINT(tfk(n)),nis,njs,nks)

        IF ( (ABS(ANINT(tfi(n))-tfi(n)) <= 1.e-5_rprec) .AND. (uu(nis,js,ks,l0(1)) > 0._rprec) ) THEN
           IF (key_alltracers) THEN
              CALL sub_quant_initside_u(NINT(tfi(n)), j0(1), k0(1), l0(1), &
                   prof, ttr(n), tf(n), sf(n), rrr)
           ELSE
              CALL sub_quant_initside_u(NINT(tfi(n)), j0(1), k0(1), l0(1), &
                   prof, ttr(n))
           ENDIF
        ENDIF

        IF ((ABS(ANINT(tfj(n))-tfj(n)) <= 1.e-5_rprec).AND.(vv(is,njs,ks,l0(1)) > 0._rprec)) THEN
           IF (key_alltracers) THEN
              CALL sub_quant_initside_v(i0(1), NINT(tfj(n)), k0(1), l0(1), &
                   prof, ttr(n), tf(n), sf(n), rrr)
           ELSE
              CALL sub_quant_initside_v(i0(1), NINT(tfj(n)), k0(1), l0(1), &
                   prof, ttr(n))
           ENDIF
        ENDIF

        IF ((ABS(ANINT(tfk(n))-tfk(n)) <= 1.e-5_rprec).AND.(ww(is,js,nks,l0(1)) < 0._rprec)) THEN
           CALL sub_quant_initside_w(i0(1), j0(1), NINT(tfk(n)), l0(1), ttr(n))
        ENDIF


        !====================================================================
        ! the initial position is written on "traj.qt" [QUALI]
        !====================================================================
     ELSEIF (TRIM(mode) == 'qualitative') THEN

        x0 = fx(gi(1),gj(1))
        y0 = fy(gi(1),gj(1))
        IF (key_roms.OR.key_symphonie) THEN
           z0 = fz(gi=gi(1), gj=gj(1), gk=gk(1), il=l0(1))
        ELSE
           z0 = fz(gi=gi(1), gj=gj(1), gk=gk(1))
        ENDIF
        write (*,*) 'initial', x0, y0, z0

        ind_traj = 1

        array_qual_traj(n,ind_traj,1) = x0
        array_qual_traj(n,ind_traj,2) = y0
        array_qual_traj(n,ind_traj,3) = z0
        array_qual_traj(n,ind_traj,4) = 0._rprec

        IF (key_alltracers) THEN

           array_qual_traj(n,ind_traj,5) = tf(n)
           array_qual_traj(n,ind_traj,6) = sf(n)
           array_qual_traj(n,ind_traj,7) = rrr

           IF (key_iU_jV_kW) THEN
              array_qual_traj(n,ind_traj,8)  = gi(1)
              array_qual_traj(n,ind_traj,9)  = gj(1)
              array_qual_traj(n,ind_traj,10) = gk(1)
           END IF

        ELSE

           IF (key_iU_jV_kW) THEN
              array_qual_traj(n,ind_traj,5) = gi(1)
              array_qual_traj(n,ind_traj,6) = gj(1)
              array_qual_traj(n,ind_traj,7) = gk(1)
           END IF

        ENDIF

     ELSE

        STOP

     ENDIF

     iswap = 0 !!NG-bug: -- BUG FIXED by this initialization --
     !! In qualitatif mode it was not possible to compute
     !! more than one particule.

     !====================================================================
     ! entry point (new mesh)
     !====================================================================
     !=====================================================================!
     ! DL1 DL1 DL1 DL1 DL1 DL1 DL1 DL1 DL1 DL1 DL1 DL1 DL1 DL1 DL1 DL1 DL1 !
     !=====================================================================!
     DL1: DO WHILE( (TRIM(mode) == 'quantitative').OR.((age < tfin).AND.(iswap /= 2)) )
        write (*,*) 'age', age, tfin, iswap
        write (*,*) 'gs', gi(1), gj(1), gk(1)
        !====================================================================
        ! linear interpolation for the local velocity
        !====================================================================
        CALL sub_reducmem_shift_or_not_ind(i0(1),j0(1),k0(1),is,js,ks)

        u = uu(is-1,js  ,ks,l0(1)) + &
             (gi(1) - REAL(i0(1)-1, kind=rprec)) * (uu(is,js,ks  ,l0(1)) - uu(is-1,js  ,ks,l0(1)))
        v = vv(is  ,js-1,ks,l0(1)) + &
             (gj(1) - REAL(j0(1)-1, kind=rprec)) * (vv(is,js,ks  ,l0(1)) - vv(is  ,js-1,ks,l0(1)))
        w = ww(is  ,js  ,ks,l0(1)) - &
             (gk(1) + REAL(k0(1)  , kind=rprec)) * (ww(is,js,ks+1,l0(1)) - ww(is  ,js  ,ks,l0(1)))
        write (*,*) u, v, w


        !====================================================================
        ! entry and exit indices (depending on the sign of the velocity)
        !====================================================================
        i1 = i0(1) - 1       !
        i2 = i0(1)           !
        IF (u < 0._rprec) THEN     ! i1 = i - (1 + SIGN(1,u)) * 1/2
           i2 = i0(1) - 1     ! i2 = i - (1 - SIGN(1,u)) * 1/2
           i1 = i0(1)         !      ????? u=0. ?????
        ENDIF                !
        !
        j1 = j0(1) - 1       !
        j2 = j0(1)           !
        IF (v < 0._rprec) THEN     ! j1 = j - (1 + SIGN(1,v)) * 1/2
           j2 = j0(1) - 1     ! j2 = j - (1 - SIGN(1,v)) * 1/2
           j1 = j0(1)         !      ????? v=0. ?????
        ENDIF                !
        !
        k1 = k0(1) + 1       !
        k2 = k0(1)           !
        IF (w < 0._rprec) THEN     ! k1 = k + (1 + SIGN(1,w)) * 1/2
           k2 = k0(1) + 1     ! k2 = k + (1 - SIGN(1,w)) * 1/2
           k1 = k0(1)         !     ????? w=0. ?????
        ENDIF                !
        write (*,*) i1, i2, j1, j2, k1, k2
        !! IF (n==35) write(*,*)'1: k1 and k2:', k1, k2

        !====================================================================
        ! times to get to each edge of the (i, j, k) grid cell:
        ! - infinite, if Uexit * Ulocal < 0
        ! - linear relationship, if Uexit = Ulocal
        ! - logarithmi! formulation, in all other cases
        !  BEWARE: use of log(u)-log(v) instead of log (u/v),
        !   better for numerical results
        !====================================================================
        CALL sub_reducmem_shift_or_not_ind(i0(1) ,j0(1) ,k0(1) ,is ,js ,ks )
        CALL sub_reducmem_shift_or_not_ind(i1,j1,k1,i1s,j1s,k1s)
        CALL sub_reducmem_shift_or_not_ind(i2,j2,k2,i2s,j2s,k2s)

        du = ( uu(i2s,js ,ks ,l0(1)) - uu(i1s,js ,ks ,l0(1)) ) * REAL(i2-i1, kind=rprec)
        dv = ( vv(is ,j2s,ks ,l0(1)) - vv(is ,j1s,ks ,l0(1)) ) * REAL(j2-j1, kind=rprec)
        dw = ( ww(is ,js ,k2s,l0(1)) - ww(is ,js ,k1s,l0(1)) ) * REAL(k1-k2, kind=rprec)
        write (*,*) 'ks', k2s, k1s
        write (*,*) du, dv, dw
        !===============
        ! Scale factors
        !===============
        IF (key_roms.OR.key_symphonie) THEN
           !- ROMS case => E3T is dependent of time -!
           tfac = e3t(is,js,ks,l0(1)) * e1t(is,js,1,1) * e2t(is,js,1,1)

        ELSE
           !- OPA case => E3T is not dependent of time -!
           IF (key_partialsteps) THEN
              tfac = e3t(is,js,ks,1) * e1t(is,js,1,1) * e2t(is,js,1,1)
           ELSE
              tfac = e3t(1,1,ks,1) * e1t(is,js,1,1) * e2t(is,js,1,1)
           ENDIF
           write (*,*) tfac

        ENDIF

        !=============================
        ! Compute time to exit a cell 
        !=============================

        IF (( u * uu(i2s,js,ks,l0(1))) <= 0._rprec) THEN
           tx = 1.e35_rprec
        !!NG: 20 may 2010 ELSEIF (du == 0._rprec) THEN
        ELSEIF (abs(du/u) <= 1.e-11_rprec) THEN
           tx = ( REAL(i2,  kind=rprec) - gi(1) ) / u
        ELSE
           tx = ( LOG(ABS(uu(i2s,js,ks,l0(1)))) - LOG(ABS(u)) ) / du
        ENDIF

        IF ((v * vv(is,j2s,ks,l0(1))) <= 0._rprec) THEN
           ty = 1.e35_rprec
        !!NG: 20 may 2010 ELSEIF (dv == 0._rprec) THEN
        ELSEIF (abs(dv/v) <= 1.e-11_rprec) THEN
           ty = ( REAL(j2, kind=rprec) - gj(1) ) / v
        ELSE
           ty = ( LOG(ABS(vv(is,j2s,ks,l0(1)))) - LOG(ABS(v)) ) / dv

        ENDIF
        write (*,*) w, ww(is,js,k2s,l0(1))
        !!NG: pour BBL
        !!OLD: IF ((ibuoy /= 0).OR.((w*ww(is,js,k2s,l0(1))) <= 0.)) THEN
        IF ((key_2dquant).OR.(ibuoy /= 0).OR.((w*ww(is,js,k2s,l0(1))) <= 0._rprec)) THEN
           tz = 1.e35_rprec
        !!NG: 20 may 2010 ELSEIF (dw == 0._rprec) THEN
        ELSEIF (abs(dw/w) <= 1.e-11_rprec) THEN
           write (*,*) k2, gk(1), w
           tz = -( REAL(k2, kind=rprec) + gk(1) ) / w
        ELSE
           tz = ( LOG(ABS(ww(is,js,k2s,l0(1)))) - LOG(ABS(w)) ) / dw
        ENDIF
        write (*,*) tx, ty, tz

        !====================================================================
        ! the true exit time is given by the shortest times (of tx, ty, tz)
        ! entry and exit indices are then updated
        !====================================================================
        ttt = MIN(tx,ty,tz)
        write (*,*) 'ttt', ttt
        stop

        !!NG: pour BBL 
        IF (.NOT.key_2dquant) THEN
           IF ( (ttt >= 1.e34_qprec) .AND. (ibuoy == 0) ) THEN
              WRITE(lun_error,7000)'      dead end',n,i0(1),j0(1),k0(1),l0(1), &
                   tfi(n),tfj(n),tfk(n),tfl(n)
              CYCLE DL0
           ENDIF
        ENDIF

        IF (ttt < 0._qprec) THEN
           IF (ABS(tfac*ttt) > 1._rprec) THEN
              WRITE(lun_error,7000)' negative time',n,i0(1),j0(1),k0(1),l0(1), &
                   tfi(n),tfj(n),tfk(n),tfl(n)
              CYCLE DL0
           ENDIF
           ttt = 0._qprec
        ENDIF

        !====================================================================
        ! do we need to update the velocity field WITHIN the gridcell ?
        !====================================================================
        iswap = 0

        !====================================================================
        ! update needed
        !====================================================================
        IF ((lmt > 1).AND.((age+ttt*tfac) > tswap)) THEN
           ttt   = ( tswap - age ) / tfac
           iswap = 1
        ENDIF

        !====================================================================
        ! end of the integration [Qualitative]
        !====================================================================
        IF ((TRIM(mode) == 'qualitative').AND.((age+ttt*tfac) > tfin)) THEN
           ttt   = ( tfin - age ) / tfac
           iswap = 2
        ENDIF

        IF (TRIM(forback) == 'forward') THEN
           hl(n) = gl(1) + ttt * tfac / tfic
        ELSE
           hl(n) = gl(1) - ttt * tfac / tfic
        ENDIF

        !====================================================================
        ! case 1: we do not reach the side of the cell (for a given direction)
        ! diagnostic of the final position
        ! - linear formulation, if Uexit = Ulocal
        ! - exponential formulation, in other cases
        !====================================================================
        IF (tx > ttt) THEN
           CALL sub_dont_reachside(hi(n), du, gi(1), u, ttt, i1, i2)
        ENDIF

        IF (ty > ttt) THEN
           CALL sub_dont_reachside(hj(n), dv, gj(1), v, ttt, j1, j2)
        ENDIF

        !!NG: BBL
        IF ((ibuoy == 0).AND.(.NOT.key_2dquant)) THEN
           IF (tz > ttt) THEN
              !! IF (n==35) write(*,*)'-k1, -k2:',-k1, -k2
              CALL sub_dont_reachside(hk(n), dw, gk(1), w, ttt, -k1, -k2)
           ENDIF
        ELSEIF((ibuoy /= 0).OR.(key_2dquant)) THEN
           hk(n) = gk(1) ! IF (ibuoy /= 0) hk(n) = gk(1)
        ENDIF

        !!NG      IF (ibuoy /= 0) hk(n) = gk(1)

        !====================================================================
        ! case 2: we DO reach the side of the cell (for a given direction)
        ! no need for time update within the gridcell
        ! - we record the exiting transport
        ! - we update the positions index (hi, hj(n), hk(n))
        !====================================================================
        IF (tx <= ttt) THEN

           hi(n)=REAL(i2,kind=rprec)

           IF (TRIM(mode) == 'quantitative') &
                CALL sub_quant_reachside_u(i2, j0(1), k0(1), l0(1), hi(n), hj(n), hk(n), hl(n), ttr(n))

           IF (i2 > i1) THEN
              i1 = i2
              i2 = i2 + 1
           ELSE
              i1 = i2
              i2 = i2 - 1
           ENDIF

        ENDIF

        IF (ty <= ttt) THEN

           hj(n) = REAL(j2,kind=rprec)

           IF (TRIM(mode) == 'quantitative') THEN
              CALL sub_quant_reachside_v(      &
                   i0(1), j2, k0(1), l0(1),    &
                   hi(n), hj(n), hk(n), hl(n), &
                   ttr(n))
           ENDIF

           IF (j2 > j1) THEN
              j1 = j2
              j2 = j2 + 1 
           ELSE
              j1 = j2
              j2 = j2 - 1
           ENDIF

        ENDIF

        IF (tz <= ttt) THEN

           IF (TRIM(mode) == 'quantitative') THEN
              CALL sub_quant_reachside_w(i0(1), j0(1), k2, l0(1), ttr(n))
           ENDIF

           !!NG: BBL
           IF (.NOT.key_2dquant) THEN
              IF (ibuoy == 0) THEN
                 hk(n) = -REAL(k2,kind=rprec)

                 IF (k2 > k1) THEN
                    k1 = k2
                    k2 = k2 + 1
                 ELSE
                    k1 = k2
                    k2 = k2 - 1
                 ENDIF

                 !! IF (n==35) write(*,*)'2: k1 and k2:', k1, k2

              ENDIF
           ENDIF

        ENDIF

        !====================================================================
        ! update of the T-grid index of the new embedding gridcell
        !====================================================================
        i0(1) = max0(i1, i2)
        j0(1) = max0(j1, j2)
        k0(1) = min0(k1, k2)
        !!NG-bug:      IF k == 0 =>  free-surface bug

        !===========================!
        !- Periodicity in OPA-ORCA -!
        !===========================!
        IF (.NOT.(key_roms.OR.key_symphonie)) THEN

           !!NG: 8 march 2013 CALL sub_orca_north_pole(i1   , j1   , j0(1)+1)
           !!NG: 8 march 2013 CALL sub_orca_north_pole(i2   , j2   , j0(1)+1)
           !!NG: 8 march 2013 CALL sub_orca_north_pole(hi(n), hj(n), j0(1)+1)
           !!NG: 8 march 2013 CALL sub_orca_north_pole(i0(1), j0(1), j0(1)+1, jperfl)

           CALL sub_orca_north_pole(i1   , j1   , j0(1))
           CALL sub_orca_north_pole(i2   , j2   , j0(1))
           CALL sub_orca_north_pole(hi(n), hj(n), j0(1))
           CALL sub_orca_north_pole(i0(1), j0(1), j0(1), jperfl)


           !!NG: 23/02/2010: BUGGED CALL sub_orca_east_west_periodic(i1   , min0(i1,i2)+1, i0(1))
           !!NG: 23/02/2010: BUGGED CALL sub_orca_east_west_periodic(i2   , min0(i1,i2)+1, i0(1))
           !!NG: 23/02/2010: BUGGED CALL sub_orca_east_west_periodic(hi(n), min0(i1,i2)+1, i0(1))
           !!NG: 23/02/2010: BUGGED CALL sub_orca_east_west_periodic(i0(1), min0(i1,i2)+1, i0(1), iperfl)

           !!NG: 23/02/2010: Replace by:
           imin_tmp = min0(i1,i2)
           imax_tmp = max0(i1,i2)
           CALL sub_orca_east_west_periodic(i1   , imin_tmp+1, imax_tmp)
           CALL sub_orca_east_west_periodic(i2   , imin_tmp+1, imax_tmp)
           CALL sub_orca_east_west_periodic(hi(n), imin_tmp+1, imax_tmp)
           CALL sub_orca_east_west_periodic(i0(1), imin_tmp+1, imax_tmp, iperfl)

        ENDIF

        !====================================================================
        ! age update
        !====================================================================
        agenew(n) = age + ttt * tfac

        IF (iswap == 1) THEN
           hl(n) = NINT(hl(n)+.5_rprec) - .5_rprec
           agenew(n) = tswap
        ENDIF

        IF (iswap == 2)  THEN
           agenew(n) = tfin
        ENDIF


        !====================================================================
        ! we stop the integration if we reach a limit of the domain (it may
        !  happen in QUALITATIVE experiments) [QUALI]
        !====================================================================
        IF (TRIM(mode) == 'qualitative') THEN

           IF (k0(1) == 0) THEN
              WRITE(lun_error,7000)'out of zdomain',n,i0(1),j0(1),k0(1),l0(1), &
                   tfi(n),tfj(n),tfk(n),tfl(n)
              CYCLE DL0
           ENDIF

           !! NG : 10 mai 2010 (BUG)
           IF ( &
                (j0(1) == dims_reg(2,1)).OR.&
                (j0(1) == dims_reg(2,2).AND.(.NOT.key_jfold))) THEN
              WRITE(lun_error,7000)'out of ydomain',n,i0(1),j0(1),k0(1),l0(1), &
                   tfi(n),tfj(n),tfk(n),tfl(n)
              CYCLE DL0
           ENDIF

           IF (.NOT.key_periodic) THEN
              IF ((i0(1) == dims_reg(1,1)).OR.(i0(1) == dims_reg(1,2))) THEN
                 WRITE(lun_error,7000)'out of xdomain',n,i0(1),j0(1),k0(1),l0(1), &
                      tfi(n),tfj(n),tfk(n),tfl(n)
                 CYCLE DL0
              ENDIF
           ENDIF

        ENDIF

        !====================================================================
        ! quantitative analysis
        !
        ! in QUANTITATIVE experiments, the 2D projections of the transport
        !  field associated with 1 particle is summed depending on the final
        !  section [QUANT]
        !====================================================================
        IF (TRIM(mode) == 'quantitative') THEN
           !------------------------------------------------------------------
           ! mtfin informs about the sections used in quantitative experiments
           !------------------------------------------------------------------
           nfin(n) = mtfin(i0(1),j0(1),k0(1))
           !------------------------------------------------
           ! test of a *transparent* hydrological criterion
           !------------------------------------------------
           IF (criter2(tfi(n),tfj(n),tfk(n),tfl(n), &
                hi(n),hj(n),hk(n),hl(n),&
                i0(1),j0(1),k0(1),l0(1))) THEN

              trueage(n) = agenew(n) / tcyc

              IF (key_nointerpolstats) THEN
                 IF (key_alltracers) THEN
                    CALL sub_reducmem_shift_or_not_ind(i0(1),j0(1),k0(1),is,js,ks)
                    t=tt(is,js,ks,l0(1))
                    s=ss(is,js,ks,l0(1))
                    r=rr(is,js,ks,l0(1))
                 ENDIF
              ELSE
                 IF (key_alltracers) THEN
                    t=zinter(tt,hi(n),hj(n),hk(n),hl(n))
                    s=zinter(ss,hi(n),hj(n),hk(n),hl(n))
                    IF (key_approximatesigma) THEN
                       r=zinter(rr,hi(n),hj(n),hk(n),hl(n))
                    ELSE
                       r=sigma(zsigma,s,t)
                    ENDIF
                 ENDIF
              ENDIF

              IF (key_alltracers) THEN
                 WRITE(lun_trans,1002)n,hi(n),hj(n),-hk(n),hl(n),ttr(n),trueage(n),t,s,r
              ELSE
                 WRITE(lun_trans,1002)n,hi(n),hj(n),-hk(n),hl(n),ttr(n),trueage(n)
              ENDIF

1002          FORMAT(i0,3(1x,f0.3),1x,f0.1,1x,f0.3,1x,f0.2, 3(1x,f0.3))
           ENDIF

           !===============================================================
           ! a final criterion may define a "final section" (criter1)
           !===============================================================
           IF ((nfin(n) > 0).OR.(criter1(tfi(n),tfj(n),tfk(n),tfl(n), &
                hi(n),hj(n),hk(n),hl(n), &
                i0(1),j0(1),k0(1),l0(1)))) THEN

              EXIT DL1

              !============================================================
              ! are we on a transparent section ?
              !============================================================
           ELSEIF (nfin(n) < 0) THEN
              IF (nfinold /= nfin(n)) THEN

                 nfinold =      nfin(n)
                 lu      = 90 - nfin(n)
                 idiru   = 'E'
                 idirv   = 'N'
                 idirw   = '^'

                 IF (u < 0) idiru = 'W'
                 IF (v < 0) idirv = 'S'
                 IF (w < 0) idirw = 'v'

                 trueage(n) = agenew(n) / tcyc

                 IF (key_nointerpolstats) THEN
                    IF (key_alltracers) THEN
                       CALL sub_reducmem_shift_or_not_ind(i0(1),j0(1),k0(1),is,js,ks)
                       t=tt(is,js,ks,l0(1))
                       s=ss(is,js,ks,l0(1))
                       r=rr(is,js,ks,l0(1))
                    ENDIF
                 ELSE
                    IF (key_alltracers) THEN
                       t=zinter(tt,hi(n),hj(n),hk(n),hl(n))
                       s=zinter(ss,hi(n),hj(n),hk(n),hl(n))
                       IF (key_approximatesigma) THEN
                          r=zinter(rr,hi(n),hj(n),hk(n),hl(n))
                       ELSE
                          r=sigma(zsigma,s,t)
                       ENDIF
                    ENDIF
                 ENDIF

                 IF (key_alltracers) THEN

                    WRITE(lu,1001)n,hi(n),hj(n),-hk(n),hl(n),ttr(n),trueage(n),idiru, &
                         idirv,idirw, t,s,r

                 ELSE

                    WRITE(lu,1001)n,hi(n),hj(n),-hk(n),hl(n),ttr(n),trueage(n),idiru, &
                         idirv,idirw

                 ENDIF

1001             FORMAT(i0,3(1x,f0.3),1x,f0.1,1x,f0.3,1x,f0.2,3(1x,a2), &
                      3(1x,f0.3))
              ENDIF

              !============================================================
              ! IF (nfin(n) == 0)
              !============================================================
           ELSE

              nfinold=0

           ENDIF


        ENDIF !=> IF (TRIM(mode) == 'quantitative') THEN
        !------------------------------
        ! end of quantitative analysis
        !------------------------------

        !!NG: Here a Bug in Roms mode when iswap /= 0, because l was modified
        !!NG: and in qualitative mode l is used to compute depth (z0) 
        !!NG: in the fz module....
        !!NG: The "time index update of the velocity field if needed" has been
        !!NG: moved after the "qualitative" test...

        !====================================================================
        ! THIS SHOULD NEVER HAPPEN !
        !====================================================================
        CALL sub_reducmem_shift_or_not_ind(i0(1),j0(1),k0(1),is,js,ks)

        IF (tmask(is,js,ks,1) == 0._rprec) THEN
           WRITE(lun_error,7000)'   coast crash',n,i0(1),j0(1),k0(1),l0(1), &
                tfi(n),tfj(n),tfk(n),tfl(n)
           WRITE(lun_coast_crash,*)' '
           WRITE(lun_coast_crash,*)'coast crash for particle # ',n
           WRITE(lun_coast_crash,*)'i j k: ',i0(1),' ',j0(1),' ',k0(1)
           WRITE(lun_coast_crash,*)'is js ks: ',is,' ',js,' ',ks
           WRITE(lun_coast_crash,*)'gi gj gk:',gi(1),' ',gj(1),' ',gk(1)
           WRITE(lun_coast_crash,*)'hi(n) hj(n) hk(n):',hi(n),' ',hj(n),' ',hk(n)
           WRITE(lun_coast_crash,*)'u: ',uu(is-1,js  ,ks  ,l0(1)),uu(is,js,ks,l0(1))
           WRITE(lun_coast_crash,*)'v: ',vv(is  ,js-1,ks  ,l0(1)),vv(is,js,ks,l0(1))
           WRITE(lun_coast_crash,*)'w: ',ww(is  ,js  ,ks+1,l0(1)),ww(is,js,ks,l0(1))
           WRITE(lun_coast_crash,*)'tx ty tz ttt: ',tx,ty,tz,ttt
           STOP
        ENDIF

        !====================================================================
        ! do we need an output of the position (qualitative experiment) ?
        ! (linear interpolation)
        ! - age: previous age
        ! - agenew: current age
        ! - tnext: output time  [QUALI]
        !====================================================================

        IF (TRIM(mode) == 'qualitative') THEN
           write (*,*) 'qualitative'
           DL2:DO WHILE (agenew(n) >= tnext)
              write(*,*) agenew(n), tnext
              ft = (tnext-age) / (agenew(n)-age)
              xt = gi(1) + (hi(n)-gi(1)) * ft
              yt = gj(1) + (hj(n)-gj(1)) * ft
              zt = gk(1) + (hk(n)-gk(1)) * ft

              IF (TRIM(forback) == 'forward') THEN
                 st = gl(1) + (tnext-age) / tfic
              ELSE
                 st = gl(1) - (tnext-age) / tfic
              ENDIF

              IF (key_periodic) THEN
                 IF (iperfl == 1) THEN
                    xt = gi(1) + (hi(n)-REAL(imt-2,kind=rprec)-gi(1)) * ft
                 ENDIF
                 IF (iperfl == 2) THEN
                    xt = gi(1) + (hi(n)+REAL(imt-2,kind=rprec)-gi(1)) * ft
                 ENDIF
              ENDIF

              IF (key_jfold) THEN
                 !-----------------------------------
                 ! extrapolation not yet implemented
                 !-----------------------------------
                 IF (jperfl == 1) THEN
                    xt = gi(1)
                    yt = gj(1)
                    zt = gk(1)
                 ENDIF
              ENDIF
              !-------------------------
              ! position output on file
              !-------------------------
              x0 = fx(xt,yt)
              y0 = fy(xt,yt)

              IF (key_roms.OR.key_symphonie) THEN
                 z0 = fz(gi=xt, gj=yt, gk=zt, il=l0(1))
              ELSE
                 z0 = fz(gi=xt, gj=yt, gk=zt)
              ENDIF

              IF (tnext <= tfin) THEN

                 ind_traj = ind_traj + 1

                 array_qual_traj(n,ind_traj,1) = x0
                 array_qual_traj(n,ind_traj,2) = y0
                 array_qual_traj(n,ind_traj,3) = z0
                 array_qual_traj(n,ind_traj,4) = tnext/tcyc

                 IF (key_alltracers) THEN

                    tf(n) = zinter(tt,xt,yt,zt,st)
                    sf(n) = zinter(ss,xt,yt,zt,st)

                    IF (key_approximatesigma) THEN

                       rrr = zinter(rr,xt,yt,zt,st)

                    ELSE

                       rrr = sigma(zsigma,sf(n),tf(n))

                    ENDIF

                    array_qual_traj(n,ind_traj,5) = tf(n)
                    array_qual_traj(n,ind_traj,6) = sf(n)
                    array_qual_traj(n,ind_traj,7) = rrr

                    IF (key_iU_jV_kW) THEN
                       array_qual_traj(n,ind_traj,8)  = xt
                       array_qual_traj(n,ind_traj,9)  = yt
                       array_qual_traj(n,ind_traj,10) = zt
                    END IF

                 ELSE

                    IF (key_iU_jV_kW) THEN
                       array_qual_traj(n,ind_traj,5) = xt
                       array_qual_traj(n,ind_traj,6) = yt
                       array_qual_traj(n,ind_traj,7) = zt
                    END IF

                 ENDIF

              ENDIF

              inext = inext + 1
              tnext = tlap * REAL(inext, kind = rprec)

           ENDDO  DL2

        ENDIF

        !====================================================================
        ! time index update of the velocity field if needed
        !====================================================================
        IF (iswap /= 0) THEN
           IF (TRIM(forback) == 'forward')  THEN 
              l0(1) = l0(1) + 1
           ELSE
              l0(1) = l0(1) - 1
           ENDIF
           IF (l0(1) < 1)   l0(1) = l0(1) + lmt
           IF (l0(1) > lmt) l0(1) = l0(1) - lmt
           tswap = tswap + tfic
        ENDIF
        write (*,*) 'l0', l0(1)

        !====================================================================
        ! swap for space and time index positions
        !====================================================================
        !! IF (n==35) write(*,*)'2', n, gk(1)

        gi(1)  = hi(n)
        gj(1)  = hj(n)
        gk(1)  = hk(n)
        gl(1)  = hl(n)
        age = agenew(n)

        !! IF (n==35) write(*,*)'3', n, gk(1)


        !=============!
        ! game over ? !
        !=============!

     ENDDO DL1
     !=====================================================================!
     ! DL1 DL1 DL1 DL1 DL1 DL1 DL1 DL1 DL1 DL1 DL1 DL1 DL1 DL1 DL1 DL1 DL1 !
     !=====================================================================!

     IF (TRIM(mode) == 'quantitative') THEN
        !==================================================================
        ! a final criterion uses "1 + nsect" as the section index
        ! we flag the existence of such a criterion
        !==================================================================
        IF (nfin(n) <= 0) THEN
           nfin(n) = nsect + 1
           icrit1  = 1
        ENDIF

        !==================================================================
        ! we distinguish meandering and recirculating particles (criter0)
        !==================================================================
        IF ((TRIM(bin) == 'nobin').AND.(nfin(n) == 1).AND. &
             criter0(hi(n),hj(n),hk(n),hl(n),i0(1),j0(1),k0(1),l0(1))) nfin(n)=0

        !==================================================================
        ! hereafter,
        ! nfin=0, if the initial section (with criter0) has been reached
        ! nfin=1, if the initial section (no criter0) or a final section has
        !  been reached
        ! nfin=1+nsect, if a final criterion (criter1) is satisfied
        !==================================================================
        sectrans(nfin(n)) = sectrans(nfin(n)) + ttr(n)

        !==================================================================
        ! there IS NOW storage of transport projections for "criter1" cases
        !  (but never try to compute 2D streamfunctions with these fields !)
        !==================================================================
        IF (.NOT.key_eco) THEN
           CALL sub_quant_store_transp(nfin(n))
        ENDIF

        !==================================================================
        ! statistics...
        !==================================================================
        CALL sub_quant_statistics(nfin(n), ttr(n), agenew(n), tcyc, &
             hi(n), hj(n), hk(n), hl(n), tfi(n), tfj(n), tfk(n), tfl(n), &
             i0(1), j0(1), k0(1), l0(1), tf(n), sf(n), rf(n), trueage(n))

        !! NG: 15_09_2008
        !! NG: Here we store the final depth which depend of time for ROMS
        !! NG: and symphonie.
        IF (key_roms.OR.key_symphonie) THEN
           prof=fz(gi=hi(n), gj=hj(n), gk=hk(n), il=1)
        ELSE
           prof=fz(gi=hi(n), gj=hj(n), gk=hk(n))
        ENDIF
        final_depth(n)=prof

        WRITE(lun_output,7000)secname(nfin(n)),n, &
             i0(1),j0(1),k0(1),l0(1), &
             tfi(n),tfj(n),tfk(n),tfl(n)
7000    FORMAT(a14,', #',i0,' en:', 1x,i0,1x,i0,1x,i0,1x,i0,',  init.=',4(1x,f0.2),': ',f0.2)

     ENDIF

     !================
     ! shoot again...
     !================

  ENDDO DL0 !=> DO n=1,ntraj
  !=====================================================================!
  ! DL0 DL0 DL0 DL0 DL0 DL0 DL0 DL0 DL0 DL0 DL0 DL0 DL0 DL0 DL0 DL0 DL0 !
  !=====================================================================!
  !======================================================================!
  !END END END - General loop on the number of trajectories - END END END!
  !======================================================================!

  WRITE(lun_standard,*)''
  WRITE(lun_standard,*)'-----------------------'
  WRITE(lun_standard,*)'= Writing Output Data ='
  WRITE(lun_standard,*)'-----------------------'

  !!===============================================================================!!
  !!============================== QUALITATIVE MODE ===============================!!
  !!===============================================================================!!
  IF (TRIM(mode) == 'qualitative' ) THEN

     !----------------!
     !- Trajectories -!
     !----------------!
     !- Netcdf -!
     !----------!

     IF (key_alltracers.AND.key_iU_jV_kW) THEN

        CALL sub_save_netcdf_trajectories(                   &
             array_qual_traj(:,:,1), array_qual_traj(:,:,2), &
             array_qual_traj(:,:,3), array_qual_traj(:,:,4), &
             array_qual_traj(:,:,5), array_qual_traj(:,:,6), &
             array_qual_traj(:,:,7), array_qual_traj(:,:,8), &
             array_qual_traj(:,:,9), array_qual_traj(:,:,10) )

     ELSEIF (key_alltracers.AND.(.NOT.key_iU_jV_kW)) THEN

        CALL sub_save_netcdf_trajectories(                   &
             array_qual_traj(:,:,1), array_qual_traj(:,:,2), &
             array_qual_traj(:,:,3), array_qual_traj(:,:,4), &
             traj_temp=array_qual_traj(:,:,5), &
             traj_salt=array_qual_traj(:,:,6), &
             traj_dens=array_qual_traj(:,:,7))
 
     ELSEIF ((.NOT.key_alltracers).AND.key_iU_jV_kW) THEN
        CALL sub_save_netcdf_trajectories(                   &
             array_qual_traj(:,:,1), array_qual_traj(:,:,2), &
             array_qual_traj(:,:,3), array_qual_traj(:,:,4), &
             traj_iU=array_qual_traj(:,:,5), &
             traj_jV=array_qual_traj(:,:,6), &
             traj_kW=array_qual_traj(:,:,7))  
     ELSE
        CALL sub_save_netcdf_trajectories(                   &
             array_qual_traj(:,:,1), array_qual_traj(:,:,2), &
             array_qual_traj(:,:,3), array_qual_traj(:,:,4)  )
     ENDIF

!!NG : 2 march 2010
!!$     IF (key_alltracers) THEN
!!$
!!$        CALL sub_save_netcdf_trajectories(                   &
!!$             array_qual_traj(:,:,1), array_qual_traj(:,:,2), &
!!$             array_qual_traj(:,:,3), array_qual_traj(:,:,4), &
!!$             array_qual_traj(:,:,5), array_qual_traj(:,:,6), &
!!$             array_qual_traj(:,:,7)                          )
!!$     ELSE
!!$        CALL sub_save_netcdf_trajectories(                   &
!!$             array_qual_traj(:,:,1), array_qual_traj(:,:,2), &
!!$             array_qual_traj(:,:,3), array_qual_traj(:,:,4)  )
!!$     ENDIF

     !---------!
     !- ASCII -!
     !---------!
     IF (key_ascii_outputs) THEN

        DO n = 1, ntraj
           DO i = 1, nb_output+1
              WRITE(lun_traj,7055)n,array_qual_traj(n,i,1:nstat)
           ENDDO
        ENDDO

     ENDIF

!!NG: 2 juin 2010: 7055 FORMAT (i0,3(1x,f0.5),1x,f0.5,1x,3(1x,f0.5))
7055 FORMAT (i0,3(1x,f0.15),1x,f0.5,1x,3(1x,f0.5))

     !!================================================================================!!
     !!============================== QUANTITATIVE MODE ===============================!!
     !!================================================================================!!
  ELSE
     !-------------------!
     !- Final positions -!
     !-------------------! 
     !---------!
     !- ASCII -!
     !---------!
     IF (key_ascii_outputs) THEN
        DO n = 1, ntraj



           IF (key_alltracers) THEN

              WRITE(lun_fin_pos,7059)hi(n),hj(n),-hk(n),hl(n), &
                   ttr(n),trueage(n),nfin(n), tf(n),sf(n),rf(n)

           ELSE

              WRITE(lun_fin_pos,7059)hi(n),hj(n),-hk(n),hl(n), &
                   ttr(n),trueage(n),nfin(n)

           ENDIF

        ENDDO

     ENDIF

7059 FORMAT(3(1x,f0.3),1x,f0.3,1x,f0.3,1x,f0.3,1x,i0,3(1x,f0.3))

     !----------------------------!
     !- statistics file in ASCII -!
     !----------------------------!
     WRITE(lun_stats,*)' '

     IF (key_unitm3) THEN
        WRITE(lun_stats,*)'total transport (in m3/s): ', trtot/(1+lmax-lmin), &
             ' ( x [1+lmax-lmin] =' ,trtot, ')'
        WRITE(lun_stats,*)'max_transport (in m3/s)  : ', max_transport
     ELSE
        WRITE(lun_stats,*)'total transport (in Sv)  : ', (trtot/(1+lmax-lmin))/1.e6_rprec , &
             ' ( x [1+lmax-lmin] =' ,trtot/1.e6_rprec, ')'
        WRITE(lun_stats,*)'max_transport (in Sv)    : ', max_transport / 1.e6_rprec
     ENDIF
     WRITE(lun_stats,*)'# particles              : '  , ntraj 

     IF ((TRIM(bin) == 'nobin').AND.(isn(-1) > 0)) THEN
        !------------------------------------------
        ! we mask uninitialized min. and max. ages
        !------------------------------------------
        DO i = 1, nstat
           IF (ABS(s1min(i,-1)) > 1.e10_rprec) s1min(i,-1) = 0._rprec
           IF (ABS(s1max(i,-1)) > 1.e10_rprec) s1max(i,-1) = 0._rprec
           IF (ABS(s0min(i,-1)) > 1.e10_rprec) s0min(i,-1) = 0._rprec
           IF (ABS(s0max(i,-1)) > 1.e10_rprec) s0max(i,-1) = 0._rprec
        END DO

        WRITE(lun_stats,*)   ' '
        WRITE(lun_stats,7057)'initial state ',' ','#',isn(-1)
        WRITE(lun_stats,7157)' stats. for:   ',(chpstat(i),i=1,nstat)
        WRITE(lun_stats,7257)' min: '        ,(s0min(i,-1),i=1,nstat)
        WRITE(lun_stats,7257)' max: '        ,(s0max(i,-1),i=1,nstat)
        WRITE(lun_stats,7257)'mean: '        ,(s0x(i,-1)/sn(-1),i=1,nstat)
        WRITE(lun_stats,7257)'std. dev.: '   ,&
             (SQRT(ABS((s0x2(i,-1)-s0x(i,-1)*s0x(i,-1)/sn(-1))/(sn(-1)-1._rprec))),i=1,nstat)
     ENDIF

     IF (TRIM(bin) /= 'nobin') THEN
        WRITE(lun_stats,*)' '
        WRITE(lun_stats,*)'no stats for initial state available since'
        WRITE(lun_stats,*)'   automatic positioning was by-passed'
     ENDIF

     WRITE(lun_stats,*)' '

     !================================
     ! 2D flow projection and z stats
     !================================
     subtime = REAL(lmt, kind=rprec) / REAL(1+lmax-lmin, kind =rprec)

     !NG Performances (sub_quant_store_transp => * r_lmt)
     r_lmt =  1._rprec / REAL(lmt, kind=rprec)

     !!NG    IF (key_eco) THEN
     !!NG      WRITE(lun_xy_zonal)((subtime*uxy(i,j,1)*r_lmt,i=1,imt),j=1,jmt)
     !!NG      WRITE(lun_xy_merid)((subtime*vxy(i,j,1)*r_lmt,i=1,imt),j=1,jmt)
     !!NG      WRITE(lun_yz_merid)((subtime*vyz(j,k,1)*r_lmt,j=1,jmt),k=1,kmt)
     !!NG      WRITE(lun_yz_verti)((subtime*wyz(j,k,1)*r_lmt,j=1,jmt),k=1,kmt)
     !!NG      WRITE(lun_xz_zonal)((subtime*uxz(i,k,1)*r_lmt,i=1,imt),k=1,kmt)
     !!NG      WRITE(lun_xz_verti)((subtime*wxz(i,k,1)*r_lmt,i=1,imt),k=1,kmt)
     !!NG      WRITE(lun_xy_stats)((uh(i,j,1)*r_lmt         ,i=1,imt),j=1,jmt)
     !!NG      WRITE(lun_xy_stats)((vh(i,j,1)*r_lmt         ,i=1,imt),j=1,jmt)
     !!NG      WRITE(lun_xy_stats)((zuh(i,j,1)*r_lmt        ,i=1,imt),j=1,jmt)
     !!NG      WRITE(lun_xy_stats)((zvh(i,j,1)*r_lmt        ,i=1,imt),j=1,jmt)
     !!NG      WRITE(lun_xy_stats)((z2uh(i,j,1)*r_lmt       ,i=1,imt),j=1,jmt)
     !!NG      WRITE(lun_xy_stats)((z2vh(i,j,1)*r_lmt       ,i=1,imt),j=1,jmt)
     !!NG
     !!NG      IF (key_alltracers) THEN
     !!NG        WRITE(lun_xy_stats)((tuh(i,j,1)*r_lmt        ,i=1,imt),j=1,jmt)
     !!NG        WRITE(lun_xy_stats)((tvh(i,j,1)*r_lmt        ,i=1,imt),j=1,jmt)
     !!NG        WRITE(lun_xy_stats)((t2uh(i,j,1)*r_lmt       ,i=1,imt),j=1,jmt)
     !!NG        WRITE(lun_xy_stats)((t2vh(i,j,1)*r_lmt       ,i=1,imt),j=1,jmt)
     !!NG        WRITE(lun_xy_stats)((suh(i,j,1)*r_lmt        ,i=1,imt),j=1,jmt)
     !!NG        WRITE(lun_xy_stats)((svh(i,j,1)*r_lmt        ,i=1,imt),j=1,jmt)
     !!NG        WRITE(lun_xy_stats)((s2uh(i,j,1)*r_lmt       ,i=1,imt),j=1,jmt)
     !!NG        WRITE(lun_xy_stats)((s2vh(i,j,1)*r_lmt       ,i=1,imt),j=1,jmt)
     !!NG        WRITE(lun_xy_stats)((ruh(i,j,1)*r_lmt        ,i=1,imt),j=1,jmt)
     !!NG        WRITE(lun_xy_stats)((rvh(i,j,1)*r_lmt        ,i=1,imt),j=1,jmt)
     !!NG        WRITE(lun_xy_stats)((r2uh(i,j,1)*r_lmt       ,i=1,imt),j=1,jmt)
     !!NG        WRITE(lun_xy_stats)((r2vh(i,j,1)*r_lmt       ,i=1,imt),j=1,jmt)
     !!NG      ENDIF
     !!NG
     !!NG    ELSE
     !!NG
     !!NG      WRITE(lun_xy_zonal)(((subtime*uxy(i,j,n)*r_lmt,i=1,imt),j=1,jmt),n=0,nsect2)
     !!NG      WRITE(lun_xy_merid)(((subtime*vxy(i,j,n)*r_lmt,i=1,imt),j=1,jmt),n=0,nsect2)
     !!NG      WRITE(lun_yz_merid)(((subtime*vyz(j,k,n)*r_lmt,j=1,jmt),k=1,kmt),n=0,nsect2)
     !!NG      WRITE(lun_yz_verti)(((subtime*wyz(j,k,n)*r_lmt,j=1,jmt),k=1,kmt),n=0,nsect2)
     !!NG      WRITE(lun_xz_zonal)(((subtime*uxz(i,k,n)*r_lmt,i=1,imt),k=1,kmt),n=0,nsect2)
     !!NG      WRITE(lun_xz_verti)(((subtime*wxz(i,k,n)*r_lmt,i=1,imt),k=1,kmt),n=0,nsect2)
     !!NG      WRITE(lun_xy_stats)(((uh(i,j,n)*r_lmt         ,i=1,imt),j=1,jmt),n=0,nsect2)
     !!NG      WRITE(lun_xy_stats)(((vh(i,j,n)*r_lmt         ,i=1,imt),j=1,jmt),n=0,nsect2)
     !!NG      WRITE(lun_xy_stats)(((zuh(i,j,n)*r_lmt        ,i=1,imt),j=1,jmt),n=0,nsect2)
     !!NG      WRITE(lun_xy_stats)(((zvh(i,j,n)*r_lmt        ,i=1,imt),j=1,jmt),n=0,nsect2)
     !!NG      WRITE(lun_xy_stats)(((z2uh(i,j,n)*r_lmt       ,i=1,imt),j=1,jmt),n=0,nsect2)
     !!NG      WRITE(lun_xy_stats)(((z2vh(i,j,n)*r_lmt       ,i=1,imt),j=1,jmt),n=0,nsect2)
     !!NG
     !!NG      IF (key_alltracers) THEN
     !!NG        WRITE(lun_xy_stats)(((tuh(i,j,n)*r_lmt        ,i=1,imt),j=1,jmt),n=0,nsect2)
     !!NG        WRITE(lun_xy_stats)(((tvh(i,j,n)*r_lmt        ,i=1,imt),j=1,jmt),n=0,nsect2)
     !!NG        WRITE(lun_xy_stats)(((t2uh(i,j,n)*r_lmt       ,i=1,imt),j=1,jmt),n=0,nsect2)
     !!NG        WRITE(lun_xy_stats)(((t2vh(i,j,n)*r_lmt       ,i=1,imt),j=1,jmt),n=0,nsect2)
     !!NG        WRITE(lun_xy_stats)(((suh(i,j,n)*r_lmt        ,i=1,imt),j=1,jmt),n=0,nsect2)
     !!NG        WRITE(lun_xy_stats)(((svh(i,j,n)*r_lmt        ,i=1,imt),j=1,jmt),n=0,nsect2)
     !!NG        WRITE(lun_xy_stats)(((s2uh(i,j,n)*r_lmt       ,i=1,imt),j=1,jmt),n=0,nsect2)
     !!NG        WRITE(lun_xy_stats)(((s2vh(i,j,n)*r_lmt       ,i=1,imt),j=1,jmt),n=0,nsect2)
     !!NG        WRITE(lun_xy_stats)(((ruh(i,j,n)*r_lmt        ,i=1,imt),j=1,jmt),n=0,nsect2)
     !!NG        WRITE(lun_xy_stats)(((rvh(i,j,n)*r_lmt        ,i=1,imt),j=1,jmt),n=0,nsect2)
     !!NG        WRITE(lun_xy_stats)(((r2uh(i,j,n)*r_lmt       ,i=1,imt),j=1,jmt),n=0,nsect2)
     !!NG        WRITE(lun_xy_stats)(((r2vh(i,j,n)*r_lmt       ,i=1,imt),j=1,jmt),n=0,nsect2)
     !!NG      ENDIF
     !!NG
     !!NG    ENDIF


     uxy(:,:,:)  = uxy(:,:,:) * subtime * r_lmt
     vxy(:,:,:)  = vxy(:,:,:) * subtime * r_lmt
     vyz(:,:,:)  = vyz(:,:,:) * subtime * r_lmt
     wyz(:,:,:)  = wyz(:,:,:) * subtime * r_lmt
     uxz(:,:,:)  = uxz(:,:,:) * subtime * r_lmt
     wxz(:,:,:)  = wxz(:,:,:) * subtime * r_lmt
     uh(:,:,:)   = uh(:,:,:)            * r_lmt
     vh(:,:,:)   = vh(:,:,:)            * r_lmt
     zuh(:,:,:)  = zuh(:,:,:)           * r_lmt
     zvh(:,:,:)  = zvh(:,:,:)           * r_lmt
     z2uh(:,:,:) = z2uh(:,:,:)          * r_lmt
     z2vh(:,:,:) = z2vh(:,:,:)          * r_lmt

     IF (key_alltracers) THEN
        tuh(:,:,:)  = tuh(:,:,:)         * r_lmt
        tvh(:,:,:)  = tvh(:,:,:)         * r_lmt
        t2uh(:,:,:) = t2uh(:,:,:)        * r_lmt
        t2vh(:,:,:) = t2vh(:,:,:)        * r_lmt
        suh(:,:,:)  = suh(:,:,:)         * r_lmt
        svh(:,:,:)  = svh(:,:,:)         * r_lmt
        s2uh(:,:,:) = s2uh(:,:,:)        * r_lmt
        s2vh(:,:,:) = s2vh(:,:,:)        * r_lmt
        ruh(:,:,:)  = ruh(:,:,:)         * r_lmt
        rvh(:,:,:)  = rvh(:,:,:)         * r_lmt
        r2uh(:,:,:) = r2uh(:,:,:)        * r_lmt
        r2vh(:,:,:) = r2vh(:,:,:)        * r_lmt
     ENDIF

     !!---------------------------!
     !!= Quantitative statistics =!
     !!---------------------------!
     trtot1 = 0._rprec

     nsect2 = nsect
     IF (icrit1 == 1) THEN
        nsect2 = 1 + nsect
     ENDIF

     DO n = 0, nsect2
        IF (key_unitm3) THEN
           WRITE(lun_stats,7357)secname(n), subtime*sectrans(n)*r_lmt, n
        ELSE
           WRITE(lun_stats,7357)secname(n), subtime*sectrans(n)/1.e6_rprec*r_lmt, n
        ENDIF
        trtot1 = trtot1 + sectrans(n)
     END DO

     IF (key_unitm3) THEN
        WRITE(lun_stats,7357)'total'      ,subtime*trtot*r_lmt
        WRITE(lun_stats,7357)'except mnds',subtime*(trtot-sectrans(0))*r_lmt
        WRITE(lun_stats,7357)'lost'       ,subtime*(trtot-trtot1)*r_lmt
     ELSE
        WRITE(lun_stats,7357)'total'      ,subtime*trtot/1.e6_rprec*r_lmt
        WRITE(lun_stats,7357)'except mnds',subtime*(trtot-sectrans(0))/1.e6_rprec*r_lmt
        WRITE(lun_stats,7357)'lost'       ,subtime*(trtot-trtot1)/1.e6_rprec*r_lmt
     ENDIF

7357 FORMAT(a15,1x,f0.4,1x,i0)

     DO n=0,nsect2
        IF (isn(n) > 0) THEN
           WRITE(lun_stats,*)' '
           WRITE(lun_stats,7057)'final state ',TRIM(secname(n)),' # ',isn(n)
           WRITE(lun_stats,7157)'stats. ini: ',(chpstat(i),i=1,nstat)
           WRITE(lun_stats,7257)' min: ',(s0min(i,n),i=1,nstat)
           WRITE(lun_stats,7257)' max: ',(s0max(i,n),i=1,nstat)
           WRITE(lun_stats,7257)'mean: ',(s0x(i,n)/sn(n),i=1,nstat)
           WRITE(lun_stats,7257)'std. dev.: ',( SQRT( &
                ABS((s0x2(i,n)-s0x(i,n)*s0x(i,n)/sn(n))/(sn(n)-1._rprec))),i=1,nstat)
           WRITE(lun_stats,7157)'stats. fin: ',(chpstat(i),i=1,nstat)
           WRITE(lun_stats,7257)' min: ',(s1min(i,n),i=1,nstat)
           WRITE(lun_stats,7257)' max: ',(s1max(i,n),i=1,nstat)
           WRITE(lun_stats,7257)'mean: ',(s1x(i,n)/sn(n),i=1,nstat)
           WRITE(lun_stats,7257)'std. dev.: ',( SQRT( &
                ABS((s1x2(i,n)-s1x(i,n)*s1x(i,n)/sn(n))/(sn(n)-1._rprec))),i=1,nstat)
        ENDIF
     END DO

     !----------------------------------------------------------------------!
     !- Define statistic variables and their attributes in quantative mode -!
     !----------------------------------------------------------------------!
     CALL sub_save_netcdf_data_init_stats()

     !------------------------------------!
     !- Save statistics in Netcdf Format -!
     !------------------------------------!
     Call sub_save_netcdf_stats()

     !----------------------------!
     !- Close Netcdf output file -!
     !----------------------------!
     WRITE(lun_standard,*)''
     WRITE(lun_standard,*)'--------------------------------'
     WRITE(lun_standard,*)'= Close Statistics NetCDF File ='
     WRITE(lun_standard,*)'--------------------------------'
     CALL sub_save_netcdf_data_close(ncid_data_stats)

  ENDIF

7057 FORMAT(a13,a15,a2,i7)
7157 FORMAT(a13,7a10)
7257 FORMAT(a13,7f10.3)

  !=====================================!
  !- Binary storage of final positions -!
  !=====================================!
  IF (key_ascii_outputs) THEN
     DO n = 1, ntraj
        WRITE(lun_final)hi(n),hj(n),-hk(n),hl(n),ttr(n),agenew(n)
        WRITE(lun_init)tfi(n),tfj(n),tfk(n),tfl(n),ttr(n),tage(n)
     ENDDO
  ENDIF

  !! NG: 15_09_2008
  IF (TRIM(mode) == 'quantitative') THEN 
     !-----------------!
     !- Initial Depth -!
     !-----------------!
     !- Netcdf -!
     !----------!

     CALL sub_save_netcdf_init_depth(init_depth(:))

  ENDIF

  !-------------------!
  !- Final positions -!
  !-------------------! 
  !- Netcdf -!
  !----------!
  IF (key_alltracers) THEN
     CALL sub_save_netcdf_final_pos( &
          hi(:),hj(:),-hk(:),hl(:) , &
          agenew(:), nfin(:)       , &
          final_depth(:)           , &
          tf(:), sf(:), rf(:)        )
  ELSE
     CALL sub_save_netcdf_final_pos( &
          hi(:),hj(:),-hk(:),hl(:) , &
          agenew(:), nfin(:),final_depth(1:ntraj) )
  ENDIF

  !-----------------------------!
  !- Close Netcdf output files -!
  !-----------------------------!
  WRITE(lun_standard,*)''
  WRITE(lun_standard,*)'-------------------------------'
  WRITE(lun_standard,*)'= Close Positions NetCDF File ='
  WRITE(lun_standard,*)'-------------------------------'
  CALL sub_save_netcdf_data_close(ncid_data_pos)

  !=====================!
  !- DEALLOCATE memory -!
  !=====================!

  WRITE(lun_standard,*)''
  WRITE(lun_standard,*)'---------------------'
  WRITE(lun_standard,*)'= Deallocate Memory ='
  WRITE(lun_standard,*)'---------------------'


  IF (ALLOCATED(i0)) THEN
     CALL sub_memory(-size(i0),'i','i0','trajec_seq')
     DEALLOCATE(i0)
  END IF
  IF (ALLOCATED(j0)) THEN
     CALL sub_memory(-size(j0),'i','j0','trajec_seq')
     DEALLOCATE(j0)
  END IF
  IF (ALLOCATED(k0)) THEN
     CALL sub_memory(-size(k0),'i','k0','trajec_seq')
     DEALLOCATE(k0)
  END IF
  IF (ALLOCATED(l0)) THEN
     CALL sub_memory(-size(l0),'i','l0','trajec_seq')
     DEALLOCATE(l0)
  END IF
  IF (ALLOCATED(nfin)) THEN
     CALL sub_memory(-size(nfin),'i','nfin','trajec_seq')
     DEALLOCATE(nfin)
  END IF

  IF (ALLOCATED(agenew)) THEN
     CALL sub_memory(-size(agenew),'r','agenew','trajec_seq')
     DEALLOCATE(agenew)
  END IF
  IF (ALLOCATED(trueage)) THEN
     CALL sub_memory(-size(trueage),'r','truage','trajec_seq')
     DEALLOCATE(trueage)
  END IF
  IF (ALLOCATED(gi)) THEN
     CALL sub_memory(-size(gi),'r','gi','trajec_seq')
     DEALLOCATE(gi)
  END IF
  IF (ALLOCATED(gj)) THEN
     CALL sub_memory(-size(gj),'r','gj','trajec_seq')
     DEALLOCATE(gj)
  END IF
  IF (ALLOCATED(gk)) THEN
     CALL sub_memory(-size(gk),'r','gk','trajec_seq')
     DEALLOCATE(gk)
  END IF
  IF (ALLOCATED(gl)) THEN
     CALL sub_memory(-size(gl),'r','gl','trajec_seq')
     DEALLOCATE(gl)
  END IF
  IF (ALLOCATED(hi)) THEN
     CALL sub_memory(-size(hi),'r','hi','trajec_seq')
     DEALLOCATE(hi)
  END IF
  IF (ALLOCATED(hj)) THEN
     CALL sub_memory(-size(hj),'r','hj','trajec_seq')
     DEALLOCATE(hj)
  END IF
  IF (ALLOCATED(hk)) THEN
     CALL sub_memory(-size(hk),'r','hk','trajec_seq')
     DEALLOCATE(hk)
  END IF
  IF (ALLOCATED(hl)) THEN
     CALL sub_memory(-size(hl),'r','hl','trajec_seq')
     DEALLOCATE(hl)
  END IF

  IF (ALLOCATED(tf)) THEN
     CALL sub_memory(-size(tf),'r','tf','trajec_seq')
     DEALLOCATE(tf)
  END IF
  IF (ALLOCATED(sf)) THEN
     CALL sub_memory(-size(sf),'r','sf','trajec_seq')
     DEALLOCATE(sf)
  END IF
  IF (ALLOCATED(rf)) THEN
     CALL sub_memory(-size(rf),'r','rf','trajec_seq')
     DEALLOCATE(rf)
  END IF

  !! NG: 15_09_2008
  IF (ALLOCATED(init_depth)) THEN
     CALL sub_memory(-size(init_depth),'i','init_depth','trajec_seq')
     DEALLOCATE(init_depth)
  END IF
  IF (ALLOCATED(final_depth)) THEN
     CALL sub_memory(-size(final_depth),'i','final_depth','trajec_seq')
     DEALLOCATE(final_depth)
  END IF

  IF (ALLOCATED(array_qual_traj)) THEN
     CALL sub_memory(-size(array_qual_traj),'r','array_qual_traj','trajec_seq')
     DEALLOCATE(array_qual_traj)
  END IF

  CALL sub_netcdf_dealloc_mem()
  CALL sub_input_data_dealloc_mem()
  CALL sub_coord_dealloc()
  CALL sub_scalef_dealloc()
  IF (key_roms)  CALL sub_h_roms_dealloc()
  CALL sub_transp_dealloc()
  IF (key_alltracers) THEN
     CALL sub_tracer_dealloc()
  ENDIF
  CALL sub_posin_dealloc()
  CALL sub_tmask_dealloc()


  IF (TRIM(mode) == 'quantitative') THEN ! [QUANT]
     CALL sub_txt_dealloc()
     CALL sub_flags_dealloc()
     CALL sub_stats_dealloc()
     CALL sub_stati_dealloc()
     CALL sub_quant_dealloc()
  ENDIF

  !=================
  ! MAIN PART - END
  !=================

END SUBROUTINE trajec
!!***
