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
MODULE mod_init_particules

  USE mod_precision
  USE mod_memory
  USE mod_namelist
  USE mod_lun
  USE mod_input_data
  USE mod_posin
  USE mod_seq

  IMPLICIT NONE

  INTEGER( kind = iprec) :: & !
       maxsect = 0        , & ! quantitative mode = max of sections
       maxsegm = 0        , & ! quantitative mode = max of segments
       ntraj   = 0        , & ! number of trajectories
       nsect   = 0        , & ! number of sections
       nfnt    = 0        , & !
       icrit1  = 0

  REAL(kind = rprec)     :: & !
       trtot   = 0._rprec  

CONTAINS

  !=========================================================================
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !========================================================================
  ! four options to initialize the positions:
  ! - qualitative: on ASCII file "initial_positions.txt"
  ! - qualitative: on BINARY file "initial.bin"
  ! - quantitative: auto positionning (from file "segments")
  ! - quantitative: from file "initial.bin"
  !=========================================================================
  !=========================================================================

  SUBROUTINE sub_init_particules_positions()

    INTEGER(kind = iprec) :: nb_dims
    INTEGER(kind = iprec), DIMENSION(1) :: dims

    !================!
    !- Declarations -!
    !================!
    INTEGER(kind=iprec) :: & !
         n               , & !
         ind             , & !
         nb_subset             , & !
         alloc_size      , & !
         ncid            , & !
         varid
    INTEGER(kind=iprec) :: is_loop

    INTEGER(kind=iprec) :: test_subset_order

    REAL(kind=rprec) ::     &  !
         trmin   =    1._rprec ! A particle is set if the transport is
    ! higher thant this value. 
    ! Reduce this value if your transport is
    ! very small.

    REAL(kind=rprec) :: maxtfl  ! used in backward mode august 2020 NG!


    INTEGER(kind=iprec), DIMENSION(:), ALLOCATABLE :: ind_subset


    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)'==============================='
    WRITE(lun_standard,*)'= INITIAL PARTICULE POSITIONS ='
    WRITE(lun_standard,*)'==============================='
    WRITE(lun_standard,*)'       (read or computed)'

    !=====================!
    !- Quantitative MODE -!
    !=====================!
    IF (TRIM(mode) == 'quantitative') THEN

      !-------------------------------------------------!
      !- Compute in sections.txt : maxsect and maxsegm -!
      !-------------------------------------------------!
      CALL sub_init_part_maxsect_maxsegm()

      !---------------------!
      !- Auto positioning - !
      !---------------------!
      IF (key_sequential) THEN

        id_comments = .TRUE. !- Write input data reading information -!

        CALL posini_seq(trmin,ntraj,nsect,trtot,nfnt)

      ELSE

        CALL posini(trmin, ntraj, nsect, trtot, nfnt)

      ENDIF

      !----------------------------------------
      ! initialization of transparent sections
      !----------------------------------------
      !!NG      IF (nfnt > 0) THEN
      !!NG        DO lu = 91, 90 + nfnt
      !!NG          WRITE(fname,'(a12,a15)')'transparent.',fntname(lu-90)
      !!NG          OPEN(lu,form='FORMATTED',file=fname)
      !!NG        END DO
      !!NG      ENDIF

    ENDIF

    !=====================================================!
    !- indices read on ASCII file (QUALITATIVE analysis) -!
    !=====================================================!
    IF ((TRIM(mode) == 'qualitative').AND.(TRIM(bin) == 'nobin')) THEN

      OPEN(unit=lun_dummy, file='initial_positions.txt', action="read")

      DO n = 1, nmax + 1

        READ(lun_dummy,*,END=7054,ERR=7054)tfi(n),tfj(n),tfk(n),tfl(n),ttr(n)

      END DO

      WRITE(lun_error,*)
      WRITE(lun_error,*)'# of positions in initial.positions.txt exceeds nmax'
      WRITE(lun_error,*)'--- STOP ARIANE ---'
      STOP

7054  ntraj = n - 1

      IF (key_ascii_outputs) THEN
        WRITE(lun_output,*)
        WRITE(lun_output,*)'======================================================'
        WRITE(lun_output,*)'ASCII positions are read on file initial_positions.txt'
        WRITE(lun_output,*)'======================================================'
        WRITE(lun_output,*)
        WRITE(lun_output,*)'# of particles that could be read: ',ntraj
      ENDIF

      IF ( ntraj == 0) THEN
        WRITE(lun_error,*)
        WRITE(lun_error,*)')-: The number of particule is 0 :-( '
        WRITE(lun_error,*)'Please verify that you have correctly set the file'
        WRITE(lun_error,*)'initial.positions.txt (param0 is obsoleted).'
        WRITE(lun_error,*)'--- STOP ARIANE ---'
      ENDIF

      IF (key_sequential) THEN

        IF (key_alltracers) THEN
          IF (key_roms) THEN
            alloc_size = 7
          ELSEIF (key_symphonie) THEN
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
          IF (key_roms) THEN
            alloc_size = 3
          ELSEIF (key_symphonie) THEN
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

!!$          If (key_alltracers) Then
!!$             IF ((key_roms).OR.(key_symphonie))THEN
!!$                alloc_size = 7
!!$             ELSE
!!$                alloc_size = 6
!!$             ENDIF
!!$          Else
!!$             alloc_size = 3
!!$          Endif

        CALL sub_seq_alloc(alloc_size)

      ENDIF

    ENDIF

    !===========================================!
    !- initial positions in file "initial.bin" -!
    !===========================================!
    IF (TRIM(bin) /= 'nobin') THEN
      !-----------------------------
      ! indices read on BINARY file
      !-----------------------------
      !!NG      WRITE(lun_output,*)
      !!NG      WRITE(lun_output,*)'==============================================='
      !!NG      WRITE(lun_output,*)'BINARY positions are read from file initial.bin'
      !!NG      WRITE(lun_output,*)'==============================================='
      !!NG
      !!NG      OPEN(UNIT=lun_initial,file='initial.bin', form='UNFORMATTED', ACTION='READ')
      !!NG
      !!NG      DO n = 1, nmax + 1
      !!NG
      !!NG        READ(UNIT=lun_initial,END=7053,ERR=7053)tfi(n),tfj(n),tfk(n),tfl(n),ttr(n)
      !!NG
      !!NG        trtot = trtot + ttr(n)
      !!NG
      !!NG      END DO
      !!NG
      !!NG      WRITE(lun_error,*)
      !!NG      WRITE(lun_error,*)'# of available positions in initial.bin exceeds nmax'
      !!NG      WRITE(lun_error,*)'--- STOP ARIANE ---'
      !!NG      STOP
      !!NG
      !!NG7053  ntraj = n - 1 

      IF (key_ascii_outputs) THEN
        WRITE(lun_output,*)
        WRITE(lun_output,*)'=================================================='
        WRITE(lun_output,*)'= Positions are read from ariane_initial.nc file ='
        WRITE(lun_output,*)'=================================================='
      ENDIF

      !------------------------!
      !- Open the netcdf file -!
      !------------------------!
      CALL sub_open_netcdf_file('.', 'ariane_initial.nc', ncid)
      WRITE(lun_standard, *)'  --- ariane_initial.nc is opened ---'

      CALL sub_read_netcdf_varid_ndims(     &
           ncid   = ncid                  , &
           nc_var = TRIM(init_final)//'_x', &
           varid  = varid                 , &
           ndims  = nb_dims                 )

      CALL sub_read_netcdf_var_dims( &
           ncid = ncid             , &
           varid = varid           , &
           ndims = nb_dims         , &
           dims  = dims(:)           )

      ntraj = dims(1)

      IF (key_ascii_outputs) THEN
        WRITE(lun_output,*)
        WRITE(lun_output,*)'# of particles that could be read: ',ntraj
      ENDIF

      IF (ntraj > nmax) THEN
        WRITE(lun_error,*)''
        WRITE(lun_error,*)'Error: mod_init_particules.f90: sub_init_particules_positions'
        WRITE(lun_error,*)'Error: nmax < ntraj !!!'
        WRITE(lun_error,*)'Error: Please, in the namelist file, increase the value of nmax'
        WRITE(lun_error,*)'Error: to verify nmax >= ntraj...'
        STOP
      ENDIF


      !-------!
      !- TFI -! (rien a voir avec la chaine TF1 ;-)
      !-------!
      CALL sub_read_netcdf_var1d( ncid, varid, tfi(1:ntraj))
      WRITE(lun_standard, *)'  - X positions are read from ', &
           &TRIM(init_final)//'_x', ' -'

      !-------!
      !- TFJ -!
      !-------!
      CALL sub_read_netcdf_varid_ndims(ncid, TRIM(init_final)//'_y', varid)

      CALL sub_read_netcdf_var1d( ncid, varid, tfj(1:ntraj))
      WRITE(lun_standard, *)'  - Y positions are read from ', &
           &TRIM(init_final)//'_y', ' -'

      !-------!
      !- TFK -!
      !-------!
      CALL sub_read_netcdf_varid_ndims(ncid, TRIM(init_final)//'_z', varid)

      CALL sub_read_netcdf_var1d( ncid, varid, tfk(1:ntraj))
      WRITE(lun_standard, *)'  - Z positions are read from ', &
           &TRIM(init_final)//'_z', ' -'

      !-------!
      !- TFL -!
      !-------!
      CALL sub_read_netcdf_varid_ndims(ncid, TRIM(init_final)//'_t', varid)

      CALL sub_read_netcdf_var1d( ncid, varid, tfl(1:ntraj))
      WRITE(lun_standard, *)'  - Time positions read from ', &
           &TRIM(init_final)//'_t', ' -'

      !-------!
      !- AGE -!
      !-------!
      IF (key_read_age) THEN
        CALL sub_read_netcdf_varid_ndims(ncid, TRIM(init_final)//'_age', varid)

        CALL sub_read_netcdf_var1d( ncid, varid, tage(1:ntraj))

        trtot = SUM(tage(1:ntraj)) / REAL(ntraj, kind = rprec)

        WRITE(lun_standard, *)'  - Age are read from ', &
             &TRIM(init_final)//'_age', ' - (mean =',  trtot, ')'
      ENDIF

      !-------!
      !- TTR -! (transport)
      !-------!
      CALL sub_read_netcdf_varid_ndims(ncid, TRIM(init_final)//'_transp', varid)

      CALL sub_read_netcdf_var1d( ncid, varid, ttr(1:ntraj))

      trtot = SUM(ttr(1:ntraj))

      WRITE(lun_standard, *)'  - Transports are read from ', &
           &TRIM(init_final)//'_transp', ' - (total = ', trtot, ')'

      IF (key_alltracers.AND.(TRIM(mode) == 'quantitative')) THEN

        !---------------!
        !- TEMPERATURE -!
        !---------------!
        CALL sub_read_netcdf_varid_ndims(ncid, TRIM(init_final)//'_temp', varid)

        CALL sub_read_netcdf_var1d( ncid, varid, init_temp(1:ntraj))
        WRITE(lun_standard, *)'  - Temperature values are read from ', &
             &TRIM(init_final)//'_temp', ' -'

        !------------!
        !- SALYNITY -!
        !------------!
        CALL sub_read_netcdf_varid_ndims(ncid, TRIM(init_final)//'_salt', varid)

        CALL sub_read_netcdf_var1d( ncid, varid, init_salt(1:ntraj))
        WRITE(lun_standard, *)'  - Salinity values are read from ', &
             &TRIM(init_final)//'_salt', ' -'

        !-----------!
        !- DENSITY -!
        !-----------!
        CALL sub_read_netcdf_varid_ndims(ncid, TRIM(init_final)//'_dens', varid)

        CALL sub_read_netcdf_var1d( ncid, varid, init_dens(1:ntraj))
        WRITE(lun_standard, *)'  - Density values are read from ', &
             &TRIM(init_final)//'_dens', ' -'

      ENDIF

      !-------------------------!
      !- Close the netcdf file -!
      !-------------------------!
      CALL sub_close_netcdf_file(ncid)
      WRITE(lun_standard, *)'  --- ariane_initial.nc is closed ---'

    ENDIF ! bin /= 'nobin'

    !=========================!
    !- Read a subset of data -!
    !=========================!
    IF (TRIM(bin) == 'subbin') THEN

      ALLOCATE(ind_subset(nmax))
      CALL sub_memory(SIZE(ind_subset),'i','ind_subset','sub_init_particules_positions')

      ind_subset(:)=0

      OPEN(lun_subset, file='subset.txt')

      WRITE(lun_output,*)'A subset of indices is reading on file subset.txt.'

      DO nb_subset = 1, ntraj+1
        READ(lun_subset,*,END=1010,ERR=1010)ind_subset(nb_subset)
      END DO

      WRITE(lun_error,*)
      WRITE(lun_error,*)'# of indices in subset exceeds ntraj'
      WRITE(lun_error,*)'--- STOP ARIANE ---'
      STOP

1010  CONTINUE

      CLOSE (lun_subset)

      nb_subset = nb_subset - 1

      IF (MAXVAL(ind_subset(:)) > ntraj) THEN
        WRITE(lun_error,*)'Error: mod_init_particules.f90'
        WRITE(lun_error,*)' An indice in the subset.txt file is'
        WRITE(lun_error,*)' higher than the number of particles !!!'
        WRITE(lun_error,*)'--- STOP ARIANE ---'
        STOP
      ENDIF

      trtot     = 0._rprec

      ! Because the indices of the subset particules increase
      ! we can not use an intermediate array.
      ! This means that ind-subset(ind) is always >= at ind.
      ! A little bit dangerous...
      test_subset_order = 0
      DO ind = 1, nb_subset
        if (ind_subset(ind) < test_subset_order) then
          WRITE(lun_error,*)
          WRITE(lun_error,*)'# In subset.txt file, particles are not sorted !!!'
          WRITE(lun_error,*)'# Particle number have to increase and never decrease !!!'
          WRITE(lun_error,*)'--- STOP ARIANE ---'
          STOP
        endif
        test_subset_order = ind_subset(ind)
        tfi(ind) = tfi(ind_subset(ind))
        tfj(ind) = tfj(ind_subset(ind))
        tfk(ind) = tfk(ind_subset(ind))
        tfl(ind) = tfl(ind_subset(ind))
        ttr(ind) = ttr(ind_subset(ind))
        tage(ind) = tage(ind_subset(ind))
        trtot    = trtot + ttr(ind)
        IF (key_alltracers.AND.(TRIM(mode) == 'quantitative')) THEN
          init_temp(ind) = init_temp(ind_subset(ind))
          init_salt(ind) = init_salt(ind_subset(ind))
          init_dens(ind) = init_dens(ind_subset(ind))
        ENDIF
      END DO


      ntraj = nb_subset

      WRITE(lun_output,*)
      WRITE(lun_output,*)'# of particles that were kept as a subset: ', ntraj

      CALL sub_memory(-SIZE(ind_subset),'i','ind_subset','sub_init_particules_positions')
      DEALLOCATE(ind_subset)


    ENDIF ! bin == subbin

    !- TO BE SURE THAT TIME INDICES ARE < lmt
    IF(TRIM(forback) == 'forward')  THEN !NG august 2020

       IF (MAXVAL(tfl) > REAL(lmt,kind = rprec)) THEN
          WRITE(lun_standard,*)'  - One or more time indices are greater than the number of records !'
          !WRITE(lun_standard,*) MAXVAL(tfl),' > ', REAL(lmt,kind = rprec) 
          IF ((key_sequential).AND.(maxcycles > 1)) THEN
             WRITE(lun_standard,*)'  - In sequential mode and maxcycle > 1'
             WRITE(lun_standard,*)'  - Calculation of the time indices modulo number of records!'
             WHERE (tfl > REAL(lmt,kind = rprec)) 
                tfl = mod(tfl,REAL(lmt,kind = rprec))
             END WHERE
          ELSEIF (.NOT.(key_sequential)) THEN
             WRITE(lun_standard,*)'  - InNo sequential mode'
             WRITE(lun_standard,*)'  - Calculation of the time indices modulo number of records!'
             WHERE (tfl > REAL(lmt,kind = rprec)) 
                tfl = mod(tfl,REAL(lmt,kind = rprec))
             END WHERE
          ELSE
             WRITE(lun_error,*)'Error: One or more time indices are greater than the number of records'
             WRITE(lun_error,*)'Error:',MAXVAL(tfl),' > ', REAL(lmt,kind = rprec) 
             WRITE(lun_error,*)'Error: key_sequential is .FALSE. or maxcycles <= 1 !!'
             WRITE(lun_error,*)'Error: ARIANE STOP...'
             STOP
          ENDIF
       ENDIF

    ELSE ! backward !NG august 2020

       IF (MAXVAL(tfl) > REAL(lmt,kind = rprec)+rHalf) THEN
          WRITE(lun_standard,*)'  - One or more time indices are greater than the number of records !'
          !WRITE(lun_standard,*) MAXVAL(tfl),' > ', REAL(lmt,kind = rprec) 
          IF ((key_sequential).AND.(maxcycles > 1)) THEN
             WRITE(lun_standard,*)'  - In sequential mode and maxcycle > 1'
             WRITE(lun_standard,*)'  - Calculation of the time indices modulo number of records!'
             WHERE (tfl > REAL(lmt,kind = rprec)) 
                tfl = mod(tfl,REAL(lmt,kind = rprec))
             END WHERE
          ELSEIF (.NOT.(key_sequential)) THEN
             WRITE(lun_standard,*)'  - In No sequential mode'
             WRITE(lun_standard,*)'  - Calculation of the time indices modulo number of records!'
             WHERE (tfl > REAL(lmt,kind = rprec)) 
                tfl = mod(tfl,REAL(lmt,kind = rprec))
             END WHERE
          ELSE
             WRITE(lun_error,*)'Error: One or more time indices are greater than the number of records'
             WRITE(lun_error,*)'Error:',MAXVAL(tfl),' > ', REAL(lmt,kind = rprec) 
             WRITE(lun_error,*)'Error: key_sequential is .FALSE. or maxcycles <= 1 !!'
             WRITE(lun_error,*)'Error: ARIANE STOP...'
             STOP
          ENDIF
       ELSEIF(MAXVAL(tfl) == REAL(lmt,kind = rprec)+rHalf) THEN
          WHERE (tfl == REAL(lmt,kind = rprec)+rHalf) 
             tfl = tfl - 0.000000000000001
          END WHERE

       ENDIF

    ENDIF

    !========================================!
    !- Write initial positions in ASCI file -!
    !========================================!
    IF ((TRIM(mode) == 'quantitative').AND.key_ascii_outputs ) THEN

      WRITE(lun_standard,*)'  - Writing initial positions...'

      DO is_loop = 1, ntraj
        WRITE(lun_init_pos,7050) &
             tfi(is_loop),tfj(is_loop),tfk(is_loop),tfl(is_loop), &
             ttr(is_loop)
      ENDDO

      CLOSE(lun_init_pos)

    ENDIF

7050 FORMAT(5(1x,f0.3))

    !---------------------------------------------!
    !- Compute lmin and lmax in qualitative mode -!
    !---------------------------------------------!
    IF (TRIM(mode) == 'qualitative') THEN

      IF(TRIM(forback) == 'forward')  THEN

        WRITE(lun_standard,*) &
             '  - Computing lmin in qualitative mode...'
        lmin = NINT(MINVAL(tfl, mask = tfl >  rZero  ))
        WRITE(lun_standard,*) '    - lmin =', lmin

        IF (lmin == 0) THEN
          WRITE(lun_error,*) &
               ' lmin = 0, there is a problem in the initial_positions.txt file'
          STOP
        ENDIF

      ELSEIF(TRIM(forback) == 'backward') THEN

        WRITE(lun_standard,*) &
             '  - Computing lmax in qualitative mode...'

          maxtfl=MAXVAL(tfl, mask = tfl > rZero)

          if (maxtfl==REAL(lmt,kind = rprec)+rHalf) then
             maxtfl = maxtfl - 0.000000000000001
          endif

          lmax = NINT(maxtfl)
          WRITE(lun_standard,*) '    - lmax =', lmax

       ENDIF

    ENDIF

    !======================================!
    !- Allocate arrays in sequential mode -!
    !======================================!
    IF ((TRIM(bin) /= 'nobin')) THEN

      IF (key_sequential) THEN

        IF (key_alltracers) THEN
          IF (key_roms) THEN
            alloc_size = 7
          ELSEIF (key_symphonie) THEN
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
          IF (key_roms) THEN
            alloc_size = 3
          ELSEIF (key_symphonie) THEN
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
!!$          If (key_alltracers) Then
!!$             IF ((key_roms).OR.(key_symphonie))THEN
!!$                alloc_size = 7
!!$             ELSE
!!$                alloc_size = 6
!!$             ENDIF
!!$          Else
!!$             alloc_size = 3
!!$          Endif

        CALL sub_seq_alloc(alloc_size)
      ENDIF

    ENDIF

    !========================================================================
    ! I/O - I/O - I/O - I/O - I/O - I/O - I/O - I/O - I/O - I/O - I/O - 
    ! [QUANT]
    !========================================================================
    IF ((TRIM(mode) == 'quantitative').and.key_ascii_outputs) THEN
      WRITE(lun_output,*)' '
      WRITE(lun_output,*)'QUANTITATIVE EXPERIMENT'
      WRITE(lun_output,*)'======================='
      WRITE(lun_output,*)'    minimum gridcell transport: ', trmin
      WRITE(lun_output,*)'  maximum individual transport: ', max_transport
      WRITE(lun_output,*)'      exact number of segments: ', maxsegm
      WRITE(lun_output,*)'     exact number of particles: ', ntraj
      WRITE(lun_output,*)'      exact number of sections: ', nsect
      WRITE(lun_output,*)'    total documented transport: ', trtot
      WRITE(lun_output,*)'number of transparent sections: ', nfnt
      WRITE(lun_output,*)' '
    ENDIF

  END SUBROUTINE sub_init_particules_positions

  !=========================================================================
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !=========================================================================

  SUBROUTINE sub_init_part_maxsect_maxsegm()

    INTEGER(kind=iprec) :: idum

    !-----------------------------------------------------------------!
    !- Open the file sections.txt and compute the number of segments -!
    !- and the number of sections.                                   -!
    !- (remember that a section could be composed by more than       -!
    !-  one segment)                                                 -!
    !-----------------------------------------------------------------!
    OPEN(unit=lun_dummy, file='sections.txt', action="read")
    DO 
      READ(lun_dummy,*, END=7777, ERR=7777) idum
      maxsegm = maxsegm + 1
      IF (idum > maxsect) THEN
        maxsect = idum
      ENDIF
    ENDDO

7777 CONTINUE

    IF (key_ascii_outputs) THEN
      WRITE(lun_output,*)'Number of segments in sections.txt file = ', maxsegm
    ENDIF

    CLOSE(lun_dummy)

  END SUBROUTINE sub_init_part_maxsect_maxsegm

END MODULE mod_init_particules