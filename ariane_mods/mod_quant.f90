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
!!****h* ariane/mod_quant
!! NAME
!!   mod_quant (mod_quant.f90 - Fortran90 module)
!!
!! USAGE
!!   Include 'USE mod_quant' in the header of your Fortran 90 source 
!!   code.
!!   Then you'll have access to the subroutine:
!!      - sub_quant_initside_u
!!      - sub_quant_initside_v
!!      - sub_quant_initside_w
!!      - sub_quant_reachside_u
!!      - sub_quant_reachside_v
!!      - sub_quant_reachside_w
!!      - sub_quant_store_transp
!!      - sub_quant_statistics
!!      - sub_quant_alloc
!!      - sub_quant_init
!!      - sub_quant_initorig
!!      - sub_quant_inittmp
!!      - sub_quant_dealloc
!!
!! FUNCTION
!!   All quantitative specific parts.
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
!!
!! USES
!!   * USE mod_precision
!!   * USE mod_namelist
!!   * USE mod_input_grid
!!   * USE mod_input_data
!!   * USE mod_flags
!!   * USE mod_stati
!!   * USE mod_stats
!!   * USE mod_lun
!!   * USE mod_txt
!!   * USE mod_fx
!!   * USE mod_fy
!!   * USE mod_fz
!!   * USE mod_criter0
!!   * USE mod_criter1
!!   * USE mod_criter2
!!   * USE mod_sigma
!!   * USE mod_zinter
!!
!! USED BY
!!   * trajec
!!
!! SOURCE
!!=========================================================================
MODULE mod_quant

  !------------------!
  ! USE ASSOCIAITION !
  !------------------!
  USE mod_precision
  USE mod_memory
  USE mod_namelist
  USE mod_input_grid
  USE mod_input_data
  USE mod_init_particules
  USE mod_flags
  USE mod_stati
  USE mod_stats
  USE mod_lun
  USE mod_txt
  USE mod_fx
  USE mod_fy
  USE mod_fz
  USE mod_criter0
  USE mod_criter1
  USE mod_criter2
  USE mod_sigma
  USE mod_zinter

  !-------------!
  ! DECLARATION !
  !-------------!
  IMPLICIT NONE

  REAL(kind=rprec), DIMENSION(:,:,:), ALLOCATABLE :: &
       uxy, & !
       vxy, & !
       uxz, & !
       wxz, & !
       vyz, & !
       wyz    !

  REAL(kind=rprec), DIMENSION(:,:)  , ALLOCATABLE :: &
       uxytmp, & !
       vxytmp, & !
       uxztmp, & !
       wxztmp, & !
       vyztmp, & !
       wyztmp    !

  REAL(kind=rprec), DIMENSION(:,:,:), ALLOCATABLE :: &
       uh  , & !
       vh  , & !
       zuh , & !
       zvh , & !
       z2uh, & !
       z2vh    !

  REAL(kind=rprec), DIMENSION(:,:)  , ALLOCATABLE :: &
       uht  , & !
       vht  , & !
       zuht , & !
       zvht , & !
       z2uht, & !
       z2vht 

  REAL(kind=rprec), DIMENSION(:,:,:), ALLOCATABLE :: & 
       ruh , & !
       rvh , & !
       r2uh, & !
       r2vh, & !
       suh , & !
       svh , & !
       s2uh, & !
       s2vh, & !
       tuh , & !
       tvh , & !
       t2uh, & !
       t2vh    !

  REAL(kind=rprec), DIMENSION(:,:)  , ALLOCATABLE :: &
       ruht , & !
       rvht , & !
       r2uht, & !
       r2vht, & !
       suht , & !
       svht , & !
       s2uht, & !
       s2vht, & !
       tuht , & !
       tvht , & !
       t2uht, & !
       t2vht    !

CONTAINS
  !!***
  !=========================================================================
  !!****f* mod_quant/sub_quant_initside_u()
  !! NAME
  !!   sub_quant_initside_u()
  !!
  !! FUNCTION
  !!  the transport on the initial section is stored in the relevant
  !!  transport arrays.
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
  !!   * INPUTS:
  !!       - ii    : i indice
  !!       - jj    : j indice
  !!       - kk    : k indice
  !!       - ll    : time indice (l)
  !!       - hi    :
  !!       - hj    :
  !!       - hk    :
  !!       - hl    :
  !!       - trans : transport
  !!       - zsigma:
  !!
  !! TODO
  !!   
  !!
  !! USED BY
  !!   * trajec
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_quant_initside_u(ii, jj, kk, ll, prof, trans, temp, salt, dens)

    !-------------!
    ! Declaration !
    !-------------!
    !- arguments -!
    INTEGER(kind=iprec), INTENT(in) :: ii, jj, kk, ll
    REAL(kind=rprec),    INTENT(in) :: prof  ! depth
    REAL(kind=rprec),    INTENT(in) :: trans ! transport

    REAL(kind=rprec), OPTIONAL, INTENT(in) :: temp    ! temperature
    REAL(kind=rprec), OPTIONAL, INTENT(in) :: salt    ! salinity
    REAL(kind=rprec), OPTIONAL, INTENT(in) :: dens   ! density

    !-------------!
    ! Code begins !
    !-------------!

    IF (key_eco) THEN

      uxy(ii,jj,1)  =  uxy(ii,jj,1) + trans
      uxz(ii,kk,1)  =  uxz(ii,kk,1) + trans
      uh(ii,jj,1)   =   uh(ii,jj,1) + trans
      zuh(ii,jj,1)  =  zuh(ii,jj,1) - prof        * trans
      z2uh(ii,jj,1) = z2uh(ii,jj,1) + prof * prof * trans

      IF (key_alltracers) THEN

        IF (PRESENT(temp).AND.PRESENT(salt).AND.PRESENT(dens)) THEN
          ruh(ii,jj,1) =  ruh(ii,jj,1)  + dens        * trans
          r2uh(ii,jj,1) = r2uh(ii,jj,1) + dens * dens * trans
          suh(ii,jj,1) =  suh(ii,jj,1)  + salt        * trans
          s2uh(ii,jj,1) = s2uh(ii,jj,1) + salt * salt * trans
          tuh(ii,jj,1) =  tuh(ii,jj,1)  + temp        * trans
          t2uh(ii,jj,1) = t2uh(ii,jj,1) + temp * temp * trans
        ELSE
          STOP
        ENDIF

!!NG: seq-roms          kr=nint((float(nrclas-1)*ri(n)+(rclasmax-rclasmin*nrclas))/drclas )
!!NG: seq-roms          if (kr <= 1) kr=1
!!NG: seq-roms          if (kr >= nrclas) kr=nrclas
!!NG: seq-roms          uxr(nint(fi(n)),kr)=uxr(nint(fi(n)),kr)+ttr(n)

      ENDIF

    ELSE

      uxytmp(ii,jj) = trans
      uxztmp(ii,kk) = trans
      uht(ii,jj)    = trans
      zuht(ii,jj)   =-prof        * trans
      z2uht(ii,jj)  = prof * prof * trans

      IF (key_alltracers) THEN

        IF (PRESENT(temp).AND.PRESENT(salt).AND.PRESENT(dens)) THEN
          ruht(ii,jj) = dens       * trans
          r2uht(ii,jj)= dens * dens * trans
          suht(ii,jj) = salt        * trans
          s2uht(ii,jj)= salt  * salt  * trans
          tuht(ii,jj) = temp        * trans
          t2uht(ii,jj)= temp  * temp  * trans
        ELSE
          STOP
        ENDIF

      ENDIF

    ENDIF

  END SUBROUTINE sub_quant_initside_u
  !!***
  !=========================================================================
  !!****f* mod_quant/sub_quant_initside_v()
  !! NAME
  !!   sub_quant_initside_v()
  !!
  !! FUNCTION
  !!  the transport on the initial section is stored in the relevant
  !!  transport arrays.
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
  !!   * INPUTS:
  !!       - ii    : i indice
  !!       - jj    : j indice
  !!       - kk    : k indice
  !!       - ll    : time indice (l)
  !!       - hi    :
  !!       - hj    :
  !!       - hk    :
  !!       - hl    :
  !!       - trans : transport
  !!       - zsigma:
  !!
  !! TODO
  !!   
  !!
  !! USED BY
  !!   * trajec
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_quant_initside_v(ii, jj, kk, ll, prof, trans, temp, salt, dens)

    !-------------!
    ! Declaration !
    !-------------!
    !- arguments -!
    INTEGER(kind=iprec), INTENT(in)  :: ii, jj, kk, ll
    REAL(kind=rprec),    INTENT(in) :: prof  ! depth
    REAL(kind=rprec),    INTENT(in) :: trans ! transport

    REAL(kind=rprec), OPTIONAL, INTENT(in) :: temp    ! temperature
    REAL(kind=rprec), OPTIONAL, INTENT(in) :: salt    ! salinity
    REAL(kind=rprec), OPTIONAL, INTENT(in) :: dens   ! density

    !-------------!
    ! Code begins !
    !-------------!

    IF (key_eco) THEN

      vxy(ii,jj,1)=vxy(ii,jj,1)+trans
      vyz(jj,kk,1)=vyz(jj,kk,1)+trans
      vh(ii,jj,1)=vh(ii,jj,1)+trans

      zvh(ii,jj,1)=zvh(ii,jj,1)   - prof        * trans
      z2vh(ii,jj,1)=z2vh(ii,jj,1) + prof * prof * trans
      IF (key_alltracers) THEN


        IF (PRESENT(temp).AND.PRESENT(salt).AND.PRESENT(dens)) THEN
          rvh(ii,jj,1)=rvh(ii,jj,1)   + dens       * trans
          r2vh(ii,jj,1)=r2vh(ii,jj,1) + dens * dens * trans
          svh(ii,jj,1)=svh(ii,jj,1)   + salt        * trans
          s2vh(ii,jj,1)=s2vh(ii,jj,1) + salt  * salt  * trans
          tvh(ii,jj,1)=tvh(ii,jj,1)   + temp        * trans
          t2vh(ii,jj,1)=t2vh(ii,jj,1) + temp  * temp  * trans
        ELSE
          STOP
        ENDIF
!!NG: seq-roms           kr=nint((float(nrclas-1)*ri(n)+(rclasmax-rclasmin*nrclas))/drclas )
!!NG: seq-roms           if (kr <= 1) kr=1
!!NG: seq-roms           if (kr >= nrclas) kr=nrclas
!!NG: seq-roms           uxr(nint(tempi(n)),kr)=uxr(nint(tempi(n)),kr)+ttr(n)
      ENDIF

    ELSE

      vxytmp(ii,jj)=trans
      vyztmp(jj,kk)=trans
      vht(ii,jj)=trans

      zvht(ii,jj)  = -prof        * trans
      z2vht(ii,jj) =  prof * prof * trans
      IF (key_alltracers) THEN

        IF (PRESENT(temp).AND.PRESENT(salt).AND.PRESENT(dens)) THEN
          rvht(ii,jj)  = dens       * trans
          r2vht(ii,jj) = dens * dens * trans
          svht(ii,jj)  = salt        * trans
          s2vht(ii,jj) = salt  * salt  * trans
          tvht(ii,jj)  = temp        * trans
          t2vht(ii,jj) = temp  * temp  * trans
        ELSE
          STOP
        ENDIF
      ENDIF

    ENDIF

  END SUBROUTINE sub_quant_initside_v
  !!***
  !=========================================================================
  !!****f* mod_quant/sub_quant_initside_w()
  !! NAME
  !!   sub_quant_initside_w()
  !!
  !! FUNCTION
  !!  the transport on the initial section is stored in the relevant
  !!  transport arrays.
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
  !!   * INPUTS:
  !!       - ii    : i indice
  !!       - jj    : j indice
  !!       - kk    : k indice
  !!       - ll    : time indice (l)
  !!       - hi    :
  !!       - hj    :
  !!       - hk    :
  !!       - hl    :
  !!       - trans : transport
  !!       - zsigma:
  !!
  !! TODO
  !!   
  !!
  !! USED BY
  !!   * trajec
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_quant_initside_w(ii, jj, kk, ll, trans)

    !-------------!
    ! Declaration !
    !-------------!
    !- arguments -!
    INTEGER(kind=iprec), INTENT(in)  :: ii, jj, kk, ll
    REAL(kind=rprec),    INTENT(in)  :: trans

    !-------------!
    ! Code begins !
    !-------------!

    IF (key_eco) THEN
      wxz(ii,kk,1)=wxz(ii,kk,1)-trans
      wyz(jj,kk,1)=wyz(jj,kk,1)-trans
    ELSE
      wxztmp(ii,kk)=-trans
      wyztmp(jj,kk)=-trans
    ENDIF

  END SUBROUTINE sub_quant_initside_w
  !!***
  !=========================================================================
  !!****f* mod_quant/sub_quant_reachside_u()
  !! NAME
  !!   sub_quant_reachside_u()
  !!
  !! FUNCTION
  !!  the transport on the initial section is stored in the relevant
  !!  transport arrays.
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
  !!   * INPUTS:
  !!       - ii    : i indice
  !!       - jj    : j indice
  !!       - kk    : k indice
  !!       - ll    : time indice (l)
  !!       - hi    :
  !!       - hj    :
  !!       - hk    :
  !!       - hl    :
  !!       - trans : transport
  !!       - zsigma:
  !!
  !! TODO
  !!   
  !!
  !! USED BY
  !!   * trajec
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_quant_reachside_u(ii, jj, kk, ll, hi, hj, hk, hl, trans)

    !-------------!
    ! Declaration !
    !-------------!
    !- arguments -!
    INTEGER(kind=iprec), INTENT(in)  :: ii, jj, kk, ll
    REAL(kind=rprec),    INTENT(in)  :: hi, hj, hk, hl, trans

    !- local variables -!
    INTEGER(kind=iprec) :: iis, jjs, kks
    REAL(kind=rprec)    :: prof !
    REAL(kind=rprec)    :: temp   !
    REAL(kind=rprec)    :: salt   !
    REAL(kind=rprec)    :: dens  !

    !-------------!
    ! Code begins !
    !-------------!
    IF (key_eco) THEN

      uh(ii,jj,1) = uh(ii,jj,1) + trans

      IF (key_roms.OR.key_mars.OR.key_symphonie) THEN
        prof=fz(gi=hi, gj=hj, gk=-ABS(hk), il=ll)
      ELSE
        prof=fz(gi=hi, gj=hj, gk=-ABS(hk))
      ENDIF

      zuh(ii,jj,1)  = zuh(ii,jj,1)  - prof        * trans
      z2uh(ii,jj,1) = z2uh(ii,jj,1) + prof * prof * trans

      IF (key_alltracers) THEN
        temp  = zinter(tt,hi,hj,hk,hl)
        salt  = zinter(ss,hi,hj,hk,hl)
        IF (key_approximatesigma) THEN
          dens = zinter(rr,hi,hj,hk,hl)
        ELSE
          dens = sigma(zsigma,salt,temp)
        ENDIF
        ruh(ii,jj,1)  = ruh(ii,jj,1)  + dens       * trans
        r2uh(ii,jj,1) = r2uh(ii,jj,1) + dens * dens * trans
        suh(ii,jj,1)  = suh(ii,jj,1)  + salt        * trans
        s2uh(ii,jj,1) = s2uh(ii,jj,1) + salt  * salt  * trans
        tuh(ii,jj,1)  = tuh(ii,jj,1)  + temp        * trans
        t2uh(ii,jj,1) = t2uh(ii,jj,1) + temp  * temp  * trans
      ENDIF

      CALL sub_reducmem_shift_or_not_ind(ii,jj,kk,iis,jjs,kks)

      IF (uu(iis,jjs,kks,ll) > 0._rprec) THEN
        uxy(ii,jj,1) = uxy(ii,jj,1) + trans
        uxz(ii,kk,1) = uxz(ii,kk,1) + trans
      ELSE
        uxy(ii,jj,1) = uxy(ii,jj,1) - trans
        uxz(ii,kk,1) = uxz(ii,kk,1) - trans
      ENDIF

    ELSE

      uht(ii,jj) = uht(ii,jj) + trans

      IF (key_roms.OR.key_mars.OR.key_symphonie) THEN
        prof=fz(gi=hi, gj=hj, gk=-ABS(hk), il=ll)
      ELSE
        prof=fz(gi=hi, gj=hj, gk=-ABS(hk))
      ENDIF

      zuht(ii,jj)  = zuht(ii,jj)  - prof        * trans
      z2uht(ii,jj) = z2uht(ii,jj) + prof * prof * trans

      IF (key_alltracers) THEN
        temp  = zinter(tt,hi,hj,hk,hl)
        salt  = zinter(ss,hi,hj,hk,hl)
        IF (key_approximatesigma) THEN
          dens = zinter(rr,hi,hj,hk,hl)
        ELSE
          dens = sigma(zsigma,salt,temp)
        ENDIF
        ruht(ii,jj)  = ruht(ii,jj)  + dens       * trans
        r2uht(ii,jj) = r2uht(ii,jj) + dens * dens * trans
        suht(ii,jj)  = suht(ii,jj)  + salt        * trans
        s2uht(ii,jj) = s2uht(ii,jj) + salt  * salt  * trans
        tuht(ii,jj)  = tuht(ii,jj)  + temp        * trans
        t2uht(ii,jj) = t2uht(ii,jj) + temp  * temp  * trans
      ENDIF

      CALL sub_reducmem_shift_or_not_ind(ii,jj,kk,iis,jjs,kks)

      IF (uu(iis,jjs,kks,ll) > 0._rprec) THEN
        uxytmp(ii,jj) = uxytmp(ii,jj) + trans
        uxztmp(ii,kk) = uxztmp(ii,kk) + trans
      ELSE
        uxytmp(ii,jj) = uxytmp(ii,jj) - trans
        uxztmp(ii,kk) = uxztmp(ii,kk) - trans
      ENDIF
    ENDIF

  END SUBROUTINE sub_quant_reachside_u

  !!***
  !=========================================================================
  !!****f* mod_quant/sub_quant_reachside_v()
  !! NAME
  !!   sub_quant_reachside_v()
  !!
  !! FUNCTION
  !!   the transport on the initial section is stored in the relevant
  !!   transport arrays.
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
  !!   * INPUTS:
  !!       - ii    : i indice
  !!       - jj    : j indice
  !!       - kk    : k indice
  !!       - ll    : time indice (l)
  !!       - hi    :
  !!       - hj    :
  !!       - hk    :
  !!       - hl    :
  !!       - trans : transport
  !!       - zsigma:
  !!
  !! TODO
  !!   
  !!
  !! USED BY
  !!   * trajec
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_quant_reachside_v(ii, jj, kk, ll, hi, hj, hk, hl, trans)

    !-------------!
    ! Declaration !
    !-------------!
    !- arguments -!
    INTEGER(kind=iprec), INTENT(in) :: ii, jj, kk, ll
    INTEGER(kind=iprec) :: iis, jjs, kks
    REAL(kind=rprec),   INTENT(in)  :: hi, hj, hk, hl, trans

    !- local variables -!
    REAL(kind=rprec)   :: prof  !
    REAL(kind=rprec)    :: temp   !
    REAL(kind=rprec)    :: salt   !
    REAL(kind=rprec)    :: dens  !

    !-------------!
    ! Code begins !
    !-------------!
    IF (key_eco) THEN
      vh(ii,jj,1) = vh(ii,jj,1) + trans
      IF (key_roms.OR.key_mars.OR.key_symphonie) THEN
        prof=fz(gi=hi, gj=hj, gk=-ABS(hk), il=ll)
      ELSE
        prof=fz(gi=hi, gj=hj, gk=-ABS(hk))
      ENDIF
      zvh(ii,jj,1)  = zvh(ii,jj,1)  - prof        * trans
      z2vh(ii,jj,1) = z2vh(ii,jj,1) + prof * prof * trans

      IF (key_alltracers) THEN
        temp  = zinter(tt,hi,hj,hk,hl)
        salt  = zinter(ss,hi,hj,hk,hl)
        IF (key_approximatesigma) THEN
          dens = zinter(rr,hi,hj,hk,hl)
        ELSE
          dens = sigma(zsigma,salt,temp)
        ENDIF
        rvh(ii,jj,1)  = rvh(ii,jj,1)  + dens       * trans
        r2vh(ii,jj,1) = r2vh(ii,jj,1) + dens * dens * trans
        svh(ii,jj,1)  = svh(ii,jj,1)  + salt        * trans
        s2vh(ii,jj,1) = s2vh(ii,jj,1) + salt  * salt  * trans
        tvh(ii,jj,1)  = tvh(ii,jj,1)  + temp        * trans
        t2vh(ii,jj,1) = t2vh(ii,jj,1) + temp  * temp  * trans
      ENDIF

      CALL sub_reducmem_shift_or_not_ind(ii,jj,kk,iis,jjs,kks)

      IF (vv(iis,jjs,kks,ll) > 0._rprec) THEN
        vxy(ii,jj,1) = vxy(ii,jj,1) + trans
        vyz(jj,kk,1) = vyz(jj,kk,1) + trans
      ELSE
        vxy(ii,jj,1) = vxy(ii,jj,1) - trans
        vyz(jj,kk,1) = vyz(jj,kk,1) - trans
      ENDIF

    ELSE

      vht(ii,jj) = vht(ii,jj) + trans
      IF (key_roms.OR.key_mars.OR.key_symphonie) THEN
        prof=fz(gi=hi, gj=hj, gk=-ABS(hk), il=ll)
      ELSE
        prof=fz(gi=hi, gj=hj, gk=-ABS(hk))
      ENDIF
      zvht(ii,jj)  = zvht(ii,jj)  - prof        * trans
      z2vht(ii,jj) = z2vht(ii,jj) + prof * prof * trans

      IF (key_alltracers) THEN
        temp  = zinter(tt,hi,hj,hk,hl)
        salt  = zinter(ss,hi,hj,hk,hl)
        IF (key_approximatesigma) THEN
          dens = zinter(rr,hi,hj,hk,hl)
        ELSE
          dens = sigma(zsigma,salt,temp)
        ENDIF
        rvht(ii,jj)  = rvht(ii,jj)  + dens       * trans
        r2vht(ii,jj) = r2vht(ii,jj) + dens * dens * trans
        svht(ii,jj)  = svht(ii,jj)  + salt        * trans
        s2vht(ii,jj) = s2vht(ii,jj) + salt  * salt  * trans
        tvht(ii,jj)  = tvht(ii,jj)  + temp        * trans
        t2vht(ii,jj) = t2vht(ii,jj) + temp  * temp  * trans
      ENDIF

      CALL sub_reducmem_shift_or_not_ind(ii,jj,kk,iis,jjs,kks)

      IF (vv(iis,jjs,kks,ll) > 0._rprec) THEN
        vxytmp(ii,jj) = vxytmp(ii,jj) + trans
        vyztmp(jj,kk) = vyztmp(jj,kk) + trans
      ELSE
        vxytmp(ii,jj) = vxytmp(ii,jj) - trans
        vyztmp(jj,kk) = vyztmp(jj,kk) - trans
      ENDIF

    ENDIF

  END SUBROUTINE sub_quant_reachside_v

  !!***
  !=========================================================================
  !!****f* mod_quant/sub_quant_reachside_w()
  !! NAME
  !!   sub_quant_reachside_w()
  !!
  !! FUNCTION
  !!   the transport on the initial section is stored in the relevant
  !!   transport arrays.
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
  !!   * INPUTS:
  !!       - ii    : i indice
  !!       - jj    : j indice
  !!       - kk    : k indice
  !!       - ll    : time indice (l)
  !!       - trans : transport
  !!
  !! TODO
  !!   
  !!
  !! USED BY
  !!   * trajec
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_quant_reachside_w(ii, jj, kk, ll, trans)

    !-------------!
    ! Declaration !
    !-------------!
    !- arguments -!
    INTEGER(kind=iprec), INTENT(in) :: ii, jj, kk, ll
    INTEGER(kind=iprec)             :: iis, jjs, kks
    REAL(kind=rprec),    INTENT(in) :: trans

    !-------------!
    ! Code begins !
    !-------------!
    CALL sub_reducmem_shift_or_not_ind(ii,jj,kk,iis,jjs,kks)

    IF (key_eco) THEN
      IF (ww(iis,jjs,kks,ll) > 0._rprec) THEN
        wxz(ii,kk,1) = wxz(ii,kk,1) + trans
        wyz(jj,kk,1) = wyz(jj,kk,1) + trans
      ELSE
        wxz(ii,kk,1) = wxz(ii,kk,1) - trans
        wyz(jj,kk,1) = wyz(jj,kk,1) - trans
      ENDIF
    ELSE

      IF (ww(iis,jjs,kks,ll) > 0._rprec) THEN
        wxztmp(ii,kk) = wxztmp(ii,kk) + trans
        wyztmp(jj,kk) = wyztmp(jj,kk) + trans
      ELSE
        wxztmp(ii,kk) = wxztmp(ii,kk) - trans
        wyztmp(jj,kk) = wyztmp(jj,kk) - trans
      ENDIF
    ENDIF

  END SUBROUTINE sub_quant_reachside_w

  !!***
  !=========================================================================
  !!****f* mod_quant/sub_quant_store_transp()
  !! NAME
  !!   sub_quant_store_transp()
  !!
  !! FUNCTION
  !!   Store transport.
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
  !!   * INPUT:
  !!       - nfin: 
  !!
  !! TODO
  !!   * Cost too much - 55% of the total CPU time !!!!!!!!
  !!     Probably a lot of memory accesses.
  !!   Gprof output:
  !!    %   cumulative   self              self     total           
  !!  time   seconds   seconds    calls   s/call   s/call  name   
  !! 55.00     61.23    61.23     3630     0.02     0.02  mod_quant_mp_sub_quant_store_transp_
  !! 22.97     86.80    25.57     3631     0.01     0.01  mod_quant_mp_sub_quant_inittmp_
  !!
  !! USED BY
  !!   * trajec
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_quant_store_transp(nfin)

    !-------------!
    ! Declaration !
    !-------------!
    !- arguments -!
    INTEGER(kind=iprec), INTENT(in) :: nfin

    !-------------!
    ! Code begins !
    !-------------!
    uxy (:,:,nfin) = uxy (:,:,nfin) + uxytmp(:,:)
    vxy (:,:,nfin) = vxy (:,:,nfin) + vxytmp(:,:)
    uh  (:,:,nfin) = uh  (:,:,nfin) + uht   (:,:)
    vh  (:,:,nfin) = vh  (:,:,nfin) + vht   (:,:)
    zuh (:,:,nfin) = zuh (:,:,nfin) + zuht  (:,:)
    zvh (:,:,nfin) = zvh (:,:,nfin) + zvht  (:,:)
    z2uh(:,:,nfin) = z2uh(:,:,nfin) + z2uht (:,:)
    z2vh(:,:,nfin) = z2vh(:,:,nfin) + z2vht (:,:)

    vyz(:,:,nfin) = vyz(:,:,nfin) + vyztmp(:,:)
    wyz(:,:,nfin) = wyz(:,:,nfin) + wyztmp(:,:)

    uxz(:,:,nfin) = uxz(:,:,nfin) + uxztmp(:,:)
    wxz(:,:,nfin) = wxz(:,:,nfin) + wxztmp(:,:)

    IF (key_alltracers) THEN
      ruh (:,:,nfin) = ruh (:,:,nfin) + ruht (:,:)
      rvh (:,:,nfin) = rvh (:,:,nfin) + rvht (:,:)
      r2uh(:,:,nfin) = r2uh(:,:,nfin) + r2uht(:,:)
      r2vh(:,:,nfin) = r2vh(:,:,nfin) + r2vht(:,:)
      suh (:,:,nfin) = suh (:,:,nfin) + suht (:,:)
      svh (:,:,nfin) = svh (:,:,nfin) + svht (:,:)
      s2uh(:,:,nfin) = s2uh(:,:,nfin) + s2uht(:,:)
      s2vh(:,:,nfin) = s2vh(:,:,nfin) + s2vht(:,:)
      tuh (:,:,nfin) = tuh (:,:,nfin) + tuht (:,:)
      tvh (:,:,nfin) = tvh (:,:,nfin) + tvht (:,:)
      t2uh(:,:,nfin) = t2uh(:,:,nfin) + t2uht(:,:)
      t2vh(:,:,nfin) = t2vh(:,:,nfin) + t2vht(:,:)
    ENDIF

  END SUBROUTINE sub_quant_store_transp

  !!***
  !=========================================================================
  !!****f* mod_quant/sub_quant_statistics()
  !! NAME
  !!   sub_quant_statistics()
  !!
  !! FUNCTION
  !!   Compute some statistics.
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
  !!   * INPUTS :
  !!       - nfin   :
  !!       - ii     :
  !!       - jj     :
  !!       - k      :
  !!       - ll     :
  !!       - trans  :
  !!       - agenew :
  !!       - tcyc   :
  !!       - hi     :
  !!       - hj     :
  !!       - hk     :
  !!       - hl     :
  !!       - fi     :
  !!       - fj     :
  !!       - fk     :
  !!       - fl     :
  !!
  !!   * OUTPUTS :
  !!       - trueage :
  !!       - temp      :
  !!       - salt      :
  !!       - rf      :
  !!
  !! TODO
  !!   
  !!
  !! USED BY
  !!   * trajec
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_quant_statistics(nfin, trans, agenew, tcyc, &
       hi, hj, hk, hl, fi, fj, fk, fl, ii, jj, k, ll,  &
       temp, salt, rf, trueage, tiseq, siseq, riseq)

    !-------------!
    ! Declaration !
    !-------------!
    !- arguments -!
    INTEGER(kind=iprec), INTENT(in)  :: nfin,ii, jj, k, ll
    REAL(kind=rprec)   , INTENT(in)  :: trans, agenew, tcyc, &
         hi, hj, hk, hl, &
         fi, fj, fk, fl
    REAL(kind=rprec)   , INTENT(out) :: trueage, temp, salt, rf
    REAL(kind=rprec)   , OPTIONAL, INTENT(IN) :: tiseq, siseq, riseq


    !- local variables -!
    INTEGER(kind=iprec) :: kk !
    INTEGER(kind=iprec) :: iis, jjs, kks !
    REAL(kind=rprec)    :: ti !
    REAL(kind=rprec)    :: si !
    REAL(kind=rprec)    :: ri !
    REAL(kind=rprec)    :: prof, profm !

    !-------------!
    ! Code begins !
    !-------------!
    isn(nfin) = isn(nfin) + 1
    sn(nfin)  = sn(nfin)  + trans
    trueage   = agenew    / tcyc

    IF (key_roms.OR.key_mars.OR.key_symphonie) THEN
      prof  = fz(gi=hi, gj=hj, gk=hk , il=ll)
      profm = fz(gi=fi, gj=fj, gk=-fk, il=ll)
    ELSE
      prof  = fz(gi=hi, gj=hj, gk=hk)
      profm = fz(gi=fi, gj=fj, gk=-fk)
    ENDIF

    IF (fx(hi,hj) <= s1min(1,nfin)) s1min(1,nfin) = fx(hi,hj)
    IF (fy(hi,hj) <= s1min(2,nfin)) s1min(2,nfin) = fy(hi,hj)
    IF (-prof     <= s1min(3,nfin)) s1min(3,nfin) = -prof
    IF (trueage   <= s1min(4,nfin)) s1min(4,nfin) = trueage

    IF (fx(hi,hj) >= s1max(1,nfin)) s1max(1,nfin) = fx(hi,hj)
    IF (fy(hi,hj) >= s1max(2,nfin)) s1max(2,nfin) = fy(hi,hj)
    IF (-prof     >= s1max(3,nfin)) s1max(3,nfin) = -prof
    IF (trueage   >= s1max(4,nfin)) s1max(4,nfin) = trueage

    s1x(1,nfin)  = s1x(1,nfin)  + fx(hi,hj)             * trans
    s1x(2,nfin)  = s1x(2,nfin)  + fy(hi,hj)             * trans
    s1x(3,nfin)  = s1x(3,nfin)  - prof                  * trans
    s1x(4,nfin)  = s1x(4,nfin)  + trueage               * trans

    s1x2(1,nfin) = s1x2(1,nfin) + fx(hi,hj) * fx(hi,hj) * trans
    s1x2(2,nfin) = s1x2(2,nfin) + fy(hi,hj) * fy(hi,hj) * trans
    s1x2(3,nfin) = s1x2(3,nfin) + prof      * prof      * trans
    s1x2(4,nfin) = s1x2(4,nfin) + trueage   * trueage   * trans

    IF (fx(fi,fj) <= s0min(1,nfin)) s0min(1,nfin) = fx(fi,fj)
    IF (fy(fi,fj) <= s0min(2,nfin)) s0min(2,nfin) = fy(fi,fj)
    IF (-profm    <= s0min(3,nfin)) s0min(3,nfin) = -profm
    IF (trueage   <= s0min(4,nfin)) s0min(4,nfin) = trueage

    IF (fx(fi,fj) >= s0max(1,nfin)) s0max(1,nfin) = fx(fi,fj)
    IF (fy(fi,fj) >= s0max(2,nfin)) s0max(2,nfin) = fy(fi,fj)
    IF (-profm    >= s0max(3,nfin)) s0max(3,nfin) = -profm
    IF (trueage   >= s0max(4,nfin)) s0max(4,nfin) = trueage

    s0x(1,nfin)  = s0x(1,nfin)  + fx(fi,fj)             * trans
    s0x(2,nfin)  = s0x(2,nfin)  + fy(fi,fj)             * trans
    s0x(3,nfin)  = s0x(3,nfin)  - profm                 * trans
    s0x(4,nfin)  = s0x(4,nfin)  + trueage               * trans

    s0x2(1,nfin) = s0x2(1,nfin) + fx(fi,fj) * fx(fi,fj) * trans
    s0x2(2,nfin) = s0x2(2,nfin) + fy(fi,fj) * fy(fi,fj) * trans
    s0x2(3,nfin) = s0x2(3,nfin) + profm     * profm     * trans
    s0x2(4,nfin) = s0x2(4,nfin) + trueage   * trueage   * trans

    kk = k
    IF (k <= 0) kk = k + 1

    IF (key_alltracers) THEN

      IF (key_nointerpolstats) THEN
        CALL sub_reducmem_shift_or_not_ind(ii,jj,kk,iis,jjs,kks) !!NG: bug correction 14/04/2009
        temp = tt(iis,jjs,kks,ll) !!NG: bug correction 14/04/2009
        salt = ss(iis,jjs,kks,ll) !!NG: bug correction 14/04/2009
        rf   = rr(iis,jjs,kks,ll) !!NG: bug correction 14/04/2009
        ti   = 0._rprec
        si   = 0._rprec
        ri   = 0._rprec
      ELSE
        temp = zinter(tt,hi,hj,hk,hl)
        salt = zinter(ss,hi,hj,hk,hl)
        IF (key_approximatesigma) THEN
          rf = zinter(rr,hi,hj,hk,hl)
        ELSE
          rf = sigma(zsigma,salt,temp)
        ENDIF
        IF (key_sequential) THEN
          IF (PRESENT(tiseq).AND.PRESENT(siseq).AND.PRESENT(riseq)) THEN
            ti = tiseq
            si = siseq
            ri = riseq
          ELSE
            STOP
          ENDIF
        ELSE
          ti = zinter(tt,fi,fj,-fk,fl)
          si = zinter(ss,fi,fj,-fk,fl)
          IF (key_approximatesigma) THEN
            ri = zinter(rr,fi,fj,-fk,fl)
          ELSE
            ri = sigma(zsigma,si,ti)
          ENDIF
        ENDIF
      ENDIF

!!NG      IF (.NOT.key_nointerpolstats) THEN
!!NG        temp = zinter(tt,hi,hj,hk,hl)
!!NG        salt = zinter(ss,hi,hj,hk,hl)
!!NG        IF (key_approximatesigma) THEN
!!NG          rf = zinter(rr,hi,hj,hk,hl)
!!NG        ELSE
!!NG          rf = sigma(zsigma,salt,temp)
!!NG        ENDIF
!!NG        ti = zinter(tt,fi,fj,-fk,fl)
!!NG        si = zinter(ss,fi,fj,-fk,fl)
!!NG        IF (key_approximatesigma) THEN
!!NG          ri = zinter(rr,fi,fj,-fk,fl)
!!NG        ELSE
!!NG          ri = sigma(zsigma,si,ti)
!!NG        ENDIF
!!NG      ELSE
!!NG        temp = tt(ii,jj,kk,ll)
!!NG        salt = ss(ii,jj,kk,ll)
!!NG        rf = rr(ii,jj,kk,ll)
!!NG        ti = 0.
!!NG        si = 0.
!!NG        ri = 0.
!!NG      ENDIF

      IF (temp <= s1min(5,nfin)) s1min(5,nfin) = temp
      IF (salt <= s1min(6,nfin)) s1min(6,nfin) = salt
      IF (rf   <= s1min(7,nfin)) s1min(7,nfin) = rf
      IF (temp >= s1max(5,nfin)) s1max(5,nfin) = temp
      IF (salt >= s1max(6,nfin)) s1max(6,nfin) = salt
      IF (rf   >= s1max(7,nfin)) s1max(7,nfin) = rf

      s1x(5,nfin)  = s1x(5,nfin)  + temp        * trans
      s1x(6,nfin)  = s1x(6,nfin)  + salt        * trans
      s1x(7,nfin)  = s1x(7,nfin)  + rf          * trans
      s1x2(5,nfin) = s1x2(5,nfin) + temp * temp * trans
      s1x2(6,nfin) = s1x2(6,nfin) + salt * salt * trans
      s1x2(7,nfin)=  s1x2(7,nfin) + rf   * rf   * trans

      IF (ti <= s0min(5,nfin)) s0min(5,nfin) = ti
      IF (si <= s0min(6,nfin)) s0min(6,nfin) = si
      IF (ri <= s0min(7,nfin)) s0min(7,nfin) = ri
      IF (ti >= s0max(5,nfin)) s0max(5,nfin) = ti
      IF (si >= s0max(6,nfin)) s0max(6,nfin) = si
      IF (ri >= s0max(7,nfin)) s0max(7,nfin) = ri

      s0x(5,nfin)  = s0x(5,nfin)  + ti      * trans
      s0x(6,nfin)  = s0x(6,nfin)  + si      * trans
      s0x(7,nfin)  = s0x(7,nfin)  + ri      * trans
      s0x2(5,nfin) = s0x2(5,nfin) + ti * ti * trans
      s0x2(6,nfin) = s0x2(6,nfin) + si * si * trans
      s0x2(7,nfin) = s0x2(7,nfin) + ri * ri * trans

    ENDIF

  END SUBROUTINE sub_quant_statistics

  !!***
  !=========================================================================
  !!****f* mod_quant/sub_quant_alloc()
  !! NAME
  !!   sub_quant_alloc()
  !!
  !! FUNCTION
  !!   Dynamic allocation of arrays.
  !!
  !! AUTHOR
  !!   * Origin  : Nicolas Grima (April-May 2005)
  !! 
  !! CREATION DATE
  !!   * April-May 2005
  !!
  !! HISTORY
  !!   Date (dd/mm/yyyy/) - Modification(s)
  !!   October 2005 - zoom
  !!
  !! ARGUMENTS
  !!   * No arguments
  !!
  !! TODO
  !!   
  !!
  !! USED BY
  !!   * trajec
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_quant_alloc(dimx, dimy, dimk)

    INTEGER(kind = iprec), INTENT(IN) :: dimx, dimy, dimk

    !-------------!
    ! Code begins !
    !-------------!
    IF (key_eco) THEN
      ALLOCATE(uxy(dimx,dimy,1))
      CALL sub_memory(size(uxy),'r','uxy','sub_quant_alloc')
      ALLOCATE(vxy(dimx,dimy,1))
      CALL sub_memory(size(vxy),'r','vxy','sub_quant_alloc')

      ALLOCATE(uxz(dimx,dimk,1))
      CALL sub_memory(size(uxz),'r','uxz','sub_quant_alloc')
      ALLOCATE(wxz(dimx,dimk,1))
      CALL sub_memory(size(wxz),'r','wxz','sub_quant_alloc')

      ALLOCATE(vyz(dimy,dimk,1))
      CALL sub_memory(size(vyz),'r','vyz','sub_quant_alloc')
      ALLOCATE(wyz(dimy,dimk,1))
      CALL sub_memory(size(vyz),'r','vyz','sub_quant_alloc')

      ALLOCATE(uh(dimx,dimy,1))
      CALL sub_memory(size(uh),'r','uh','sub_quant_alloc')
      ALLOCATE(vh(dimx,dimy,1))
      CALL sub_memory(size(vh),'r','vh','sub_quant_alloc')

      ALLOCATE(zuh(dimx,dimy,1))
      CALL sub_memory(size(zuh),'r','zuh','sub_quant_alloc')
      ALLOCATE(zvh(dimx,dimy,1))
      CALL sub_memory(size(zvh),'r','zvh','sub_quant_alloc')

      ALLOCATE(z2uh(dimx,dimy,1))
      CALL sub_memory(size(z2uh),'r','z2uh','sub_quant_alloc')
      ALLOCATE(z2vh(dimx,dimy,1))
      CALL sub_memory(size(z2vh),'r','z2vh','sub_quant_alloc')

      IF (key_alltracers) THEN
        ALLOCATE(ruh(dimx,dimy,1))
        CALL sub_memory(size(ruh),'r','ruh','sub_quant_alloc')
        ALLOCATE(rvh(dimx,dimy,1))
        CALL sub_memory(size(rvh),'r','rvh','sub_quant_alloc')
        ALLOCATE(r2uh(dimx,dimy,1))
        CALL sub_memory(size(r2uh),'r','r2uh','sub_quant_alloc')
        ALLOCATE(r2vh(dimx,dimy,1))
        CALL sub_memory(size(r2vh),'r','r2vh','sub_quant_alloc')
        ALLOCATE(suh(dimx,dimy,1))
        CALL sub_memory(size(suh),'r','suh','sub_quant_alloc')
        ALLOCATE(svh(dimx,dimy,1))
        CALL sub_memory(size(svh),'r','svh','sub_quant_alloc')
        ALLOCATE(s2uh(dimx,dimy,1))
        CALL sub_memory(size(s2uh),'r','s2uh','sub_quant_alloc')
        ALLOCATE(s2vh(dimx,dimy,1))
        CALL sub_memory(size(s2vh),'r','s2vh','sub_quant_alloc')
        ALLOCATE(tuh(dimx,dimy,1))
        CALL sub_memory(size(tuh),'r','tuh','sub_quant_alloc')
        ALLOCATE(tvh(dimx,dimy,1))
        CALL sub_memory(size(tvh),'r','tvh','sub_quant_alloc')
        ALLOCATE(t2uh(dimx,dimy,1))
        CALL sub_memory(size(t2uh),'r','t2uh','sub_quant_alloc')
        ALLOCATE(t2vh(dimx,dimy,1))
        CALL sub_memory(size(t2vh),'r','t2vh','sub_quant_alloc')
      ENDIF
    ELSE
      ALLOCATE(uxy(dimx,dimy,0:maxsect+1))
      CALL sub_memory(size(uxy),'r','uxy','sub_quant_alloc')
      ALLOCATE(vxy(dimx,dimy,0:maxsect+1))
      CALL sub_memory(size(vxy),'r','vxy','sub_quant_alloc')
      ALLOCATE(uxytmp(dimx,dimy))
      CALL sub_memory(size(uxytmp),'r','uxytmp','sub_quant_alloc')
      ALLOCATE(vxytmp(dimx,dimy))
      CALL sub_memory(size(vxytmp),'r','vxytmp','sub_quant_alloc')

      ALLOCATE(uxz(dimx,dimk,0:maxsect+1))
      CALL sub_memory(size(uxz),'r','uxz','sub_quant_alloc')
      ALLOCATE(wxz(dimx,dimk,0:maxsect+1))
      CALL sub_memory(size(wxz),'r','wxz','sub_quant_alloc')
      ALLOCATE(uxztmp(dimx,dimk))
      CALL sub_memory(size(uxztmp),'r','uxztmp','sub_quant_alloc')
      ALLOCATE(wxztmp(dimx,dimk))
      CALL sub_memory(size(wxztmp),'r','wxztmp','sub_quant_alloc')

      ALLOCATE(vyz(dimy,dimk,0:maxsect+1))
      CALL sub_memory(size(vyz),'r','vyz','sub_quant_alloc')
      ALLOCATE(wyz(dimy,dimk,0:maxsect+1))
      CALL sub_memory(size(wyz),'r','wyz','sub_quant_alloc')
      ALLOCATE(vyztmp(dimy,dimk))
      CALL sub_memory(size(vyztmp),'r','vyztmp','sub_quant_alloc')
      ALLOCATE(wyztmp(dimy,dimk))
      CALL sub_memory(size(wyztmp),'r','wyztmp','sub_quant_alloc')

      ALLOCATE(uh(dimx,dimy,0:maxsect+1))
      CALL sub_memory(size(uh),'r','uh','sub_quant_alloc')
      ALLOCATE(vh(dimx,dimy,0:maxsect+1))
      CALL sub_memory(size(vh),'r','vh','sub_quant_alloc')
      ALLOCATE(uht(dimx,dimy))
      CALL sub_memory(size(uht),'r','uht','sub_quant_alloc')
      ALLOCATE(vht(dimx,dimy))
      CALL sub_memory(size(vht),'r','vht','sub_quant_alloc')


      ALLOCATE(zuh(dimx,dimy,0:maxsect+1))
      CALL sub_memory(size(zuh),'r','zuh','sub_quant_alloc')
      ALLOCATE(zvh(dimx,dimy,0:maxsect+1))
      CALL sub_memory(size(zvh),'r','zvh','sub_quant_alloc')
      ALLOCATE(z2uh(dimx,dimy,0:maxsect+1))
      CALL sub_memory(size(z2uh),'r','z2uh','sub_quant_alloc')
      ALLOCATE(z2vh(dimx,dimy,0:maxsect+1))
      CALL sub_memory(size(z2vh),'r','z2vh','sub_quant_alloc')
      ALLOCATE(zuht(dimx,dimy))
      CALL sub_memory(size(zuht),'r','zuht','sub_quant_alloc')
      ALLOCATE(zvht(dimx,dimy))
      CALL sub_memory(size(zvht),'r','zvht','sub_quant_alloc')
      ALLOCATE(z2uht(dimx,dimy))
      CALL sub_memory(size(z2uht),'r','z2uht','sub_quant_alloc')
      ALLOCATE(z2vht(dimx,dimy))
      CALL sub_memory(size(z2vht),'r','z2vht','sub_quant_alloc')

      IF (key_alltracers) THEN

        ALLOCATE(ruh(dimx,dimy,0:maxsect+1))
        CALL sub_memory(size(ruh),'r','ruh','sub_quant_alloc')
        ALLOCATE(rvh(dimx,dimy,0:maxsect+1))
        CALL sub_memory(size(rvh),'r','rvh','sub_quant_alloc')
        ALLOCATE(r2uh(dimx,dimy,0:maxsect+1))
        CALL sub_memory(size(r2uh),'r','r2uh','sub_quant_alloc')
        ALLOCATE(r2vh(dimx,dimy,0:maxsect+1))
        CALL sub_memory(size(r2vh),'r','r2vh','sub_quant_alloc')
        ALLOCATE(suh(dimx,dimy,0:maxsect+1))
        CALL sub_memory(size(suh),'r','suh','sub_quant_alloc')
        ALLOCATE(svh(dimx,dimy,0:maxsect+1))
        CALL sub_memory(size(svh),'r','svh','sub_quant_alloc')
        ALLOCATE(s2uh(dimx,dimy,0:maxsect+1))
        CALL sub_memory(size(s2uh),'r','s2uh','sub_quant_alloc')
        ALLOCATE(s2vh(dimx,dimy,0:maxsect+1))
        CALL sub_memory(size(s2vh),'r','s2vh','sub_quant_alloc')
        ALLOCATE(tuh(dimx,dimy,0:maxsect+1))
        CALL sub_memory(size(tuh),'r','tuh','sub_quant_alloc')
        ALLOCATE(tvh(dimx,dimy,0:maxsect+1))
        CALL sub_memory(size(tvh),'r','tvh','sub_quant_alloc')
        ALLOCATE(t2uh(dimx,dimy,0:maxsect+1))
        CALL sub_memory(size(t2uh),'r','t2uh','sub_quant_alloc')
        ALLOCATE(t2vh(dimx,dimy,0:maxsect+1))
        CALL sub_memory(size(t2vh),'r','t2vh','sub_quant_alloc')
        ALLOCATE(ruht(dimx,dimy))
        CALL sub_memory(size(ruht),'r','ruht','sub_quant_alloc')
        ALLOCATE(rvht(dimx,dimy))
        CALL sub_memory(size(rvht),'r','rvht','sub_quant_alloc')
        ALLOCATE(r2uht(dimx,dimy))
        CALL sub_memory(size(r2uht),'r','r2uht','sub_quant_alloc')
        ALLOCATE(r2vht(dimx,dimy))
        CALL sub_memory(size(r2vht),'r','r2vht','sub_quant_alloc')
        ALLOCATE(suht(dimx,dimy))
        CALL sub_memory(size(suht),'r','suht','sub_quant_alloc')
        ALLOCATE(svht(dimx,dimy))
        CALL sub_memory(size(svht),'r','svht','sub_quant_alloc')
        ALLOCATE(s2uht(dimx,dimy))
        CALL sub_memory(size(s2uht),'r','s2uht','sub_quant_alloc')
        ALLOCATE(s2vht(dimx,dimy))
        CALL sub_memory(size(s2vht),'r','s2vht','sub_quant_alloc')
        ALLOCATE(tuht(dimx,dimy))
        CALL sub_memory(size(tuht),'r','tuht','sub_quant_alloc')
        ALLOCATE(tvht(dimx,dimy))
        CALL sub_memory(size(tvht),'r','tvht','sub_quant_alloc')
        ALLOCATE(t2uht(dimx,dimy))
        CALL sub_memory(size(t2uht),'r','t2uht','sub_quant_alloc')
        ALLOCATE(t2vht(dimx,dimy))
        CALL sub_memory(size(t2vht),'r','t2vht','sub_quant_alloc')
      ENDIF
    ENDIF
  END SUBROUTINE sub_quant_alloc

  !!***
  !=========================================================================
  !!****f* mod_quant/sub_quant_init()
  !! NAME
  !!   sub_quant_init()
  !!
  !! FUNCTION
  !!   Array initalizations.
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
  !!    Cost too much - 23% of the total CPU time !!!!!!
  !!
  !! USED BY
  !!   * trajec
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_quant_init()

    !-------------!
    ! Code begins !
    !-------------!

    CALL sub_quant_initorig()

    IF (.NOT.key_eco) THEN
      CALL sub_quant_inittmp()
    ENDIF

  END SUBROUTINE sub_quant_init

  !!***
  !=========================================================================
  !!****f* mod_quant/sub_quant_initorig()
  !! NAME
  !!   sub_quant_initorig()
  !!
  !! FUNCTION
  !!   First part of array initializations.
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
  !!   * sub_quant_init (mod_quant.f90)
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_quant_initorig()

    !-------------!
    ! Code begins !
    !-------------!
    uxy(:,:,:) = 0._rprec 
    vxy(:,:,:) = 0._rprec

    uxz(:,:,:) = 0._rprec
    wxz(:,:,:) = 0._rprec

    vyz(:,:,:) = 0._rprec
    wyz(:,:,:) = 0._rprec

    uh(:,:,:) = 0._rprec
    vh(:,:,:) = 0._rprec

    zuh(:,:,:)  = 0._rprec
    zvh(:,:,:)  = 0._rprec

    z2uh(:,:,:) = 0._rprec
    z2vh(:,:,:) = 0._rprec

    IF (key_alltracers) THEN

      ruh(:,:,:)  = 0._rprec
      rvh(:,:,:)  = 0._rprec
      r2uh(:,:,:) = 0._rprec
      r2vh(:,:,:) = 0._rprec
      suh(:,:,:)  = 0._rprec
      svh(:,:,:)  = 0._rprec
      s2uh(:,:,:) = 0._rprec
      s2vh(:,:,:) = 0._rprec
      tuh(:,:,:)  = 0._rprec
      tvh(:,:,:)  = 0._rprec
      t2uh(:,:,:) = 0._rprec
      t2vh(:,:,:) = 0._rprec
    ENDIF
  END SUBROUTINE sub_quant_initorig

  !!***
  !=========================================================================
  !!****f* mod_quant/sub_quant_inittmp()
  !! NAME
  !!   sub_quant_inittmp()
  !!
  !! FUNCTION
  !!   Second part of array initializations.
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
  !!   *
  !!
  !! USED BY
  !!   * sub_quant_init (mod_quant.f90)
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_quant_inittmp()

    !-------------!
    ! Code begins !
    !-------------!
!!NG:09/01/2007    IF (.NOT.key_eco) THEN

    uxytmp(:,:) = 0._rprec !
    vxytmp(:,:) = 0._rprec !

    uxztmp(:,:) = 0._rprec !
    wxztmp(:,:) = 0._rprec !

    vyztmp(:,:) = 0._rprec !
    wyztmp(:,:) = 0._rprec !

    uht(:,:)    = 0._rprec !
    vht(:,:)    = 0._rprec !
    zuht(:,:)   = 0._rprec !
    zvht(:,:)   = 0._rprec !
    z2uht(:,:)  = 0._rprec !
    z2vht(:,:)  = 0._rprec !

    IF (key_alltracers) THEN
      ruht(:,:)   = 0._rprec !
      rvht(:,:)   = 0._rprec !
      r2uht(:,:)  = 0._rprec !
      r2vht(:,:)  = 0._rprec !
      suht(:,:)   = 0._rprec !
      svht(:,:)   = 0._rprec !
      s2uht(:,:)  = 0._rprec !
      s2vht(:,:)  = 0._rprec !
      tuht(:,:)   = 0._rprec !
      tvht(:,:)   = 0._rprec !
      t2uht(:,:)  = 0._rprec !
      t2vht(:,:)  = 0._rprec !
    ENDIF

!!NG:09/01/2007    ENDIF

  END SUBROUTINE sub_quant_inittmp

  !!***
  !=========================================================================
  !!****f* mod_quant/sub_quant_dealloc()
  !! NAME
  !!   sub_quant_dealloc()
  !!
  !! FUNCTION
  !!   Deallocate dynamic array memory.
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
  !!   * No argument.
  !!
  !! TODO
  !!   
  !!
  !! USED BY
  !!   * trajec
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_quant_dealloc()

    !-------------!
    ! Code begins !
    !-------------!
    IF (key_eco) THEN
       CALL sub_memory(-size(uxy),'r','uxy','sub_quant_dealloc')
       DEALLOCATE(uxy)
       CALL sub_memory(-size(vxy),'r','vxy','sub_quant_dealloc')
       DEALLOCATE(vxy)

       CALL sub_memory(-size(uxz),'r','uxz','sub_quant_dealloc')
       DEALLOCATE(uxz)
       CALL sub_memory(-size(wxz),'r','wxz','sub_quant_dealloc')
       DEALLOCATE(wxz)

       CALL sub_memory(-size(vyz),'r','vyz','sub_quant_dealloc')
       DEALLOCATE(vyz)
       CALL sub_memory(-size(wyz),'r','wyz','sub_quant_dealloc')
       DEALLOCATE(wyz)

       CALL sub_memory(-size(uh),'r','uh','sub_quant_dealloc')
       DEALLOCATE(uh)
       CALL sub_memory(-size(vh),'r','vh','sub_quant_dealloc')
       DEALLOCATE(vh)
       CALL sub_memory(-size(zuh),'r','zuh','sub_quant_dealloc')
       DEALLOCATE(zuh)
       CALL sub_memory(-size(zvh),'r','zvh','sub_quant_dealloc')
       DEALLOCATE(zvh)
       CALL sub_memory(-size(z2uh),'r','z2uh','sub_quant_dealloc')
       DEALLOCATE(z2uh)
       CALL sub_memory(-size(z2vh),'r','z2vh','sub_quant_dealloc')
       DEALLOCATE(z2vh)


       IF (key_alltracers) THEN
          CALL sub_memory(-size(ruh),'r','ruh','sub_quant_dealloc')
          DEALLOCATE(ruh)
          CALL sub_memory(-size(rvh),'r','rvh','sub_quant_dealloc')
          DEALLOCATE(rvh)
          CALL sub_memory(-size(r2uh),'r','r2uh','sub_quant_dealloc')
          DEALLOCATE(r2uh)
          CALL sub_memory(-size(r2vh),'r','r2vh','sub_quant_dealloc')
          DEALLOCATE(r2vh)
          CALL sub_memory(-size(suh),'r','suh','sub_quant_dealloc')
          DEALLOCATE(suh)
          CALL sub_memory(-size(svh),'r','svh','sub_quant_dealloc')
          DEALLOCATE(svh)
          CALL sub_memory(-size(s2uh),'r','s2uh','sub_quant_dealloc')
          DEALLOCATE(s2uh)
          CALL sub_memory(-size(s2vh),'r','s2vh','sub_quant_dealloc')
          DEALLOCATE(s2vh)
          CALL sub_memory(-size(tuh),'r','tuh','sub_quant_dealloc')
          DEALLOCATE(tuh)
          CALL sub_memory(-size(tvh),'r','tvh','sub_quant_dealloc')
          DEALLOCATE(tvh)
          CALL sub_memory(-size(t2uh),'r','t2uh','sub_quant_dealloc')
          DEALLOCATE(t2uh)
          CALL sub_memory(-size(t2vh),'r','t2vh','sub_quant_dealloc')
          DEALLOCATE(t2vh)
       ENDIF
    ELSE

       CALL sub_memory(-size(uxy),'r','uxy','sub_quant_dealloc')
       DEALLOCATE(uxy)
       CALL sub_memory(-size(vxy),'r','vxy','sub_quant_dealloc')
       DEALLOCATE(vxy)
       CALL sub_memory(-size(uxytmp),'r','uxytmp','sub_quant_dealloc')
       DEALLOCATE(uxytmp)
       CALL sub_memory(-size(vxytmp),'r','vxytmp','sub_quant_dealloc')
       DEALLOCATE(vxytmp)

       CALL sub_memory(-size(uxz),'r','uxz','sub_quant_dealloc')
       DEALLOCATE(uxz)
       CALL sub_memory(-size(wxz),'r','wxz','sub_quant_dealloc')
       DEALLOCATE(wxz)
       CALL sub_memory(-size(uxztmp),'r','uxztmp','sub_quant_dealloc')
       DEALLOCATE(uxztmp)
       CALL sub_memory(-size(wxztmp),'r','wxztmp','sub_quant_dealloc')
       DEALLOCATE(wxztmp)

       CALL sub_memory(-size(vyz),'r','vyz','sub_quant_dealloc')
       DEALLOCATE(vyz)
       CALL sub_memory(-size(wyz),'r','wyz','sub_quant_dealloc')
       DEALLOCATE(wyz)
       CALL sub_memory(-size(vyztmp),'r','vyztmp','sub_quant_dealloc')
       DEALLOCATE(vyztmp)
       CALL sub_memory(-size(wyztmp),'r','wyztmp','sub_quant_dealloc')
       DEALLOCATE(wyztmp)

       CALL sub_memory(-size(uh),'r','uh','sub_quant_dealloc')
       DEALLOCATE(uh)
       CALL sub_memory(-size(vh),'r','vh','sub_quant_dealloc')
       DEALLOCATE(vh)
       CALL sub_memory(-size(zuh),'r','zuh','sub_quant_dealloc')
       DEALLOCATE(zuh)
       CALL sub_memory(-size(zvh),'r','zvh','sub_quant_dealloc')
       DEALLOCATE(zvh)
       CALL sub_memory(-size(z2uh),'r','z2uh','sub_quant_dealloc')
       DEALLOCATE(z2uh)
       CALL sub_memory(-size(z2vh),'r','z2vh','sub_quant_dealloc')
       DEALLOCATE(z2vh)

       CALL sub_memory(-size(uht),'r','uht','sub_quant_dealloc')
       DEALLOCATE(uht)
       CALL sub_memory(-size(vht),'r','vht','sub_quant_dealloc')
       DEALLOCATE(vht)
       CALL sub_memory(-size(zuht),'r','zuht','sub_quant_dealloc')
       DEALLOCATE(zuht)
       CALL sub_memory(-size(zvht),'r','zvht','sub_quant_dealloc')
       DEALLOCATE(zvht)
       CALL sub_memory(-size(z2uht),'r','z2uht','sub_quant_dealloc')
       DEALLOCATE(z2uht)
       CALL sub_memory(-size(z2vht),'r','z2vht','sub_quant_dealloc')
       DEALLOCATE(z2vht)

       IF (key_alltracers) THEN
          CALL sub_memory(-size(ruh),'r','ruh','sub_quant_dealloc')
          DEALLOCATE(ruh)
       CALL sub_memory(-size(rvh),'r','rvh','sub_quant_dealloc')
          DEALLOCATE(rvh)
       CALL sub_memory(-size(r2uh),'r','r2uh','sub_quant_dealloc')
          DEALLOCATE(r2uh)
       CALL sub_memory(-size(r2vh),'r','r2vh','sub_quant_dealloc')
          DEALLOCATE(r2vh)
       CALL sub_memory(-size(suh),'r','suh','sub_quant_dealloc')
          DEALLOCATE(suh)
          CALL sub_memory(-size(svh),'r','svh','sub_quant_dealloc')
          DEALLOCATE(svh)
           CALL sub_memory(-size(s2uh),'r','s2uh','sub_quant_dealloc')
          DEALLOCATE(s2uh)
          CALL sub_memory(-size(s2vh),'r','s2vh','sub_quant_dealloc')
          DEALLOCATE(s2vh)
          CALL sub_memory(-size(tuh),'r','tuh','sub_quant_dealloc')
          DEALLOCATE(tuh)
          CALL sub_memory(-size(tvh),'r','tvh','sub_quant_dealloc')
          DEALLOCATE(tvh)
          CALL sub_memory(-size(t2uh),'r','t2uh','sub_quant_dealloc')
          DEALLOCATE(t2uh)
          CALL sub_memory(-size(t2vh),'r','t2vh','sub_quant_dealloc')
          DEALLOCATE(t2vh)

          CALL sub_memory(-size(ruht),'r','ruht','sub_quant_dealloc')
          DEALLOCATE(ruht)
          CALL sub_memory(-size(rvht),'r','rvht','sub_quant_dealloc')
          DEALLOCATE(rvht)
          CALL sub_memory(-size(r2uht),'r','r2uht','sub_quant_dealloc')
          DEALLOCATE(r2uht)
          CALL sub_memory(-size(r2vht),'r','r2vht','sub_quant_dealloc')
          DEALLOCATE(r2vht)
          CALL sub_memory(-size(suht),'r','suht','sub_quant_dealloc')
          DEALLOCATE(suht)
          CALL sub_memory(-size(svht),'r','svht','sub_quant_dealloc')
          DEALLOCATE(svht)
          CALL sub_memory(-size(s2uht),'r','s2uht','sub_quant_dealloc')
          DEALLOCATE(s2uht)
          CALL sub_memory(-size(s2vht),'r','s2vht','sub_quant_dealloc')
          DEALLOCATE(s2vht)
          CALL sub_memory(-size(tuht),'r','tuht','sub_quant_dealloc')
          DEALLOCATE(tuht)
          CALL sub_memory(-size(tvht),'r','tvht','sub_quant_dealloc')
          DEALLOCATE(tvht)
          CALL sub_memory(-size(t2uht),'r','t2uht','sub_quant_dealloc')
          DEALLOCATE(t2uht)
          CALL sub_memory(-size(t2vht),'r','t2vht','sub_quant_dealloc')
          DEALLOCATE(t2vht)
       ENDIF
    ENDIF

  END SUBROUTINE sub_quant_dealloc

  !!***
END MODULE mod_quant
