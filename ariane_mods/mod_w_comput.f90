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
!!****h* ariane/mod_w_comput
!! NAME
!!   mod_w_comput (mod_w_comput.f90 - Fortran90 module)
!!
!! USAGE
!!   Include 'USE mod_w_comput' in the header of your Fortran 90 source 
!!   code.
!!   Then you'll have access to the subroutine:
!!      - sub_w_comput
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
!!
!! USES
!!   * USE mod_precision
!!   * USE mod_namelist
!!
!! USED BY
!!   * mod_input_data.f90
!!
!! SOURCE
!!=========================================================================
MODULE mod_w_comput

  !------------------!
  ! USE ASSOCIAITION !
  !------------------!
  USE mod_precision
  USE mod_namelist
  USE mod_input_grid

  !-------------!
  ! DECLARATION !
  !-------------!
  IMPLICIT NONE

CONTAINS
  !!***
  !=========================================================================
  !!****f* mod_w_comput/sub_w_comput()
  !! NAME
  !!   sub_w_comput()
  !!
  !! FUNCTION
  !!   Compute vertical current from zonal and meridional currents.
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
  !!      - uu:
  !!      - vv:
  !!
  !!   * OUTPUTS:
  !!      - ww:
  !!
  !! TODO
  !!   
  !!
  !! USED BY
  !!   * mod_input_data.f90
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_w_comput(uu,vv,ww)

    !-------------!
    ! Declaration !
    !-------------!
    !- arguments -!
    REAL(kind = rprec), DIMENSION(:,:,:,:), INTENT(in)  :: uu, vv
    REAL(kind = rprec), DIMENSION(:,:,:,:), INTENT(out) :: ww

    !- local variables -!
    REAL(kind = rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: div
    REAL(kind = rprec), DIMENSION(:,:,:), ALLOCATABLE :: stretch

    INTEGER(kind = iprec) :: isize
    INTEGER(kind = iprec) :: jsize
    INTEGER(kind = iprec) :: ksize
    INTEGER(kind = iprec) :: lsize
    INTEGER(kind=iprec) :: ii
    INTEGER(kind=iprec) :: jj
    INTEGER(kind=iprec) :: kk
    INTEGER(kind=iprec) :: ll
    !-------------!
    ! Code begins !
    !-------------!
    !- intialization -!
    isize= SIZE(uu, dim=1)
    jsize= SIZE(uu, dim=2)
    ksize= SIZE(uu, dim=3)
    lsize= SIZE(uu, dim=4)


    ALLOCATE(div(1:isize,1:jsize,1:ksize,1:lsize))
    call sub_memory(size(div),'r','div','sub_comput_w')
    ALLOCATE(stretch(1:isize,1:jsize,1:lsize))
    call sub_memory(size(stretch),'r','stretch','sub_comput_w')

    !=============================================
    ! the 3D velocity field MUST be NON-DIVERGENT unless key_vvl
    !=============================================
    ww(:,:,:,:)=0._rprec

    IF (key_2D) THEN
      ww(:,:,2,:)=-EPSILON(1.0_rprec)
    ELSE
      div(2:isize,2:jsize,:,:) = uu(1:isize-1,2:jsize,:,:) - uu(2:isize,2:jsize,:,:) &
           + vv(2:isize,1:jsize-1,:,:) - vv(2:isize,2:jsize,:,:)
      if (key_vvl) then
         stretch = SUM(div, DIM=3)
         DO ll = 1,lsize
            DO jj = 2,jsize
               DO ii = 2,isize
                  DO kk = NINT(mbathy(ii,jj,1,1)),1,-1
                     ww(ii,jj,kk,ll) = ww(ii,jj,kk+1,ll) + div(ii,jj,kk,ll) &
                          - stretch(ii,jj,ll) * e3t0(ii,jj,kk,ll)/totaldepth(ii,jj,1,1)
                  ENDDO
                  ww(ii,jj,1,ll) = 0.
               ENDDO
            ENDDO
         ENDDO
      else
      DO ll = 1,lsize
        DO kk = ksize-1,1,-1
          DO jj = 2,jsize
            DO ii = 2,isize
              ww(ii,jj,kk,ll) = ww(ii,jj,kk+1,ll) + div(ii,jj,kk,ll)
            END DO
          END DO
        END DO
     END DO
     endif
    ENDIF

    !!NG     DO ll = 1,lsize
    !!NG       DO kk = 1,ksize-1
    !!NG         DO jj = 2,jsize
    !!NG           DO ii = 2,isize
    !!NG             ww(ii,jj,kk+1,ll)=ww(ii,jj,kk,ll) &
    !!NG                  - uu(ii-1,jj  ,kk,ll) + uu(ii,jj,kk,ll) &
    !!NG                  - vv(ii  ,jj-1,kk,ll) + vv(ii,jj,kk,ll)
    !!NG           END DO
    !!NG         END DO
    !!NG       END DO
    !!NG     END DO

    IF (key_periodic) THEN
      ww(1,:,:,:)=ww(isize-1,:,:,:)
    ENDIF

    DEALLOCATE(div)
    call sub_memory(-size(div),'r','div','sub_comput_w')
    call sub_memory(-size(stretch),'r','stretch','sub_comput_w')

  END SUBROUTINE sub_w_comput
  !!***
END MODULE mod_w_comput