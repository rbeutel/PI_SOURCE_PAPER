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
!!****h* ariane/mod_fz
!! NAME
!!   mod_fz (mod_fz.f90 - Fortran90 module)
!!
!! USAGE
!!   Include 'USE mod_fz' in the header of your Fortran 90 source 
!!   code.
!!   Then you'll have access to the subroutine:
!!      - fz
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
!!   * USE mod_fz
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
!!   * mod_fx
!!   * mod_fy
!!
!! USES
!!   * USE mod_precision
!!   * USE mod_input_grid
!!
!! USED BY
!!   * posini
!!   * trajec
!!   * 
!!
!! SOURCE
!!=========================================================================
MODULE mod_fz

  USE mod_precision
  USE mod_input_grid
  USE mod_reducmem

  IMPLICIT NONE

CONTAINS
  !!***
  !=========================================================================
  !!****f* mod_fz/fz()
  !! NAME
  !!   fz()
  !!
  !! FUNCTION
  !!   fz computes the ...
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
  !!   * gk : k position.
  !!
  !! TODO
  !!   * add comments
  !!
  !! USED BY
  !!   * posini
  !!   * trajec
  !!   * mod_quant
  !!
  !! SOURCE
  !!=======================================================================
  FUNCTION fz(gi, gj, gk, il)

    !-------------!
    ! Declaration !
    !-------------!

    !- arguments -!
    INTEGER(kind = iprec), INTENT(in), OPTIONAL :: il
    REAL(kind=rprec), INTENT(in) :: gi, gj
    REAL(kind=rprec), INTENT(in)  :: gk
    REAL(kind=rprec)              :: fz

    !- local variables -!
    INTEGER(kind=iprec) :: i1, i2, i1s, i2s
    INTEGER(kind=iprec) :: j1, j2, j1s, j2s
    INTEGER(kind=iprec) :: k1, k2, k1s, k2s
    REAL(kind=rprec)    :: a,b,c

    !-------------!
    ! Code begins !
    !-------------!
    !- Compute k1 and k2 -!
    c  = -gk - AINT(-gk, kind = rprec)
    k1 = INT(-gk, kind = iprec)


    !!NG: 26_02_2008 ===> Error if k1=kmt <====
    k2 = k1 + 1_iprec
    !!NG: 26_02_2008 Add this line to solve this problem
    IF (k2.GT.kmt) k2 = kmt

    !---------------------------------------------------!
    !                 Here fz is computed               !
    !---------------------------------------------------!
    IF (key_roms.OR.key_mars.OR.key_symphonie) THEN

      IF (PRESENT(il)) THEN

        !- compute i1 and i2 -!
        i1 = INT(gi+0.5_rprec, kind = iprec)
        i2 = i1 + 1_iprec
        a  = gi - ( REAL(i1, kind = rprec) - .5_rprec)
        IF (i2.GT.imt) i2 = i1

        !- compute j1 and j2 -!
        j1 = INT(gj+0.5_rprec, kind = iprec)
        j2 = j1 + 1_iprec
        b  = gj - ( REAL(j1, kind = rprec) - .5_rprec)
        IF (j2.GT.jmt) j2 = j1

        !- -!
        CALL sub_reducmem_shift_or_not_ind(i1,j1,k1,i1s,j1s,k1s)
        CALL sub_reducmem_shift_or_not_ind(i2,j2,k2,i2s,j2s,k2s)

        !- Compute fz (ROMS and Symphonie case) -!
        fz= &
             (1._rprec-a)*(1._rprec-b)*(-(1._rprec-c)*zz_ww(i1s,j1s,k1s,il) - c*zz_ww(i1s,j1s,k2s,il))+ &
             (1._rprec-a)*          b *(-(1._rprec-c)*zz_ww(i1s,j2s,k1s,il) - c*zz_ww(i1s,j2s,k2s,il))+ &
             a           *      (1.-b)*(-(1._rprec-c)*zz_ww(i2s,j1s,k1s,il) - c*zz_ww(i2s,j1s,k2s,il))+ &
             a           *          b *(-(1._rprec-c)*zz_ww(i2s,j2s,k1s,il) - c*zz_ww(i2s,j2s,k2s,il))
         
      ELSE
        WRITE(lun_error,*)' mod_fz: fz: there''s a bug in arguments...'
        STOP
      ENDIF

    ELSE
      i1 = INT(gi+0.5_rprec, kind = iprec)
      j1 = INT(gj+0.5_rprec, kind = iprec)
      CALL sub_reducmem_shift_or_not_ind(i1,j1,k1,i1s,j1s,k1s)
      CALL sub_reducmem_shift_or_not_ind(i1,j1,k2,i1s,j1s,k2s)

      !- Compute fz (default case = OPA) -!
      ! write (lun_standard, *) zz_ww(i1s,j1s,k1s,1), zz_ww(i1s,j1s,k2s,1), k1s, k2s
      fz= -(1._rprec - c) * zz_ww(i1s,j1s,k1s,1) - c * zz_ww(i1s,j1s,k2s,1)
   ENDIF
    !---------------------------------------------------!
    RETURN

  END FUNCTION fz
  !!***
END MODULE mod_fz
