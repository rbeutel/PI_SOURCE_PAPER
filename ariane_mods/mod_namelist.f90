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
!============================================================================
! #    #  ####  #####      #    #   ##   #    # ###### #      #  ####  #####
! ##  ## #    # #    #     ##   #  #  #  ##  ## #      #      # #        #
! # ## # #    # #    #     # #  # #    # # ## # #####  #      #  ####    #
! #    # #    # #    #     #  # # ###### #    # #      #      #      #   #
! #    # #    # #    #     #   ## #    # #    # #      #      # #    #   #
! #    #  ####  ##### #### #    # #    # #    # ###### ###### #  ####    #
!============================================================================
!!****h* ariane/mod_namelist
!! NAME
!!   mod_namelist (mod_namelist.f90 - Fortran90 module)
!!
!! USAGE
!!   Include 'USE mod_namelist' in the header of your Fortran 90 source 
!!   code.
!!   Then you'll have access to the subroutine:
!!     * sub_read_namelist()
!!
!! FUNCTION
!!   DECLARE global key variables, parameters and some data and
!!   READ their values in a namelist.
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
!!   1/ Key variables, parameters and some data are declared in global mode.
!!   2/ They are initialized from a namelist file (see sub_read_namelist).
!!
!! EXAMPLES
!!   * trajec.f90:  USE mod_namelist
!!                  CALL sub_read_namelist()
!!
!! NOTES
!!   ROBODoc header style.
!!
!! TODO
!!   * Delete some variables need only in the obsolete single precision 
!!     mode (ex.: key_doublelogexp).
!!   * The namelist file should contain only common Quantitative and 
!!     Qualitative parameters.
!!   * Create specific files and subroutines to read these files for 
!!     the Quantitative and Qualitative parameters 
!!     (namequant and namequal - sub_read_namequant and sub_read_namequal).
!!
!! PORTABILITY
!!         Machine-OS    - Fortran90/95 compiler
!!   * i686-pc-linux-gnu -         ifort
!!
!! SEE ALSO
!!   * namelist file
!!
!! USES
!!   * USE mod_precision
!!
!! USED BY
!!   * trajec.f90:  USE mod_namelist
!!                  CALL sub_read_namelist()
!!   * posini.f90:  USE mod_namelist
!!   * grid42.f90:  USE mod_namelist
!!   * speed40.f90:  USE mod_namelist
!!   * mod_fx.f90:  USE mod_namelist
!!   * mod_fy.f90:  USE mod_namelist
!!   * mod_reducmem.f90:  USE mod_namelist
!!   * mod_quant.f90:  USE mod_namelist
!!   * mod_tmask.f90:  USE mod_namelist
!!   * mod_tracer.f90:  USE mod_namelist
!!   * mod_txt.f90:  USE mod_namelist
!!   * mod_zinter.f90:  USE mod_namelist
!!
!! SOURCE
!!=========================================================================
MODULE mod_namelist

  !-----------------!
  ! USE ASSOCIATION !
  !-----------------!
  USE mod_precision 
  USE mod_cst
  USE mod_lun

  !-------------!
  ! DECLARATION !
  !-------------!
  IMPLICIT NONE
  !-----------------------------------------------------------------------
  ! Default values for integer, real character and logical namelist variables
!!$  INTEGER(kind=iprec), PARAMETER :: i_defp = 1
!!$  INTEGER(kind=iprec), PARAMETER :: i_defn = -1
!!$  INTEGER(kind=iprec), PARAMETER :: i_high = 9999999
!!$  REAL(kind=rprec)   , PARAMETER :: r_def  = 0._rprec
!!$  CHARACTER(len = 4) , PARAMETER :: c_def  = 'NONE'
!!$  LOGICAL            , PARAMETER :: l_def  = .FALSE.

  REAL(kind = rprec)             :: zsigma_buffer

  !--------------------!
  !- Public variables -!
  !--------------------!
  INTEGER(kind = iprec) :: &!
       nstat             = i_defn ! dimension of the statistical arrays (key_alltracers)

  !-------------!
  ! ARIANE item !
  !-------------!
  LOGICAL               ::             & !
       key_roms             = l_def  , & !# Read ROMS data format (default OPA format)
       key_mars             = l_def  , & !# Read MARS3D data format
       key_write_transport  = l_def  , & !# Write transport in a Netcdf File !
       key_symphonie        = l_def  , & !# Read SYMPHONIE data format (default OPA format)
       key_B2C_grid         = l_def  , & !# Read data on B grid and interpol them on C grid
       key_sequential       = l_def  , & !# read input data sequentially.
       key_alltracers       = l_def  , & !# activation de T, S, RAU
       key_ascii_outputs    = .TRUE. , & !# to activate outputs in ascii files like before.
       key_iU_jV_kW         = l_def  , & !# Qualitative = add i on grid U, j on V and k on W
       key_2D               = l_def  , & !# Ariaen 2D (Warning: Work in Progress january 2019!)
       key_read_age         = l_def      !# keep or not the age of the particule when
  !# using bin or subbin options.
  CHARACTER(len = 12)   ::             & !
       mode                 = c_def  , & !# qualitative or quantitative mode
       forback              = c_def  , & !# forward or backward
       bin                  = c_def  , & !# nobin, bin or subbin
       init_final           = c_def      !# if bin /= nobin => init or final
  INTEGER(kind = iprec) ::             & !
       nmax                 = i_defn , & !# Maximun number of particles
       it_ind               = i_defp     !# temporal dimension when interpolation (seq)
  REAL(kind = rprec)    ::             & !
       tunit                = r_def      !# Convenient unit time
  INTEGER(kind = iprec) ::             & !
       ntfic                = i_defn     !# Sampling time (in number of tnunit)
  REAL(kind = rprec)    ::             & !
       tcyc                 = r_def      !# NEW : unit to print ages
  LOGICAL               ::             & !
       key_approximatesigma = l_def      !# linear interpolation of the density
  LOGICAL               ::             & !
       key_computesigma     = l_def      !# Compute the density from the
  !# salinity and temperature
  REAL(kind=rprec)      ::             & !
       zsigma               = 1000._rprec !# density used in mod_rhostp.f90
  LOGICAL               ::             & !
       memory_log           = .TRUE.     ! Active memory log
  LOGICAL               ::             & !
       output_netcdf_large_file = l_def  ! to save Netcdf Large file
  LOGICAL               ::             & !
       small_virtual_memory     = l_def  ! to reduce virtual memory but can increase the writting time (DON'T USE IT - IT IS NOT OPERATIONAL)
  LOGICAL               ::             & !
       roms_global_attribute    = .TRUE. ! read hc and so on as global attributes !

  NAMELIST/ARIANE/           & !
       key_roms            , & !
       roms_global_attribute , & !
       key_mars            , & !
       key_write_transport , & !
       key_symphonie       , & !
       key_B2C_grid        , & !
       key_sequential      , & !
       key_alltracers      , & ! 
       key_ascii_outputs   , & !
       key_iU_jV_kW        , & !
       key_2D              , & !   Work in Progress
       mode                , & ! 
       forback             , & ! 
       bin                 , & !
       init_final          , & !
       key_read_age        , & !
       nmax                , & ! 
       tunit               , & ! 
       ntfic               , & ! 
       tcyc                , & !
       key_approximatesigma, & !
       key_computesigma    , & !
       zsigma              , & !
       memory_log          , & !
       output_netcdf_large_file, & !
       small_virtual_memory    !

  !-----------------!
  ! SEQUENTIAL item !
  !-----------------!
  INTEGER(kind = iprec) :: & !
       maxcycles = i_defp

  LOGICAL               :: & !
       key_interp_temporal  = l_def  !# Temporal interpolation IN TEST

  NAMELIST/SEQUENTIAL/     & !
       key_interp_temporal , & !
       maxcycles

  !-------------------!
  ! QUANTITATIVE item !
  !-------------------!
  LOGICAL             ::             & !
       key_2dquant         = l_def , & !# qunatitative experiment in 2D (same immersion)
       key_eco             = l_def , & !# supprime des calculs lourds
       key_reducmem        = l_def , & !# Reduce memory in reading only the selected region
       key_unitm3          = l_def , & !# m3/s au lieu de Sv dans stats.qt
       key_nointerpolstats = l_def     !# statsistiques non interpolees
  REAL(kind = rprec)    ::           & !
       max_transport       = -rOne     !# Maximum transport value (m^3/s)
  INTEGER(kind = iprec) ::           & !
       lmin                = i_defn, & !
       lmax                = i_defn

  NAMELIST/QUANTITATIVE/             & !
       key_2dquant                 , & !
       key_eco                     , & !
       key_reducmem                , & !
       key_unitm3                  , & !
       key_nointerpolstats         , & !
       max_transport               , & !
       lmin                        , & !
       lmax                            !

  !------------------!
  ! QUALITATIVE item !
  !------------------!
  REAL(kind = rprec)    ::   & !
       delta_t    = r_def      !# Define a convenient unit time (in seconds)
  INTEGER(kind = iprec) ::   & !
       frequency  = i_defn , & !# Frequency of output
       nb_output  = i_high     !# Maximum number of output trajectories
  !!   - used in qualitative mode
  !!   - but must have very high value in quant. mode
  LOGICAL               ::   & !
       key_region  = l_def , & !# To select a geographic region (reduce memory)
       mask        = l_def     !# Suppress land point if it is true

  NAMELIST/QUALITATIVE/      & !
       delta_t             , & !
       frequency           , & !
       nb_output           , & !
       key_region          , & !
       mask                    !

  !---------------!
  ! OPAPARAM item !
  !---------------!
  INTEGER(kind=iprec) ::           & !
       imt              = i_defn , & ! Dimension in x (longitude).
       jmt              = i_defn , & ! Dimension in y (latitude).
       kmt              = i_defn , & ! Dimension in z (dept).
       lmt              = i_defn     ! Time dimension (time).

  REAL(kind = rprec) :: &
       epr_coef         = r_def      ! Coefficient to tranfert E-P-R in m/s 
  ! by default E-P-R is in kg/m2/s
  ! a coefficeint of 1.e-3 transfom the
  ! kg/m2/s in m/s

  CHARACTER(len=1)    :: & !
       pivot            = 'T' ! Pivot point "T" or "F"
  !                           ! T-point pivot: ORCA 4, 2, 025
  !                           ! F-point pivot: ORCA 05

  CHARACTER(len=5)    :: & ! 'E-P-R' or 'zero' or ''
       w_surf_option    = ''

  LOGICAL             ::           & !
       key_periodic     = l_def  , & !# domaine periodique en i=imt
       key_jfold        = l_def  , & !# repliement de la grille en j=jmt
       key_computew     = l_def  , & !# recalcul du transport vertical
       key_partialsteps = l_def  , & !# partial steps
       key_vvl          = l_def  , & !# stretch grid with ssh.
       key_sigma        = l_def      !

  !!NG: modified 27_02_2008
  !!NG: WARNING - I think there is a bug here, because key_sigma and zsigma should
  !!NG: WARNING - be also available for ROMS... Works only if zsigma = 0.
  !!NG: WARNING - key_sigma should be renamed key_computesigma     !
  !!NG: WARNING - and key_approximatesigma -> key_interpolatesigma !
  !!NG: modified 27_02_2008

  NAMELIST/OPAPARAM/         & !
       imt                 , & !
       jmt                 , & !
       kmt                 , & !
       lmt                 , & !
       key_computew        , & !
       key_partialsteps    , & !
       key_vvl             , & !
       key_jfold           , & !
       pivot               , & !
       key_periodic        , & !
       key_sigma           , & !
       zsigma              , & !
       w_surf_option       , & !
       epr_coef          

  !-----------------!
  ! ROMSPARAMS item !
  !-----------------!
  INTEGER(kind=iprec) ::  & !
       xi_rho  = i_defn , & ! Dimension in x (longitude).
       eta_rho = i_defn , & ! Dimension in y (latitude).
       s_w     = i_defn , & ! Dimension in z (dept).
       time    = i_defn     ! Time dimension (time).

  NAMELIST/ROMSPARAM/     & !
       xi_rho           , & !
       eta_rho          , & !
       s_w              , & !
       time                 !

  !-------------------!
  ! MARSPARAMS item !
  !-------------------!
  INTEGER(kind=iprec) ::  & !
       x_t = i_defn , & ! Dimension in x (longitude).
       y_t = i_defn , & ! Dimension in y (latitude).
       sigma_t   = i_defn !, & ! Dimension in z (dept).
  !       time    = i_defn     ! Time dimension (time).

  NAMELIST/MARSPARAM/   & !
       x_t        , & !
       y_t        , & !
       sigma_t          , & !
       time               !

  !----------------------!
  ! SYMPHONIEPARAMS item !
  !----------------------!
  INTEGER(kind=iprec) ::  & !
       x_dim = i_defn , & ! Dimension in x (longitude).
       y_dim = i_defn , & ! Dimension in y (latitude).
       z_dim = i_defn     ! Dimension in z (dept).
  !!       time  = i_defn     ! Time dimension (time). (already been assigned)

  NAMELIST/SYMPHONIEPARAM/ & !
       x_dim             , & !
       y_dim             , & !
       z_dim             , & !
       time                  !

  !----------!
  ! B2C item !
  !----------!
  LOGICAL               ::          & !
       key_B2C_save_data = l_def  , & !
       key_read_w        = l_def  , & ! # Vertical velocity is available or not
       periodic_lon      = l_def  , & ! # periodicity in longitude
       periodic_lat      = l_def  , &
       key_add_bottom    = l_def      ! add a level at the bottom : ex. for MOM

  INTEGER(kind = iprec) ::    &
       nb_dim_lon   = i_defn, &
       nb_dim_lat   = i_defn, &
       nb_dim_depth = i_defn, &
       nb_dim_time  = i_defn  

  CHARACTER(len = 5)   :: B2C_grid_Z_or_Sigma

  !!  REAL(kind = rprec)    ::

  NAMELIST/B2C/             & !
       B2C_grid_Z_or_Sigma, & !
       nb_dim_lon         , & !
       nb_dim_lat         , & !
       nb_dim_depth       , & !
       nb_dim_time        , & !
       key_add_bottom     , & !
       key_partialsteps   , & !
       key_B2C_save_data  , & !
       key_read_w         , & !
       periodic_lon       , & !
       periodic_lat   

  !------------------------------------------------------------------------
  INTEGER(kind = iprec) :: &
       ind0_zo = i_defn, indn_zo = i_defn, maxsize_zo = i_defn, &
       ind0_me = i_defn, indn_me = i_defn, maxsize_me = i_defn, &
       ind0_ve = i_defn, indn_ve = i_defn, maxsize_ve = i_defn, &
       ind0_ssh = i_defn, indn_ssh = i_defn, maxsize_ssh = i_defn, &
       ind0_te = i_defn, indn_te = i_defn, maxsize_te = i_defn, &
       ind0_sa = i_defn, indn_sa = i_defn, maxsize_sa = i_defn, &
       ind0_de = i_defn, indn_de = i_defn, maxsize_de = i_defn, &
       ind0_ze = i_defn, indn_ze = i_defn, maxsize_ze = i_defn, &
       ind0_sse = i_defn, indn_sse = i_defn, maxsize_sse = i_defn, &
       ind0_zt = i_defn, indn_zt = i_defn, maxsize_zt = i_defn, &
       ind0_ep = i_defn, indn_ep = i_defn, maxsize_ep = i_defn

  CHARACTER(len = 64) :: &
       c_prefix_zo    = c_def, c_suffix_zo    = c_def, nc_var_zo = c_def, &
       nc_var_eivu    = c_def, nc_att_mask_zo = c_def                   , &
       c_prefix_me    = c_def, c_suffix_me    = c_def, nc_var_me = c_def, &
       nc_var_eivv    = c_def, nc_att_mask_me = c_def                   , &
       c_prefix_ve    = c_def, c_suffix_ve    = c_def, nc_var_ve = c_def, &
       nc_var_eivw    = c_def, nc_att_mask_ve = c_def                   , &
       c_prefix_ssh    = c_def, c_suffix_ssh    = c_def, nc_var_ssh = c_def, &
       nc_att_mask_ssh = c_def                                           , &
       c_prefix_te    = c_def, c_suffix_te    = c_def, nc_var_te = c_def, &
       nc_att_mask_te = c_def                                           , &
       c_prefix_sa    = c_def, c_suffix_sa    = c_def, nc_var_sa = c_def, &
       nc_att_mask_sa = c_def                                           , &
       c_prefix_de    = c_def, c_suffix_de    = c_def, nc_var_de = c_def, &
       nc_att_mask_de = c_def                                           , &
       c_prefix_ze    = c_def, c_suffix_ze    = c_def, nc_var_ze = c_def, &
       nc_att_mask_ze = c_def                                           , &
       c_prefix_sse    = c_def, c_suffix_sse    = c_def, nc_var_sse = c_def, &
       nc_att_mask_sse = c_def                                          , &
       c_prefix_zt    = c_def, c_suffix_zt    = c_def, nc_var_zt = c_def, &
       nc_att_mask_zt = c_def                                          , &
       c_prefix_ep    = c_def, c_suffix_ep    = c_def, nc_var_ep = c_def, &
       nc_att_mask_ep = c_def

  CHARACTER(len = 128) :: &
       c_dir_zo = c_def, c_dir_me = c_def, c_dir_ve = c_def, &
       c_dir_ssh = c_def, &
       c_dir_te = c_def, c_dir_sa = c_def, c_dir_de = c_def, &
       c_dir_ze = c_def, c_dir_sse = c_def, c_dir_ep = c_def, &
       c_dir_zt = c_def

  NAMELIST/ZONALCRT/ &
       c_dir_zo, c_prefix_zo, ind0_zo, indn_zo, maxsize_zo, c_suffix_zo, &
       nc_var_zo, nc_var_eivu, nc_att_mask_zo
  NAMELIST/MERIDCRT/ &
       c_dir_me, c_prefix_me, ind0_me, indn_me, maxsize_me, c_suffix_me, &
       nc_var_me, nc_var_eivv, nc_att_mask_me
  NAMELIST/VERTICRT/ &
       c_dir_ve, c_prefix_ve, ind0_ve, indn_ve, maxsize_ve, c_suffix_ve, &
       nc_var_ve, nc_var_eivw, nc_att_mask_ve
  NAMELIST/SSH/ &
       c_dir_ssh, c_prefix_ssh, ind0_ssh, indn_ssh, maxsize_ssh, c_suffix_ssh, &
       nc_var_ssh, nc_att_mask_ssh
  NAMELIST/E_P_R/ &
       c_dir_ep, c_prefix_ep, ind0_ep, indn_ep, maxsize_ep, c_suffix_ep, &
       nc_var_ep, nc_att_mask_ep
  NAMELIST/TEMPERAT/ &
       c_dir_te, c_prefix_te, ind0_te, indn_te, maxsize_te, c_suffix_te, &
       nc_var_te, nc_att_mask_te
  NAMELIST/SALINITY/ &
       c_dir_sa, c_prefix_sa, ind0_sa, indn_sa, maxsize_sa, c_suffix_sa, &
       nc_var_sa, nc_att_mask_sa
  NAMELIST/DENSITY/ &
       c_dir_de, c_prefix_de, ind0_de, indn_de, maxsize_de, c_suffix_de, &
       nc_var_de, nc_att_mask_de
  NAMELIST/ZETA/ &
       c_dir_ze, c_prefix_ze, ind0_ze, indn_ze, maxsize_ze, c_suffix_ze, &
       nc_var_ze, nc_att_mask_ze
  NAMELIST/SSE/ &
       c_dir_sse, c_prefix_sse, ind0_sse, indn_sse, maxsize_sse, c_suffix_sse, &
       nc_var_sse, nc_att_mask_sse
  NAMELIST/ZZ_TT_SIGMA/&
       c_dir_zt, c_prefix_zt, ind0_zt, indn_zt, maxsize_zt, c_suffix_zt, &
       nc_var_zt, nc_att_mask_zt

  !------------------------------------------------------------------------
  CHARACTER(len = 128) :: dir_glbatt = c_def, fn_glbatt = c_def
  CHARACTER(len = 64) :: &
       nc_glbatt_hc   = c_def , &
       nc_glbatt_sc_w = c_def , &
       nc_glbatt_Cs_w = c_def
  NAMELIST/GLOBALATT/         &
       dir_glbatt, fn_glbatt, &
       nc_glbatt_hc         , &
       nc_glbatt_sc_w       , &
       nc_glbatt_Cs_w
  !------------------------------------------------------------------------
  CHARACTER(len = 128) :: dir_mesh = c_def, fn_mesh = c_def
  CHARACTER(len = 64) :: &
       nc_var_xx_tt = c_def, nc_var_xx_uu = c_def, &
       nc_var_xx_vv = c_def, nc_var_xx_ff = c_def, &
       nc_var_yy_tt = c_def, nc_var_yy_uu = c_def, &
       nc_var_yy_vv = c_def, nc_var_yy_ff = c_def, &
       nc_var_zz_tt = c_def, nc_var_zz_ww = c_def, &
       nc_var_e1f   = c_def, nc_var_e2f   = c_def, &
       nc_var_e2u   = c_def, nc_var_e1v   = c_def, &
       nc_var_e1t   = c_def, nc_var_e2t   = c_def, &
       nc_var_e3t   = c_def, nc_var_e3f   = c_def, &
       nc_var_tmask = c_def                      , &
       mesh_type    = c_def, nc_var_e3t2d = c_def, &
       nc_var_e3tz  = c_def, nc_var_mbathy= c_def   !! Nemov3 and vvl
  REAL(kind = rprec) :: nc_mask_val = r_def

  NAMELIST/MESH/ &
       dir_mesh   , fn_mesh      , &
       nc_var_xx_tt, nc_var_xx_uu, nc_var_xx_vv, nc_var_xx_ff, &
       nc_var_yy_tt, nc_var_yy_uu, nc_var_yy_vv, nc_var_yy_ff, &
       nc_var_zz_tt, nc_var_zz_ww, &
       nc_var_e2u  , nc_var_e1v  , &
       nc_var_e1t  , nc_var_e2t  , nc_var_e3t, &
       nc_var_tmask, nc_mask_val , &
       mesh_type, nc_var_e3t2d, nc_var_e3tz, nc_var_mbathy

  CHARACTER(len = 128) :: dir_grd_roms = c_def, fn_grd_roms = c_def
  CHARACTER(len = 64) :: &
       nc_var_lon_rho_roms = c_def, nc_var_lon_u_roms    = c_def, &
       nc_var_lat_rho_roms = c_def, nc_var_lat_v_roms    = c_def, &
       nc_var_pm_roms      = c_def, nc_var_pn_roms       = c_def, &
       nc_var_h_roms       = c_def, nc_var_mask_rho_roms = c_def

  NAMELIST/GRDROMS/ &
       dir_grd_roms       , fn_grd_roms      , &
       nc_var_lon_rho_roms, nc_var_lon_u_roms, &
       nc_var_lat_rho_roms, nc_var_lat_v_roms, &
       nc_var_pm_roms     , nc_var_pn_roms   , &
       nc_var_h_roms      , nc_var_mask_rho_roms

  CHARACTER(len = 128) :: dir_grd_mars = c_def, fn_grd_mars = c_def
  CHARACTER(len = 64) :: &
       nc_var_lon_t_mars   = c_def, nc_var_lon_u_mars  = c_def, &
       nc_var_lat_t_mars   = c_def, nc_var_lat_v_mars  = c_def, &
       nc_var_e1t_mars     = c_def, nc_var_e2t_mars    = c_def, &
       nc_var_e2u_mars     = c_def, nc_var_e1v_mars    = c_def, &
       nc_var_hc_mars      = c_def,                             &
       nc_var_Cs_w_mars    = c_def, nc_var_sc_w_mars   = c_def, &
       nc_var_bathy_t_mars = c_def, nc_var_mask_t_mars = c_def

  NAMELIST/GRDMARS/ &
       dir_grd_mars       , fn_grd_mars      , &
       nc_var_lon_t_mars  , nc_var_lon_u_mars, &
       nc_var_lat_t_mars  , nc_var_lat_v_mars, &
       nc_var_e1t_mars    , nc_var_e2t_mars  , &
       nc_var_e2u_mars    , nc_var_e1v_mars  , &
       nc_var_hc_mars     ,                    &
       nc_var_Cs_w_mars   , nc_var_sc_w_mars , &
       nc_var_bathy_t_mars, nc_var_mask_t_mars

  REAL(kind=rprec)     :: cst_scale_factor = r_def  !! in meter
  CHARACTER(len = 128) :: dir_grd_symp = c_def, fn_grd_symp = c_def
  CHARACTER(len = 64 ) ::                                         &
       nc_var_lon_t_symp   = c_def, nc_var_lon_u_symp    = c_def, &
       nc_var_lat_t_symp   = c_def, nc_var_lat_v_symp    = c_def, &
       nc_var_depth_t_symp = c_def 

  NAMELIST/GRDSYMPHONIE/                       &
       dir_grd_symp       , fn_grd_symp      , &
       cst_scale_factor                      , &
       nc_var_lon_t_symp  , nc_var_lon_u_symp, &
       nc_var_lat_t_symp  , nc_var_lat_v_symp, &
       nc_var_depth_t_symp   

  !-----------------!
  !- B2CGRIDZ item -!
  !-----------------!
  CHARACTER(len = 128) ::          &
       dir_B2C_grid       = c_def, &
       file_name_B2C_grid = c_def

  NAMELIST/B2CGRIDZ/       &
       dir_B2C_grid      , &
       file_name_B2C_grid, &
       nc_var_xx_tt, nc_var_xx_uu, &
       nc_var_yy_tt, nc_var_yy_vv, &
       nc_var_zz_ww      , &
       nc_var_e1f        , & !# This is for interpolation from B to C grid
       nc_var_e2f        , & !# This is for interpolation from B to C grid
       nc_var_e3f        , & !# This is for interpolation from B to C grid
       nc_var_e2u        , &
       nc_var_e1v        , &
       nc_var_e1t  , nc_var_e2t  , nc_var_e3t, &
       nc_var_tmask      , &
       nc_mask_val

CONTAINS
  !!***
  !=========================================================================
  !!****f* mod_namelist/sub_read_namelist()
  !! NAME
  !!   sub_read_namelist()
  !!
  !! FUNCTION
  !!   READ the values of global key variables, parameters and some data
  !!   from a namelist file.
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
  !!   No arguments.
  !!
  !! TODO
  !!   Add the possibility to specify n optionnal Logical Unit number 
  !!   (argument).
  !!
  !! USED BY
  !!   * trajec.f90: CALL sub_read_namelist()
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_read_namelist(l_info)

    LOGICAL, OPTIONAL, INTENT(IN) :: l_info

    LOGICAL :: l_write

    IF (PRESENT(l_info)) THEN
      IF (l_info) THEN
        l_write=.TRUE.
      ELSE
        l_write=.FALSE.
      ENDIF
    ELSE
      l_write=.FALSE.
    ENDIF

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)'============'
    WRITE(lun_standard,*)'= NAMELIST ='
    WRITE(lun_standard,*)'============'

    OPEN(unit=lun_nml, File='namelist', ACTION='READ')
    WRITE(lun_standard,*)' --- Successful Opening ---'

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' - Reading ARIANE item:'
    REWIND(unit=lun_nml)
    READ(unit=lun_nml, nml=ARIANE)
    IF (l_write) THEN
      WRITE(lun_standard,*)'   - key_2D (WIP)        = ', key_2D
      WRITE(lun_standard,*)'   - key_roms            = ', key_roms
      IF (key_roms) THEN
        WRITE(lun_standard,*)'   - roms_global_attribute= ', roms_global_attribute
      ENDIF
      WRITE(lun_standard,*)'   - key_mars            = ', key_mars
      WRITE(lun_standard,*)'   - key_write_transport = ', key_write_transport
      WRITE(lun_standard,*)'   - key_symphonie       = ', key_symphonie
      WRITE(lun_standard,*)'   - key_B2C_grid        = ', key_B2C_grid
      WRITE(lun_standard,*)'   - key_sequential      = ', key_sequential
      WRITE(lun_standard,*)'   - key_alltracers      = ', key_alltracers
      WRITE(lun_standard,*)'   - key_ascii_outputs   = ', key_ascii_outputs
      WRITE(lun_standard,*)'   - key_iU_jV_kW        = ', key_iU_jV_kW 
      WRITE(lun_standard,*)'   - mode                = ', TRIM(mode)
      WRITE(lun_standard,*)'   - forback             = ', TRIM(forback)
      WRITE(lun_standard,*)'   - bin                 = ', TRIM(bin)
      IF (TRIM(bin) /= 'nobin') THEN
        WRITE(lun_standard,*)'   - init_final          = ', TRIM(init_final)
        WRITE(lun_standard,*)'   - key_read_age        = ', key_read_age
      ENDIF
      WRITE(lun_standard,*)'   - nmax                = ', nmax
      WRITE(lun_standard,*)'   - tunit               = ', tunit
      WRITE(lun_standard,*)'   - ntfic               = ', ntfic
      IF ( tcyc == r_def) THEN
        WRITE(lun_standard,*)'   - tcyc                = "value will be computed later"'
      ELSE
        WRITE(lun_standard,*)'   - tcyc                =', tcyc
      ENDIF

      !! NG: Bug fixed 15 10 2008 between key_computesigma and key_approximatesigma
      IF (key_roms.OR.key_mars) THEN
        key_computesigma=.TRUE.
      ENDIF
      IF (key_computesigma) THEN
        key_approximatesigma=.FALSE.
      ELSE
        key_approximatesigma=.TRUE.
      ENDIF

      WRITE(lun_standard,*)'   - key_approximatesigma (OBSOLETE and FORCED)=', key_approximatesigma
      WRITE(lun_standard,*)'   - key_computesigma    =', key_computesigma
      IF (key_roms.OR.key_mars) THEN
        key_computesigma=.TRUE.
      ENDIF
      IF (key_computesigma) THEN
        WRITE(lun_standard,*)'   - zsigma              =', zsigma
      ENDIF

      WRITE(lun_standard,*)'   - memory_log          =', memory_log
      WRITE(lun_standard,*)'   - output_netcdf_large_file=', &
           output_netcdf_large_file
    ENDIF

    small_virtual_memory     = l_def !! DON'T USE IT - NOT OPERATIONAL

    !-------------------!
    !- Sequential Mode -!
    !-------------------!
    IF (key_sequential) THEN
      CALL sub_read_namelist_sequential(l_write)
    ENDIF

    !-------------!
    !- Test mode -!
    !-------------!
    IF  ((TRIM(mode) == 'QUANTITATIVE').OR.(TRIM(mode) == 'quantitative')) THEN
      mode='quantitative'
      IF (key_iU_jV_kW) THEN
        key_iU_jV_kW = .FALSE.
        WRITE(lun_standard,*)'   - change key_iU_jV_kW to .False. '
      END IF
    ELSEIF ((TRIM(mode) == 'QUALITATIVE').OR.(TRIM(mode) == 'qualitative')) THEN
      mode='qualitative'
    ELSE
      WRITE(lun_error,*) 'string must be: quantitative OR qualitative'
      STOP
    ENDIF

    !-----------!
    !- forback -!
    !-----------!
    IF ((TRIM(forback) == 'FORWARD').OR.(TRIM(forback) == 'forward')) THEN
      forback='forward'
    ELSEIF ((TRIM(forback) == 'BACKWARD').OR.(TRIM(forback) == 'backward')) THEN
      forback='backward'
    ELSE
      WRITE(lun_error,*)'string must be: forward OR backward'
      STOP
    ENDIF

    !-------!
    !- bin -!
    !-------!
    IF ((TRIM(bin) == 'BIN').OR.(TRIM(bin) == 'bin')) THEN
      bin='bin'
    ELSEIF ((TRIM(bin) == 'SUBBIN').OR.(TRIM(bin) == 'subbin')) THEN
      bin='subbin'
    ELSEIF ((TRIM(bin) == 'NOBIN').OR.(TRIM(bin) == 'nobin')) THEN
      bin='nobin'
    ELSE
      WRITE(lun_error,*)'string must be: nobin OR bin OR subbin'
      STOP
    ENDIF

    !--------------!
    !- init_final -!
    !--------------!
    IF (TRIM(bin) /= 'nobin') THEN
      IF ((TRIM(init_final) == 'INIT').OR.(TRIM(init_final) == 'init')) THEN
        init_final = 'init'
      ELSEIF  ((TRIM(init_final) == 'FINAL').OR.(TRIM(init_final) == 'final')) THEN
        init_final = 'final'
      ELSE
        WRITE(lun_error,*)'IF bin is different than ''nobin'', &
             &init_final have to be set.'
        WRITE(lun_error,*)' and init_final string must be: init or final.'
        STOP
      ENDIF
    ENDIF

    !--------------------!
    !- Key_computesigma -!
    !--------------------!
    IF (key_computesigma) THEN
      IF (zsigma == r_def) THEN
        WRITE(lun_standard,*)''
        WRITE(lun_standard,*) &
             '==========='
        WRITE(lun_standard,*) &
             '= WARNING = key_computesigma=.TRUE. and zsigma is ', r_def
        WRITE(lun_standard,*) &
             '==========='
      ENDIF
      zsigma_buffer = zsigma
    ENDIF

    !--------------------------------------------------!
    !- Initialize allocatable array dimensions (stat) -!
    !--------------------------------------------------!
    IF (TRIM(mode) == 'qualitative') THEN

      IF (key_alltracers) THEN
        nstat = 7
      ELSE
        nstat = 4
      END IF

      IF (key_iU_jV_kW) THEN
        nstat = nstat + 3
      ENDIF

    ELSEIF (TRIM(mode) == 'quantitative') THEN

      IF (key_alltracers) THEN
        nstat = 7
      ELSE
        nstat = 4
      END IF

    ELSE

      WRITE(lun_error,*)' Error:mod_namelist:  STOP (nstat)'
      STOP

    ENDIF

    !-------------------------------!
    !- Read OPA or ROMS parameters -!
    !-------------------------------!
    IF (key_roms) THEN
      CALL sub_read_namelist_roms(l_write)
    ELSEIF (key_mars) THEN
      CALL sub_read_namelist_mars(l_write)
    ELSEIF (key_symphonie) THEN
      CALL sub_read_namelist_symphonie(l_write)
    ELSEIF (key_B2C_grid) THEN
      CALL sub_read_namelist_B2C(l_write)
    ELSE
      CALL sub_read_namelist_opa(l_write)
    ENDIF

    !-----------------------------------------------!
    !- Read QUANTITATIVE or QUALITATIVE parameters -!
    !-----------------------------------------------!
    IF (TRIM(mode) == 'quantitative') THEN
      CALL sub_read_namelist_quantitative(l_write)
    ELSEIF (TRIM(mode) == 'qualitative') THEN
      CALL sub_read_namelist_qualitative(l_write)
    ELSE
      WRITE(lun_error,*)'mod_namelist.f90: sub_read_namelist:'
      WRITE(lun_error,*)'Problem with index mode in ARIANE item'
      WRITE(lun_error,*)'-- WE STOP --'
      STOP
    ENDIF

    !-----------------------------------!
    !- Print tcyc if computed from lmt -!
    !-----------------------------------!
    IF (tcyc == r_def) THEN
      tcyc = tunit * REAL(ntfic, kind = rprec) * REAL(lmt, kind = rprec)
      IF (l_write) THEN
        WRITE(lun_standard,*)''
        WRITE(lun_standard,*)' - More:'
        WRITE(lun_standard,*)'     - tcyc                =', tcyc
      ENDIF
    ENDIF

    !---------------------------!
    !- Close the namelist file -!
    !---------------------------!
    CLOSE(unit=lun_nml)

    !- Comments -!
    IF (l_write) THEN
      IF (TRIM(mode) == 'qualitative') THEN 
        WRITE(lun_standard,*)''
        WRITE(lun_standard,*)'                          -'
        WRITE(lun_standard,*)'                         ---'
        WRITE(lun_standard,*)'                       -------'
        WRITE(lun_standard,*)'                     -----------'
        WRITE(lun_standard,*)'                -<>- QUALITATIVE -<>-'
        WRITE(lun_standard,*)'                     -----------'
        WRITE(lun_standard,*)'                       -------'
        WRITE(lun_standard,*)'                         ---'
        WRITE(lun_standard,*)'                          -'
      ELSEIF (TRIM(mode) == 'quantitative') THEN
        WRITE(lun_standard,*)''
        WRITE(lun_standard,*)'                          -'
        WRITE(lun_standard,*)'                         ---'
        WRITE(lun_standard,*)'                       -------'
        WRITE(lun_standard,*)'                     -----------'
        WRITE(lun_standard,*)'                -<>- QUANTITATIVE -<>-'
        WRITE(lun_standard,*)'                     -----------'
        WRITE(lun_standard,*)'                       -------'
        WRITE(lun_standard,*)'                         ---'
        WRITE(lun_standard,*)'                          -'
      ELSE
        STOP ! just for fun because not necessary
      ENDIF
    ENDIF

    key_sigma = key_computesigma

  END SUBROUTINE sub_read_namelist

  !!***
  !=========================================================================
  !!****f* mod_namelist/sub_read_namelist_quantitative()
  !! NAME
  !!   sub_read_namelist_quantitative()
  !!
  !! FUNCTION
  !!   READ  a namelist file.
  !!
  !! AUTHOR
  !!   * Origin  : Nicolas Grima
  !! 
  !! CREATION DATE
  !!   * January 2006
  !!
  !! HISTORY
  !!   Date (dd/mm/yyyy/) - Modification(s)
  !!
  !! ARGUMENTS
  !!
  !! TODO
  !!
  !! USED BY
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_read_namelist_sequential(l_write)

    LOGICAL, INTENT(IN) :: l_write

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' - Reading SEQUENTIAL item:'
    REWIND(unit=lun_nml)
    READ(unit=lun_nml, nml=SEQUENTIAL)

    IF (l_write) THEN
      WRITE(lun_standard,*)'   - key_interp_temporal =', key_interp_temporal
      WRITE(lun_standard,*)'   - maxcycles           =', maxcycles
    ENDIF

    !--------------------------!
    !- Temporal interpolation -!
    !--------------------------!
    IF (key_interp_temporal) THEN
      it_ind = 2
    ELSE
      it_ind = 1
    ENDIF

  END SUBROUTINE sub_read_namelist_sequential

  !!***
  !=========================================================================
  !!****f* mod_namelist/sub_read_namelist_quantitative()
  !! NAME
  !!   sub_read_namelist_quantitative()
  !!
  !! FUNCTION
  !!   READ  a namelist file.
  !!
  !! AUTHOR
  !!   * Origin  : Nicolas Grima
  !! 
  !! CREATION DATE
  !!   * January 2006
  !!
  !! HISTORY
  !!   Date (dd/mm/yyyy/) - Modification(s)
  !!
  !! ARGUMENTS
  !!
  !! TODO
  !!
  !! USED BY
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_read_namelist_quantitative(l_write)

    LOGICAL, INTENT(IN) :: l_write

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' - Reading QUANTITATIVE item:'
    REWIND(unit=lun_nml)
    READ(unit=lun_nml, nml=QUANTITATIVE)

    IF (lmin == i_defn) THEN 
      lmin = 1
    ELSEIF ( (lmin <= 0).OR.(lmin > lmt) ) THEN
      WRITE(lun_standard,*)''
      WRITE(lun_standard,*)' Namelist: ERROR - lmin value is invalid!'
      STOP
    ENDIF

    IF (lmax == i_defn) THEN
      lmax = lmt
    ELSEIF ( (lmax < lmin).OR.(lmax > lmt) ) THEN
      WRITE(lun_standard,*)''
      WRITE(lun_standard,*)' Namelist: ERROR - lmax value is invalid!'
      STOP
    ENDIF

    IF (key_sequential) THEN
      key_eco = .TRUE.
      IF (l_write) THEN
        WRITE(lun_standard,*)''
        WRITE(lun_standard,*)' Namelist: INFO - In sequential + quantitative modes &
             & KEY_ECO is forced to .TRUE..'
      ENDIF
    ENDIF

    IF (l_write) THEN
      WRITE(lun_standard,*)'   - key_2dquant         =', key_2dquant
      WRITE(lun_standard,*)'   - key_eco             =', key_eco
      WRITE(lun_standard,*)'   - key_reducmem        =', key_reducmem
      WRITE(lun_standard,*)'   - key_unitm3          =', key_unitm3
      WRITE(lun_standard,*)'   - key_nointerpolstats =', key_nointerpolstats
      WRITE(lun_standard,*)'   - max_transport       =', max_transport
      WRITE(lun_standard,*)'   - lmin                =', lmin
      WRITE(lun_standard,*)'   - lmax                =', lmax
    ENDIF

    !! BLINDAGE !!
    IF (max_transport < rZero) THEN
      WRITE(lun_error,*)''
      WRITE(lun_error,*)'Error: namelist file: max_transport has to be > 0 !'
      WRITE(lun_error,*)'Error: Stop Ariane'
      STOP
    ENDIF

  END SUBROUTINE sub_read_namelist_quantitative
  !!***
  !=========================================================================
  !!****f* mod_namelist/sub_read_namelist_qualitative()
  !! NAME
  !!   sub_read_namelist_qualitative()
  !!
  !! FUNCTION
  !!   READ  a namelist file.
  !!
  !! AUTHOR
  !!   * Origin  : Nicolas Grima
  !! 
  !! CREATION DATE
  !!   * January 2006
  !!
  !! HISTORY
  !!   Date (dd/mm/yyyy/) - Modification(s)
  !!
  !! ARGUMENTS
  !!
  !! TODO
  !!
  !! USED BY
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_read_namelist_qualitative(l_write)

    LOGICAL, INTENT(IN) :: l_write

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' - Reading QUALITATIVE item:'
    REWIND(unit=lun_nml)
    READ(unit=lun_nml, nml=QUALITATIVE)

    IF (delta_t == r_def) THEN
      delta_t = tunit
    ENDIF

    IF (frequency == i_defn) THEN
      frequency = ntfic
    ENDIF

    IF (nb_output == i_high) THEN
      nb_output = (INT(tunit) * ntfic * lmt) / (INT(delta_t) * frequency)
    ENDIF

    IF (l_write) THEN
      WRITE(lun_standard,*)'   - delta_t             =', delta_t
      WRITE(lun_standard,*)'   - frequency           =', frequency
      WRITE(lun_standard,*)'   - nb_output           =', nb_output
      WRITE(lun_standard,*)'   - key_region          =', key_region
      WRITE(lun_standard,*)'   - mask                =', mask
    ENDIF

  END SUBROUTINE sub_read_namelist_qualitative
  !!***
  !=========================================================================
  !!****f* mod_namelist/sub_read_namelist_opa()
  !! NAME
  !!   sub_read_namelist_opa()
  !!
  !! FUNCTION
  !!   READ the values of global key variables, parameters and some data
  !!   from a namelist file.
  !!
  !! AUTHOR
  !!   * Origin  : Nicolas Grima (November 2005)
  !! 
  !! CREATION DATE
  !!   * November 2005
  !!
  !! HISTORY
  !!   Date (dd/mm/yyyy/) - Modification(s)
  !!
  !! ARGUMENTS
  !!   No arguments.
  !!
  !! TODO
  !!   Add the possibility to specify n optionnal Logical Unit number 
  !!   (argument).
  !!
  !! USED BY
  !!   * trajec.f90: CALL sub_read_namelist()
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_read_namelist_opa(l_write)

    LOGICAL, INTENT(IN) :: l_write

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' - Reading OPAPARAM item:'
    REWIND(unit=lun_nml)
    READ(unit=lun_nml, nml=OPAPARAM)

    IF    ((TRIM(pivot) == 't').OR.(TRIM(pivot) == 'T')) THEN
      pivot = 'T'
    ELSEIF((TRIM(pivot) == 'f').OR.(TRIM(pivot) == 'F')) THEN
      pivot = 'F'
    ELSE
      WRITE(lun_error,*)''
      WRITE(lun_error,*)'Error: namelist/mod_namelist.f90'
      WRITE(lun_error,*)'       pivot value must be T or F'
      WRITE(lun_error,*)'       pivot= ', TRIM(pivot)
      STOP
    ENDIF

    IF (l_write) THEN
      WRITE(lun_standard,*)'   - imt                 =', imt
      WRITE(lun_standard,*)'   - jmt                 =', jmt
      WRITE(lun_standard,*)'   - kmt                 =', kmt
      WRITE(lun_standard,*)'   - lmt                 =', lmt
      WRITE(lun_standard,*)'   - key_computew        =', key_computew
      WRITE(lun_standard,*)'   - key_vvl             =', key_vvl
      WRITE(lun_standard,*)'   - key_partialsteps    =', key_partialsteps
      WRITE(lun_standard,*)'   - key_jfold           =', key_jfold
      WRITE(lun_standard,*)'   - pivot               =', pivot
      WRITE(lun_standard,*)'   - key_periodic        =', key_periodic
      WRITE(lun_standard,*)'   - w_surf_option       =', w_surf_option
      IF (TRIM(w_surf_option) == 'E-P-R') THEN
        WRITE(lun_standard,*)'   - epr_coef            =', epr_coef
      ENDIF
    ENDIF

    IF ((TRIM(w_surf_option) /= 'E-P-R').AND.&
         (TRIM(w_surf_option) /= 'zero').AND.&
         (TRIM(w_surf_option) /= '')) THEN
      WRITE(lun_error,*)'ERROR: w_sur_option has to be ''E-P-R'' or ''zero'' or '''''
      STOP
    ENDIF

    IF ((TRIM(w_surf_option) == 'E-P-R')) THEN
      IF (epr_coef == r_def) THEN
        WRITE(lun_error,*)'ERROR: If w_surf_option = ''E-P-R'' then epr_coef /= 0.'
        STOP
      ENDIF
    ENDIF

    IF ((key_sigma).AND.(key_computesigma)) THEN
      IF (l_write) THEN
        WRITE(lun_standard,*)'==========='
        WRITE(lun_standard,*) &
             '= WARNING = Key_computesigma and zsigma are defined in ARIANE item.'
        WRITE(lun_standard,*) &
             '            The next two values are NOT USED.'
        WRITE(lun_standard,*)'==========='
        WRITE(lun_standard,*) &
             '   - key_sigma (OBSOLETE and NOT used)= ', key_sigma
        WRITE(lun_standard,*) &
             '   - zsigma    (OBSOLETE and NOT used)= ', zsigma
        WRITE(lun_standard,*)' The values used are :'
        WRITE(lun_standard,*)'   - key_computesigma   =', key_computesigma
      ENDIF
      zsigma = zsigma_buffer
      IF (l_write) THEN
        WRITE(lun_standard,*)'   - zsigma      =', zsigma 
      ENDIF
    ELSEIF ((key_sigma).AND.(.NOT.key_computesigma)) THEN
      IF (l_write) THEN
        WRITE(lun_standard,*)'==========='
        WRITE(lun_standard,*) &
             '= WARNING = Ariane will use these values but in futur'
        WRITE(lun_standard,*) &
             '=========== please define key_computesigma and &
             &zsigma in ARIANE item.'
      ENDIF
      key_computesigma=key_sigma
      IF (l_write) THEN
        WRITE(lun_standard,*) &
             '   - key_sigma is OBSOLETE Ariane use key_computesigma'
        WRITE(lun_standard,*)'   - key_computesigma      =', key_computesigma

        IF (zsigma /= r_def) THEN
          WRITE(lun_standard,*) &
               '   - zsigma  (OBSOLETE in OPAPARAM)= ', zsigma
        ELSE
          WRITE(lun_standard,*) &
               '= WARNING = key_computesigma=.TRUE. and zsigma is ', r_def
        ENDIF
      ENDIF
    ENDIF

    !! NG: 16_10_2008 Bug fixed between key_computesigma and key_approximatesigma
    IF (key_computesigma) THEN
      key_approximatesigma=.FALSE.
    ELSE
      key_approximatesigma=.TRUE.
    ENDIF

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' - Reading ZONALCRT item:'
    REWIND(unit=lun_nml)
    READ(unit=lun_nml, nml=ZONALCRT)
    IF (l_write) THEN
      WRITE(lun_standard,*)'   - c_dir_zo            = ', TRIM(c_dir_zo)
      WRITE(lun_standard,*)'   - c_prefix_zo         = ', TRIM(c_prefix_zo)
      WRITE(lun_standard,*)'   - ind0_zo             = ', ind0_zo
      WRITE(lun_standard,*)'   - indn_zo             = ', indn_zo
      WRITE(lun_standard,*)'   - maxsize_zo          = ', maxsize_zo
      WRITE(lun_standard,*)'   - c_suffix_zo         = ', TRIM(c_suffix_zo)
      WRITE(lun_standard,*)'   - nc_var_zo           = ', TRIM(nc_var_zo)
      WRITE(lun_standard,*)'   - nc_var_eivu         = ', TRIM(nc_var_eivu)
      WRITE(lun_standard,*)'   - nc_att_mask_zo      = ', TRIM(nc_att_mask_zo)
    ENDIF

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' - Reading MERIDCRT item:'
    REWIND(unit=lun_nml)
    READ(unit=lun_nml, nml=MERIDCRT)
    IF (l_write) THEN
      WRITE(lun_standard,*)'   - c_dir_me            = ', TRIM(c_dir_me)
      WRITE(lun_standard,*)'   - c_prefix_me         = ', TRIM(c_prefix_me)
      WRITE(lun_standard,*)'   - ind0_me             = ', ind0_me
      WRITE(lun_standard,*)'   - indn_me             = ', indn_me
      WRITE(lun_standard,*)'   - maxsize_me          = ', maxsize_me
      WRITE(lun_standard,*)'   - c_suffix_me         = ', TRIM(c_suffix_me)
      WRITE(lun_standard,*)'   - nc_var_me           = ', TRIM(nc_var_me) 
      WRITE(lun_standard,*)'   - nc_var_eivv         = ', TRIM(nc_var_eivv)
      WRITE(lun_standard,*)'   - nc_att_mask_me      = ', TRIM(nc_att_mask_me)
    ENDIF

    IF (TRIM(w_surf_option) == 'E-P-R') THEN
      WRITE(lun_standard,*)''
      WRITE(lun_standard,*)' - Reading E-P-R item:'
      REWIND(unit=lun_nml)
      READ(unit=lun_nml, nml=E_P_R)
      IF (l_write) THEN
        WRITE(lun_standard,*)'   - c_dir_ep            = ', TRIM(c_dir_ep)
        WRITE(lun_standard,*)'   - c_prefix_ep         = ', TRIM(c_prefix_ep)
        WRITE(lun_standard,*)'   - ind0_ep             = ', ind0_ep
        WRITE(lun_standard,*)'   - indn_ep             = ', indn_ep
        WRITE(lun_standard,*)'   - maxsize_ep          = ', maxsize_ep
        WRITE(lun_standard,*)'   - c_suffix_ep         = ', TRIM(c_suffix_ep)
        WRITE(lun_standard,*)'   - nc_var_ep           = ', TRIM(nc_var_ep) 
        WRITE(lun_standard,*)'   - nc_att_mask_ep      = ', TRIM(nc_att_mask_ep)
      ENDIF
    ENDIF

    IF (.NOT.(key_computew)) THEN
      key_read_w=.TRUE.
      WRITE(lun_standard,*)''
      WRITE(lun_standard,*)' - Reading VERTICRT item:'
      REWIND(unit=lun_nml)
      READ(unit=lun_nml, nml=VERTICRT)
      IF (l_write) THEN
        WRITE(lun_standard,*)'   - c_dir_ve            = ', TRIM(c_dir_ve)
        WRITE(lun_standard,*)'   - c_prefix_ve         = ', TRIM(c_prefix_ve)
        WRITE(lun_standard,*)'   - ind0_ve             = ', ind0_ve
        WRITE(lun_standard,*)'   - indn_ve             = ', indn_ve
        WRITE(lun_standard,*)'   - maxsize_ve          = ', maxsize_ve
        WRITE(lun_standard,*)'   - c_suffix_ve         = ', TRIM(c_suffix_ve)
        WRITE(lun_standard,*)'   - nc_var_ve           = ', TRIM(nc_var_ve)
        WRITE(lun_standard,*)'   - nc_var_eivw         = ', TRIM(nc_var_eivw)
        WRITE(lun_standard,*)'   - nc_att_mask_ve      = ', TRIM(nc_att_mask_ve)
      ENDIF
    ENDIF

    if (key_vvl) then
       WRITE(lun_standard,*)''
       WRITE(lun_standard,*)' - Reading SSH item:'
       REWIND(unit=lun_nml)
       READ(unit=lun_nml, nml=SSH)
       IF (l_write) THEN
          WRITE(lun_standard,*)'   - c_dir_ssh            = ', TRIM(c_dir_ssh)
          WRITE(lun_standard,*)'   - c_prefix_ssh         = ', TRIM(c_prefix_ssh)
          WRITE(lun_standard,*)'   - ind0_ssh             = ', ind0_ssh
          WRITE(lun_standard,*)'   - indn_ssh             = ', indn_ssh
          WRITE(lun_standard,*)'   - maxsize_ssh          = ', maxsize_ssh
          WRITE(lun_standard,*)'   - c_suffix_ssh         = ', TRIM(c_suffix_ssh)
          WRITE(lun_standard,*)'   - nc_var_ssh           = ', TRIM(nc_var_ssh)
          WRITE(lun_standard,*)'   - nc_att_mask_ssh      = ', TRIM(nc_att_mask_ssh)
       ENDIF
    endif

    IF (key_alltracers) THEN

      WRITE(lun_standard,*)''
      WRITE(lun_standard,*)' - Reading TEMPERAT item:'
      REWIND(unit=lun_nml)
      READ(unit=lun_nml, nml=TEMPERAT)
      IF (l_write) THEN
        WRITE(lun_standard,*)'   - c_dir_te            = ', TRIM(c_dir_te)
        WRITE(lun_standard,*)'   - c_prefix_te         = ', TRIM(c_prefix_te)
        WRITE(lun_standard,*)'   - ind0_te             = ', ind0_te
        WRITE(lun_standard,*)'   - indn_te             = ', indn_te
        WRITE(lun_standard,*)'   - maxsize_te          = ', maxsize_te
        WRITE(lun_standard,*)'   - c_suffix_te         = ', TRIM(c_suffix_te)
        WRITE(lun_standard,*)'   - nc_var_te           = ', TRIM(nc_var_te) 
        WRITE(lun_standard,*)'   - nc_att_mask_te      = ', TRIM(nc_att_mask_te)
      ENDIF

      WRITE(lun_standard,*)''
      WRITE(lun_standard,*)' - Reading SALINITY item:'
      REWIND(unit=lun_nml)
      READ(unit=lun_nml, nml=SALINITY)
      IF (l_write) THEN
        WRITE(lun_standard,*)'   - c_dir_sa            = ', TRIM(c_dir_sa)
        WRITE(lun_standard,*)'   - c_prefix_sa         = ', TRIM(c_prefix_sa)
        WRITE(lun_standard,*)'   - ind0_sa             = ', ind0_sa
        WRITE(lun_standard,*)'   - indn_sa             = ', indn_sa
        WRITE(lun_standard,*)'   - maxsize_sa          = ', maxsize_sa
        WRITE(lun_standard,*)'   - c_suffix_sa         = ', TRIM(c_suffix_sa)
        WRITE(lun_standard,*)'   - nc_var_sa           = ', TRIM(nc_var_sa) 
        WRITE(lun_standard,*)'   - nc_att_mask_sa      = ', TRIM(nc_att_mask_sa)
      ENDIF

      IF (.NOT.(key_computesigma)) THEN
        WRITE(lun_standard,*)''
        WRITE(lun_standard,*)' - Reading DENSITY item:'
        REWIND(unit=lun_nml)
        READ(unit=lun_nml, nml=DENSITY)
        IF (l_write) THEN
          WRITE(lun_standard,*)'   - c_dir_de            = ', TRIM(c_dir_de)
          WRITE(lun_standard,*)'   - c_prefix_de         = ', TRIM(c_prefix_de)
          WRITE(lun_standard,*)'   - ind0_de             = ', ind0_de
          WRITE(lun_standard,*)'   - indn_de             = ', indn_de
          WRITE(lun_standard,*)'   - maxsize_de          = ', maxsize_de
          WRITE(lun_standard,*)'   - c_suffix_de         = ', TRIM(c_suffix_de)
          WRITE(lun_standard,*)'   - nc_var_de           = ', TRIM(nc_var_de) 
          WRITE(lun_standard,*)'   - nc_att_mask_de      = ', TRIM(nc_att_mask_de)
        ENDIF
      ENDIF

    ENDIF

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' - Reading MESH item:'
    REWIND(unit=lun_nml)
    READ(unit=lun_nml, nml=MESH)
    IF (l_write) THEN
      WRITE(lun_standard,*)'   - dir_mesh            = ', TRIM(dir_mesh)
      WRITE(lun_standard,*)'   - fn_mesh             = ', TRIM(fn_mesh)
      WRITE(lun_standard,*)'   - nc_var_xx_tt        = ', TRIM(nc_var_xx_tt)
      WRITE(lun_standard,*)'   - nc_var_xx_uu        = ', TRIM(nc_var_xx_uu)
      IF (TRIM(nc_var_xx_vv) /= c_def)  &
           WRITE(lun_standard,*)'   - nc_var_xx_vv        = OBSOLETE'
      IF (TRIM(nc_var_xx_ff) /= c_def)  &
           WRITE(lun_standard,*)'   - nc_var_xx_ff        = OBSOLETE'
      WRITE(lun_standard,*)'   - nc_var_yy_tt        = ', TRIM(nc_var_yy_tt)
      IF (TRIM(nc_var_yy_uu) /= c_def)  &
           WRITE(lun_standard,*)'   - nc_var_yy_uu        = OBSOLETE'
      WRITE(lun_standard,*)'   - nc_var_yy_vv        = ', TRIM(nc_var_yy_vv)
      IF (TRIM(nc_var_yy_ff) /= c_def)  &
           WRITE(lun_standard,*)'   - nc_var_yy_ff        = OBSOLETE'
      IF (TRIM(nc_var_zz_tt) /= c_def)  &
           WRITE(lun_standard,*)'   - nc_var_zz_tt        = OBSOLETE'
      WRITE(lun_standard,*)'   - nc_var_zz_ww        = ', TRIM(nc_var_zz_ww)
      WRITE(lun_standard,*)'   - nc_var_e2u          = ', TRIM(nc_var_e2u)
      WRITE(lun_standard,*)'   - nc_var_e1v          = ', TRIM(nc_var_e1v)
      WRITE(lun_standard,*)'   - nc_var_e1t          = ', TRIM(nc_var_e1t)
      WRITE(lun_standard,*)'   - nc_var_e2t          = ', TRIM(nc_var_e2t)
      IF (TRIM(mesh_type) == c_def) THEN
         WRITE(lun_standard,*)'   - nc_var_e3t          = ', TRIM(nc_var_e3t)
         if (key_vvl) then
            WRITE(lun_standard,*)'   - nc_var_mbathy       = ', TRIM(nc_var_mbathy)
         endif
      ELSEIF(TRIM(mesh_type) == 'nemov3') THEN
        WRITE(lun_standard,*)'   - nc_var_e3t2d        = ', TRIM(nc_var_e3t2D)
        WRITE(lun_standard,*)'   - nc_var_e3tz         = ', TRIM(nc_var_e3tz)
        WRITE(lun_standard,*)'   - nc_var_mbathy       = ', TRIM(nc_var_mbathy)
      ELSE
        STOP
      ENDIF
      WRITE(lun_standard,*)'   - nc_var_tmask        = ', TRIM(nc_var_tmask)
      WRITE(lun_standard,*)'   - nc_mask_val         = ', nc_mask_val
    ENDIF

  END SUBROUTINE sub_read_namelist_opa
  !!***
  !=========================================================================
  !!****f* mod_namelist/sub_read_namelist_roms()
  !! NAME
  !!   sub_read_namelist()
  !!
  !! FUNCTION
  !!   READ the values of global key variables, parameters and some data
  !!   from a namelist file.
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
  !!   No arguments.
  !!
  !! TODO
  !!   Add the possibility to specify n optionnal Logical Unit number 
  !!   (argument).
  !!
  !! USED BY
  !!   * trajec.f90: CALL sub_read_namelist()
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_read_namelist_roms(l_write)

    LOGICAL, INTENT(IN) :: l_write


    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' - Reading ROMSPARAM:'
    REWIND(unit=lun_nml)
    READ(unit=lun_nml, nml=ROMSPARAM)
    imt = xi_rho
    jmt = eta_rho
    kmt = s_w
    lmt = time
    IF (l_write) THEN
      WRITE(lun_standard,*)'   - xi_rho              =', imt
      WRITE(lun_standard,*)'   - eta_rho             =', jmt
      WRITE(lun_standard,*)'   - s_w                 =', kmt
      WRITE(lun_standard,*)'   - time                =', lmt
    ENDIF

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' - Reading ZONALCRT:'
    REWIND(unit=lun_nml)
    READ(unit=lun_nml, nml=ZONALCRT)
    IF (l_write) THEN
      WRITE(lun_standard,*)'   - c_dir_zo            = ', TRIM(c_dir_zo)
      WRITE(lun_standard,*)'   - c_prefix_zo         = ', TRIM(c_prefix_zo)
      WRITE(lun_standard,*)'   - ind0_zo             = ', ind0_zo
      WRITE(lun_standard,*)'   - indn_zo             = ', indn_zo
      WRITE(lun_standard,*)'   - maxsize_zo          = ', maxsize_zo
      WRITE(lun_standard,*)'   - c_suffix_zo         = ', TRIM(c_suffix_zo)
      WRITE(lun_standard,*)'   - nc_var_zo           = ', TRIM(nc_var_zo)
      WRITE(lun_standard,*)'   - nc_att_mask_zo      = ', TRIM(nc_att_mask_zo)
    ENDIF

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' - Reading MERIDCRT:'
    REWIND(unit=lun_nml)
    READ(unit=lun_nml, nml=MERIDCRT)
    IF (l_write) THEN
      WRITE(lun_standard,*)'   - c_dir_me            = ', TRIM(c_dir_me)
      WRITE(lun_standard,*)'   - c_prefix_me         = ', TRIM(c_prefix_me)
      WRITE(lun_standard,*)'   - ind0_me             = ', ind0_me
      WRITE(lun_standard,*)'   - indn_me             = ', indn_me
      WRITE(lun_standard,*)'   - maxsize_me          = ', maxsize_me
      WRITE(lun_standard,*)'   - c_suffix_me         = ', TRIM(c_suffix_me)
      WRITE(lun_standard,*)'   - nc_var_me           = ', TRIM(nc_var_me) 
      WRITE(lun_standard,*)'   - nc_att_mask_me      = ', TRIM(nc_att_mask_me)
    ENDIF

    IF (key_alltracers) THEN

      WRITE(lun_standard,*)''
      WRITE(lun_standard,*)' - Reading TEMPERAT:'
      REWIND(unit=lun_nml)
      READ(unit=lun_nml, nml=TEMPERAT)
      IF (l_write) THEN
        WRITE(lun_standard,*)'   - c_dir_te            = ', TRIM(c_dir_te)
        WRITE(lun_standard,*)'   - c_prefix_te         = ', TRIM(c_prefix_te)
        WRITE(lun_standard,*)'   - ind0_te             = ', ind0_te
        WRITE(lun_standard,*)'   - indn_te             = ', indn_te
        WRITE(lun_standard,*)'   - maxsize_te          = ', maxsize_te
        WRITE(lun_standard,*)'   - c_suffix_te         = ', TRIM(c_suffix_te)
        WRITE(lun_standard,*)'   - nc_var_te           = ', TRIM(nc_var_te) 
        WRITE(lun_standard,*)'   - nc_att_mask_te      = ', TRIM(nc_att_mask_te)
      ENDIF

      WRITE(lun_standard,*)''
      WRITE(lun_standard,*)' - Reading SALINITY:'
      REWIND(unit=lun_nml)
      READ(unit=lun_nml, nml=SALINITY)
      IF (l_write) THEN
        WRITE(lun_standard,*)'   - c_dir_sa            = ', TRIM(c_dir_sa)
        WRITE(lun_standard,*)'   - c_prefix_sa         = ', TRIM(c_prefix_sa)
        WRITE(lun_standard,*)'   - ind0_sa             = ', ind0_sa
        WRITE(lun_standard,*)'   - indn_sa             = ', indn_sa
        WRITE(lun_standard,*)'   - maxsize_sa          = ', maxsize_sa
        WRITE(lun_standard,*)'   - c_suffix_sa         = ', TRIM(c_suffix_sa)
        WRITE(lun_standard,*)'   - nc_var_sa           = ', TRIM(nc_var_sa) 
        WRITE(lun_standard,*)'   - nc_att_mask_sa      = ', TRIM(nc_att_mask_sa)
      ENDIF

    ENDIF

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' - Reading ZETA:'
    REWIND(unit=lun_nml)
    READ(unit=lun_nml, nml=ZETA)
    IF (l_write) THEN
      WRITE(lun_standard,*)'   - c_dir_ze            = ', TRIM(c_dir_ze)
      WRITE(lun_standard,*)'   - c_prefix_ze         = ', TRIM(c_prefix_ze)
      WRITE(lun_standard,*)'   - ind0_ze             = ', ind0_ze
      WRITE(lun_standard,*)'   - indn_ze             = ', indn_ze
      WRITE(lun_standard,*)'   - maxsize_ze          = ', maxsize_ze
      WRITE(lun_standard,*)'   - c_suffix_ze         = ', TRIM(c_suffix_ze)
      WRITE(lun_standard,*)'   - nc_var_ze           = ', TRIM(nc_var_ze) 
      WRITE(lun_standard,*)'   - nc_att_mask_ze      = ', TRIM(nc_att_mask_ze)
    ENDIF

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' - Reading GLOBALATT:'
    REWIND(unit=lun_nml)
    READ(unit=lun_nml, nml=GLOBALATT)
    IF (l_write) THEN
      WRITE(lun_standard,*)'   - dir_glbatt          = ', TRIM(dir_glbatt)
      WRITE(lun_standard,*)'   - fn_glbatt           = ', TRIM(fn_glbatt)
      WRITE(lun_standard,*)'   - nc_glbatt_hc        = ', TRIM(nc_glbatt_hc)
      WRITE(lun_standard,*)'   - nc_glbatt_sc_w      = ', TRIM(nc_glbatt_sc_w)
      WRITE(lun_standard,*)'   - nc_glbatt_Cs_w      = ', TRIM(nc_glbatt_Cs_w)
    ENDIF

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' - Reading GRDROMS:'
    REWIND(unit=lun_nml)
    READ(unit=lun_nml, nml=GRDROMS)
    IF (l_write) THEN
      WRITE(lun_standard,*)'   - dir_grd_roms        = ', TRIM(dir_grd_roms)
      WRITE(lun_standard,*)'   - fn_grd_roms         = ', TRIM(fn_grd_roms)
      WRITE(lun_standard,*)'   - nc_var_lon_rho_roms = ', TRIM(nc_var_lon_rho_roms)
      WRITE(lun_standard,*)'   - nc_var_lon_u_roms   = ', TRIM(nc_var_lon_u_roms)
      WRITE(lun_standard,*)'   - nc_var_lat_rho_roms = ', TRIM(nc_var_lat_rho_roms)
      WRITE(lun_standard,*)'   - nc_var_lat_v_roms   = ', TRIM(nc_var_lat_v_roms)
      WRITE(lun_standard,*)'   - nc_var_pm_roms      = ', TRIM(nc_var_pm_roms)
      WRITE(lun_standard,*)'   - nc_var_pn_roms      = ', TRIM(nc_var_pn_roms)
      WRITE(lun_standard,*)'   - nc_var_h_roms       = ', TRIM(nc_var_h_roms)
      WRITE(lun_standard,*)'   - nc_var_mask_rho_roms    = ', TRIM(nc_var_mask_rho_roms)
    ENDIF

  END SUBROUTINE sub_read_namelist_roms

  !!***
  !=========================================================================
  !!****f* mod_namelist/sub_read_namelist_mars()
  !! NAME
  !!   sub_read_namelist()
  !!
  !! FUNCTION
  !!   READ the values of global key variables, parameters and some data
  !!   from a namelist file.
  !!
  !! AUTHOR
  !!   * Origin  : Nicolas Grima (June 2018)
  !! 
  !! CREATION DATE
  !!   * April 2007
  !!
  !! HISTORY
  !!   Date (dd/mm/yyyy/) - Modification(s)
  !!
  !! ARGUMENTS
  !!   No arguments.
  !!
  !! TODO
  !!   Add the possibility to specify n optionnal Logical Unit number 
  !!   (argument).
  !!
  !! USED BY
  !!   * trajec.f90: CALL sub_read_namelist()
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_read_namelist_mars(l_write)

    LOGICAL, INTENT(IN) :: l_write

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' - Reading MARSPARAM:'
    REWIND(unit=lun_nml)
    READ(unit=lun_nml, nml=MARSPARAM)
    imt = x_t
    jmt = y_t
    kmt = sigma_t
    lmt = time
    IF (l_write) THEN
      WRITE(lun_standard,*)'   - x_t               =', imt
      WRITE(lun_standard,*)'   - y_t               =', jmt
      WRITE(lun_standard,*)'   - sigma_t           =', kmt
      WRITE(lun_standard,*)'   - time              =', lmt
    ENDIF

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' - Reading ZONALCRT:'
    REWIND(unit=lun_nml)
    READ(unit=lun_nml, nml=ZONALCRT)
    IF (l_write) THEN
      WRITE(lun_standard,*)'   - c_dir_zo            = ', TRIM(c_dir_zo)
      WRITE(lun_standard,*)'   - c_prefix_zo         = ', TRIM(c_prefix_zo)
      WRITE(lun_standard,*)'   - ind0_zo             = ', ind0_zo
      WRITE(lun_standard,*)'   - indn_zo             = ', indn_zo
      WRITE(lun_standard,*)'   - maxsize_zo          = ', maxsize_zo
      WRITE(lun_standard,*)'   - c_suffix_zo         = ', TRIM(c_suffix_zo)
      WRITE(lun_standard,*)'   - nc_var_zo           = ', TRIM(nc_var_zo)
      WRITE(lun_standard,*)'   - nc_var_eivu         = ', TRIM(nc_var_eivu)
      WRITE(lun_standard,*)'   - nc_att_mask_zo      = ', TRIM(nc_att_mask_zo)
    ENDIF

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' - Reading MERIDCRT:'
    REWIND(unit=lun_nml)
    READ(unit=lun_nml, nml=MERIDCRT)
    IF (l_write) THEN
      WRITE(lun_standard,*)'   - c_dir_me            = ', TRIM(c_dir_me)
      WRITE(lun_standard,*)'   - c_prefix_me         = ', TRIM(c_prefix_me)
      WRITE(lun_standard,*)'   - ind0_me             = ', ind0_me
      WRITE(lun_standard,*)'   - indn_me             = ', indn_me
      WRITE(lun_standard,*)'   - maxsize_me          = ', maxsize_me
      WRITE(lun_standard,*)'   - c_suffix_me         = ', TRIM(c_suffix_me)
      WRITE(lun_standard,*)'   - nc_var_me           = ', TRIM(nc_var_me) 
      WRITE(lun_standard,*)'   - nc_att_mask_me      = ', TRIM(nc_att_mask_me)
    ENDIF

    IF (key_alltracers) THEN

      WRITE(lun_standard,*)''
      WRITE(lun_standard,*)' - Reading TEMPERAT:'
      REWIND(unit=lun_nml)
      READ(unit=lun_nml, nml=TEMPERAT)
      IF (l_write) THEN
        WRITE(lun_standard,*)'   - c_dir_te            = ', TRIM(c_dir_te)
        WRITE(lun_standard,*)'   - c_prefix_te         = ', TRIM(c_prefix_te)
        WRITE(lun_standard,*)'   - ind0_te             = ', ind0_te
        WRITE(lun_standard,*)'   - indn_te             = ', indn_te
        WRITE(lun_standard,*)'   - maxsize_te          = ', maxsize_te
        WRITE(lun_standard,*)'   - c_suffix_te         = ', TRIM(c_suffix_te)
        WRITE(lun_standard,*)'   - nc_var_te           = ', TRIM(nc_var_te) 
        WRITE(lun_standard,*)'   - nc_att_mask_te      = ', TRIM(nc_att_mask_te)
      ENDIF

      WRITE(lun_standard,*)''
      WRITE(lun_standard,*)' - Reading SALINITY:'
      REWIND(unit=lun_nml)
      READ(unit=lun_nml, nml=SALINITY)
      IF (l_write) THEN
        WRITE(lun_standard,*)'   - c_dir_sa            = ', TRIM(c_dir_sa)
        WRITE(lun_standard,*)'   - c_prefix_sa         = ', TRIM(c_prefix_sa)
        WRITE(lun_standard,*)'   - ind0_sa             = ', ind0_sa
        WRITE(lun_standard,*)'   - indn_sa             = ', indn_sa
        WRITE(lun_standard,*)'   - maxsize_sa          = ', maxsize_sa
        WRITE(lun_standard,*)'   - c_suffix_sa         = ', TRIM(c_suffix_sa)
        WRITE(lun_standard,*)'   - nc_var_sa           = ', TRIM(nc_var_sa) 
        WRITE(lun_standard,*)'   - nc_att_mask_sa      = ', TRIM(nc_att_mask_sa)
      ENDIF

      IF (.NOT.(key_computesigma)) THEN
        WRITE(lun_standard,*)''
        WRITE(lun_standard,*)' - Reading DENSITY item:'
        REWIND(unit=lun_nml)
        READ(unit=lun_nml, nml=DENSITY)
        IF (l_write) THEN
          WRITE(lun_standard,*)'   - c_dir_de            = ', TRIM(c_dir_de)
          WRITE(lun_standard,*)'   - c_prefix_de         = ', TRIM(c_prefix_de)
          WRITE(lun_standard,*)'   - ind0_de             = ', ind0_de
          WRITE(lun_standard,*)'   - indn_de             = ', indn_de
          WRITE(lun_standard,*)'   - maxsize_de          = ', maxsize_de
          WRITE(lun_standard,*)'   - c_suffix_de         = ', TRIM(c_suffix_de)
          WRITE(lun_standard,*)'   - nc_var_de           = ', TRIM(nc_var_de) 
          WRITE(lun_standard,*)'   - nc_att_mask_de      = ', TRIM(nc_att_mask_de)
        ENDIF
      ENDIF

    ENDIF

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' - Reading ZETA:'
    REWIND(unit=lun_nml)
    READ(unit=lun_nml, nml=ZETA)
    IF (l_write) THEN
      WRITE(lun_standard,*)'   - c_dir_ze            = ', TRIM(c_dir_ze)
      WRITE(lun_standard,*)'   - c_prefix_ze         = ', TRIM(c_prefix_ze)
      WRITE(lun_standard,*)'   - ind0_ze             = ', ind0_ze
      WRITE(lun_standard,*)'   - indn_ze             = ', indn_ze
      WRITE(lun_standard,*)'   - maxsize_ze          = ', maxsize_ze
      WRITE(lun_standard,*)'   - c_suffix_ze         = ', TRIM(c_suffix_ze)
      WRITE(lun_standard,*)'   - nc_var_ze           = ', TRIM(nc_var_ze) 
      WRITE(lun_standard,*)'   - nc_att_mask_ze      = ', TRIM(nc_att_mask_ze)
    ENDIF

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' - Reading GRDMARS:'
    REWIND(unit=lun_nml)
    READ(unit=lun_nml, nml=GRDMARS)
    IF (l_write) THEN
      WRITE(lun_standard,*)'   - dir_grd_mars         = ', TRIM(dir_grd_mars)
      WRITE(lun_standard,*)'   - fn_grd_mars          = ', TRIM(fn_grd_mars)
      WRITE(lun_standard,*)'   - nc_var_lon_t_mars    = ', TRIM(nc_var_lon_t_mars)
      WRITE(lun_standard,*)'   - nc_var_lon_u_mars    = ', TRIM(nc_var_lon_u_mars)
      WRITE(lun_standard,*)'   - nc_var_lat_t_mars    = ', TRIM(nc_var_lat_t_mars)
      WRITE(lun_standard,*)'   - nc_var_lat_v_mars    = ', TRIM(nc_var_lat_v_mars)
      WRITE(lun_standard,*)'   - nc_var_e1t_mars      = ', TRIM(nc_var_e1t_mars)
      WRITE(lun_standard,*)'   - nc_var_e2t_mars      = ', TRIM(nc_var_e2t_mars)
      WRITE(lun_standard,*)'   - nc_var_e2u_mars      = ', TRIM(nc_var_e2u_mars)
      WRITE(lun_standard,*)'   - nc_var_e1v_mars      = ', TRIM(nc_var_e1v_mars)
      WRITE(lun_standard,*)'   - nc_var_bathy_t_mars  = ', TRIM(nc_var_bathy_t_mars)
      WRITE(lun_standard,*)'   - nc_var_hc_mars       = ', TRIM(nc_var_hc_mars)
      WRITE(lun_standard,*)'   - nc_var_Cs_w_mars     = ', TRIM(nc_var_Cs_w_mars)
      WRITE(lun_standard,*)'   - nc_var_sc_w_mars     = ', TRIM(nc_var_sc_w_mars)
      WRITE(lun_standard,*)'   - nc_var_mask_t_mars   = ', TRIM(nc_var_mask_t_mars)
    ENDIF

  END SUBROUTINE sub_read_namelist_mars

  !!***
  !=========================================================================
  !!****f* mod_namelist/sub_read_namelist_symphonie()
  !! NAME
  !!   sub_read_namelist()
  !!
  !! FUNCTION
  !!   READ the values of global key variables, parameters and some data
  !!   from a namelist file.
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
  !!   No arguments.
  !!
  !! TODO
  !!   Add the possibility to specify n optionnal Logical Unit number 
  !!   (argument).
  !!
  !! USED BY
  !!   * trajec.f90: CALL sub_read_namelist()
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_read_namelist_symphonie(l_write)

    LOGICAL, INTENT(IN) :: l_write

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' - Reading SYMPHONIEPARAM:'
    REWIND(unit=lun_nml)
    READ(unit=lun_nml, nml=SYMPHONIEPARAM)
    imt = x_dim
    jmt = y_dim
    kmt = z_dim + 1
    lmt = time
    IF (l_write) THEN
      WRITE(lun_standard,*)'   - x_dim               =', imt
      WRITE(lun_standard,*)'   - y_dim               =', jmt
      WRITE(lun_standard,*)'   - z_dim + 1           =', kmt
      WRITE(lun_standard,*)'   - time                =', lmt
    ENDIF

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' - Reading ZONALCRT:'
    REWIND(unit=lun_nml)
    READ(unit=lun_nml, nml=ZONALCRT)
    IF (l_write) THEN
      WRITE(lun_standard,*)'   - c_dir_zo            = ', TRIM(c_dir_zo)
      WRITE(lun_standard,*)'   - c_prefix_zo         = ', TRIM(c_prefix_zo)
      WRITE(lun_standard,*)'   - ind0_zo             = ', ind0_zo
      WRITE(lun_standard,*)'   - indn_zo             = ', indn_zo
      WRITE(lun_standard,*)'   - maxsize_zo          = ', maxsize_zo
      WRITE(lun_standard,*)'   - c_suffix_zo         = ', TRIM(c_suffix_zo)
      WRITE(lun_standard,*)'   - nc_var_zo           = ', TRIM(nc_var_zo)
      WRITE(lun_standard,*)'   - nc_var_eivu         = ', TRIM(nc_var_eivu)
      WRITE(lun_standard,*)'   - nc_att_mask_zo      = ', TRIM(nc_att_mask_zo)
    ENDIF

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' - Reading MERIDCRT:'
    REWIND(unit=lun_nml)
    READ(unit=lun_nml, nml=MERIDCRT)
    IF (l_write) THEN
      WRITE(lun_standard,*)'   - c_dir_me            = ', TRIM(c_dir_me)
      WRITE(lun_standard,*)'   - c_prefix_me         = ', TRIM(c_prefix_me)
      WRITE(lun_standard,*)'   - ind0_me             = ', ind0_me
      WRITE(lun_standard,*)'   - indn_me             = ', indn_me
      WRITE(lun_standard,*)'   - maxsize_me          = ', maxsize_me
      WRITE(lun_standard,*)'   - c_suffix_me         = ', TRIM(c_suffix_me)
      WRITE(lun_standard,*)'   - nc_var_me           = ', TRIM(nc_var_me) 
      WRITE(lun_standard,*)'   - nc_att_mask_me      = ', TRIM(nc_att_mask_me)
    ENDIF

    IF (key_alltracers) THEN

      WRITE(lun_standard,*)''
      WRITE(lun_standard,*)' - Reading TEMPERAT:'
      REWIND(unit=lun_nml)
      READ(unit=lun_nml, nml=TEMPERAT)
      IF (l_write) THEN
        WRITE(lun_standard,*)'   - c_dir_te            = ', TRIM(c_dir_te)
        WRITE(lun_standard,*)'   - c_prefix_te         = ', TRIM(c_prefix_te)
        WRITE(lun_standard,*)'   - ind0_te             = ', ind0_te
        WRITE(lun_standard,*)'   - indn_te             = ', indn_te
        WRITE(lun_standard,*)'   - maxsize_te          = ', maxsize_te
        WRITE(lun_standard,*)'   - c_suffix_te         = ', TRIM(c_suffix_te)
        WRITE(lun_standard,*)'   - nc_var_te           = ', TRIM(nc_var_te) 
        WRITE(lun_standard,*)'   - nc_att_mask_te      = ', TRIM(nc_att_mask_te)
      ENDIF

      WRITE(lun_standard,*)''
      WRITE(lun_standard,*)' - Reading SALINITY:'
      REWIND(unit=lun_nml)
      READ(unit=lun_nml, nml=SALINITY)
      IF (l_write) THEN
        WRITE(lun_standard,*)'   - c_dir_sa            = ', TRIM(c_dir_sa)
        WRITE(lun_standard,*)'   - c_prefix_sa         = ', TRIM(c_prefix_sa)
        WRITE(lun_standard,*)'   - ind0_sa             = ', ind0_sa
        WRITE(lun_standard,*)'   - indn_sa             = ', indn_sa
        WRITE(lun_standard,*)'   - maxsize_sa          = ', maxsize_sa
        WRITE(lun_standard,*)'   - c_suffix_sa         = ', TRIM(c_suffix_sa)
        WRITE(lun_standard,*)'   - nc_var_sa           = ', TRIM(nc_var_sa) 
        WRITE(lun_standard,*)'   - nc_att_mask_sa      = ', TRIM(nc_att_mask_sa)
      ENDIF

      IF (.NOT.(key_computesigma)) THEN
        WRITE(lun_standard,*)''
        WRITE(lun_standard,*)' - Reading DENSITY item:'
        REWIND(unit=lun_nml)
        READ(unit=lun_nml, nml=DENSITY)
        IF (l_write) THEN
          WRITE(lun_standard,*)'   - c_dir_de            = ', TRIM(c_dir_de)
          WRITE(lun_standard,*)'   - c_prefix_de         = ', TRIM(c_prefix_de)
          WRITE(lun_standard,*)'   - ind0_de             = ', ind0_de
          WRITE(lun_standard,*)'   - indn_de             = ', indn_de
          WRITE(lun_standard,*)'   - maxsize_de          = ', maxsize_de
          WRITE(lun_standard,*)'   - c_suffix_de         = ', TRIM(c_suffix_de)
          WRITE(lun_standard,*)'   - nc_var_de           = ', TRIM(nc_var_de) 
          WRITE(lun_standard,*)'   - nc_att_mask_de      = ', TRIM(nc_att_mask_de)
        ENDIF
      ENDIF

    ENDIF

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' - Reading SSE:'
    REWIND(unit=lun_nml)
    READ(unit=lun_nml, nml=SSE)
    IF (l_write) THEN
      WRITE(lun_standard,*)'   - c_dir_sse            = ', TRIM(c_dir_sse)
      WRITE(lun_standard,*)'   - c_prefix_sse         = ', TRIM(c_prefix_sse)
      WRITE(lun_standard,*)'   - ind0_sse             = ', ind0_sse
      WRITE(lun_standard,*)'   - indn_sse             = ', indn_sse
      WRITE(lun_standard,*)'   - maxsize_sse          = ', maxsize_sse
      WRITE(lun_standard,*)'   - c_suffix_sse         = ', TRIM(c_suffix_sse)
      WRITE(lun_standard,*)'   - nc_var_sse           = ', TRIM(nc_var_sse) 
      WRITE(lun_standard,*)'   - nc_att_mask_sse      = ', TRIM(nc_att_mask_sse)
    ENDIF

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' - Reading ZZ_TT_SIGMA:'
    REWIND(unit=lun_nml)
    READ(unit=lun_nml, nml=ZZ_TT_SIGMA)
    IF (l_write) THEN
      WRITE(lun_standard,*)'   - c_dir_zt            = ', TRIM(c_dir_zt)
      WRITE(lun_standard,*)'   - c_prefix_zt         = ', TRIM(c_prefix_zt)
      WRITE(lun_standard,*)'   - ind0_zt             = ', ind0_zt
      WRITE(lun_standard,*)'   - indn_zt             = ', indn_zt
      WRITE(lun_standard,*)'   - maxsize_zt          = ', maxsize_zt
      WRITE(lun_standard,*)'   - c_suffix_zt         = ', TRIM(c_suffix_zt)
      WRITE(lun_standard,*)'   - nc_var_zt           = ', TRIM(nc_var_zt) 
      WRITE(lun_standard,*)'   - nc_att_mask_zt      = ', TRIM(nc_att_mask_zt)
    ENDIF

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' - Reading GRDSYMPHONIE:'
    REWIND(unit=lun_nml)
    READ(unit=lun_nml, nml=GRDSYMPHONIE)
    IF (l_write) THEN
      WRITE(lun_standard,*)'   - dir_grd_symp         = ', TRIM(dir_grd_symp)
      WRITE(lun_standard,*)'   - fn_grd_symp          = ', TRIM(fn_grd_symp)
      WRITE(lun_standard,*)'   - cst_scale_factor     = ', cst_scale_factor
      WRITE(lun_standard,*)'   - nc_var_lon_t_symp    = ', TRIM(nc_var_lon_t_symp)
      WRITE(lun_standard,*)'   - nc_var_lon_u_symp    = ', TRIM(nc_var_lon_u_symp)
      WRITE(lun_standard,*)'   - nc_var_lat_t_symp    = ', TRIM(nc_var_lat_t_symp)
      WRITE(lun_standard,*)'   - nc_var_lat_v_symp    = ', TRIM(nc_var_lat_v_symp)
      WRITE(lun_standard,*)'   - nc_var_depth_t_symp  = ', TRIM(nc_var_depth_t_symp)
    ENDIF

  END SUBROUTINE sub_read_namelist_symphonie

  !!***
  !=========================================================================
  !!****f* mod_namelist/sub_read_namelist_B2C()
  !! NAME
  !!   sub_read_namelist_B2C()
  !!
  !! FUNCTION
  !!   READ the values of global key variables, parameters and some data
  !!   from a namelist file.
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
  !!   No arguments.
  !!
  !! TODO
  !!   Add the possibility to specify n optionnal Logical Unit number 
  !!   (argument).
  !!
  !! USED BY
  !!   * trajec.f90: CALL sub_read_namelist()
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_read_namelist_B2C(l_write)

    LOGICAL, INTENT(IN) :: l_write

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' - Reading B2C:'
    REWIND(unit=lun_nml)
    READ(unit=lun_nml, nml=B2C)

    imt = nb_dim_lon
    jmt = nb_dim_lat
    IF (key_add_bottom) THEN
      kmt = nb_dim_depth + 1
    ELSE
      kmt = nb_dim_depth
    ENDIF
    lmt = nb_dim_time

    key_read_w=.FALSE.

    IF (l_write) THEN
      WRITE(lun_standard,*)'   - B2C_grid_Z_or_Sigma = ', TRIM(B2C_grid_Z_or_Sigma)
      WRITE(lun_standard,*)'   - nb_dim_lon          = ', imt
      WRITE(lun_standard,*)'   - nb_dim_lat          = ', jmt
      WRITE(lun_standard,*)'   - nb_dim_depth        = ', kmt
      WRITE(lun_standard,*)'   - nb_dim_time         = ', lmt
      WRITE(lun_standard,*)'   - key_add_bottom      = ', key_add_bottom
      WRITE(lun_standard,*)'   - key_partialsteps    = ', key_partialsteps
      WRITE(lun_standard,*)'   - key_B2C_save_data   = ', key_B2C_save_data
      WRITE(lun_standard,*)'   - key_read_w (FORCED to FALSE) = ', key_read_w
      WRITE(lun_standard,*)'   - periodic_lon        = ', periodic_lon
      WRITE(lun_standard,*)'   - periodic_lat        = ', periodic_lat
    ENDIF

    IF ((TRIM(B2C_grid_Z_or_Sigma) == 'Z').OR.&
         (TRIM(B2C_grid_Z_or_Sigma) == 'z')) THEN

      B2C_grid_Z_or_Sigma ='Z'
      CALL sub_read_namelist_B2C_gridz(l_write)

    ELSEIF ((TRIM(B2C_grid_Z_or_Sigma) == 'S').OR.&
         (TRIM(B2C_grid_Z_or_Sigma) == 'Sigma').OR.&
         (TRIM(B2C_grid_Z_or_Sigma) == 'sigma').OR.&
         (TRIM(B2C_grid_Z_or_Sigma) == 'SIGMA')) THEN

      B2C_grid_Z_or_Sigma ='Sigma'
      WRITE(lun_standard,*)'STOP ==== B2C SIGMA GRID NOT YET IMPLEMENTED ==== STOP '
      STOP

    ELSE

      WRITE(lun_error,*)'ERROR === B2C_grid_Z_or_Sigma is not correct === ERROR'
      WRITE(lun_error,*)'B2C_grid_Z_or_Sigma should be: Z or Sigma'
      STOP

    ENDIF

    IF (periodic_lon) THEN

      WRITE(lun_standard,*)' '
      WRITE(lun_standard,*)'STOP ==== periodic_lon IS NOT YET IMPLEMENTED ==== STOP '
      WRITE(lun_standard,*)' '
      STOP

    ENDIF

    IF (periodic_lat) THEN

      WRITE(lun_standard,*)' '
      WRITE(lun_standard,*)'STOP ==== periodic_lat IS NOT YET IMPLEMENTED ==== STOP '
      WRITE(lun_standard,*)' '
      STOP

    ENDIF

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' - Reading ZONALCRT item:'
    REWIND(unit=lun_nml)
    READ(unit=lun_nml, nml=ZONALCRT)
    IF (l_write) THEN
      WRITE(lun_standard,*)'   - c_dir_zo            = ', TRIM(c_dir_zo)
      WRITE(lun_standard,*)'   - c_prefix_zo         = ', TRIM(c_prefix_zo)
      WRITE(lun_standard,*)'   - ind0_zo             = ', ind0_zo
      WRITE(lun_standard,*)'   - indn_zo             = ', indn_zo
      WRITE(lun_standard,*)'   - maxsize_zo          = ', maxsize_zo
      WRITE(lun_standard,*)'   - c_suffix_zo         = ', TRIM(c_suffix_zo)
      WRITE(lun_standard,*)'   - nc_var_zo           = ', TRIM(nc_var_zo)
      WRITE(lun_standard,*)'   - nc_var_eivu         = ', TRIM(nc_var_eivu)
      WRITE(lun_standard,*)'   - nc_att_mask_zo      = ', TRIM(nc_att_mask_zo)
    ENDIF

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' - Reading MERIDCRT item:'
    REWIND(unit=lun_nml)
    READ(unit=lun_nml, nml=MERIDCRT)
    IF (l_write) THEN
      WRITE(lun_standard,*)'   - c_dir_me            = ', TRIM(c_dir_me)
      WRITE(lun_standard,*)'   - c_prefix_me         = ', TRIM(c_prefix_me)
      WRITE(lun_standard,*)'   - ind0_me             = ', ind0_me
      WRITE(lun_standard,*)'   - indn_me             = ', indn_me
      WRITE(lun_standard,*)'   - maxsize_me          = ', maxsize_me
      WRITE(lun_standard,*)'   - c_suffix_me         = ', TRIM(c_suffix_me)
      WRITE(lun_standard,*)'   - nc_var_me           = ', TRIM(nc_var_me) 
      WRITE(lun_standard,*)'   - nc_var_eivv         = ', TRIM(nc_var_eivv)
      WRITE(lun_standard,*)'   - nc_att_mask_me      = ', TRIM(nc_att_mask_me)
    ENDIF

    IF (key_read_w) THEN
      WRITE(lun_standard,*)''
      WRITE(lun_standard,*)' - Reading VERTICRT item:'
      REWIND(unit=lun_nml)
      READ(unit=lun_nml, nml=VERTICRT)
      IF (l_write) THEN
        WRITE(lun_standard,*)'   - c_dir_ve            = ', TRIM(c_dir_ve)
        WRITE(lun_standard,*)'   - c_prefix_ve         = ', TRIM(c_prefix_ve)
        WRITE(lun_standard,*)'   - ind0_ve             = ', ind0_ve
        WRITE(lun_standard,*)'   - indn_ve             = ', indn_ve
        WRITE(lun_standard,*)'   - maxsize_ve          = ', maxsize_ve
        WRITE(lun_standard,*)'   - c_suffix_ve         = ', TRIM(c_suffix_ve)
        WRITE(lun_standard,*)'   - nc_var_ve           = ', TRIM(nc_var_ve)
        WRITE(lun_standard,*)'   - nc_var_eivw         = ', TRIM(nc_var_eivw)
        WRITE(lun_standard,*)'   - nc_att_mask_ve      = ', TRIM(nc_att_mask_ve)
      ENDIF
    ENDIF


    IF (key_alltracers) THEN

      WRITE(lun_standard,*)''
      WRITE(lun_standard,*)' - Reading TEMPERAT item:'
      REWIND(unit=lun_nml)
      READ(unit=lun_nml, nml=TEMPERAT)
      IF (l_write) THEN
        WRITE(lun_standard,*)'   - c_dir_te            = ', TRIM(c_dir_te)
        WRITE(lun_standard,*)'   - c_prefix_te         = ', TRIM(c_prefix_te)
        WRITE(lun_standard,*)'   - ind0_te             = ', ind0_te
        WRITE(lun_standard,*)'   - indn_te             = ', indn_te
        WRITE(lun_standard,*)'   - maxsize_te          = ', maxsize_te
        WRITE(lun_standard,*)'   - c_suffix_te         = ', TRIM(c_suffix_te)
        WRITE(lun_standard,*)'   - nc_var_te           = ', TRIM(nc_var_te) 
        WRITE(lun_standard,*)'   - nc_att_mask_te      = ', TRIM(nc_att_mask_te)
      ENDIF

      WRITE(lun_standard,*)''
      WRITE(lun_standard,*)' - Reading SALINITY item:'
      REWIND(unit=lun_nml)
      READ(unit=lun_nml, nml=SALINITY)
      IF (l_write) THEN
        WRITE(lun_standard,*)'   - c_dir_sa            = ', TRIM(c_dir_sa)
        WRITE(lun_standard,*)'   - c_prefix_sa         = ', TRIM(c_prefix_sa)
        WRITE(lun_standard,*)'   - ind0_sa             = ', ind0_sa
        WRITE(lun_standard,*)'   - indn_sa             = ', indn_sa
        WRITE(lun_standard,*)'   - maxsize_sa          = ', maxsize_sa
        WRITE(lun_standard,*)'   - c_suffix_sa         = ', TRIM(c_suffix_sa)
        WRITE(lun_standard,*)'   - nc_var_sa           = ', TRIM(nc_var_sa) 
        WRITE(lun_standard,*)'   - nc_att_mask_sa      = ', TRIM(nc_att_mask_sa)
      ENDIF

      IF (.NOT.(key_computesigma)) THEN
        WRITE(lun_standard,*)''
        WRITE(lun_standard,*)' - Reading DENSITY item:'
        REWIND(unit=lun_nml)
        READ(unit=lun_nml, nml=DENSITY)
        IF (l_write) THEN
          WRITE(lun_standard,*)'   - c_dir_de            = ', TRIM(c_dir_de)
          WRITE(lun_standard,*)'   - c_prefix_de         = ', TRIM(c_prefix_de)
          WRITE(lun_standard,*)'   - ind0_de             = ', ind0_de
          WRITE(lun_standard,*)'   - indn_de             = ', indn_de
          WRITE(lun_standard,*)'   - maxsize_de          = ', maxsize_de
          WRITE(lun_standard,*)'   - c_suffix_de         = ', TRIM(c_suffix_de)
          WRITE(lun_standard,*)'   - nc_var_de           = ', TRIM(nc_var_de) 
          WRITE(lun_standard,*)'   - nc_att_mask_de      = ', TRIM(nc_att_mask_de)
        ENDIF
      ENDIF

    ENDIF

  END SUBROUTINE sub_read_namelist_B2C

  !!***
  !=========================================================================
  !!****f* mod_namelist/sub_read_namelist_B2C_gridz()
  !! NAME
  !!   sub_read_namelist_B2C_gridz()
  !!
  !! FUNCTION
  !!   READ the values of global key variables, parameters and some data
  !!   from a namelist file.
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
  !!   No arguments.
  !!
  !! TODO
  !!   Add the possibility to specify n optionnal Logical Unit number 
  !!   (argument).
  !!
  !! USED BY
  !!   * trajec.f90: CALL sub_read_namelist()
  !!
  !! SOURCE
  !!=======================================================================
  SUBROUTINE sub_read_namelist_B2C_gridz(l_write)

    LOGICAL, INTENT(IN) :: l_write

    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)' - Reading B2CGRIDZ:'
    REWIND(unit=lun_nml)
    READ(unit=lun_nml, nml=B2CGRIDZ)

    IF (l_write) THEN
      WRITE(lun_standard,*)'   - dir_B2C_grid        = ', TRIM(dir_B2C_grid)
      WRITE(lun_standard,*)'   - file_name_B2C_grid  = ', TRIM(file_name_B2C_grid)
      WRITE(lun_standard,*)'   - nc_var_xx_tt        = ', TRIM(nc_var_xx_tt)
      WRITE(lun_standard,*)'   - nc_var_xx_uu        = ', TRIM(nc_var_xx_uu)
      WRITE(lun_standard,*)'   - nc_var_yy_tt        = ', TRIM(nc_var_yy_tt)
      WRITE(lun_standard,*)'   - nc_var_yy_vv        = ', TRIM(nc_var_yy_vv)
      WRITE(lun_standard,*)'   - nc_var_zz_ww        = ', TRIM(nc_var_zz_ww)
      WRITE(lun_standard,*)'   - nc_var_e1f          = ', TRIM(nc_var_e1f)
      WRITE(lun_standard,*)'   - nc_var_e2f          = ', TRIM(nc_var_e2f)
      WRITE(lun_standard,*)'   - nc_var_e3f          = ', TRIM(nc_var_e3f)
      WRITE(lun_standard,*)'   - nc_var_e2u          = ', TRIM(nc_var_e2u)
      WRITE(lun_standard,*)'   - nc_var_e1v          = ', TRIM(nc_var_e1v)
      WRITE(lun_standard,*)'   - nc_var_e1t          = ', TRIM(nc_var_e1t)
      WRITE(lun_standard,*)'   - nc_var_e2t          = ', TRIM(nc_var_e2t)
      WRITE(lun_standard,*)'   - nc_var_e3t          = ', TRIM(nc_var_e3t)
      WRITE(lun_standard,*)'   - nc_var_tmask        = ', TRIM(nc_var_tmask)
      WRITE(lun_standard,*)'   - nc_mask_val         = ', nc_mask_val
    ENDIF

  END SUBROUTINE sub_read_namelist_B2C_gridz
  !!***
END MODULE mod_namelist
