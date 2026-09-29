!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify it    !!
!! under the terms of the GNU General Public License as published by !!
!! the Free Software Foundation, either version 3 of the License, or !!
!! (at your option) any later version.                               !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS. if not, see: http://www.gnu.org/licenses/gpl.txt  !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module two_stream_gray_rad_mod

! ==================================================================================
! ==================================================================================

   use fms_mod,               only: open_file, check_nml_error, &
                                    mpp_pe, close_file, error_mesg, &
                                    NOTE, FATAL,  uppercase

   use constants_mod,         only: stefan, cp_air, grav, pstd_mks, pstd_mks_earth, seconds_per_sol, orbital_period, radius, rdgas

   use    diag_manager_mod,   only: register_diag_field, send_data

   use    time_manager_mod,   only: time_type, length_of_year, length_of_day, &
                                    operator(+), operator(-), operator(/=), get_time

   use astronomy_mod,         only: astronomy_init, diurnal_solar

   use interpolator_mod,      only: interpolate_type, interpolator_init, interpolator, ZERO, interpolator_end

!==================================================================================
implicit none
private
!==================================================================================

! version information

character(len=128) :: version='$Id: two_stream_gray_rad.F90,v 1.1.2.1 2016/02/23$'
character(len=128) :: tag='Two-stream gray atmosphere'

! Version Details
! [2016/01/28] <James Penn>:  GFDL Astronomy module used to calculate insolation
!   - Diurnal & Seasonal cycle implemented in astronomy_mod.
!   - set `do_seasonal` to .true. in namelist to activate.
!   - default (`do_seasonal = .false.`) remains as permanent global equinox.
!   - tidal locking can be achieved with by setting `2*pi / orbital_period = omega`.
! [2016/02/23] <James Penn>:  Window parameterisation in the longwave stream
!   - reference: Ruth Geen etal, GRL 2016 (supp. information).
!   - Set `wv_exponent` and `solar_exponent` to -1 in namelist to activate.
!   - Moisture and CO2 feedback on SW and LW controlled by `ir_tau...` constants.
!
! Generic condensate cloud scheme (single species, grey)
!   The default condensation curve is MgSiO3 (Roman & Rauscher 2019, ApJ 872, 1).
!   Set the cld_* namelist values to model a different cloud species.
!
!   Condensation curve:
!   - set `do_cld_condensation` to .true. to activate.
!   - condensation temperature from the general log-linear form
!         log10(P_bar) * cld_p_coeff = cld_A - cld_B/T - cld_mh_coeff * log10(cld_mh)
!     =>  Tcond = cld_B / (cld_A - cld_p_coeff*log10(P_bar) - cld_mh_coeff*log10(cld_mh))
!     Defaults (cld_A = 13.43, cld_B = 28665, cld_p_coeff = cld_mh_coeff = 1)
!     give the MgSiO3 curve. Fits of the form
!         1e4/Tcond = a - b*log10(P_bar) - c*[M/H]
!     map onto this with cld_B = 1e4, cld_A = a, cld_p_coeff = b, cld_mh_coeff = c.
!   - `cld_mh` is the non-log metallicity relative to solar (1 = solar, 10 = 10x solar).
!   - `cld_species` is a label only, used in log messages and diagnostic long names.
!   - outputs cld_tcond (K), cld_cond_flag (1 where T < Tcond, 0 elsewhere) and
!     cld_flag (the cloud actually seen by the radiation, after the
!     vertical-extent limit below).
!
!   Cloud vertical extent:
!   - `cld_extent` = 'none' (default), 'thick' or 'thin'. Cloud base is the
!     deepest (highest-pressure) level in a column with T < Tcond.
!       'none'  : every cell with T < Tcond is cloudy.
!       'thick' : cloud from the base up to `cld_thick_ptop` (Pa, default
!                 100 Pa = 1 mbar).
!       'thin'  : cloud from the base up `cld_thin_nscale` pressure scale
!                 heights (default 1), i.e. up to p_base*exp(-cld_thin_nscale).
!     Cells between the base and the limit are cloudy only where T < Tcond.
!     The limit is applied to cld_flag, so it affects the TOA albedo and the
!     LW and SW opacities. Outputs cld_p_cbase and cld_p_ctop (Pa, 0 if cloud-free).
!
!   Cloud-dependent TOA albedo:
!   - set `do_cld_toa_albedo` to .true. to set TOA albedo = `cld_toa_albedo`
!     (default 0.6) in any column with at least one cloudy level (cld_flag).
!     Requires do_cld_condensation = .true. Outputs cld_toa_albedo.
!
!   Cloud optical depth (Roman & Rauscher 2019, ApJ 872, 1, Eqs. 1-2):
!   - per cloudy cell (cld_flag), mass of condensible gas and optical depth
!         m_g   = (dP/g) * cld_x_gas * (cld_mu_gas / mu_atm)
!         tau   = 3 * m_g * Q * cld_f_cond / (4 * cld_r_part * cld_rho_part)
!     with Q = cld_Q_lw_eff (longwave) or cld_Q_sw_eff (shortwave).
!     User inputs: cld_x_gas (mole fraction of the condensible gas), cld_mu_gas
!     (its molecular weight, g/mol), cld_Q_lw_eff, cld_Q_sw_eff, cld_f_cond,
!     cld_r_part, cld_rho_part. mu_atm is derived from the model's specific
!     gas constant (mu_atm = universal_gas_const / rdgas), so it matches
!     whatever working fluid this experiment's constants_nml has set.
!     Outputs cld_mgas (kg/m2) when either opacity is active.
!
!   Cloud longwave opacity:
!   - set `do_cld_lw_opacity` to .true. to add tau_lw_cloud to cloudy cells,
!     so tau_lw = tau_lw_atm + tau_lw_cloud there. Requires
!     do_cld_condensation = .true., rad_scheme = "frierson", and
!     do_normal_integration_method = .true. Outputs cld_tau_lw.
!
!   Cloud shortwave opacity:
!   - set `do_cld_sw_opacity` to .true. to add tau_sw_cloud to cloudy cells,
!     so tau_sw = tau_sw_atm + tau_sw_cloud. Because the Frierson SW stream
!     is computed from the cumulative optical depth at half levels
!     (sw_down = (1-A_toa) * S * exp(-sw_tau)), the per-layer cloud optical
!     depth is summed downward from the model top and added to sw_tau.
!     The grey SW stream is pure attenuation (no scattering), so cloud SW
!     extinction is deposited as heating in the cloudy layers; reflection by
!     the cloud is represented only by do_cld_toa_albedo. Requires
!     do_cld_condensation = .true. and rad_scheme = "frierson" or "byrne".
!     Outputs cld_tau_sw.
!==================================================================================

! public interfaces

public :: two_stream_gray_rad_init, two_stream_gray_rad_down, two_stream_gray_rad_up, two_stream_gray_rad_end
!==================================================================================


! module variables
logical :: initialized     = .false.
real    :: solar_constant  = 1360.0
real    :: del_sol         = 1.4
real    :: del_sw          = 0.0
real    :: ir_tau_eq       = 6.0
real    :: ir_tau_pole     = 1.5
real    :: atm_abs         = 0.0
real    :: cloud_albedo    = 0.0
real    :: odp = 1.0  ! odp = optical depth parameter. Used to simulate GHG concentration, taken from Tapio FMS, added by Nathanael Wong
real    :: sw_diff         = 0.0
real    :: linear_tau      = 0.1
real    :: wv_exponent     = 4.0
real    :: solar_exponent  = 4.0
logical :: do_seasonal     = .false.
integer :: solday          = -10  !s Day of year to run perpetually if do_seasonal=True and solday>0
real    :: equinox_day     = 0.75 !s Fraction of year [0,1] where NH autumn equinox occurs (only really useful if calendar has defined months).
logical :: use_time_average_coszen = .false. !s if .true., then time-averaging is done on coszen so that insolation doesn't depend on timestep
real    :: dt_rad_avg     = -1
logical :: do_tl = .false.
logical :: do_closein = .false. 
logical :: do_normal_integration_method = .true. 
real :: R_stellar = 476700000. 
real :: d_stellar = 1117512000.0

! ---------------------------------------------------------------------------------
! Generic condensate cloud scheme (default condensation curve = MgSiO3, R&R 2019)
! ---------------------------------------------------------------------------------
character(len=32) :: cld_species = 'MgSiO3'  ! label only (log messages, diag long names)

! condensation curve: Tcond = cld_B / (cld_A - cld_p_coeff*log10(P_bar) - cld_mh_coeff*log10(cld_mh))
logical :: do_cld_condensation = .false.
real    :: cld_mh       = 100.0      ! non-log metallicity relative to solar (1 = solar, 10 = 10x solar)
real    :: cld_A        = 13.43    ! intercept (MgSiO3)
real    :: cld_B        = 28665.0  ! slope in 1/T, K (MgSiO3)
real    :: cld_p_coeff  = 1.0      ! coefficient on log10(P_bar)
real    :: cld_mh_coeff = 1.0      ! coefficient on log10(cld_mh)
real, parameter :: pa_per_bar = 1.0e5

! cloud vertical extent (applied to cld_flag, so it limits the TOA albedo
! and the LW and SW opacities)
!   'none'  : every cell with T < Tcond is cloudy
!   'thick' : from the deepest T < Tcond level up to cld_thick_ptop
!   'thin'  : from the deepest T < Tcond level up cld_thin_nscale pressure
!             scale heights, i.e. to p_base * exp(-cld_thin_nscale)
character(len=16) :: cld_extent = 'none'
real    :: cld_thick_ptop  = 100.0   ! Pa (1 mbar), cloud top for 'thick'
real    :: cld_thin_nscale = 1.0     ! scale heights above cloud base for 'thin'
integer, parameter :: CLD_EXT_NONE = 0, CLD_EXT_THICK = 1, CLD_EXT_THIN = 2
integer :: cld_ext = CLD_EXT_NONE

! cloud-dependent TOA albedo
logical :: do_cld_toa_albedo = .false.
real    :: cld_toa_albedo    = 0.6   ! TOA albedo in columns with any cloudy level

! cloud optical depth (Roman & Rauscher 2019, ApJ 872, 1, eqs. 1-2)
logical :: do_cld_lw_opacity = .false.
logical :: do_cld_sw_opacity = .false.
real    :: cld_x_gas    = 0.0     ! mole fraction of the condensible gas (chi_g)
real    :: cld_mu_gas   = 100.39  ! molecular weight of the condensible species, g/mol (MgSiO3)
real    :: cld_Q_lw_eff = 0.01    ! particle extinction efficiency, longwave (Q_e,LW)
real    :: cld_Q_sw_eff = 0.01    ! particle extinction efficiency, shortwave (Q_e,SW)
real    :: cld_f_cond   = 0.1     ! fraction of the condensible gas that forms condensate
real    :: cld_r_part   = 0.2e-6  ! particle radius, m
real    :: cld_rho_part = 3200.0  ! particle density, kg/m^3 (MgSiO3, R&R 2019 Table 2)

! universal gas constant, J/(kmol K); mean atmospheric molecular weight is
! derived from this and the model's specific gas constant (rdgas), so it is
! consistent with whatever this experiment's constants_nml has set.
real, parameter :: universal_gas_const = 8314.4621

real, allocatable, dimension(:,:,:) :: tcond_cld, cld_flag, cld_cond_flag
real, allocatable, dimension(:,:,:) :: mgas_cld, tau_lw_cld, tau_sw_cld
real, allocatable, dimension(:,:)   :: p_cbase_cld, p_ctop_cld   ! Pa, 0 if cloud-free
real, allocatable, dimension(:,:)   :: toa_albedo

character(len=32) :: rad_scheme = 'frierson'

integer, parameter :: B_GEEN = 1,  B_FRIERSON = 2, &
                      B_BYRNE = 3, B_SCHNEIDER_LIU=4
integer, private :: sw_scheme = B_FRIERSON
integer, private :: lw_scheme = B_FRIERSON

! constants for RG radiation version
real    :: ir_tau_co2_win  = 0.2150
real    :: ir_tau_wv_win1  = 147.11
real    :: ir_tau_wv_win2  = 1.0814e4
real    :: ir_tau_co2      = 0.1
real    :: ir_tau_wv1      = 23.8
real    :: ir_tau_wv2      = 254.0
real    :: window          = 0.3732
real    :: carbon_conc     = 360.0

! constants for SCHNEIDER_LIU radiation version
real    :: single_albedo      = 0.8
real    :: back_scatter       = 0.398
real    :: lw_tau_0_gp        = 80.0
real    :: sw_tau_0_gp        = 3.0
real    :: lw_tau_exponent_gp = 2.0
real    :: sw_tau_exponent_gp = 1.0
real    :: diabatic_acce = 1.0
real,save :: gp_albedo, Ga_asym, g_asym

! parameters for Byrne and OGorman radiation scheme
real :: bog_a = 0.8678
real :: bog_b = 1997.9
real :: bog_mu = 1.0

real, allocatable, dimension(:,:)   :: insolation, p2, lw_tau_0, sw_tau_0 !s albedo now defined in mixed_layer_init
real, allocatable, dimension(:,:)   :: lat_tl, lat_tl_insol
real, allocatable, dimension(:,:)   :: b_surf, b_surf_gp
real, allocatable, dimension(:,:,:) :: b, tdt_rad, tdt_solar
real, allocatable, dimension(:,:,:) :: lw_up, lw_down, lw_flux, sw_up, sw_down, sw_flux, rad_flux
real, allocatable, dimension(:,:,:) :: lw_tau, sw_tau, lw_dtrans
real, allocatable, dimension(:,:)   :: olr, net_lw_surf, toa_sw_in, coszen, fracsun

! window parameterisation (RG, 2015)
real, allocatable, dimension(:,:,:) :: lw_up_win, lw_down_win, lw_dtrans_win
real, allocatable, dimension(:,:,:) :: b_win, sw_dtrans
real, allocatable, dimension(:,:)   :: sw_wv, del_sol_tau, sw_tau_k, lw_del_tau, lw_del_tau_win

real, save :: pi, deg_to_rad , rad_to_deg

!extras for reading in co2 concentration
logical                             :: do_read_co2=.false.
type(interpolate_type),save         :: co2_interp           ! use external file for co2
character(len=256)                  :: co2_file='co2'       !  file name of co2 file to read
character(len=256)                  :: co2_variable_name='co2'       !  file name of co2 file to read

logical :: do_toa_albedo = .false.

namelist /two_stream_gray_rad_nml/ solar_constant, del_sol, &
           ir_tau_eq, ir_tau_pole, odp, atm_abs, cloud_albedo, &
           sw_diff, linear_tau, del_sw, wv_exponent, &
           solar_exponent, do_seasonal, &
           ir_tau_co2_win, ir_tau_wv_win1, ir_tau_wv_win2, &
           ir_tau_co2, ir_tau_wv1, ir_tau_wv2, &
           window, carbon_conc, rad_scheme, &
           do_read_co2, co2_file, co2_variable_name, solday, equinox_day, bog_a, bog_b, bog_mu, &
           use_time_average_coszen, dt_rad_avg, &
           diabatic_acce, &
           do_toa_albedo, do_tl, do_closein, R_stellar, d_stellar, do_normal_integration_method, &
           cld_species, &
           do_cld_condensation, cld_mh, cld_A, cld_B, cld_p_coeff, cld_mh_coeff, &
           cld_extent, cld_thick_ptop, cld_thin_nscale, &
           do_cld_toa_albedo, cld_toa_albedo, &
           do_cld_lw_opacity, do_cld_sw_opacity, cld_x_gas, cld_mu_gas, &
           cld_Q_lw_eff, cld_Q_sw_eff, cld_f_cond, cld_r_part, cld_rho_part
!==================================================================================
!-------------------- diagnostics fields -------------------------------

integer :: id_olr, id_swdn_sfc, id_swdn_toa, id_net_lw_surf, id_lwdn_sfc, id_lwup_sfc, &
           id_tdt_rad, id_tdt_solar, id_flux_rad, id_flux_lw, id_flux_sw, id_coszen, id_fracsun, &
           id_lw_dtrans, id_lw_dtrans_win, id_sw_dtrans, id_co2, id_nudge_flux, id_tau, &
           id_tcond_cld, id_cld_flag, id_toa_albedo, id_mgas_cld, id_tau_lw_cld, id_tau_sw_cld, &
           id_cld_cond_flag, id_p_cbase_cld, id_p_ctop_cld

character(len=10), parameter :: mod_name = 'two_stream'

real :: missing_value = -999.


contains


! ==================================================================================
! ==================================================================================


subroutine two_stream_gray_rad_init(is, ie, js, je, num_levels, axes, Time, lonb, latb, dt_real)

!-------------------------------------------------------------------------------------
integer, intent(in), dimension(4) :: axes
type(time_type), intent(in)       :: Time
integer, intent(in)               :: is, ie, js, je, num_levels
real ,dimension(:,:),intent(in),optional :: lonb,latb !s Changed to 2d arrays as 2013 interpolator expects this.
real, intent(in)                  :: dt_real !s atmospheric timestep, used for radiation averaging
!-------------------------------------------------------------------------------------
integer, dimension(3) :: half = (/1,2,4/)
integer :: ierr, io, unit
character(len=64) :: sp
!-----------------------------------------------------------------------------------------
! read namelist and copy to logfile

unit = open_file ('input.nml', action='read')
ierr=1
do while (ierr /= 0)
   read  (unit, nml=two_stream_gray_rad_nml, iostat=io, end=10)
   ierr = check_nml_error (io, 'two_stream_gray_rad_nml')
enddo
10 call close_file (unit)

unit = open_file ('logfile.out', action='append')
if ( mpp_pe() == 0 ) then
  write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
  write (unit, nml=two_stream_gray_rad_nml)
endif
call close_file (unit)

pi         = 4. * atan(1.)
deg_to_rad = pi/180.
rad_to_deg = 180./pi

sp = trim(cld_species)

call astronomy_init

if(dt_rad_avg .le. 0) dt_rad_avg = dt_real !s if dt_rad_avg is set to a value in nml then it will be used instead of dt_real

if (do_cld_condensation) then
  if (cld_mh <= 0.0) then
    call error_mesg('two_stream_gray_rad','cld_mh must be positive (non-log metallicity).', FATAL)
  endif
  if (cld_B <= 0.0) then
    call error_mesg('two_stream_gray_rad','cld_B must be positive.', FATAL)
  endif
  call error_mesg('two_stream_gray_rad','Condensate cloud scheme active for species: '//trim(sp), NOTE)
endif

! cloud vertical extent
select case (uppercase(trim(cld_extent)))
case ('NONE')
  cld_ext = CLD_EXT_NONE
case ('THICK')
  cld_ext = CLD_EXT_THICK
case ('THIN')
  cld_ext = CLD_EXT_THIN
case default
  call error_mesg('two_stream_gray_rad','cld_extent must be "none", "thick" or "thin".', FATAL)
end select
if (cld_ext /= CLD_EXT_NONE .and. .not. do_cld_condensation) then
  call error_mesg('two_stream_gray_rad','cld_extent /= "none" requires do_cld_condensation = .true.', FATAL)
endif
if (cld_thick_ptop <= 0.0 .or. cld_thin_nscale <= 0.0) then
  call error_mesg('two_stream_gray_rad','cld_thick_ptop and cld_thin_nscale must be positive.', FATAL)
endif
if (do_cld_condensation) then
  call error_mesg('two_stream_gray_rad','Cloud vertical extent: '//trim(cld_extent), NOTE)
endif

if (do_cld_toa_albedo .and. .not. do_cld_condensation) then
  call error_mesg('two_stream_gray_rad','do_cld_toa_albedo requires do_cld_condensation = .true.', FATAL)
endif
if (cld_toa_albedo < 0.0 .or. cld_toa_albedo >= 1.0) then
  call error_mesg('two_stream_gray_rad','cld_toa_albedo must be in [0,1).', FATAL)
endif

! cloud longwave opacity: scheme-specific requirements
if (do_cld_lw_opacity) then
  if (.not. do_cld_condensation) then
    call error_mesg('two_stream_gray_rad','do_cld_lw_opacity requires do_cld_condensation = .true.', FATAL)
  endif
  if (uppercase(trim(rad_scheme)) /= 'FRIERSON') then
    call error_mesg('two_stream_gray_rad','do_cld_lw_opacity is currently only implemented for rad_scheme = "frierson".', FATAL)
  endif
  if (.not. do_normal_integration_method) then
    call error_mesg('two_stream_gray_rad','do_cld_lw_opacity currently requires do_normal_integration_method = .true.', FATAL)
  endif
endif

! cloud shortwave opacity: scheme-specific requirements
if (do_cld_sw_opacity) then
  if (.not. do_cld_condensation) then
    call error_mesg('two_stream_gray_rad','do_cld_sw_opacity requires do_cld_condensation = .true.', FATAL)
  endif
  if (uppercase(trim(rad_scheme)) /= 'FRIERSON' .and. uppercase(trim(rad_scheme)) /= 'BYRNE') then
    call error_mesg('two_stream_gray_rad','do_cld_sw_opacity is currently only implemented for rad_scheme = "frierson" or "byrne".', FATAL)
  endif
endif

! cloud optical depth: shared parameter checks
if (do_cld_lw_opacity .or. do_cld_sw_opacity) then
  if (cld_r_part <= 0.0 .or. cld_rho_part <= 0.0 .or. cld_mu_gas <= 0.0) then
    call error_mesg('two_stream_gray_rad','cld_r_part, cld_rho_part, and cld_mu_gas must be positive.', FATAL)
  endif
  if (cld_Q_lw_eff < 0.0 .or. cld_Q_sw_eff < 0.0) then
    call error_mesg('two_stream_gray_rad','cld_Q_lw_eff and cld_Q_sw_eff must be non-negative.', FATAL)
  endif
  if (cld_x_gas < 0.0 .or. cld_f_cond < 0.0 .or. cld_f_cond > 1.0) then
    call error_mesg('two_stream_gray_rad','cld_x_gas must be non-negative and cld_f_cond must be in [0,1].', FATAL)
  endif
  if (rdgas <= 0.0) then
    call error_mesg('two_stream_gray_rad','rdgas must be positive to derive the mean atmospheric molecular weight for the cloud opacity.', FATAL)
  endif
endif

if(do_read_co2)then
   call interpolator_init (co2_interp, trim(co2_file)//'.nc', lonb, latb, data_out_of_bounds=(/ZERO/))
endif


if(uppercase(trim(rad_scheme)) == 'GEEN') then
  lw_scheme = B_GEEN
  call error_mesg('two_stream_gray_rad','Using Geen (2015) radiation scheme.', NOTE)
else if(uppercase(trim(rad_scheme)) == 'FRIERSON') then
  lw_scheme = B_FRIERSON
  !call error_mesg('two_stream_gray_rad','Using "'//trim(rad_scheme)//'" radiation scheme.', NOTE)
  call error_mesg('two_stream_gray_rad','Using Frierson (2006) radiation scheme.', NOTE)
else if(uppercase(trim(rad_scheme)) == 'BYRNE') then
  lw_scheme = B_BYRNE
  call error_mesg('two_stream_gray_rad','Using Byrne & OGorman (2013) radiation scheme.', NOTE)
else if(uppercase(trim(rad_scheme)) == 'SCHNEIDER') then
  lw_scheme = B_SCHNEIDER_LIU
  call error_mesg('two_stream_gray_rad','Using Schneider & Liu (2009) radiation scheme for GIANT PLANETS.', NOTE)
else
  call error_mesg('two_stream_gray_rad','"'//trim(rad_scheme)//'"'//' is not a valid radiation scheme.', FATAL)
endif
sw_scheme = lw_scheme


if(lw_scheme == B_SCHNEIDER_LIU) then
	g_asym          = 1 - 2.*back_scatter;                 ! asymmetry factor
	gp_albedo     = ( sqrt(1. - g_asym*single_albedo) - sqrt(1. - single_albedo) )        &
           / ( sqrt(1. - g_asym*single_albedo) + sqrt(1. - single_albedo) );
	Ga_asym         = 2.*sqrt( (1. - single_albedo) * (1. - g_asym*single_albedo) );
endif

if ((lw_scheme == B_BYRNE).or.(lw_scheme == B_GEEN)) then
    if (pstd_mks/=pstd_mks_earth) call error_mesg('two_stream_gray_rad','Pstd_mks and pstd_mks_earth are not the same in the this run, but lw scheme will use pstd_mks_earth because abs coeffs in Byrne and Geen schemes are non-dimensionalized by Earth surface pressure.', NOTE)
endif



initialized = .true.

allocate (b                (ie-is+1, je-js+1, num_levels))
allocate (tdt_rad          (ie-is+1, je-js+1, num_levels))
allocate (tdt_solar        (ie-is+1, je-js+1, num_levels))

allocate (lw_dtrans        (ie-is+1, je-js+1, num_levels))
allocate (lw_tau           (ie-is+1, je-js+1, num_levels+1))
allocate (sw_tau           (ie-is+1, je-js+1, num_levels+1))
allocate (lw_up            (ie-is+1, je-js+1, num_levels+1))
allocate (lw_down          (ie-is+1, je-js+1, num_levels+1))
allocate (lw_flux          (ie-is+1, je-js+1, num_levels+1))
allocate (sw_up            (ie-is+1, je-js+1, num_levels+1))
allocate (sw_down          (ie-is+1, je-js+1, num_levels+1))
allocate (sw_flux          (ie-is+1, je-js+1, num_levels+1))
allocate (rad_flux         (ie-is+1, je-js+1, num_levels+1))

allocate (b_surf           (ie-is+1, je-js+1))
allocate (lw_tau_0         (ie-is+1, je-js+1))
allocate (sw_tau_0         (ie-is+1, je-js+1))
allocate (olr              (ie-is+1, je-js+1))
allocate (net_lw_surf      (ie-is+1, je-js+1))
allocate (toa_sw_in        (ie-is+1, je-js+1))
allocate (lat_tl           (ie-is+1, je-js+1))
allocate (lat_tl_insol     (ie-is+1, je-js+1))

allocate (insolation       (ie-is+1, je-js+1))
allocate (p2               (ie-is+1, je-js+1))
allocate (coszen           (ie-is+1, je-js+1))
allocate (fracsun          (ie-is+1, je-js+1)) !jp from astronomy.f90 : fraction of sun on surface

if (do_cld_condensation) then
  allocate (tcond_cld      (ie-is+1, je-js+1, num_levels))
  allocate (cld_flag       (ie-is+1, je-js+1, num_levels))
  allocate (cld_cond_flag  (ie-is+1, je-js+1, num_levels))
  allocate (p_cbase_cld    (ie-is+1, je-js+1))
  allocate (p_ctop_cld     (ie-is+1, je-js+1))
endif

if (do_cld_lw_opacity .or. do_cld_sw_opacity) then
  allocate (mgas_cld    (ie-is+1, je-js+1, num_levels))
endif
if (do_cld_lw_opacity) then
  allocate (tau_lw_cld  (ie-is+1, je-js+1, num_levels))
endif
if (do_cld_sw_opacity) then
  allocate (tau_sw_cld  (ie-is+1, je-js+1, num_levels))
endif

allocate (toa_albedo (ie-is+1, je-js+1))

!Allocate RG variables only if option selected.
select case(lw_scheme)
case(B_GEEN)
  allocate (b_win            (ie-is+1, je-js+1, num_levels))
  allocate (lw_dtrans_win    (ie-is+1, je-js+1, num_levels))
  allocate (lw_up_win        (ie-is+1, je-js+1, num_levels+1))
  allocate (lw_down_win      (ie-is+1, je-js+1, num_levels+1))
  allocate (lw_del_tau       (ie-is+1, je-js+1))
  allocate (lw_del_tau_win   (ie-is+1, je-js+1))

  allocate (sw_dtrans        (ie-is+1, je-js+1, num_levels))
  allocate (sw_wv            (ie-is+1, je-js+1))
  allocate (del_sol_tau      (ie-is+1, je-js+1))
  allocate (sw_tau_k         (ie-is+1, je-js+1))

case(B_BYRNE)
  allocate (lw_del_tau       (ie-is+1, je-js+1))

case(B_SCHNEIDER_LIU)
	allocate (b_surf_gp      (ie-is+1, je-js+1))
end select


!-----------------------------------------------------------------------
!------------ initialize diagnostic fields ---------------

    id_olr = &
    register_diag_field ( mod_name, 'olr', axes(1:2), Time, &
               'outgoing longwave radiation', &
               'W/m2', missing_value=missing_value               )
    id_swdn_sfc = &
    register_diag_field ( mod_name, 'swdn_sfc', axes(1:2), Time, &
               'Absorbed SW at surface', &
               'W/m2', missing_value=missing_value               )
    id_swdn_toa = &
    register_diag_field ( mod_name, 'swdn_toa', axes(1:2), Time, &
               'SW flux down at TOA', &
               'W/m2', missing_value=missing_value               )
    id_lwup_sfc = &
    register_diag_field ( mod_name, 'lwup_sfc', axes(1:2), Time, &
               'LW flux up at surface', &
               'W/m2', missing_value=missing_value               )

    id_lwdn_sfc = &
    register_diag_field ( mod_name, 'lwdn_sfc', axes(1:2), Time, &
               'LW flux down at surface', &
               'W/m2', missing_value=missing_value               )

    id_net_lw_surf = &
    register_diag_field ( mod_name, 'net_lw_surf', axes(1:2), Time, &
               'Net upward LW flux at surface', &
               'W/m2', missing_value=missing_value               )

    id_tdt_rad = &
        register_diag_field ( mod_name, 'tdt_rad', axes(1:3), Time, &
               'Temperature tendency due to radiation', &
               'K/s', missing_value=missing_value               )

    id_tdt_solar = &
        register_diag_field ( mod_name, 'tdt_solar', axes(1:3), Time, &
               'Temperature tendency due to solar radiation', &
               'K/s', missing_value=missing_value               )

    id_flux_rad = &
        register_diag_field ( mod_name, 'flux_rad', axes(half), Time, &
               'Total radiative flux (positive up)', &
               'W/m^2', missing_value=missing_value               )
    id_flux_lw = &
        register_diag_field ( mod_name, 'flux_lw', axes(half), Time, &
               'Net longwave radiative flux (positive up)', &
               'W/m^2', missing_value=missing_value               )
    id_flux_sw = &
        register_diag_field ( mod_name, 'flux_sw', axes(half), Time, &
               'Net shortwave radiative flux (positive up)', &
               'W/m^2', missing_value=missing_value               )

    id_coszen  = &
               register_diag_field ( mod_name, 'coszen', axes(1:2), Time, &
                 'cosine of zenith angle', &
                 'none', missing_value=missing_value      )
                 
    id_tau     = &
               register_diag_field ( mod_name, 'tau', axes(1:2), Time, &
                 'surface lw optical depth', &
                 'none', missing_value=missing_value      )             
                 
    id_fracsun  = &
               register_diag_field ( mod_name, 'fracsun', axes(1:2), Time, &
                 'daylight fraction of time interval', &
                 'none', missing_value=missing_value      )

    id_co2  = &
               register_diag_field ( mod_name, 'co2', Time, &
                 'co2 concentration', &
                 'ppmv', missing_value=missing_value      )

  if (lw_scheme.eq.B_GEEN) then
    id_lw_dtrans_win  = &
               register_diag_field ( mod_name, 'lw_dtrans_win', axes(1:3), Time, &
                 'LW window transmission', &
                 'none', missing_value=missing_value      )
  endif
  if (sw_scheme.eq.B_GEEN) then
    id_sw_dtrans  = &
                 register_diag_field ( mod_name, 'sw_dtrans', axes(1:3), Time, &
                   'SW transmission', &
                   'none', missing_value=missing_value      )
  endif
    id_lw_dtrans  = &
               register_diag_field ( mod_name, 'lw_dtrans', axes(1:3), Time, &
                 'LW transmission (non window)', &
                 'none', missing_value=missing_value      )
   id_nudge_flux = &
   register_diag_field ( mod_name, 'nudge_flux', axes(1:2), Time, &
               'Sea ice nudging flux', &
               'W/m2', missing_value=missing_value               )

    id_tcond_cld = &
    register_diag_field ( mod_name, 'cld_tcond', axes(1:3), Time, &
               trim(sp)//' condensation temperature', &
               'K', missing_value=missing_value               )

    id_cld_cond_flag = &
    register_diag_field ( mod_name, 'cld_cond_flag', axes(1:3), Time, &
               trim(sp)//' condensation flag before vertical-extent limit (1 where T < Tcond)', &
               'none', missing_value=missing_value            )

    id_cld_flag = &
    register_diag_field ( mod_name, 'cld_flag', axes(1:3), Time, &
               trim(sp)//' cloud flag after vertical-extent limit (used by radiation)', &
               'none', missing_value=missing_value            )

    id_p_cbase_cld = &
    register_diag_field ( mod_name, 'cld_p_cbase', axes(1:2), Time, &
               trim(sp)//' cloud base pressure (0 if cloud-free)', &
               'Pa', missing_value=missing_value              )

    id_p_ctop_cld = &
    register_diag_field ( mod_name, 'cld_p_ctop', axes(1:2), Time, &
               trim(sp)//' cloud top pressure (0 if cloud-free)', &
               'Pa', missing_value=missing_value              )

    id_toa_albedo = &
    register_diag_field ( mod_name, 'cld_toa_albedo', axes(1:2), Time, &
               'TOA albedo applied to insolation', &
               'none', missing_value=missing_value            )

    id_mgas_cld = &
    register_diag_field ( mod_name, 'cld_mgas', axes(1:3), Time, &
               'Mass of condensible gas ('//trim(sp)//' cloud) per unit area in cloudy cells', &
               'kg/m2', missing_value=missing_value            )

    id_tau_lw_cld = &
    register_diag_field ( mod_name, 'cld_tau_lw', axes(1:3), Time, &
               trim(sp)//' cloud longwave optical depth added in cloudy cells', &
               'none', missing_value=missing_value            )

    id_tau_sw_cld = &
    register_diag_field ( mod_name, 'cld_tau_sw', axes(1:3), Time, &
               trim(sp)//' cloud shortwave optical depth added in cloudy cells (per layer)', &
               'none', missing_value=missing_value            )
return
end subroutine two_stream_gray_rad_init

! ==================================================================================

subroutine two_stream_gray_rad_down (is, js, Time_diag, lat, lon, p_half, t,         &
                           net_surf_sw_down, surf_lw_down, albedo, q, const_correct, nudge_out) ! NTL ICE CHANGE 

! Begin the radiation calculation by computing downward fluxes.
! This part of the calculation does not depend on the surface temperature.

integer, intent(in)                 :: is, js
type(time_type), intent(in)         :: Time_diag
real, intent(in), dimension(:,:)    :: lat, lon, albedo
real, intent(out), dimension(:,:)   :: net_surf_sw_down
real, intent(out), dimension(:,:)   :: surf_lw_down
real, intent(in), dimension(:,:,:)  :: t, q,  p_half
real, intent(in), dimension(:,:) :: const_correct, nudge_out
integer :: i, j, k, n, dyofyr

integer :: seconds, year_in_s, days
real :: r_seconds, frac_of_day, frac_of_year, gmt, time_since_ae, rrsun, day_in_s, r_solday, r_total_seconds, r_days, r_dt_rad_avg, dt_rad_radians
logical :: used

! all for closein insolation 
real :: ci_lat_f, ci_lat_d
real, dimension(size(q,1), size(q,2)) :: ci_rho_dist, ci_J_f, ci_alpha, ci_beta, ci_phi, & 
                                         ci_t, ci_l, ci_s, ci_delta1, ci_delta2, ci_J_p

! cumulative (TOA-downward) cloud SW optical depth
real, dimension(size(q,1), size(q,2)) :: tau_sw_cld_cum

real ,dimension(size(q,1),size(q,2),size(q,3)) :: co2f

n = size(t,3)



! albedo(:,:) = albedo_value !s albedo now set in mixed_layer_init.

! =================================================================================
! SHORTWAVE RADIATION

! insolation at TOA
if (do_seasonal) then
  ! Seasonal Cycle: Use astronomical parameters to calculate insolation
  call get_time(Time_diag, seconds, days)
  call get_time(length_of_year(), year_in_s)
  r_seconds = real(seconds)
  day_in_s = length_of_day()
  frac_of_day = r_seconds / day_in_s

  if(solday .ge. 0) then
      r_solday=real(solday)
      frac_of_year = (r_solday*day_in_s) / year_in_s
  else
      r_days=real(days)
      r_total_seconds=r_seconds+(r_days*day_in_s)
      frac_of_year = r_total_seconds / year_in_s
  endif

  gmt = abs(mod(frac_of_day, 1.0)) * 2.0 * pi

  time_since_ae = modulo(frac_of_year-equinox_day, 1.0) * 2.0 * pi

  if(use_time_average_coszen) then

     r_dt_rad_avg=real(dt_rad_avg)
     dt_rad_radians = (r_dt_rad_avg/day_in_s)*2.0*pi

     call diurnal_solar(lat, lon, gmt, time_since_ae, coszen, fracsun, rrsun, dt_rad_radians)
  else
     call diurnal_solar(lat, lon, gmt, time_since_ae, coszen, fracsun, rrsun)
  end if

     insolation = solar_constant * coszen
     
  if (do_toa_albedo) then 
      p2          = (1. - 3.*sin(lat)**2)/4.
      insolation  = insolation * (0.75 + 0.15 * 2. * p2)
  endif 

     insolation = insolation 
     
else if (do_tl) then 


  if (do_closein) then 
      
      lat_tl_insol = asin(cos(lat)*cos(lon)) + pi/2 
      
      ci_lat_f = acos((radius+R_stellar) / d_stellar)
      ci_lat_d = acos((radius-R_stellar) / d_stellar) 
      
      ci_rho_dist = sqrt(d_stellar**2 + radius**2 - 2*d_stellar*radius*cos(lat_tl_insol))
      
      ! J_f
      ci_J_f = (d_stellar * cos(lat_tl_insol) - radius)*R_stellar**2 / ci_rho_dist**3 
      
      ! J_p
      ci_alpha = acos((d_stellar*cos(lat_tl_insol) - radius) / ci_rho_dist) 
      ci_beta = pi/2 - ci_alpha 
      ci_phi = asin(R_stellar / ci_rho_dist)
      ci_t = ci_rho_dist * cos(ci_phi)
      ci_l = ci_t * sin(ci_phi)
      ci_s = ci_t * sin(ci_phi) * cos(ci_alpha) 
      ci_delta1 = acos(cos(ci_phi)/cos(ci_beta)) 
      ci_delta2 = acos(tan(ci_beta)/tan(ci_phi)) 
      ci_J_p = ci_l * ci_s / pi / ci_t**2 * (pi - ci_delta2 + sin(ci_delta2)*cos(ci_delta2)) + 1 / pi * (ci_delta1 - sin(ci_delta1)*cos(ci_delta1))
      
      where (lat_tl_insol .le. ci_lat_f)
          insolation = solar_constant * (d_stellar / R_stellar)**2 * ci_J_f 
      elsewhere ((lat_tl_insol .gt. ci_lat_f) .and. (lat_tl_insol .lt. ci_lat_d)) 
          insolation = solar_constant * (d_stellar / R_stellar)**2 * ci_J_p 
      elsewhere  (lat_tl_insol .ge. ci_lat_d) 
          insolation = 0. 
      endwhere 
      
      
      
      
  else

      insolation = solar_constant * cos(lat) * cos(lon)
  
      where (insolation .le. 0. )
          insolation = 0. 
      endwhere   
  endif 
  
  
  coszen = insolation / solar_constant 
  

else if (sw_scheme==B_SCHNEIDER_LIU) then
  insolation = (solar_constant/pi)*cos(lat)
else
  ! Default: Averaged Earth insolation at all longitudes
  p2          = (1. - 3.*sin(lat)**2)/4.
  insolation  = 0.25 * solar_constant * (1.0 + del_sol * p2 + del_sw * sin(lat))
  if (do_toa_albedo) then 
      insolation  = insolation * (0.75 + 0.15 * 2. * p2)
  endif 
  insolation = insolation 
end if

! =================================================================================
! Condensate cloud flags (incl. vertical-extent limit), cloud optical depths and
! cloud-dependent TOA albedo (must precede the SW calculation). cld_flag here is
! already limited by cld_extent, so the albedo and the LW/SW opacities all see
! the limited cloud.
if (do_cld_condensation) call cld_condensation(p_half, t, Time_diag)
if (do_cld_lw_opacity .or. do_cld_sw_opacity) call cld_opacity(p_half, Time_diag)

toa_albedo = cloud_albedo   ! background value everywhere
if (do_cld_toa_albedo) then
  where (any(cld_flag > 0.5, dim=3)) toa_albedo = cld_toa_albedo
endif

select case(sw_scheme)
case(B_GEEN)
  ! RG scheme: optical depth a function of wv and co2
  sw_tau_k = 0.
  do k = 1, n
    sw_wv = sw_tau_k + 0.5194
    sw_wv = exp( 0.01887 / (sw_tau_k + 0.009522)                      &
                 + 1.603 / ( sw_wv*sw_wv ) )
    del_sol_tau(:,:) = ( 0.0596 + 0.0029 * log(carbon_conc/360.)       &
                                + sw_wv(:,:) * q(:,:,k) )             &
                     * ( p_half(:,:,k+1) - p_half(:,:,k) ) / p_half(:,:,n+1)
    sw_dtrans(:,:,k) = exp( - del_sol_tau(:,:) )
    sw_tau_k = sw_tau_k + del_sol_tau(:,:)
  end do

  ! compute downward shortwave flux
  sw_down(:,:,1) = insolation(:,:)
  do k = 1, n
     sw_down(:,:,k+1)   = sw_down(:,:,k) * sw_dtrans(:,:,k)
  end do

case(B_FRIERSON, B_BYRNE)
  ! Default: Frierson handling of SW radiation
  ! SW optical thickness
  sw_tau_0    = (1.0 - sw_diff*sin(lat)**2)*atm_abs

  ! compute optical depths for each model level
  do k = 1, n+1
    sw_tau(:,:,k) = sw_tau_0 * (p_half(:,:,k)/pstd_mks)**solar_exponent
  end do

  ! add condensate cloud shortwave optical depth: tau_sw_total = tau_sw_atm + tau_sw_cloud.
  ! sw_tau is cumulative from the model top at half levels, and half level k+1
  ! is the bottom of layer k, so the per-layer cloud optical depth is summed
  ! downward. tau_sw_cld is zero outside cloudy cells; sw_tau(:,:,1) (TOA) is unchanged.
  if (do_cld_sw_opacity) then
    tau_sw_cld_cum = 0.0
    do k = 1, n
      tau_sw_cld_cum   = tau_sw_cld_cum + tau_sw_cld(:,:,k)
      sw_tau(:,:,k+1)  = sw_tau(:,:,k+1) + tau_sw_cld_cum
    end do
  endif

  ! compute downward shortwave flux

  do k = 1, n+1
     sw_down(:,:,k)   = (1.0 - toa_albedo(:,:)) * insolation(:,:) * exp(-sw_tau(:,:,k))
  end do

case(B_SCHNEIDER_LIU)
  ! Schneider & Liu 2009 Giant planet scheme
  ! SW optical thickness

  ! compute optical depths for each model level
  do k = 1, n+1
    sw_tau(:,:,k) = sw_tau_0_gp * (p_half(:,:,k)/pstd_mks)**sw_tau_exponent_gp
  end do

  ! compute downward shortwave flux
  do k = 1, n+1
     sw_down(:,:,k)   = insolation(:,:) * (1.0-gp_albedo)* exp(- Ga_asym * sw_tau(:,:,k))
  end do

case default
 call error_mesg('two_stream_gray_rad','invalid radiation scheme',FATAL)

end select

! =================================================================================
! LONGWAVE RADIATION

! longwave source function
b = stefan*t**4

if(do_read_co2)then
  call interpolator( co2_interp, Time_diag, p_half, co2f, trim(co2_variable_name))
  carbon_conc = maxval(co2f) !Needs maxval just because co2f is a 3d array of a constant, so maxval is just a way to pick out one number
endif


select case(lw_scheme)
case(B_GEEN)
  ! split LW in 2 bands: water-vapour window and remaining = non-window
  ! ref: Ruth Geen etal, GRL 2016 (supp. information).
  lw_dtrans_win = 1.
  do k = 1, n
    lw_del_tau    = ( ir_tau_co2 + 0.2023 * log(carbon_conc/360.)                  &
                    + ir_tau_wv1*log(ir_tau_wv2*q(:,:,k) + 1) )                    &
               * ( p_half(:,:,k+1)-p_half(:,:,k) ) / pstd_mks_earth
    lw_dtrans(:,:,k) = exp( - lw_del_tau )
    lw_del_tau_win   = ( ir_tau_co2_win + 0.0954 * log(carbon_conc/360.)           &
                                     + ir_tau_wv_win1*q(:,:,k)                 &
                                     + ir_tau_wv_win2*q(:,:,k)*q(:,:,k) )      &
                  * ( p_half(:,:,k+1)-p_half(:,:,k) ) / pstd_mks_earth
    lw_dtrans_win(:,:,k) = exp( - lw_del_tau_win )
  end do

  ! also compute downward longwave flux for window
  ! Allocate a fraction of the longwave spectrum as window radiation
  b_win = window*b
  b = (1.0 - window)*b
  lw_down_win(:,:,1) = 0.0
  lw_down(:,:,1) = 0.0
  do k = 1,n
    lw_down(:,:,k+1) = lw_down(:,:,k)*lw_dtrans(:,:,k) + b(:,:,k)*(1. - lw_dtrans(:,:,k))
    lw_down_win(:,:,k+1) = lw_down_win(:,:,k)*lw_dtrans_win(:,:,k)              &
                      + b_win(:,:,k)*(1.0 - lw_dtrans_win(:,:,k))
  end do
  lw_down = lw_down + lw_down_win

case(B_BYRNE)
  ! dtau/ds = a*mu + b*q
  ! ref: Byrne, M. P. & O'Gorman, P. A.
  !      Land–ocean warming contrast over a wide range of climates:
  !      Convective quasi-equilibrium theory and idealized simulations.
  !      J. Climate 26, 4000–4106 (2013).

  do k = 1, n
    lw_del_tau    = (bog_a*bog_mu + 0.17 * log(carbon_conc/360.)  + bog_b*q(:,:,k)) * (( p_half(:,:,k+1)-p_half(:,:,k) ) / pstd_mks_earth)
    lw_dtrans(:,:,k) = exp( - lw_del_tau )

  end do

  ! compute downward longwave flux by integrating downward
  lw_down(:,:,1)      = 0.
  do k = 1, n
     lw_down(:,:,k+1) = lw_down(:,:,k)*lw_dtrans(:,:,k) + b(:,:,k)*(1. - lw_dtrans(:,:,k))
  end do

case(B_FRIERSON)
  ! longwave optical thickness function of latitude and pressure
  if (do_tl) then 
      lat_tl = pi/2 - asin(cos(lat)*cos(lon))
      lw_tau_0 = ir_tau_eq + (ir_tau_pole - ir_tau_eq)*sin(lat_tl)**2
      where (lat_tl .ge. (pi/2)) 
          lw_tau_0 = ir_tau_pole
      endwhere 
  else 
      lw_tau_0 = ir_tau_eq + (ir_tau_pole - ir_tau_eq)*sin(lat)**2
  endif 
  lw_tau_0 = lw_tau_0 * odp ! scale by optical depth parameter - default 1

  ! compute optical depths for each model level
  do k = 1, n+1
  lw_tau(:,:,k) = lw_tau_0 * ( linear_tau * p_half(:,:,k)/pstd_mks     &
       + (1.0 - linear_tau) * (p_half(:,:,k)/pstd_mks)**wv_exponent )
  end do

  ! longwave differential transmissivity
  do k = 1, n
     lw_dtrans(:,:,k) = exp( -(lw_tau(:,:,k+1) - lw_tau(:,:,k)) )
  end do

  ! add condensate cloud longwave optical depth in cloudy cells: tau_total = tau_atm + tau_cloud,
  ! i.e. lw_dtrans_total = lw_dtrans_atm * exp(-tau_cloud)
  if (do_cld_lw_opacity) then
    do k = 1, n
      where (cld_flag(:,:,k) > 0.5)
        lw_dtrans(:,:,k) = lw_dtrans(:,:,k) * exp( -tau_lw_cld(:,:,k) )
      endwhere
    end do
  endif

  ! compute downward longwave flux by integrating downward
  lw_down(:,:,1)      = 0.
  if (do_normal_integration_method) then
      do k = 1, n
         lw_down(:,:,k+1) = lw_down(:,:,k)*lw_dtrans(:,:,k) + b(:,:,k)*(1. - lw_dtrans(:,:,k))
      end do
         
  else
      do k = 1, n
         lw_down(:,:,k+1) = 2.*b(:,:,k) * (lw_tau(:,:,k+1) - lw_tau(:,:,k))/(2.+ (lw_tau(:,:,k+1) - lw_tau(:,:,k))) + &
                            lw_down(:,:,k) * (2.- (lw_tau(:,:,k+1) - lw_tau(:,:,k)))/(2.+ (lw_tau(:,:,k+1) - lw_tau(:,:,k)))
      end do
  endif

case(B_SCHNEIDER_LIU)

  ! compute optical depths for each model level
  do k = 1, n+1
  lw_tau(:,:,k) = lw_tau_0_gp * (p_half(:,:,k)/pstd_mks)**lw_tau_exponent_gp
  end do

  ! longwave differential transmissivity
  do k = 1, n
     lw_dtrans(:,:,k) = exp( -(lw_tau(:,:,k+1) - lw_tau(:,:,k)) )
  end do

  ! compute downward longwave flux by integrating downward
  lw_down(:,:,1)      = 0.
  do k = 1, n
     lw_down(:,:,k+1) = lw_down(:,:,k)*lw_dtrans(:,:,k) + b(:,:,k)*(1. - lw_dtrans(:,:,k))
  end do


case default
 call error_mesg('two_stream_gray_rad','invalid radiation scheme',FATAL)

end select

! =================================================================================
surf_lw_down     = lw_down(:, :, n+1) 
toa_sw_in        = sw_down(:, :, 1)
net_surf_sw_down = sw_down(:, :, n+1) * (1. - albedo) 
! =================================================================================

if(lw_scheme.eq.B_SCHNEIDER_LIU) then
	b_surf_gp=surf_lw_down(:,:)+net_surf_sw_down(:,:)
endif


!------- downward lw flux surface -------
if ( id_lwdn_sfc > 0 ) then
   used = send_data ( id_lwdn_sfc, surf_lw_down, Time_diag)
endif

!surf_lw_down = surf_lw_down + const_correct + nudge_out ! NTL UPDATE NUDGING 

!------- incoming sw flux toa -------
if ( id_swdn_toa > 0 ) then
   used = send_data ( id_swdn_toa, toa_sw_in, Time_diag)
endif

! NTL NUDGING
if( id_nudge_flux > 0) then 
   used = send_data ( id_nudge_flux, const_correct+nudge_out, Time_diag)
endif 


!------- downward sw flux surface -------
if ( id_swdn_sfc > 0 ) then
   used = send_data ( id_swdn_sfc, net_surf_sw_down, Time_diag)
endif

!------- cosine of zenith angle ------------
if ( id_coszen > 0 ) then
   used = send_data ( id_coszen, coszen, Time_diag)
endif

if (id_tau > 0) then 
   used = send_data ( id_tau, lw_tau_0, Time_diag)
endif

!------- carbon dioxide concentration ------------
if ( id_co2 > 0 ) then
   used = send_data ( id_co2, carbon_conc, Time_diag)
endif

!------- TOA albedo applied to insolation ------------
if ( id_toa_albedo > 0 ) then
   used = send_data ( id_toa_albedo, toa_albedo, Time_diag)
endif

return
end subroutine two_stream_gray_rad_down

! ==================================================================================

subroutine two_stream_gray_rad_up (is, js, Time_diag, lat, p_half, t_surf, t, tdt, albedo)

! Now complete the radiation calculation by computing the upward and net fluxes.

integer, intent(in)                 :: is, js
type(time_type), intent(in)         :: Time_diag
real, intent(in) , dimension(:,:)   :: lat, albedo
real, intent(in) , dimension(:,:)   :: t_surf
real, intent(in) , dimension(:,:,:) :: t, p_half
real, intent(inout), dimension(:,:,:) :: tdt


integer :: i, j, k, n

logical :: used

n = size(t,3)
b_surf            = stefan*t_surf**4

select case(lw_scheme)
case(B_GEEN)
  ! integrate upward, including window contribution
  lw_up(:,:,n+1)     = b_surf*(1-window)
  lw_up_win(:,:,n+1) = b_surf*window

  do k = n, 1, -1
     lw_up(:,:,k)   = lw_up(:,:,k+1)*lw_dtrans(:,:,k) + b(:,:,k)*(1.0 - lw_dtrans(:,:,k))
     lw_up_win(:,:,k)   = lw_up_win(:,:,k+1)*lw_dtrans_win(:,:,k) + b_win(:,:,k)*(1.0 - lw_dtrans_win(:,:,k))
  end do
  lw_up = lw_up + lw_up_win

case(B_FRIERSON, B_BYRNE)
  ! compute upward longwave flux by integrating upward
  lw_up(:,:,n+1)    = b_surf
  if (do_normal_integration_method) then
      do k = n, 1, -1
         lw_up(:,:,k)   = lw_up(:,:,k+1)*lw_dtrans(:,:,k) + b(:,:,k)*(1.0 - lw_dtrans(:,:,k))
      end do
         
  else
      do k = n, 1, -1
        lw_up(:,:,k) = 2.*b(:,:,k) * -1.*(lw_tau(:,:,k) - lw_tau(:,:,k+1))/ (2.- (lw_tau(:,:,k) - lw_tau(:,:,k+1))) &
                        + lw_up(:,:,k+1)* (2. + (lw_tau(:,:,k) - lw_tau(:,:,k+1)))/(2. - (lw_tau(:,:,k) - lw_tau(:,:,k+1)))      
      end do                        
  endif

case(B_SCHNEIDER_LIU)
  ! compute upward longwave flux by integrating upward

  lw_up(:,:,n+1)    = b_surf_gp
  do k = n, 1, -1
     lw_up(:,:,k)   = lw_up(:,:,k+1)*lw_dtrans(:,:,k) + b(:,:,k)*(1.0 - lw_dtrans(:,:,k))
  end do

case default
 call error_mesg('two_stream_gray_rad','invalid radiation scheme',FATAL)

end select

! compute upward shortwave flux (here taken to be constant)
do k = 1, n+1
   sw_up(:,:,k)   = albedo(:,:) * sw_down(:,:,n+1)
end do

! net fluxes (positive up)
lw_flux  = lw_up - lw_down
sw_flux  = sw_up - sw_down
rad_flux = lw_flux + sw_flux

do k = 1, n
   tdt_rad(:,:,k)   = diabatic_acce * ( rad_flux(:,:,k+1) - rad_flux(:,:,k) )  &
        * grav/( cp_air*(p_half(:,:,k+1) - p_half(:,:,k)) )

   tdt_solar(:,:,k) = ( sw_flux(:,:,k+1) - sw_flux(:,:,k) )  &
        * grav/( cp_air*(p_half(:,:,k+1) - p_half(:,:,k)) )

   tdt(:,:,k) = tdt(:,:,k) + tdt_rad(:,:,k)
end do

olr         = lw_up(:,:,1)
net_lw_surf = lw_flux(:, :, n+1)

!------- outgoing lw flux toa (olr) -------
if ( id_olr > 0 ) then
   used = send_data ( id_olr, olr, Time_diag)
endif
!------- upward lw flux surface -------
if ( id_lwup_sfc > 0 ) then
   used = send_data ( id_lwup_sfc, b_surf, Time_diag)
endif
!------- net upward lw flux surface -------
if ( id_net_lw_surf > 0 ) then
   used = send_data ( id_net_lw_surf, net_lw_surf, Time_diag)
endif
!------- temperature tendency due to radiation ------------
if ( id_tdt_rad > 0 ) then
   used = send_data ( id_tdt_rad, tdt_rad, Time_diag)
endif
!------- temperature tendency due to solar radiation ------------
if ( id_tdt_solar > 0 ) then
   used = send_data ( id_tdt_solar, tdt_solar, Time_diag)
endif
!------- total radiative flux (at half levels) -----------
if ( id_flux_rad > 0 ) then
   used = send_data ( id_flux_rad, rad_flux, Time_diag)
endif
!------- longwave radiative flux (at half levels) --------
if ( id_flux_lw > 0 ) then
   used = send_data ( id_flux_lw, lw_flux, Time_diag)
endif
if ( id_flux_sw > 0 ) then
   used = send_data ( id_flux_sw, sw_flux, Time_diag)
endif
!------- radiative transmissivities (at half levels) --------
if ( id_lw_dtrans > 0 ) then
   used = send_data ( id_lw_dtrans, lw_dtrans, Time_diag)
endif
if ( id_lw_dtrans_win > 0 ) then
   used = send_data ( id_lw_dtrans_win, lw_dtrans_win, Time_diag)
endif
if ( id_sw_dtrans > 0 ) then
   used = send_data ( id_sw_dtrans, sw_dtrans, Time_diag)
endif

return
end subroutine two_stream_gray_rad_up

! ==================================================================================

subroutine cld_condensation(p_half, t, Time_diag)

! 1) Computes the condensation temperature of the cloud species at each full
!    level and flags cells where T < Tcond (cld_cond_flag), using
!
!      Tcond = cld_B / (cld_A - cld_p_coeff*log10(P_bar) - cld_mh_coeff*log10(cld_mh))
!
!    With the defaults (cld_A = 13.43, cld_B = 28665, cld_p_coeff = cld_mh_coeff = 1)
!    this is the MgSiO3 curve  log10(P_bar) = 13.43 - 28665/T - log10(mh).
!
! 2) Applies the vertical-extent limit to give cld_flag, which is what the TOA
!    albedo and LW/SW opacities use. Cloud base = deepest level with T < Tcond
!    (k = 1 is the model top). Cells with full-level pressure below p_lim are
!    cleared:
!      'none' : p_lim = 0 (no ceiling)
!      'thick': p_lim = cld_thick_ptop
!      'thin' : p_lim = p_base * exp(-cld_thin_nscale)
!    Cells between the base and p_lim are cloudy only where T < Tcond.

real, intent(in), dimension(:,:,:) :: p_half, t
type(time_type), intent(in)        :: Time_diag

integer :: i, j, k, n, kbase, ktop
real    :: log_mh, p_lim
real, dimension(size(t,1), size(t,2))            :: p_bar, denom
real, dimension(size(t,1), size(t,2), size(t,3)) :: p_full
logical :: used

n      = size(t,3)
log_mh = log10(cld_mh)

! ---- 1) raw condensation flag ----
do k = 1, n
  ! full-level pressure = mean of bounding half levels (Pa); converted to bar for Tcond
  p_full(:,:,k) = 0.5 * (p_half(:,:,k) + p_half(:,:,k+1))
  p_bar = p_full(:,:,k) / pa_per_bar
  denom = cld_A - cld_p_coeff*log10(p_bar) - cld_mh_coeff*log_mh

  if (any(denom <= 0.0)) then
    call error_mesg('two_stream_gray_rad', &
         trim(cld_species)//' condensation: pressure and metallicity give non-positive Tcond denominator.', FATAL)
  endif

  tcond_cld(:,:,k) = cld_B / denom

  where (t(:,:,k) < tcond_cld(:,:,k))
    cld_cond_flag(:,:,k) = 1.0
  elsewhere
    cld_cond_flag(:,:,k) = 0.0
  endwhere
end do

! ---- 2) vertical extent ----
cld_flag    = cld_cond_flag
p_cbase_cld = 0.0
p_ctop_cld  = 0.0

do j = 1, size(t,2)
  do i = 1, size(t,1)

    ! cloud base: deepest (highest-pressure) condensing level
    kbase = 0
    do k = n, 1, -1
      if (cld_cond_flag(i,j,k) > 0.5) then
        kbase = k
        exit
      endif
    end do
    if (kbase == 0) cycle   ! no condensation in this column

    select case (cld_ext)
    case (CLD_EXT_THICK)
      p_lim = cld_thick_ptop
    case (CLD_EXT_THIN)
      p_lim = p_full(i,j,kbase) * exp(-cld_thin_nscale)
    case default
      p_lim = 0.0             ! 'none': no ceiling
    end select

    where (p_full(i,j,:) < p_lim) cld_flag(i,j,:) = 0.0

    ! effective cloud top after truncation
    ktop = 0
    do k = 1, kbase
      if (cld_flag(i,j,k) > 0.5) then
        ktop = k
        exit
      endif
    end do
    if (ktop == 0) cycle   ! e.g. 'thick' with base already above cld_thick_ptop

    p_cbase_cld(i,j) = p_full(i,j,kbase)
    p_ctop_cld(i,j)  = p_full(i,j,ktop)
  end do
end do

if ( id_tcond_cld > 0 ) then
   used = send_data ( id_tcond_cld, tcond_cld, Time_diag)
endif
if ( id_cld_cond_flag > 0 ) then
   used = send_data ( id_cld_cond_flag, cld_cond_flag, Time_diag)
endif
if ( id_cld_flag > 0 ) then
   used = send_data ( id_cld_flag, cld_flag, Time_diag)
endif
if ( id_p_cbase_cld > 0 ) then
   used = send_data ( id_p_cbase_cld, p_cbase_cld, Time_diag)
endif
if ( id_p_ctop_cld > 0 ) then
   used = send_data ( id_p_ctop_cld, p_ctop_cld, Time_diag)
endif

end subroutine cld_condensation

! ==================================================================================

subroutine cld_opacity(p_half, Time_diag)

! Computes the condensate cloud optical depths added in cloudy cells
! (cld_flag, i.e. T < Tcond and within the cld_extent limit), following
! Roman & Rauscher (2019), ApJ 872, 1, Eqs. 1-2:
!
!   m_g   = (dP/g) * x_gas * (mu_gas / mu_atm)         [Eq. 1, mass of condensible
!                                                       gas per unit area in the layer]
!   tau_a = 3 * m_g * Q_eff * f / (4 * r * rho)        [Eq. 2, aerosol optical depth]
!
! evaluated separately for the longwave (Q_eff = cld_Q_lw_eff -> tau_lw_cld)
! and the shortwave (Q_eff = cld_Q_sw_eff -> tau_sw_cld). Both are per-layer
! optical depths at full levels.
!
! mu_gas (cld_mu_gas) is the condensible gas species and mu_atm is the bulk
! atmosphere, derived from the model's own specific gas constant via
! mu_atm = universal_gas_const / rdgas. mgas_cld (and hence both optical
! depths) is zero in non-cloudy cells.

real, intent(in), dimension(:,:,:) :: p_half
type(time_type), intent(in)        :: Time_diag

integer :: k, n
real    :: mu_atm, geom_fac
real, dimension(size(p_half,1), size(p_half,2)) :: delta_p
logical :: used

n        = size(p_half,3) - 1
mu_atm   = universal_gas_const / rdgas   ! mean molecular weight of the bulk atmosphere
geom_fac = 3.0 * cld_f_cond / (4.0 * cld_r_part * cld_rho_part)   ! m^2/kg per unit Q_eff

do k = 1, n
  delta_p = p_half(:,:,k+1) - p_half(:,:,k)

  where (cld_flag(:,:,k) > 0.5)
    mgas_cld(:,:,k) = (delta_p/grav) * cld_x_gas * (cld_mu_gas/mu_atm)
  elsewhere
    mgas_cld(:,:,k) = 0.0
  endwhere
end do

if (do_cld_lw_opacity) tau_lw_cld = cld_Q_lw_eff * geom_fac * mgas_cld
if (do_cld_sw_opacity) tau_sw_cld = cld_Q_sw_eff * geom_fac * mgas_cld

if ( id_mgas_cld > 0 ) then
   used = send_data ( id_mgas_cld, mgas_cld, Time_diag)
endif
if ( do_cld_lw_opacity .and. id_tau_lw_cld > 0 ) then
   used = send_data ( id_tau_lw_cld, tau_lw_cld, Time_diag)
endif
if ( do_cld_sw_opacity .and. id_tau_sw_cld > 0 ) then
   used = send_data ( id_tau_sw_cld, tau_sw_cld, Time_diag)
endif

end subroutine cld_opacity

! ==================================================================================


subroutine two_stream_gray_rad_end

deallocate (b, tdt_rad, tdt_solar)
deallocate (lw_up, lw_down, lw_flux, sw_up, sw_down, sw_flux, rad_flux)
deallocate (b_surf, olr, net_lw_surf, toa_sw_in, lw_tau_0, sw_tau_0)
deallocate (lw_dtrans, lw_tau, sw_tau)
deallocate (insolation, p2) !s albedo

! deallocate RG variables
if (lw_scheme.eq.B_GEEN) then
  deallocate (b_win, lw_dtrans_win, lw_up_win, lw_down_win, lw_del_tau, lw_del_tau_win)
endif
if (sw_scheme.eq.B_GEEN) then
  deallocate (sw_dtrans, sw_wv, del_sol_tau, sw_tau_k)
endif
if (lw_scheme.eq.B_BYRNE) then
  deallocate (lw_del_tau)
endif
if(lw_scheme.eq.B_SCHNEIDER_LIU) then
	deallocate (b_surf_gp)
endif

if (do_cld_condensation) deallocate (tcond_cld, cld_flag, cld_cond_flag, p_cbase_cld, p_ctop_cld)
if (do_cld_lw_opacity .or. do_cld_sw_opacity) deallocate (mgas_cld)
if (do_cld_lw_opacity) deallocate (tau_lw_cld)
if (do_cld_sw_opacity) deallocate (tau_sw_cld)

deallocate (toa_albedo)

if(do_read_co2)call interpolator_end(co2_interp)

end subroutine two_stream_gray_rad_end

! ==================================================================================

end module two_stream_gray_rad_mod