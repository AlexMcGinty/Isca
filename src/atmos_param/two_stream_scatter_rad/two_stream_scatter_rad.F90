!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Two stream grey radiative transfer code with scattering solved    !!
!! following Pierrehumbert (2010), Principles of Planerary Climate.  !!
!! Chapter 5, Page 353.                                              !!
!!                                                                   !!
!! Forumulates radiative transfer equations as matrix problem and    !!
!! solves via LU decomposition.                                      !!
!!                                                                   !!
!! Code that carries out LU decomposition and solves for a generic   !!
!! system Ax = B was downloaded at the following link:               !!
!! http://jean-pierre.moreau.pagesperso-orange.fr/f_matrices.html    !!
!!                                                                   !!
!! Please direct any questions to Neil Lewis, Univ. Oxford, at       !! 
!! neil.lewis <at> physics.ox.ac.uk                                  !!
!!                                                                   !!
!! Distributed under GNU general public license, available here:     !!
!! http://www.gnu.org/licenses/gpl.txt                               !!
!!                                                                   !!
!! Modifications:                                                    !!
!!  - Tidally locked insolation (do_tl) and close-in finite stellar  !!
!!    disc (do_closein), ported from two_stream_gray_rad.            !!
!!  - Direct beam attenuated exactly, exp(-|dtau|/mu), instead of    !!
!!    the Pade factor (which goes negative when |dtau| > 2 mu).      !!
!!  - Direct-beam scattering source in the solver uses the actual    !!
!!    beam loss per layer, so scattered energy is conserved exactly. !!
!!  - IAMAX initialises its output index.                            !!
!!  - v1.3: solver builds the band matrix directly (no N x N arrays),!!
!!    optional radiation timestep dt_rad, timing clock.              !!
!!  - v1.4: coszen_fixed (fixed direct-beam mu, default branch).     !!
!!  - v1.5: book-consistent normalisation and asymmetric scattering: !!
!!    * half_tau_diffuse = .false.: the diffuse equations use the    !!
!!      full layer optical depth, so gamma, gammaprime and gamma+-   !!
!!      follow Pierrehumbert Eq. 5.30 exactly. .true. reproduces     !!
!!      v1.4 (diffuse equations saw dtau/2, gamma+- did not).        !!
!!    * two_stream_closure: 'hemi', 'quadrature', 'eddington' set    !!
!!      gamma, gammaprime, ghat_factor from Table 5.3 / Eq. 5.21     !!
!!      (and force half_tau_diffuse = .false.); 'custom' (default)   !!
!!      uses the namelist values.                                    !!
!!        closure      gamma       gammaprime  ghat_factor           !!
!!        hemi         1           1           3/2                   !!
!!        quadrature   sqrt(3)/2   sqrt(3)/2   1                     !!
!!        eddington    3/4         1           1                     !!
!!      ghat = ghat_factor * g_asym, g_asym = Mie asymmetry factor.  !!
!!      Eddington couples the streams at omega = 0 (gamma2 = -1/4),  !!
!!      so it does not reduce to the no-scattering scheme.           !!
!!    * do_delta_scaling: delta-Eddington (Joseph et al. 1976),      !!
!!      f = g^2, tau' = (1 - w f) tau, w' = (1 - f) w / (1 - w f),   !!
!!      g' = g / (1 + g); the direct beam is attenuated with tau'.   !!
!!    * one-time NOTE if gamma+ < 0 (unphysical upward split).       !!
!!  - v1.6: separate SW and LW closures:                            !!
!!    sw_two_stream_closure, lw_two_stream_closure ('' = inherit     !!
!!    two_stream_closure). Recommended: SW 'eddington' with          !!
!!    do_delta_scaling (forward-peaked cloud scattering), LW 'hemi'  !!
!!    (thermal IR; keeps the validated clear-sky LW).                !!
!!    One-time NOTE if ghat >= 1 in a band (two-stream breaks down:  !!
!!    1 - ghat*omega <= 0; use delta scaling or another closure).    !!
!!    LW opacity mapping: with half_tau_diffuse = .false. and hemi,  !!
!!    the diffuse coefficient is 2 (diffusivity 2), so a Frierson-   !!
!!    style optical depth tau_F (diffusivity folded in) maps to      !!
!!    lw_abs_a = (tau_F / 2) * g / p_ref.                            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module two_stream_scatter_rad_mod

! ==================================================================================
! ==================================================================================

   use mpp_mod,               only: mpp_clock_id, mpp_clock_begin, mpp_clock_end

   use fms_mod,               only: open_file, check_nml_error, &
                                    mpp_pe, close_file, error_mesg, &
                                    NOTE, FATAL,  uppercase

   use constants_mod,         only: stefan, cp_air, grav, pstd_mks, pstd_mks_earth, &
                                    seconds_per_sol, orbital_period, radius

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

character(len=128) :: version='$Id: two_stream_scatter_rad.F90,v 1.6 2026/11/12$'
character(len=128) :: tag='Two-stream semi-gray atmosphere with scattering'


! public interfaces

public :: two_stream_scatter_rad_init, two_stream_scatter_rad, two_stream_scatter_rad_end
!==================================================================================


! module variables
logical :: initialized     = .false.

! astronomical parameters
real    :: solar_constant  = 1360.0
real    :: del_sol         = 1.4
real    :: del_sw          = 0.0
logical :: do_seasonal     = .false.
integer :: solday          = -10  !s Day of year to run perpetually if do_seasonal=True and solday>0
real    :: equinox_day     = 0.75 !s Fraction of year [0,1] where NH autumn equinox occurs (only really useful if calendar has defined months).
logical :: use_time_average_coszen = .false. !s if .true., then time-averaging is done on coszen so that insolation doesn't depend on timestep
real    :: dt_rad_avg     = -1
real :: diabatic_acce = 1.0 ! artefact from gp scheme in two_stream_grey_rad 
! fixed cosine of zenith angle for the direct beam (default insolation branch only).
! <= 0 (default): coszen = insolation / solar_constant (= 0.25 for del_sol = del_sw = 0).
! Lets a column use a physical mu (e.g. 0.5 dayside mean, 1 substellar) without
! changing the incident flux, which is still set by solar_constant/4.
real    :: coszen_fixed = -1.0

! radiation timestep (s). Radiative fluxes and heating are recomputed every
! dt_rad seconds and held fixed in between. <= 0 (default): every timestep.
! Must be a multiple of the model timestep.
real    :: dt_rad = 0.0
integer :: rad_step_interval = 1, rad_step_counter = 0
integer :: id_clock_rad

! tidally locked insolation (point source: substellar point at lat = 0, lon = 0;
! do_closein: substellar point at lon = 180 deg)
logical :: do_tl      = .false.
logical :: do_closein = .false.      ! finite stellar disc (requires do_tl)
real    :: R_stellar  = 476700000.   ! m
real    :: d_stellar  = 1117512000.0 ! m, star-planet separation

! co2 
real :: carbon_conc = 360.0 ! ppmv

character(len=32) :: sw_optical_depth = 'generic'
integer :: sw_flag 
integer, parameter :: sw_TRANSPARENT = 1
integer, parameter :: sw_GENERIC = 2

! parameters for lw optical depths etc (defaults to BOG... i.e. no scattering)
real :: lw_sca_a = 0.0
real :: lw_sca_b = 0.0
real :: lw_sca_c = 0.0 
real :: lw_abs_a = 0.1627 * 9.80 / 101325.0 ! dtau = a * dsigma   --->   dtau = 1/g * c * dp   where c = ga/p_s 
real :: lw_abs_b = 1997.9 * 9.80 / 101325.0 
real :: lw_abs_c = 0.17   * 9.80 / 101325.0

! parameters for sw optical depths (defaults to BOG... i.e. no scattering or absorption)
real :: sw_sca_a = 0.0
real :: sw_sca_b = 0.0
real :: sw_sca_c = 0.0
real :: sw_abs_a = 0.0
real :: sw_abs_b = 0.0 

real :: gamma = 1.0 ! associated with closure for integrating over all angles 
real :: gammaprime = 1.0 ! " " " (see discussion in Pierrehumbert, 2010, Chapter 5)
real :: g_asym = 0.0 ! asymmetry factor (Mie g, Pierrehumbert Eq. 5.17); used for SW and LW

! v1.5: closure, normalisation and asymmetric scattering
character(len=32) :: two_stream_closure = 'custom' ! 'custom' | 'hemi' | 'quadrature' | 'eddington'
real    :: ghat_factor      = 1.5      ! ghat = ghat_factor * g_asym (3/2 hemi, 1 Eddington/quadrature)
logical :: half_tau_diffuse = .true.   ! .true.: v1.4 behaviour (diffuse eqs see dtau/2); .false.: Eq. 5.30
logical :: do_delta_scaling = .false.  ! delta-Eddington scaling of tau, omega, g
real    :: tau_fac          = 0.5      ! set in init: 0.5 if half_tau_diffuse, else 1
logical :: warned_gammaplus = .false.

! v1.6: per-band closures ('' = use two_stream_closure). Resolved in init into the
! band coefficients below; 'custom' takes gamma, gammaprime, ghat_factor from the namelist.
character(len=32) :: sw_two_stream_closure = ''
character(len=32) :: lw_two_stream_closure = ''
real    :: sw_gamma, sw_gammaprime, sw_ghat_factor
real    :: lw_gamma, lw_gammaprime, lw_ghat_factor
logical :: warned_ghat_sw = .false., warned_ghat_lw = .false.


real, allocatable, dimension(:,:)   :: insolation, p2,  sw_tau_0 !s albedo now defined in mixed_layer_init
real, allocatable, dimension(:,:)   :: pi_B_surf
real, allocatable, dimension(:,:,:) :: pi_B, tdt_rad, tdt_solar
real, allocatable, dimension(:,:,:) :: lw_up, lw_down, lw_flux, &
                                       sw_up, sw_down, sw_down_direct, sw_flux, &
                                       rad_flux
real, allocatable, dimension(:,:,:) :: lw_tau, sw_tau, lw_del_tau, sw_del_tau
real, allocatable, dimension(:,:,:) :: lw_ss_albedo, sw_ss_albedo
real, allocatable, dimension(:,:,:) :: lw_abs_coeff, lw_scatter_coeff, &
                                       sw_abs_coeff, sw_scatter_coeff
real, allocatable, dimension(:,:,:) :: sw_gammab, sw_gammaone, sw_gammatwo, &
                                       sw_gammaminus, sw_gammaplus 
real, allocatable, dimension(:,:,:) :: lw_gammab, lw_gammaone, lw_gammatwo, &
                                       lw_gammaminus, lw_gammaplus 
real, allocatable, dimension(:,:)   :: olr, net_lw_surf, toa_sw_in, &
                                       toa_sw_out, coszen, fracsun


real, save :: pi, deg_to_rad , rad_to_deg

!extras for reading in co2 concentration
logical                             :: do_read_co2=.false.
type(interpolate_type),save         :: co2_interp           ! use external file for co2
character(len=256)                  :: co2_file='co2'       !  file name of co2 file to read
character(len=256)                  :: co2_variable_name='co2'       !  file name of co2 file to read


namelist/two_stream_scatter_rad_nml/ solar_constant, del_sol, del_sw, &
           do_seasonal, solday, equinox_day,  &
           use_time_average_coszen, dt_rad_avg,&
           diabatic_acce,& !Schneider Liu values 
           lw_abs_a, lw_abs_b, lw_abs_c, &
           lw_sca_a, lw_sca_b, lw_sca_c, &
           sw_abs_a, sw_abs_b, &
           sw_sca_a, sw_sca_b, sw_sca_c, &
           gamma, gammaprime, g_asym, &
           do_read_co2, co2_file, co2_variable_name, carbon_conc, &
           sw_optical_depth, &
           do_tl, do_closein, R_stellar, d_stellar, dt_rad, coszen_fixed, &
           two_stream_closure, ghat_factor, half_tau_diffuse, do_delta_scaling, &
           sw_two_stream_closure, lw_two_stream_closure


!==================================================================================
!-------------------- diagnostics fields -------------------------------

integer :: id_olr, id_swdn_sfc, id_swdn_toa, id_swup_toa, &
           id_net_lw_surf, id_lwdn_sfc, &
           id_lwup_sfc, id_tdt_rad, id_tdt_solar, id_flux_rad, id_flux_lw, &
           id_flux_sw, id_coszen, id_fracsun, &
           id_flux_lw_up, id_flux_lw_down, id_flux_sw_up, id_flux_sw_down, &
           id_flux_sw_down_direct, id_co2

character(len=18), parameter :: mod_name = 'two_stream_scatter'

real :: missing_value = -999.


contains


! ==================================================================================
! ==================================================================================


subroutine two_stream_scatter_rad_init(is, ie, js, je, num_levels, axes, Time, lonb, latb, dt_real)

!-------------------------------------------------------------------------------------
integer, intent(in), dimension(4) :: axes
type(time_type), intent(in)       :: Time
integer, intent(in)               :: is, ie, js, je, num_levels
real ,dimension(:,:),intent(in),optional :: lonb,latb !s Changed to 2d arrays as 2013 interpolator expects this.
real, intent(in)                  :: dt_real !s atmospheric timestep, used for radiation averaging
!-------------------------------------------------------------------------------------
integer, dimension(3) :: half = (/1,2,4/)
integer :: ierr, io, unit
logical :: sw_preset, lw_preset
!-----------------------------------------------------------------------------------------
! read namelist and copy to logfile

unit = open_file ('input.nml', action='read')
ierr=1
do while (ierr /= 0)
   read  (unit, nml=two_stream_scatter_rad_nml, iostat=io, end=10)
   ierr = check_nml_error (io, 'two_stream_scatter_rad_nml')
enddo
10 call close_file (unit)

unit = open_file ('logfile.out', action='append')
if ( mpp_pe() == 0 ) then
  write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
  write (unit, nml=two_stream_scatter_rad_nml)
endif
call close_file (unit)

pi         = 4. * atan(1.)
deg_to_rad = pi/180.
rad_to_deg = 180./pi

call astronomy_init

if (do_closein .and. .not. do_tl) then
  call error_mesg('two_stream_scatter_rad','do_closein requires do_tl = .true.', FATAL)
endif
if (do_tl .and. do_seasonal) then
  call error_mesg('two_stream_scatter_rad','do_tl and do_seasonal are mutually exclusive.', FATAL)
endif
if (coszen_fixed > 1.0) then
  call error_mesg('two_stream_scatter_rad','coszen_fixed must be <= 1.', FATAL)
endif
if (coszen_fixed > 0.0 .and. (do_tl .or. do_seasonal)) then
  call error_mesg('two_stream_scatter_rad', &
       'coszen_fixed only applies to the default insolation branch; ignored with do_tl/do_seasonal.', NOTE)
endif
if (len_trim(sw_two_stream_closure) == 0) sw_two_stream_closure = two_stream_closure
if (len_trim(lw_two_stream_closure) == 0) lw_two_stream_closure = two_stream_closure
call set_closure(sw_two_stream_closure, sw_gamma, sw_gammaprime, sw_ghat_factor, sw_preset)
call set_closure(lw_two_stream_closure, lw_gamma, lw_gammaprime, lw_ghat_factor, lw_preset)
if ((sw_preset .or. lw_preset) .and. half_tau_diffuse) then
  half_tau_diffuse = .false.
  call error_mesg('two_stream_scatter_rad', &
       'closure preset given: half_tau_diffuse set to .false. (Eq. 5.30 normalisation).', NOTE)
endif
if (half_tau_diffuse) then
  tau_fac = 0.5
else
  tau_fac = 1.0
endif
if (ghat_factor <= 0.0) call error_mesg('two_stream_scatter_rad','ghat_factor must be > 0.', FATAL)
if (mpp_pe() == 0) then
  write(*,'(a,3f9.5)') ' two_stream_scatter_rad: SW closure '//trim(sw_two_stream_closure)// &
       '  gamma, gammaprime, ghat_factor =', sw_gamma, sw_gammaprime, sw_ghat_factor
  write(*,'(a,3f9.5)') ' two_stream_scatter_rad: LW closure '//trim(lw_two_stream_closure)// &
       '  gamma, gammaprime, ghat_factor =', lw_gamma, lw_gammaprime, lw_ghat_factor
  write(*,'(a,l2,a,l2)') ' two_stream_scatter_rad: half_tau_diffuse =', half_tau_diffuse, &
       '   do_delta_scaling =', do_delta_scaling
endif
if (g_asym <= -1.0 .or. g_asym >= 1.0) call error_mesg('two_stream_scatter_rad','g_asym must lie in (-1, 1).', FATAL)

if(dt_rad_avg .le. 0) dt_rad_avg = dt_real !s if dt_rad_avg is set to a value in nml then it will be used instead of dt_real

if (dt_rad > 0.0) then
  rad_step_interval = nint(dt_rad / dt_real)
  if (rad_step_interval < 1 .or. abs(rad_step_interval*dt_real - dt_rad) > 1.e-6*dt_rad) then
    call error_mesg('two_stream_scatter_rad','dt_rad must be a positive multiple of the model timestep.', FATAL)
  endif
else
  rad_step_interval = 1
endif
rad_step_counter = 0
id_clock_rad = mpp_clock_id('Two-stream scatter radiation')

if(do_read_co2)then
   call interpolator_init (co2_interp, trim(co2_file)//'.nc', lonb, latb, data_out_of_bounds=(/ZERO/))
endif

if (uppercase(trim(sw_optical_depth)) == 'TRANSPARENT') then 
  sw_flag = sw_TRANSPARENT
elseif (uppercase(trim(sw_optical_depth)) == 'GENERIC') then
  sw_flag = sw_GENERIC
else
  call error_mesg('two_stream_scatter_rad','"'//trim(sw_optical_depth)//'"'//' is not a valid sw parametrisation.', FATAL)
endif

initialized = .true.

allocate (sw_tau_0         (ie-is+1, je-js+1))
allocate (insolation       (ie-is+1, je-js+1))
allocate (p2               (ie-is+1, je-js+1))

allocate (pi_B_surf        (ie-is+1, je-js+1))

allocate (pi_B             (ie-is+1, je-js+1, num_levels))
allocate (tdt_rad          (ie-is+1, je-js+1, num_levels))
allocate (tdt_solar        (ie-is+1, je-js+1, num_levels))

allocate (lw_up            (ie-is+1, je-js+1, num_levels+1))
allocate (lw_down          (ie-is+1, je-js+1, num_levels+1))
allocate (lw_flux          (ie-is+1, je-js+1, num_levels+1))
allocate (sw_up            (ie-is+1, je-js+1, num_levels+1))
allocate (sw_down          (ie-is+1, je-js+1, num_levels+1))
allocate (sw_down_direct   (ie-is+1, je-js+1, num_levels+1))
allocate (sw_flux          (ie-is+1, je-js+1, num_levels+1))
allocate (rad_flux         (ie-is+1, je-js+1, num_levels+1))

allocate (lw_tau           (ie-is+1, je-js+1, num_levels+1))
allocate (sw_tau           (ie-is+1, je-js+1, num_levels+1))
allocate (lw_del_tau       (ie-is+1, je-js+1, num_levels))
allocate (sw_del_tau       (ie-is+1, je-js+1, num_levels))

allocate (lw_ss_albedo     (ie-is+1, je-js+1, num_levels))
allocate (sw_ss_albedo     (ie-is+1, je-js+1, num_levels))

allocate (lw_abs_coeff     (ie-is+1, je-js+1, num_levels))
allocate (sw_abs_coeff     (ie-is+1, je-js+1, num_levels))
allocate (lw_scatter_coeff (ie-is+1, je-js+1, num_levels))
allocate (sw_scatter_coeff (ie-is+1, je-js+1, num_levels))

allocate (sw_gammab        (ie-is+1, je-js+1, num_levels))
allocate (sw_gammaminus    (ie-is+1, je-js+1, num_levels))
allocate (sw_gammaplus     (ie-is+1, je-js+1, num_levels))
allocate (sw_gammaone      (ie-is+1, je-js+1, num_levels))
allocate (sw_gammatwo      (ie-is+1, je-js+1, num_levels))

allocate (lw_gammab        (ie-is+1, je-js+1, num_levels))
allocate (lw_gammaminus    (ie-is+1, je-js+1, num_levels))
allocate (lw_gammaplus     (ie-is+1, je-js+1, num_levels))
allocate (lw_gammaone      (ie-is+1, je-js+1, num_levels))
allocate (lw_gammatwo      (ie-is+1, je-js+1, num_levels))

allocate (olr              (ie-is+1, je-js+1))
allocate (net_lw_surf      (ie-is+1, je-js+1))
allocate (toa_sw_in        (ie-is+1, je-js+1))
allocate (toa_sw_out       (ie-is+1, je-js+1))
allocate (coszen           (ie-is+1, je-js+1))
allocate (fracsun          (ie-is+1, je-js+1)) !jp from astronomy.f90 : fraction of sun on surface


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
    id_swup_toa = &
    register_diag_field ( mod_name, 'swup_toa', axes(1:2), Time, &
                'SW flux up at TOA', &
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

    id_flux_lw_up = &
        register_diag_field ( mod_name, 'flux_lw_up', axes(half), Time, &
               'Upward longwave radiative flux (positive up)', &
               'W/m^2', missing_value=missing_value               )
    id_flux_lw_down = &
        register_diag_field ( mod_name, 'flux_lw_down', axes(half), Time, &
               'Downward longwave radiative flux (positive down)', &
               'W/m^2', missing_value=missing_value               )
    id_flux_sw_up = &
        register_diag_field ( mod_name, 'flux_sw_up', axes(half), Time, &
               'Upward shortwave radiative flux (positive up)', &
               'W/m^2', missing_value=missing_value               )
    id_flux_sw_down = &
        register_diag_field ( mod_name, 'flux_sw_down', axes(half), Time, &
               'Downward (diffuse) shortwave radiative flux (positive down)', &
               'W/m^2', missing_value=missing_value               )
    id_flux_sw_down_direct = &
        register_diag_field ( mod_name, 'flux_sw_down_direct', axes(half), Time, &
               'Downward direct shortwave radiative flux (positive down)', &
               'W/m^2', missing_value=missing_value               )

    id_coszen  = &
               register_diag_field ( mod_name, 'coszen', axes(1:2), Time, &
                 'cosine of zenith angle', &
                 'none', missing_value=missing_value      )
    id_fracsun  = &
               register_diag_field ( mod_name, 'fracsun', axes(1:2), Time, &
                 'daylight fraction of time interval', &
                 'none', missing_value=missing_value      )

    id_co2  = &
               register_diag_field ( mod_name, 'co2', Time, &
                 'co2 concentration', &
                 'ppmv', missing_value=missing_value      )

return
end subroutine two_stream_scatter_rad_init

! ==================================================================================

subroutine two_stream_scatter_rad (is, js, ie, je, Time_diag, lat, lon, &
                                   p_half, p_full, &
                                   t_surf, t, tdt, net_surf_sw_down, surf_lw_down, albedo, q)


integer, intent(in)                 :: is, js, ie, je
type(time_type), intent(in)         :: Time_diag
real, intent(in), dimension(:,:)    :: lat, lon, albedo
real, intent(out), dimension(:,:)   :: net_surf_sw_down
real, intent(out), dimension(:,:)   :: surf_lw_down
real, intent(in), dimension(:,:,:)  :: t, q,  p_half, p_full
real, intent(in), dimension(:,:)    :: t_surf
real, intent(inout), dimension(:,:,:) :: tdt
integer :: i, j, k, n, dyofyr

integer :: seconds, year_in_s, days
real :: r_seconds, frac_of_day, frac_of_year, gmt, time_since_ae, rrsun, day_in_s, r_solday, r_total_seconds, r_days, r_dt_rad_avg, dt_rad_radians
logical :: used

real ,dimension(size(q,1),size(q,2),size(q,3)) :: co2f
real ,dimension(size(q,1),size(q,2),size(q,3)) :: ext_tmp   ! extinction coefficient (abs + sca)
real :: g_eff, ghat                                          ! delta-scaled g, ghat = ghat_factor*g_eff

! close-in insolation
real :: ci_lat_f, ci_lat_d
real, dimension(size(q,1), size(q,2)) :: lat_tl_insol, ci_rho_dist, ci_J_f, ci_alpha, ci_beta, &
                                         ci_phi, ci_t, ci_l, ci_s, ci_delta1, ci_delta2, ci_J_p

n = size(t,3)

! Radiation timestep: on non-radiation steps reuse the stored heating rate and
! surface fluxes (module arrays from the last full calculation).
rad_step_counter = rad_step_counter + 1
if (rad_step_counter > 1 .and. mod(rad_step_counter-1, rad_step_interval) /= 0) then
  tdt = tdt + tdt_rad
  surf_lw_down     = lw_down(:, :, n+1)
  net_surf_sw_down = (sw_down_direct(:,:,n+1)+sw_down(:, :, n+1)) * (1. - albedo)
  return
endif

call mpp_clock_begin(id_clock_rad)

if(do_read_co2)then
  call interpolator( co2_interp, Time_diag, p_half, co2f, trim(co2_variable_name))
  carbon_conc = maxval(co2f) !Needs maxval just because co2f is a 3d array of a constant, so maxval is just a way to pick out one number
endif

! =================================================================================
! COMPUTE INSOLATION 

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

else if (do_tl) then
  ! Tidally locked: substellar point fixed at lat = 0, lon = 0

  if (do_closein) then
    ! finite stellar disc (close-in planet)
    lat_tl_insol = asin(cos(lat)*cos(lon)) + pi/2

    ci_lat_f = acos((radius + R_stellar) / d_stellar)
    ci_lat_d = acos((radius - R_stellar) / d_stellar)

    ci_rho_dist = sqrt(d_stellar**2 + radius**2 - 2*d_stellar*radius*cos(lat_tl_insol))

    ! J_f: fully illuminated region
    ci_J_f = (d_stellar*cos(lat_tl_insol) - radius) * R_stellar**2 / ci_rho_dist**3

    ! J_p: partially illuminated region
    ci_alpha  = acos((d_stellar*cos(lat_tl_insol) - radius) / ci_rho_dist)
    ci_beta   = pi/2 - ci_alpha
    ci_phi    = asin(R_stellar / ci_rho_dist)
    ci_t      = ci_rho_dist * cos(ci_phi)
    ci_l      = ci_t * sin(ci_phi)
    ci_s      = ci_t * sin(ci_phi) * cos(ci_alpha)
    ci_delta1 = acos(cos(ci_phi)/cos(ci_beta))
    ci_delta2 = acos(tan(ci_beta)/tan(ci_phi))
    ci_J_p    = ci_l*ci_s/pi/ci_t**2 * (pi - ci_delta2 + sin(ci_delta2)*cos(ci_delta2)) &
              + 1/pi * (ci_delta1 - sin(ci_delta1)*cos(ci_delta1))

    where (lat_tl_insol .le. ci_lat_f)
      insolation = solar_constant * (d_stellar/R_stellar)**2 * ci_J_f
    elsewhere ((lat_tl_insol .gt. ci_lat_f) .and. (lat_tl_insol .lt. ci_lat_d))
      insolation = solar_constant * (d_stellar/R_stellar)**2 * ci_J_p
    elsewhere
      insolation = 0.
    endwhere

    ! effective mu for the direct-beam slant path; clipped, since the close-in
    ! flux can exceed S*cos(theta) and remains non-zero beyond the geometric terminator
    coszen = min(max(insolation / solar_constant, 0.0), 1.0)

  else
    ! point-source star
    insolation = max(solar_constant * cos(lat) * cos(lon), 0.0)
    coszen     = insolation / solar_constant       ! = cos(lat)cos(lon) on the dayside: exact
  endif
else
  ! Default: Averaged Earth insolation at all longitudes
  p2          = (1. - 3.*sin(lat)**2)/4.
  insolation  = 0.25 * solar_constant * (1.0 + del_sol * p2 + del_sw * sin(lat))
  if (coszen_fixed > 0.0) then
   coszen = coszen_fixed                      ! slant path set independently of flux
  else
   coszen = insolation / solar_constant
  endif
end if


! =================================================================================
! COMPUTE SHORTWAVE RADIATION 

! zero fluxes before starting 
sw_down_direct = 0.0
sw_up          = 0.0
sw_down        = 0.0
sw_scatter_coeff = 0.0
sw_abs_coeff     = 0.0
sw_ss_albedo     = 0.0 


select case (sw_flag)
case(sw_TRANSPARENT)
  ! just reflect at surface (direct down is reflected into diffuse up)
  do k = 1, n+1
    sw_down_direct(:,:,k) = insolation(:,:)
  enddo 
  do k = 1, n+1
    sw_up(:,:,k) = albedo(:,:) * sw_down_direct(:,:,n+1)
  enddo 

case(sw_GENERIC)
  ! Absorption and scattering in shortwave applied to both the upward and 
  ! downward beams. Total SW flux is made up of three components: a downward 
  ! direct component, and upward and downward diffuse components. Direct 
  ! radiation impinging on the surface is reflected into the upward diffuse 
  ! beam. 
  ! 
  ! First order differencing is applied to derivatives w.r.t tau. 
  !
  ! Once optical depths have been defined, the problem is solved as a matrix 
  ! equation of the form Ax = B, where x represents the upward and downward 
  ! fluxes (combined into a single vector). A is a banded matrix, and the problem 
  ! is solved by LU decomposition. This is carried out by the routine 
  ! two_stream_solver. 
  !
  ! Method follows Pierrehumbert (2010): Princples of Planetary Climate. Chapter 
  ! 5, Page 353. 

  
  ! Compute optical properties of atmosphere 
  ! Compute scattering coefficient (currently no contribution from carbon_conc, easily added though...)
  sw_scatter_coeff  = sw_sca_a + sw_sca_b * q 
  ! Compute absorption coefficient 
  sw_abs_coeff = sw_abs_a + sw_abs_b*q 
  ! Compute single scattering albedo 
  ext_tmp = sw_abs_coeff + sw_scatter_coeff
  where(ext_tmp .ne. 0.0) 
    sw_ss_albedo = sw_scatter_coeff / ext_tmp
  endwhere 
  ! delta-Eddington scaling (forward peak treated as unscattered)
  call delta_scale(ext_tmp, sw_ss_albedo, g_eff)
  ghat = sw_ghat_factor * g_eff
  if (.not. warned_ghat_sw .and. ghat >= 1.0 .and. any(sw_ss_albedo > 0.0)) then
    warned_ghat_sw = .true.
    call error_mesg('two_stream_scatter_rad','SW ghat >= 1: two-stream breaks down (1 - ghat*omega <= 0); '// &
         'use do_delta_scaling or another sw_two_stream_closure.', NOTE)
  endif
  ! Compute SW del_tau (negative: tau = 0 at surface, increasing upward)
  sw_del_tau =  1 / grav * ext_tmp * & 
               (p_half(:,:,1:n) - p_half(:,:,2:n+1) )


  ! Compute SW optical depth 
  ! define optical depth to be zero at surface, following Pierrehumbert(2010)
  sw_tau(:,:,n+1) = 0.0 
  do k = 1, n
    sw_tau(:,:,n+1-k) = sw_tau(:,:,n+2-k) - sw_del_tau(:,:,n+1-k)
  enddo 

  !! gammas (as in Pierrehumbert, 2010)  
  ! Pierrehumbert Eq. 5.30 with ghat = sw_ghat_factor * g (custom, 3/2 reproduces v1.4)
  sw_gammaone = sw_gamma *(1.0 - ghat*sw_ss_albedo) + sw_gammaprime*(1-sw_ss_albedo)
  sw_gammatwo = sw_gamma *(1.0 - ghat*sw_ss_albedo) - sw_gammaprime*(1-sw_ss_albedo)
  sw_gammab = sw_gammaone - sw_gammatwo 
  do k = 1, n
    sw_gammaplus(:,:,k)  = 0.5*sw_ss_albedo(:,:,k) - sw_gamma*sw_ss_albedo(:,:,k)*ghat*coszen(:,:) 
    sw_gammaminus(:,:,k) = 0.5*sw_ss_albedo(:,:,k) + sw_gamma*sw_ss_albedo(:,:,k)*ghat*coszen(:,:) 
  enddo 
  if (.not. warned_gammaplus .and. any(sw_gammaplus < -1.e-12)) then
    warned_gammaplus = .true.
    call error_mesg('two_stream_scatter_rad', &
         'gamma+ < 0 somewhere (negative upward scattering of the direct beam): '// &
         'reduce g_asym or set do_delta_scaling = .true.', NOTE)
  endif

  ! Compute direct downward radiation: exact attenuation exp(-|dtau|/mu).
  ! sw_del_tau < 0, so exp(sw_del_tau/coszen) = exp(-|dtau|/mu). The exponent
  ! is always <= 0, so small coszen drives the beam to zero (no overflow, never negative).
  sw_down_direct(:,:,1) = insolation(:,:)
  do k = 1, n
    where (coszen > 0.0)
      sw_down_direct(:,:,k+1) = sw_down_direct(:,:,k) * exp(sw_del_tau(:,:,k) / coszen(:,:))
    elsewhere
      sw_down_direct(:,:,k+1) = 0.0
    endwhere
    ! flush a negligible beam to zero (avoids slow denormal arithmetic downstream)
    where (sw_down_direct(:,:,k+1) < 1.e-20*insolation) sw_down_direct(:,:,k+1) = 0.0
  enddo
  
  !!!! Zero blackbody radiation for shortwave !!!!
  pi_B_surf = 0.0
  pi_B = 0.0


  do i = 1, (ie-is+1)
    do j = 1, (je-js+1)
      if (coszen(i,j) .le. 0.0) then 
        sw_up(i,j,:) = 0.0
        sw_down(i,j,:) = 0.0
      else 
      call two_stream_solver(n, albedo(i,j), coszen(i,j), pi_B_surf(i,j), &
                           pi_B(i,j,:), sw_down_direct(i,j,:), &
                           sw_del_tau(i,j,:)*tau_fac, sw_gammaone(i,j,:), &
                           sw_gammatwo(i,j,:), sw_gammab(i,j,:), &
                           sw_gammaplus(i,j,:), sw_gammaminus(i,j,:), &
                           sw_up(i,j,:), sw_down(i,j,:)) ! zero contribution from in-atmosphere emission (see above)
                           ! tau_fac = 1 (Eq. 5.30) or 0.5 (v1.4 legacy, half_tau_diffuse); diffuse streams only
      endif 
    
    enddo
  enddo 

end select 


! =================================================================================
! LONGWAVE RADIATION
! Supports absorption and scattering in the longwave. 
! 
! First order differencing is applied to derivatives w.r.t tau. 
!
! Once optical depths have been defined, the problem is solved as a matrix 
! equation of the form Ax = B, where x represents the upward and downward 
! fluxes (combined into a single vector). A is a banded matrix, and the problem 
! is solved by LU decomposition. This is carried out by the routine 
! two_stream_solver. 
!
! Method follows Pierrehumbert (2010): Princples of Planetary Climate. Chapter 
! 5, Page 353. 

! zero fluxes before starting 
lw_up            = 0.0
lw_down          = 0.0 
lw_scatter_coeff = 0.0
lw_abs_coeff     = 0.0
lw_ss_albedo     = 0.0

! Compute optical properties of atmosphere 
! Compute scattering coefficient (currently no contribution from carbon_conc, easily added though...)
lw_scatter_coeff  = lw_sca_a + lw_sca_b * q 
! Compute absorption coefficient 
lw_abs_coeff = lw_abs_a + lw_abs_b * q + lw_abs_c * log(carbon_conc / 360.)
! Compute single scattering albedo 
ext_tmp = lw_abs_coeff + lw_scatter_coeff
where (ext_tmp .ne. 0.0) 
  lw_ss_albedo = lw_scatter_coeff / ext_tmp
endwhere 
! delta-Eddington scaling (no effect when lw_ss_albedo = 0)
call delta_scale(ext_tmp, lw_ss_albedo, g_eff)
ghat = lw_ghat_factor * g_eff
if (.not. warned_ghat_lw .and. ghat >= 1.0 .and. any(lw_ss_albedo > 0.0)) then
  warned_ghat_lw = .true.
  call error_mesg('two_stream_scatter_rad','LW ghat >= 1: two-stream breaks down (1 - ghat*omega <= 0); '// &
       'use do_delta_scaling or another lw_two_stream_closure.', NOTE)
endif

! Compute lw del_tau 
lw_del_tau = 1 / grav * ext_tmp * & 
             (p_half(:,:,1:n) - p_half(:,:,2:n+1))

! Compute lw optical depth 
!define optical depth to be zero at surface, following Pierrehumbert (2010)
lw_tau(:,:,n+1) = 0.0 
do k = 1, n
  lw_tau(:,:,n+1-k) = lw_tau(:,:,n+2-k) - lw_del_tau(:,:,n+1-k)
enddo 

!! gammas, as defined in Pierrehumbert (2010)
lw_gammaone = lw_gamma *(1.0 - ghat*lw_ss_albedo) + lw_gammaprime*(1-lw_ss_albedo)
lw_gammatwo = lw_gamma *(1.0 - ghat*lw_ss_albedo) - lw_gammaprime*(1-lw_ss_albedo)
lw_gammab = lw_gammaone - lw_gammatwo 
! no direct beam in the LW, so gamma+- only multiply zero; kept for the solver interface
do k = 1, n
  lw_gammaplus(:,:,k)  = 0.5*lw_ss_albedo(:,:,k)
  lw_gammaminus(:,:,k) = 0.5*lw_ss_albedo(:,:,k)
enddo 

! Compute black body radiation for longwave 
pi_B_surf = stefan * t_surf * t_surf * t_surf * t_surf
pi_B = stefan * t * t * t * t

do i = 1, (ie-is+1)
  do j = 1, (je-js+1)
  call two_stream_solver(n, 0.0, 1.0, pi_B_surf(i,j), &
                         pi_B(i,j,:), sw_down_direct(i,j,:)*0.0, &
                         lw_del_tau(i,j,:)*tau_fac, lw_gammaone(i,j,:), &
                         lw_gammatwo(i,j,:), lw_gammab(i,j,:), &
                         lw_gammaplus(i,j,:), lw_gammaminus(i,j,:), &
                         lw_up(i,j,:), lw_down(i,j,:))
  ! zero contribution from direct solar beam 
  ! tau_fac = 1 (Eq. 5.30) or 0.5 (v1.4 legacy, half_tau_diffuse)
  enddo
enddo 


! =================================================================================
! Compute net fluxes (positive up) and heating rates 

lw_flux  = lw_up - lw_down
sw_flux  = sw_up - sw_down - sw_down_direct
rad_flux = lw_flux + sw_flux

do k = 1, n
   tdt_rad(:,:,k)   = diabatic_acce * ( rad_flux(:,:,k+1) - rad_flux(:,:,k) )  &
        * grav/( cp_air*(p_half(:,:,k+1) - p_half(:,:,k)) )

   tdt_solar(:,:,k) = ( sw_flux(:,:,k+1) - sw_flux(:,:,k) )  &
        * grav/( cp_air*(p_half(:,:,k+1) - p_half(:,:,k)) )

   tdt(:,:,k) = tdt(:,:,k) + tdt_rad(:,:,k)
end do

! =================================================================================
surf_lw_down     = lw_down(:, :, n+1)
toa_sw_in        = sw_down_direct(:, :, 1)
toa_sw_out       = sw_up(:,:,1)
net_surf_sw_down = (sw_down_direct(:,:,n+1)+sw_down(:, :, n+1)) * (1. - albedo)
olr         = lw_up(:,:,1)
net_lw_surf = lw_flux(:, :, n+1)
! =================================================================================


!------- downward lw flux surface -------
if ( id_lwdn_sfc > 0 ) then
  used = send_data ( id_lwdn_sfc, surf_lw_down, Time_diag)
endif
!------- incoming sw flux toa -------
if ( id_swdn_toa > 0 ) then
  used = send_data ( id_swdn_toa, toa_sw_in, Time_diag)
endif
!------- outgoing sw flux toa -------
if ( id_swup_toa > 0 ) then
  used = send_data ( id_swup_toa, toa_sw_out, Time_diag)
endif
!------- downward sw flux surface -------
if ( id_swdn_sfc > 0 ) then
  used = send_data ( id_swdn_sfc, net_surf_sw_down, Time_diag)
endif

!------- cosine of zenith angle ------------
if ( id_coszen > 0 ) then
  used = send_data ( id_coszen, coszen, Time_diag)
endif

!------- carbon dioxide concentration ------------
if ( id_co2 > 0 ) then
  used = send_data ( id_co2, carbon_conc, Time_diag)
endif

!------- outgoing lw flux toa (olr) -------
if ( id_olr > 0 ) then
   used = send_data ( id_olr, olr, Time_diag)
endif
!------- upward lw flux surface -------
if ( id_lwup_sfc > 0 ) then
   used = send_data ( id_lwup_sfc, pi_B_surf, Time_diag)
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
! Upward and downward fluxes 
if ( id_flux_lw_up > 0 ) then
  used = send_data ( id_flux_lw_up, lw_up, Time_diag)
endif
if ( id_flux_lw_down > 0 ) then
  used = send_data ( id_flux_lw_down, lw_down, Time_diag)
endif
if ( id_flux_sw_up > 0 ) then
  used = send_data ( id_flux_sw_up, sw_up, Time_diag)
endif
if ( id_flux_sw_down > 0 ) then
  used = send_data ( id_flux_sw_down, sw_down, Time_diag)
endif
if ( id_flux_sw_down_direct > 0 ) then
  used = send_data ( id_flux_sw_down_direct, sw_down_direct, Time_diag)
endif

call mpp_clock_end(id_clock_rad)

return
end subroutine two_stream_scatter_rad

! ==================================================================================


subroutine two_stream_scatter_rad_end

  deallocate (sw_tau_0         )
  deallocate (insolation       )
  deallocate (p2               )
  
  deallocate (pi_B_surf        )
  
  deallocate (pi_B             )
  deallocate (tdt_rad          )
  deallocate (tdt_solar        )
  
  deallocate (lw_up            )
  deallocate (lw_down          )
  deallocate (lw_flux          )
  deallocate (sw_up            )
  deallocate (sw_down          )
  deallocate (sw_down_direct   )
  deallocate (sw_flux          )
  deallocate (rad_flux         )
  
  deallocate (lw_tau           )
  deallocate (sw_tau           )
  deallocate (lw_del_tau       )
  deallocate (sw_del_tau       )
  
  deallocate (lw_ss_albedo     )
  deallocate (sw_ss_albedo     )
  
  deallocate (lw_abs_coeff     )
  deallocate (sw_abs_coeff     )
  deallocate (lw_scatter_coeff )
  deallocate (sw_scatter_coeff )
  
  deallocate (sw_gammab        )
  deallocate (sw_gammaminus    )
  deallocate (sw_gammaplus     )
  deallocate (sw_gammaone      )
  deallocate (sw_gammatwo      )
  
  deallocate (lw_gammab        )
  deallocate (lw_gammaminus    )
  deallocate (lw_gammaplus     )
  deallocate (lw_gammaone      )
  deallocate (lw_gammatwo      )
  
  deallocate (olr              )
  deallocate (net_lw_surf      )
  deallocate (toa_sw_in        )
  deallocate (toa_sw_out       )
  deallocate (coszen           )
  deallocate (fracsun          ) 

if(do_read_co2)call interpolator_end(co2_interp)

end subroutine two_stream_scatter_rad_end

! ==================================================================================

subroutine set_closure(name, g, gp, gh, preset)
! Resolve a closure name into (gamma, gammaprime, ghat_factor); Pierrehumbert
! Table 5.3 and Eq. 5.21. 'custom' uses the namelist gamma, gammaprime, ghat_factor.
character(len=*), intent(in) :: name
real, intent(out)    :: g, gp, gh
logical, intent(out) :: preset
preset = .true.
select case (uppercase(trim(name)))
case ('CUSTOM')
  g = gamma;             gp = gammaprime;            gh = ghat_factor;  preset = .false.
case ('HEMI', 'HEMI-ISOTROPIC', 'HEMI_ISOTROPIC')
  g = 1.0;               gp = 1.0;                   gh = 1.5
case ('QUADRATURE')
  g = 0.5*sqrt(3.0);     gp = 0.5*sqrt(3.0);         gh = 1.0
case ('EDDINGTON')
  g = 0.75;              gp = 1.0;                   gh = 1.0
case default
  call error_mesg('two_stream_scatter_rad','"'//trim(name)// &
       '" is not a valid closure (custom, hemi, quadrature, eddington).', FATAL)
end select
end subroutine set_closure

! ==================================================================================

subroutine delta_scale(ext, ssalb, g_out)
! delta-Eddington scaling (Joseph, Wiscombe & Weinman 1976), applied in place when
! do_delta_scaling: f = g^2, ext' = (1 - w f) ext, w' = (1 - f) w / (1 - w f),
! g' = (g - f)/(1 - f) = g/(1 + g). Without scaling, returns g_out = g_asym unchanged.
real, intent(inout), dimension(:,:,:) :: ext, ssalb
real, intent(out) :: g_out
real :: f
if (do_delta_scaling .and. g_asym /= 0.0) then
  f = g_asym * g_asym
  where (ssalb > 0.0)
    ext   = ext * (1.0 - ssalb*f)
    ssalb = (1.0 - f) * ssalb / (1.0 - ssalb*f)
  endwhere
  g_out = g_asym / (1.0 + g_asym)
else
  g_out = g_asym
endif
end subroutine delta_scale

! ==================================================================================

subroutine two_stream_solver(nlevels, alpha_g, coszen, pi_B_surf, pi_B,  I_direct, del_tau, gammaone, & 
  gammatwo, gammab, gammaplus, gammaminus, I_up, I_down)

! del_tau is the layer optical depth times tau_fac (1 for Eq. 5.30, 0.5 in v1.4
! legacy mode; negative in this code's convention). I_direct is the direct-beam flux at half levels (zero for LW).
! coszen is retained in the interface (gammaplus/gammaminus depend on it via g*mu)
! but is no longer used in the source term.

integer, intent(in) :: nlevels 
real, intent(in)              :: alpha_g, coszen,  pi_B_surf  ! pi_B_surf zero if shortwave
real, intent(in), dimension(nlevels) :: pi_B,  del_tau, gammaone, gammatwo, &
         gammab, gammaplus, gammaminus 
         !pi_B zero if SW, I_direct zero if longwave  
real, intent(in), dimension(nlevels+1) :: I_direct


real, intent(out), dimension(nlevels+1) :: I_up, I_down  

integer :: N  ! N = 2*(nlevels + 1)
integer, parameter :: MU = 3, ML = 3 ! number of upper and lower diagonals 
! band storage only: LINPACK layout needs 2*ML+MU+1 rows (rows 1..ML for pivot fill-in)
integer, parameter :: LDB = 2*ML + MU + 1, MB = ML + MU + 1   ! MB = row of the main diagonal
real, dimension(LDB, 2*(nlevels+1)) :: B
real, dimension(2*(nlevels+1)) :: I_comb, source_term
integer :: i, lcounter1, lcounter2, IND
integer, dimension(2*(nlevels+1)) :: IPVT

N = 2*(nlevels + 1)


! MAKE MATRIX OF ABSORPTION / SCATTERING COEFFICIENTS - band diagonal matrix 
! entries stored directly in band form: full-matrix A(i,j) lives at B(i-j+MB, j)
B = 0.0
! fill main diagonal 
lcounter1 = 1
lcounter2 = 1
do i = 1, N
if (i.eq.2) then 
B((i)-(i)+MB, i) = 1 ! upper bc
elseif (i.eq.N-1) then 
B((i)-(i)+MB, i) = 1
else
if (mod(i,2).eq.1) then 
B((i)-(i)+MB, i) = -1 + gammaone(lcounter1) / 2 * del_tau(lcounter1)
lcounter1 = lcounter1 + 1
else 
B((i)-(i)+MB, i) = 1 - gammaone(lcounter2) / 2 * del_tau(lcounter2)
lcounter2 = lcounter2 + 1
endif 
endif 
enddo 

! fill first upper and lower diagonal 
lcounter1 = 1
lcounter2 = 1
do i = 1, N-1
if (mod(i,2).eq.1) then
if (i.eq.N-1) then 
B((i)-(i+1)+MB, i+1) = - alpha_g ! lower bc
else 
B((i)-(i+1)+MB, i+1) = -gammatwo(lcounter1) / 2 * del_tau(lcounter1)
lcounter1 = lcounter1 + 1
endif 
if (i.ne.1) then ! for i = 1, A(2,1) = 0 for upper bc
B((i+1)-(i)+MB, i) = gammatwo(lcounter2) / 2 * del_tau(lcounter2)
lcounter2 = lcounter2 + 1
endif 
endif 
enddo 

! fill second upper and lower diagonal 
lcounter1 = 1
lcounter2 = 1
do i = 1, N-2
if (mod(i,2).eq.1) then 
B((i)-(i+2)+MB, i+2) = 1+ gammaone(lcounter1) / 2 * del_tau(lcounter1)
lcounter1 = lcounter1 + 1
else 
B((i+2)-(i)+MB, i) = -1 - gammaone(lcounter2) / 2 * del_tau(lcounter2)
lcounter2 = lcounter2 + 1
endif 
enddo 

!fill third upper and lower diagonal 
lcounter1 = 1
lcounter2 = 1
do i = 1, N-3
if (mod(i,2).eq.1) then 
B((i)-(i+3)+MB, i+3) = - gammatwo(lcounter1) / 2 * del_tau(lcounter1)
lcounter1 = lcounter1 + 1
B((i+3)-(i)+MB, i) = gammatwo(lcounter2) / 2 * del_tau(lcounter2)
lcounter2 = lcounter2 + 1
endif 
enddo 


! Make source term
! Direct-beam scattering source uses the actual beam loss in each layer,
! dF(k) = I_direct(k) - I_direct(k+1) >= 0, split up/down by gammaplus/gammaminus
! (gammaplus + gammaminus = omega), so scattered energy is conserved exactly.
! Upward-equation rows (odd) take -gammaplus*dF; downward rows (even) take +gammaminus*dF.
source_term = 0.0
source_term(1) = del_tau(1) * gammab(1) * pi_B(1) &
               - gammaplus(1) * (I_direct(1) - I_direct(2))
lcounter1 = 2
lcounter2 = 1
do i = 3, N-2
  if (mod(i,2).eq.1) then
    source_term(i) = del_tau(lcounter1) * gammab(lcounter1) * pi_B(lcounter1) &
                   - gammaplus(lcounter1) * (I_direct(lcounter1) - I_direct(lcounter1+1))
    lcounter1 = lcounter1 + 1
  else
    source_term(i) = -del_tau(lcounter2) * gammab(lcounter2) * pi_B(lcounter2) &
                   + gammaminus(lcounter2) * (I_direct(lcounter2) - I_direct(lcounter2+1))
    lcounter2 = lcounter2 + 1
  endif
enddo
source_term(N-1) = alpha_g * I_direct(nlevels+1) + (1 - alpha_g) * pi_B_surf
source_term(N)   = -del_tau(nlevels) * gammab(nlevels) * pi_B(nlevels) &
                 + gammaminus(nlevels) * (I_direct(nlevels) - I_direct(nlevels+1))


! call LU factorization  
CALL NSBFAC(B,LDB,N,ML,MU,IPVT,IND)
! solve AX = B for X 
if (IND .ne. 0) then 
write(*,*) 'MATRIX IS ILL-CONDITIONED'
endif 
CALL NSBSLV(B,LDB,N,ML,MU,IPVT,source_term,I_comb)

lcounter1 = 1
lcounter2 = 1
do i = 1, N
if (mod(i,2).eq.1) then 
I_up(lcounter1) = I_comb(i)
lcounter1 = lcounter1 + 1
else 
I_down(lcounter2) = I_comb(i)
lcounter2 = lcounter2 + 1
endif 
enddo 


end subroutine two_stream_solver

SUBROUTINE NSBFAC(B,LDB,N,ML,MU,IPVT,IND)
!-------------------------------------------------------------------
!     LU factorization of a band matrix (non symmetric) with partial
!     pivoting.
!     INPUTS:
!     B  : banded matrix. The correspondance between full matrix
!     A(i,j) and band matrix B(k,l) is given by following sequence:
!     m=ml+mu+1
!     do j=1,n
!       i1=max(1,j-mu)
!       i2=min(n,j+ml)
!       do i=i1,i2
!         k=i-j+m
!         b(k,j)=a(i,j)
!       enddo
!     enddo
!     LDB : 1st dimension of B in calling program (ldb.ge.2ml+mu+1)
!     N   : size of B                             -----------------
!     ML  : number of lower diagonals
!     MU  : number of upper diagonals
!     OUTPUTS:
!     B   : banded matrix storing the LU elements, the ML first lines
!           of which are used for pivoting.
!     IPVT: integer vector of size N storing the pivoting indices.
!     IND : flag = 0,  B is non singular
!                = k,  B may be singular
!-------------------------------------------------------------------
INTEGER :: N, LDB, ML, MU, IND, OUT
INTEGER :: I, M, J0, J1, JZ, I0, J, JU, K, KP1, L, LM, MM, n_m_one 
REAL :: B(LDB,N),T
INTEGER IPVT(*)
M=ML+MU+1
IND=0

J0=MU+2
J1=MIN(N,M)-1
IF(J1.GE.J0) THEN
DO JZ=J0,J1
I0=M+1-JZ
DO I=I0,ML
B(I,JZ)=0.D0
ENDDO
ENDDO
ENDIF
JZ=J1
JU=0

n_m_one=N-1
IF(n_m_one.GE.1) THEN
DO K=1,n_m_one
KP1=K+1
JZ=JZ+1
IF(JZ.LE.N.AND.ML.GE.1) THEN
DO I=1,ML
B(I,JZ)=0.D0
ENDDO
ENDIF

LM=MIN(ML,N-K)
call IAMAX(B(M,K),LM+1, OUT)
L=OUT+M-1
IPVT(K)=L+K-M
IF(B(L,K).EQ.0.D0) GO TO 10
IF(L.NE.M) THEN
T=B(L,K)
B(L,K)=B(M,K)
B(M,K)=T
ENDIF
T=-1.D0/B(M,K)
CALL SCALE(LM,T,B(M+1,K))

JU=MIN(MAX(JU,MU+IPVT(K)),N)
MM=M
IF(JU.GE.KP1) THEN
DO J=KP1,JU
L=L-1
MM=MM-1
T=B(L,J)
IF(L.NE.MM) THEN
B(L,J)=B(MM,J)
B(MM,J)=T
ENDIF
CALL DAXPY(LM,T,B(M+1,K),B(MM+1,J))
ENDDO
ENDIF
GO TO 20
10     ind=k
20     CONTINUE
ENDDO
ENDIF
IPVT(N)=N
IF(B(M,N).EQ.0.D0) IND=N
END SUBROUTINE NSBFAC

SUBROUTINE DAXPY(N,A,X,Y)
INTEGER N, I
REAL X(*),Y(*),A
DO I=1,N
Y(I)=Y(I)+A*X(I)
ENDDO
END SUBROUTINE DAXPY 

SUBROUTINE IAMAX(A,N,FINAL)
INTEGER N, I, FINAL 
REAL A(*),T
T=0.d0
FINAL=1   ! defined result even if all entries are zero
DO I=1,N
IF(ABS(A(I)).GT.T) THEN
t=abs(A(I))
FINAL=I
ENDIF
enddo
END SUBROUTINE IAMAX 

SUBROUTINE SCALE(N,T,A)
INTEGER N, I
REAL A(*),T
DO I=1,N
A(I)=T*A(I)
ENDDO
END SUBROUTINE SCALE 

SUBROUTINE NSBSLV(A,LDA,N,ML,MU,IPVT,b,x)
!--------------------------------------------------------------------
!     Solve banded linear system Ax = b
!     INPUTS:
!     A   : banded matrix as output of LU factorization by NSBFAC
!           (see storing mode in NSBFAC subroutine).
!     LDA : 1st dimension of A in calling program (lda.ge.2ml+mu+1)
!     N   : order of A                             -----------------
!     ML  : number of lower diagonals
!     MU  : number of upper diagonals
!     IPVT: integer vector of size N storing the pivoting indices
!           as output of NSBFAC.
!     b   : second member vector
!     OUTPUT:
!     x   : solution vector
!---------------------------------------------------------------------

INTEGER :: N, ML, MU, LDA
REAL :: A(LDA,*),B(*),X(*),T
INTEGER IPVT(*)
INTEGER :: I, K, L, M, LM, LA, LB, n_m_one 
DO I=1,N
X(I)=B(I)
ENDDO
M=ML+MU+1
n_m_one=N-1
!     solve L*y = b
IF(ML.NE.0.AND.n_m_one.GE.1) THEN
DO K=1,n_m_one
LM=MIN(ML,N-K)
L=IPVT(K)
T=X(L)
IF(L.NE.K) THEN
X(L)=X(K)
X(K)=T
ENDIF
CALL DAXPY(LM,T,A(M+1,K),X(K+1))
ENDDO
ENDIF
!     solve U*y = x
DO K=N,1,-1
X(K)=X(K)/A(M,K)
LM=MIN(K,M)-1
LA=M-LM
LB=K-LM
T=-X(K)
CALL DAXPY(LM,T,A(LA,K),X(LB))
ENDDO
END SUBROUTINE NSBSLV

end module two_stream_scatter_rad_mod