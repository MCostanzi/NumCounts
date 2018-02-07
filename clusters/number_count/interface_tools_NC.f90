MODULE interface_tools_NC
use cosmosis_modules
IMPLICIT none
    integer, parameter :: dl=8
    
    ! =============== DEFINE GLOBAL VARIABLE ==================
	!	SOME USEFUL CONSTANTS =================================
	double precision :: pi_value=3.141592653589793d0
	double precision :: onethr=1.d0/3.d0
	double precision :: c_light=2997.92458d0 ![km/s/100]
	double precision :: rho_c=2.7751428946d11 ! =3/8*100 2 /pi/G*Mpc/M_sun * 1e6 [in units M_sun / Mpc^3 at z=0]
	double precision :: delta_cr = 1.686d0
	! =========================================================
	! Global variables for integration ========================
	double precision :: rsmooth, z_red 
	!$OMP THREADPRIVATE(rsmooth, z_red)
    ! the above line is to use openmp since these are global variables,
    ! but with values which are specific for each thread
	! =========================================================
    ! Cosmological parameter ==================================
    double precision :: w_0,w_a, Omega_dm, Omega_m, Omega_k, Omega_v, Omega_nu, A_s
    ! =========================================================
    ! Array to allocate the effective area and z ==============
    double precision, dimension(:), allocatable :: z_eff_area,eff_area!,P_b_ztrue_array
    integer line_counter ! To interpolate effective area
    ! =========================================================
    double precision Delta_Omega ! Survey Area in Steradian
    double precision cos_theta_s(1) ! Apretura azimutale della windows function = 1 - Delta_Omega/2pi
    integer, parameter :: num_of_elle=150 ! Number of multipoles computed to compute the Window Function
    double precision :: P_l(num_of_elle+1) ! Where to store the Legendre Polynomials value for the window function
    double precision red1,red2
    integer ell,n_z1,n_z2
    integer, parameter :: num_of_kappas=200 ! number of k values used to compute the integral of the Window Function
    double precision kappa
    double precision :: k_array(num_of_kappas) ! where to store k-value for interpolation
    double precision, dimension(:,:,:), allocatable :: int_vol_bessel ! where to store \int dVol*Bessel
    double precision, dimension(:,:,:), allocatable :: sum_ell_int_vol_bessel ! sum over ell (\int dVol*Bessel)/Vol*K_ell for interpolation
    ! =========================================================
    ! Some useful Global variables ============================
    double precision LnLambda_min,LnLambda_max,z_min, z_max ! extrema for integration
    double precision mass, redshift ! Global variables for integration
    integer n_L1,n_L2 ! Global variables for integration
    !$OMP THREADPRIVATE(mass,redshift,z_min,z_max,LnLambda_min,LnLambda_max,n_L1,n_L2)
    !$OMP THREADPRIVATE(kappa,ell)
    !$OMP THREADPRIVATE(n_z1,n_z2,red1,red2)
    
    ! the above line is to use openmp since these are global variables,
    ! but with values which are specific for each thread
    integer task ! switch to compute different function in compute_number_count.f90
    ! task = 1 : d n(Lambda,z)/dz = int_0^inf dLn(M) M n(M,z) Phi(M) -> for usual number count
    ! task = 2 : d n(Lambda,z)/dz = int_0^inf dLn(M) M n(M,z) Phi(M) * M -> for mean mass in the richness bin
    ! task = 3 : d n(Lambda,z)/dz = int_0^inf dLn(M) M n(M,z) Phi(M) * b(M,z) -> for the sample variance
    !$OMP THREADPRIVATE(task)
    
    logical use_WL ! Flag used in compute_NC
    logical use_NC ! Flag used in compute_NC
    logical use_Plobltr ! Flag used in compute_NC
    logical full_NC_cov ! Flag used in compute_NC
    logical use_only_Poisson ! Flag used in compute_NC
    logical use_HOD_LM ! Flag used in compute_NC
    logical only_intr_scat ! Flag used in compute_NC
    logical use_powM_scatter ! Flag used in compute_NC
    logical use_wind_func ! Flag used in compute_NC
    logical LMrel_is_mean ! Flag used in compute_NC
    logical use_skewnorm ! Flag used in compute_NC
    logical use_BuzzardHMFcorr ! Flag used in compute_NC
    logical use_GaussPrior_HMF ! Flag used in compute_NC
    logical add_photoz_scat ! Flag used in compute_NC
    logical use_eff_area ! Flag used in compute_NC
    logical use_HMF_b_corr ! Flag used in compute_NC
    ! =========================================================
    ! NUISANCE PARAMETERS FOR THE HMF =========================
    double precision :: slope, inter
    ! NUISANCE PARAMETERS FOR THE RICHNESS-MASS RELATION ======
    double precision :: delta_a0,delta_b0,delta_c0,Log10Mminsat,sig_intr
    double precision :: Log10Mmax ! log10 upper mass limit for integration
    ! =========================================================
    ! NUISANCE PARAMETERS FOR THE PROJECTION EFFECTS ==========
    double precision :: frac_p,epsi
    ! NUISANCE PARAMETERS FOR THE BARYON CORRECTION ===========
    double precision :: p_a,p_c,p_d
    ! =========================================================
    ! PARAMETER TO SCALE THE COV MATRIX =======================
    double precision :: NC_COV_scale,WL_COV_scale
    ! =========================================================


    type ini_settings
        integer :: feedback
        double precision :: Delta_mean
        double precision :: Delta_crit
        integer :: num_redshift_bin
        integer :: num_lambda_bin
        logical :: use_WL_mass
        logical :: use_NC_data
        character(len=300) :: file_data_NC
        character(len=300) :: file_cov_misc
        character(len=300) :: file_data_WL
        character(len=300) :: file_cov_WL
        character(len=300) :: file_cov_HMFnuis
        character(len=300) :: file_dlogMdOm
        character(len=10) :: PDF_lob_lin_version
    end type

    type number_counts
        double precision, dimension(:,:), allocatable :: n_Li_zj
        double precision, dimension(:,:), allocatable :: n_Li_zj_var
        double precision, dimension(:,:), allocatable :: n_Li_zj_data
        double precision, dimension(:,:), allocatable :: mean_m_Li_zj
        double precision, dimension(:,:), allocatable :: mean_m_Li_zj_data
        double precision, dimension(:,:), allocatable :: cov_misc_i_j
        double precision, dimension(:), allocatable :: data_i
        double precision, dimension(:,:,:,:), allocatable :: cov_Li_zi_Lj_zj
        double precision, dimension(:,:), allocatable :: cov_i_j
        double precision, dimension(:), allocatable :: LnLambda_min_array
        double precision, dimension(:), allocatable :: LnLambda_max_array
        double precision, dimension(:), allocatable :: z_min_array
        double precision, dimension(:), allocatable :: z_max_array
        double precision, dimension(:,:), allocatable :: dlogMdOmega
        integer num_of_L_bin,num_of_z_bin
        double precision sigma8_NC
    end type


    type pk_settings
        real(dl), dimension(:), allocatable :: redshifts
        real(dl), dimension(:), allocatable :: kh
        real(dl), dimension(:,:), allocatable :: matpower
        real(dl), dimension(:,:), allocatable :: ddmat
        integer num_z, num_k
    end type

    contains

function load_matter_power(block, PK) result(status)
    use cosmosis_modules
    integer(cosmosis_block) :: block
    integer(cosmosis_status) :: status
    type(pk_settings) :: PK
    real(dl), allocatable, dimension(:) :: k, z
    real(dl), allocatable, dimension(:,:) :: P

    !Get the data columns from the fits data
    status = 0

    !Load k, z, P
    status = datablock_get_double_grid(block, matter_power_lin_section, &
    "K_H", k, "Z", z, "P_K", P)

    if (status .ne. 0) then
        write(*,*) "Could not find K_H, Z, or P_K in block"
        return
    endif

    !Fill in data structure
    PK%num_k = size(k)
    PK%num_z = size(z)
    call allocate_matterpower(PK)
    PK%kh = k
    PK%redshifts = z
    PK%matpower = P


    !Clean up
    deallocate(k, z, P)
end function


subroutine allocate_matterpower(PK)
    type(pk_settings) :: PK
    allocate(PK%redshifts(PK%num_z))
    allocate(PK%kh(PK%num_k))
    allocate(PK%matpower(PK%num_k,PK%num_z))
    allocate(PK%ddmat(PK%num_k,PK%num_z))
end subroutine

subroutine deallocate_matterpower(PK)
    type(pk_settings) :: PK
    deallocate(PK%redshifts)
    deallocate(PK%kh)
    deallocate(PK%matpower)
    deallocate(PK%ddmat)
end subroutine

subroutine allocate_NC(NC,settings)
    type(number_counts) :: NC
    type(ini_settings) :: settings
    allocate(NC%n_Li_zj(NC%num_of_L_bin,NC%num_of_z_bin))
    allocate(NC%n_Li_zj_var(NC%num_of_L_bin,NC%num_of_z_bin))
    allocate(NC%cov_misc_i_j(NC%num_of_L_bin*NC%num_of_z_bin,NC%num_of_L_bin*NC%num_of_z_bin))
    allocate(NC%n_Li_zj_data(NC%num_of_L_bin,NC%num_of_z_bin))
    allocate(NC%mean_m_Li_zj_data(NC%num_of_L_bin,NC%num_of_z_bin))
    allocate(NC%mean_m_Li_zj(NC%num_of_L_bin,NC%num_of_z_bin))
    allocate(NC%dlogMdOmega(NC%num_of_L_bin,NC%num_of_z_bin))
    if (settings%use_WL_mass .and. settings%use_NC_data) then
        allocate(NC%data_i(NC%num_of_L_bin*NC%num_of_z_bin*2))
        allocate(NC%cov_i_j(NC%num_of_L_bin*NC%num_of_z_bin*2,NC%num_of_L_bin*NC%num_of_z_bin*2))
    else
        allocate(NC%data_i(NC%num_of_L_bin*NC%num_of_z_bin))
        allocate(NC%cov_i_j(NC%num_of_L_bin*NC%num_of_z_bin,NC%num_of_L_bin*NC%num_of_z_bin))
    end if
    allocate(NC%cov_Li_zi_Lj_zj(NC%num_of_L_bin,NC%num_of_z_bin,NC%num_of_L_bin,NC%num_of_z_bin))
    allocate(NC%LnLambda_min_array(NC%num_of_L_bin))
    allocate(NC%LnLambda_max_array(NC%num_of_L_bin))
    allocate(NC%z_min_array(NC%num_of_z_bin))
    allocate(NC%z_max_array(NC%num_of_z_bin))
end subroutine

subroutine deallocate_NC(NC)
    type(number_counts) :: NC
    deallocate(NC%n_Li_zj)
    deallocate(NC%n_Li_zj_var)
    deallocate(NC%cov_misc_i_j)
    deallocate(NC%n_Li_zj_data)
    deallocate(NC%mean_m_Li_zj)
    deallocate(NC%mean_m_Li_zj_data)
    deallocate(NC%data_i)
    deallocate(NC%cov_Li_zi_Lj_zj)
    deallocate(NC%cov_i_j)
    deallocate(NC%LnLambda_min_array)
    deallocate(NC%LnLambda_max_array)
    deallocate(NC%z_min_array)
    deallocate(NC%z_max_array)
    deallocate(NC%dlogMdOmega)
end subroutine



END MODULE interface_tools_NC
