
function setup(options) result(result)
	USE cosmosis_modules
	USE interface_tools_NC
	USE read_data
	implicit none
	integer(cosmosis_block), value :: options
	integer(cosmosis_status) :: status
	type(ini_settings), pointer :: settings
	type(c_ptr) :: result
	allocate(settings)
	settings%feedback = 0
	settings%Delta_mean = 0.d0
	settings%Delta_crit = 0.d0

	status = 0
	status = status + datablock_get_int_default(options, option_section, "feedback", 0, settings%feedback)
	status = status + datablock_get_double_default(options, option_section, "Delta_mean", 0.d0, settings%Delta_mean)
	status = status + datablock_get_double_default(options, option_section, "Delta_crit", 0.d0, settings%Delta_crit)
	status = status + datablock_get_double_default(options, option_section, "Delta_Omg", 3.1415d0, Delta_Omega)
	status = status + datablock_get_double_default(options, option_section, "Log10Mmax", 16.6d0, Log10Mmax)
	status = status + datablock_get_int_default(options, option_section, "num_redshift_bin", 1, settings%num_redshift_bin)
	status = status + datablock_get_int_default(options, option_section, "num_lambda_bin", 1, settings%num_lambda_bin)
	status = status + datablock_get_string(options, option_section, "PDF_lob_lin_version", settings%PDF_lob_lin_version)
	status = status + datablock_get_logical_default(options, option_section, "use_NC_data", .false., settings%use_NC_data)
	status = status + datablock_get_logical_default(options, option_section, "use_WL_mass", .false., settings%use_WL_mass)
	status = status + datablock_get_logical_default(options, option_section, "use_P_lob_ltr", .true., use_Plobltr)
	status = status + datablock_get_logical_default(options, option_section, "use_HOD_LM_rel", .true., use_HOD_LM)
	status = status + datablock_get_logical_default(options, option_section, "use_only_intr_scat", .true., only_intr_scat)
	status = status + datablock_get_logical_default(options, option_section, "use_powlaw_mass_scat", .true., use_powM_scatter)
	status = status + datablock_get_logical_default(options, option_section, "LMrel_is_mean", .true., LMrel_is_mean)
	status = status + datablock_get_logical_default(options, option_section, "use_skewnorm", .true., use_skewnorm)
	status = status + datablock_get_logical_default(options, option_section, "use_BuzzardHMFcorr", .false., use_BuzzardHMFcorr)
	status = status + datablock_get_logical_default(options, option_section, "add_photoz_scat", .true., add_photoz_scat)
	status = status + datablock_get_logical_default(options, option_section, "use_eff_area", .true., use_eff_area)
	status = status + datablock_get_logical_default(options, option_section, "use_HMF_b_corr", .false., use_HMF_b_corr)
	status = status + datablock_get_string(options, option_section, "data_file_bin", settings%file_data_NC)

	full_NC_cov= .false.
	use_only_Poisson= .true.
	use_wind_func= .false.


	if (status .ne. 0) then 
		write(*,*) "Failed setup of Number Counts!", status
		stop
	endif
	result = c_loc(settings)
	if (settings%Delta_mean > 1.d-3 .and. settings%Delta_crit > 1.d-3 .or. settings%Delta_mean < 1.d-3 .and. settings%Delta_crit < 1.d-3) then
		write(*,*) "Warning: both Delta_mean and Delta_crit larger or equal to zero;"
		write(*,*) "Define just one according to your mass definition"
		stop
	end if
	if ( .not. settings%use_NC_data .and. .not. settings%use_WL_mass) then
		write(*,*) "Warning: you have to use at least one data set;"
		write(*,*) "Set to TRUE one of the two flgas use_NC_data / use_WL_mass"
		stop
	end if
	
	! Read the Fitting Parameters for the Tinker HMF ========
	call init_tinker_param
	! READ TABLES FOR INT P(L_ob|L_tr) 
	if (use_Plobltr) call read_int_P_lob_ltr(settings%num_lambda_bin,settings%PDF_lob_lin_version)
	! READ TABLES FOR SKEW-GAUSS PARAMETERS (<lambda_sat>,sig_intr) FOR P(l_tr|M)
	if (use_skewnorm) call read_skew_gauss_table()


end function setup

function execute(block, config) result(status)

	use interface_tools_NC
	use cosmosis_modules
	use compute_number_count
	implicit none
	integer(cosmosis_block), value :: block
	integer(cosmosis_status) :: status
	type(c_ptr), value :: config
	type(ini_settings), pointer :: settings
	type(pk_settings) :: PK
	type(number_counts) :: NC

	integer nz,nL,nz2,nL2,index_i,index_j
	double precision total_lnlike
	double precision M1,sig_intr_Mmax,sig_intr_Mmin,Mpivot
	logical compute_like ! Logical parameter to check model priors
	! HMF Nuisance parameters Gauss priors
	double precision mean_slope,mean_inter,sig_slope,sig_inter 
	! Richness-mass relation parameters Gauss priors
	double precision mean_delta_a0,mean_delta_b0,mean_delta_c0,sig_delta_a0,sig_delta_b0,sig_delta_c0,mean_sig_intr,sig_sig_intr,mean_epsi,sig_epsi
	! Pure fraction parameter Gauss priors
	double precision mean_p_frac,sig_p_frac ! OLD MODEL; NOT USED ANYMORE
	double precision, dimension(2,2) :: cov_hmf
	double precision, dimension(2) :: data_hmf
	double precision lnlike_hmf

	status = 0
	call c_f_pointer(config, settings)


	! Get Cosmological Parameter Values
	status = status + datablock_get_double(block, cosmological_parameters_section, "OMEGA_M", Omega_dm) ! NOTE "Omega_m" here is actually Omega_b+Omega_cdm
	status = status + datablock_get_double(block, cosmological_parameters_section, "OMEGA_LAMBDA", Omega_v)
	status = status + datablock_get_double_default(block,cosmological_parameters_section, "OMEGA_NU", 0.0D0, Omega_nu)
	status = status + datablock_get_double_default(block, cosmological_parameters_section, "OMEGA_K", 0.0D0, Omega_k)
	status = status + datablock_get_double_default(block, cosmological_parameters_section, "W", -1.0D0, w_0)
	status = status + datablock_get_double_default(block, cosmological_parameters_section, "WA", 0.0D0, w_a)
	status = status + datablock_get_double_default(block, cosmological_parameters_section, "A_s", 0.0D0, A_s)
	if (A_s <= 0.d0) then
		    status = status + datablock_get_double(block, cosmological_parameters_section, "LnA_s",     A_s)
		    A_s=dexp(A_s)/1.0d10
	end if
	! Get Nuisance Parameters for HMF
	status = status + datablock_get_double_default(block, number_count_section,"HMF_SLOPE", 0.0D0, slope)
	status = status + datablock_get_double_default(block, number_count_section,"HMF_INTER", 1.0D0, inter)
	! Get Nuisance Parameters for Richness-mass relation
	status = status + datablock_get_double_default(block, number_count_section,"RM_SLOPE", 0.0D0, delta_a0)
	status = status + datablock_get_double_default(block, number_count_section,"RM_INTER", 0.0D0, delta_b0)
	status = status + datablock_get_double_default(block, number_count_section,"RM_Z_EVOL", 0.0D0, delta_c0)
	status = status + datablock_get_double_default(block, number_count_section,"LOG10MMINSAT", 5.0D12, Log10Mminsat)
	status = status + datablock_get_double_default(block, number_count_section,"RM_INTR_SIG", 0.20D0,sig_intr)
	status = status + datablock_get_double_default(block, number_count_section,"EPSILON", 0.0D0, epsi)
	! Get Nuisance Parameters for Projection effects
	status = status + datablock_get_double_default(block, number_count_section,"P_FRACT", 1.0D0, frac_p) ! OLD NOT USED
	! Get Nuisance Parameters for HMF baryon correction
	status = status + datablock_get_double_default(block, number_count_section,"P_A", -0.0872D0, p_a)
	status = status + datablock_get_double_default(block, number_count_section,"P_C", 13.6339D0, p_c)
	status = status + datablock_get_double_default(block, number_count_section,"P_D", 0.3509D0, p_d)


	! Re-define parameters according to the usal definition
	Omega_m = Omega_dm + Omega_nu ! Omega_m = Omega_cdm + Omega_baryon + Omega_nu

	!load in the matter power spectrum
	status = load_matter_power(block,PK)
	if (status .ne. 0) then
		write(*,*) "Could not load matter power"
		status=3
		return 
	endif

	use_WL=settings%use_WL_mass ! Global variable used in compute_NC
	use_NC=settings%use_NC_data ! Global variable used in compute_NC


!	! SET THE NUMBER OF LAMBDA AND REDSHIFT BIN TO COMPUTE
	NC%num_of_L_bin= settings%num_lambda_bin
	NC%num_of_z_bin= settings%num_redshift_bin
	call allocate_NC(NC,settings)

	call read_bins(NC,settings%file_data_NC)

    ! =========================================================
	Delta_c=settings%Delta_crit ! Set Delta_crit > 0 in the .ini file for halo masses defined wrt the critical density
	Delta_m=settings%Delta_mean ! Set Delta_mean > 0 in the .ini file for halo masses defined wrt the mean background density
	Delta_m_omp=Delta_m ! Variable set to avid problem with OMP

	if (NC%z_max_array(NC%num_of_z_bin)>PK%redshifts(PK%num_z)) then
	    write(*,*) "WARNING out of redshift range computed by CAMB; increse z_max" 
	    return
	end if

	call compute_NC(PK,NC) ! compute n(Lambda_i,z_j))

    if (settings%feedback>=1) then
    write(*,*) 'sigma8',NC%sigma8_NC,'Omega_m',Omega_m,'LnA_s*10^10',log(A_s*1.d10),'Omega_nu',Omega_nu
    write(*,*) 'HMF slope',slope,'HMF amp',inter
    write(*,*) 'RMR SLOPE',delta_a0,'RMR LOG10MPIVOT/AMP',delta_b0
    write(*,*) 'Log10Mminsat',Log10Mminsat,'sig_intr',sig_intr
    write(*,*) 'scatter EVOL',delta_c0,'z EVOL',epsi
    if (use_HMF_b_corr) write(*,*) 'p_a',p_a,'p_c',p_c,'p_d',p_d
    write(*,*) '==================================='
    write(*,*) 'Z LOW Z HI Lambda LOW Lambda HI NC THEO   LogM_WL THEO'
	do nz=1,NC%num_of_z_bin ! Start loop over redshift bin
	    do nL=1,NC%num_of_L_bin ! Start loop over Lambda bin
	        index_i=(nz-1)*NC%num_of_L_bin+nL
	        write(*,'(6f16.4)') NC%z_min_array(nz),NC%z_max_array(nz),exp(NC%LnLambda_min_array(nL)),exp(NC%LnLambda_max_array(nL)),NC%n_Li_zj(nL,nz),NC%mean_m_Li_zj(nL,nz)
	    end do
	end do
	end if


	! NB: number_count_section has to be added to ~/cosmosis/datablock/cosmosis_f90/cosmosis_section_names.f90
	! Add manually the following line at the end of cosmosis/datablock/section_names.txt
	!    #number_count
	!    number_count
	status = datablock_put_double_grid(block,number_count_section,"Lambda",NC%LnLambda_min_array,"z",NC%z_min_array, "n_Lambda_z",NC%n_Li_zj)


	call deallocate_matterpower(PK)
	call deallocate_NC(NC)

end function


function cleanup(config) result(status)
	use interface_tools_NC
	use cosmosis_modules
	type(c_ptr), value :: config
	type(ini_settings), pointer :: settings
	integer(cosmosis_status) :: status

	!Free memory allocated in the setup function
	call c_f_pointer(config, settings)
	deallocate(settings)

	status = 0

end function cleanup

