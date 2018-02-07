
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
	if (settings%use_NC_data) status = status + datablock_get_string(options, option_section, "data_file_NC", settings%file_data_NC)
	if (settings%use_WL_mass) status = status + datablock_get_string(options, option_section, "data_file_WL", settings%file_data_WL)
	if (settings%use_WL_mass) status = status + datablock_get_string(options, option_section, "data_file_COV_WL", settings%file_cov_WL)
	if (settings%use_WL_mass) status = status + datablock_get_string(options, option_section, "dlogMdOm_file", settings%file_dlogMdOm)
	if (settings%use_NC_data) status = status + datablock_get_string(options, option_section, "data_file_COV_MISC", settings%file_cov_misc)
	if (settings%use_NC_data) status = status + datablock_get_double_default(options, option_section, "scale_NC_COV", 1.d0, NC_COV_scale)
	if (settings%use_WL_mass) status = status + datablock_get_double_default(options, option_section, "scale_WL_COV", 1.d0, WL_COV_scale)
	status = status + datablock_get_logical_default(options, option_section, "use_P_lob_ltr", .true., use_Plobltr)
	status = status + datablock_get_logical_default(options, option_section, "use_full_NC_cov", .true., full_NC_cov)
	status = status + datablock_get_logical_default(options, option_section, "use_only_Poisson_cov", .true., use_only_Poisson)
	status = status + datablock_get_logical_default(options, option_section, "use_HOD_LM_rel", .true., use_HOD_LM)
	status = status + datablock_get_logical_default(options, option_section, "use_only_intr_scat", .true., only_intr_scat)
	status = status + datablock_get_logical_default(options, option_section, "use_powlaw_mass_scat", .true., use_powM_scatter)
	status = status + datablock_get_logical_default(options, option_section, "use_wind_func", .true., use_wind_func)
	status = status + datablock_get_logical_default(options, option_section, "LMrel_is_mean", .true., LMrel_is_mean)
	status = status + datablock_get_logical_default(options, option_section, "use_skewnorm", .true., use_skewnorm)
	status = status + datablock_get_logical_default(options, option_section, "use_BuzzardHMFcorr", .false., use_BuzzardHMFcorr)
	status = status + datablock_get_logical_default(options, option_section, "use_GaussPrior_HMF", .true., use_GaussPrior_HMF)
	status = status + datablock_get_logical_default(options, option_section, "add_photoz_scat", .true., add_photoz_scat)
	if (use_GaussPrior_HMF) status = status + datablock_get_string(options, option_section, "data_file_COV_HMFnuis", settings%file_cov_HMFnuis)
	status = status + datablock_get_logical_default(options, option_section, "use_eff_area", .true., use_eff_area)
	status = status + datablock_get_logical_default(options, option_section, "use_HMF_b_corr", .false., use_HMF_b_corr)


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
	use MatrixUtils ! perche' venga compilato devo modificare il makefile di camb, e uncomment the EXTQLCS = MatrixUtils line 
	! Forse devo fare la gabola poi di ricommentare la linea del makefile di camb
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
	! Get parameters for Gaussian priors for Nuisance parameters
	status = status + datablock_get_double_default(block, number_count_section,"MEAN_RM_SLOPE", 0.0D0, mean_delta_a0)
	status = status + datablock_get_double_default(block, number_count_section,"MEAN_RM_INTER", 0.0D0, mean_delta_b0)
	status = status + datablock_get_double_default(block, number_count_section,"MEAN_RM_Z_EVOL", 0.0D0, mean_delta_c0)
	status = status + datablock_get_double_default(block, number_count_section,"MEAN_RM_INTR_SIG", 0.0D0,mean_sig_intr)
	status = status + datablock_get_double_default(block, number_count_section,"MEAN_HMF_SLOPE", 0.0D0, mean_slope)
	status = status + datablock_get_double_default(block, number_count_section,"MEAN_HMF_INTER", 0.0D0, mean_inter)
	status = status + datablock_get_double_default(block, number_count_section,"MEAN_P_FRACT", 0.0D0,mean_p_frac)  ! OLD NOT USED
	status = status + datablock_get_double_default(block, number_count_section,"MEAN_EPSILON", 0.0D0, mean_epsi)
	status = status + datablock_get_double_default(block, number_count_section,"SIG_RM_SLOPE", 0.0D0, sig_delta_a0)
	status = status + datablock_get_double_default(block, number_count_section,"SIG_RM_INTER", 0.0D0, sig_delta_b0)
	status = status + datablock_get_double_default(block, number_count_section,"SIG_RM_Z_EVOL", 0.0D0, sig_delta_c0)
	status = status + datablock_get_double_default(block, number_count_section,"SIG_RM_INTR_SIG", 0.0D0,sig_sig_intr)
	status = status + datablock_get_double_default(block, number_count_section,"SIG_HMF_SLOPE", 0.0D0, sig_slope)
	status = status + datablock_get_double_default(block, number_count_section,"SIG_HMF_INTER", 0.0D0, sig_inter)
	status = status + datablock_get_double_default(block, number_count_section,"SIG_P_FRACT", 0.0D0, sig_p_frac)  ! OLD NOT USED
	status = status + datablock_get_double_default(block, number_count_section,"SIG_EPSILON", 0.0D0, sig_epsi)

	! Re-define parameters according to the usal definition
	Omega_m = Omega_dm + Omega_nu ! Omega_m = Omega_cdm + Omega_baryon + Omega_nu

    compute_like= .True. ! Initialize flag

    ! CHECK PRIORS RANGE ======================================
    ! Check sigma_intr priors
    if (use_HOD_LM) then
        Mpivot=10.d0**delta_b0
    else ! use power law lmabda-mass relation
        Mpivot=10.d0**14.344d0
    end if
    if (use_powM_scatter) then
        if (delta_c0<=0.d0) then ! If slope is negative
            sig_intr_Mmax=sig_intr*(10.**Log10Mminsat/Mpivot)**delta_c0
            sig_intr_Mmin=sig_intr*(10.**Log10Mmax/Mpivot)**delta_c0
        else ! If slope is positive
            sig_intr_Mmin=sig_intr*(10.**Log10Mminsat/Mpivot)**delta_c0
            sig_intr_Mmax=sig_intr*(10.**Log10Mmax/Mpivot)**delta_c0
        end if
    else
        if (delta_c0<=0.d0) then ! If slope is negative
            sig_intr_Mmax=sig_intr+log(10.**Log10Mminsat/Mpivot)*delta_c0
            sig_intr_Mmin=sig_intr+log(10.**Log10Mmax/Mpivot)*delta_c0
        else ! If slope is positive
            sig_intr_Mmin=sig_intr+log(10.**Log10Mminsat/Mpivot)*delta_c0
            sig_intr_Mmax=sig_intr+log(10.**Log10Mmax/Mpivot)*delta_c0
        end if
    end if
!    write(*,*) sig_intr_Mmin,sig_intr_Mmax
    if ( sig_intr_Mmin>0.05d0 .and. sig_intr_Mmax<2.0d0) then
            compute_like= .True.  ! CHECK sigma_intr CONDITION
    else
            compute_like= .False.
    end if

    ! If use_HOD_LM_rel check the M1 priors
    if( use_HOD_LM .and. compute_like) then
        M1=10.d0**Log10Mminsat + 10.d0**delta_b0 ! M1= Mass at which <l_sat|M1> = 1
        if (M1 >= 10.d0*(10.d0**Log10Mminsat) .and. M1<= 30.d0*(10.d0**Log10Mminsat)) then
            compute_like= .True.
        else
            compute_like= .False.
        end if
    end if

    ! If use hard prior on HMF nuisance parameters check them
    if (.not. use_GaussPrior_HMF .and. compute_like) then 
        if (slope>0.017d0 .and. slope<0.08d0 .and. inter>0.935d0+1.5d0*slope .and. inter<0.98d0+1.5d0*slope) then
            compute_like= .True.
        else
            compute_like= .False.
        end if
    end if
    ! END CHECK PRIORS RANGE ==================================

    if (compute_like) then

	!load in the matter power spectrum
	status = load_matter_power(block,PK)
	if (status .ne. 0) then
		write(*,*) "Could not load matter power"
		status=3
		return 
	endif



!	! SET THE NUMBER OF LAMBDA AND REDSHIFT BIN TO COMPUTE
	NC%num_of_L_bin= settings%num_lambda_bin
	NC%num_of_z_bin= settings%num_redshift_bin
	call allocate_NC(NC,settings)


    ! =========================================================
	Delta_c=settings%Delta_crit ! Set Delta_crit > 0 in the .ini file for halo masses defined wrt the critical density
	Delta_m=settings%Delta_mean ! Set Delta_mean > 0 in the .ini file for halo masses defined wrt the mean background density
	Delta_m_omp=Delta_m ! Variable set to avid problem with OMP


	NC%cov_i_j=0.d0
	NC%data_i=0.d0
	total_lnlike=0.d0

	use_WL=settings%use_WL_mass ! Global variable used in compute_NC
	use_NC=settings%use_NC_data ! Global variable used in compute_NC
	if (settings%use_NC_data) call read_NC_data(NC,settings%file_data_NC) ! read NC data
	if (settings%use_NC_data) call read_COV_NCMISS_data(NC,settings%file_cov_misc) ! read COV NC MISCENTERING
	if (settings%use_WL_mass) call read_WL_data(NC,settings%file_data_WL) ! read WL-masses(Omega_m=0.30) data
	if (settings%use_WL_mass) call read_COV_WL_data(NC,settings%file_cov_WL) ! read COV of WL data
	if (NC%z_max_array(NC%num_of_z_bin)>PK%redshifts(PK%num_z)) then
	    write(*,*) "WARNING out of redshift range computed by CAMB; increse z_max" 
	    return
	end if

	call compute_NC(PK,NC) ! compute n(Lambda_i,z_j))
	if (settings%use_WL_mass) then
	! Reascle Log10M(Omega_m=0.3) for Omega_m =================
	call read_dlogMdOm(NC,settings%file_dlogMdOm,NC%dlogMdOmega)
	do nz=1,NC%num_of_z_bin
	    do nL=1,NC%num_of_L_bin
	        NC%mean_m_Li_zj_data(nL,nz)=NC%mean_m_Li_zj_data(nL,nz)+NC%dlogMdOmega(nL,nz)*(Omega_m-0.3d0)
!	        write(*,*) NC%dlogMdOmega(nL,nz)
	    end do
	end do
	end if
    ! =========================================================


!	! Sorting the covariance and data array conviniently to compute loglike
	do nz=1,NC%num_of_z_bin ! Start loop over redshift bin
	    do nL=1,NC%num_of_L_bin ! Start loop over Lambda bin
	        index_i=(nz-1)*NC%num_of_L_bin+nL
	        if (settings%use_NC_data) NC%data_i(index_i)= NC%n_Li_zj(nL,nz)-NC%n_Li_zj_data(nL,nz) ! NC data
	        if (.not. settings%use_NC_data .and. settings%use_WL_mass) NC%data_i(index_i)= NC%mean_m_Li_zj(nL,nz)-NC%mean_m_Li_zj_data(nL,nz) ! WL data
	        if (settings%use_WL_mass .and. settings%use_NC_data) NC%data_i(index_i+NC%num_of_L_bin*NC%num_of_z_bin)= NC%mean_m_Li_zj(nL,nz)-NC%mean_m_Li_zj_data(nL,nz) ! WL data
	        do nz2=1,NC%num_of_z_bin
	            do nL2=1,NC%num_of_L_bin
	                index_j=(nz2-1)*NC%num_of_L_bin+nL2
                    if (settings%use_NC_data) then
	                    if (index_j>=index_i) NC%cov_i_j(index_i,index_j)= (NC%cov_Li_zi_Lj_zj(nL,nz,nL2,nz2)+NC%cov_misc_i_j(index_i,index_j))*NC_COV_scale ! add miscentering component and scale COV
	                    if (index_j<index_i)  NC%cov_i_j(index_i,index_j)= NC%cov_i_j(index_j,index_i) ! matrice simmetrica
	                end if
	            end do
	        end do
	    end do
	end do

    if (settings%feedback>=1) then
    write(*,*) 'sigma8',NC%sigma8_NC,'Omega_m',Omega_m,'LnA_s*10^10',log(A_s*1.d10),'Omega_nu',Omega_nu
    write(*,*) 'HMF slope',slope,'HMF amp',inter
    write(*,*) 'RMR SLOPE',delta_a0,'RMR LOG10MPIVOT/AMP',delta_b0
    write(*,*) 'Log10Mminsat',Log10Mminsat,'sig_intr',sig_intr
    write(*,*) 'scatter EVOL',delta_c0,'z EVOL',epsi
    if (use_HMF_b_corr) write(*,*) 'p_a',p_a,'p_c',p_c,'p_d',p_d
    write(*,*) '==================================='
    write(*,*) 'NC THEO    NC DATA   DIFF(NC)   LogM_WL THEO    LogM_WL DATA   DIFF(M_WL)'
	do nz=1,NC%num_of_z_bin ! Start loop over redshift bin
	    do nL=1,NC%num_of_L_bin ! Start loop over Lambda bin
	        index_i=(nz-1)*NC%num_of_L_bin+nL
	        write(*,*) NC%n_Li_zj(nL,nz),NC%n_Li_zj_data(nL,nz),NC%n_Li_zj(nL,nz)-NC%n_Li_zj_data(nL,nz),NC%mean_m_Li_zj(nL,nz),NC%mean_m_Li_zj_data(nL,nz),NC%mean_m_Li_zj(nL,nz)-NC%mean_m_Li_zj_data(nL,nz)!NC%cov_i_j(index_i+NC%num_of_L_bin*NC%num_of_z_bin,index_i+NC%num_of_L_bin*NC%num_of_z_bin)
	    end do
	end do
	end if


	total_lnlike= - Matrix_GaussianLogLikeDouble(NC%cov_i_j,NC%data_i) ! non considera il termine -0.5ln(2*pi)

	! Gaussian prior on the nuisance parameter for the HMF ====
	! From the comparison of Buzzard HMF / Tinker HMF =========
	! computed for the fiducial values ========================
	if (use_GaussPrior_HMF) then
	    call read_COV_HMF_nuis(settings%file_cov_HMFnuis,cov_hmf)
	    data_hmf=[slope-mean_slope,inter-mean_inter] ! My Cov
	    lnlike_hmf= -Matrix_GaussianLogLikeDouble(cov_hmf ,data_hmf)
	    total_lnlike=total_lnlike+lnlike_hmf
	end if
	! =========================================================

	! Gaussian prior on Richness-mass relation parameters =====

	total_lnlike=total_lnlike+LnGauss_prior(delta_a0,mean_delta_a0,sig_delta_a0)+LnGauss_prior(delta_b0,mean_delta_b0,sig_delta_b0)+LnGauss_prior(sig_intr,mean_sig_intr,sig_sig_intr)
	! =========================================================


!	! Add here other Gaussian prior if needed =================
!	total_lnlike=total_lnlike+LnGauss_prior(new_param,mean_new_param,sig_new_param)
!	! =========================================================


	status = datablock_put_double(block, likelihoods_section, "NC_LIKE", total_lnlike)
	status = status + datablock_put_double(block,number_count_section,"sigma8_NC",NC%sigma8_NC)

	if(settings%feedback >0) then
	    write(*,*) 'NC_LNLIKE: ',total_lnlike
	end if

	call deallocate_matterpower(PK)
	call deallocate_NC(NC)

	else ! IF PARAMETERS DOES NOT SATISFY THE PRIOR CONDITION
	total_lnlike=-1.0d50
	status = datablock_put_double(block, likelihoods_section, "NC_LIKE", total_lnlike)
	status = status + datablock_put_double(block,number_count_section,"sigma8_NC",1000.d0)
	if(settings%feedback >0) then
	    write(*,*) 'NC_LNLIKE: ',total_lnlike
	end if
	end if

end function execute


function cleanup(config) result(status)
	use interface_tools_NC
	use cosmosis_modules
	USE utils_NC
	type(c_ptr), value :: config
	type(ini_settings), pointer :: settings
	integer(cosmosis_status) :: status

	!Free memory allocated in the setup function
	call c_f_pointer(config, settings)
	deallocate(settings)

	if (use_Plobltr) deallocate(int_P_lob_ltr)
	if (use_Plobltr) deallocate(z_int_P_lob)
	if (use_Plobltr) deallocate(l_int_P_lob)
	if (use_skewnorm) deallocate(skew_table)
	if (use_skewnorm) deallocate(sig_skew_table)
	if (use_skewnorm) deallocate(l_sat_grid)
	if (use_skewnorm) deallocate(sig_intr_grid)
	if (use_skewnorm) deallocate(ddskew_table)
	if (use_skewnorm) deallocate(ddsig_skew_table)

	status = 0

end function cleanup
