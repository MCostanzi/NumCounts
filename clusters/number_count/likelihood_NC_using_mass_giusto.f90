
function setup(options) result(result)
	USE cosmosis_modules
	USE interface_tools_NC
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
	status = status + datablock_get_double_default(options, option_section, "Delta_Omg", 3.1415d0, settings%Delta_Omg)
	status = status + datablock_get_int_default(options, option_section, "num_redshift_bin", 1, settings%num_redshift_bin)
	status = status + datablock_get_int_default(options, option_section, "num_lambda_bin", 1, settings%num_lambda_bin)
	status = status + datablock_get_double_default(options, option_section, "minimum_redshift", 0.d0, settings%minimum_redshift)
	status = status + datablock_get_double_default(options, option_section, "maximum_redshift", 1.d0, settings%maximum_redshift)
	status = status + datablock_get_double_default(options, option_section, "minimum_lambda", 5.d0, settings%minimum_lambda)
	status = status + datablock_get_double_default(options, option_section, "maximum_lambda", 125.d0, settings%maximum_lambda)
	status = status + datablock_get_logical_default(options, option_section, "use_NC_data", .false., settings%use_NC_data)
	if (settings%use_NC_data) status = status + datablock_get_string(options, option_section, "data_file_NC", settings%file_data_NC)
	status = status + datablock_get_logical_default(options, option_section, "use_WL_mass", .false., settings%use_WL_mass)
	if (settings%use_WL_mass) status = status + datablock_get_string(options, option_section, "data_file_WL", settings%file_data_WL)
	if (settings%use_WL_mass) status = status + datablock_get_string(options, option_section, "data_file_COV_WL", settings%file_cov_WL)
	status = status + datablock_get_logical_default(options, option_section, "blind_result", .false., settings%blind)
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
end function setup

function execute(block, config) result(status)

	use interface_tools_NC
	use cosmosis_modules
	use compute_number_count
	use compute_Cl_cluster
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
	double precision lnlike,total_lnlike
	! HMF Nuisance Likelihood parameters
	double precision mean_slope,mean_inter,sig_slope,sig_inter 
	! Ln(Lambda(M,z)) Nuisance Likelihood parameters
	double precision mean_delta_a0,mean_delta_b0,mean_delta_c0,sig_delta_a0,sig_delta_b0,sig_delta_c0,mean_sig_intr,sig_sig_intr
	! Ln(Lambda(M,z)) Nuisance Likelihood parameters
	double precision mean_p_frac,mean_epsi,sig_p_frac,sig_epsi

	! Parameter for blinding result
	! delta_blind = displacements for cosmological parameters
	! delta_blind = [dOmega_m, dOmega_b, dOmega_k, dA_s, dn_s, dh_fid, d_m_nu, dw, dw_a]
    integer,         parameter                    :: REAL_KIND = 4
    integer,         parameter                    :: UNIT = 10
    integer,         parameter                    :: SAMPLE_LENGTH = 9
    real(REAL_KIND), dimension(0:SAMPLE_LENGTH-1) :: delta_blind
    integer                                       :: i_blind


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
	! Get Nuisance Parameters for HMF
	status = status + datablock_get_double_default(block, number_count_section,"HMF_SLOPE", 0.0D0, slope)
	status = status + datablock_get_double_default(block, number_count_section,"HMF_INTER", 1.0D0, inter)
	! Get Nuisance Parameters for Richness-mass relation
	status = status + datablock_get_double_default(block, number_count_section,"RM_SLOPE", 0.0D0, delta_a0)
	status = status + datablock_get_double_default(block, number_count_section,"RM_INTER", 0.0D0, delta_b0)
	status = status + datablock_get_double_default(block, number_count_section,"RM_Z_EVOL", 0.0D0, delta_c0)
	status = status + datablock_get_double_default(block, number_count_section,"RM_INTR_SIG", 0.20D0,sig_intr)
	! Get Nuisance Parameters for Projection effects
	status = status + datablock_get_double_default(block, number_count_section,"P_FRACT", 1.0D0, frac_p)
	status = status + datablock_get_double_default(block, number_count_section,"EPSILON", 0.0D0, epsi)
	! Get parameters for Gaussian priors for Nuisance parameters
	status = status + datablock_get_double_default(block, number_count_section,"MEAN_RM_SLOPE", 0.0D0, mean_delta_a0)
	status = status + datablock_get_double_default(block, number_count_section,"MEAN_RM_INTER", 0.0D0, mean_delta_b0)
	status = status + datablock_get_double_default(block, number_count_section,"MEAN_RM_Z_EVOL", 0.0D0, mean_delta_c0)
	status = status + datablock_get_double_default(block, number_count_section,"MEAN_RM_INTR_SIG", 0.0D0,mean_sig_intr)
	status = status + datablock_get_double_default(block, number_count_section,"MEAN_HMF_SLOPE", 0.0D0, mean_slope)
	status = status + datablock_get_double_default(block, number_count_section,"MEAN_HMF_INTER", 0.0D0, mean_inter)
	status = status + datablock_get_double_default(block, number_count_section,"MEAN_P_FRACT", 0.0D0,mean_p_frac)
	status = status + datablock_get_double_default(block, number_count_section,"MEAN_EPSILON", 0.0D0, mean_epsi)
	status = status + datablock_get_double_default(block, number_count_section,"SIG_RM_SLOPE", 0.0D0, sig_delta_a0)
	status = status + datablock_get_double_default(block, number_count_section,"SIG_RM_INTER", 0.0D0, sig_delta_b0)
	status = status + datablock_get_double_default(block, number_count_section,"SIG_RM_Z_EVOL", 0.0D0, sig_delta_c0)
	status = status + datablock_get_double_default(block, number_count_section,"SIG_RM_INTR_SIG", 0.0D0,sig_sig_intr)
	status = status + datablock_get_double_default(block, number_count_section,"SIG_HMF_SLOPE", 0.0D0, sig_slope)
	status = status + datablock_get_double_default(block, number_count_section,"SIG_HMF_INTER", 0.0D0, sig_inter)
	status = status + datablock_get_double_default(block, number_count_section,"SIG_P_FRACT", 0.0D0, sig_p_frac)
	status = status + datablock_get_double_default(block, number_count_section,"SIG_EPSILON", 0.0D0, sig_epsi)

	! Re-define parameters according to the usal definition
	Omega_m = Omega_dm + Omega_nu ! Omega_m = Omega_cdm + Omega_baryon + Omega_nu

    ! Blind the cosmological parameters
    if (settings%blind) then
        open(UNIT, file="/home/costanzi/blind_displacement.dat", form='unformatted',access='direct', recl=4)
        do i_blind = 0,SAMPLE_LENGTH-1
            read(UNIT, rec=i_blind+1) delta_blind(i_blind)
!           write(*,*) delta_blind(i_blind)
        end do
        Omega_dm = Omega_dm + delta_blind(0)
        Omega_k = Omega_k + delta_blind(2)
        Omega_nu = Omega_nu + delta_blind(6)
        Omega_m = Omega_dm + Omega_nu
        Omega_v = 1 - Omega_m - Omega_k
        w_0 = w_0 + delta_blind(7)
        w_a = w_a + delta_blind(8)
        A_s = sqrt((A_s)/(A_s + delta_blind(3))) ! used only to rescale NC%sigma8(A_s+dA_s,Om+dOm) -> NC%sigma8(A_s,Om+dOm)
!        write(*,*) delta_blind
!        write(*,*) Omega_dm,Omega_k,Omega_nu,Omega_m,Omega_v,w_0,w_a
    end if
    if (.not. settings%blind) A_s = 1.d0

	!load in the matter power spectrum
	status = load_matter_power(block,PK)
	if (status .ne. 0) then
		write(*,*) "Could not load matter power"
		status=3
		return 
	endif

	! SURVEY AREA IN STERADIAN
	Delta_Omega=settings%Delta_Omg
!	! READ EFFECTIVE AREA
	call read_effective_area()
	! SET THE NUMBER OF LAMBDA AND REDSHIFT BIN TO COMPUTE
	NC%num_of_L_bin= settings%num_lambda_bin
	NC%num_of_z_bin= settings%num_redshift_bin
	call allocate_NC(NC,settings)
    ! =========================================================
	Delta_c=settings%Delta_crit ! Set Delta_crit > 0 in the .ini file for halo masses defined wrt the critical density
	Delta_m=settings%Delta_mean ! Set Delta_mean > 0 in the .ini file for halo masses defined wrt the mean background density
	Delta_m_omp=Delta_m ! Variable set to avid problem with OMP

	total_lnlike=0.d0
	NC%cov_i_j=0.d0
	NC%data_i=0.d0

!	call compute_Cl(PK,cosmo,NC)
!	stop 24
	use_WL=settings%use_WL_mass ! Global variable used in compute_NC
	use_NC=settings%use_NC_data ! Global variable used in compute_NC
	if (settings%use_NC_data) call read_NC_data(NC,settings%file_data_NC) ! read the data
!	if (settings%use_WL_mass) call read_WL_data(NC,settings%file_data_WL) ! read the data
!	if (settings%use_WL_mass) call read_COV_WL_data(NC,settings%file_cov_WL) ! read COV of WL data
!	if (NC%z_max_array(NC%num_of_z_bin)>PK%redshifts(PK%num_z)) then
!	    write(*,*) "WARNING out of redshift range computed by CAMB; increse z_max" 
!	    return
!	end if
	call compute_NC(PK,NC) ! compute n(Lambda_i,z_j))

!	! Reascle WL data and Cov matrix for Omega_m ==============
!	if (settings%use_WL_mass) NC%mean_m_Li_zj_data=NC%mean_m_Li_zj_data*10.d0**(-0.39d0*(Omega_m-0.3d0))
!	if (settings%use_NC_data .and. settings%use_WL_mass) NC%cov_i_j(NC%num_of_L_bin*NC%num_of_z_bin+1:,NC%num_of_L_bin*NC%num_of_z_bin+1:)=NC%cov_i_j(NC%num_of_L_bin*NC%num_of_z_bin+1:,NC%num_of_L_bin*NC%num_of_z_bin+1:)*10.d0**(-0.39d0*(Omega_m-0.3d0)*2.d0)
!	if (.not. settings%use_NC_data .and. settings%use_WL_mass) NC%cov_i_j=NC%cov_i_j*10.d0**(-0.39d0*(Omega_m-0.3d0)*2.d0)
!!	write(*,*) log10( NC%mean_m_Li_zj_data),log10(abs(NC%mean_m_Li_zj_data-NC%mean_m_Li_zj))
!    ! =========================================================

!	! Sorting the covariance and data array conviniently to compute loglike
	do nz=1,NC%num_of_z_bin ! Start loop over redshift bin
	    do nL=1,NC%num_of_L_bin ! Start loop over Lambda bin
	        index_i=(nz-1)*NC%num_of_L_bin+nL
	        if (settings%use_NC_data) NC%data_i(index_i)= NC%n_Li_zj(nL,nz)-NC%n_Li_zj_data(nL,nz) ! NC data
!!	        write(*,'(8e16.8)') nz*1.d0,nL*1.d0,(8.d13*10.d0**((nL-1)*(0.2d0))+8.d13*10.d0**((nL)*(0.2d0)))*0.5d0,NC%data_i(index_i)/NC%n_Li_zj_data(nL,nz)
	        if (.not. settings%use_NC_data .and. settings%use_WL_mass) NC%data_i(index_i)= NC%mean_m_Li_zj(nL,nz)-NC%mean_m_Li_zj_data(nL,nz) ! WL data
	        if (settings%use_WL_mass .and. settings%use_NC_data) NC%data_i(index_i+NC%num_of_L_bin*NC%num_of_z_bin)= NC%mean_m_Li_zj(nL,nz)-NC%mean_m_Li_zj_data(nL,nz) ! WL data
!	        NC%data_i(index_i)= NC%mean_m_Li_zj(nL,nz)-NC%mean_m_Li_zj_data(nL,nz) ! WL data only
	        do nz2=1,NC%num_of_z_bin
	            do nL2=1,NC%num_of_L_bin
	                index_j=(nz2-1)*NC%num_of_L_bin+nL2
!	                if (nz2>=nz .and. nL2>=nL) NC%cov_i_j(index_i,index_j)= NC%cov_Li_zi_Lj_zj(nL,nz,nL2,nz2)
!	                if (nz2<=nz .and. nL2<nL)  NC%cov_i_j(index_i,index_j)= NC%cov_Li_zi_Lj_zj(nL2,nz2,nL,nz) ! matrice simmetrica
                    if (settings%use_NC_data) then
	                    if (index_j>=index_i) NC%cov_i_j(index_i,index_j)= NC%cov_Li_zi_Lj_zj(nL,nz,nL2,nz2)
	                    if (index_j<index_i)  NC%cov_i_j(index_i,index_j)= NC%cov_i_j(index_j,index_i) ! matrice simmetrica
	                end if
!	                NC%cov_i_j(index_i,index_j)= NC%cov_Li_zi_Lj_zj(nL,nz,nL2,nz2)
!	                NC%cov_i_j(index_j,index_i)= NC%cov_i_j(index_i,index_j) ! matrice simmetrica
!	                if (nz2==nz .and. nL2==nL) then ! start test if the diagonal matrix works
!	                    NC%cov_i_j(index_i,index_j)= NC%cov_Li_zi_Lj_zj(nL,nz,nL2,nz2)
!	                if (nz2>=nz .and. nL2>=nL .and. NC%cov_Li_zi_Lj_zj(nL,nz,nL2,nz2) .gt. 1.d-4) write(*,'(10e16.8)') nz*1.d0,nL*1.d0,nz2*1.d0,nL2*1.d0,NC%cov_Li_zi_Lj_zj(nL,nz,nL2,nz2)!,-(NC%data_i(index_i))**2.d0/NC%cov_i_j(index_i,index_j)/2.d0 - 0.5*log(2*pi_value*NC%cov_i_j(index_i,index_j))
!	                else 
!	                    NC%cov_i_j(index_i,index_j)=0.d0
!	                end if ! end test if the diagonal matrix works
                    ! Only WL data
!                    if (.not. settings%use_NC_data .and. settings%use_WL_mass) then
!                        if (nz2==nz .and. nL2==nL) NC%cov_i_j(index_i,index_j)=(NC%mean_m_Li_zj_data(nL,nz)*0.04d0)**2.d0 ! WL Cov Matrix
!                    end if
!                    ! NC and WL data
!                    if (settings%use_WL_mass .and. settings%use_NC_data) then
!                        if (nz2==nz .and. nL2==nL) NC%cov_i_j(index_i+NC%num_of_L_bin*NC%num_of_z_bin,index_j+NC%num_of_L_bin*NC%num_of_z_bin)=(NC%mean_m_Li_zj_data(nL,nz)*0.04d0)**2.d0 ! WL Cov Matrix
!                    end if
!                    if (nz2==nz .and. nL2==nL) NC%cov_i_j(index_i,index_j)=(NC%mean_m_Li_zj_data(nL,nz)*0.04d0)**2.d0 ! WL Cov Matrix only
	            end do
	        end do
	    end do
	end do
	
!	open(90,file='outputs_random/Cov_ij_4lbin_2zbin_sub_vol_kcut_SDSS_30_bestfit.txt')
!	do nz=1,NC%num_of_z_bin*NC%num_of_L_bin*2
!	    write(90,'(200e16.8)')NC%cov_i_j(nz,:)
!	end do
!	close(90)
	
	
	total_lnlike= - Matrix_GaussianLogLikeDouble(NC%cov_i_j,NC%data_i) ! non considera il termine -0.5ln(2*pi)
	
	! Gaussian prior on the nuisance parameter for the HMF ====
	! From the comparison of Buzzard HMF / Tinker HMF =========
	! computed for the fiducial values ========================

	total_lnlike=total_lnlike+LnGauss_prior(slope,mean_slope,sig_slope)+LnGauss_prior(inter,mean_inter,sig_inter)
	! =========================================================

!	! Gaussian prior on nuisance parameter for lnlambda(M,z) ==

	total_lnlike=total_lnlike+LnGauss_prior(delta_a0,mean_delta_a0,sig_delta_a0)+LnGauss_prior(delta_b0,mean_delta_b0,sig_delta_b0)+LnGauss_prior(sig_intr,mean_sig_intr,sig_sig_intr)
!	! =========================================================

!	! Gaussian prior on nuisance parameter for Projections ====

	total_lnlike=total_lnlike+LnGauss_prior(frac_p,mean_p_frac,sig_p_frac)+LnGauss_prior(epsi,mean_epsi,sig_epsi)
!	! =========================================================

!	! Add here other Gaussian prior if needed =================
!	total_lnlike=total_lnlike+LnGauss_prior(new_param,mean_new_param,sig_new_param)
!	! =========================================================

	status = datablock_put_double(block, likelihoods_section, "NC_LIKE", total_lnlike)
	status = status + datablock_put_double(block,number_count_section,"sigma8_NC",NC%sigma8_NC)

	if(settings%feedback >0) then
	    write(*,*) 'NC_LNLIKE: ',total_lnlike
	end if

	deallocate(z_eff_area)
	deallocate(eff_area)
	call deallocate_matterpower(PK)
	call deallocate_NC(NC)
!    STOP 1024
end function execute


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
