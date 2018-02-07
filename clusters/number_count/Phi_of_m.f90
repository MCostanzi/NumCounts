module Phi_of_m

USE interface_tools_NC
USE quadpack ! if I want to use qags integration
USE utils_NC
implicit none



contains

	! =========================================================
	! Probability of having a halo of mass M with =============
	! Lambda MOCK within a given range \Delta_Lambda ==========
	function Phi_M(halo_mass)
	double precision halo_mass, x_min, x_max, sigma2_lnL, Phi_M
	double precision sig2_lnl_z,sig2_intr,Mpivot


	if (use_Plobltr) then
	    ! For P(Lambda_ob|M,z) = \int dlambda_tr P(lambda_ob| ===== 
	    ! lambda_tr,z_tr) P(lambda_tr|M,z_tr)
	    call bilin_interp(z_int_P_lobMz,M_int_P_lobMz,int_P_lob_M_z(n_L1,:,:),redshift,log(halo_mass),Phi_M)
	else
!	    ! For P(Lambda|M) Lognormal : =============================
!	    ! Phi_M = 1/2 * {erfc[x(Lambda_min)]-erfc[x(Lambda_max)]} =
!	    ! w\ x =Ln(Lambda)-<Ln(Lambda)|(M,z)>/sqrt(2*sigma_lnL^2) =
!	    ! w/ sigma_lnL^2 = 1/exp(< ln lambda|M,z> ) + \sigma^2 ====
	    if (use_HOD_LM) then
	        Mpivot=10.d0**delta_b0
	    else ! use power law lmabda-mass relation
!	        Mpivot=10.d0**14.35217d0 ! Test DES BUZZARD
	        Mpivot=10.d0**14.344d0 ! SDSS DATA
	    end if
	    if (use_powM_scatter) then
	        sig2_intr=(sig_intr*(halo_mass/Mpivot)**delta_c0)**2.d0
	    else
	        sig2_intr=(sig_intr+ delta_c0*log(halo_mass/Mpivot))**2.d0
	    end if
	    if (only_intr_scat) then
	        sigma2_lnL =sig2_intr
	    else
	        sigma2_lnL =sig2_intr + (exp(LnLambda_M(halo_mass,redshift))-1.d0)*exp(-LnLambda_M(halo_mass,redshift))**2.d0
	        if (exp(LnLambda_M(halo_mass,redshift))< 1.d0) sigma2_lnL =sig2_intr
	    end if
	    if (use_skewnorm) then  ! If P(lambda_tr|M,z_tr)=SkewNormal
	        mass=halo_mass ! mass global variable OMP THREADPRIVATE used in P_l_tr_M_z
	        call trapzd_int_1args(P_l_tr_M_z,exp(LnLambda_min),exp(LnLambda_max),redshift,Phi_M,int(min(exp(LnLambda_max)-exp(LnLambda_min),50.d0)))
	    else ! If P(lambda_tr|M,z_tr)=LogNormal
	        if (LMrel_is_mean) then ! if <ln(Lambda(M))>=ln<Lambda(M)>-0.5*Sigma is the mean of the log-norm distriubtion
	            x_min = (LnLambda_min-LnLambda_M(halo_mass,redshift)+0.5d0*sigma2_lnL)/sqrt(2.d0*sigma2_lnL)
	            x_max = (LnLambda_max-LnLambda_M(halo_mass,redshift)+0.5d0*sigma2_lnL)/sqrt(2.d0*sigma2_lnL)
	        else ! if <ln(Lambda(M))> is the mean of the log-norm distriubtion
	            x_min = (LnLambda_min-LnLambda_M(halo_mass,redshift))/sqrt(2.d0*sigma2_lnL)
	            x_max = (LnLambda_max-LnLambda_M(halo_mass,redshift))/sqrt(2.d0*sigma2_lnL)
	        end if
	        Phi_M = 0.5d0*(derfc(x_min)-derfc(x_max))
	    end if
!	!	write(*,'(10e16.8)') redshift,log10(halo_mass),sigma2_lnL,erfc(x_min),erfc(x_max),Phi_M
	end if
	return
	end function Phi_M
	! =========================================================

	! =========================================================
	function dP_l_ob_M_z(lambda_tr)
	double precision lambda_tr,dP_l_ob_M_z,dummy
	double precision int_P_lob_ltr_ztr
	
	!Interpolate for \int_Delta_lambda_i d\lambda_obs P(\lambda_obs|\lambda_tr,z_tr)
!	call ESLGQ(z_int_P_lob,l_int_P_lob,int_P_lob_ltr(n_L1,:,:),line_counter_z,line_counter_l,redshift,lambda_tr,int_P_lob_ltr_ztr)
!	write(*,*) n_L1,redshift,lambda_tr,int_P_lob_ltr_ztr
	call bilin_interp(z_int_P_lob,l_int_P_lob,int_P_lob_ltr(n_L1,:,:),redshift,lambda_tr,int_P_lob_ltr_ztr)
	if (int_P_lob_ltr_ztr<0.d0) int_P_lob_ltr_ztr=0.d0
	if (int_P_lob_ltr_ztr>1.d0) int_P_lob_ltr_ztr=1.d0

!	call bilin_interp(z_int_P_lob,l_int_P_lob,int_P_lob_ltr(4,:,:),0.75d0,1.d0,int_P_lob_ltr_ztr)
!    write(*,*) int_P_lob_ltr_ztr
!	call bilin_interp(z_int_P_lob,l_int_P_lob,int_P_lob_ltr(4,:,:),0.75d0,100.d0,int_P_lob_ltr_ztr)
!    write(*,*) int_P_lob_ltr_ztr
!    stop 24

	dP_l_ob_M_z=P_l_tr_M_z(lambda_tr,redshift)*int_P_lob_ltr_ztr

	return
	end function dP_l_ob_M_z
	! =========================================================


!	! P(lambda_tr|M,z_tr)=Normal ==============================
!	function P_l_tr_M_z(lambda_tr,z_reds)
!	double precision lambda_tr,z_reds, sigma2_L,P_l_tr_M_z
!	double precision sig2_intr,Lambda_M

!    Lambda_M=exp(LnLambda_M(mass,z_reds)) ! <lambda_tr_NOPRJ|M,z> mass\redshift = global variables
!    sig2_intr=sig_intr**2.d0
!    sigma2_L =sig2_intr*(Lambda_M)**2.d0 + (Lambda_M-1.d0) ! Questo qllo giusto in generale
!    if (Lambda_M<= 0.0d0) sigma2_L =sig2_intr*0.001d0**2.d0 ! in teoria sarebbe zero
!    P_l_tr_M_z = exp(-0.5d0*(lambda_tr-Lambda_M)**2.d0/sigma2_L)/sqrt(2.d0*pi_value*sigma2_L)
!	return
!	end function P_l_tr_M_z
!	! =========================================================

!	! P(lambda_tr|M,z_tr)=Log-Normal ==========================
!	function P_l_tr_M_z(lambda_tr,z_reds)
!	double precision lambda_tr,z_reds, sigma2_lnL,P_l_tr_M_z
!	double precision sig2_intr,Lambda_M

!    Lambda_M=exp(LnLambda_M(mass,z_reds)) ! <lambda_tr_NOPRJ|M,z> mass\redshift = global variables
!    sig2_intr=sig_intr**2.d0
!    sigma2_lnL =sig2_intr + (Lambda_M-1.d0)/(Lambda_M)**2.d0
!    if (Lambda_M< 1.d0) sigma2_lnL =sig2_intr
!    P_l_tr_M_z = exp(-0.5d0*(log(lambda_tr)-(LnLambda_M(mass,z_reds)-0.5d0*sigma2_lnL))**2.d0/sigma2_lnL)/sqrt(2.d0*pi_value*sigma2_lnL)/lambda_tr
!	return
!	end function P_l_tr_M_z
!	! =========================================================

	! P(lambda_tr|M,z_tr)=Skew-Normal =========================
	function P_l_tr_M_z(lambda_tr,z_reds)
	double precision l_sat_M,sig_skew,skew,lambda_tr,z_reds, P_l_tr_M_z,sig2_ltr,sig_intr_M,Mpivot
	double precision term1,term2

	if (use_HOD_LM) then
	    Mpivot=10.d0**delta_b0
	else ! use power law lmabda-mass relation
!	    Mpivot=10.d0**14.35217d0 ! Test DES BUZZARD
	    Mpivot=10.d0**14.344d0 ! SDSS DATA
	end if
	if (use_powM_scatter) then
	    sig_intr_M=sig_intr*(mass/Mpivot)**delta_c0
	else
	    sig_intr_M=sig_intr +delta_c0*log(mass/Mpivot)
	end if
	
	if (use_skewnorm) then ! USE P(lambda_tr|M,z_tr)=Skew-Normal
	    l_sat_M=exp(LnLambda_M(mass,z_reds))-1.d0 ! mass = global variables
	    skew=interp2d(l_sat_M,sig_intr_M,n_lsat_grid,l_sat_grid,n_sig_grid,sig_intr_grid,skew_table,ddskew_table)
	    sig_skew=interp2d(l_sat_M,sig_intr_M,n_lsat_grid,l_sat_grid,n_sig_grid,sig_intr_grid,sig_skew_table,ddsig_skew_table)
    !	call bilin_interp(l_sat_grid,sig_intr_grid,sig_skew_table,l_sat_M,sig_intr,sig_skew)
    !	call bilin_interp(l_sat_grid,sig_intr_grid,skew_table,l_sat_M,sig_intr,skew)
	    sig2_ltr=sig_skew**2.d0

	    term1=dexp(-0.5*(lambda_tr-l_sat_M)**2./sig2_ltr)
	    if (skew<=0.d0) then ! the skew-norm is Gaussian
	        term2=1.d0
	    else
	        term2=derfc(-skew*(lambda_tr-l_sat_M)/sqrt(2.d0*sig2_ltr))
	    end if
        P_l_tr_M_z = term1*term2/sqrt(2.d0*pi_value*sig2_ltr)
    !	write(*,*) lambda_tr,l_sat_M,sig_intr_M,skew,sig_skew,P_l_tr_M_z
    !	stop 24
	else ! USE P(lambda_tr|M,z_tr)=Log-Normal
	    l_sat_M=exp(LnLambda_M(mass,z_reds)) ! is actually <lambda_true> not <lambda_sat>
	    if (only_intr_scat) then
	        sig2_ltr =sig_intr_M**2.d0
	    else
	        sig2_ltr =sig_intr_M**2.d0 + (l_sat_M-1.d0)/l_sat_M**2.d0
	        if (l_sat_M< 1.d0) sig2_ltr = sig_intr_M**2.d0
	    end if
	    if (LMrel_is_mean) then ! if <ln(Lambda(M))>=ln<Lambda(M)>-0.5*Sigma is the mean of the log-norm distriubtion
	        P_l_tr_M_z = exp(-0.5d0*(log(lambda_tr)-(LnLambda_M(mass,z_reds)-0.5d0*sig2_ltr))**2.d0/sig2_ltr)/sqrt(2.d0*pi_value*sig2_ltr)/lambda_tr
	    else ! if <ln(Lambda(M))> is the mean of the log-norm distriubtion
	        P_l_tr_M_z = exp(-0.5d0*(log(lambda_tr)-LnLambda_M(mass,z_reds))**2.d0/sig2_ltr)/sqrt(2.d0*pi_value*sig2_ltr)/lambda_tr
	    end if
	end if
	return
	end function P_l_tr_M_z
	! =========================================================

	! Observable-Mass relation from HOD model =================
	function LnLambda_M(halo_mass,z_reds)
	double precision halo_mass,Mpivot, LnLambda_M,z_reds


	if (use_HOD_LM) then
	    Mpivot=10.d0**delta_b0
	    if (halo_mass>10.d0**Log10Mminsat) then ! To avoid numerical error for halo_mass==10.d0**Log10Mminsat
	        LnLambda_M=log(1.d0 + ((halo_mass-10.d0**Log10Mminsat)/Mpivot)**delta_a0*((1.d0+z_reds)/(1.d0+0.2d0))**epsi)
	    else ! by construction this should not happen
	        LnLambda_M=log(1.d0)
	    end if
	else ! use power law lmabda-mass relation
!	    Mpivot=10.d0**14.35217d0 ! Test DES BUZZARD
	    Mpivot=10.d0**14.344d0 ! SDSS DATA
	    LnLambda_M=delta_b0 + log(halo_mass/Mpivot)*delta_a0 + log((1.d0+z_reds)/(1.d0+0.2d0))*epsi
	end if
	return
	end function LnLambda_M
	! =========================================================

end module Phi_of_m
