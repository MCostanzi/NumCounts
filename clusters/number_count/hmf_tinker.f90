	module hmf_tinker

	use quadpack ! if I want to use qags integration
	USE interface_tools_NC
	use utils_NC

	implicit none

	contains


	! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	! Get Fitting Paramter for Tinker HMF
	subroutine Get_fit_param(z_redshift)
	implicit none
	double precision z_redshift, a_scale_fac, om_mat_at_z
    ! INTERPOLATE FOR Delta_m=Delta_c/Omega_m(z) ===============
    if (Delta_c .gt. 1.d-3) then ! i.e. if Delta_c is defined
        a_scale_fac=1.d0/(1.d0+z_redshift) !scale factor
        om_mat_at_z=Omega_m*(1.d0+z_redshift)**3.d0/e(a_scale_fac) !mean matter density at scale factor a
        Delta_m=Delta_c/om_mat_at_z 
!        write(*,*) 'Delta_c: ',Delta_c ,'and Delta_m', Delta_m, 'at z=',z_redshift
    end if

    ! Cubic spline interpolation for the desired delta_m =======
    if (Delta_m == 200.d0) then ! to avoid interpolation problem at the edge
        aatin0=0.186d0
        atin0=1.47d0
        btin0=2.57d0
        ctin0 =1.19d0
    else
        aatin0=0.26d0 ! A0 value for Delta_m larger then 1600 fixed to 0.26
        if (Delta_m .lt. 1600.d0) aatin0=Cubic_spline_interp(Delta_m,aa_tin0,ddaa_tin0)
        atin0 =Cubic_spline_interp(Delta_m,a_tin0,dda_tin0)
        btin0 =Cubic_spline_interp(Delta_m,b_tin0,ddb_tin0)
        ctin0 =Cubic_spline_interp(Delta_m,c_tin0,ddc_tin0)
    end if

	end subroutine Get_fit_param

	! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	! Halo Mass Function dn/dM= f(sigma(M,z))/(4/3pi R^3) * d Ln(sigma(M,z))^-1/dM
	function dndM(halo_mass,z_redshift)
	implicit none
	double precision z_redshift,Radius,halo_mass,sigma2_Rz,dndM
	double precision s_fit,int_fit

	Radius=M_to_R(halo_mass)

	sigma2_Rz=sigma2R_at_z(Radius,z_redshift) ! from interpolation
	dndM = fnu_nu(sigma2_Rz,z_redshift)/(4.d0*pi_value*Radius**3.d0/3.d0)*dlog_sigdm(Radius,z_redshift,sigma2_Rz)/halo_mass
	dndM = dndM * (slope*(log10(halo_mass)-log10(64929581370445.531d0))+inter) ! Nuisance HMF pivot central point
	if (use_BuzzardHMFcorr) then ! If use Mock catalog from Buzzard you can use the actual correction to the Tinker HMF
	    if (z_redshift<0.1d0) then
	        s_fit=-0.0551919385627d0
	        int_fit=1.01860574955d0
	    else if (z_redshift>=0.1d0 .and. z_redshift<0.2d0) then
	        s_fit=-0.0201402191055d0
	        int_fit=1.0712385319d0
	    else if (z_redshift>=0.2d0 .and. z_redshift<0.3d0) then
	        s_fit=-0.059003735557d0
	        int_fit=1.01455431653d0
	    else if (z_redshift>=0.3d0 .and. z_redshift<0.4d0) then
	        s_fit=0.0179784529536d0
	        int_fit=1.04963055699d0
	    else if (z_redshift>=0.4d0 .and. z_redshift<0.5d0) then
	        s_fit=0.0640008689536d0
	        int_fit=1.03237512997d0
	    else
	        s_fit=0.0973768656564d0
	        int_fit=1.05784457948d0
	    end if
	    dndM = dndM * (s_fit*(log10(halo_mass)-log10(64929581370445.531d0))+int_fit)
	end if
	
	return
	end function dndM
	
	! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	! Tinker et al. 2008 mass function
	function fnu_nu(x,y) !x=sigma^2 y=redshift
	implicit double precision (a-h,o-z)
	implicit integer (i-n)
	
	aatin= aatin0*(1.d0+y)**(-0.14d0)
	atin  = atin0 * (1.d0 + y)**(-0.06d0)
	bbtin = 0.75d0/(log10(Delta_m/75.d0))
	bbbtin = -1.d0* bbtin**(1.2d0)
	b4tin = 10.d0**bbbtin
	btin = btin0 * (1.d0 + y)**(-b4tin)
	z=1.d0/sqrt(x)
	fnu_nu=aatin*((btin*z)**(atin)+1.d0)*exp(-ctin0*z**2.d0)

	return
	end function fnu_nu

	! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	!d ln(sig^-1)/ d M = - d sig^2/ d R * 1/ sig^2 * R / (6*M)
	function dlog_sigdm(x,z_redshift,sigma2_r_z) !x=rsigma (R relativo a sigma(R)),sigma2_r_z=sigma^2(R,z)
	implicit double precision (a-h,o-z)
	implicit integer (i-n)

	dlog_sigdm=dsigma2RdR_at_z(x,z_redshift)*x/sigma2_r_z/6.d0 ! 1/M lo metto direttamente quando computo la hmf
	! NB dsigma2RdR_at_z= -dsig2dr , obtained from interpolation
!	write(*,*) 'dlogsigmadm',dsig2dr(x,z_redshift),sigma2_r_z
	return
	end function dlog_sigdm
	
	! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	!d sig^2(R) / d R =1/(2*Pi*Pi) \int (dk k^3 P(k) * 2 * W(kR) * dW/dr)
	! SE INTEGRO CON QROMB SEMBREREBBE DARE LEGGERI PROBLEMI NUMERICI
	function dsig2dr(x,z_redshift) !x=rsigma (R relativo a sigma(R))
    implicit double precision (a-h,l-z)
	implicit integer (i-k)
	
	epsabs = 0.0E+00
    epsrel = 0.001E+00
    key = 6
	
	rsmooth=x ! Assign the Global Variable
	z_red=z_redshift ! Assign the Global Variable

	Ln_k_min=log(0.0001d0)
	Ln_k_max=log(200.d0)


    call qag(d_dsig2dr, Ln_k_min,Ln_k_max, epsabs, epsrel, key, ss, abserr, keval, ier )
	dsig2dr=ss
	return
	end function dsig2dr


	function d_dsig2dr(x) !integranda di d sigma^2 / d R; x=wave number
	implicit double precision (a-h,o-z)
	implicit integer (i-n)
	double precision f_k, R_k, dwdr 
	f_k=exp(x)
	R_k=f_k*rsmooth
	if (R_k < 5e-4) then
	    dwdr = 0.d0
	else
	    dwdr=3.d0 * (sin(R_k) * (R_k**2.d0 - 3.d0) + 3.d0 * R_k * cos(R_k)) / (R_k**4.d0) !derivata top_hat rispetto R
	endif
	d_dsig2dr=(f_k**4.d0)*Pk_at_z(f_k,z_red)*filtTH(R_k)*dwdr/pi_value/pi_value
	return
	end function d_dsig2dr

	! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	! sigma^2(R,z) = 1/2pi^2 int_0^inf dk k^2 P(k,z) W(kR)^2 = 1/2pi^2 int_-inf^inf dLn(k) k^3 P(k,z) W(kR)^2
!	! INTEGRATION WITH QAGS SEEMS TO BE THE MOST ACCURATE AND FAST
	function sigma_r(r,z_redshift) 
    implicit double precision (a-h,l-z)
	implicit integer (i-k)
	
	epsabs = 0.0E+00
    epsrel = 0.001E+00
    key = 6

    ! Extrema of the logarithmic integral
	Ln_k_min=log(0.0001d0)
	Ln_k_max=log(200.d0)
	! Assign Global Variables
	rsmooth=r ! smoothing radius
	z_red=z_redshift ! redshift considered
	

    call qag(dsigma_r, Ln_k_min,Ln_k_max, epsabs, epsrel, key, ss, abserr, keval, ier )
!    call qtrap_1e7(dsigma_r,Ln_k_min,Ln_k_max,sigma_r) ! da un po' di problemi numerici per alcuni valori di R
	sigma_r=ss/2.d0 /pi_value /pi_value
!	write(*,*) "sigmasq",sigma_r
	return
	end function sigma_r

	function dsigma_r(x) ! d sigma^2(R,z) = k^3 P(k,z) W(kR)^2 ; input: x= Ln(k)
	implicit double precision (a-h,l-z)
	implicit integer (i-k)

	wave_num=exp(x) ! wave_num = k
	rk=rsmooth*wave_num
	dsigma_r=wave_num**3.d0 * Pk_at_z(wave_num,z_red)*filtTH(rk)**2.d0
	return
	end function dsigma_r

!	! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	function filtTH(rk) ! Top-Hat filter k-space
	implicit none
	double precision rk,filtTH
    if (rk < 1.0e-6) then
        filtTH = 1.d0 !to avoid precision problems
	else 
	    filtTH=3.d0*(sin(rk)-rk*cos(rk))/rk**3.d0
	endif
	return
	end function filtTH

	! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	
	function M_to_R(halo_mass) ! conversion from mass to radius: input units [M_sun/h]
	implicit none
	double precision halo_mass, M_to_R
	
	M_to_R = (halo_mass*3.d0/(4.d0*pi_value*rho_c*Omega_dm))**onethr ! [Mpc/h] ; NB neutrino component neglected
	return
	end function M_to_R
	! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	function R_to_M(Radius) ! conversion from radius to mass: input units [Mpc/h]
	implicit none
	double precision  Radius, R_to_M
	
	R_to_M = 4.d0*pi_value*Radius**3.d0/3.d0*rho_c*Omega_dm ! [M_sun/h] ; NB neutrino component neglected
	return
	end function R_to_M
	! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	
	! :::::: E(z)^2 :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !E(z)^2=O_m * (1+z)^3 + (1- O_m - O_l) * (1+z)^2 + O_l * (1+z)^3 * exp [-3 int_1^a eosDE / a da]=
    !      = O_m / a^3 + (1- O_m - O_l)/a^2 + O_l/a^3 * exp [-3 int_1^a eosDE / a da]
	function e(a) ! input: scale factor
	implicit double precision (a-h,o-z)
	implicit integer (i-n)
!	e=a*omem+a**4*omev !for Omega_m + Omega_l = 1.0 & w_0=-1 & w_a=0
!	e=omem/(a**3)+(1-omem-omev)/a**2+omev !for w_0=-1 & w_a=0
	espo=-3.d0*(1.d0+w_0+w_a) ! per w_0 /= 1 & w_a /= 0
	e=Omega_m/(a**3.d0)+(Omega_k)/(a**2.d0)+Omega_v*a**espo*exp(3.d0*w_a*(a-1.d0)) ! per w_0 /= 1 & w_a /= 0

	return
	end function e
	!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

	! Bias Tinker et al. 2010 ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	function bias_tin(halo_mass,z_redshift) !x=delta_c/sigma
	implicit double precision (a-h,o-z)
	implicit integer (i-n)

	Radius=M_to_R(halo_mass)
	sigma2_Rz=sigma2R_at_z(Radius,z_redshift) ! from interpolation
	fy=delta_cr/sqrt(sigma2_Rz) ! delta_cr/sigma(R,z)

	y=log10(Delta_m) !log10(Delta)
	flarge_A = 1.d0+0.24d0*y*exp(-(4.d0/y)**4.d0)
	fsmall_A = 0.44d0*y-0.88d0
	flarge_B = 0.183d0
	fsmall_B = 1.5d0
	flarge_C = 0.019d0+0.107d0*y+0.19d0*exp(-(4.d0/y)**4.d0)
	fsmall_C = 2.4d0
	bias_tin=1.d0-flarge_A*fy**fsmall_A/(fy**fsmall_A+delta_cr**fsmall_A)+flarge_B*fy**fsmall_B+flarge_C*fy**fsmall_C
	return
	end function bias_tin

	! HMF BARYON CORRECTION FROM VELLISCIG14 1402.4461 ::::::::
	! (dn/dlogM)^DM+b / (dn/dlogM)^DM = 10^[A + B/ ::::::::::::
	! {1+exp[(-log(M_DM)+C)/D]}] ::::::::::::::::::::::::::::::
	function hmf_baryon_corr(halo_mass,z_redshift)
	double precision z_redshift,halo_mass
	double precision hmf_baryon_corr,p_b

	p_b=-p_a ! This follows from the fit in the Velliscig paper
	! this also ensure hmf_baryon_corr -> 0 for log10M>16.0

	hmf_baryon_corr=p_a + p_b / (1.d0 + dexp(-(log10(halo_mass)+p_c)/p_d))
	hmf_baryon_corr=10.d0**hmf_baryon_corr
	return
	end function hmf_baryon_corr

	end module hmf_tinker
