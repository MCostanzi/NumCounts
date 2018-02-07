

    Module compute_dVdzdOmega
    use hmf_tinker
    implicit none

    contains
	! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	! ======== Comoving Volume between zmin and zmax per UNIT SOLID ANGLE =====================
	function com_vol(zmin,zmax)
	implicit double precision (a-z)
	double precision com_vol
	integer keval, ier
	
	epsabs = 0.0E-00
    epsrel = 0.001E+00

	call qags (derivs_com_dis, zmin,zmax, epsabs, epsrel, com_vol, abserr, keval, ier )

	return
	end function com_vol

	! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	!====== COMOVING VOLUME ELEMENT PER UNIT SOLID ANGLE AND REDSHIFT =========================
! d V / d Omega dz= (c/H_0)^3 * 1/|Omega_k| [sinn(sqrt(|Omega_k|) * int_0^z 1/E(z) dz)]^2 /E(z)
	function derivs_com_dis(z)
	implicit double precision (a-z)
	double precision integrale
	integer keval, ier
	
	epsabs = 0.0E-00
    epsrel = 0.001E+00
!    omem=0.3d0
!    omev=0.7d0
!	derivs_com_dis=(omem*(1+z)**3+1-omem)**(-0.5) ! for w_0=-1 & w_a=0 & O_k=0
!	call qromb(f_inv_e,0.d0,z,integrale) ! sembrerebbe indifferente qromb o qtrap
!	call qtrap(f_inv_e,0.d0,z,integrale)
	call qags (f_inv_e, 0.d0,z, epsabs, epsrel, integrale, abserr, keval, ier )


	if (abs(Omega_k) < 1.d-6) then ! Omega_k=0
		derivs_com_dis=(c_light)**3.d0*f_inv_e(z)*integrale**2.d0 !(Mpc/h)^3
	else if (Omega_k <= 1.d-6) then ! Omega_0>1
		derivs_com_dis=(c_light)**3.d0*f_inv_e(z)*(sin(sqrt(-Omega_k)*integrale))**2.d0/(-Omega_k) !(Mpc/h)^3
	else ! Omega_0<1
		derivs_com_dis=(c_light)**3.d0*f_inv_e(z)*(sinh(sqrt(Omega_k)*integrale))**2.d0/(Omega_k) !(Mpc/h)^3
	end if

	return
	end function derivs_com_dis

	! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	! Wrapper to integrate E(z)^-1 using qags

	function f_inv_e(z) !funzione inversa di E(z) per usare qromb/qtrap etc etc
	double precision a_scalef, f_inv_e,z
	
	a_scalef=1.d0/(1.d0+z)
	f_inv_e=1.d0/sqrt(e(a_scalef))
	
	return
	end function f_inv_e

	! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	
	! TRANSVERSE COMOVING DISTANCE: D(z)= c/H_0* int_0^z dz'/E(z')
    function dcom(z)
    implicit double precision (a-z)
    double precision integrale
	integer keval, ier
	
	epsabs = 0.0E-00
    epsrel = 0.001E+00
    
!    call qromb(f_inv_e,0.d0,z,integrale) !int_0^z dz'/E(z') ! sembrerebbe indifferente qromb o qtrap
	call qags (f_inv_e, 0.d0,z, epsabs, epsrel, integrale, abserr, keval, ier ) !int_0^z dz'/E(z')
!    call qtrap(f_inv_e,0.d0,z,integrale) !int_0^z dz'/E(z')

	if (abs(Omega_k) < 1.d-6) then ! Omega_k=0
		dcom=c_light*integrale!Mpc/h
	else if (Omega_k <= 1.d-6) then ! Omega_0>1
		dcom=c_light*sin(sqrt(-Omega_k)*integrale)/sqrt(-Omega_k) !Mpc/h
	else ! Omega_0<1
		dcom=c_light*sinh(sqrt(Omega_k)*integrale)/sqrt(Omega_k) !Mpc/h
	end if

	return
	end function dcom
	! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

	! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	
!	! EFFECTIVE AREA in STERADIAN for SDSS OLD ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!    function eff_area_z(zs)
!    double precision eff_area_z,zs
!    double precision poly_coeff_vol(5)

!	if (zs<0.31d0) then
!	    poly_coeff_vol=[-3.83606104d1,   2.34041353d1,  -4.11057491d0,  -1.90910475d-02, 3.0877461d0]
!	else if (zs>=0.31d0 .and. zs<0.4d0) then
!	    poly_coeff_vol=[ 100.48392173d0,  -85.05121543d0  , -0.4508528d0   , 13.05021727d0  ,  0.6336882d0]
!	else
!	    poly_coeff_vol=[0.0d0, -33.10516966d0,  52.58469869d0, -28.45746784d0,   8.00190717d0]
!	end if
!	
!	eff_area_z= poly_coeff_vol(1)*zs**4.d0+poly_coeff_vol(2)*zs**3.d0+poly_coeff_vol(3)*zs**2.d0+poly_coeff_vol(4)*zs +poly_coeff_vol(5)

!	return
!	end function eff_area_z
!	! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	
	! EFFECTIVE AREA in STERADIAN for SDSS  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	! FIT FROM dr8_rand_zmask2_redmapper_v5.10_randcat_z0.05-0.60_lgt005_lam020-140_vl05_area.fit :::::::::
    function eff_area_z(zs)
    double precision eff_area_z,z_pivot
    double precision, intent(in) :: zs ! uso intent(in) perche non si modifichi nell'ouput'
    double precision poly_coeff_vol(12)
    integer i,poly_deg

	poly_coeff_vol=[ -1.14293122d05,   5.96846869d04,   9.24239180d03,  -2.23118813d03, -4.52580713d03,   1.18404878d03,   1.27951911d02,  -5.05716847d01, 1.01744577d00,  -3.11253383d-01,   5.48481084d-03,   3.12629987d00]
	
	poly_deg=12
	z_pivot=0.2d0
	eff_area_z=0.d0
	do i=1,poly_deg
	    eff_area_z=eff_area_z+ poly_coeff_vol(i)*(zs-z_pivot)**(poly_deg-i)
	end do

	return
	end function eff_area_z
	! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



!	! EFFECTIVE AREA in STERADIAN for DES Y1 MOCK :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!    function eff_area_z(zs)
!    double precision eff_area_z
!    double precision, intent(in) :: zs ! uso intent(in) perche non si modifichi nell'ouput'
!    double precision poly_coeff_vol(7)

!	if (zs<0.525d0) then
!	    poly_coeff_vol=[0.d0, 0.d0, 0.d0,0.d0,0.d0,0.d0, 0.33748012d0]
!	else if (zs>=0.525d0 .and. zs<0.695d0) then
!	    poly_coeff_vol=[ -6.88583364d4,   9.63306282d3,   5.22302755d3,   3.59004532d2, -2.35436263d1,  -3.88872811d0,   1.80278818d-1]
!	    zs=zs-0.635d0
!	else
!	    poly_coeff_vol=[0.d0, 0.d0, 0.d0,0.d0,0.d0,0.d0,   0.0001863d0]
!	end if
!	
!	eff_area_z= poly_coeff_vol(1)*zs**6.d0+poly_coeff_vol(2)*zs**5.d0+poly_coeff_vol(3)*zs**4.d0+poly_coeff_vol(4)*zs**3.d0+poly_coeff_vol(5)*zs**2.d0+poly_coeff_vol(6)*zs +poly_coeff_vol(7)

!	return
!	end function eff_area_z
!	! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


	end module compute_dVdzdOmega
