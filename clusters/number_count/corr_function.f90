	module corr_function
	use quadpack ! if I want to use qags integration
	USE interface_tools_NC
	use utils_NC
	implicit none
	
	! CORRELATION FUNCTION QUANTITIES =========================
	double precision, dimension(:), allocatable :: Ln_dist_of_xi !dove allocare i Ln(dist) a cui e' calcolato xi(R,z)
	double precision, dimension(:,:), allocatable :: out_xi !dove allocare i Ln(Xi) per diversi R e z
	double precision, dimension(:,:), allocatable :: ddout_xi !dove allocare i dd Ln(Xi) per diversi R e z; serve per xi_at_zbin
	integer :: n_log_R_step ! number of points at which Xi and dXidr are computed

	contains

    ! ALLOCATE ARRAYS TO STORE xi(R,z) and derivative ===
    ! this is do to speed up the code =========================
    subroutine allocate_xi
    implicit none
    double precision log_R_min,log_R_max,d_logr
    integer nr,nz
    
    n_log_R_step=100 ! number of points where compute xi(R,z) and dxi(R,z)/dR

    allocate(Ln_dist_of_xi(n_log_R_step))
    allocate(out_xi(n_log_R_step,num_of_z)) ! uso lo stesso numero di redshift del P(k)
    allocate(ddout_xi(n_log_R_step,num_of_z)) ! uso lo stesso numero di redshift del P(k)

    log_R_min=log(0.1d0)!(5.d12*3.d0/(4.d0*pi_value*rho_c*Omega_dm*Delta_m))**onethr) ! this must be equal or lower than the radius of the smallest halo
    log_R_max=log(100.d0) ! this must be equal or larger than dcom(z_max)-dcom(z_min)
    d_logr=(log_R_max-log_R_min)/(n_log_R_step-1.d0) 
    ! Compute and store xi(R,z) and dxi(R,z)/dR ===
    !$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(STATIC), PRIVATE(nr,nz)
    do nr=1,n_log_R_step
        Ln_dist_of_xi(nr)=log_R_min+d_logr*(nr-1)
        do nz=1,num_of_z
            out_xi(nr,nz)=log(calc_corr_at_R(exp(Ln_dist_of_xi(nr)),redshifts_pk(nz)))
!            out_xi(nr,nz)=log(corr_at_R_z(exp(Ln_dist_of_xi(nr)),redshifts_pk(nz)))
!            write(*,*) exp(Ln_dist_of_xi(nr)),exp(out_xi(nr,nz))
        end do
    end do ! end redshift loop
    !$OMP END PARALLEL DO
!    stop(24)
    do nz=1,num_of_z
        call spline(Ln_dist_of_xi,out_xi(:,nz),n_log_R_step,cllo,clhi,ddout_xi(:,nz))  ! dd Ln(Xi) for interpolation
    end do 
    end subroutine allocate_xi

	!==========================================================

	subroutine deallocate_xi

    deallocate(Ln_dist_of_xi)
    deallocate(out_xi) ! uso lo stesso numero di redshift del P(k)
    deallocate(ddout_xi) ! uso lo stesso numero di redshift del P(k)

    end subroutine deallocate_xi
	!==========================================================


	! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ! Get Xi(R,z) at particular R by interpolation
    ! NB the redshift of the Xi(R,z) is given by index_z
	function XiR_at_zbin(radius,index_z) result(XiR_interp)

    double precision radius
    integer index_z
    double precision logradius
    integer klo,khi
    double precision XiR_interp, dp
    double precision ho,a0,b0
    double precision, dimension(2) :: Xi_r, dd_Xi_r
    integer, save :: i_last = 1

    logradius = log(radius)
    if (logradius < Ln_dist_of_xi(1)) then
        Xi_r=out_xi(1:2,index_z)
        dp = (Xi_r(2)-Xi_r(1))/(Ln_dist_of_xi(2)-Ln_dist_of_xi(1))
        XiR_interp = Xi_r(1) + dp*(logradius-Ln_dist_of_xi(1))
    else if (logradius > Ln_dist_of_xi(n_log_R_step)) then
        !Do dodgy linear extrapolation on assumption accuracy of result won't matter
        Xi_r=out_xi(n_log_R_step-1:n_log_R_step,index_z)
        dp = (Xi_r(2)-Xi_r(1))/(Ln_dist_of_xi(n_log_R_step)-Ln_dist_of_xi(n_log_R_step-1))
        XiR_interp = Xi_r(2) + dp*(logradius-Ln_dist_of_xi(n_log_R_step))
    else
        ! cubic spline interpolation
        klo=min(i_last,n_log_R_step)
        do while (Ln_dist_of_xi(klo) > logradius)
            klo=klo-1
        end do
        do while (Ln_dist_of_xi(klo+1)< logradius)
            klo=klo+1
        end do
        i_last =klo
        khi=klo+1

        Xi_r=out_xi(klo:khi,index_z)
        dd_Xi_r = ddout_xi(klo:khi,index_z)

        ho=Ln_dist_of_xi(khi)-Ln_dist_of_xi(klo)
        a0=(Ln_dist_of_xi(khi)-logradius)/ho
        b0=1-a0

        XiR_interp = a0*Xi_r(1)+b0*Xi_r(2)+((a0**3-a0)*dd_Xi_r(1) &
        + (b0**3-b0)*dd_Xi_r(2))*ho**2/6
    end if

    XiR_interp = exp(XiR_interp)
    XiR_interp = max(1.d-50,XiR_interp)
    end function XiR_at_zbin
	! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !Get Xi(R,z) at particular R and z  by interpolation
    ! This function call XiR_at_zbin to interpolate for R
    function XiR_at_z(radius, z)


    double precision radius, z, XiR_at_z
    integer zlo, zhi, iz, itf
    double precision ho,a0,b0
    double precision, dimension(4) :: Xi_r, dd_Xi_r, zvec
    integer, save :: zi_last = 1

    if(z>redshifts_pk(num_of_z) .or. z<redshifts_pk(1)) then
        write (*,*) ' z out of bounds in XiR_at_z (',z,') max redshift is', redshifts_pk(num_of_z), 'min redshift is', redshifts_pk(1)
        STOP 24
    end if

    if (redshifts_pk(num_of_z) - z < 1.d-4) then ! if z = z_max
        XiR_at_z=XiR_at_zbin(radius,num_of_z)
    else
    zlo=min(zi_last,num_of_z)
    do while (redshifts_pk(zlo) > z)
        zlo=zlo-1
    end do
    do while (redshifts_pk(zlo+1)< z)
        zlo=zlo+1
    end do
    zi_last=zlo
    zhi=zlo+1

    if(zlo==1)then
        iz = 2
        zvec(2:4)=redshifts_pk(zlo:zhi+1)
        do itf=zlo, zhi+1
            Xi_r(iz) = log(XiR_at_zbin(radius,itf))
            iz=iz+1
        end do
        call spline_double(zvec(2:4),Xi_r(2:4),3,dd_Xi_r(2:4)) 
    else
        iz = 1
        zvec(:)=redshifts_pk(zlo-1:zhi+1)
        do itf=zlo-1, zhi+1
            Xi_r(iz) = log(XiR_at_zbin(radius,itf))
            iz=iz+1
        end do
        call spline_double(zvec,Xi_r,4,dd_Xi_r)
    end if

    ho=zvec(3)-zvec(2)
    a0=(zvec(3)-z)/ho
    b0=(z-zvec(2))/ho

    XiR_at_z = a0*Xi_r(2)+b0*Xi_r(3)+((a0**3-a0)*dd_Xi_r(2) &
    +(b0**3-b0)*dd_Xi_r(3))*ho**2/6

    XiR_at_z = exp(XiR_at_z)
    XiR_at_z = max(1.d-50,XiR_at_z)
    end if

    end function XiR_at_z
	! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



	!Correlation function ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	function corr_at_R_z(dist,z_redshift)
	implicit none
	double precision corr_at_R_z,z_redshift,dist,Ln_k_min,Ln_k_max
	double precision epsabs, epsrel, abserr
	integer neval, ier, key
	epsabs = 0.0E+00
	epsrel = 0.001E+00
    key = 6
    ! Extrema of the logarithmic integral
	Ln_k_min=log(0.0001d0)
	Ln_k_max=log(200.d0)

	z_global=z_redshift ! Assign global variable
	rsmooth=dist  ! Assign global variable
!	open(64,file='/home/costanzi/cosmosis_ltsp01/cosmosis_sdss_v2/outputs_random/pk_lin_xi_test2.dat')
	call qag(d_xi_R_z, Ln_k_min,Ln_k_max, epsabs, epsrel, key, corr_at_R_z, abserr, neval, ier )
	corr_at_R_z=corr_at_R_z/2.d0/pi_value/pi_value
	return
	end function corr_at_R_z
	
	function d_xi_R_z(lnk)
	implicit none
	double precision lnk,k_num,d_xi_R_z

	k_num=dexp(lnk)
    if (use_NLPk) then
	    d_xi_R_z= k_num**3.d0*NLPk_at_z(k_num,z_global)*dsin(k_num*rsmooth)/(k_num*rsmooth)
	else
	    d_xi_R_z= k_num**3.d0*Pk_at_z(k_num,z_global)*dsin(k_num*rsmooth)/(k_num*rsmooth)
	end if
!	write(64,*) k_num,Pk_at_z(k_num,z_global),dsin(k_num*rsmooth)/(k_num*rsmooth)
	return
	end function d_xi_R_z



!  FROM: https://github.com/tmcclintock/FastCorr/blob/master/fastcorr.c
!  This is the routine where xi(R) is actually calculated.
!  dist: tangential radial distance - units in either Mpc/h or Mpc
!  N: number of roots of j_0 to evaluate
!  h: step size for the quadrature routine
!  sinx: precalculated sin(x)
!  dpsi: precalculated function dpsi(x) (see Ogata et al. 2005)
 
	function calc_corr_at_R(dist,z_redshift,N_root,h_size)
	implicit none
	double precision f,sums,k_r,psi,PIsinht,dpsi,t,hsize
	double precision z_redshift,dist,calc_corr_at_R
	integer i,Nroot
	integer, optional :: N_root
	double precision, optional :: h_size

	if (present(N_root)) then
	    Nroot=N_root
	else
	    Nroot=600
	end if

	if (present(h_size)) then
	    hsize=h_size
	else
	    hsize=0.002d0
	end if

!open(64,file='/home/costanzi/cosmosis_ltsp01/cosmosis_sdss_v2/outputs_random/pk_lin_xi_test.dat')
	sums = 0.d0
	do i=0,Nroot
	    t=hsize*(i+1.d0)
	    psi = t*dtanh(dsinh(t)*0.5d0*pi_value)
	    k_r=psi*pi_value/hsize
	    PIsinht = pi_value*dsinh(t)
	    dpsi= (pi_value*t*dcosh(t)+dsinh(PIsinht))/(1.d0+dcosh(PIsinht));
	    if (use_NLPk) then
	        f = k_r*NLPk_at_z(k_r/dist,z_redshift)
	    else
	        f = k_r*Pk_at_z(k_r/dist,z_redshift)
	    end if
	    sums = sums + f*dsin(k_r)*dpsi
!	    write(64,*) k_r/dist,Pk_at_z(k_r/dist,z_redshift)
	end do
	calc_corr_at_R = sums/(dist**3.d0*pi_value)
	return
	end function calc_corr_at_R




	end module corr_function
