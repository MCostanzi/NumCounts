module utils_NC

	USE interface_tools_NC
	implicit none

	! POWER SPECTRUM QUANTITIES ===============================
	double precision, dimension(:), allocatable :: redshifts_pk !dove allocare i redshift a cui e' calcolato il power spectrum
	double precision, dimension(:), allocatable :: Ln_k_of_pk !dove allocare i Ln(k-wavenumber) a cui e' calcolato il power spectrum
	double precision, dimension(:,:), allocatable :: outpower !dove allocare i Ln(matter power spectrum) per diversi k e z
	double precision, dimension(:,:), allocatable :: ddoutpower !dove allocare i dd Ln(matter power spectrum) per diversi k e z; serve per Pk_at_zbin
	integer :: num_of_z ! number of redshift at which P(k) is computed
	integer :: num_of_k ! number of k at which P(k) is computed
	! =========================================================
	! ARRAYS FOR TINKER HMF PARAMETERS=========================
	integer, parameter :: num_row=9 ! numero line file
	double precision :: Log_delta_input(num_row),aa_tin0(num_row),a_tin0(num_row),b_tin0(num_row),c_tin0(num_row) ! parametri letti
	double precision :: ddaa_tin0(num_row),dda_tin0(num_row),ddb_tin0(num_row),ddc_tin0(num_row) ! parametri letti
	double precision :: Delta_c, Delta_m, Delta_m_omp ! Reference overdensity used to define the halo masses
	double precision :: aatin0,atin0,btin0,ctin0 ! parametri fittati
	!$OMP THREADPRIVATE(aatin0,atin0,btin0,ctin0,Delta_m)
	!==========================================================
	! SIGMA^2(R,z) and dSIGMA^2(R,z)/dM =======================
	double precision, dimension(:), allocatable :: Ln_R_array ! array for sigma and dsigmadr interpolation
	double precision, dimension(:,:), allocatable :: Ln_sigma2Rz
	double precision, dimension(:,:), allocatable :: dd_Ln_sigma2Rz ! second derivative for interpolation
	double precision, dimension(:,:), allocatable :: Ln_dsigma2Rz
	double precision, dimension(:,:), allocatable :: dd_Ln_dsigma2Rz ! second derivative for interpolation
	integer :: n_log_M_step ! number of points at which sigma and dsigmadr are computed
	!==========================================================
    ! Array to allocate the \int_Delta_lambda_i d\lambda_obs ==
    !P(\lambda_obs|\lambda_tr,z_tr) ===========================
    double precision, dimension(:,:,:), allocatable :: int_P_lob_ltr ! used for interpolation
    double precision, dimension(:), allocatable :: z_int_P_lob !redshifts in cui e' calcolato int_P_lob_ltr
    double precision, dimension(:), allocatable :: l_int_P_lob !lambdas in cui e' calcolato int_P_lob_ltr
    integer line_counter_z,line_counter_l ! Number of redshift/lambdas in the table, used to interpolate
    !==========================================================
    ! Array to allocate \int d\lambda_in P(\lambda_tr|M,z) ====
    ! I(\Delta_lambda_i,\lambda_tr,z) =========================
    double precision, dimension(:,:,:), allocatable :: int_P_lob_M_z ! used for interpolation
    double precision, dimension(:), allocatable :: z_int_P_lobMz !redshifts in cui e' calcolato
    double precision, dimension(:), allocatable :: M_int_P_lobMz !lambdas in cui e' calcolato
    integer nz_step,nlogM_step ! Number of redshift/lambdas in the table, used to interpolate
    !==========================================================
    !==========================================================
    ! Array to allocate SKEW-Gauss parameters =================
    ! I(\Delta_lambda_i,\lambda_tr,z) =========================
    double precision, dimension(:,:), allocatable :: skew_table ! used for interpolation
    double precision, dimension(:,:), allocatable :: sig_skew_table ! used for interpolation
    double precision, dimension(:,:), allocatable :: ddskew_table ! used for interpolation
    double precision, dimension(:,:), allocatable :: ddsig_skew_table ! used for interpolation
    double precision, dimension(:), allocatable :: sig_intr_grid !sig_intr in cui e' calcolato
    double precision, dimension(:), allocatable :: l_sat_grid !lambda_sat in cui e' calcolato
    integer n_lsat_grid, n_sig_grid
    !==========================================================
	double precision, parameter :: cllo=1.d30,clhi=1.d30 ! spline interpolation parameters

    contains

    ! ALLOCATE POWER SPECTRUM QUANTITIES ======================
    subroutine allocate_pow_spec(PK)
    type(pk_settings) :: PK
    integer nz

    num_of_z=PK%num_z
    num_of_k=PK%num_k
    
    allocate(redshifts_pk(num_of_z))
    allocate(Ln_k_of_pk(num_of_k))
    allocate(outpower(num_of_k,num_of_z))
    allocate(ddoutpower(num_of_k,num_of_z)) ! second derivative of Pk per interpolare con Pk_at_zbin
    
    Ln_k_of_pk=log(PK%kh) ! Log(kh) for interpolation
    redshifts_pk=PK%redshifts ! power spectrum redshifts
    outpower=log(PK%matpower) ! Ln(matter power spectrum) for interpolation
    do nz=1,PK%num_z
        call spline(Ln_k_of_pk,outpower(:,nz),num_of_k,cllo,clhi,ddoutpower(:,nz))  ! dd Ln(matter power spectrum) for interpolation
    end do 
    
!    open(96,file='cosmo1_output/ps_camb_z0.0.dat.MatteoC')
!    do nz=1,PK%num_k
!        write(96,*) PK%kh(nz),PK%matpower(nz,1)
!    end do
!    close(96)
!    STOP 1024
    ! =========================================================
    end subroutine allocate_pow_spec

	subroutine deallocate_pow_spec
	
	deallocate(redshifts_pk)
    deallocate(Ln_k_of_pk)
    deallocate(outpower)
    deallocate(ddoutpower)
    
    end subroutine deallocate_pow_spec
    
	subroutine allocate_Ln_sigma

    allocate(Ln_R_array(n_log_M_step))
    allocate(Ln_sigma2Rz(n_log_M_step,num_of_z))
    allocate(dd_Ln_sigma2Rz(n_log_M_step,num_of_z))
    allocate(Ln_dsigma2Rz(n_log_M_step,num_of_z))
    allocate(dd_Ln_dsigma2Rz(n_log_M_step,num_of_z))

    end subroutine allocate_Ln_sigma
    
	subroutine deallocate_Ln_sigma
	deallocate(Ln_R_array)
    deallocate(Ln_sigma2Rz)
    deallocate(dd_Ln_sigma2Rz)
    deallocate(Ln_dsigma2Rz)
    deallocate(dd_Ln_dsigma2Rz)
    end subroutine deallocate_Ln_sigma

	! ::::::::::::::: INTERPOLATION FUNCTIONS ::::::::::::::::::::::::::::::::::::::::::::::::::
	! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ! Get matter power spectrum at particular k/h by interpolation, from power_spec.f90 (Cosmomc Dic13)
    ! NB the redshift of the Pk is given by index_z
	function Pk_at_zbin(kh,index_z) result(pk_interp)

    double precision kh
    integer index_z
    double precision logk
    integer klo,khi
    double precision pk_interp, dp
    double precision ho,a0,b0
    double precision, dimension(2) :: matpower, ddmat
    integer, save :: i_last = 1

    logk = log(kh)
    if (logk < Ln_k_of_pk(1)) then
        matpower=outpower(1:2,index_z)
        dp = (matpower(2)-matpower(1))/(Ln_k_of_pk(2)-Ln_k_of_pk(1))
        pk_interp = matpower(1) + dp*(logk-Ln_k_of_pk(1))
    else if (logk > Ln_k_of_pk(num_of_k)) then
        !Do dodgy linear extrapolation on assumption accuracy of result won't matter
        matpower=outpower(num_of_k-1:num_of_k,index_z)
        dp = (matpower(2)-matpower(1))/(Ln_k_of_pk(num_of_k)-Ln_k_of_pk(num_of_k-1))
        pk_interp = matpower(2) + dp*(logk-Ln_k_of_pk(num_of_k))
    else
        ! cubic spline interpolation
        klo=min(i_last,num_of_k)
        do while (Ln_k_of_pk(klo) > logk)
            klo=klo-1
        end do
        do while (Ln_k_of_pk(klo+1)< logk)
            klo=klo+1
        end do
        i_last =klo
        khi=klo+1

        matpower=outpower(klo:khi,index_z)
        ddmat = ddoutpower(klo:khi,index_z)

        ho=Ln_k_of_pk(khi)-Ln_k_of_pk(klo)
        a0=(Ln_k_of_pk(khi)-logk)/ho
        b0=1-a0

        pk_interp = a0*matpower(1)+b0*matpower(2)+((a0**3-a0)*ddmat(1) &
        + (b0**3-b0)*ddmat(2))*ho**2/6
    end if

    pk_interp = exp(max(-50.d0,pk_interp)) 

    end function Pk_at_zbin
	! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !Get matter power spectrum at particular k/h and z  by interpolation
    ! This function call Pk_at_zbin to interpolate for k
    function Pk_at_z(kh, z)


    double precision kh, z, Pk_at_z
    integer zlo, zhi, iz, itf
    double precision ho,a0,b0
    double precision, dimension(4) :: matpower, ddmat, zvec
    integer, save :: zi_last = 1

    if(z>redshifts_pk(num_of_z) .or. z<redshifts_pk(1)) then
        write (*,*) ' z out of bounds in Pk_at_z (',z,') max redshift is', redshifts_pk(num_of_z), 'min redshift is', redshifts_pk(1)
        STOP 24
    end if

    if (redshifts_pk(num_of_z) - z < 1.d-4) then ! if z = z_max
        Pk_at_z=Pk_at_zbin(kh,num_of_z)
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
            matpower(iz) = log(Pk_at_zbin(kh,itf))
            iz=iz+1
        end do
        call spline_double(zvec(2:4),matpower(2:4),3,ddmat(2:4))
    elseif(zhi==num_of_z)then
        iz = 1
        zvec(1:3)=redshifts_pk(zlo-1:zhi)
        do itf=zlo-1, zhi
            matpower(iz) = log(Pk_at_zbin(kh,itf))
            iz=iz+1
        end do
        call spline_double(zvec(1:3),matpower(1:3),3,ddmat(1:3))
    else
        iz = 1
        zvec(:)=redshifts_pk(zlo-1:zhi+1)
        do itf=zlo-1, zhi+1
            matpower(iz) = log(Pk_at_zbin(kh,itf))
            iz=iz+1
        end do
        call spline_double(zvec,matpower,4,ddmat)
    end if

    ho=zvec(3)-zvec(2)
    a0=(zvec(3)-z)/ho
    b0=(z-zvec(2))/ho

    Pk_at_z = a0*matpower(2)+b0*matpower(3)+((a0**3-a0)*ddmat(2) &
    +(b0**3-b0)*ddmat(3))*ho**2/6

    Pk_at_z = exp(max(-50.d0,Pk_at_z))
    end if

    end function Pk_at_z

	! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ! Get sigma^2(R,z) at particular R by interpolation
    ! NB the redshift of the sigma^2(R,z) is given by index_z
	function sigma2R_at_zbin(radius,index_z) result(sigma2R_interp)

    double precision radius
    integer index_z
    double precision logradius
    integer klo,khi
    double precision sigma2R_interp, dp
    double precision ho,a0,b0
    double precision, dimension(2) :: sigmasq_r, dd_sigmasq_r
    integer, save :: i_last = 1

    logradius = log(radius)
    if (logradius < Ln_R_array(1)) then
        sigmasq_r=Ln_sigma2Rz(1:2,index_z)
        dp = (sigmasq_r(2)-sigmasq_r(1))/(Ln_R_array(2)-Ln_R_array(1))
        sigma2R_interp = sigmasq_r(1) + dp*(logradius-Ln_R_array(1))
    else if (logradius > Ln_R_array(n_log_M_step)) then
        !Do dodgy linear extrapolation on assumption accuracy of result won't matter
        sigmasq_r=Ln_sigma2Rz(n_log_M_step-1:n_log_M_step,index_z)
        dp = (sigmasq_r(2)-sigmasq_r(1))/(Ln_R_array(n_log_M_step)-Ln_R_array(n_log_M_step-1))
        sigma2R_interp = sigmasq_r(2) + dp*(logradius-Ln_R_array(n_log_M_step))
    else
        ! cubic spline interpolation
        klo=min(i_last,n_log_M_step)
        do while (Ln_R_array(klo) > logradius)
            klo=klo-1
        end do
        do while (Ln_R_array(klo+1)< logradius)
            klo=klo+1
        end do
        i_last =klo
        khi=klo+1

        sigmasq_r=Ln_sigma2Rz(klo:khi,index_z)
        dd_sigmasq_r = dd_Ln_sigma2Rz(klo:khi,index_z)

        ho=Ln_R_array(khi)-Ln_R_array(klo)
        a0=(Ln_R_array(khi)-logradius)/ho
        b0=1-a0

        sigma2R_interp = a0*sigmasq_r(1)+b0*sigmasq_r(2)+((a0**3-a0)*dd_sigmasq_r(1) &
        + (b0**3-b0)*dd_sigmasq_r(2))*ho**2/6
    end if

    sigma2R_interp = exp(sigma2R_interp)
    sigma2R_interp = max(1.d-50,sigma2R_interp)
    end function sigma2R_at_zbin
	! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !Get sigma^2(R,z) at particular R and z  by interpolation
    ! This function call sigma2R_at_zbin to interpolate for R
    function sigma2R_at_z(radius, z)


    double precision radius, z, sigma2R_at_z
    integer zlo, zhi, iz, itf
    double precision ho,a0,b0
    double precision, dimension(4) :: sigmasq_r, dd_sigmasq_r, zvec
    integer, save :: zi_last = 1

    if(z>redshifts_pk(num_of_z) .or. z<redshifts_pk(1)) then
        write (*,*) ' z out of bounds in sigma2R_at_z (',z,') max redshift is', redshifts_pk(num_of_z), 'min redshift is', redshifts_pk(1)
        STOP 24
    end if

    if (redshifts_pk(num_of_z) - z < 1.d-4) then ! if z = z_max
        sigma2R_at_z=sigma2R_at_zbin(radius,num_of_z)
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
            sigmasq_r(iz) = log(sigma2R_at_zbin(radius,itf))
            iz=iz+1
        end do
        call spline_double(zvec(2:4),sigmasq_r(2:4),3,dd_sigmasq_r(2:4))
    elseif(zhi==num_of_z)then
        iz = 1
        zvec(1:3)=redshifts_pk(zlo-1:zhi)
        do itf=zlo-1, zhi
            sigmasq_r(iz) = log(sigma2R_at_zbin(radius,itf))
            iz=iz+1
        end do
        call spline_double(zvec(1:3),sigmasq_r(1:3),3,dd_sigmasq_r(1:3))
    else
        iz = 1
        zvec(:)=redshifts_pk(zlo-1:zhi+1)
        do itf=zlo-1, zhi+1
            sigmasq_r(iz) = log(sigma2R_at_zbin(radius,itf))
            iz=iz+1
        end do
        call spline_double(zvec,sigmasq_r,4,dd_sigmasq_r)
    end if

    ho=zvec(3)-zvec(2)
    a0=(zvec(3)-z)/ho
    b0=(z-zvec(2))/ho

    sigma2R_at_z = a0*sigmasq_r(2)+b0*sigmasq_r(3)+((a0**3-a0)*dd_sigmasq_r(2) &
    +(b0**3-b0)*dd_sigmasq_r(3))*ho**2/6

    sigma2R_at_z = exp(sigma2R_at_z)
    sigma2R_at_z = max(1.d-50,sigma2R_at_z)
    end if

    end function sigma2R_at_z
	! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ! Get d sigma^2(R,z)/d R at particular R by interpolation
    ! NB the redshift of the dsigma^2(R,z)/dR is given by index_z
	function dsigma2RdR_at_zbin(radius,index_z) result(dsigma2RdR_interp)

    double precision radius
    integer index_z
    double precision logradius
    integer klo,khi
    double precision dsigma2RdR_interp, dp
    double precision ho,a0,b0
    double precision, dimension(2) :: sigmasq_r, dd_sigmasq_r
    integer, save :: i_last = 1

    logradius = log(radius)
    if (logradius < Ln_R_array(1)) then
        sigmasq_r=Ln_dsigma2Rz(1:2,index_z)
        dp = (sigmasq_r(2)-sigmasq_r(1))/(Ln_R_array(2)-Ln_R_array(1))
        dsigma2RdR_interp = sigmasq_r(1) + dp*(logradius-Ln_R_array(1))
    else if (logradius > Ln_R_array(n_log_M_step)) then
        !Do dodgy linear extrapolation on assumption accuracy of result won't matter
        sigmasq_r=Ln_dsigma2Rz(n_log_M_step-1:n_log_M_step,index_z)
        dp = (sigmasq_r(2)-sigmasq_r(1))/(Ln_R_array(n_log_M_step)-Ln_R_array(n_log_M_step-1))
        dsigma2RdR_interp = sigmasq_r(2) + dp*(logradius-Ln_R_array(n_log_M_step))
    else
        ! cubic spline interpolation
        klo=min(i_last,n_log_M_step)
        do while (Ln_R_array(klo) > logradius)
            klo=klo-1
        end do
        do while (Ln_R_array(klo+1)< logradius)
            klo=klo+1
        end do
        i_last =klo
        khi=klo+1

        sigmasq_r=Ln_dsigma2Rz(klo:khi,index_z)
        dd_sigmasq_r = dd_Ln_dsigma2Rz(klo:khi,index_z)

        ho=Ln_R_array(khi)-Ln_R_array(klo)
        a0=(Ln_R_array(khi)-logradius)/ho
        b0=1-a0

        dsigma2RdR_interp = a0*sigmasq_r(1)+b0*sigmasq_r(2)+((a0**3-a0)*dd_sigmasq_r(1) &
        + (b0**3-b0)*dd_sigmasq_r(2))*ho**2/6
    end if

    dsigma2RdR_interp = exp(dsigma2RdR_interp)
    dsigma2RdR_interp = max(1.d-50,dsigma2RdR_interp)
    end function dsigma2RdR_at_zbin
	! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !Get dsigma^2(R,z)/dR at particular R and z  by interpolation
    ! This function call dsigma2RdR_at_zbin to interpolate for R
    function dsigma2RdR_at_z(radius, z)


    double precision radius, z, dsigma2RdR_at_z
    integer zlo, zhi, iz, itf
    double precision ho,a0,b0
    double precision, dimension(4) :: sigmasq_r, dd_sigmasq_r, zvec
    integer, save :: zi_last = 1

    if(z>redshifts_pk(num_of_z) .or. z<redshifts_pk(1)) then
        write (*,*) ' z out of bounds in dsigma2RdR_at_z (',z,') max redshift is', redshifts_pk(num_of_z), 'min redshift is', redshifts_pk(1)
        STOP 24
    end if

    if (redshifts_pk(num_of_z) - z < 1.d-4) then ! if z = z_max
        dsigma2RdR_at_z=dsigma2RdR_at_zbin(radius,num_of_z)
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
            sigmasq_r(iz) = log(dsigma2RdR_at_zbin(radius,itf))
            iz=iz+1
        end do
        call spline_double(zvec(2:4),sigmasq_r(2:4),3,dd_sigmasq_r(2:4))
    elseif(zhi==num_of_z)then
        iz = 1
        zvec(1:3)=redshifts_pk(zlo-1:zhi)
        do itf=zlo-1, zhi
            sigmasq_r(iz) = log(dsigma2RdR_at_zbin(radius,itf))
            iz=iz+1
        end do
        call spline_double(zvec(1:3),sigmasq_r(1:3),3,dd_sigmasq_r(1:3))
    else
        iz = 1
        zvec(:)=redshifts_pk(zlo-1:zhi+1)
        do itf=zlo-1, zhi+1
            sigmasq_r(iz) = log(dsigma2RdR_at_zbin(radius,itf))
            iz=iz+1
        end do
        call spline_double(zvec,sigmasq_r,4,dd_sigmasq_r)
    end if

    ho=zvec(3)-zvec(2)
    a0=(zvec(3)-z)/ho
    b0=(z-zvec(2))/ho

    dsigma2RdR_at_z = a0*sigmasq_r(2)+b0*sigmasq_r(3)+((a0**3-a0)*dd_sigmasq_r(2) &
    +(b0**3-b0)*dd_sigmasq_r(3))*ho**2/6

    dsigma2RdR_at_z = exp(dsigma2RdR_at_z)
    dsigma2RdR_at_z = max(1.d-50,dsigma2RdR_at_z)
    end if

    end function dsigma2RdR_at_z

	! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	! Cubic spline interpolation for Tinker HMF paramter
	function Cubic_spline_interp(Delta,tin_param_array,ddtin_param_array)
    double precision Delta,Log_Delta
    integer llo,lhi
    double precision Cubic_spline_interp
    double precision ho,a0,b0
    double precision, dimension(9) :: tin_param_array, ddtin_param_array
    integer, save :: i_last = 1

    Log_Delta = log10(Delta)

    if (Log_Delta < Log_delta_input(1) .or. Log_Delta > Log_delta_input(num_row)) then
        write (*,*) 'Cubic_spline_interp: value Delta out of range ', Delta
        stop
    end if

    llo=1
    do while (Log_delta_input(llo+1) < Log_Delta)  
        llo = llo + 1
    end do

    lhi=llo+1
    ho=Log_delta_input(lhi)-Log_delta_input(llo)
    a0=(Log_delta_input(lhi)-Log_Delta)/ho
    b0=(Log_Delta-Log_delta_input(llo))/ho
    
    Cubic_spline_interp = a0*tin_param_array(llo)+ b0*tin_param_array(lhi)+((a0**3-a0)* ddtin_param_array(llo) &
    +(b0**3-b0)*ddtin_param_array(lhi))*ho**2/6.

    return
    end function Cubic_spline_interp
	! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! SUBROUTINE FOR INTERPOLATION ============================
	SUBROUTINE ENLGR(X,Y,N,T,Z)
	implicit double precision (a-h,o-z)
	implicit integer (i-n)
	dimension X(N),Y(N)
	DOUBLE PRECISION X,Y,T,Z,S

	Z=0.0
	IF (N.LE.0) RETURN
	IF (N.EQ.1) THEN
	  Z=Y(1)
	  RETURN
	END IF
	IF (N.EQ.2) THEN
	  Z=(Y(1)*(T-X(2))-Y(2)*(T-X(1)))/(X(1)-X(2))
	  RETURN
	END IF
	I=1
10	IF (X(I).LT.T) THEN
	  I=I+1
	  IF (I.LE.N) GOTO 10
	END IF
	K=I-4
	IF (K.LT.1) K=1
	M=I+3
	IF (M.GT.N) M=N
	DO 30 I=K,M
	  S=1.0
	  DO 20 J=K,M
	    IF (J.NE.I) THEN
	      S=S*(T-X(J))/(X(I)-X(J))
	    END IF
20	  CONTINUE
	  Z=Z+S*Y(I)
30	CONTINUE
	RETURN
	END SUBROUTINE ENLGR

    ! SUBROUTINE FOR 2D INTERPOLATION =========================
	SUBROUTINE ESLGQ(X,Y,Z,N,M,U,V,W)
	implicit integer (i-n)
	DIMENSION X(N),Y(M),Z(N,M),B(10)
	DOUBLE PRECISION X,Y,Z,U,V,W,B,HH
	IF (U.LE.X(1)) THEN
	  IP=1
	  IPP=4
	ELSE IF (U.GE.X(N)) THEN
	  IP=N-3
	  IPP=N
	ELSE
	  I=1
	  J=N
10	  IF (IABS(I-J).NE.1) THEN
	    L=(I+J)/2
	    IF (U.LT.X(L)) THEN
	      J=L
	    ELSE
	      I=L
	    END IF
	    GOTO 10
	  END IF
	  IP=I-3
	  IPP=I+4
	END IF
	IF (IP.LT.1) IP=1
	IF (IPP.GT.N) IPP=N
	IF (V.LE.Y(1)) THEN
	  IQ=1
	  IQQ=4
	ELSE IF (V.GE.Y(M)) THEN
	  IQ=M-3
	  IQQ=M
	ELSE
	  I=1
	  J=M
20	  IF (IABS(J-I).NE.1) THEN
	    L=(I+J)/2
	    IF (V.LT.Y(L)) THEN
	      J=L
	    ELSE
	      I=L
	    END IF
	    GOTO 20
	  END IF
	  IQ=I-3
	  IQQ=I+4
	END IF
	IF (IQ.LT.1) IQ=1
	IF (IQQ.GT.M) IQQ=M
	DO 50 I=IP,IPP
	  B(I-IP+1)=0.0
	  DO 40 J=IQ,IQQ
	    HH=Z(I,J)
	    DO 30 K=IQ,IQQ
	      IF (K.NE.J) THEN
	        HH=HH*(V-Y(K))/(Y(J)-Y(K))
	      END IF
30	    CONTINUE
	    B(I-IP+1)=B(I-IP+1)+HH
40	  CONTINUE
50	CONTINUE
	W=0.0
	DO 70 I=IP,IPP
	  HH=B(I-IP+1)
	  DO 60 J=IP,IPP
	    IF (J.NE.I) THEN
	      HH=HH*(U-X(J))/(X(I)-X(J))
	    END IF
60	  CONTINUE
	  W=W+HH
70	CONTINUE
	RETURN
	END SUBROUTINE ESLGQ

    ! SUBROUTINE FOR SPILNE INTERPOLATION
	  subroutine spline_double(x,y,n,d2)
      integer, intent(in) :: n
      integer, parameter :: dp=KIND(1.d0)
      real(dp), intent(in) :: x(n), y(n)
      real(dp), intent(out) :: d2(n)
      real(dp), dimension(:), allocatable :: u
      integer i
      real(dp) xp,sig,xxdiv,d1l,d1r

      allocate(u(1:n-1))

      d2(1)=0._dp
      u(1)=0._dp

      d1r= (y(2)-y(1))/(x(2)-x(1))
      do i=2,n-1
        d1l=d1r
        d1r=(y(i+1)-y(i))/(x(i+1)-x(i))
        xxdiv=1._dp/(x(i+1)-x(i-1))
        sig=(x(i)-x(i-1))*xxdiv
        xp=1._dp/(sig*d2(i-1)+2._dp)
        d2(i)=(sig-1._dp)*xp
        u(i)=(6._dp*(d1r-d1l)*xxdiv-sig*u(i-1))*xp
      end do

      d2(n)=0._dp
      do i=n-1,1,-1
        d2(i)=d2(i)*d2(i+1)+u(i)
      end do

      deallocate(u)
    end subroutine spline_double
	
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! calculates array of second derivatives used by cubic spline
! interpolation. y2 is array of second derivatives, yp1 and ypn are first
! derivatives at end points.

!Thanks Martin Reinecke
    subroutine spline(x,y,n,d11,d1n,d2)
      integer, parameter :: dl = KIND(1.d0)
      integer, parameter :: sp = KIND(1.0)
      integer, intent(in) :: n
      real(dl), intent(in) :: x(n), y(n), d11, d1n
      real(dl), intent(out) :: d2(n)
      integer i
      real(dl) xp,qn,sig,un,xxdiv,u(n-1),d1l,d1r

      d1r= (y(2)-y(1))/(x(2)-x(1))
      if (d11>.99e30_dl) then
        d2(1)=0._dl
        u(1)=0._dl
      else
        d2(1)=-0.5_dl
        u(1)=(3._dl/(x(2)-x(1)))*(d1r-d11)
      endif

      do i=2,n-1
        d1l=d1r
        d1r=(y(i+1)-y(i))/(x(i+1)-x(i))
        xxdiv=1._dl/(x(i+1)-x(i-1))
        sig=(x(i)-x(i-1))*xxdiv
        xp=1._dl/(sig*d2(i-1)+2._dl)

        d2(i)=(sig-1._dl)*xp

        u(i)=(6._dl*(d1r-d1l)*xxdiv-sig*u(i-1))*xp
      end do
      d1l=d1r

      if (d1n>.99e30_dl) then
        qn=0._dl
        un=0._dl
      else
        qn=0.5_dl
        un=(3._dl/(x(n)-x(n-1)))*(d1n-d1l)
      endif

      d2(n)=(un-qn*u(n-1))/(qn*d2(n-1)+1._dl)
      do i=n-1,1,-1
        d2(i)=d2(i)*d2(i+1)+u(i)
      end do
    end subroutine spline 

     SUBROUTINE spline_deriv(x,y,y2,y1,n)
     !Get derivative y1 given array of x, y and y''
      implicit none
      integer, parameter :: dl = KIND(1.d0)
      integer, parameter :: sp = KIND(1.0)
      INTEGER, intent(in) :: n
      real(dl), intent(in) :: x(n), y(n), y2(n)
      real(dl), intent(out) :: y1(n)
      INTEGER i
      real(dl) dx

      do i=1, n-1
           
         dx = (x(i+1) - x(i))
         y1(i) = (y(i+1) - y(i))/dx - dx*(2*y2(i) + y2(i+1))/6
      end do
       dx = x(n) - x(n-1)
       y1(n) = (y(n) - y(n-1))/dx + dx* ( y2(i-1)  + 2*y2(i) )/6

      END SUBROUTINE spline_deriv
      
      
         SUBROUTINE BJL(L,X,JL)
        !!== MODIFIED SUBROUTINE FOR SPHERICAL BESSEL FUNCTIONS.                       ==!!
        !!== CORRECTED THE SMALL BUGS IN PACKAGE CMBFAST&CAMB(for l=4,5, x~0.001-0.002)==!! 
        !!== CORRECTED THE SIGN OF J_L(X) FOR X<0 CASE                                 ==!!
        !!== WORKS FASTER AND MORE ACCURATE FOR LOW L, X<<L, AND L<<X cases            ==!! 
        !!== zqhuang@astro.utoronto.ca                                                 ==!!
        IMPLICIT NONE
        INTEGER L
        real(dl) X,JL
        real(dl) AX,AX2
        real(dl),PARAMETER::LN2=0.6931471805599453094D0
        real(dl),PARAMETER::ONEMLN2=0.30685281944005469058277D0
        real(dl),PARAMETER::PID2=1.5707963267948966192313217D0
        real(dl),PARAMETER::PID4=0.78539816339744830961566084582D0
        real(dl),parameter::ROOTPI12 = 21.269446210866192327578D0
        real(dl),parameter::GAMMA1 =   2.6789385347077476336556D0 !/* Gamma function of 1/3 */
        real(dl),parameter::GAMMA2 =   1.3541179394264004169452D0 !/* Gamma function of 2/3 */
        real(dl),PARAMETER::PI=3.141592653589793238463D0
        real(dl) NU,NU2,BETA,BETA2,COSB
        real(dl) sx,sx2
        real(dl) cotb,cot3b,cot6b,secb,sec2b
        real(dl) trigarg,expterm,L3

        IF(L.LT.0)THEN
            write(*,*) 'Can not evaluate Spherical Bessel Function with index l<0'
            STOP
        ENDIF
        AX=DABS(X)
        AX2=AX**2
        IF(L.LT.7)THEN
            IF(L.EQ.0)THEN
                IF(AX.LT.1.D-1)THEN
                    JL=1.D0-AX2/6.D0*(1.D0-AX2/20.D0)
                ELSE
                    JL=DSIN(AX)/AX
                ENDIF

            ELSEIF(L.EQ.1)THEN
                IF(AX.LT.2.D-1)THEN
                    JL=AX/3.D0*(1.D0-AX2/10.D0*(1.D0-AX2/28.D0))
                ELSE
                    JL=(DSIN(AX)/AX-DCOS(AX))/AX
                ENDIF
            ELSEIF(L.EQ.2)THEN
                IF(AX.LT.3.D-1)THEN
                    JL=AX2/15.D0*(1.D0-AX2/14.D0*(1.D0-AX2/36.D0))
                ELSE
                    JL=(-3.0D0*DCOS(AX)/AX-DSIN(AX)*(1.D0-3.D0/AX2))/AX
                ENDIF
            ELSEIF(L.EQ.3)THEN
                IF(AX.LT.4.D-1)THEN
                    JL=AX*AX2/105.D0*(1.D0-AX2/18.D0*(1.D0-AX2/44.D0))
                ELSE
                    JL=(DCOS(AX)*(1.D0-15.D0/AX2)-DSIN(AX)*(6.D0-15.D0/AX2)/AX)/AX
                ENDIF
            ELSEIF(L.EQ.4)THEN
                IF(AX.LT.6.D-1)THEN
                    JL=AX2**2/945.D0*(1.D0-AX2/22.D0*(1.D0-AX2/52.D0))
                ELSE
                    JL=(DSIN(AX)*(1.D0-(45.D0-105.D0/AX2)/AX2)+DCOS(AX)*(10.D0-105.D0/AX2)/AX)/AX
                ENDIF
            ELSEIF(L.EQ.5)THEN
                IF(AX.LT.1.D0)THEN
                    JL=AX2**2*AX/10395.D0*(1.D0-AX2/26.D0*(1.D0-AX2/60.D0))
                ELSE
                    JL=(DSIN(AX)*(15.D0-(420.D0-945.D0/AX2)/AX2)/AX-DCOS(AX)*(1.D0-(105.D0-945.0d0/AX2)/AX2))/AX
                ENDIF
            ELSE
                IF(AX.LT.1.D0)THEN
                    JL=AX2**3/135135.D0*(1.D0-AX2/30.D0*(1.D0-AX2/68.D0))
                ELSE
                    JL=(DSIN(AX)*(-1.D0+(210.D0-(4725.D0-10395.D0/AX2)/AX2)/AX2)+ &
                        DCOS(AX)*(-21.D0+(1260.D0-10395.D0/AX2)/AX2)/AX)/AX
                ENDIF
            ENDIF
        ELSE
            NU=0.5D0+L
            NU2=NU**2
            IF(AX.LT.1.D-40)THEN
                JL=0.D0
            ELSEIF((AX2/L).LT.5.D-1)THEN
                JL=DEXP(L*DLOG(AX/NU)-LN2+NU*ONEMLN2-(1.D0-(1.D0-3.5D0/NU2)/NU2/30.D0)/12.D0/NU) &
                   /NU*(1.D0-AX2/(4.D0*NU+4.D0)*(1.D0-AX2/(8.D0*NU+16.D0)*(1.D0-AX2/(12.D0*NU+36.D0))))
            ELSEIF((real(L,dl)**2/AX).LT.5.D-1)THEN
                BETA=AX-PID2*(L+1)
                JL=(DCOS(BETA)*(1.D0-(NU2-0.25D0)*(NU2-2.25D0)/8.D0/AX2*(1.D0-(NU2-6.25)*(NU2-12.25D0)/48.D0/AX2)) &
                   -DSIN(BETA)*(NU2-0.25D0)/2.D0/AX* (1.D0-(NU2-2.25D0)*(NU2-6.25D0)/24.D0/AX2*(1.D0-(NU2-12.25)* &
                       (NU2-20.25)/80.D0/AX2)) )/AX   
            ELSE
                L3=NU**0.325
                IF(AX .LT. NU-1.31*L3) then
                    COSB=NU/AX
                    SX = DSQRT(NU2-AX2)
                    COTB=NU/SX
                    SECB=AX/NU
                    BETA=DLOG(COSB+SX/AX)
                    COT3B=COTB**3
                    COT6B=COT3B**2
                    SEC2B=SECB**2
                    EXPTERM=( (2.D0+3.D0*SEC2B)*COT3B/24.D0 &
                       - ( (4.D0+SEC2B)*SEC2B*COT6B/16.D0 &
                       + ((16.D0-(1512.D0+(3654.D0+375.D0*SEC2B)*SEC2B)*SEC2B)*COT3B/5760.D0 &
                       + (32.D0+(288.D0+(232.D0+13.D0*SEC2B)*SEC2B)*SEC2B)*SEC2B*COT6B/128.D0/NU)*COT6B/NU) &
                       /NU)/NU
                    JL=DSQRT(COTB*COSB)/(2.D0*NU)*DEXP(-NU*BETA+NU/COTB-EXPTERM)

                !          /**************** Region 2: x >> l ****************/

                ELSEIF (AX .GT. NU+1.48*L3) then
                    COSB=NU/AX
                    SX=DSQRT(AX2-NU2)
                    COTB=NU/SX
                    SECB=AX/NU
                    BETA=DACOS(COSB)
                    COT3B=COTB**3
                    COT6B=COT3B**2
                    SEC2B=SECB**2
                    TRIGARG=NU/COTB-NU*BETA-PID4 &
                           -((2.0+3.0*SEC2B)*COT3B/24.D0  &
                           +(16.D0-(1512.D0+(3654.D0+375.D0*SEC2B)*SEC2B)*SEC2B)*COT3B*COT6B/5760.D0/NU2)/NU
                    EXPTERM=( (4.D0+sec2b)*sec2b*cot6b/16.D0 &
                           -(32.D0+(288.D0+(232.D0+13.D0*SEC2B)*SEC2B)*SEC2B)*SEC2B*COT6B**2/128.D0/NU2)/NU2
                    JL=DSQRT(COTB*COSB)/NU*DEXP(-EXPTERM)*DCOS(TRIGARG)

                !          /***************** Region 3: x near l ****************/

                ELSE
                    BETA=AX-NU
                    BETA2=BETA**2
                    SX=6.D0/AX
                    SX2=SX**2
                    SECB=SX**0.3333333333333333d0
                    SEC2B=SECB**2
                    JL=( GAMMA1*SECB + BETA*GAMMA2*SEC2B &
                          -(BETA2/18.D0-1.D0/45.D0)*BETA*SX*SECB*GAMMA1 &
                          -((BETA2-1.D0)*BETA2/36.D0+1.D0/420.D0)*SX*SEC2B*GAMMA2   &
                          +(((BETA2/1620.D0-7.D0/3240.D0)*BETA2+1.D0/648.D0)*BETA2-1.D0/8100.D0)*SX2*SECB*GAMMA1 &
                          +(((BETA2/4536.D0-1.D0/810.D0)*BETA2+19.D0/11340.D0)*BETA2-13.D0/28350.D0)*BETA*SX2*SEC2B*GAMMA2 &
                          -((((BETA2/349920.D0-1.D0/29160.D0)*BETA2+71.D0/583200.D0)*BETA2-121.D0/874800.D0)* &
                           BETA2+7939.D0/224532000.D0)*BETA*SX2*SX*SECB*GAMMA1)*DSQRT(SX)/ROOTPI12
                ENDIF
            ENDIF
        ENDIF
        IF(X.LT.0.AND.MOD(L,2).NE.0)JL=-JL
        END SUBROUTINE BJL    

    ! SUBROUTINE FOR INTEGRATION
    SUBROUTINE qtrap(func,a,b,s)
    INTEGER JMAX
    double precision a,b,func,s,EPS
    EXTERNAL func
    PARAMETER (EPS=1.e-4, JMAX=20)
!    USES trapzd
!    Returns as s the integral of the function func from a to b. The parameters EPS can be set
!    to the desired fractional accuracy and JMAX so that 2 to the power JMAX-1 is the maximum
!    allowed number of steps. Integration is performed by the trapezoidal rule.
    INTEGER j
    double precision olds
    olds=-1.e30
    do 11 j=1,JMAX
        call trapzd(func,a,b,s,j)
        if (j.gt.5) then
            if (abs(s-olds).lt.EPS*abs(olds).or. (s.eq.0..and.olds.eq.0.)) return
        endif
    olds=s
11      continue
    print*, 'too many steps in qtrap'
    END SUBROUTINE qtrap

	SUBROUTINE trapzd(func,a,b,s,n)
	INTEGER n
	double precision a,b,s,func
	EXTERNAL func
	INTEGER it,j
	double precision del,sum,tnm,x
	if (n.eq.1) then
	  s=0.5*(b-a)*(func(a)+func(b))
	else
	  it=2**(n-2)
	  tnm=it
	  del=(b-a)/tnm
	  x=a+0.5*del
	  sum=0.
	  do 11 j=1,it
	    sum=sum+func(x)
	    x=x+del
11	  continue
	  s=0.5*(s+(b-a)*sum/tnm)
	endif
	return
	END SUBROUTINE trapzd
	
	

!	
	FUNCTION assert_eq(n1,n2,string)
	CHARACTER(LEN=*), INTENT(IN) :: string
	INTEGER, INTENT(IN) :: n1,n2
	INTEGER :: assert_eq
	if (n1 == n2) then
		assert_eq=n1
	else
		write (*,*) 'nrerror: an assert_eq failed with this tag:', &
			string
		STOP 'program terminated by assert_eq2'
	end if
	END FUNCTION assert_eq
	
	FUNCTION iminloc(arr)
	INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
	INTEGER, PARAMETER :: DP = KIND(1.0D0)
	REAL(DP), DIMENSION(:), INTENT(IN) :: arr
	INTEGER(I4B), DIMENSION(1) :: imin
	INTEGER(I4B) :: iminloc
	imin=minloc(arr(:))
	iminloc=imin(1)
	END FUNCTION iminloc

	! SUBROUTINE FOR BILINIAR INTERPOLATION ===================
	SUBROUTINE bilin_interp(x1a,x2a,ya,x1,x2,y)
	    IMPLICIT NONE
	    INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
	    INTEGER, PARAMETER :: DP = KIND(1.0D0)
	    REAL(DP), DIMENSION(:), INTENT(IN) :: x1a,x2a
	    REAL(DP), DIMENSION(:,:), INTENT(IN) :: ya
	    REAL(DP), INTENT(IN) :: x1,x2
	    REAL(DP), INTENT(OUT) :: y
	    INTEGER(I4B) :: j,k,m,ndum
	    DOUBLE PRECISION t,u,xx1,xx2
	    m=assert_eq(size(x1a),size(ya,1),'polin2: m')
	    ndum=assert_eq(size(x2a),size(ya,2),'polin2: ndum')
	    xx1=x1
	    xx2=x2
	    if (x1>=x1a(m)) xx1=x1a(m) ! no estrapolation
	    if (x1<=x1a(1)) xx1=x1a(1) ! no estrapolation
	    if (x2>=x2a(ndum)) xx2=x2a(ndum) ! no estrapolation
	    if (x2<=x2a(1)) xx2=x2a(1) ! no estrapolation
	    j=iminloc(abs(xx1-x1a))
	    if (xx1-x1a(j)<0.d0) j=j-1 ! if x1a(j) is the upper limit shift
	    k=iminloc(abs(xx2-x2a))
	    if (xx2-x2a(k)<0.d0) k=k-1 ! if x2a(k) is the upper limit shift
    !	write(*,*) x1a(j),x1,x1a(j+1)
    !	write(*,*) x2a(k),x1,x2a(k+1)
	    t=(xx1-x1a(j))/(x1a(j+1)-x1a(j))
	    u=(xx2-x2a(k))/(x2a(k+1)-x2a(k))
	    y=(1.d0-t)*(1.d0-u)*ya(j,k)+t*(1.d0-u)*ya(j+1,k)+t*u*ya(j+1,k+1)+(1.d0-t)*u*ya(j,k+1)
	    return
	END SUBROUTINE bilin_interp

	! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ! Interpolation 2D-grid: in one dimension spline interpolation in the other linear :::::::::
    ! z(x,y): x axes spline interpolation =====================
    ! INPUT:
    ! x: where to compute the grid along x-axes ===============
    ! x_grid: list of x values at which z(x,y) is computed ====
    ! x_grid_size: size of the x-grid =========================
    ! z_grid_aty: values of z(x) at fixed y ===================
    ! ddz_grid_aty: d^2 z(x,y) / dx^2 at fixed y ==============
	function spline_interp(x,x_grid_size,x_grid,z_grid_aty,ddz_grid_aty) result(z)

    implicit none
    double precision x
    integer x_grid_size
    double precision, intent(in) :: x_grid(x_grid_size), z_grid_aty(x_grid_size), ddz_grid_aty(x_grid_size)
    integer xlo,xhi
    double precision z, dz
    double precision ho,a0,b0
    double precision, dimension(2) :: f_x1x2, ddf_x1x2
    integer, save :: i_last = 1

    if (x < x_grid(1)) then
        f_x1x2=z_grid_aty(1:2)
        dz = (f_x1x2(2)-f_x1x2(1))/(x_grid(2)-x_grid(1))
        z = f_x1x2(1) + dz*(x-x_grid(1))
    else if (x > x_grid(x_grid_size)) then
        !Do dodgy linear extrapolation on assumption accuracy of result won't matter
        f_x1x2=z_grid_aty(x_grid_size-1:x_grid_size)
        dz = (f_x1x2(2)-f_x1x2(1))/(x_grid(x_grid_size)-x_grid(x_grid_size-1))
        z = f_x1x2(2) + dz*(x-x_grid(x_grid_size))
    else
        ! cubic spline interpolation
        xlo=min(i_last,x_grid_size)
        do while (x_grid(xlo) > x)
            xlo=xlo-1
        end do
        do while (x_grid(xlo+1)< x)
            xlo=xlo+1
        end do
        i_last =xlo
        xhi=xlo+1

        f_x1x2=z_grid_aty(xlo:xhi)
        ddf_x1x2 = ddz_grid_aty(xlo:xhi)

        ho=x_grid(xhi)-x_grid(xlo)
        a0=(x_grid(xhi)-x)/ho
        b0=1-a0

        z = a0*f_x1x2(1)+b0*f_x1x2(2)+((a0**3-a0)*ddf_x1x2(1) &
        + (b0**3-b0)*ddf_x1x2(2))*ho**2/6
    end if

    end function spline_interp
	! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ! Interpolation 2D-grid: in one dimension spline interpolation in the other linear :::::::::
    function interp2d(x,y,x_grid_size,x_grid,y_grid_size,y_grid,z_grid,ddz_grid) result(z)


    double precision x,y,z
    integer x_grid_size,y_grid_size
    double precision, intent(in) :: x_grid(x_grid_size), y_grid(y_grid_size), z_grid(x_grid_size,y_grid_size), ddz_grid(x_grid_size,y_grid_size)
    integer ylo, yhi, iy, itf
    double precision ho,a0,b0
    double precision, dimension(4) :: f_y, ddf_y, yvec
    integer, save :: yi_last = 1

!    if(y>y_grid(y_grid_size) .or. y<y_grid(1)) then
!        write (*,*) 'interp2d ERROR: y (',y,') out interpolation grid'
!        STOP 24
!    end if
    if(y>y_grid(y_grid_size)) z=spline_interp(x,x_grid_size,x_grid,z_grid(:,y_grid_size),ddz_grid(:,y_grid_size)) ! if y > y_max -> y==y_max
    if (y<y_grid(1)) z=spline_interp(x,x_grid_size,x_grid,z_grid(:,1),ddz_grid(:,1)) ! if y < y_min -> y==y_min

    if (y_grid(y_grid_size) - y < 1.d-4) then ! if y = y_max
        z=spline_interp(x,x_grid_size,x_grid,z_grid(:,y_grid_size),ddz_grid(:,y_grid_size))
    else
        ylo=min(yi_last,y_grid_size)
        do while (y_grid(ylo) > y)
            ylo=ylo-1
        end do
        do while (y_grid(ylo+1)< y)
            ylo=ylo+1
        end do
        yi_last=ylo
        yhi=ylo+1

        if(ylo==1)then
            iy = 2
            yvec(2:4)=y_grid(ylo:yhi+1)
            do itf=ylo, yhi+1
                f_y(iy) = spline_interp(x,x_grid_size,x_grid,z_grid(:,itf),ddz_grid(:,itf))
                iy=iy+1
            end do
            call spline_double(yvec(2:4),f_y(2:4),3,ddf_y(2:4))
        elseif (yhi==y_grid_size)then
            iy = 1
            yvec(1:3)=y_grid(ylo-1:yhi)
            do itf=ylo-1, yhi
                f_y(iy) = spline_interp(x,x_grid_size,x_grid,z_grid(:,itf),ddz_grid(:,itf))
                iy=iy+1
            end do
            call spline_double(yvec(1:3),f_y(1:3),3,ddf_y(1:3))
        else
            iy = 1
            yvec(:)=y_grid(ylo-1:yhi+1)
            do itf=ylo-1, yhi+1
                f_y(iy) = spline_interp(x,x_grid_size,x_grid,z_grid(:,itf),ddz_grid(:,itf))
                iy=iy+1
            end do
            call spline_double(yvec,f_y,4,ddf_y)
        end if

        ho=yvec(3)-yvec(2)
        a0=(yvec(3)-y)/ho
        b0=(y-yvec(2))/ho

        z = a0*f_y(2)+b0*f_y(3)+((a0**3-a0)*ddf_y(2) &
        +(b0**3-b0)*ddf_y(3))*ho**2/6

    end if

    end function interp2d

	! SUBROUTINE FOR TRAPEZOIDAL INTEGRATION ==================
	SUBROUTINE trapzd_int(func,a,b,s,n)
	    INTEGER n
	    double precision a,b,s,func
	    EXTERNAL func
	    INTEGER j
	    double precision del,tnm,x
	    if (n.eq.1) then
	        s=0.5*(b-a)*(func(a)+func(b))
	    else
	        tnm=n
	        del=(b-a)/tnm
	        x=a
	        s=0.d0
	        do j=1,n
	            s=s+0.5d0*del*(func(x)+func(x+del))
	            x=x+del
	        end do
    !	    write(*,*) a,b,func(a),func(b),n,redshift,s,j,x
	    endif
	    return
	END SUBROUTINE trapzd_int


	! SUBROUTINE FOR TRAPEZOIDAL INTEGRATION ==================
	! FOR FUNCTION WHICH HAS 2 ARGUMENTS AS INPUT =============

	SUBROUTINE trapzd_int_2args(func,a,b,args,s,n)
	    INTEGER n ! number of suddivision of the interval
	    double precision, intent(in) :: args(2)
	    double precision a,b,s,func
	    EXTERNAL func
	    INTEGER j
	    double precision del,tnm,x
	    if (n.eq.1) then
	        s=0.5*(b-a)*(func(a,args(1),args(2))+func(b,args(1),args(2)))
	    else
	        tnm=n
	        del=(b-a)/tnm
	        x=a
	        s=0.d0
	        do j=1,n
	            s=s+0.5d0*del*(func(x,args(1),args(2))+func(x+del,args(1),args(2)))
	            x=x+del
	        end do
    !	    write(*,*) a,b,func(a),func(b),n,redshift,s,j,x
	    endif
	    return
	END SUBROUTINE trapzd_int_2args
	
	! SUBROUTINE FOR TRAPEZOIDAL INTEGRATION ==================
	! FOR FUNCTION WHICH HAS 1 ARGUMENT AS INPUT ==============
	
	SUBROUTINE trapzd_int_1args(func,a,b,arg,s,n)
	    INTEGER n ! number of suddivision of the interval
	    double precision, intent(in) :: arg
	    double precision a,b,s,func
	    EXTERNAL func
	    INTEGER j
	    double precision del,tnm,x
	    if (n.eq.1) then
	        s=0.5*(b-a)*(func(a,arg)+func(b,arg))
	    else
	        tnm=n
	        del=(b-a)/tnm
	        x=a
	        s=0.d0
	        do j=1,n
	            s=s+0.5d0*del*(func(x,arg)+func(x+del,arg))
	            x=x+del
	        end do
    !	    write(*,*) a,b,func(a),func(b),n,redshift,s,j,x
	    endif
	    return
	END SUBROUTINE trapzd_int_1args
    
    
    
end module utils_NC
