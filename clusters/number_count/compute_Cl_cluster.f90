! Computing halo number density in a given redshift and observable bin:


    MODULE compute_Cl_cluster
    USE legendre_polynomial
    USE hmf_tinker
    USE compute_dVdzdOmega
    USE interface_tools_NC
    USE compute_number_count
    USE Phi_of_m
    IMPLICIT none
    integer, parameter :: num_of_elle_Cl=200 ! Number of multipoles to compute for C_ell
    double precision, dimension(:,:,:), allocatable :: int_bessel_Pk ! where to store \int dk k^2*Bessel_z1*Bessel_z2*Pk_z1_z2
    ! =========================================================
    ! Some useful Global variables ============================
    double precision raggio1,raggio2 ! For integration Bessel_l_z1*Bessel_l_z2*P(k,z1,z2)
    double precision rez1,rez2 ! For integration Bessel_l_z1*Bessel_l_z2*P(k,z1,z2)
    integer ell_Cl ! For integration Bessel_l_z1*Bessel_l_z2*P(k,z1,z2)
    !$OMP THREADPRIVATE(raggio1,raggio2,rez1,rez2,ell_Cl)
    double precision kappa_value,z_min_Cl,z_max_Cl
    contains



    SUBROUTINE compute_Cl(PK,NC)
    implicit none
    type(pk_settings) :: PK
    type(number_counts) :: NC
    integer nm, nz, nL, ell ! counter
    integer nz1, nz2, nL1, nL2 ! counter
    double precision T0, T1 ! computation time
    integer now0(3),now1(3)
    double precision d_logM,log_M_max,log_M_min,log_mass
    double precision epsabs, epsrel, abserr ! for integral for j_l(rk)*r^2
    integer neval, ier, key ! for integral for j_l(rk)*r^2
    epsabs = 0.0E+00 ! for integral for j_l(rk)*r^2
    epsrel = 0.0001E+00 ! for integral for j_l(rk)*r^2
    key = 6 ! for integral for j_l(rk)*r^2

    call CPU_time(T0)
    call itime(now0)

	! Read the Fitting Parameters for the Tinker HMF ========
	call init_tinker_param

	! If Delta_m is defined get fitting parameters for Tinker hmf for Delta_m
	! otherwise the fitting parameters are computed each time n(M,z) is called
	! for a different redshift for Delta_m=Delta_c/Omega_m(z)
    if (Delta_m .gt. 1.d-3) call Get_fit_param(0.d0)

    ! =========================================================
    ! ALLOCATE POWER SPECTRUM QUANTITIES ======================
    call allocate_pow_spec(PK)
    ! =========================================================
    ! ALLOCATE ARRAYS TO STORE \SIGMA^2(R,z) and derivative ===
    ! this is do to speed up the code =========================
    n_log_M_step=100 ! number of points where compute \sigma^2(R,z) and d\sigma^2(R,z)/dR
    call allocate_Ln_sigma

    log_M_min=log(1.d12) ! this must be equal or lower than log_m_lo in dndz(z_reds)
    log_M_max=log(1.d17) ! this must be equal or larger than log_m_up in dndz(z_reds)
    d_logM=(log_M_max-log_M_min)/(n_log_M_step-1.d0) 
    ! Compute and store \sigma^2(R,z) and d\sigma^2(R,z)/dR ===
    !$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(STATIC), PRIVATE(nm,nz,log_mass)
    do nm=1,n_log_M_step
        log_mass=log_M_min+d_logM*(nm-1)
        Ln_R_array(nm)=log(M_to_R(exp(log_mass)))
        do nz=1,num_of_z
            Ln_sigma2Rz(nm,nz)=log(sigma_r(M_to_R(exp(log_mass)),redshifts_pk(nz)))
            Ln_dsigma2Rz(nm,nz)=log(-dsig2dr(M_to_R(exp(log_mass)),redshifts_pk(nz))) ! NB disgma2/dR is negative
        end do
!        write(*,*) nm,log_mass
    end do ! end redshift loop
    !$OMP END PARALLEL DO
    ! Compute second derivative for interpolation =============
    do nz=1,PK%num_z
        call spline(Ln_R_array,Ln_sigma2Rz(:,nz),n_log_M_step,cllo,clhi,dd_Ln_sigma2Rz(:,nz))  ! dd Ln(sigma2(R,z))
        call spline(Ln_R_array,Ln_dsigma2Rz(:,nz),n_log_M_step,cllo,clhi,dd_Ln_dsigma2Rz(:,nz))  ! dd Ln(sigma2(R,z))
    end do

!    ! =========================================================
!    ! Compute integral of k^2 j_l(kr(z1))*j_l(kr(z2))* ========
!    ! sqrt(P(k,z1)*P(k,z2)) ===================================
!    allocate(int_bessel_Pk(NC%num_of_z_bin,NC%num_of_z_bin,num_of_elle_Cl))
!    !$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(STATIC), PRIVATE(nz1,nz2,nL)
!    do nz1=1,NC%num_of_z_bin
!        rez1=(NC%z_min_array(nz1)+NC%z_max_array(nz1))*0.5d0
!        raggio1=dcom(rez1)
!        do nz2=nz1,NC%num_of_z_bin
!            rez2=(NC%z_min_array(nz2)+NC%z_max_array(nz2))*0.5d0
!            raggio2=dcom(rez2)
!            do nL=2,num_of_elle_Cl
!                ell_Cl=nL
!                call qag(d_bessel_Pk,log(1.d-4),log(1.d3), epsabs, epsrel, key, int_bessel_Pk(nz1,nz2,ell_Cl), abserr, neval, ier )
!!                write(*,*) rez1,rez2,ell_Cl,int_bessel_Pk(nz1,nz2,ell_Cl),d_bessel_Pk(log(0.01d0))
!!                call BJL(ell_Cl,0.01d0*raggio1,raggio2)
!!                write(*,'(6e16.8)') 0.01d0*raggio1,raggio2**2.d0*Pk_at_z(0.01d0,rez1),d_bessel_Pk(log(0.01d0)),Pk_at_z(0.01d0,rez1)
!            end do
!        end do
!    end do
!    !$OMP END PARALLEL DO
!    ! =========================================================
!!!    stop 24

    call CPU_time(T1)
    call itime(now1)
    write(*,*)'computation time: ',T1-T0,'time end ',now1(2),'m ',now1(3),'s; time start',now0(2),'m',now0(3),'s'


!    Delta_Omega = NC%D_Omega ! Survey Area in steradians
    
!    open(64,file='Cl_output/bjl_kr_and_Pkz.dat') ! Test Buzzard
!    do ell=2,num_of_elle_Cl
!        do nz1=1,NC%num_of_z_bin ! Start loop over redshift bin 1
!            call BJL(ell,1.d-4*dcom(NC%z_min_array(nz1)),raggio1)
!            call BJL(ell,1.d1*dcom(NC%z_min_array(nz1)),raggio2)
!            write(64,'(10e16.8)') ell*1.d0,NC%z_min_array(nz1),1.d-4*dcom(NC%z_min_array(nz1)),raggio1,sqrt(Pk_at_z(1.d-4,NC%z_min_array(nz1))),1.d1*dcom(NC%z_min_array(nz1)),raggio2,sqrt(Pk_at_z(1.d1,NC%z_min_array(nz1)))
!        end do
!        write(64,*) ''
!        write(64,*) ''
!    end do
!    stop 24
    open(80,file='outputs_random/Cl_b1.dat') ! Test Buzzard
    NC%n_Li_zj=0.d0
    NC%n_Li_zj_var=0.d0

    do nz1=1,NC%num_of_z_bin ! Start loop over redshift bin 1
        do nL1=1,NC%num_of_L_bin ! Start loop over Lambda bin 1
!!            ! se uso lambda
            task=1 ! NC
            NC%n_Li_zj(nL1,nz1)=n_Lambda_z(nL1,nz1,NC)
!!            task=2 ! mean mass
!!            NC%mean_m_Li_zj(nL1,nz1)=n_Lambda_z(nL1,nz1,NC)/NC%n_Li_zj(nL1,nz1) 
!            task=3 ! b(M,z)*N_ij
!            NC%n_Li_zj_var(nL1,nz1)=n_Lambda_z(nL1,nz1,NC)
!            write(*,'(5e16.8)') NC%z_min_array(nz1),NC%z_max_array(nz1), 2.d13*10.d0**((nL1-1)*(0.2d0)),2.d13*10.d0**((nL1)*(0.2d0)),NC%n_Li_zj(nL1,nz1)
        end do
    end do
!    stop 24
    task=1 ! in verita non necessario perche dopo non chiamo piu' n_lambda_z
    do nz1=1,NC%num_of_z_bin ! Start loop over redshift bin 1
!        do nz2=nz1,NC%num_of_z_bin ! Start loop over redshift bin 2
            do nL1=1,NC%num_of_L_bin ! Start loop over Lambda bin 1
!                do nL2=nL1,NC%num_of_L_bin ! Start loop over Lambda bin 2
                    do ell=2,num_of_elle_Cl
!                        write(*,'(5e16.8)') ell*1.d0,2.d0/pi_value*NC%n_Li_zj_var(nL1,nz1)*NC%n_Li_zj_var(nL1,nz1)*int_bessel_Pk(nz1,nz1,ell)/NC%n_Li_zj(nL1,nz1)/NC%n_Li_zj(nL1,nz1),C_ell(ell,nz1,nL1,NC)/NC%n_Li_zj(nL1,nz1)**2.d0
                        write(80,'(5e16.8)') ell*1.d0,C_ell_limb(ell,nz1,nL1,NC)/(NC%n_Li_zj(nL1,nz1)/Delta_Omega)**2.d0!,C_ell(ell,nz1,nL1,NC)/(NC%n_Li_zj(nL1,nz1)/Delta_Omega)**2.d0
!                        write(80,'(5e16.8)') ell*1.d0,C_ell_limb(ell,nz1,nL1,NC)/(NC%n_Li_zj(nL1,nz1))**2.d0*ell*(ell+1.d0)/2.d0/pi_value
                    end do
!                end do
            end do
!        end do
    end do


    call CPU_time(T1)
    call itime(now1)
    write(*,*)'computation time: ',T1-T0,'time end ',now1(2),'m ',now1(3),'s; time start',now0(2),'m',now0(3),'s'
    NC%sigma8_NC=sqrt(sigma2R_at_z(8.d0,0.d0))
    write(*,*) 'sigma8',NC%sigma8_NC,'Omega_m',Omega_m
    write(*,*) 'slope',slope,'inter',inter
    write(*,*) 'delta_a0',delta_a0,'delta_b0',delta_b0
    write(*,*) 'frac_p',frac_p,'epsi',epsi
    ! =========================================================

!    deallocate(int_bessel_Pk)

    call deallocate_pow_spec
    call deallocate_Ln_sigma

    END SUBROUTINE compute_Cl

	! Integrand of ============================================
	! \int dLnk k^3 j_l(r1*k)*j_l(r2*k) P(k,z1,z2) ============
	function d_bessel_Pk(lnk)
	double precision d_bessel_Pk,kerre1,kerre2,bess_ell_kr1,bess_ell_kr2,lnk
	integer elle
	
	elle=ell_Cl
	kerre1=exp(lnk)*raggio1
	kerre2=exp(lnk)*raggio2
	call BJL(elle,kerre1,bess_ell_kr1)
	call BJL(elle,kerre2,bess_ell_kr2)
	d_bessel_Pk=exp(lnk)**3.d0*bess_ell_kr1*bess_ell_kr2*sqrt(Pk_at_z(exp(lnk),rez1)*Pk_at_z(exp(lnk),rez2))
!	d_bessel_Pk=bess_ell_kr1**2.d0*sqrt(Pk_at_z(exp(lnk),rez1)*Pk_at_z(exp(lnk),rez2))
	end function d_bessel_Pk
    ! =========================================================

    ! C_ell for Limber approximation
    function C_ell_limb(ell,n_z,n_L,NC)
    implicit none
    type(number_counts) :: NC
    double precision C_ell_limb,z_lo,z_up,z_reds
    integer ell,n_z,n_L
    double precision epsabs, epsrel, abserr
    integer neval, ier
    epsabs = 0.0E+00
    epsrel = 0.001E+00

    z_lo=NC%z_min_array(n_z)
    z_up=NC%z_max_array(n_z)

    ! Delta Ln_Lambda
    LnLambda_min = NC%LnLambda_min_array(n_L) ! Global Variable used in Phi_M(halo_mass)
    LnLambda_max = NC%LnLambda_max_array(n_L) ! Global Variable used in Phi_M(halo_mass)
    
    ell_Cl= ell! Assign the Global Variable
    call qags (d_Cl_limb, z_lo,z_up, epsabs, epsrel, C_ell_limb, abserr, neval, ier )
    return
!    C_ell_limb=C_ell_limb*4.d0*pi_value
    end function C_ell_limb
    ! =========================================================

    ! integrand of C_ell for Limber approximation
    function d_Cl_limb(z_red)
    double precision d_Cl_limb,z_red,int_nbP
    
!    d_Cl=f_inv_e(z_red)*c_light*derivs_com_dis(z_red)**2.d0/dcom(z_red)**2.d0*b_eff_at_z(z_red)*Pk_at_z((ell_Cl+0.5d0)/dcom(z_red),z_red)
    task=1
    int_nbP=dndz(z_red)
    task=1
    d_Cl_limb=int_nbP**2.d0*Pk_at_z((ell_Cl+0.5d0)/dcom(z_red),z_red)*derivs_com_dis(z_red)!*Delta_Omega
    end function d_Cl_limb

    ! Effective bias in Richness bin ==========================
    ! LnLambda_min LnLambda_max and z =========================
    ! b_eff(\Delta Lambda_i,z) = (int_0^inf dM n(M,z)* b(M,z) =
    ! * Phi(M)) / (int_0^inf dM n(M,z)*Phi(M)) ================
    function b_eff_at_z(zed)
    implicit none
    double precision b_eff_at_z,zed
    double precision int_nbP,int_nP

    task=1
    int_nP=dndz(zed)
    task=3
    int_nbP=dndz(zed)
    b_eff_at_z=int_nbP/int_nP
    task=1
    return
    end function b_eff_at_z
    ! =========================================================


    ! C_ell*N(Delta_z)^2 ======================================
    function C_ell(ell,n_z,n_L,NC)
    implicit none
    type(number_counts) :: NC
    double precision C_ell,lnk_lo,lnk_up
    integer ell,n_z,n_L
    double precision epsabs, epsrel, abserr
    integer neval, ier
    epsabs = 0.0E+00
    epsrel = 0.001E+00

    lnk_lo=log(1.d-4)
    lnk_up=log(10.d0)
    
    z_min_Cl=NC%z_min_array(n_z) ! Global Variable used in W_k
    z_max_Cl=NC%z_max_array(n_z) ! Global Variable used in W_k

    ! Delta Ln_Lambda
    LnLambda_min = NC%LnLambda_min_array(n_L) ! Global Variable used in Phi_M(halo_mass)
    LnLambda_max = NC%LnLambda_max_array(n_L) ! Global Variable used in Phi_M(halo_mass)
    
    ell_Cl= ell! Assign the Global Variable
    call qags (d_Cl, lnk_lo, lnk_up, epsabs, epsrel, C_ell, abserr, neval, ier )
    C_ell=C_ell*2.d0/pi_value
    return
    end function C_ell
    ! =========================================================

    ! integrand of C_ell ======================================
    ! dln k k^3 [W(k)_i]^2 ====================================
    function d_Cl(ln_k)
    double precision d_Cl,ln_k

    d_Cl=exp(ln_k)**3.d0*W_k(ln_k)**2.d0
    end function d_Cl


    !  ========================================================
    ! W(k) = \int_Delta_z dz dV/dz/dOmega * Delta_Omega * =====
    ! [/int dM n(M,z) b(M,z) \int_Delta_lambda d lambda * =====
    ! P(lambda|M,z)] * j_l(kr(z)) * sqrt(P(k,z)) ==============
    function W_k(ln_k)
    implicit none 
    double precision W_k,ln_k
    double precision epsabs, epsrel, abserr
    integer neval, ier, key
    epsabs = 0.0E+00
    epsrel = 0.01E+00
    key = 5 ! for integral for j_l(rk)*r^2


    kappa_value=exp(ln_k) ! Assign the Global Variable

    call qag (integrand_W_k, z_min_Cl,z_max_Cl, epsabs, epsrel, key, W_k, abserr, neval, ier )
    return
    end function W_k
    ! =========================================================

	! integrand of W(k)_i =====================================
	FUNCTION integrand_W_k(z_reds)
    double precision z_reds,int_nbP,bessel_l_kr,kerre,integrand_W_k
    task=3
    int_nbP=dndz(z_red)
    task=1
    
    kerre=kappa_value*dcom(z_reds)
    call BJL(ell_Cl,kerre,bessel_l_kr)
    integrand_W_k = int_nbP*derivs_com_dis(z_red)*bessel_l_kr*sqrt(Pk_at_z(kappa_value,z_reds))!*Delta_Omega
!    write(*,'(2e18.9)') exp(Ln_halo_mass),Phi_M(exp(Ln_halo_mass))
    return
	END FUNCTION integrand_W_k
	! =========================================================







END MODULE compute_Cl_cluster

