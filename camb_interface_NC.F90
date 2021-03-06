module camb_interface_tools
	use camb
	use cosmosis_modules
	implicit none

	integer :: standard_lmax = 1200
	integer, parameter :: CAMB_MODE_ALL = 1
	integer, parameter :: CAMB_MODE_CMB = 2
	integer, parameter :: CAMB_MODE_BG  = 3
	integer, parameter :: CAMB_MODE_THERMAL  = 4

	real(8) :: linear_zmin=0.0, linear_zmax=4.0
	integer :: linear_nz = 401
	real(8) :: linear_kmax=50.0 ! mod_costanzi

	integer :: k_eta_max_scalar = 2400
	logical :: do_lensing, do_nonlinear, do_tensors
	logical :: do_nonlinear_Pk ! mod_costanzi
	integer :: sterile_neutrino = 0
	real(dl) :: delta_neff = 0.0
	real(dl) :: sterile_mass_fraction = 0.0
	real(dl) :: cmb_output_scale = 7.4311e12

	real(dl), parameter :: default_yhe = 0.24
	real(dl), parameter :: default_cs2de = 1.0
	real(dl), parameter :: default_r = 0.0
	real(dl), parameter :: default_nrun = 0.0
	real(dl), parameter :: default_w = -1.0
	real(dl), parameter :: default_wa = 0.0
	real(dl), parameter :: default_pivot_scalar = 0.05
	integer,  parameter :: default_massive_nu = 0
	integer,  parameter :: default_sterile_neutrinos = 0


	contains


	function camb_comoving_sound_horizon() result(rsdrag)
		use ModelParams
		use Precision
		use ThermoData, only : z_drag
		implicit none
		real(dl) ::  adrag, atol, rsdrag
		real(dl), external :: rombint
		integer error

		adrag = 1.0d0/(1.0d0+z_drag)
		atol = 1e-6
		rsdrag = rombint(dsound_da,1d-8,adrag,atol)
	end function camb_comoving_sound_horizon


	function camb_shift_parameter(params) result(shift_parameter)
		type(cambparams) :: params
		real(dl) :: omega_m, ombh2, omdmh2, zstar, shift_parameter
		real(dl), parameter :: c_km_per_s = 299792.458

		omega_m = params%omegac + params%omegab + params%omegan

         ombh2 = CP%omegab*(CP%h0/100.0d0)**2
         omdmh2 = (CP%omegac+CP%omegan)*(CP%h0/100.0d0)**2

    !From Hu & Sugiyama (via modules.f90)
		zstar =  1048*(1+0.00124*ombh2**(-0.738))*(1+ &
			(0.0783*ombh2**(-0.238)/(1+39.5*ombh2**0.763)) * &
			(omdmh2+ombh2)**(0.560/(1+21.1*ombh2**1.81)))

		shift_parameter = sqrt(omega_m) * params%H0 / c_km_per_s * &
		&   (1+zstar)*AngularDiameterDistance(zstar)

	end function

	function camb_initial_setup(block, mode, fixed_mode) result(status)
		integer default_lmax
		integer(c_size_t) :: block
		integer status
		character(64) :: mode_name=""
		integer :: mode
		integer, optional :: fixed_mode
		integer::  use_tabulated_w_int
		default_lmax = standard_lmax
		status=0
		! There are currently three camb modes - "background", "cmb", and "all"
		! This code may get called with a fixed mode, or with not in which case
		! we read from file
		! First in the fixed mode case, we just use that as the output mode
		if (present(fixed_mode)) then
			mode=fixed_mode
		else
			!Otherwise read from ini file

			status = datablock_get_string(block, option_section, "mode", mode_name)
			if (trim(mode_name) == "background") then
				mode=CAMB_MODE_BG
			else if (trim(mode_name) == "cmb") then
				mode=CAMB_MODE_CMB
			else if (trim(mode_name) == "all") then
				mode=CAMB_MODE_ALL
			else if (trim(mode_name) == "thermal") then
				mode=CAMB_MODE_THERMAL
			else
				write(*,*) "You need to specify a mode to use the camb module you chose."
				write(*,*) "In the camb section of your ini file, please specify one of:"
				write(*,*) "mode=background  ; For background quantities like D_A(z) only"
				write(*,*) "mode=cmb         ; For background + cmb power spectra"
				write(*,*) "mode=all         ; For background + cmb + linear matter power spectra"
				write(*,*) "mode=thermal     ; For background + thermal history params"
				write(*,*) ""
				write(*,*) "We found error status: ", status
				write(*,*) "And mode=", mode_name
				write(*,*) "Quitting now."
				stop 1
			endif
		endif

		status = 0

		!We do not use the CMB lmax if only using the background mode
		if (mode .ne. CAMB_MODE_BG) then
			status = status + datablock_get_int_default(block, option_section, "lmax", default_lmax, standard_lmax)
			status = status + datablock_get_int_default(block, option_section, "k_eta_max_scalar", 2*standard_lmax, k_eta_max_scalar)
		endif

		!We can always set an optional feedback level,
		!which defaults to zero (silent)
		status = status + datablock_get_int_default(block, option_section, "feedback", 0, FeedbackLevel)
		status = status + datablock_get_logical_default(block, option_section, "use_tabulated_w", .false., use_tabulated_w)
		status = status + datablock_get_logical_default(block, option_section, "do_tensors", .false., do_tensors)

		if (mode == CAMB_MODE_ALL) then
			status = status + datablock_get_double_default(block, option_section,"zmin", linear_zmin, linear_zmin)
			status = status + datablock_get_double_default(block, option_section,"zmax", linear_zmax, linear_zmax)
			status = status + datablock_get_int_default(block, option_section,"nz", linear_nz, linear_nz)
			status = status + datablock_get_double_default(block, option_section,"kmax", linear_kmax, linear_kmax) ! mod_costanzi
		endif

		status = status + datablock_get_logical_default(block, option_section, "do_nonlinear", .false. , do_nonlinear)
		status = status + datablock_get_logical_default(block, option_section, "do_lensing", .false. , do_lensing)
		status = status + datablock_get_logical_default(block, option_section, "do_nonlinear_Pk", .false., do_nonlinear_Pk) ! mod_costanzi

		!Error check
		if (status .ne. 0) then
			write(*,*) "Problem setting some options for camb. Status code =  ", status
			return
		endif


 		if (do_lensing) then
 			status = status + datablock_get_string(block, option_section, "high_ell_template", highL_unlensed_cl_template)
 			if ((status .ne. 0 ) .or. trim(highL_unlensed_cl_template)=="") then
 				status = 1
 				write(*,*) "If you set do_lensing=1 then you also need to set"
 				write(*,*) "the parameter high_ell_template to the complete path"
 				write(*,*) "to the file HighLExtrapTemplate_lenspotentialCls.dat"
 				write(*,*) "which comes with CAMB - i.e. in your ini file camb section, put:"
 				write(*,*) "high_ell_template = /path/to/cosmosis/src/standard-library/boltzmann/camb/HighLExtrapTemplate_lenspotentialCls.dat"
 			elseif (.not. FileExists(trim(highL_unlensed_cl_template))) then
 				status = 2
 				write(*,*) "You set the parameter high_ell_template in the ini file to the value:"
 				write(*,*) trim(highL_unlensed_cl_template)
 				write(*,*) "But I could not find a file there.  You need to include the full"
 				write(*,*) "path to that file as that parameter, i.e.:"
 				write(*,*) "high_ell_template = /path/to/cosmosis/src/standard-library/boltzmann/camb/HighLExtrapTemplate_lenspotentialCls.dat"
 			endif
 		endif
 		


		!If noisy, report relevant params
		if (FeedbackLevel .gt. 0) then
			write(*,*) "camb mode  = ", mode
			if (mode .ne. CAMB_MODE_BG) write(*,*) "camb cmb_lmax = ", standard_lmax
			write(*,*) "camb FeedbackLevel = ", FeedbackLevel
			if (status .ne. 0) write(*,*) "Setup status: ", status
		endif
	end function camb_initial_setup

	function camb_interface_set_params(block, params, background_only) result(status)
		integer (c_int) :: status
		integer (c_size_t) :: block
		logical, optional :: background_only
		logical :: perturbations
		type(CambParams) :: params
		integer :: sterile_neutrino
		real(8), dimension(:), allocatable :: w_array, a_array
		character(*), parameter :: cosmo = cosmological_parameters_section
		

		perturbations = .true.
		if (present(background_only)) perturbations = .not. background_only

	
		call CAMB_SetDefParams(params)
		status = 0

        status = status + datablock_get_double(block, cosmo, "omega_b", params%omegab)
        status = status + datablock_get_double(block, cosmo, "omega_c", params%omegac)
        status = status + datablock_get_double(block, cosmo, "omega_lambda", params%omegav)
        status = status + datablock_get_double(block, cosmo, "omega_nu", params%omegan)
		status = status + datablock_get_double(block, cosmo, "omega_k", params%omegak)
		status = status + datablock_get_double(block, cosmo, "hubble", params%H0)
		
		if (perturbations) then
			status = status + datablock_get_double(block, cosmo, "n_s",     params%initpower%an(1))
            status = status + datablock_get_double_default(block, cosmo, "k_s", default_pivot_scalar, params%initpower%k_0_scalar)
			status = status + datablock_get_double_default(block, cosmo, "A_s",0.d0,params%initpower%ScalarPowerAmp(1))
			if (params%initpower%ScalarPowerAmp(1) <= 0.d0) then ! mod_costanzi vary in the chain Ln(10^10A_s)
			    status = status + datablock_get_double(block, cosmo, "LnA_s",     params%initpower%ScalarPowerAmp(1)) ! mod_costanzi
			    params%initpower%ScalarPowerAmp(1)=dexp(params%initpower%ScalarPowerAmp(1))/1.0d10 ! mod_costanzi
			end if ! mod_costanzi
			status = status + datablock_get_double(block, cosmo, "tau", params%Reion%optical_depth)
			status = status + datablock_get_double_default(block, cosmo, "R_T", default_r, params%initpower%rat(1))

			status = status + datablock_get_double_default(block, cosmo, "n_run", default_nrun, params%initpower%n_run(1))
			if (params%initpower%rat(1) .ne. 0) then
				status = status + datablock_get_double(block, cosmo, "n_T", params%initpower%ant(1))
			endif
		endif

		!Neutrinos

		status = status + datablock_get_double_default(block, cosmo, "cs2_de", default_cs2de, cs2_lam)
		status = status + datablock_get_double_default(block, cosmo, "yhe", default_yhe, params%yhe)

		status = status + datablock_get_double_default(block, cosmo, "massless_nu", params%Num_Nu_massless, params%Num_Nu_massless) ! mod_costanzi this parameter has to be read even if omega_nu=0
		if (params%omegan .ne. 0) then
			status = status + datablock_get_int_default(block, cosmo, "sterile_neutrino", default_sterile_neutrinos, sterile_neutrino)
!			status = status + datablock_get_double_default(block, cosmo, "massless_nu", params%Num_Nu_massless, params%Num_Nu_massless) ! mod_costanzi
			status = status + datablock_get_int_default(block, cosmo, "massive_nu", default_massive_nu, params%Num_Nu_massive)

			!  We have coded for two massive neturino scenarios so far:
			!  sterile neutrinos, and a single massive neutrino.
			params%nu_mass_numbers = 0.d0 ! mod_costanzi
			if (sterile_neutrino > 0) then
				status = status + datablock_get_double(block, cosmo, "delta_neff", delta_neff)
				status = status + datablock_get_double(block, cosmo, "sterile_mass_fraction", sterile_mass_fraction)
				params%share_delta_neff = .false.
				params%Num_Nu_massless = 2.0307
				params%Nu_mass_eigenstates = 2
				params%Num_Nu_massive = 2
				params%nu_mass_degeneracies(1) = 1.0153
				params%nu_mass_degeneracies(2) = delta_neff
				params%nu_mass_fractions(1) = (1.0 - sterile_mass_fraction) 
				params%nu_mass_fractions(2) = sterile_mass_fraction
			elseif (params%Num_Nu_massive == 1) then
				params%Nu_mass_eigenstates = 1
				params%nu_mass_numbers(1) = 1
				params%Nu_mass_fractions(1) = 1.0
				params%share_delta_neff = .true.
			! start mod_costanzi
			elseif (params%Num_Nu_massive == 3) then 
				params%Nu_mass_eigenstates = 1
				params%nu_mass_numbers(1) = 3 !array of the integer number of physical neutrinos per eigenstate
				params%Nu_mass_fractions(1) = 1.0
				params%share_delta_neff = .true.
			! end mod_costanzi
			elseif (params%Num_Nu_massive == 0) then
				write(*,*) 'You need massive_nu>0 to have any omega_nu!=0'
				status=1
				return
			else
				stop "Sorry - we have not coded up the neutrino scenario your parameters implied"
			endif
		endif


		call setcgammappf()

		! tabulated dark energy EoS
		if (use_tabulated_w) then
			status = status + datablock_get_double_array_1d(block, de_equation_of_state_section, "w", w_array, nw_ppf)
			status = status + datablock_get_double_array_1d(block, de_equation_of_state_section, "a", a_array, nw_ppf)
			if (nw_ppf .gt. nwmax) then
				write(*,*) "The size of the w(a) table was too large ", nw_ppf, nwmax
				status=nw_ppf
				return
			endif
			w_ppf(1:nw_ppf) = w_array(1:nw_ppf)
			a_ppf(1:nw_ppf) = dlog(a_array(1:nw_ppf))  !a is stored as log(a)
			deallocate(w_array, a_array)
			call setddwa()
			call interpolrde()
		else
			status = status + datablock_get_double_default(block, cosmo, "w", -1.0D0, w_lam)
			status = status + datablock_get_double_default(block, cosmo, "wa",  0.0D0, wa_ppf)
!			write(*,*) 'w_0',w_lam,'w_a',wa_ppf
			if (w_lam+wa_ppf .gt. 0) then
				write(*,*) "Unphysical w_0 + w_a = ", w_lam, " + ", wa_ppf, " = ", w_lam+wa_ppf, " > 0"
				status = 1
			endif
		endif	

		params%wantTransfer = .true.
		params%transfer%high_precision = .TRUE. ! mod_costanzi
		params%transfer%kmax = linear_kmax !50.0 ! mod_costanzi
		params%wantTensors = (params%initpower%rat(1) .ne. 0.0) .or. do_tensors
!		write(*,*) 'params%transfer%k_per_logint',params%transfer%k_per_logint ! mod_costanzi
!		write(*,*) 'params%transfer%high_precision',params%transfer%high_precision ! mod_costanzi
!		write(*,*) 'params%transfer%kmax',params%transfer%kmax ! mod_costanzi
        params%Max_l=standard_lmax
        params%Max_eta_k=2*standard_lmax

        params%DoLensing = do_lensing

		!Set nonlinear behaviour
		if (do_nonlinear) then
		    if (do_nonlinear_Pk) then ! mod_costanzi
		        params%NonLinear=3 ! mod_costanzi
		    else ! mod_costanzi
			    params%NonLinear=2 ! mod_costanzi
			end if ! mod_costanzi
!			params%NonLinear=2 ! mod_costanzi Original
		else
			params%NonLinear=0
		endif

	
	
		!Some extras and modifications 
		params%want_zdrag = .true.
		params%want_zstar = .true.
		params%reion%use_optical_depth = .true.
		params%reion%delta_redshift = 0.5

		use_spline_template=params%DoLensing
		params%AccurateReionization = .true.
        params%Transfer%PK_num_redshifts = 1
        params%Transfer%PK_redshifts = 0

	end function
	
	function camb_interface_setup_zrange(params) result(status)
		integer(cosmosis_status) :: status
		type(CambParams) :: params
		real(8) :: zmin, zmax, dz
		integer nz, i

		zmin = linear_zmin
		zmax = linear_zmax
		nz = linear_nz

		dz=(zmax-zmin)/(nz-1.0)
		params%transfer%num_redshifts =  nz
        params%Transfer%PK_num_redshifts = nz

		if (nz .gt. max_transfer_redshifts) then
			write(*,*) "Requested too many redshifts for CAMB to handle: ", nz, " = (", zmax, " - ", zmin, ") / ", dz, " + 1"
			status = 1
		endif
		
        do i=1,params%transfer%num_redshifts
			params%transfer%redshifts(nz-i+1)  = zmin + dz*(i-1)
	        params%transfer%pk_redshifts(nz-i+1)  = zmin + dz*(i-1)
    	enddo


    	call Transfer_SortAndIndexRedshifts(params%transfer)
		status = 0
	end function



	function camb_interface_save_cls(block) result(status)
	
		integer (cosmosis_block) :: block
		integer (cosmosis_status) :: status
	
		integer, parameter :: input_set = 1
		real  :: cls(2:standard_lmax,1:4)
		real(8)  :: cls_double(2:standard_lmax,1:4), cls_phi(2:standard_lmax)
		integer  :: ell(2:standard_lmax), l
		logical, parameter :: switch_polarization_convention = .false.	
	
		status = 0
		call CAMB_GetCls(cls, standard_lmax, input_set, switch_polarization_convention)
		cls_double(:,1:4) = cls * 7.4311e12  !cmb output scale
	    do l=2,standard_lmax
			ell(l) = l
		enddo

		if (do_lensing) then
		    do l=2,standard_lmax
				cls_phi(l) = Cl_scalar(l, input_set,  C_phi) * (l+1.0) / ((l*1.0)**3 * twopi)
			enddo
		endif
	
		status = status + datablock_put_int_array_1d(block, cmb_cl_section, "ELL", ell)
		status = status + datablock_put_double_array_1d(block, cmb_cl_section, "TT", cls_double(:,1))
		status = status + datablock_put_double_array_1d(block, cmb_cl_section, "EE", cls_double(:,2))
		status = status + datablock_put_double_array_1d(block, cmb_cl_section, "BB", cls_double(:,3))
		status = status + datablock_put_double_array_1d(block, cmb_cl_section, "TE", cls_double(:,4))
		if (do_lensing) then
			status = status + datablock_put_double_array_1d(block, cmb_cl_section, "PP", cls_phi)
		endif
	
		if (status .ne. 0) then
			write(*,*) "Failed to save cmb!."
			return
		endif
	end function

	function camb_interface_save_sigma8(block) result(status)
		!Save sigma8 at z=0 to the cosmological parameters section of the file
		integer (cosmosis_block) :: block
		integer (cosmosis_status) :: status
		real(8) :: sigma8
		real(8), parameter :: radius8 = 8.0_8
		integer nz

		!Ask camb for sigma8
		status = 0
		sigma8=0.0
		call Transfer_Get_sigma8(MT,radius8)
		
		!It gives us the array sigma8(z).
		!We want the entry for z=0
		nz = CP%Transfer%num_redshifts
		sigma8 = MT%sigma_8(nz,1)

		!Save sigma8
		status = status + datablock_put_double(block, cosmological_parameters_section, "SIGMA_8", sigma8)
		return
	end function
	
	function camb_interface_save_transfer(block) result(status)
		integer (cosmosis_block) :: block
		integer (cosmosis_status) :: status
		Type(MatterPowerData) :: PK
		integer nz, nk, iz, ik
		real(8), allocatable, dimension(:) :: k, z
		real(8), allocatable, dimension(:,:) :: P, T
		real(8), allocatable, dimension(:,:) :: P_NL ! mod_costanzi

		call Transfer_GetMatterPowerData(MT, PK, 1,var1 = transfer_nonu,var2=transfer_nonu) ! mod_costanzi
		! In this way we compute P_cdm+baryon; change var1=var2 for different power spectra

		nz = CP%Transfer%num_redshifts
		nk = MT%num_q_trans

		allocate(k(nk))
		allocate(z(nz))
		allocate(P(nk,nz))
		allocate(T(nk,nz))

		do ik=1,nk
			k(ik) = MT%TransferData(Transfer_kh,ik,1)
		enddo

		do iz=1,nz
			z(iz) = CP%Transfer%Redshifts(nz-iz+1)
		enddo

		do ik=1,nk
			do iz=1,nz
			    P(ik,iz) = exp(PK%matpower(ik,nz-iz+1)) ! mod_costanzi to speed up a bit we don t need to interpolate again
!				P(ik,iz) = MatterPowerData_k(PK, k(ik), nz-iz+1) ! mod_costanzi original
				T(ik,iz) = MT%TransferData(Transfer_cdm,ik,nz-iz+1)
!				write(*,*) P(ik,iz),exp(PK%matpower(ik,nz-iz+1))
			enddo
		enddo
		
!		! start mod_costanzi
!		open(64,file='ps_k_z.dat')
!		do iz=1,nz
!		    do ik=1,nk
!		        write(64,'(3e16.8)') z(iz),k(ik),P(ik,iz)
!		    end do
!    		write(64,*) ''
!    		write(64,*) ''
!		end do
!		! end mod_costanzi

		status = datablock_put_double_grid(block, matter_power_lin_section, &
        	"k_h", k, "z", z, "P_k", P)

		if (status .ne. 0) then
			write(*,*) "Failed to save transfer function in CAMB."
		endif

		status = datablock_put_double_grid(block, linear_cdm_transfer_section, &
        	"k_h", k, "z", z, "delta_cdm", T) ! mod_costanzi prima c era P al posto di T

		if (status .ne. 0) then
			write(*,*) "Failed to save matter power in CAMB."
		endif

!		! start mod_costanzi
		if (do_nonlinear_Pk) then
		    allocate(P_NL(nk,nz))
		    call MatterPowerdata_MakeNonlinear(PK)
			do ik=1,nk
			    do iz=1,nz
			    P_NL(ik,iz) = exp(PK%matpower(ik,nz-iz+1))
			    end do
			end do

		    status = datablock_put_double_grid(block, matter_power_nl_section, &
        	"k_h", k, "z", z, "PNL_k", P_NL)

    		if (status .ne. 0) then
    			write(*,*) "Failed to save transfer function in CAMB."
    		endif
    		deallocate(P_NL)
		end if
!		! end mod_costanzi

		deallocate(k, z, P, T)

	end function

	
	function camb_interface_save_da(params, block, save_density, save_thermal) result(status)
		integer (cosmosis_block) :: block
		type(CambParams) :: params
		integer (c_int) :: status
		logical, optional :: save_density, save_thermal
		logical :: density, thermal
		real(8), dimension(:), allocatable :: distance, z, rho
		character(*), parameter :: dist = distances_section
		integer nz, i
		

		! Rho as given by the code is actually 8 * pi * G * rho / c**2 , and it is measured in (Mpc)**-2
		! There really isn't a sensible set of units to do this in, so let's just use kg/m**3
		! c**2 / (8 pi G) = 5.35895884e24 kg/m
		! 1 Mpc = 3.08568025e24 m
		real(8), parameter :: mpc_in_m = 3.08568025e22
		real(8), parameter :: c2_8piG_kgm = 5.35895884e25
		real(8), parameter :: rho_units = c2_8piG_kgm / (mpc_in_m**2)
		real(8) :: shift, rs_zdrag

		density = .true.
		if (present(save_density)) density = save_density
		thermal = .true.
		if (present(save_thermal)) thermal = save_thermal

		status = 0


		nz = params%transfer%num_redshifts
		allocate(distance(nz))
		allocate(z(nz))
		!if (density) allocate(rho(nz))

		do i=1,nz
			z(i) = params%transfer%redshifts(i)
			distance(i) = AngularDiameterDistance(z(i))
			!if (density) rho(i) = MT%TransferData(Transfer_rho_tot,1,i) * rho_units
		enddo
		

		shift = camb_shift_parameter(params)
		status = status + datablock_put_double(block, dist, "CMBSHIFT", shift)


		if (thermal) then
			status = status + datablock_put_double(block, dist, &
				"AGE", ThermoDerivedParams( derived_Age ))
			status = status + datablock_put_metadata(block, dist, "AGE", "unit", "Gyr")

			status = status + datablock_put_double(block, dist, &
				"RS_ZDRAG", ThermoDerivedParams( derived_rdrag ))
			status = status + datablock_put_metadata(block, dist, "RS_ZDRAG", "unit", "Mpc")

			!There is an 
			status = status + datablock_put_double(block, dist, &
				"THETASTAR", ThermoDerivedParams( derived_thetastar ))
			status = status + datablock_put_metadata(block, dist, "THETASTAR", "unit", "100 radian")

			status = status + datablock_put_double(block, dist, &
				"ZDRAG", ThermoDerivedParams( derived_zdrag ))


			status = status + datablock_put_double(block, dist, &
				"K_D", ThermoDerivedParams( derived_kD ))
			status = status + datablock_put_metadata(block, dist, "K_D", "unit", "1/Mpc")

			status = status + datablock_put_double(block, dist, &
				"THETA_D", ThermoDerivedParams( derived_thetaD ))
			status = status + datablock_put_metadata(block, dist, "THETA_D", "unit", "100 radian")

			status = status + datablock_put_double(block, dist, &
				"Z_EQUALITY", ThermoDerivedParams( derived_zEQ ))

			status = status + datablock_put_double(block, dist, &
				"K_EQUALITY", ThermoDerivedParams( derived_keq ))
			status = status + datablock_put_metadata(block, dist, "K_EQUALITY", "unit", "1/Mpc")


			status = status + datablock_put_double(block, dist, &
				"THETA_EQUALITY", ThermoDerivedParams( derived_thetaEQ ))
			status = status + datablock_put_metadata(block, dist, "THETA_EQUALITY", "unit", "100 radian")

			status = status + datablock_put_double(block, dist, &
				"THETA_RS_EQUALITY", ThermoDerivedParams( derived_theta_rs_EQ ))
			status = status + datablock_put_metadata(block, dist, "THETA_RS_EQUALITY", "unit", "100 radian")

			status = status + datablock_put_double(block, dist, &
				"DA_STAR", ThermoDerivedParams( derived_DAstar ))
			status = status + datablock_put_metadata(block, dist, "DA_STAR", "unit", "Gpc")

			status = status + datablock_put_double(block, dist, &
				"R_STAR", ThermoDerivedParams( derived_rstar ))
			status = status + datablock_put_metadata(block, dist, "R_STAR", "unit", "Mpc")

			status = status + datablock_put_double(block, dist, &
				"ZSTAR", ThermoDerivedParams( derived_zstar ))

			status = status + datablock_put_double(block, dist, &
				"CHISTAR", ComovingRadialDistance(ThermoDerivedParams( derived_zstar )))
			status = status + datablock_put_metadata(block, dist, "CHISTAR", "unit", "Myr")
		else
			status = status + datablock_put_double(block, dist, &
				"AGE", DeltaPhysicalTimeGyr(0.0_dl,1.0_dl))
			status = status + datablock_put_metadata(block, dist, "AGE", "unit", "Gyr")
		endif


		status = status + datablock_put_double_array_1d(block, dist, "Z", z)
		status = status + datablock_put_double_array_1d(block, dist, "D_A", distance)
		status = status + datablock_put_metadata(block, dist, "D_A", "unit", "Mpc")

		distance = distance * (1+z) !Convert to D_M
		status = status + datablock_put_double_array_1d(block, dist, "D_M", distance)
		status = status + datablock_put_metadata(block, dist, "D_M", "unit", "Mpc")

		distance = distance * (1+z) !Convert to D_L
		status = status + datablock_put_double_array_1d(block, dist, "D_L", distance)
		status = status + datablock_put_metadata(block, dist, "D_L", "unit", "Mpc")

		distance = 5*log10(distance)+25 !Convert to distance modulus
		! The distance is already the dimensionful one, so we do not
		! multiply be c/H0
		status = status + datablock_put_double_array_1d(block, dist, "MU", distance)

		! Save H(z)
		do i=1,nz
			distance(i) = HofZ(z(i))
		enddo
		status = status + datablock_put_double_array_1d(block, dist, "H", distance)
		status = status + datablock_put_metadata(block, dist, "H", "unit", "Mpc/c")

		!if (density) then
		!	status = status + datablock_put_double_array_1d(block, dist, "RHO", rho)
		!	status = status + datablock_put_metadata(block, dist, "RHO", "unit", "KG/M^3")
		!endif


		status = status + datablock_put_int(block, dist, "NZ", nz)

		!And finally save a
		z = 1.0/(1+z)
		status = status + datablock_put_double_array_1d(block, dist, "A", z)


		if (status .ne. 0) then
			write(*,*) "Failed to write redshift-distance column data in block section."
		endif
		
		deallocate(distance)
		deallocate(z)
		!if (density) deallocate(rho)
		
	end function
	
end module camb_interface_tools



