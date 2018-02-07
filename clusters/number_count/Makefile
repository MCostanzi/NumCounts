include ${COSMOSIS_SRC_DIR}/config/compilers.mk

all: likelihood_NC.so number_count_module.so

test:
	@echo "Alas,  Tinker MF  has no tests"

clean:
	rm -f *.so *.o *.mod
	rm -rf *.so.dSYM/

export RECOMBINATION ?= recfast
CAMBDIR=${COSMOSIS_SRC_DIR}/cosmosis-standard-library/boltzmann/camb/camb_Jan15
CAMBLIB=libcamb_$(RECOMBINATION).a

LAPACKL =  -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a -Wl,--end-group -L $(INTEL_HOME)/linux/compiler/lib/intel64 -liomp5 -lpthread -lm

OBJS=quadpack.f90 legendre_polynomial.f90 interface_tools_NC.f90 utils_NC.f90 read_data.f90 Phi_of_m.f90 hmf_tinker.f90 compute_dVdzdOmega.f90 compute_number_count.f90 likelihood_NC.f90

likelihood_NC.so: $(OBJS)
	$(FC) $(FFLAGS) -shared $(OBJS) -o $@ $(INTERFACE) ${COSMOSIS_SRC_DIR}/cosmosis/datablock/libcosmosis_fortran.so -I$(CAMBDIR)  -L$(CAMBDIR) -lcamb_$(RECOMBINATION) $(LAPACKL)

OBJS2=quadpack.f90 legendre_polynomial.f90 interface_tools_NC.f90 utils_NC.f90 read_data.f90 Phi_of_m.f90 hmf_tinker.f90 compute_dVdzdOmega.f90 compute_number_count.f90 number_count_module.f90

number_count_module.so: $(OBJS2)
	$(FC) $(FFLAGS) -shared $(OBJS2) -o $@ $(INTERFACE) ${COSMOSIS_SRC_DIR}/cosmosis/datablock/libcosmosis_fortran.so





