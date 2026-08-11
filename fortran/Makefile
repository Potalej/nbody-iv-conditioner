FC = gfortran
FFLAGS = -O3
INSTALL_DIR ?= ../lib

SRCS = utils.f90 conditioners.f90
LIB_FLAGS = -I../lib -L../lib -lconditioners_mod -lutils

lib:
	mkdir -p ${INSTALL_DIR}

	$(FC) $(FFLAGS) -c utils.f90 -J${INSTALL_DIR}
	ar rcs libutils.a utils.o

	$(FC) $(FFLAGS) -lutils -c conditioners.f90 -J${INSTALL_DIR}
	ar rcs libconditioners_mod.a conditioners.o

	$(FC) $(FFLAGS) -c ${SRCS} -J${INSTALL_DIR}
	
ifneq (${INSTALL_DIR},./)
	mv *.o *.a ${INSTALL_DIR}/
endif

test:
	$(FC) $(FFLAGS) -o test ./tests/* ${LIB_FLAGS}

clean:
	rm -rf conditioners test lib/ ../lib *.o *.a *.mod