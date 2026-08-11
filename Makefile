FC = gfortran
FFLAGS = -O3
INSTALL_DIR ?= ../lib

.PHONY: all lib test clean

all: lib

lib: $(INSTALL_DIR)/libconditioners.a $(INSTALL_DIR)/libutils.a

$(INSTALL_DIR):
	mkdir -p $@

utils.o: utils.f90 | $(INSTALL_DIR)
	$(FC) $(FFLAGS) -c $< -J$(INSTALL_DIR)

$(INSTALL_DIR)/libutils.a: utils.o
	ar rcs $@ $<

conditioners.o: conditioners.f90 utils.o | $(INSTALL_DIR)
	$(FC) $(FFLAGS) -c $< -I$(INSTALL_DIR) -J$(INSTALL_DIR)

$(INSTALL_DIR)/libconditioners.a: conditioners.o
	ar rcs $@ $<

test: lib
	$(FC) $(FFLAGS) -o test1 ./tests/test1.f90 \
		-I$(INSTALL_DIR) \
		-L$(INSTALL_DIR) \
		-lconditioners -lutils

clean:
	rm -rf test1 test2 test3 *.o *.a *.mod
	rm -rf $(INSTALL_DIR)